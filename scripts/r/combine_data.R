####### 
## This script does the following:
##
##  1. USFS Data + Models: Read/Process, create log-log OLS re-fits in USFS data
##  2. Inventory Data + Model fits: Load data from the _processed folder
##  3. Create combined data + compute coverage: 
##          - create combined data set
##          - create training/holdout from USFS
##          - compute coverage for all 4 combinations
##          - add results to _processed
##
## Note that the USFS model fit steps happen twice, once regionally and once pooling all of CA.

#####################
### 1. USFS Data
#####################

###########################################################################
## INPUT   data/_raw/TS3_Raw_tree_data.csv
##         data/_raw/TS1_Regional_information.csv
## OUTPUT  data/_processed/usfs_prepped.rds
##           $trees       one row per tree, incl. region + holdout randomization
##           $eq_regional one log-log equation per species x region
##           $eq_pooled   one log-log equation per species (regions pooled)
##           $params      the settings used, for provenance
## =============================================================================

library(tidyverse)

## ---- settings ---------------------------------------------------------------
## Overridable by a driver script: each is only set if not already defined.
if (!exists("SEED"))          SEED          <- 212
if (!exists("PCT_TRAIN"))     PCT_TRAIN     <- 0.90  # fraction used to FIT equations
if (!exists("MIN_N_PER_SPP")) MIN_N_PER_SPP <- 10    # min train trees to fit a species
if (!exists("CA_REGIONS"))    CA_REGIONS    <- c("InlEmp", "InlVal", "NoCalC",
                                                 "SacVal", "SoCalC")

CM_PER_IN <- 2.54
M_TO_FT   <- 3.281

dir.create("data/_processed", showWarnings = FALSE, recursive = TRUE)

## ---- genus helper -----------------------------------------------------------
## Matches the inventory pipeline's handling of hybrid names ("x Cupressocyparis").
get_genus <- function(name) {
  name <- str_squish(name)
  if_else(str_starts(name, regex("x\\s", ignore_case = TRUE)),
          word(str_remove(name, regex("^x\\s+", ignore_case = TRUE)), 1),
          word(name, 1))
}

## ---- read + clean -----------------------------------------------------------
## check.names = TRUE turns "DBH (cm)" -> DBH..cm. and "TreeHt (m)" -> TreeHt..m.
regions <- read.csv("data/_raw/TS1_Regional_information.csv", check.names = TRUE) |>
  transmute(region      = Region.Code,
            region_name = Region.Name,
            city        = City,
            state       = State) |>
  as_tibble()

ts3 <- read.csv("data/_raw/TS3_Raw_tree_data.csv", check.names = TRUE) |>
  as_tibble() |>
  filter(Region %in% CA_REGIONS,
         !is.na(DBH..cm.), !is.na(TreeHt..m.),
         DBH..cm. > 0, TreeHt..m. > 0) |>
  rename(dbh_cm = DBH..cm., height_m = TreeHt..m., region = Region) |>
  mutate(
    ScientificName = str_squish(ScientificName),
    ## NOTE: TreeID is only unique WITHIN a region (3941 distinct TreeID vs 4049
    ## distinct region+TreeID in the CA subset), so the tree key is composite.
    tree_id = paste(region, TreeID, sep = "-")
  )

## ---- one row per tree -------------------------------------------------------
## The raw file has a few repeated measures; keep one so a tree cannot land in
## both the training and the holdout partition. (~3% of rows.)
set.seed(SEED)
ts3 <- ts3 |>
  group_by(tree_id) |>
  slice_sample(n = 1) |>
  ungroup()

## ---- holdout randomization, stratified within species -----------------------
## u_holdout is the draw itself (kept so the split is reproducible/auditable);
## split is the resulting label. Ranking within species gives an exact
## PCT_TRAIN proportion per species rather than a binomial one.
set.seed(SEED)
usfs_trees <- ts3 |>
  mutate(u_holdout = runif(n())) |>
  mutate(
    split = if_else(rank(u_holdout, ties.method = "first") <= PCT_TRAIN * n(),
                    "train", "holdout"),
    .by = ScientificName
  ) |>
  left_join(regions, by = "region") |>
  transmute(
    tree_id,
    ScientificName,
    spp   = ScientificName,          # name the inventory model expects
    genus = get_genus(ScientificName),
    region, region_name, city, state,
    dbh_cm, height_m,
    dbh_inch       = dbh_cm   / CM_PER_IN,
    tree_height_ft = height_m * M_TO_FT,
    u_holdout,
    split = factor(split, levels = c("train", "holdout")),
    data_source = "USFS"
  )

## ---- fit the log-log equations on TRAIN only --------------------------------
fit_loglog <- function(df) {
  df |>
    mutate(ok = map_lgl(data, ~ nrow(.x) >= MIN_N_PER_SPP)) |>
    filter(ok) |>
    mutate(m = map(data, ~ lm(log(height_m) ~ log(dbh_cm), data = .x))) |>
    mutate(
      n                   = map_int(data, nrow),
      intercept           = map_dbl(m, ~ coef(.x)[1]),
      slope               = map_dbl(m, ~ coef(.x)[2]),
      intercept_se        = map_dbl(m, ~ sqrt(vcov(.x)[1, 1])),
      slope_se            = map_dbl(m, ~ sqrt(vcov(.x)[2, 2])),
      intercept_slope_cov = map_dbl(m, ~ vcov(.x)[1, 2]),
      sigma               = map_dbl(m, ~ summary(.x)$sigma)
    ) |>
    select(-data, -m, -ok)
}

usfs_train <- usfs_trees |> filter(split == "train")

## one equation per species x region
usfs_eq_regional <- usfs_train |>
  nest(.by = c(region, ScientificName)) |>
  fit_loglog() |>
  mutate(spp = ScientificName, .after = ScientificName)

## one equation per species, regions pooled (this is what the plots use)
usfs_eq_pooled <- usfs_train |>
  nest(.by = ScientificName) |>
  fit_loglog() |>
  mutate(spp = ScientificName, region = NA_character_, .after = ScientificName)

## ---- save -------------------------------------------------------------------
usfs <- list(
  trees       = usfs_trees,
  eq_regional = usfs_eq_regional,
  eq_pooled   = usfs_eq_pooled,
  params      = list(seed = SEED, pct_train = PCT_TRAIN,
                     min_n_per_spp = MIN_N_PER_SPP, ca_regions = CA_REGIONS)
)

saveRDS(usfs, "data/_processed/usfs_prepped.rds")

message(sprintf(
  "01_prep_usfs: %d trees (%d train / %d holdout), %d species, %d pooled eqs, %d regional eqs",
  nrow(usfs_trees),
  sum(usfs_trees$split == "train"), sum(usfs_trees$split == "holdout"),
  n_distinct(usfs_trees$ScientificName),
  nrow(usfs_eq_pooled), nrow(usfs_eq_regional)
))

#####################
### 2. Inventory
#####################
## INPUT   data/_processed/wca-partitions.RData     (wca_train, wca_val)
##         rslt/dbh-ht/_dbh-ht-model-obj.rds      ## update RUN_ID if new run
##         rslt/dbh-ht/dbh-ht-coefs.xlsx          ## update RUN_ID if new run
## OUTPUT  data/_processed/inventory_prepped.rds
##           $trees      one row per tree, with split
##           $eq         species-level coefficient table (used for species overlap)
##           $sigma      population residual sd
##           $model_path path to the merMod, re-read by the third section to combine below.
##           $params
## =============================================================================



## ---- settings ---------------------------------------------------------------
if (!exists("SEED")) SEED <- 212

## dbh-ht results were renamed from the dated run-id directory to rslt/dbh-ht/
RUN_ID        <- "dbh-ht"

PARTITION_PATH <- "data/_processed/wca-partitions.RData"


MODEL_PATH <- file.path("rslt", RUN_ID, "_dbh-ht-model-obj.rds") ## got from Cami, did not rerun- but IDS match
COEFS_PATH <- file.path("rslt", RUN_ID, "dbh-ht-coefs.xlsx")

CM_PER_IN <- 2.54
M_TO_FT   <- 3.281

dir.create("data/_processed", showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(PARTITION_PATH), file.exists(MODEL_PATH),
          file.exists(COEFS_PATH))

## ---- genus helper (same rule as 01) -----------------------------------------
get_genus <- function(name) {
  name <- str_squish(name)
  if_else(str_starts(name, regex("x\\s", ignore_case = TRUE)),
          word(str_remove(name, regex("^x\\s+", ignore_case = TRUE)), 1),
          word(name, 1))
}

## ---- read partitions and tag the split --------------------------------------
load(PARTITION_PATH)   # wca_train, wca_val

inv_raw <- bind_rows(
  wca_train |> mutate(split = "train"),  ## nomenclature used in 2.6, 3.4
  wca_val   |> mutate(split = "holdout")
)
rm(wca_train, wca_val)

## ---- clean + one row per tree -----------------------------------------------
## The inventory has repeated measures per tree_id across survey dates; keep one
## row per tree so a tree contributes a single point to coverage.
set.seed(SEED)
inv_trees <- inv_raw |>
  filter(dbh_height_is_anom == FALSE,
         !is.na(dbh), !is.na(height), dbh > 0, height > 0) |>
  group_by(tree_id) |>
  slice_sample(n = 1) |>
  ungroup() |>
  transmute(
    tree_id        = as.character(tree_id),
    ScientificName = spp,
    spp,
    genus          = coalesce(genus, get_genus(spp)),
    ## the inventory has no USFS region; lat/lon are kept in case you ever want
    ## to assign one, but region is deliberately NA here. Could replace with region 
    ## via  spatial join if we had region polygons!
    region         = NA_character_,
    region_name    = NA_character_,
    city           = NA_character_,
    state          = "CA",
    any_of(c("lat", "lon")),   # absent in the public release; passthrough only
    dbh_cm         = dbh,
    height_m       = height,
    dbh_inch       = dbh    / CM_PER_IN,
    tree_height_ft = height * M_TO_FT,
    u_holdout      = NA_real_,   # split is inherited, not drawn here
    split          = factor(split, levels = c("train", "holdout")),
    data_source    = "Inventory"
  )

## ---- fitted model + extracted coefficients ----------------------------------
inv_model <- readRDS(MODEL_PATH)

population <- readxl::read_xlsx(COEFS_PATH, sheet = "population")
inv_sigma  <- population$estimate[population$term == "sigma"]

inv_eq <- readxl::read_xlsx(COEFS_PATH, sheet = "species") |>
  rename(ScientificName = spp) |>
  mutate(spp = ScientificName, .after = ScientificName)

## ---- save -------------------------------------------------------------------
## The merMod itself is not copied into the rds (it is large and re-read
## from MODEL_PATH), only the path is recorded.
inv <- list(
  trees      = inv_trees,
  eq         = inv_eq,
  sigma      = inv_sigma,
  model_path = MODEL_PATH,
  params     = list(seed = SEED, partition_path = PARTITION_PATH,
                    coefs_path = COEFS_PATH, run_id = RUN_ID)
)

saveRDS(inv, "data/_processed/inventory_prepped.rds")

message(sprintf(
  "02_prep_inventory: %d trees (%d train / %d holdout), %d species, %d species eqs",
  nrow(inv_trees),
  sum(inv_trees$split == "train"), sum(inv_trees$split == "holdout"),
  n_distinct(inv_trees$ScientificName), nrow(inv_eq)
))

#####################
### 3. Combine + test/train split + compute coverage
#####################
## Produces TWO combined tables:
##   combined_regional  USFS equations kept region-by-region. A USFS tree is
##                      scored by its own region's equation; an inventory tree
##                      (no region) is scored by every regional equation, so it
##                      appears once per region. These are used to generate the 
##                      region plots in the zipfile in appendix.
##   combined_pooled    one equation per species (regions pooled) -> exactly one
##                      row per tree per equation source. This is the table the
##                      plots in the paper use, and the only one where the joint
##                      "status" (contained by both / one / neither) is defined.
##
## Region columns in both tables:
##   region_data  region the TREE came from      (NA for inventory)
##   region_eq    region the EQUATION came from  (NA for pooled/inventory eqs)
##   region       coalesce(region_eq, region_data) -- the operative region
##
##
## INPUT   data/_processed/usfs_prepped.rds
##         data/_processed/inventory_prepped.rds
## OUTPUT  data/_processed/combined_regional.rds
##         data/_processed/combined_pooled.rds
##         data/_processed/line_df.rds
##         data/_processed/common_species.rds
## =============================================================================

## merTools is deliberately NOT attached: it pulls in MASS, whose select()
## masks dplyr::select(). It is called namespaced below.

## ---- settings ---------------------------------------------------------------
TCRIT     <- qt(0.975, df = Inf)   # 1.96, for the USFS delta-method PIs
N_SIMS    <- 1000                  # merTools sims for the tree-level PIs
N_SIMS_LINE <- 10000               # merTools sims for the "smooth" plotting lines
GRID_N    <- 200                   # points per species in the plotting lines
SEED      <- 22626

## Which partitions to carry into the combined tables. "holdout" alone is much
## faster (predictInterval is the bottleneck); keeping both makes the `split`
## column informative so downstream code can choose.
SPLITS_TO_SCORE    <- c("train", "holdout")
## Cap on inventory TRAINING trees per species (the inventory train partition is
## ~1.7M rows). Holdout rows are never subsampled. Set to Inf to keep everything.
INV_MAX_TRAIN_PER_SPP <- 200

inv_model <- readRDS(inv$model_path)

usfs_trees <- usfs$trees |> rename(region_data = region)
inv_trees  <- inv$trees  |> rename(region_data = region)

usfs_eq_regional <- usfs$eq_regional |> rename(region_eq = region)
usfs_eq_pooled   <- usfs$eq_pooled   |> rename(region_eq = region)

## ---- restrict to splits + subsample inventory training ----------------------
set.seed(SEED)
usfs_trees <- usfs_trees |> filter(split %in% SPLITS_TO_SCORE)
inv_trees <- inv_trees |>
  filter(split %in% SPLITS_TO_SCORE) |>
  (\(d) bind_rows(
    d |> filter(split == "holdout"),
    d |> filter(split == "train") |>
      slice_sample(n = min(INV_MAX_TRAIN_PER_SPP, .Machine$integer.max),
                   by = ScientificName)
  ))()

## ---- species present in all four objects ------------------------------------
common_species <- Reduce(intersect, list(
  unique(usfs_trees$ScientificName),
  unique(inv_trees$ScientificName),
  unique(usfs_eq_pooled$ScientificName),
  unique(inv$eq$ScientificName)
))

usfs_trees <- usfs_trees |> filter(ScientificName %in% common_species)
inv_trees  <- inv_trees  |> filter(ScientificName %in% common_species)

message(sprintf("03_combine: %d species in common; %d USFS trees, %d inventory trees",
                length(common_species), nrow(usfs_trees), nrow(inv_trees)))

## ---- the two scoring engines ------------------------------------------------

## USFS equations: delta-method PI on the log-log line, back-transformed to feet.
score_usfs <- function(trees, eq, by) {
  trees |>
    inner_join(eq, by = by, relationship = "many-to-many") |>
    mutate(
      log_dbh  = log(dbh_inch * CM_PER_IN),
      log_fit  = intercept + slope * log_dbh,
      line_var = intercept_se^2 + 2 * log_dbh * intercept_slope_cov +
        log_dbh^2 * slope_se^2,
      log_se   = sqrt(sigma^2 + line_var),
      fit = exp(log_fit)                 * M_TO_FT,
      lo  = exp(log_fit - TCRIT * log_se) * M_TO_FT,
      hi  = exp(log_fit + TCRIT * log_se) * M_TO_FT,
      eq_source = "USFS"
    ) |>
    select(-log_dbh, -log_fit, -line_var, -log_se,
           -intercept, -slope, -intercept_se, -slope_se,
           -intercept_slope_cov, -sigma, -n)
}

## Inventory equations: merTools::predictInterval on the fitted lmer.
## Trees the model never saw (species or genus) are dropped.
score_inventory <- function(trees, model, n_sims = N_SIMS, seed = SEED) {
  model_spp <- levels(model@frame$spp)
  model_gen <- if ("genus" %in% names(model@frame)) levels(model@frame$genus) else NULL
  
  nd <- trees |> filter(spp %in% model_spp)
  if (!is.null(model_gen)) nd <- nd |> filter(genus %in% model_gen)
  nd <- nd |> mutate(dbh = dbh_inch * CM_PER_IN)
  
  if (nrow(nd) == 0) stop("no rows left for the inventory model")
  
  pi <- merTools::predictInterval(
    merMod = model, newdata = as.data.frame(nd), which = "full", level = 0.95,
    n.sims = n_sims, stat = "median", type = "linear.prediction",
    include.resid.var = TRUE, seed = seed
  )
  
  nd |>
    bind_cols(pi) |>
    mutate(
      fit = exp(fit) * M_TO_FT,
      lo  = exp(lwr) * M_TO_FT,
      hi  = exp(upr) * M_TO_FT,
      region_eq = NA_character_,   # the inventory model has no regional structure
      eq_source = "Inventory"
    ) |>
    select(-lwr, -upr, -dbh)
}

## ---- inventory-equation blocks (computed once, reused by both tables) --------
usfs_inv <- score_inventory(usfs_trees, inv_model)   # USFS data,      Inventory eq
inv_inv  <- score_inventory(inv_trees,  inv_model)   # Inventory data, Inventory eq

## ---- USFS-equation blocks ---------------------------------------------------
## regional: USFS tree matched to its OWN region's equation
usfs_usfs_reg <- score_usfs(
  usfs_trees, usfs_eq_regional,
  by = c("ScientificName", "spp", "region_data" = "region_eq")
) |>
  mutate(region_eq = region_data)

## regional: inventory trees have no region, so they are scored by every
## regional equation -> one row per tree per region
inv_usfs_reg <- score_usfs(
  inv_trees, usfs_eq_regional, by = c("ScientificName", "spp")
)

## pooled: one equation per species (usfs_eq_pooled carries region_eq = NA, so
## the region columns still line up when the blocks are stacked)
usfs_usfs_pool <- score_usfs(usfs_trees, usfs_eq_pooled,
                             by = c("ScientificName", "spp"))
inv_usfs_pool  <- score_usfs(inv_trees,  usfs_eq_pooled,
                             by = c("ScientificName", "spp"))

## ---- assemble ---------------------------------------------------------------
finalise <- function(df) {
  df |>
    filter(tree_height_ft > 0, dbh_inch > 0) |>
    mutate(
      region      = coalesce(region_eq, region_data),
      contained   = tree_height_ft >= lo & tree_height_ft <= hi,
      data_source = factor(data_source, levels = c("Inventory", "USFS")),
      eq_source   = factor(eq_source,   levels = c("Inventory", "USFS")),
      combo       = paste(as.character(data_source), as.character(eq_source), sep = "_"),
      data_model  = factor(
        paste0("Data: ", data_source, ", Model: ", eq_source),
        levels = c("Data: USFS, Model: USFS",
                   "Data: Inventory, Model: USFS",
                   "Data: USFS, Model: Inventory",
                   "Data: Inventory, Model: Inventory")
      ),
      ## a row is out-of-sample if the tree was held out of its own source's fit
      ## OR if it is being scored by the other source's equation entirely
      is_holdout = split == "holdout" | data_source != eq_source
    ) |>
    select(combo, data_model, data_source, eq_source,
           ScientificName, spp, genus, tree_id,
           region, region_data, region_eq,
           dbh_inch, tree_height_ft, fit, lo, hi, contained,
           u_holdout, split, is_holdout)
}

combined_regional <- bind_rows(
  usfs_usfs_reg, inv_usfs_reg, usfs_inv, inv_inv
) |> finalise()

combined_pooled <- bind_rows(
  usfs_usfs_pool, inv_usfs_pool, usfs_inv, inv_inv
) |>
  finalise() |>
  ## joint status: for a given tree, is it inside the USFS PI, the inventory PI,
  ## both, or neither? Only well defined in the pooled table (1 row per tree
  ## per equation source).
  mutate(
    in_USFS      = any(contained[eq_source == "USFS"]),
    in_Inventory = any(contained[eq_source == "Inventory"]),
    .by = c(data_source, tree_id)
  ) |>
  mutate(
    status = factor(
      case_when(
        in_USFS &  in_Inventory ~ "both",
        in_USFS & !in_Inventory ~ "usfs_only",
        !in_USFS &  in_Inventory ~ "inventory_only",
        !in_USFS & !in_Inventory ~ "neither"
      ),
      levels = c("both", "usfs_only", "inventory_only", "neither")
    )
  )

## ---- smooth plotting lines (pooled equations only) --------------------------
## USFS line: deterministic from the pooled coefficients
usfs_line <- combined_pooled |>
  filter(eq_source == "USFS") |>
  summarise(lo = min(dbh_inch), hi = max(dbh_inch), .by = ScientificName) |>
  reframe(dbh_inch = seq(lo, hi, length.out = GRID_N), .by = ScientificName) |>
  inner_join(usfs_eq_pooled, by = "ScientificName") |>
  mutate(
    log_dbh  = log(dbh_inch * CM_PER_IN),
    log_fit  = intercept + slope * log_dbh,
    line_var = intercept_se^2 + 2 * log_dbh * intercept_slope_cov +
      log_dbh^2 * slope_se^2,
    log_se   = sqrt(sigma^2 + line_var),
    fit = exp(log_fit)                 * M_TO_FT,
    lo  = exp(log_fit - TCRIT * log_se) * M_TO_FT,
    hi  = exp(log_fit + TCRIT * log_se) * M_TO_FT,
    eq_source = "USFS"
  ) |>
  select(ScientificName, dbh_inch, fit, lo, hi, eq_source)

## Inventory line: predictInterval on the same grid
inv_grid <- combined_pooled |>
  filter(eq_source == "Inventory") |>
  summarise(lo = min(dbh_inch), hi = max(dbh_inch), .by = ScientificName) |>
  reframe(dbh_inch = seq(lo, hi, length.out = GRID_N), .by = ScientificName) |>
  mutate(spp   = ScientificName,
         genus = word(ScientificName, 1),
         dbh   = dbh_inch * CM_PER_IN)

inv_pi <- merTools::predictInterval(
  merMod = inv_model, newdata = as.data.frame(inv_grid), which = "full",
  level = 0.95, n.sims = N_SIMS_LINE, stat = "median",
  type = "linear.prediction", include.resid.var = TRUE, seed = SEED
)

inv_line <- inv_grid |>
  bind_cols(inv_pi) |>
  mutate(fit = exp(fit) * M_TO_FT,
         lo  = exp(lwr) * M_TO_FT,
         hi  = exp(upr) * M_TO_FT,
         eq_source = "Inventory") |>
  select(ScientificName, dbh_inch, fit, lo, hi, eq_source)

line_df <- bind_rows(usfs_line, inv_line)

## ---- save -------------------------------------------------------------------
saveRDS(common_species,    "data/_processed/common_species.rds")
saveRDS(combined_regional, "data/_processed/combined_regional.rds")
saveRDS(combined_pooled,   "data/_processed/combined_pooled.rds")
saveRDS(line_df,           "data/_processed/line_df.rds")

message(sprintf(
  "03_combine: regional = %d rows, pooled = %d rows, line_df = %d rows",
  nrow(combined_regional), nrow(combined_pooled), nrow(line_df)
))

message("combine_data.R: data combined and prepped.")
