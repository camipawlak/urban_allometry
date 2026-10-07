## =============================================================================
## regional_plots_by_spp.R
##
## Per-species scatter + line plots, faceted by USFS region.
##
## Facet logic:
##   - a USFS tree belongs to ONE region and appears in that facet only
##   - an inventory tree has no region, so it is duplicated into every facet
##   - the USFS line differs per facet (one regional equation each)
##   - the inventory line is identical in every facet (the lmer has no regional
##     structure) and is drawn repeated, as a fixed reference
##
## The number of facets for a species = the number of regions in which that
## species has a USFS regional equation, i.e. where it cleared MIN_N_PER_SPP
## in combine_data.r. That is usually FEWER than 5, and for some species is 1.
##
##
## INPUT   data/_processed/combined_regional.rds
##         data/_processed/usfs_prepped.rds      (regional coefficients + trees)
##         data/_processed/line_df.rds           (pooled inventory line)
##         data/_processed/common_species.rds
## OUTPUT  _figures/by-region/our-refits/*.png
##         _figures/by-region/our-refits-and-mcpherson/*.png
## =============================================================================

library(tidyverse)
library(ggnewscale)

## ---- settings ---------------------------------------------------------------
CM_PER_IN <- 2.54
M_TO_FT   <- 3.281
TCRIT     <- qt(0.975, df = Inf)
GRID_N    <- 200

## Only draw a facet if the region has at least this many holdout USFS trees.
## Set to 0 to keep every region that has an equation.
MIN_USFS_HOLDOUT_PER_REGION <- 1

dir.create("_figures/by-region/our-refits", showWarnings = FALSE, recursive = TRUE)

## ---- theming ----------------------------------------------------------------
src_colors  <- c(Inventory = "#1f78b4", USFS = "#e31a1c")
src_shapes  <- c(Inventory = 16,        USFS = 17)
src_lty     <- c(Inventory = 1,         USFS = 2)
status_cols <- c(both           = "grey",
                 neither        = "black",
                 inventory_only = "#1f78b4",
                 usfs_only      = "#e31a1c")

## ---- data -------------------------------------------------------------------
combined  <- readRDS("data/_processed/combined_regional.rds")
usfs      <- readRDS("data/_processed/usfs_prepped.rds")
line_pool <- readRDS("data/_processed/line_df.rds")

usfs_eq_regional <- usfs$eq_regional        # has a `region` column

holdout <- combined |> filter(split == "holdout")

## =============================================================================
## Assign every row to a facet region
## =============================================================================
## Three of the four combos already carry a region; only Inventory x Inventory
## is region-free and has to be fanned out.

## regions in which each species actually has a USFS equation
spp_regions <- usfs_eq_regional |> distinct(ScientificName, region)

with_facet <- bind_rows(
  ## USFS data (either equation): the tree's own region
  holdout |>
    filter(data_source == "USFS") |>
    mutate(facet_region = region_data),
  
  ## Inventory data, USFS equation: already one row per regional equation
  holdout |>
    filter(data_source == "Inventory", eq_source == "USFS") |>
    mutate(facet_region = region_eq),
  
  ## Inventory data, Inventory equation: no region at all -> duplicate into
  ## every region where this species has a USFS equation
  holdout |>
    filter(data_source == "Inventory", eq_source == "Inventory") |>
    select(-region, -region_eq) |>
    inner_join(spp_regions, by = "ScientificName",
               relationship = "many-to-many") |>
    mutate(facet_region = region, region_eq = NA_character_)
) |>
  filter(!is.na(facet_region))

## =============================================================================
## Joint status, recomputed WITHIN each facet
## =============================================================================
## The pooled table's `status` is not valid here: a tree can sit inside one
## region's interval and outside another's. Within a single facet, though, each
## tree has exactly one USFS row and one inventory row, so it is well defined
## again -- just region-specific. Recompute rather than reuse.
with_facet <- with_facet |>
  mutate(
    in_USFS      = any(contained[eq_source == "USFS"]),
    in_Inventory = any(contained[eq_source == "Inventory"]),
    .by = c(facet_region, data_source, tree_id)
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

## =============================================================================
## Lines, per species x region
## =============================================================================
## dbh range to draw over: the span of points actually in that facet
facet_range <- with_facet |>
  summarise(dbh_lo = min(dbh_inch), dbh_hi = max(dbh_inch),
            .by = c(ScientificName, facet_region))

## USFS: one regional equation per facet, delta-method PI
usfs_line_reg <- facet_range |>
  reframe(dbh_inch = seq(dbh_lo, dbh_hi, length.out = GRID_N),
          .by = c(ScientificName, facet_region)) |>
  inner_join(usfs_eq_regional,
             by = c("ScientificName", "facet_region" = "region")) |>
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
  select(ScientificName, facet_region, dbh_inch, fit, lo, hi, eq_source)

## Inventory: the pooled line, repeated in every facet. No refit needed -- the
## lmer has no regional term, so the regional line IS the pooled line.
inv_line_reg <- line_pool |>
  filter(eq_source == "Inventory") |>
  inner_join(distinct(facet_range, ScientificName, facet_region),
             by = "ScientificName", relationship = "many-to-many") |>
  select(ScientificName, facet_region, dbh_inch, fit, lo, hi, eq_source)

line_reg <- bind_rows(usfs_line_reg, inv_line_reg)

## =============================================================================
## Regional CV coverage for the annotation labels
## =============================================================================
## Mirrors the METHODOLOGY of 06_coverage_cv.R, but per region -- so these
## numbers will NOT match the (pooled) CV dotplot in value, only in method:
##
##   (A) USFS data      x USFS regional eq   k-fold CV. Folds stratified within
##                                           species x region; the regional
##                                           equation is refit per fold and the
##                                           held-out fold is scored.
##   (B) Inventory data x USFS regional eq   k-fold CV. Regional eqs refit per
##                                           fold; ALL inventory points scored
##                                           against every region's equation,
##                                           containment averaged per point
##                                           across folds first.
##   (C) USFS data      x Inventory eq       SINGLE FIT (fixed lmer, already
##                                           out-of-sample), all USFS points,
##                                           summarised per region.
##   (D) Inventory data x Inventory eq       SINGLE FIT, holdout only. The
##                                           inventory eq and points carry no
##                                           region, so this is ONE number per
##                                           species, repeated in every facet.
##
## As in 06, k comes from the holdout fraction. Same seed constant as 06, but
## folds are drawn within species x region, so fold membership differs too.
SEED_CV <- 22626

common_species <- readRDS("data/_processed/common_species.rds")
PCT_TRAIN      <- usfs$params$pct_train
MIN_N_PER_SPP  <- usfs$params$min_n_per_spp

stopifnot(PCT_TRAIN > 0, PCT_TRAIN < 1)
K <- round(1 / (1 - PCT_TRAIN))
if (K < 2) stop("PCT_TRAIN = ", PCT_TRAIN, " implies k = ", K,
                "; need a holdout fraction of at most 0.5")

message(sprintf("regional_plots_by_spp: regional coverage CV, k = %d folds", K))

## fold on ALL USFS trees, stratified within species x region
set.seed(SEED_CV)
usfs_folded <- usfs$trees |>
  filter(ScientificName %in% common_species) |>
  mutate(fold = sample(rep_len(1:K, n())), .by = c(ScientificName, region)) |>
  select(ScientificName, region, tree_id, dbh_cm, height_m,
         dbh_inch, tree_height_ft, fold)

## inventory target points for (B): one row per tree (region-free)
inv_pts_cv <- combined |>
  filter(data_source == "Inventory", eq_source == "USFS") |>
  distinct(ScientificName, tree_id, dbh_inch, tree_height_ft)

## fit per species x region log-log OLS on whatever is handed in
fit_fold_eq_reg <- function(df) {
  df |>
    nest(.by = c(ScientificName, region)) |>
    mutate(ok = map_lgl(data, ~ nrow(.x) >= MIN_N_PER_SPP)) |>
    filter(ok) |>
    mutate(m = map(data, ~ lm(log(height_m) ~ log(dbh_cm), data = .x))) |>
    transmute(
      ScientificName, region,
      intercept           = map_dbl(m, ~ coef(.x)[1]),
      slope               = map_dbl(m, ~ coef(.x)[2]),
      intercept_se        = map_dbl(m, ~ sqrt(vcov(.x)[1, 1])),
      slope_se            = map_dbl(m, ~ sqrt(vcov(.x)[2, 2])),
      intercept_slope_cov = map_dbl(m, ~ vcov(.x)[1, 2]),
      sigma               = map_dbl(m, ~ summary(.x)$sigma)
    )
}

## score points against a regional coefficient table (delta-method PI, as 06)
score_usfs_eq_reg <- function(pts, eq, by) {
  pts |>
    inner_join(eq, by = by, relationship = "many-to-many") |>
    mutate(
      log_dbh = log(dbh_inch * CM_PER_IN),
      log_fit = intercept + slope * log_dbh,
      log_se  = sqrt(sigma^2 + intercept_se^2 +
                       2 * log_dbh * intercept_slope_cov +
                       log_dbh^2 * slope_se^2),
      lo = exp(log_fit - TCRIT * log_se) * M_TO_FT,
      hi = exp(log_fit + TCRIT * log_se) * M_TO_FT,
      contained = tree_height_ft >= lo & tree_height_ft <= hi
    )
}

## (A) USFS data x USFS regional eq -- genuine k-fold CV
cv_usfs_usfs_reg <- map_dfr(1:K, function(f) {
  eq <- usfs_folded |> filter(fold != f) |> fit_fold_eq_reg()
  usfs_folded |>
    filter(fold == f) |>
    select(ScientificName, region, tree_id, dbh_inch, tree_height_ft) |>
    score_usfs_eq_reg(eq, by = c("ScientificName", "region")) |>
    select(ScientificName, region, tree_id, contained)
}) |>
  summarise(coverage = mean(contained), n = n(),
            .by = c(ScientificName, region)) |>
  mutate(eq_source = "USFS", data_source = "USFS", estimator = "k-fold CV")

## (B) Inventory data x USFS regional eq -- eq refit per fold, all points
## scored against every region's equation; average per point across folds
## first so a point is not weighted by how many folds its equation survived
cv_inv_usfs_reg <- map_dfr(1:K, function(f) {
  eq <- usfs_folded |> filter(fold != f) |> fit_fold_eq_reg()
  inv_pts_cv |>
    score_usfs_eq_reg(eq, by = "ScientificName") |>
    select(ScientificName, region, tree_id, contained)
}) |>
  summarise(contained = mean(contained),
            .by = c(ScientificName, region, tree_id)) |>
  summarise(coverage = mean(contained), n = n(),
            .by = c(ScientificName, region)) |>
  mutate(eq_source = "USFS", data_source = "Inventory",
         estimator = "k-fold CV")

## (C) USFS data x Inventory eq -- fixed model, already out-of-sample;
## all USFS points (not just holdout), per region, mirroring 06
sf_usfs_inv_reg <- combined |>
  filter(data_source == "USFS", eq_source == "Inventory") |>
  summarise(coverage = mean(contained), n = n(),
            .by = c(ScientificName, region_data)) |>
  rename(region = region_data) |>
  mutate(eq_source = "Inventory", data_source = "USFS",
         estimator = "single fit")

## (D) Inventory data x Inventory eq -- holdout only, region-free; fan the
## one species-level number into every region with a USFS equation
sf_inv_inv_reg <- combined |>
  filter(data_source == "Inventory", eq_source == "Inventory",
         split == "holdout") |>
  summarise(coverage = mean(contained), n = n(), .by = ScientificName) |>
  inner_join(spp_regions, by = "ScientificName",
             relationship = "many-to-many") |>
  mutate(eq_source = "Inventory", data_source = "Inventory",
         estimator = "single fit (region-free)")

cov_cv_reg <- bind_rows(cv_usfs_usfs_reg, cv_inv_usfs_reg,
                        sf_usfs_inv_reg, sf_inv_inv_reg)

saveRDS(cov_cv_reg, "data/_processed/coverage_cv_regional.rds")

## label-ready version
cov_cv_reg_pct <- cov_cv_reg |>
  mutate(pct = paste0(round(100 * coverage), "%")) |>
  select(ScientificName, region, eq_source, data_source, pct)

## =============================================================================
## Plot function
## =============================================================================
species_plot_region <- function(focus, dat = with_facet, lines = line_reg,
                                min_usfs = MIN_USFS_HOLDOUT_PER_REGION) {
  
  pts <- dat |> filter(ScientificName %in% focus)
  if (nrow(pts) == 0) stop("no holdout rows for ", focus)
  
  ## drop facets with too few USFS holdout trees to say anything
  keep <- pts |>
    summarise(n_usfs = n_distinct(tree_id[data_source == "USFS"]),
              .by = facet_region) |>
    filter(n_usfs >= min_usfs) |>
    pull(facet_region)
  if (length(keep) == 0) stop("no region clears min_usfs for ", focus)
  
  pts <- pts |> filter(facet_region %in% keep)
  lin <- lines |> filter(ScientificName %in% focus, facet_region %in% keep)
  
  ## coverage labels: regional CV table (see block above) -- CV for the USFS
  ## equation, single fit for the inventory equation, per facet
  cov_lab <- cov_cv_reg_pct |>
    filter(ScientificName %in% focus, region %in% keep) |>
    select(facet_region = region, eq_source, data_source, pct) |>
    complete(facet_region = keep, eq_source = c("USFS", "Inventory"),
             data_source = c("USFS", "Inventory")) |>
    mutate(pct = replace_na(pct, "n/a")) |>
    pivot_wider(names_from = c(eq_source, data_source),
                values_from = pct, names_sep = "_") |>
    mutate(
      lab_usfs = sprintf("USFS eq: \u25b2 %s  \u25cf %s",
                         USFS_USFS, USFS_Inventory),
      lab_inv  = sprintf("Inv. eq: \u25b2 %s  \u25cf %s",
                         Inventory_USFS, Inventory_Inventory)
    )
  
  ## one marker per tree per facet
  pts_1 <- pts |>
    distinct(facet_region, data_source, tree_id,
             dbh_inch, tree_height_ft, status)
  
  n_facets <- length(keep)
  
  ggplot(pts_1, aes(dbh_inch, tree_height_ft)) +
    geom_point(aes(color = status, fill = status, shape = data_source),
               alpha = 1, size = 1.4) +
    scale_fill_manual(values = status_cols, name = "Contained in",
                      guide = guide_legend(override.aes = list(shape = 22, size = 4)),
                      drop = FALSE) +
    scale_color_manual(values = status_cols, name = "Contained in", drop = FALSE) +
    scale_shape_manual(values = src_shapes, name = "Data Source", drop = FALSE) +
    new_scale_color() +
    geom_line(data = lin, aes(y = fit, color = eq_source, linetype = eq_source),
              linewidth = 0.8) +
    geom_line(data = lin, aes(y = lo,  color = eq_source, linetype = eq_source),
              linewidth = 0.8) +
    geom_line(data = lin, aes(y = hi,  color = eq_source, linetype = eq_source),
              linewidth = 0.8) +
    scale_linetype_manual(values = src_lty, name = "Equation Source") +
    scale_color_manual(values = src_colors, name = "Equation Source") +
    ## shared scales across facets: the whole point is comparing regions
    facet_wrap(~ facet_region, nrow = 1) +
    labs(x = "DBH (in)", y = "Tree Height (ft)",
         title = focus,
         subtitle = paste0("Holdout only. USFS equation is region-specific; ",
                           "inventory equation is the same in every panel.\n",
                           "Coverage: USFS eq = regional k-fold CV; Inv. eq = ",
                           "single fit (all USFS points per region; inventory ",
                           "holdout, region-free).")) +
    theme_minimal() +
    theme(legend.position       = "bottom",
          legend.title.position = "top",
          plot.title            = element_text(face = "italic"),
          strip.text            = element_text(face = "bold")) +
    guides(
      linetype = guide_legend(ncol = 1, title.position = "top",
                              keywidth = unit(1.5, "cm")),
      shape    = guide_legend(ncol = 1, title.position = "top"),
      color    = guide_legend(ncol = 2, title.position = "top")
    ) +
    geom_text(data = cov_lab, aes(x = -Inf, y = Inf, label = lab_usfs),
              hjust = -0.05, vjust = 1.6, size = 2.6,
              color = src_colors["USFS"], inherit.aes = FALSE) +
    geom_text(data = cov_lab, aes(x = -Inf, y = Inf, label = lab_inv),
              hjust = -0.05, vjust = 3.0, size = 2.6,
              color = src_colors["Inventory"], inherit.aes = FALSE)
}

## =============================================================================
## Run over every species
## =============================================================================
all_spp <- with_facet |>
  distinct(ScientificName) |>
  semi_join(distinct(line_reg, ScientificName), by = "ScientificName") |>
  pull(ScientificName) |>
  sort()

safe_plot <- possibly(species_plot_region, otherwise = NULL)
plots <- set_names(map(all_spp, safe_plot), all_spp)

failed <- names(plots)[map_lgl(plots, is.null)]
if (length(failed))
  message("skipped ", length(failed), ": ", paste(failed, collapse = ", "))
plots <- compact(plots)

## width scales with facet count so panels stay readable
facet_counts <- with_facet |>
  summarise(k = n_distinct(facet_region), .by = ScientificName) |>
  deframe()

## one PNG per species
iwalk(plots, function(p, sp) {
  k  <- max(1, facet_counts[[sp]])
  fn <- file.path("_figures/by-region/our-refits",
                  paste0(str_replace_all(sp, "[^A-Za-z0-9]+", "_"), ".png"))
  ggsave(fn, p, width = min(4 + 3 * k, 20), height = 5.5, dpi = 300)
})

message("regional_plots_by_spp: wrote ", length(plots), " plots")

########################################
### With original McPherson equations too
########################################

## ---- settings ---------------------------------------------------------------
TS6_PATH  <- "data/_raw/TS6_Growth_coefficients.csv"
CM_PER_IN <- 2.54
M_TO_FT   <- 3.281
TCRIT     <- qt(0.975, df = Inf)
GRID_N    <- 200
MIN_USFS_HOLDOUT_PER_REGION <- 1

## pipeline region label -> TS6 region code. Identity by default; add entries
## like `"Inland Empire" = "InlEmp"` if your labels differ.
REGION_MAP <- identity

dir.create("_figures/by-region/our-refits-and-mcpherson", showWarnings = FALSE, recursive = TRUE)

## ---- theming ----------------------------------------------------------------
src_colors  <- c(Inventory = "#1f78b4", USFS = "#e31a1c", `TS6 published` = "#33a02c")
src_shapes  <- c(Inventory = 16,        USFS = 17)
src_lty     <- c(Inventory = 1,         USFS = 2,         `TS6 published` = 4)
status_cols <- c(both           = "grey",
                 neither        = "black",
                 inventory_only = "#1f78b4",
                 usfs_only      = "#e31a1c")

## =============================================================================
## TS6 coefficients: read + prediction functions
## =============================================================================
ts6_raw <- read_csv(TS6_PATH, show_col_types = FALSE) |>
  rename_with(str_trim)   # "Predicts component " has a trailing space

ts6_ht <- ts6_raw |>
  filter(`Independent variable` == "dbh", `Predicts component` == "tree ht") |>
  transmute(
    region_ts6     = Region,
    ScientificName = `Scientific Name`,
    eqname         = EqName,
    wt             = `Model weight`,
    a, b, c, d, e,
    sigma    = Sigma,
    apps_min = `Apps min`,   # in predicted units (m of height)
    apps_max = `Apps max`
  )

## WLS weight evaluated at x (dbh, cm). Residual variance at x = sigma^2 / w.
ts6_weight <- function(wt, x) {
  case_when(
    wt %in% c("1")                    ~ 1,
    wt %in% c("1/sqrt(dbh)")          ~ 1 / sqrt(x),
    wt %in% c("1/dbh")                ~ 1 / x,
    wt %in% c("1/dbh^2", "1/(dbh)^2") ~ 1 / x^2,
    .default = 1                      # age weights never occur for dbh models
  )
}

## Vectorised fit + approximate 95% PI. x in cm; returns meters.
ts6_predict <- function(eqname, wt, a, b, c, d, e, sigma, x) {
  w    <- ts6_weight(wt, x)
  llx  <- log(log(x + 1))
  poly <- coalesce(a, 0) + coalesce(b, 0) * x + coalesce(c, 0) * x^2 +
    coalesce(d, 0) * x^3 + coalesce(e, 0) * x^4
  is_log <- str_starts(eqname, "loglog") | str_starts(eqname, "expow")
  mse    <- if_else(is_log, c, NA_real_)          # c column stores MSE here
  logfit <- case_when(
    str_starts(eqname, "loglog") ~ a + b * llx + mse / (2 * w),
    str_starts(eqname, "expow")  ~ a + b * x   + mse / (2 * w),
    .default = NA_real_
  )
  log_se <- sqrt(mse / w)
  raw_se <- coalesce(sigma, 0) / sqrt(w)
  tibble(
    fit_m = if_else(is_log, exp(logfit), poly),
    lo_m  = if_else(is_log, exp(logfit - TCRIT * log_se), poly - TCRIT * raw_se),
    hi_m  = if_else(is_log, exp(logfit + TCRIT * log_se), poly + TCRIT * raw_se)
  )
}

## =============================================================================
## Rebuild the faceted holdout table (same logic as 05_plots_by_region.R)
## =============================================================================
combined  <- readRDS("data/_processed/combined_regional.rds")
usfs      <- readRDS("data/_processed/usfs_prepped.rds")
line_pool <- readRDS("data/_processed/line_df.rds")

## Refit-equation coverage labels come from the REGIONAL CV table computed and
## saved by 05_plots_by_region.R: same methodology as the CV dotplot (USFS eq
## refit per fold; inventory eq single fit), but per region -- so the values
## will not match the pooled dotplot. The TS6 label below is unchanged: the
## published equation is a fixed single fit and is not cross-validated.
## NOTE: 05 must be run before this script.
CVREG_PATH <- "data/_processed/coverage_cv_regional.rds"
if (!file.exists(CVREG_PATH))
  stop("missing ", CVREG_PATH,
       " -- run 05_plots_by_region.R first (it computes the regional CV)")
cov_cv_reg_pct <- readRDS(CVREG_PATH) |>
  mutate(pct = paste0(round(100 * coverage), "%")) |>
  select(ScientificName, region, eq_source, data_source, pct)

usfs_eq_regional <- usfs$eq_regional
holdout <- combined |> filter(split == "holdout")

spp_regions <- usfs_eq_regional |> distinct(ScientificName, region)

with_facet <- bind_rows(
  holdout |> filter(data_source == "USFS") |>
    mutate(facet_region = region_data),
  holdout |> filter(data_source == "Inventory", eq_source == "USFS") |>
    mutate(facet_region = region_eq),
  holdout |> filter(data_source == "Inventory", eq_source == "Inventory") |>
    select(-region, -region_eq) |>
    inner_join(spp_regions, by = "ScientificName",
               relationship = "many-to-many") |>
    mutate(facet_region = region, region_eq = NA_character_)
) |>
  filter(!is.na(facet_region))

with_facet <- with_facet |>
  mutate(
    in_USFS      = any(contained[eq_source == "USFS"]),
    in_Inventory = any(contained[eq_source == "Inventory"]),
    .by = c(facet_region, data_source, tree_id)
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

## =============================================================================
## Lines
## =============================================================================
facet_range <- with_facet |>
  summarise(dbh_lo = min(dbh_inch), dbh_hi = max(dbh_inch),
            .by = c(ScientificName, facet_region))

## --- refit USFS regional line (as in 05) -------------------------------------
usfs_line_reg <- facet_range |>
  reframe(dbh_inch = seq(dbh_lo, dbh_hi, length.out = GRID_N),
          .by = c(ScientificName, facet_region)) |>
  inner_join(usfs_eq_regional,
             by = c("ScientificName", "facet_region" = "region")) |>
  mutate(
    log_dbh  = log(dbh_inch * CM_PER_IN),
    log_fit  = intercept + slope * log_dbh,
    line_var = intercept_se^2 + 2 * log_dbh * intercept_slope_cov +
      log_dbh^2 * slope_se^2,
    log_se   = sqrt(sigma^2 + line_var),
    fit = exp(log_fit)                  * M_TO_FT,
    lo  = exp(log_fit - TCRIT * log_se) * M_TO_FT,
    hi  = exp(log_fit + TCRIT * log_se) * M_TO_FT,
    eq_source = "USFS", extrapolated = FALSE
  ) |>
  select(ScientificName, facet_region, dbh_inch, fit, lo, hi,
         eq_source, extrapolated)

## --- pooled inventory line (as in 05) ----------------------------------------
inv_line_reg <- line_pool |>
  filter(eq_source == "Inventory") |>
  inner_join(distinct(facet_range, ScientificName, facet_region),
             by = "ScientificName", relationship = "many-to-many") |>
  mutate(extrapolated = FALSE) |>
  select(ScientificName, facet_region, dbh_inch, fit, lo, hi,
         eq_source, extrapolated)

## --- NEW: published TS6 line --------------------------------------------------
ts6_line_reg <- facet_range |>
  mutate(region_ts6 = REGION_MAP(facet_region)) |>
  inner_join(ts6_ht, by = c("ScientificName", "region_ts6"),
             relationship = "many-to-one") |>
  reframe(
    dbh_inch = seq(dbh_lo, dbh_hi, length.out = GRID_N),
    pick(everything()),
    .by = c(ScientificName, facet_region)
  ) |>
  mutate(
    dbh_cm = dbh_inch * CM_PER_IN,
    ts6_predict(eqname, wt, a, b, c, d, e, sigma, dbh_cm),
    fit = fit_m * M_TO_FT,
    lo  = lo_m  * M_TO_FT,
    hi  = hi_m  * M_TO_FT,
    ## outside the published application range (defined on predicted height)
    extrapolated = fit_m < apps_min | fit_m > apps_max,
    eq_source = "TS6 published"
  ) |>
  select(ScientificName, facet_region, dbh_inch, fit, lo, hi,
         eq_source, extrapolated)

matched <- distinct(ts6_line_reg, ScientificName, facet_region)
message("TS6 dbh->height match: ", nrow(matched), " species x region combos (",
        n_distinct(matched$ScientificName), " species)")
unmatched <- anti_join(distinct(facet_range, ScientificName, facet_region),
                       matched, by = c("ScientificName", "facet_region"))
if (nrow(unmatched))
  message("no TS6 height equation for ", nrow(unmatched),
          " species x region combos (name or region-code mismatch?)")

line_reg <- bind_rows(usfs_line_reg, inv_line_reg, ts6_line_reg)


#### rectangle showing USFS limits:
## published recommended limits for y, and the dbh implied by the fitted model for x

## y-bounds: published application range, in ft
apps_rect <- facet_range |>
  mutate(region_ts6 = REGION_MAP(facet_region)) |>
  inner_join(ts6_ht, by = c("ScientificName", "region_ts6"),
             relationship = "many-to-one") |>
  transmute(ScientificName, facet_region,
            ymin = apps_min * M_TO_FT,
            ymax = apps_max * M_TO_FT)

## x-bounds: dbh span over which the TS6 fit is actually inside that range
apps_x <- ts6_line_reg |>
  filter(!extrapolated) |>
  summarise(xmin = min(dbh_inch), xmax = max(dbh_inch),
            .by = c(ScientificName, facet_region))

apps_rect <- apps_rect |>
  left_join(apps_x, by = c("ScientificName", "facet_region")) |>
  mutate(xmin = coalesce(xmin, -Inf),
         xmax = coalesce(xmax, Inf))

## =============================================================================
## Coverage of the TS6 approximate PI, per facet x data source
## =============================================================================
ts6_cov <- with_facet |>
  distinct(ScientificName, facet_region, data_source, tree_id,
           dbh_inch, tree_height_ft) |>
  mutate(region_ts6 = REGION_MAP(facet_region)) |>
  inner_join(ts6_ht, by = c("ScientificName", "region_ts6"),
             relationship = "many-to-one") |>
  mutate(
    dbh_cm = dbh_inch * CM_PER_IN,
    ts6_predict(eqname, wt, a, b, c, d, e, sigma, dbh_cm),
    contained_ts6 = tree_height_ft >= lo_m * M_TO_FT &
      tree_height_ft <= hi_m * M_TO_FT
  ) |>
  summarise(cov = mean(contained_ts6), n = n(),
            .by = c(ScientificName, facet_region, data_source))

write_csv(ts6_cov, "_figures/species_usfs_eq_coverage.csv")

## =============================================================================
## Plot function (same as pooled_plots_by_spp, plus the TS6 layer + third coverage line)
## =============================================================================
species_plot_ts6 <- function(focus, dat = with_facet, lines = line_reg,
                             cov_ts6 = ts6_cov, apps = apps_rect,
                             min_usfs = MIN_USFS_HOLDOUT_PER_REGION) {
  
  pts <- dat |> filter(ScientificName %in% focus)
  if (nrow(pts) == 0) stop("no holdout rows for ", focus)
  
  keep <- pts |>
    summarise(n_usfs = n_distinct(tree_id[data_source == "USFS"]),
              .by = facet_region) |>
    filter(n_usfs >= min_usfs) |>
    pull(facet_region)
  if (length(keep) == 0) stop("no region clears min_usfs for ", focus)
  
  pts <- pts |> filter(facet_region %in% keep)
  lin <- lines |> filter(ScientificName %in% focus, facet_region %in% keep)
  rect <- apps |> filter(ScientificName %in% focus, facet_region %in% keep)
  
  ## coverage labels: refit equations from the regional CV table (05) +
  ## TS6 (approx PI, per facet, computed here)
  cov_lab <- cov_cv_reg_pct |>
    filter(ScientificName %in% focus, region %in% keep) |>
    select(facet_region = region, eq_source, data_source, pct) |>
    complete(facet_region = keep, eq_source = c("USFS", "Inventory"),
             data_source = c("USFS", "Inventory")) |>
    mutate(pct = replace_na(pct, "n/a")) |>
    pivot_wider(names_from = c(eq_source, data_source),
                values_from = pct, names_sep = "_") |>
    mutate(
      lab_usfs = sprintf("USFS eq: \u25b2 %s  \u25cf %s",
                         USFS_USFS, USFS_Inventory),
      lab_inv  = sprintf("Inv. eq: \u25b2 %s  \u25cf %s",
                         Inventory_USFS, Inventory_Inventory)
    )
  
  cov_lab_ts6 <- cov_ts6 |>
    filter(ScientificName %in% focus, facet_region %in% keep) |>
    mutate(pct = paste0(round(100 * cov), "%")) |>
    select(facet_region, data_source, pct) |>
    complete(facet_region = keep, data_source = c("USFS", "Inventory")) |>
    mutate(pct = replace_na(pct, "n/a")) |>
    pivot_wider(names_from = data_source, values_from = pct) |>
    mutate(lab_ts6 = sprintf("TS6 eq*: \u25b2 %s  \u25cf %s", USFS, Inventory))
  
  pts_1 <- pts |>
    distinct(facet_region, data_source, tree_id,
             dbh_inch, tree_height_ft, status)
  
  ggplot(pts_1, aes(dbh_inch, tree_height_ft)) +
    ## application-range band: fill mapped through alpha so it gets a legend key
    geom_rect(data = rect,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  alpha = "rec"),
              inherit.aes = FALSE, fill = src_colors[["TS6 published"]]) +
    ## separate outline layer so alpha doesn't wash out the border
    geom_rect(data = rect,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA,
              color = src_colors[["TS6 published"]],
              linetype = 3, linewidth = 0.3) +
    geom_point(aes(color = status, fill = status, shape = data_source),
               alpha = 1, size = 1.4,
               show.legend = c(color = FALSE, fill = TRUE, shape = TRUE)) +
    scale_fill_manual(values = status_cols, name = "Contained in",
                      drop = FALSE,
                      guide = guide_legend(order = 4, nrow = 2,
                                           override.aes = list(shape = 22,
                                                               size = 4))) +
    scale_color_manual(values = status_cols, name = "Contained in",
                       drop = FALSE, guide = "none") +
    scale_shape_manual(values = src_shapes, name = "Data Source", drop = FALSE,
                       guide = guide_legend(order = 1, ncol = 1,
                                            title.position = "top")) +
    scale_alpha_manual(
      values = c("rec" = 0.10),
      name = "",
      labels = c("rec" = "USFS usage"),
      guide  = guide_legend(order = 5, title.position = "top",
                            override.aes = list(fill = src_colors[["TS6 published"]],
                                                alpha = 0.35))
    ) +
    new_scale_color() +
    ## fits: solid within TS6 application range, dotted where extrapolated
    geom_line(data = lin |> filter(!extrapolated),
              aes(y = fit, color = eq_source, linetype = eq_source),
              linewidth = 0.8) +
    geom_line(data = lin |> filter(extrapolated),
              aes(y = fit, color = eq_source, group = eq_source),
              linetype = 3, linewidth = 0.6, show.legend = FALSE) +
    ## PI bounds (TS6 bounds are residual-only)
    geom_line(data = lin, aes(y = lo, color = eq_source, linetype = eq_source),
              linewidth = 0.5, alpha = 0.8, show.legend = FALSE) +
    geom_line(data = lin, aes(y = hi, color = eq_source, linetype = eq_source),
              linewidth = 0.5, alpha = 0.8, show.legend = FALSE) +
    scale_linetype_manual(values = src_lty, name = "Equation Source",
                          guide = guide_legend(order = 2, ncol = 1,
                                               title.position = "top",
                                               keywidth = unit(1.1, "cm"))) +
    scale_color_manual(values = src_colors, name = "Equation Source",
                       guide = guide_legend(order = 2, ncol = 1,
                                            title.position = "top",
                                            keywidth = unit(1.1, "cm"))) +
    facet_wrap(~ facet_region, nrow = 1) +
    labs(x = "DBH (in)", y = "Tree Height (ft)",
         title = focus) +
    theme_minimal() +
    theme(legend.position       = "bottom",
          legend.box.just       = "left",
          legend.title.position = "top",
          plot.title            = element_text(face = "italic"),
          strip.text            = element_text(face = "bold")) +
## coverage stats
    geom_text(data = cov_lab, aes(x = -Inf, y = Inf, label = lab_usfs),
              hjust = -0.05, vjust = 1.6, size = 3.5,
              color = src_colors[["USFS"]], inherit.aes = FALSE) +
    geom_text(data = cov_lab, aes(x = -Inf, y = Inf, label = lab_inv),
              hjust = -0.05, vjust = 3.0, size = 3.5,
              color = src_colors[["Inventory"]], inherit.aes = FALSE) +
    geom_text(data = cov_lab_ts6, aes(x = -Inf, y = Inf, label = lab_ts6),
              hjust = -0.05, vjust = 4.4, size = 3.5,
              color = src_colors[["TS6 published"]], inherit.aes = FALSE)
}

## =============================================================================
## Run over every species that has at least one TS6 height equation
## =============================================================================
all_spp <- with_facet |>
  distinct(ScientificName) |>
  semi_join(matched, by = "ScientificName") |>
  pull(ScientificName) |>
  sort()

safe_plot <- possibly(species_plot_ts6, otherwise = NULL)
plots <- set_names(map(all_spp, safe_plot), all_spp)

failed <- names(plots)[map_lgl(plots, is.null)]
if (length(failed))
  message("skipped ", length(failed), ": ", paste(failed, collapse = ", "))
plots <- compact(plots)

facet_counts <- with_facet |>
  summarise(k = n_distinct(facet_region), .by = ScientificName) |>
  deframe()


iwalk(plots, function(p, sp) {
  k  <- max(1, facet_counts[[sp]])
  fn <- file.path("_figures/by-region/our-refits-and-mcpherson",
                  paste0(str_replace_all(sp, "[^A-Za-z0-9]+", "_"), ".png"))
  ggsave(fn, p, width = min(4 + 3 * k, 20), height = 5.8, dpi = 300)
})

message("regional_plots_by_spp: wrote ", length(plots), " plots")

