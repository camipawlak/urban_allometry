## =============================================================================
## pooled_plots.R
##
## Coverage dotplots: cross-validated and holdout only. k is derived from the holdout fraction
## used in 01_prep_usfs.R:  k = round(1 / (1 - PCT_TRAIN)).
## PCT_TRAIN = 0.9 -> k = 10;  0.8 -> 5;  2/3 -> 3.
##
## Attaches its own packages so it can be run standalone. ggnewscale supplies
## new_scale_color(), used by species_plot(); without it every per-species plot
## fails, and because those calls are wrapped in possibly() the failure surfaces
## only as "skipped N" rather than an error.
##
## IMPORTANT -- only two of the four combos are genuinely cross-validated:
##
##   (A) USFS data     x USFS eq        CV. Refit per fold, score the held-out
##                                      fold. This is the real thing.
##   (B) Inventory data x USFS eq       CV. Refit per fold, score ALL inventory
##                                      points each time, then average. The
##                                      inventory points were never in any USFS
##                                      fit, so this is out-of-sample either
##                                      way; the folds measure how much the
##                                      USFS equation itself wobbles.
##   (C) USFS data     x Inventory eq   SINGLE FIT. The lmer is fixed and USFS
##                                      data never entered it, so every point is
##                                      already out-of-sample. Nothing to CV.
##   (D) Inventory data x Inventory eq  SINGLE FIT, holdout partition only.
##                                      Cross-validating this would mean
##                                      refitting the lmer k times on ~1.7M
##                                      rows, which this project does not do.
##
## So (A) and (B) are k-fold estimates and (C) and (D) are single-split
## estimates sitting in the same panel. They are comparable in intent but not
## in variance -- the CV points are more stable. The plot says so.
##
## FOR CV Coverage
##
## INPUT   data/_processed/usfs_prepped.rds
##         data/_processed/combined_pooled.rds
##         data/_processed/common_species.rds
## OUTPUT  _figures/coverage-cv-by-species.png
##         data/_processed/coverage_cv.rds
##
## FOR holdout only:
## Holdout-only _figures, built from the POOLED combined table:
##   coverage-by-species.png   prediction-interval coverage by species,
##                                faceted by equation source
##
## Run with working directory = usfs_comparison/
##
## INPUT   data/_processed/combined_pooled.rds
##         data/_processed/line_df.rds
## OUTPUT  _figures/*.png
##        _figures/coverage-by-species.png
## =============================================================================
## ---- settings ---------------------------------------------------------------
library(tidyverse)
library(ggnewscale)
source("scripts/r/fig-theme.R")

CM_PER_IN <- 2.54
M_TO_FT   <- 3.281
TCRIT     <- qt(0.975, df = Inf)
SEED_CV   <- 22626

dir.create("_figures", showWarnings = FALSE)

## ---- theming ----------------------------------------------------------------
src_colors <- c(Inventory = "#1f78b4", USFS = "#e31a1c")
src_shapes <- c(Inventory = 16,        USFS = 17)

## ---- data -------------------------------------------------------------------
usfs           <- readRDS("data/_processed/usfs_prepped.rds")
combined       <- readRDS("data/_processed/combined_pooled.rds")
common_species <- readRDS("data/_processed/common_species.rds")

PCT_TRAIN     <- usfs$params$pct_train
MIN_N_PER_SPP <- usfs$params$min_n_per_spp

## k from the holdout fraction
stopifnot(PCT_TRAIN > 0, PCT_TRAIN < 1)
K <- round(1 / (1 - PCT_TRAIN))
if (K < 2) stop("PCT_TRAIN = ", PCT_TRAIN, " implies k = ", K,
                "; need a holdout fraction of at most 0.5")

message(sprintf("06_coverage_cv: PCT_TRAIN = %.3f -> k = %d folds", PCT_TRAIN, K))

## ---- fold assignment --------------------------------------------------------
## Fold on ALL USFS trees, train and holdout alike -- CV replaces the single
## split rather than nesting inside it. Stratified within species so each fold
## has a comparable number of trees per species.
set.seed(20000)
usfs_folded <- usfs$trees |>
  filter(ScientificName %in% common_species) |>
  mutate(fold = sample(rep_len(1:K, n())), .by = ScientificName) |>
  select(ScientificName, spp, tree_id, dbh_cm, height_m,
         dbh_inch, tree_height_ft, fold)

## inventory target points for (B): one row per tree
inv_pts <- combined |>
  filter(combo == "Inventory_USFS") |>
  distinct(ScientificName, tree_id, dbh_inch, tree_height_ft)

## ---- helpers ----------------------------------------------------------------
## fit per-species log-log OLS on whatever is handed in
fit_fold_eq <- function(df) {
  df |>
    nest(.by = ScientificName) |>
    mutate(ok = map_lgl(data, ~ nrow(.x) >= MIN_N_PER_SPP)) |>
    filter(ok) |>
    mutate(m = map(data, ~ lm(log(height_m) ~ log(dbh_cm), data = .x))) |>
    transmute(
      ScientificName,
      intercept           = map_dbl(m, ~ coef(.x)[1]),
      slope               = map_dbl(m, ~ coef(.x)[2]),
      intercept_se        = map_dbl(m, ~ sqrt(vcov(.x)[1, 1])),
      slope_se            = map_dbl(m, ~ sqrt(vcov(.x)[2, 2])),
      intercept_slope_cov = map_dbl(m, ~ vcov(.x)[1, 2]),
      sigma               = map_dbl(m, ~ summary(.x)$sigma)
    )
}

## score points against a coefficient table (delta-method PI)
score_usfs_eq <- function(pts, eq) {
  pts |>
    inner_join(eq, by = "ScientificName") |>
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

## =============================================================================
## (A) USFS data x USFS equation -- genuine k-fold CV
## =============================================================================
cv_usfs_usfs <- map_dfr(1:K, function(f) {
  eq <- usfs_folded |> filter(fold != f) |> fit_fold_eq()
  usfs_folded |>
    filter(fold == f) |>
    select(ScientificName, tree_id, dbh_inch, tree_height_ft) |>
    score_usfs_eq(eq) |>
    select(ScientificName, tree_id, contained)
}) |>
  summarise(coverage = mean(contained), n = n(), .by = ScientificName) |>
  mutate(eq_source = "USFS", data_source = "USFS", estimator = "k-fold CV")

## =============================================================================
## (B) Inventory data x USFS equation -- eq refit per fold, all points scored
## =============================================================================
## Average containment per POINT across folds first, then average over points,
## so a species is not weighted by how many folds it survived.
cv_inv_usfs <- map_dfr(1:K, function(f) {
  eq <- usfs_folded |> filter(fold != f) |> fit_fold_eq()
  inv_pts |>
    score_usfs_eq(eq) |>
    select(ScientificName, tree_id, contained)
}) |>
  summarise(contained = mean(contained), .by = c(ScientificName, tree_id)) |>
  summarise(coverage = mean(contained), n = n(), .by = ScientificName) |>
  mutate(eq_source = "USFS", data_source = "Inventory", estimator = "k-fold CV")

## =============================================================================
## (C) USFS data x Inventory equation -- fixed model, already out-of-sample
## =============================================================================
sf_usfs_inv <- combined |>
  filter(combo == "USFS_Inventory") |>
  summarise(coverage = mean(contained), n = n(), .by = ScientificName) |>
  mutate(eq_source = "Inventory", data_source = "USFS",
         estimator = "single fit")

## =============================================================================
## (D) Inventory data x Inventory equation -- holdout partition only
## =============================================================================
## Restricted to split == "holdout": the training rows are in-sample for the
## lmer and would inflate coverage.
sf_inv_inv <- combined |>
  filter(combo == "Inventory_Inventory", split == "holdout") |>
  summarise(coverage = mean(contained), n = n(), .by = ScientificName) |>
  mutate(eq_source = "Inventory", data_source = "Inventory",
         estimator = "single fit")

## ---- assemble ---------------------------------------------------------------
cov_cv <- bind_rows(cv_usfs_usfs, cv_inv_usfs, sf_usfs_inv, sf_inv_inv) |>
  mutate(
    ScientificName = fct_reorder(ScientificName, coverage, .fun = mean),
    eq_source      = factor(eq_source, levels = c("USFS", "Inventory")),
    data_source    = factor(data_source, levels = c("Inventory", "USFS")),
    facet_lab      = paste0("Model: ", eq_source)
  )

saveRDS(cov_cv, "data/_processed/coverage_cv.rds")

## ---- plot -------------------------------------------------------------------
p_cv <- ggplot(cov_cv, aes(x = coverage, y = ScientificName)) +
  geom_vline(xintercept = 0.95, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_point(aes(color = data_source, shape = data_source), size = 3, alpha = 0.9) +
  scale_color_manual(values = src_colors, name = "Data Source") +
  scale_shape_manual(values = src_shapes, name = "Data Source") +
  scale_x_continuous(labels = scales::label_percent(scale = 1, accuracy = 1)) +
  facet_wrap(~ facet_lab) +
  labs(
    x = "Coverage (points inside 95% PI)", y = NULL,
    title = "Prediction-interval coverage by species",
  ) +
  theme_minimal() +
  theme(legend.position       = "top",
        legend.title.position = "left",
        plot.title = element_text(size = 17),
        panel.grid.major.y    = element_line(color = "grey92"),
        strip.text.x.top = element_text(size = 14),
        axis.text.y           = element_text(size = 11, face = "italic"),
        axis.text.x = element_text(size = 13),
        panel.spacing.x = unit(1.2, "lines")) +
  coord_cartesian(xlim = c(0.4, 1))

ggsave("_figures/coverage-cv-by-species.png", p_cv,
       width = 8, height = 9, dpi = 300)

## ---- console summary --------------------------------------------------------
cov_cv |>
  summarise(species = n(),
            median_coverage = median(coverage),
            .by = c(eq_source, data_source, estimator)) |>
  arrange(eq_source, data_source) |>
  print()

message("wrote _figures/coverage-cv-by-species.png")

#####
## Plot versions without cross validation
###

## ---- settings ---------------------------------------------------------------
MIN_N_HOLDOUT     <- 0   # drop species with fewer than this many holdout trees

dir.create("_figures", showWarnings = FALSE)

## ---- theming (was R/plot_theming.R) -----------------------------------------
src_colors  <- c(Inventory = "#1f78b4", USFS = "#e31a1c")
src_shapes  <- c(Inventory = 16,        USFS = 17)
src_lty     <- c(Inventory = 1,         USFS = 2)
status_cols <- c(both           = "grey",
                 neither        = "black",
                 inventory_only = "#1f78b4",
                 usfs_only      = "#e31a1c")

## ---- data -------------------------------------------------------------------
combined <- readRDS("data/_processed/combined_pooled.rds")
line_df  <- readRDS("data/_processed/line_df.rds")

## Coverage numbers for the per-species annotations come from the CV dotplot
## so the labels on these plots match that figure
## exactly. NOTE: 06 must therefore be run before this script.
cov_cv <- readRDS("data/_processed/coverage_cv.rds") |>
  mutate(across(c(ScientificName, eq_source, data_source), as.character),
         pct = round(100 * coverage)) |>
  select(ScientificName, eq_source, data_source, pct)

## HOLDOUT ONLY: trees that were held out of the fit for their own source.
holdout <- combined |> filter(split == "holdout")

## =============================================================================
## 1. Prediction-interval coverage by species
## =============================================================================
cov_df <- holdout |>
  summarise(coverage = mean(contained),
            n        = n(),
            .by = c(ScientificName, eq_source, data_source)) |>
  filter(n >= MIN_N_HOLDOUT) |>
  mutate(
    ScientificName = fct_reorder(ScientificName, coverage, .fun = mean),
    eq_source      = factor(eq_source, levels = c("USFS", "Inventory")),
    facet_lab      = paste0("Model: ", eq_source)
  )

p_coverage <- ggplot(cov_df, aes(x = coverage, y = ScientificName)) +
  geom_vline(xintercept = 0.95, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_point(aes(color = data_source, shape = data_source), size = 3, alpha = 0.9) +
  scale_color_manual(values = src_colors, name = "Data Source") +
  scale_shape_manual(values = src_shapes, name = "Data Source") +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  facet_wrap(~ facet_lab) +
  labs(x = "Coverage (points inside 95% PI)", y = NULL,
       title = "Prediction-interval coverage by species", 
       subtitle = "No cross-validation, so lots of 100% coverage for small holdout samples.") +
  theme_minimal() +
  theme(legend.position       = "top",
        legend.title.position = "left",
        plot.title = element_text(size = 17),
        panel.grid.major.y    = element_line(color = "grey92"),
        strip.text.x.top = element_text(size = 14),
        axis.text.y           = element_text(size = 11, face = "italic"),
        axis.text.x = element_text(size = 13),
        panel.spacing.x = unit(1.2, "lines")) + 
  coord_cartesian(xlim = c(0.4, 1))

ggsave("_figures/coverage-by-species.png", p_coverage,
       width = 8, height = 9, dpi = 300)

## =============================================================================
## 2-3. Per-species scatter + fitted line and PI
## =============================================================================
species_plot <- function(focus, dat = holdout, lines = line_df) {
  
  pts <- dat |> filter(ScientificName %in% focus)
  if (nrow(pts) == 0) stop("no holdout rows for ", focus)
  
  ## coverage annotation: one number per (equation source x data source),
  ## taken from coverage_cv.rds so it matches the CV dotplot (06). USFS-eq
  ## values are k-fold CV; inventory-eq values are single-fit (see 06 header).
  cov_lab <- cov_cv |>
    filter(ScientificName %in% focus) |>
    complete(ScientificName = focus,
             eq_source   = c("USFS", "Inventory"),
             data_source = c("USFS", "Inventory")) |>
    pivot_wider(names_from = c(eq_source, data_source),
                values_from = pct, names_sep = "_") |>
    mutate(
      hdr_usfs  = "USFS Eq. Coverage (CV):",
      body_usfs = sprintf("\u2009\u25b2 USFS DATA = %s\n  \u25cf \u2009  INV. DATA = %s",
                          ifelse(is.na(USFS_USFS),      "n/a", paste0(USFS_USFS, "%")),
                          ifelse(is.na(USFS_Inventory), "n/a", paste0(USFS_Inventory, "%"))),
      hdr_inv   = "Inv. Eq. Coverage:",
      body_inv  = sprintf("\u2009\u25b2 USFS DATA = %s\n  \u25cf \u2009 INV. DATA = %s",
                          ifelse(is.na(Inventory_USFS),      "n/a", paste0(Inventory_USFS, "%")),
                          ifelse(is.na(Inventory_Inventory), "n/a", paste0(Inventory_Inventory, "%")))
    )
  
  ## one marker per tree (the table holds one row per equation source, but the
  ## point itself and its joint status are identical across them)
  pts_1 <- pts |>
    distinct(ScientificName, data_source, tree_id,
             dbh_inch, tree_height_ft, status)
  
  lin <- lines |> filter(ScientificName %in% focus)
  
  ggplot(pts_1, aes(dbh_inch, tree_height_ft)) +
    geom_point(aes(color = status, fill = status, shape = data_source),
               alpha = 1, size = 1.5) +
    scale_fill_manual(values = status_cols, name = "Contained in",
                      guide = guide_legend(override.aes = list(shape = 22, size = 4)),
                      drop = FALSE) +
    scale_color_manual(values = status_cols, name = "Contained in", drop = FALSE) +
    scale_shape_manual(values = src_shapes, name = "Data Source", drop = FALSE) +
    new_scale_color() +
    geom_line(data = lin, aes(y = fit, color = eq_source, linetype = eq_source),
              linewidth = 0.9) +
    geom_line(data = lin, aes(y = lo,  color = eq_source, linetype = eq_source),
              linewidth = 0.9) +
    geom_line(data = lin, aes(y = hi,  color = eq_source, linetype = eq_source),
              linewidth = 0.9) +
    scale_linetype_manual(values = src_lty, name = "Equation Source") +
    scale_color_manual(values = src_colors, name = "Equation Source") +
    facet_wrap(~ ScientificName, scales = "free") +
    labs(x = "DBH (in)", y = "Tree Height (ft)") +
    theme_minimal() +
    theme(legend.position = "bottom", 
          legend.title.position = "top",
          strip.text = element_text(face = "italic", size = 14)) +
    guides(
      linetype = guide_legend(ncol = 1, title.position = "top",
                              keywidth = unit(1.5, "cm")),
      shape    = guide_legend(ncol = 1, title.position = "top"),
      color    = guide_legend(ncol = 2, title.position = "top")
    ) +
    ## USFS coverage block
    geom_text(data = cov_lab, aes(x = -Inf, y = Inf, label = hdr_usfs),
              hjust = -.25, vjust = 1.4, size = 3, color = src_colors["USFS"],
              fontface = "bold", inherit.aes = FALSE) +
    geom_text(data = cov_lab, aes(x = -Inf, y = Inf, label = body_usfs),
              hjust = -.35, vjust = 1.8, size = 3, color = "black",
              lineheight = 0.95, inherit.aes = FALSE) +
    ## Inventory coverage block, pushed below the USFS one
    geom_text(data = cov_lab, aes(x = -Inf, y = Inf, label = hdr_inv),
              hjust = -.25, vjust = 6.6, size = 3, color = src_colors["Inventory"],
              fontface = "bold", inherit.aes = FALSE) +
    geom_text(data = cov_lab, aes(x = -Inf, y = Inf, label = body_inv),
              hjust = -.35, vjust = 4, size = 3, color = "black",
              lineheight = 0.95, inherit.aes = FALSE)
}


## ---- all species ------------------------------------------------------------
## Only species that have both holdout points and fitted lines will plot.
all_spp <- intersect(unique(holdout$ScientificName),
                     unique(line_df$ScientificName)) |> sort()

## safe() so one bad species doesn't kill the whole run
safe_plot <- possibly(species_plot, otherwise = NULL)
plots <- set_names(map(all_spp, safe_plot), all_spp)

failed <- names(plots)[map_lgl(plots, is.null)]
if (length(failed))
  message("skipped ", length(failed), ": ", paste(failed, collapse = ", "))
plots <- compact(plots)


## (b) one PNG per species
dir.create("_figures/species", showWarnings = FALSE)
iwalk(plots, function(p, sp) {
  fn <- file.path("_figures/species",
                  paste0(str_replace_all(sp, "[^A-Za-z0-9]+", "_"), ".png"))
  save_fig(p, fn, height = 4.4)
})

message("all-species: wrote ", length(plots), " plots")


message("pooled_plots: wrote _figures/coverage-by-species.png")
