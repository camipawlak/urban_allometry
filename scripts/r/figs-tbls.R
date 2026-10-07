library(tidyverse)
library(lme4)
library(writexl)
library(patchwork)
library(cowplot)
source("scripts/r/profile-signature.R")
source("scripts/r/fig-theme.R")

rslt_dir <- "rslt/"

# helper: extract a tidy parameter table from a fitted lmer model
# profile_ci: matrix from confint(profile(..., which = "theta_")), one row per
#             VarCorr entry plus the residual SD; NULL when no current profile exists
#
# Variance and correlation rows carry the profile-likelihood bounds directly. Those
# bounds are asymmetric about the estimate and are the interval to report; `se` on
# those rows is only the symmetric width the interval implies, kept for continuity.
# Fixed effects have no profile here (profile() was run with which = "theta_"), so
# they carry Wald bounds. `ci_method` records which applies to each row.
extract_params <- function(fit, model_label, profile_ci = NULL) {

  # fixed effects with Wald SEs and Wald intervals
  fe <- summary(fit)$coefficients
  fixed_tbl <- tibble(
    param_type = "fixed",
    group      = NA_character_,
    term       = rownames(fe),
    estimate   = fe[, "Estimate"],
    se         = fe[, "Std. Error"],
    lower      = fe[, "Estimate"] - 1.96 * fe[, "Std. Error"],
    upper      = fe[, "Estimate"] + 1.96 * fe[, "Std. Error"],
    ci_method  = "Wald"
  )

  # variance and correlation parameters from VarCorr
  vc_df <- as.data.frame(VarCorr(fit))

  # Map profile CI rows onto vc_df rows.
  #
  # lme4 profiles on the sd/correlation scale by default (prof.scale = "sdcor"),
  # emitting one row per VarCorr entry, in VarCorr's own order, with the residual
  # SD last. So the mapping is positional against vc_df.
  #
  # Do NOT reconstruct Cholesky positions here. Some lme4 versions emit the
  # variance parameters in Cholesky order (SD, correlation, SD) instead of sdcor
  # order (SD, SD, correlation); assuming the former silently pairs correlation
  # bounds with an SD row, producing negative bounds on a standard deviation and
  # intervals that do not contain their estimate. The check at the end of this
  # block exists to catch exactly that, whatever the ordering turns out to be.
  se_vc <- rep(NA_real_, nrow(vc_df))
  lo_vc <- rep(NA_real_, nrow(vc_df))
  hi_vc <- rep(NA_real_, nrow(vc_df))
  ci_vc <- rep(NA_character_, nrow(vc_df))

  if (!is.null(profile_ci)) {
    pnames  <- rownames(profile_ci)
    # ".sigma" with signames = TRUE (the default), "sigma" with signames = FALSE
    sig_idx <- grep("^[.]?sigma$", pnames)
    re_idx  <- setdiff(seq_along(pnames), sig_idx)

    vc_re  <- which(vc_df$grp != "Residual")
    vc_res <- which(vc_df$grp == "Residual")

    if (length(re_idx) != length(vc_re))
      stop(sprintf(paste("profile has %d random-effect rows but the model has %d",
                         "variance parameters -- cannot map them safely"),
                   length(re_idx), length(vc_re)))

    lo_vc[vc_re] <- profile_ci[re_idx, "2.5 %"]
    hi_vc[vc_re] <- profile_ci[re_idx, "97.5 %"]
    se_vc[vc_re] <- (hi_vc[vc_re] - lo_vc[vc_re]) / (2 * 1.96)
    ci_vc[vc_re] <- "profile"

    # residual SD: no SE reported (a symmetric SE misrepresents a scale parameter),
    # but the profile interval itself is carried through.
    if (length(sig_idx) == 1L && length(vc_res) == 1L) {
      lo_vc[vc_res] <- profile_ci[sig_idx, "2.5 %"]
      hi_vc[vc_res] <- profile_ci[sig_idx, "97.5 %"]
      ci_vc[vc_res] <- "profile"
    }

    # every interval must bracket its own estimate; if not, the rows are misaligned
    bad <- which(!is.na(lo_vc) &
                 (vc_df$sdcor < lo_vc - 1e-8 | vc_df$sdcor > hi_vc + 1e-8))
    if (length(bad) > 0)
      stop(sprintf(paste("profile bounds do not bracket the estimate for %s",
                         "(est %.4f not in [%.4f, %.4f]) -- profile rows are",
                         "misaligned with VarCorr"),
                   model_label, vc_df$sdcor[bad[1]], lo_vc[bad[1]], hi_vc[bad[1]]))
  }

  var_idx <- which(is.na(vc_df$var2))
  var_tbl <- tibble(
    param_type = "variance (SD)",
    group      = vc_df$grp[var_idx],
    term       = coalesce(vc_df$var1[var_idx], ""),
    estimate   = vc_df$sdcor[var_idx],
    se         = se_vc[var_idx],
    lower      = lo_vc[var_idx],
    upper      = hi_vc[var_idx],
    ci_method  = ci_vc[var_idx]
  )

  cor_idx <- which(!is.na(vc_df$var2))
  cor_tbl <- tibble(
    param_type = "correlation",
    group      = vc_df$grp[cor_idx],
    term       = paste0(vc_df$var1[cor_idx], " × ", vc_df$var2[cor_idx]),
    estimate   = vc_df$sdcor[cor_idx],
    se         = se_vc[cor_idx],
    lower      = lo_vc[cor_idx],
    upper      = hi_vc[cor_idx],
    ci_method  = ci_vc[cor_idx]
  )

  bind_rows(fixed_tbl, var_tbl, cor_tbl) |>
    mutate(model = model_label, .before = everything())
}

## load model objects -----------------------------------------------------------

# A stored profile is only valid for the fit it was computed from, so check the
# signature it carries rather than trusting timestamps (see profile-signature.R).
# Returning NULL makes extract_params() leave the variance-parameter bounds and
# ci_method as NA, rather than silently pairing new estimates with old intervals.
read_profile_if_current <- function(label, fit, ci_path) {
  if (!file.exists(ci_path)) {
    warning(sprintf("%s: no profile found -- variance CIs will be empty.", label),
            call. = FALSE)
    return(NULL)
  }
  ci <- read_rds(ci_path)
  if (!profile_matches_fit(ci, fit)) {
    warning(sprintf(paste("%s: stored profile was computed from a different fit;",
                          "ignoring it -- variance CIs will be empty. Re-run",
                          "scripts/r/profile.R."), label),
            call. = FALSE)
    return(NULL)
  }
  ci
}

message("Loading dbh-ht model...")
fit_dbh_ht     <- read_rds("rslt/dbh-ht/_dbh-ht-model-obj.rds")
ci_dbh_ht      <- read_profile_if_current("dbh-ht", fit_dbh_ht,
                                          "rslt/dbh-ht/ci_dbh_ht_profile.rds")

message("Loading age-dbh model...")
fit_age_dbh    <- read_rds("rslt/age-dbh/_dbh-age-model-obj.rds")
ci_age_dbh     <- read_profile_if_current("age-dbh", fit_age_dbh,
                                          "rslt/age-dbh/ci_age_dbh_profile.rds")

message("Loading dbh-canopy model...")
fit_dbh_canopy <- read_rds("rslt/dbh-canopy/_dbh-canopy-model-obj.rds")
ci_dbh_canopy  <- read_profile_if_current("dbh-canopy", fit_dbh_canopy,
                                          "rslt/dbh-canopy/ci_dbh_canopy_profile.rds")

message("Loading age-ht model...")
fit_age_ht     <- read_rds("rslt/age-ht/_age-ht-model-obj.rds")
ci_age_ht      <- read_profile_if_current("age-ht", fit_age_ht,
                                          "rslt/age-ht/ci_age_ht_profile.rds")

## extract parameter tables -----------------------------------------------------

params_dbh_ht     <- extract_params(fit_dbh_ht,     "dbh-ht",     ci_dbh_ht)
params_age_dbh    <- extract_params(fit_age_dbh,    "age-dbh",    ci_age_dbh)
params_dbh_canopy <- extract_params(fit_dbh_canopy, "dbh-canopy", ci_dbh_canopy)
params_age_ht     <- extract_params(fit_age_ht,     "age-ht",     ci_age_ht)

params_all <- bind_rows(params_dbh_ht, params_age_dbh, params_dbh_canopy, params_age_ht)

## print ------------------------------------------------------------------------

params_all

## export -----------------------------------------------------------------------

writexl::write_xlsx(
  list(
    all        = params_all,
    dbh_ht     = params_dbh_ht,
    age_dbh    = params_age_dbh,
    dbh_canopy = params_dbh_canopy,
    age_ht     = params_age_ht
  ),
  path = file.path(rslt_dir, "param-tables.xlsx")
)

write_csv(params_all, file.path(rslt_dir, "param-tables.csv"))

message("Exported to rslt/param-tables.xlsx and rslt/param-tables.csv")

## allometric curves figure -----------------------------------------------------

# load partitions and coefficient tables
load("data/_processed/wca-partitions.RData")

coefs_dbh_ht     <- readRDS("rslt/dbh-ht/dbh-ht-coefs.rds")
coefs_age_dbh    <- readRDS("rslt/age-dbh/dbh-age-coefs.rds")
coefs_dbh_canopy <- readRDS("rslt/dbh-canopy/dbh-canopy-coefs.rds")
coefs_age_ht     <- readRDS("rslt/age-ht/age-ht-coefs.rds")

# model-specific training subsets
dbh_ht_data <- wca_train |> filter(dbh_height_is_anom == FALSE)

age_dbh_data <- wca_train |>
  filter(age_dbh_height_is_anom == FALSE, !is.na(age)) |>
  add_count(spp, name = "n_spp") |>
  filter(n_spp >= 30) |>
  select(-n_spp)

dbh_canopy_data <- wca_train |>
  filter(dbh_height_canopy_width_is_anom == FALSE) |>
  drop_na(dbh, canopy_width_avg_m) |>
  add_count(spp, name = "n_spp") |>
  filter(n_spp >= 30) |>
  select(-n_spp) |>
  rename(canopy = canopy_width_avg_m)

age_ht_data <- wca_train |>
  filter(age_dbh_height_is_anom == FALSE, !is.na(age), !is.na(height)) |>
  add_count(spp, name = "n_spp") |>
  filter(n_spp >= 30) |>
  select(-n_spp)

# generate back-transformed curves from coefficient tables
make_curves <- function(spp_coefs, pop_coefs, x_range, n = 200) {
  x_seq <- exp(seq(log(x_range[1]), log(x_range[2]), length.out = n))

  spp_curves <- spp_coefs |>
    select(spp, intercept, slope) |>
    rowwise() |>
    reframe(x = x_seq, y = exp(intercept) * x^slope, spp = spp)

  pop_int   <- pop_coefs$estimate[pop_coefs$term == "intercept"]
  pop_slope <- pop_coefs$estimate[pop_coefs$term == "slope"]
  pop_curve <- tibble(x = x_seq, y = exp(pop_int) * x_seq^pop_slope)

  list(spp = spp_curves, pop = pop_curve)
}

curves_dbh_ht     <- make_curves(coefs_dbh_ht$species,     coefs_dbh_ht$population,     range(dbh_ht_data$dbh))
curves_age_dbh    <- make_curves(coefs_age_dbh$species,    coefs_age_dbh$population,    range(age_dbh_data$age))
curves_dbh_canopy <- make_curves(coefs_dbh_canopy$species, coefs_dbh_canopy$population, range(dbh_canopy_data$dbh))
curves_age_ht     <- make_curves(coefs_age_ht$species,     coefs_age_ht$population,     range(age_ht_data$age))

# panel builder
make_panel <- function(curves, x_lab, y_lab) {
  ggplot() +
    geom_line(
      data = curves$spp,
      aes(x = x, y = y, group = spp, color = "Species"),
      linewidth = 0.35, alpha = 0.6
    ) +
    geom_line(
      data = curves$pop,
      aes(x = x, y = y, color = "Population"),
      linewidth = 1.2
    ) +
    scale_color_manual(
      name = NULL,
      values = c("Species"    = unname(fig_colors["species"]),
                 "Population" = unname(fig_colors["population"]))
    ) +
    labs(x = x_lab, y = y_lab) +
    theme_allometry()
}

fig_allometric_curves <-
  (make_panel(curves_dbh_ht,     "DBH (cm)",    "Height (m)") |
   make_panel(curves_dbh_canopy, "DBH (cm)",    "Canopy width (m)") |
   make_panel(curves_age_dbh,    "Age (years)", "DBH (cm)") |
   make_panel(curves_age_ht,     "Age (years)", "Height (m)")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# four panels across the full text width
save_fig(fig_allometric_curves, "_img/fig-allometric-curves.png", height = 2.5)

message("Saved _img/fig-allometric-curves.png")
