library(lme4)
source("scripts/r/profile-signature.R")

# Compute likelihood profile CIs for the variance-covariance parameters (theta_)
# of each fitted lmer model and save them alongside the model in rslt/.
#
# Profiling is expensive: profile() refits the model repeatedly while stepping each
# variance parameter across its likelihood surface, so cost scales with the fit. The
# dbh-ht model (~720k rows, nested genus/spp effects) took over four hours at lme4's
# defaults on one core; the others take about a minute each.
#
# Two things keep that in check:
#
#   1. A run is skipped when the stored profile is already newer than the fitted
#      object it came from. Set PROFILE_FORCE <- TRUE to reprofile regardless.
#   2. Profiling is parallelised across parameters and the sampling is coarsened.
#      Parallelism is free accuracy-wise; the tolerance and evaluation-point
#      settings trade a little endpoint precision for a lot of time. All are
#      overridable by defining the variable before sourcing this script.
#
# Parallelism is per-parameter, so there is nothing to gain from more workers than
# the model has variance parameters (seven for dbh-ht). Each worker holds its own
# copy of the model, so memory, not core count, is usually the binding constraint.

if (!exists("PROFILE_FORCE"))  PROFILE_FORCE  <- FALSE
if (!exists("PROFILE_NCPUS"))  PROFILE_NCPUS  <- max(1L, min(6L, parallel::detectCores() - 1L))
if (!exists("PROFILE_MAXPTS")) PROFILE_MAXPTS <- 50      # lme4 default 100
if (!exists("PROFILE_DEVTOL")) PROFILE_DEVTOL <- 1e-7    # lme4 default 1e-9

# forking is unavailable on Windows; fall back to serial there
if (!exists("PROFILE_PARALLEL"))
  PROFILE_PARALLEL <- if (.Platform$OS.type == "windows") "no" else "multicore"
if (PROFILE_PARALLEL == "no") PROFILE_NCPUS <- 1L

# Which models to profile. Restrict this to run dbh-ht on its own in a long
# background job while the quick models are profiled in the foreground, e.g.
#   PROFILE_MODELS <- "dbh-ht"; source("scripts/r/profile.R")
if (!exists("PROFILE_MODELS"))
  PROFILE_MODELS <- c("age-dbh", "dbh-canopy", "age-ht", "dbh-ht")

profile_model <- function(label, fit_path, ci_path) {
  if (!label %in% PROFILE_MODELS) return(invisible(NA))

  if (!file.exists(fit_path))
    stop("missing fitted model: ", fit_path, " -- run the model script first")

  fit <- readRDS(fit_path)

  if (!PROFILE_FORCE && file.exists(ci_path) &&
      profile_matches_fit(readRDS(ci_path), fit)) {
    message(sprintf("Skipping %s: stored profile matches this fit.", label))
    return(invisible(FALSE))
  }

  message(sprintf("Profiling %s model (%s, ncpus = %d)...",
                  label, PROFILE_PARALLEL, PROFILE_NCPUS))
  t0  <- Sys.time()
  pr  <- profile(fit,
                 which    = "theta_",
                 maxpts   = PROFILE_MAXPTS,
                 devtol   = PROFILE_DEVTOL,
                 parallel = PROFILE_PARALLEL,
                 ncpus    = PROFILE_NCPUS)
  ci  <- confint(pr)
  attr(ci, "model_signature") <- model_signature(fit)
  saveRDS(ci, ci_path)
  message(sprintf("  done (%.1f min).",
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  invisible(TRUE)
}

profile_model("age-dbh",    "rslt/age-dbh/_dbh-age-model-obj.rds",       "rslt/age-dbh/ci_age_dbh_profile.rds")
profile_model("dbh-canopy", "rslt/dbh-canopy/_dbh-canopy-model-obj.rds", "rslt/dbh-canopy/ci_dbh_canopy_profile.rds")
profile_model("age-ht",     "rslt/age-ht/_age-ht-model-obj.rds",         "rslt/age-ht/ci_age_ht_profile.rds")
profile_model("dbh-ht",     "rslt/dbh-ht/_dbh-ht-model-obj.rds",         "rslt/dbh-ht/ci_dbh_ht_profile.rds")
