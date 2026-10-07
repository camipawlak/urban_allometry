## Shared by profile.R and figs-tbls.R.
##
## A stored profile is valid for the fit it was computed from. Deciding that by
## file mtime is unreliable: the fitted model objects are gitignored (they are
## large), so a fresh clone has the profile but not the model, and any re-run of a
## model script gives the object a newer mtime than the committed profile even
## when the refit is bit-identical. Profiling dbh-ht takes about a day, so a
## needless rerun is expensive.
##
## Compare the fit's identity instead. These five values determine the likelihood
## surface the profile was computed over; if they all match, the stored profile
## still applies no matter what the timestamps say.
model_signature <- function(fit) {
  list(
    theta   = unname(lme4::getME(fit, "theta")),
    beta    = unname(lme4::fixef(fit)),
    sigma   = unname(stats::sigma(fit)),
    nobs    = stats::nobs(fit),
    formula = paste(deparse(stats::formula(fit)), collapse = " ")
  )
}

## TRUE when `ci` carries a signature matching `fit`. Profiles written before
## signatures existed carry none; those return FALSE so they are recomputed once
## and gain a signature.
profile_matches_fit <- function(ci, fit, tol = 1e-8) {
  sig <- attr(ci, "model_signature")
  if (is.null(sig)) return(FALSE)
  isTRUE(all.equal(sig, model_signature(fit), tolerance = tol))
}
