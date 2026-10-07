library(tidyverse)
library(lme4)
library(modelr)
library(rlang)
library(ggforce)
library(writexl)
# The diagnostic plots below are for interactive inspection and are never saved.
# Under Rscript they would be auto-printed to a stray Rplots.pdf in the repo root,
# so send them to a null device unless the session is interactive. ggsave() opens
# its own device and is unaffected.
if (!interactive()) pdf(NULL)

data_dir <- 'data/_processed/'
plt_dir <- '_img/age-dbh/'
rslt_dir <- 'rslt/age-dbh/'
fs::dir_create(plt_dir)
fs::dir_create(rslt_dir)

## data viz --------------------------------------------------------------------

# load data partitions
load('data/_processed/wca-partitions.RData')

# Collapse to one record per tree before any modelling. The partition file retains
# repeated measures so that dbh-growth-rates.R can compute growth from them; those
# repeated measures must not enter a model fit. The seed is fixed so that every model
# script draws the same record per tree, and so that runs are reproducible.
set.seed(22626)
wca_train <- wca_train |> group_by(tree_id) |> slice_sample(n = 1) |> ungroup()
wca_val   <- wca_val   |> group_by(tree_id) |> slice_sample(n = 1) |> ungroup()

# drop trees failing the growth-rate quality screen (see make-partitions.R),
# anomalies, records missing age, and species with fewer than min_n observations
min_n <- 30
dbh_age_data <- wca_train |>
  filter(growth_qc_pass, age_dbh_height_is_anom == FALSE, !is.na(age)) |>
  add_count(spp, name = "n_spp") |>
  filter(n_spp >= min_n) |>
  select(-n_spp)

# inspect data for randomly chosen species (model specification checking)
dbh_age_data |>
  distinct(spp) |>
  slice_sample(n = 12) %>%
  semi_join(dbh_age_data, ., by = 'spp') |>
  ggplot(aes(x = log(age), y = log(dbh))) +
  facet_wrap(~spp) +
  geom_point(alpha = 0.2) +
  theme_bw()

# generate scatterplots for all species in series of 3x3 panels (alphabetical order)
plot_df <- dbh_age_data |>
  mutate(spp = factor(spp, levels = sort(unique(spp))))

n_pages <- plot_df |>
  distinct(spp) |>
  summarise(n = n()) |>
  pull(n) |>
  (\(n) ceiling(n / 9))()

base <- ggplot(plot_df, aes(x = log(age), y = log(height))) +
  geom_point(alpha = 0.3) +
  theme_bw()

spp_levels <- levels(plot_df$spp)

for (pg in seq_len(n_pages)) {
  i1 <- (pg - 1) * 9 + 1
  i2 <- min(pg * 9, length(spp_levels))
  spp_page <- spp_levels[i1:i2]
  
  a <- stringr::str_to_lower(stringr::str_sub(as.character(spp_page[1]), 1, 3))
  b <- stringr::str_to_lower(stringr::str_sub(as.character(spp_page[length(spp_page)]), 1, 3))
  
  p <- base + ggforce::facet_wrap_paginate(~ spp, nrow = 3, ncol = 3, page = pg)
  
  ggsave(
    filename = file.path(plt_dir, sprintf("dbh_age_spp_%s-%s_%03d.png", a, b, pg)),
    plot = p,
    width = 6, height = 6, units = "in", dpi = 300, bg = "white"
  )
  print(pg)
}

## mixed model fitting ---------------------------------------------------------

# check for within-genus replication
dbh_age_data |>
  distinct(genus, spp) |>
  count(genus) |> 
  arrange(desc(n)) |>
  filter(n > 2)

# fit random slope and intercept model
fit.lme <- lmer(log(dbh) ~ log(age) + (1 + log(age) | spp), 
                data = dbh_age_data, 
                control = lmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 2e6)))

# inspect
summary(fit.lme)

# export model object
write_rds(fit.lme, file = paste0(rslt_dir, '_dbh-age-model-obj.rds'))

## diagnostics -----------------------------------------------------------------

# scatterplots for residual diagnostics by species (same format as above)
diag_df <- dbh_age_data |>
  bind_cols(
    resid = residuals(fit.lme, scaled = TRUE),
    pred  = predict(fit.lme)
  ) |>
  group_by(spp) |>
  mutate(
    pred.scaled = as.numeric(scale(pred))
  ) |>
  ungroup() |>
  mutate(spp = factor(spp, levels = sort(unique(spp))))

n_pages <- diag_df |> 
  distinct(spp) |>
  summarise(n = n()) |>
  pull(n) |>
  (\(n) ceiling(n / 9))() # 3x3 => 9 spp per page

base_plot <- ggplot(diag_df, aes(x = pred.scaled, y = resid)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0) +
  labs(x = "fitted values", y = "residuals") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 6)
  )

spp_levels <- levels(diag_df$spp)
for (pg in seq_len(n_pages)) {
  i1 <- (pg - 1) * 9 + 1
  i2 <- min(pg * 9, length(spp_levels))
  spp_page <- spp_levels[i1:i2]
  
  a <- stringr::str_to_lower(stringr::str_sub(as.character(spp_page[1]), 1, 3))
  b <- stringr::str_to_lower(stringr::str_sub(as.character(spp_page[length(spp_page)]), 1, 3))
  
  p <- base_plot +
    ggforce::facet_wrap_paginate(~ spp, scales = "free", nrow = 3, ncol = 3, page = pg)
  
  ggsave(
    filename = file.path(plt_dir, sprintf("dbh-age-lmer-diagnostics_spp_%s-%s_%03d.png", a, b, pg)),
    plot = p,
    width = 6, height = 6, units = "in", dpi = 250, bg = "white"
  )
  
  print(pg)
}

## species-level coefficient blups + SEs --------------------------------------

# species-level blups (slope + intercept)
sp_blups <- coef(fit.lme)$spp %>%
  as.data.frame() %>%
  rownames_to_column("spp") %>%
  rename(intercept = `(Intercept)`, slope = `log(age)`)

# full covariance of blups/blues, conditional on variance parameter estimates
V <- as.matrix(vcov(fit.lme, full = TRUE))

# bookkeeping for random effects
p     <- length(fixef(fit.lme))   # number of fixed effects
Gp    <- getME(fit.lme, "Gp")     # block pointers into b (0-based)
cnms  <- getME(fit.lme, "cnms")   # term names per RE grouping factor
flist <- getME(fit.lme, "flist")  # grouping factors (levels/order)

# fixed-effect index in V
i_beta <- function(term) {
  i <- match(term, colnames(V))
  stopifnot(!is.na(i))
  i
}

# random-effect index in V
i_re <- function(grp, lvl, term) {
  j <- match(grp, names(cnms))
  stopifnot(!is.na(j))
  
  levs <- levels(flist[[j]])
  il <- match(lvl, levs)
  it <- match(term, cnms[[j]])
  stopifnot(!is.na(il), !is.na(it))
  
  nc <- length(cnms[[j]])
  start <- p + Gp[j] + 1
  start + (il - 1) * nc + (it - 1)
}

# SE for sum of components
se_sum <- function(idxs) {
  v <- V[idxs, idxs, drop = FALSE]
  sqrt(sum(v))
}

# covariance between two sets of components
cov_sum <- function(idxs_a, idxs_b) {
  sum(V[idxs_a, idxs_b, drop = FALSE])
}

# species-level SEs and intercept-slope covariance (beta + b_spp)
grp_name <- "spp"
spp_levels_vec <- levels(flist[[grp_name]])
sp_full_se <- tibble(
  spp                 = spp_levels_vec,
  intercept_se        = sapply(spp_levels_vec, \(s)
                               se_sum(c(i_beta("(Intercept)"),
                                        i_re(grp_name, s, "(Intercept)")))),
  slope_se            = sapply(spp_levels_vec, \(s)
                               se_sum(c(i_beta("log(age)"),
                                        i_re(grp_name, s, "log(age)")))),
  intercept_slope_cov = sapply(spp_levels_vec, \(s)
                               cov_sum(c(i_beta("(Intercept)"),
                                         i_re(grp_name, s, "(Intercept)")),
                                       c(i_beta("log(age)"),
                                         i_re(grp_name, s, "log(age)"))))
)

# join with point estimates
sp_blups_se <- sp_blups %>%
  left_join(select(sp_full_se, spp, intercept_se, slope_se, intercept_slope_cov),
            by = "spp")

# population-level estimates with SEs
## SE propagation ignores Cov(intercept, slope); sigma stored for prediction intervals
blues_tbl <- tibble(
  term = c("intercept", "slope", "sigma"),
  estimate = c(unname(fixef(fit.lme)["(Intercept)"]),
               unname(fixef(fit.lme)["log(age)"]),
               sigma(fit.lme)),
  se = c(sqrt(V[i_beta("(Intercept)"), i_beta("(Intercept)")]),
         sqrt(V[i_beta("log(age)"),    i_beta("log(age)")]),
         NA_real_)
)

# coefficient list
coef_list <- list(
  species = sp_blups_se,
  population = blues_tbl
)

# export as excel
writexl::write_xlsx(
  coef_list,
  path = file.path(rslt_dir, "dbh-age-coefs.xlsx")
)

# export as rds
write_rds(coef_list, file = paste0(rslt_dir, 'dbh-age-coefs.rds'))

## model visualizations (rates) ------------------------------------------------

# species-level scaling curve relative to population-level scaling curve
b0 <- blues_tbl$estimate[blues_tbl$term == "intercept"]
b1 <- blues_tbl$estimate[blues_tbl$term == "slope"]

age_small <- quantile(dbh_age_data$age, 0.25, na.rm = TRUE)
age_large <- quantile(dbh_age_data$age, 0.75, na.rm = TRUE)

s_dev <- sp_blups_se %>%
  transmute(
    spp,
    d_small = (intercept + slope * log(age_small)) - (b0 + b1 * log(age_small)),
    d_large = (intercept + slope * log(age_large)) - (b0 + b1 * log(age_large))
  ) %>%
  mutate(
    score = sqrt(d_small^2 + d_large^2),
    label = if_else(rank(-score) <= 3, spp, NA_character_)
  )

ggplot(s_dev, aes(x = d_small, y = d_large)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(alpha = 0.5, size = 1) +
  ggrepel::geom_text_repel(aes(label = label), na.rm = TRUE, size = 3) +
  coord_equal() +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.1, color = "black")
  ) +
  labs(x = paste0("log-DBH deviation at age ", round(age_small, 1)),
       y = paste0("log-DBH deviation at age ", round(age_large, 1)))

ggsave(
  filename = paste0(plt_dir, "dbh-age-lmer-spp-scaling-deviation.png"),
  width = 3, height = 3, dpi = 300
)

## model visualizations (species) ----------------------------------------------

# data in model scale
df_plot <- dbh_age_data %>%
  mutate(x = log(age), y = log(dbh))

# fixed (population) coefficients
b0 <- blues_tbl$estimate[blues_tbl$term == "intercept"]
b1 <- blues_tbl$estimate[blues_tbl$term == "slope"]

# population line (global x-grid)
x_grid <- modelr::seq_range(df_plot$x, n = 200)
pred_pop <- tibble(x = x_grid, fit = b0 + b1 * x)

# species coefficients + SEs
sb <- sp_blups_se %>%
  transmute(
    spp,
    b0s       = intercept,
    b1s       = slope,
    se_b0s    = intercept_se,
    se_b1s    = slope_se,
    cov_b0b1s = intercept_slope_cov
  )

# x-grid per species (within-species range) + predicted fit + SE band
x_rng_spp <- df_plot %>%
  group_by(spp) %>%
  summarise(
    xmin = quantile(x, 0, na.rm = TRUE),
    xmax = quantile(x, 1, na.rm = TRUE),
    .groups = "drop"
  )

pred_spp <- x_rng_spp %>%
  left_join(sb, by = "spp") %>%
  filter(!is.na(b0s), !is.na(b1s)) %>%
  rowwise() %>%
  mutate(x = list(seq(xmin, xmax, length.out = 100))) %>%
  tidyr::unnest(x) %>%
  mutate(
    fit = b0s + b1s * x,
    se  = sqrt(se_b0s^2 + 2 * x * cov_b0b1s + (x^2) * se_b1s^2),
    lo  = fit - 1.96 * se,
    hi  = fit + 1.96 * se
  ) %>%
  ungroup()

# paginate 3x3 (9 spp per page) + save one PNG per page (alphabetical)
df_plot <- df_plot %>%
  mutate(spp = factor(spp, levels = sort(unique(spp))))

pred_spp <- pred_spp %>%
  mutate(spp = factor(spp, levels = levels(df_plot$spp)))

spp_levels <- levels(df_plot$spp)
n_pages <- ceiling(length(spp_levels) / 9)

base_plot <- ggplot(df_plot, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 0.4, color = "grey") +
  geom_line(
    data = pred_pop,
    aes(y = fit, color = "Population"),
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  geom_ribbon(
    data = pred_spp,
    aes(y = fit, ymin = lo, ymax = hi),
    alpha = 0.2,
    fill = "blue"
  ) +
  geom_line(
    data = pred_spp,
    aes(y = fit, color = "Species"),
    linewidth = 0.5,
    alpha = 0.8
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Species" = "blue", "Population" = "black")
  ) +
  labs(x = "log(Age)", y = "log(DBH)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7),
    panel.grid.major = element_line(linewidth = 0.1),
    legend.position = "top"
  )

for (pg in seq_len(n_pages)) {
  i1 <- (pg - 1) * 9 + 1
  i2 <- min(pg * 9, length(spp_levels))
  spp_page <- spp_levels[i1:i2]
  
  a <- stringr::str_to_lower(stringr::str_sub(as.character(spp_page[1]), 1, 3))
  b <- stringr::str_to_lower(stringr::str_sub(as.character(spp_page[length(spp_page)]), 1, 3))
  
  p <- base_plot +
    ggforce::facet_wrap_paginate(~ spp, nrow = 3, ncol = 3, page = pg)
  
  ggsave(
    filename = file.path(plt_dir, sprintf("dbh-age-lmer-spp_%s-%s_%03d.png", a, b, pg)),
    plot = p,
    width = 6, height = 6, units = "in", dpi = 250, bg = "white"
  )
  print(pg)
}

## comparison with mcpherson model selection strategy --------------------------

# helper to compute aic based on implied likelihood (not on model scale)
aic_adj <- function(m) {
  mf <- model.frame(m)
  y_model <- model.response(mf)
  resp <- deparse(formula(m)[[2]])
  
  if (grepl("^\\s*log10\\s*\\(\\s*dbh\\s*\\)", resp)) {
    y_orig <- 10^(y_model)
    jac <- 2 * sum(log(y_orig)) + 2 * length(y_orig) * log(10)
  } else if (grepl("^\\s*log\\s*\\(\\s*dbh\\s*\\)", resp)) {
    y_orig <- exp(y_model)
    jac <- 2 * sum(log(y_orig))
  } else {
    y_orig <- y_model
    jac <- 0
  }
  
  w <- weights(m)
  if (is.null(w)) w <- rep(1, length(y_model))
  r <- resid(m)
  n <- length(r)
  rss_w <- sum(w * r^2)
  sigma2 <- rss_w / n
  k <- length(coef(m)) + 1  # betas + σ²
  
  neg2logL <- n * (1 + log(2 * pi)) + n * log(sigma2) - sum(log(w)) + jac
  2 * k + neg2logL
}

# no weights
weight_specs <- list(
  w1.none = ~ 1
)

# model specifications (age -> dbh; dbh is response)
model_specs <- list(
  lin    = dbh ~ age,
  exp    = log(dbh) ~ age,
  ln     = dbh ~ log(age),
  log10  = dbh ~ log(age, base = 10),
  loglog = log(dbh) ~ log(age)
)

# fit all model × weight combinations for each species
wfits <- dbh_age_data |>
  nest(.by = spp) |>
  crossing(model = names(model_specs),
           wname = names(weight_specs)) |>
  mutate(
    data = map2(data, wname, \(df, w_name) {
      df %>% mutate(.w = eval_tidy(f_rhs(weight_specs[[w_name]]), data = df))
    }),
    fit = map2(data, model, \(df, mod_name) {
      lm(model_specs[[mod_name]], data = df, weights = .w)
    }),
    aic = map_dbl(fit, aic_adj)
  )

# choose min aic
best_wfits <- wfits |>
  group_by(spp, model) |>
  slice_min(aic, n = 1) |>
  group_by(spp) |>
  slice_min(aic, n = 1)

# tabulate by model type
best_wfits |>
  ungroup() |>
  count(model) |>
  mutate(pct = n/sum(n)) |>
  arrange(desc(pct))

# check whether log-log models give comparable predictions where exp was selected
spp_exp <- best_wfits |>
  filter(model == 'exp') |>
  pull(spp)
rmspe_exp_ll <- wfits |>
  filter(spp %in% spp_exp,
         model %in% c('exp', 'loglog')) |>
  select(spp, model, fit) |>
  spread(model, fit) |>
  mutate(test = map(spp, 
                    ~filter(wca_val, spp == .x, age_dbh_height_is_anom == F) |> 
                      select(age, dbh) |>
                      drop_na(age, dbh)),
         exp.pred = map2(exp, test, predict),
         loglog.pred = map2(loglog, test, predict),
         log.dbh = map(test, ~ log(.x$dbh))) |>
  select(spp, exp.pred, loglog.pred, log.dbh) |>
  unnest(everything()) |>
  mutate(exp.res = exp.pred - log.dbh,
         loglog.res = loglog.pred - log.dbh) |>
  group_by(spp) |>
  summarize(across(ends_with('.res'), ~sqrt(mean(.x^2)))) |>
  mutate(delta = loglog.res - exp.res)

# Export the exp-vs-loglog RMSPE comparison. This is a different comparison from
# dbh-ht's rmspe file: that one scores the model against the McPherson equations,
# whereas this scores two functional forms of our own fit against each other on
# the validation partition. exp.res / loglog.res are already RMSPE, not residuals.
write_csv(
  rmspe_exp_ll |> rename(exp_rmspe = exp.res, loglog_rmspe = loglog.res),
  paste0(rslt_dir, 'age-dbh-rmspe-exp-vs-loglog.csv')
)

# looks mostly comparable
rmspe_exp_ll |>
  arrange(desc(delta)) |>
  ggplot(aes(x = delta, y = fct_reorder(spp, delta))) +
  geom_point() +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = 'dashed') +
  theme_bw() +
  labs(x = '% change RMSPE (ll:exp)',
       y = NULL) 

# pull exceptions and plot; most still look nonlinear, so go with log-log
spp_check <- rmspe_exp_ll |>
  filter(delta > 0.1) |>
  pull(spp)
dbh_age_data |>
  filter(spp %in% spp_check) |>
  ggplot(aes(x = (age), y = log(dbh))) +
  facet_wrap(~spp, nrow = 5, ncol = 5) +
  geom_point() +
  theme_bw()

# helper to back-transform dbh from model scale
bt_dbh <- function(fit, eta) {
  resp <- deparse(formula(fit)[[2]])
  if (grepl("^\\s*log10\\s*\\(\\s*dbh\\s*\\)", resp)) {
    10^(eta)
  } else if (grepl("^\\s*log\\s*\\(\\s*dbh\\s*\\)", resp)) {
    exp(eta)
  } else {
    eta
  }
}

# restrict validation set to species present in training
train_spp <- dbh_age_data %>%
  distinct(spp) %>%
  pull(spp)

val_filt <- wca_val %>%
  filter(
    spp %in% train_spp,
    age_dbh_height_is_anom == FALSE,
    !is.na(age)
  )

# oos predictions from individual models (McPherson-style)
pred_df1 <- best_wfits %>%
  filter(spp %in% train_spp) %>%
  left_join(
    val_filt %>% nest(.key = "test", .by = spp),
    join_by(spp)
  ) %>%
  mutate(
    pred.raw = map2(fit, test, predict),
    pred = map2(fit, pred.raw, bt_dbh)
  ) %>%
  select(spp, test, pred) %>%
  unnest(everything()) %>%
  mutate(method = "mcpherson") %>%
  select(spp, tree_id, age, dbh, method, pred)

# predictions from mixed model
pred_df2 <- val_filt %>%
  add_predictions(fit.lme) %>%
  mutate(
    method = "lmer",
    pred = exp(pred)   # because model is log(dbh) ~ log(age) + ...
  ) %>%
  select(spp, tree_id, age, dbh, method, pred)

# bind
pred_df <- bind_rows(pred_df1, pred_df2)

# plot RMSPE comparison
top_n <- 8
plot_df <- pred_df |>
  group_by(method, spp) |>
  summarize(mspe = mean((dbh - pred)^2), .groups = "drop") |>
  tidyr::pivot_wider(names_from = method, values_from = mspe) |>
  mutate(
    diff = abs(log10(lmer) - log10(mcpherson)),
    label = if_else(rank(-diff) <= top_n, spp, NA_character_)
  )

ggplot(plot_df, aes(x = lmer, y = mcpherson)) +
  geom_point(alpha = 0.75, aes(color = factor(!is.na(label)))) +
  geom_abline(slope = 1, linetype = "dashed") +
  ggrepel::geom_text_repel(
    aes(label = label),
    na.rm = TRUE, size = 3,
    hjust = 1, direction = "y", nudge_x = -0.05
  ) +
  scale_y_log10() +
  scale_x_log10() +
  guides(color = "none") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.1, color = "black")
  ) +
  labs(x = "mixed model RMSPE", y = "species model RMSPE")

paste0(plt_dir, "dbh-age-mcpherson-lmer-rmspe-comparison.png") |>
  ggsave(width = 3, height = 3, units = "in", dpi = 300)

## species-level plots (back-transformed, one png per species) -----------------

spp_plt_dir <- '_figures/age-dbh/'
# Clear the directory first: this loop writes one PNG per fitted species, and a
# run that fits fewer species than the last one would otherwise leave the old
# files behind, indistinguishable by name from current output.
fs::dir_delete(spp_plt_dir[fs::dir_exists(spp_plt_dir)])
fs::dir_create(spp_plt_dir)

sigma_val <- blues_tbl$estimate[blues_tbl$term == "sigma"]

for (sp in sort(unique(dbh_age_data$spp))) {
  sp_data <- dbh_age_data |> filter(spp == sp)
  n_obs   <- nrow(sp_data)
  sp_coef <- sp_blups_se |> filter(spp == sp)
  if (nrow(sp_coef) == 0) next
  
  x_seq   <- seq(min(sp_data$age), max(sp_data$age), length.out = 200)
  log_x   <- log(x_seq)
  fit_log <- sp_coef$intercept + sp_coef$slope * log_x
  se_ci   <- sqrt(sp_coef$intercept_se^2 + 2 * log_x * sp_coef$intercept_slope_cov + log_x^2 * sp_coef$slope_se^2)
  se_pi   <- sqrt(sp_coef$intercept_se^2 + 2 * log_x * sp_coef$intercept_slope_cov + log_x^2 * sp_coef$slope_se^2 + sigma_val^2)
  
  pred_df <- tibble(
    x     = x_seq,
    fit   = exp(fit_log),
    ci_lo = exp(fit_log - 1.96 * se_ci),
    ci_hi = exp(fit_log + 1.96 * se_ci),
    pi_lo = exp(fit_log - 1.96 * se_pi),
    pi_hi = exp(fit_log + 1.96 * se_pi)
  )
  
  eq_label <- sprintf("dbh = %.3f \u00d7 age^%.3f",
                      exp(sp_coef$intercept), sp_coef$slope)
  
  p <- ggplot(sp_data, aes(x = age, y = dbh)) +
    geom_point(position = position_jitter(width = 0.5, height = 0, seed = 42),
               color = "grey70", alpha = 0.4, size = 0.8) +
    geom_ribbon(data = pred_df,
                aes(x = x, ymin = ci_lo, ymax = ci_hi),
                inherit.aes = FALSE,
                fill = "#2E8B57", alpha = 0.2) +
    geom_line(data = pred_df, aes(x = x, y = pi_lo),
              inherit.aes = FALSE,
              color = "#2E8B57", linetype = "dashed", linewidth = 0.5) +
    geom_line(data = pred_df, aes(x = x, y = pi_hi),
              inherit.aes = FALSE,
              color = "#2E8B57", linetype = "dashed", linewidth = 0.5) +
    geom_line(data = pred_df, aes(x = x, y = fit),
              inherit.aes = FALSE,
              color = "#145A32", linewidth = 0.8) +
    annotate("text", x = -Inf, y = Inf,
             label = eq_label,
             hjust = -0.05, vjust = 1.5,
             size = 3, color = "grey30") +
    labs(title = sprintf("%s  (n = %d)", sp, n_obs),
         x = "Age (years)", y = "DBH (cm)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.1))
  
  sp_safe <- str_replace_all(sp, "[^A-Za-z0-9]+", "_") |> str_remove("^_+|_+$")
  ggsave(
    filename = file.path(spp_plt_dir, sprintf("%s.png", sp_safe)),
    plot = p,
    width = 4, height = 4, units = "in", dpi = 200, bg = "white"
  )
}