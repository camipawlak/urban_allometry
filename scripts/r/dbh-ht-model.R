library(tidyverse)
library(lme4)
library(modelr)
library(rlang)
library(ggforce)
library(writexl)
source("scripts/r/fig-theme.R")
# The diagnostic plots below are for interactive inspection and are never saved.
# Under Rscript they would be auto-printed to a stray Rplots.pdf in the repo root,
# so send them to a null device unless the session is interactive. ggsave() opens
# its own device and is unaffected.
if (!interactive()) pdf(NULL)

data_dir <- 'data/_processed/'
plt_dir <- '_img/dbh-ht/'
rslt_dir <- 'rslt/dbh-ht/'
fs::dir_create(plt_dir)
fs::dir_create(rslt_dir)

## data viz --------------------------------------------------------------------
# load data partitions
load('data/_processed/wca-partitions.RData')

# NOTE: genus is now resolved in make-partitions.R, where the column is created.

# Collapse to one record per tree before any modelling. The partition file retains
# repeated measures so that dbh-growth-rates.R can compute growth from them; those
# repeated measures must not enter a model fit. The seed is fixed so that every model
# script draws the same record per tree, and so that runs are reproducible.
set.seed(22626)
wca_train <- wca_train |> group_by(tree_id) |> slice_sample(n = 1) |> ungroup()
wca_val   <- wca_val   |> group_by(tree_id) |> slice_sample(n = 1) |> ungroup()

# drop anomalies
dbh_ht_data <- wca_train |>
  filter(dbh_height_is_anom == F)

# inspect data for randomly chosen species (model specification checking)
dbh_ht_data |>
  distinct(spp) |>
  slice_sample(n = 12) %>%
  semi_join(dbh_ht_data, ., by = 'spp') |>
  ggplot(aes(x = log(dbh), y = log(height))) +
  facet_wrap(~spp) +
  geom_point(alpha = 0.2) +
  theme_bw()

# generate scatterplots for all species in series of 3x3 panels (alphabetical order)
plot_df <- dbh_ht_data |>
  mutate(spp = factor(spp, levels = sort(unique(spp))))

n_pages <- plot_df |>
  distinct(spp) |>
  summarise(n = n()) |>
  pull(n) |>
  (\(n) ceiling(n / 9))()

base <- ggplot(plot_df, aes(x = log(dbh), y = log(height))) +
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
    filename = file.path(plt_dir, sprintf("dbh_ht_spp_%s-%s_%03d.png", a, b, pg)),
    plot = p,
    width = 6, height = 6, units = "in", dpi = 300, bg = "white"
  )
  print(pg)
}


## mixed model fitting ---------------------------------------------------------

# fit random slope and intercept model
fit.lme <- lmer(log(height) ~ log(dbh) + (1 + log(dbh) | genus/spp), 
                data = dbh_ht_data, 
                control = lmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 2e6)))

# inspect
summary(fit.lme)

# export model object
write_rds(fit.lme, file = paste0(rslt_dir, '_dbh-ht-model-obj.rds'))

## diagnostics -----------------------------------------------------------------

# scatterplots for residual diagnostics by species (same format as above)
diag_df <- dbh_ht_data |>
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
    filename = file.path(plt_dir, sprintf("dbh-ht-lmer-diagnostics_spp_%s-%s_%03d.png", a, b, pg)),
    plot = p,
    width = 6, height = 6, units = "in", dpi = 250, bg = "white"
  )
  
  print(pg)
}

## species-level and genus-level coefficient blups -----------------------------

# species-level blups (slope + intercept)
sp_blups <- coef(fit.lme)$`spp:genus` %>%
  as.data.frame() %>%
  rownames_to_column("spp_genus") %>%
  separate(spp_genus, into = c("spp", "genus"), sep = ":", remove = FALSE) %>%
  rename(intercept = `(Intercept)`, slope = `log(dbh)`)

# genus-level blups (slope + intercept)
genus_blups <- coef(fit.lme)$genus %>%
  as.data.frame() %>%
  rownames_to_column("genus") %>%
  rename(intercept = `(Intercept)`, slope = `log(dbh)`)

# full covariance of blups/blues, conditional on variance parameter estimates
V <- as.matrix(vcov(fit.lme, full = TRUE))

# bookkeeping for random effects
p     <- length(fixef(fit.lme))                 # number of fixed effects
Gp    <- getME(fit.lme, "Gp")                   # block pointers into b (0-based)
cnms  <- getME(fit.lme, "cnms")                 # term names per RE grouping factor
flist <- getME(fit.lme, "flist")                # grouping factors (levels/order)

# fixed-effect index in V
i_beta <- function(term) {
  i <- match(term, colnames(V))
  stopifnot(!is.na(i))
  i
}

# random-effect index in V
i_re <- function(grp, lvl, term) {
  j  <- match(grp, names(cnms))                 # which RE grouping factor
  stopifnot(!is.na(j))
  
  levs <- levels(flist[[j]])
  il   <- match(lvl, levs)                      # which level within grp
  it   <- match(term, cnms[[j]])                # which term within grp
  stopifnot(!is.na(il), !is.na(it))
  
  nc    <- length(cnms[[j]])                    # terms per level
  start <- p + Gp[j] + 1                        # convert 0-based pointer to 1-based V index
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

# genus-level SEs and intercept-slope covariance (beta + b_genus)
genus_levels <- levels(flist[["genus"]])
genus_full_se <- tibble(
  genus               = genus_levels,
  intercept_se        = sapply(genus_levels, \(g)
                               se_sum(c(i_beta("(Intercept)"),
                                        i_re("genus", g, "(Intercept)")))),
  slope_se            = sapply(genus_levels, \(g)
                               se_sum(c(i_beta("log(dbh)"),
                                        i_re("genus", g, "log(dbh)")))),
  intercept_slope_cov = sapply(genus_levels, \(g)
                               cov_sum(c(i_beta("(Intercept)"),
                                         i_re("genus", g, "(Intercept)")),
                                       c(i_beta("log(dbh)"),
                                         i_re("genus", g, "log(dbh)"))))
)

# genera with at least three species
keep_genus <- dbh_ht_data %>%
  distinct(genus, spp) %>%
  count(genus, name = "n_spp") %>%
  filter(n_spp >= 3) %>%
  pull(genus)

# join with point estimates and drop genera with fewer than three species
genus_blups_se <- genus_blups %>%
  left_join(genus_full_se, by = "genus") |>
  filter(genus %in% keep_genus)

# species-level SEs (beta + b_genus + b_spp:genus)
int_name <- "spp:genus"
int_levels <- levels(flist[[int_name]])
get_genus_from_int <- function(gs) {
  parts <- strsplit(gs, ":", fixed = TRUE)[[1]]
  if (int_name == "spp:genus") parts <- rev(parts)
  parts[1]
}

sp_full_se <- tibble(
  spp_genus           = int_levels,
  genus               = vapply(int_levels, get_genus_from_int, character(1)),
  intercept_se        = sapply(int_levels, \(gs) {
    g <- get_genus_from_int(gs)
    se_sum(c(i_beta("(Intercept)"),
             i_re("genus", g, "(Intercept)"),
             i_re(int_name, gs, "(Intercept)")))
  }),
  slope_se            = sapply(int_levels, \(gs) {
    g <- get_genus_from_int(gs)
    se_sum(c(i_beta("log(dbh)"),
             i_re("genus", g, "log(dbh)"),
             i_re(int_name, gs, "log(dbh)")))
  }),
  intercept_slope_cov = sapply(int_levels, \(gs) {
    g <- get_genus_from_int(gs)
    cov_sum(c(i_beta("(Intercept)"),
              i_re("genus", g, "(Intercept)"),
              i_re(int_name, gs, "(Intercept)")),
            c(i_beta("log(dbh)"),
              i_re("genus", g, "log(dbh)"),
              i_re(int_name, gs, "log(dbh)")))
  })
)

# join with point estimates
sp_blups_se <- sp_blups %>%
  left_join(select(sp_full_se, spp_genus, intercept_se, slope_se, intercept_slope_cov),
            by = "spp_genus")

# population-level estimates with SEs
## SE propagation ignores Cov(intercept, slope); sigma stored for prediction intervals
blues_tbl <- tibble(
  term = c("intercept", "slope", "sigma"),
  estimate = c(unname(fixef(fit.lme)["(Intercept)"]),
               unname(fixef(fit.lme)["log(dbh)"]),
               sigma(fit.lme)),
  se = c(sqrt(V[i_beta("(Intercept)"), i_beta("(Intercept)")]),
         sqrt(V[i_beta("log(dbh)"),    i_beta("log(dbh)")]),
         NA_real_)
)

# coefficient list
coef_list <- list(
  species = sp_blups_se,
  genus = genus_blups_se,
  population = blues_tbl
)

# export as excel
writexl::write_xlsx(
  coef_list,
  path = file.path(rslt_dir, "dbh-ht-coefs.xlsx")
)

# export as rds
write_rds(coef_list, file = paste0(rslt_dir, 'dbh-ht-coefs.rds'))

## model visualizations (rates) ------------------------------------------------

# genus-level scaling rates
g_rates <- genus_blups_se %>%
  transmute(
    genus,
    slope,
    slope_se
  ) %>%
  mutate(genus = fct_reorder(genus, slope))

ggplot(g_rates, aes(x = slope, y = genus)) +
  geom_errorbarh(aes(xmin = slope - 2*slope_se,
                     xmax = slope + 2*slope_se),
                 height = 0) +
  geom_point(size = 1.6) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.1, color = "black")) +
  labs(x = "Scaling rate (DBH \u2192 height)", y = NULL)

ggsave(filename = paste0(plt_dir, 'dbh-ht-lmer-genus-scaling-rates.png'),
       width = 3, height = 3, dpi = 300)

# species-level scaling curve relative to population-level scaling curve
b0 <- blues_tbl$estimate[blues_tbl$term == "intercept"]
b1 <- blues_tbl$estimate[blues_tbl$term == "slope"]
dbh_small <- quantile(dbh_ht_data$dbh, 0.25, na.rm = TRUE)
dbh_large <- quantile(dbh_ht_data$dbh, 0.75, na.rm = TRUE)

s_dev <- sp_blups_se %>%
  transmute(
    spp,
    d_small = (intercept + slope * log(dbh_small)) - (b0 + b1 * log(dbh_small)),
    d_large = (intercept + slope * log(dbh_large)) - (b0 + b1 * log(dbh_large))
  ) %>%
  mutate(
    score = sqrt(d_small^2 + d_large^2),
    label = if_else(rank(-score) <= 3, spp, NA_character_) # label more or fewer points
  )

ggplot(s_dev, aes(x = d_small, y = d_large)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(alpha = 0.5, size = 1) +
  ggrepel::geom_text_repel(aes(label = label), na.rm = TRUE, size = 3) +
  coord_equal() +
  theme_allometry() +
  labs(x = paste0("log-height deviation at DBH ", round(dbh_small, 1)),
       y = paste0("log-height deviation at DBH ", round(dbh_large, 1)))

# square scatter, single column
save_fig(last_plot(), paste0(plt_dir, 'dbh-ht-lmer-spp-scaling-rates.png'),
         width = FIG_WIDTH_HALF, height = 3.5)

## model visualizations (genus) ------------------------------------------------

# data in model scale
df_plot <- dbh_ht_data %>%
  mutate(x = log(dbh), y = log(height)) |>
  filter(genus %in% keep_genus)

# fixed (population) coefficients
b0 <- blues_tbl$estimate[blues_tbl$term == "intercept"]
b1 <- blues_tbl$estimate[blues_tbl$term == "slope"]

# genus coefficients + SEs
gb <- genus_blups_se %>%
  rename(
    b0g       = intercept,
    b1g       = slope,
    se_b0g    = intercept_se,
    se_b1g    = slope_se,
    cov_b0b1g = intercept_slope_cov
  )

# x-grid per genus (use within-genus range)
x_rng <- df_plot %>%
  group_by(genus) %>%
  summarise(xmin = quantile(x, 0, na.rm = TRUE),
            xmax = quantile(x, 1, na.rm = TRUE),
            .groups = "drop")

pred_genus <- x_rng %>%
  left_join(gb, by = "genus") %>%
  rowwise() %>%
  mutate(x = list(seq(xmin, xmax, length.out = 100))) %>%
  unnest(x) %>%
  mutate(
    fit = b0g + b1g * x,
    se  = sqrt(se_b0g^2 + 2 * x * cov_b0b1g + (x^2) * se_b1g^2),
    lo  = fit - 1.96 * se,
    hi  = fit + 1.96 * se
  ) %>%
  ungroup()

# population line on same x-range as overall data
x_grid <- seq_range(df_plot$x, n = 200)
pred_pop <- tibble(x = x_grid, fit = b0 + b1 * x)

base_plot <- ggplot(df_plot, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 0.4, color = "grey") +
  geom_line(
    data = pred_pop,
    aes(y = fit, color = "Population"),
    linewidth = 0.3,
    linetype = 'dashed'
  ) +
  geom_ribbon(
    data = pred_genus,
    aes(y = fit, ymin = lo, ymax = hi),
    alpha = 0.2,
    fill = "blue"
  ) +
  geom_line(
    data = pred_genus,
    aes(y = fit, color = "Genus"),
    linewidth = 0.5,
    alpha = 0.8
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Genus" = "blue", "Population" = "black")
  ) +
  labs(x = "log(DBH)", y = "log(Height)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7),
    panel.grid.major = element_line(linewidth = 0.1),
    legend.position = 'top'
  )

# # for testing
# base_plot +
#   facet_wrap_paginate(~ genus, nrow = 3, ncol = 3, page = 1)

# paginate 3x3 (9 genera per page) + save one PNG per page
n_pages <- df_plot %>%
  distinct(genus) %>%
  summarise(n = n()) %>%
  pull(n) %>%
  (\(n) ceiling(n / 9))()

genus_levels <- levels(factor(df_plot$genus, levels = sort(unique(df_plot$genus))))

for (pg in seq_len(n_pages)) {
  i1 <- (pg - 1) * 9 + 1
  i2 <- min(pg * 9, length(genus_levels))
  genus_page <- genus_levels[i1:i2]
  
  a <- stringr::str_to_lower(stringr::str_sub(as.character(genus_page[1]), 1, 3))
  b <- stringr::str_to_lower(stringr::str_sub(as.character(genus_page[length(genus_page)]), 1, 3))
  
  p <- base_plot +
    ggforce::facet_wrap_paginate(~ genus, nrow = 3, ncol = 3, page = pg)
  
  ggsave(
    filename = file.path(plt_dir, sprintf("dbh-ht-lmer-genus_%s-%s_%03d.png", a, b, pg)),
    plot = p,
    width = 5, height = 6, units = "in", dpi = 250, bg = "white"
  )
  print(pg)
}

## model visualizations (species) ----------------------------------------------

# data in model scale
df_plot <- dbh_ht_data %>% mutate(x = log(dbh), y = log(height))

# fixed (population) coefficients
b0 <- blues_tbl$estimate[blues_tbl$term == "intercept"]
b1 <- blues_tbl$estimate[blues_tbl$term == "slope"]

# population line (global x-grid)
x_grid_pop <- seq_range(df_plot$x, n = 200)
pred_pop <- tibble(x = x_grid_pop, fit = b0 + b1 * x)

# genus coefficients (expects: genus, intercept, slope)
gb <- genus_blups_se %>% transmute(genus, g_int = intercept, g_slope = slope)

# species coefficients (use existing sp_blups; expects: genus, spp, intercept, slope)
has_se <- all(c("intercept_se", "slope_se") %in% names(sp_blups_se))

# genus replication
genus_sizes <- df_plot %>%
  distinct(genus, spp) %>%
  count(genus, name = "n_spp")

g3 <- genus_sizes %>% filter(n_spp >= 3) %>% arrange(genus) %>% pull(genus)
g_lt3 <- genus_sizes %>% filter(n_spp < 3) %>% arrange(genus) %>% pull(genus)

# genera with >= 3 species: paginate within genus (facets = spp)
for (g in g3) {
  df_g <- df_plot %>% filter(genus == g)
  
  spp_levels <- df_g %>% distinct(spp) %>% arrange(spp) %>% pull(spp)
  n_pages <- ceiling(length(spp_levels) / 9)
  
  # genus line (within-genus x-range)
  x_grid_g <- seq_range(df_g$x, n = 200)
  pred_genus_g <- gb %>%
    filter(genus == g) %>%
    crossing(x = x_grid_g) %>%
    mutate(fit = g_int + g_slope * x)
  
  # species predictions (+ ribbon if SEs present), within genus
  sp_g <- sp_blups_se %>% filter(genus == g)
  
  x_rng_spp <- df_g %>%
    group_by(spp) %>%
    summarise(xmin = quantile(x, 0, na.rm = TRUE),
              xmax = quantile(x, 1, na.rm = TRUE),
              .groups = "drop")
  
  pred_spp_g <- x_rng_spp %>%
    left_join(sp_g, by = "spp") %>%
    filter(!is.na(intercept), !is.na(slope)) %>%
    rowwise() %>%
    mutate(x = list(seq(xmin, xmax, length.out = 100))) %>%
    unnest(x) %>%
    mutate(fit = intercept + slope * x) %>%
    ungroup()
  
  if (has_se) {
    pred_spp_g <- pred_spp_g %>%
      mutate(
        se = sqrt(intercept_se^2 + 2 * x * intercept_slope_cov + (x^2) * slope_se^2),
        lo = fit - 1.96 * se,
        hi = fit + 1.96 * se
      )
  }
  
  # facet order
  df_g <- df_g %>% mutate(spp = factor(spp, levels = spp_levels))
  pred_spp_g <- pred_spp_g %>% mutate(spp = factor(spp, levels = spp_levels))
  
  g_safe <- str_replace_all(g, "[^A-Za-z0-9]+", "_") %>% str_remove("^_+|_+$")
  
  base_plot_g <- ggplot(df_g, aes(x = x, y = y)) +
    geom_point(alpha = 0.6, size = 0.4, color = "grey") +
    geom_line(data = pred_pop,
              aes(y = fit, color = "Population", linetype = "Population"),
              linewidth = 0.3) +
    geom_line(data = pred_genus_g,
              aes(y = fit, color = "Genus", linetype = "Genus"),
              linewidth = 0.5, alpha = 0.8) +
    { if (has_se)
      geom_ribbon(data = pred_spp_g,
                  aes(y = fit, ymin = lo, ymax = hi),
                  alpha = 0.2, fill = "blue") } +
    geom_line(data = pred_spp_g,
              aes(y = fit, color = "Species", linetype = "Species"),
              linewidth = 0.6, alpha = 0.9) +
    scale_color_manual(name = NULL,
                       values = c("Population" = "black", "Genus" = "blue", "Species" = "blue")) +
    scale_linetype_manual(name = NULL,
                          values = c("Population" = "dashed", "Genus" = "dotted", "Species" = "solid")) +
    labs(x = "log(DBH)", y = "log(Height)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(size = 7),
          panel.grid.major = element_line(linewidth = 0.1),
          legend.position = "top")
  
  # # for testing
  # base_plot_g + facet_wrap_paginate(~ spp, nrow = 3, ncol = 3, page = 1)
  
  for (pg in seq_len(n_pages)) {
    p <- base_plot_g + facet_wrap_paginate(~ spp, nrow = 3, ncol = 3, page = pg)
    
    ggsave(
      filename = file.path(plt_dir, sprintf("dbh-ht-lmer-spp-%s_%03d.png", g_safe, pg)),
      plot = p,
      width = 6, height = 6, units = "in", dpi = 250, bg = "white"
    )
    print(paste("Saved:", g, "page", pg))
  }
}

# all other genera (< 3 species): paginate by species globally (alphabetical)
df_s <- df_plot %>% filter(genus %in% g_lt3)

# alphabetical species list (across the <3-species genera)
spp_levels_s <- df_s %>% distinct(spp) %>% arrange(spp) %>% pull(spp)
df_s <- df_s %>% mutate(spp = factor(spp, levels = spp_levels_s))

# species predictions (+ ribbon), using sp_blups
sp_s <- sp_blups_se %>% filter(genus %in% g_lt3)

x_rng_s <- df_s %>%
  group_by(genus, spp) %>%
  summarise(xmin = quantile(x, 0, na.rm = TRUE),
            xmax = quantile(x, 1, na.rm = TRUE),
            .groups = "drop")

pred_spp_s <- x_rng_s %>%
  left_join(sp_s, by = c("genus", "spp")) %>%
  filter(!is.na(intercept), !is.na(slope)) %>%
  rowwise() %>%
  mutate(x = list(seq(xmin, xmax, length.out = 100))) %>%
  unnest(x) %>%
  mutate(fit = intercept + slope * x) %>%
  ungroup()

if (has_se) {
  pred_spp_s <- pred_spp_s %>%
    mutate(
      se = sqrt(intercept_se^2 + 2 * x * intercept_slope_cov + (x^2) * slope_se^2),
      lo = fit - 1.96 * se,
      hi = fit + 1.96 * se
    )
}

# genus line per panel: join genus coef onto each species' genus and generate within-species x-range
pred_genus_s <- x_rng_s %>%
  left_join(gb, by = "genus") %>%
  rowwise() %>%
  mutate(x = list(seq(xmin, xmax, length.out = 100))) %>%
  unnest(x) %>%
  mutate(fit = g_int + g_slope * x) %>%
  ungroup()

pred_spp_s <- pred_spp_s %>% mutate(spp = factor(spp, levels = spp_levels_s))
n_pages_s <- ceiling(length(spp_levels_s) / 9)

base_plot_s <- ggplot(df_s, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 0.4, color = "grey") +
  geom_line(data = pred_pop,
            aes(y = fit, color = "Population", linetype = "Population"),
            linewidth = 0.3) +
  { if (has_se)
    geom_ribbon(data = pred_spp_s,
                aes(y = fit, ymin = lo, ymax = hi),
                alpha = 0.2, fill = "blue") } +
  geom_line(data = pred_spp_s,
            aes(y = fit, color = "Species", linetype = "Species"),
            linewidth = 0.6, alpha = 0.9) +
  scale_color_manual(name = NULL,
                     values = c("Population" = "black", "Species" = "blue")) +
  scale_linetype_manual(name = NULL,
                        values = c("Population" = "dashed", "Species" = "solid")) +
  labs(x = "log(DBH)", y = "log(Height)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 7),
        panel.grid.major = element_line(linewidth = 0.1),
        legend.position = "top")

for (pg in seq_len(n_pages_s)) {
  i1 <- (pg - 1) * 9 + 1
  i2 <- min(pg * 9, length(spp_levels_s))
  spp_page <- spp_levels_s[i1:i2]
  
  a <- str_to_lower(str_sub(as.character(spp_page[1]), 1, 3))
  b <- str_to_lower(str_sub(as.character(spp_page[length(spp_page)]), 1, 3))
  
  p <- base_plot_s + facet_wrap_paginate(~ spp, nrow = 3, ncol = 3, page = pg)
  
  ggsave(
    filename = file.path(plt_dir, sprintf("dbh-ht-lmer-spp_singles_%s-%s_%03d.png", a, b, pg)),
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
  
  # Identify the transformation and compute Jacobian term
  if (grepl("^\\s*log10\\s*\\(\\s*height\\s*\\)", resp)) {
    # Model fit on log10(height): y_model = log10(height)
    y_orig <- 10^(y_model)
    jac <- 2 * sum(log(y_orig)) + 2 * length(y_orig) * log(10)
  } else if (grepl("^\\s*log\\s*\\(\\s*height\\s*\\)", resp)) {
    # Model fit on natural log(height)
    y_orig <- exp(y_model)
    jac <- 2 * sum(log(y_orig))
  } else {
    # Model fit on original scale
    y_orig <- y_model
    jac <- 0
  }
  
  # Weighted least squares pieces on model scale
  w <- weights(m)
  if (is.null(w)) w <- rep(1, length(y_model))
  r <- resid(m)
  n <- length(r)
  rss_w <- sum(w * r^2)
  sigma2 <- rss_w / n
  k <- length(coef(m)) + 1  # betas + σ²
  
  # -2 log L on original scale, with weight + Jacobian corrections
  neg2logL <- n * (1 + log(2 * pi)) + n * log(sigma2) - sum(log(w)) + jac
  aic_orig <- 2 * k + neg2logL
  aic_orig
}

# no weights
weight_specs <- list(
  w1.none = ~ 1
)

# model specifications (COPIED FROM PY SCRIPT)
model_specs <- list(
  lin    = height ~ dbh,
  exp    = log(height) ~ dbh,
  ln     = height ~ log(dbh),
  log10  = height ~ log(dbh, base = 10),
  loglog = log(height) ~ log(dbh)
)

# fit all model × weight combinations for each species
wfits <- dbh_ht_data |>
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

# helper to back-transform height from model scale
bt_height <- function(fit, eta) {
  resp <- deparse(formula(fit)[[2]])
  if (grepl("^\\s*log10\\s*\\(\\s*height\\s*\\)", resp)) {
    10^(eta)
  } else if (grepl("^\\s*log\\s*\\(\\s*height\\s*\\)", resp)) {
    exp(eta)
  } else {
    eta
  }
}

# oos predictions from individual models
pred_df1 <- best_wfits |>
  left_join(wca_val |> filter(dbh_height_is_anom == F) |> nest(.key = 'test', .by = spp),
            join_by(spp)) |>
  mutate(pred.raw = map2(fit, test, predict),
         pred = map2(fit, pred.raw, bt_height)) |>
  select(spp, test, pred) |>
  unnest(everything()) |>
  mutate(method = 'mcpherson') |>
  select(spp, tree_id, dbh, height, method, pred)

# predictions from mixed model
pred_df2 <- wca_val |>
  filter(dbh_height_is_anom == F) |>
  add_predictions(fit.lme) |>
  mutate(method = 'lmer',
         pred = exp(pred)) |>
  select(spp, tree_id, dbh, height, method, pred)

# bind
pred_df <- bind_rows(pred_df1, pred_df2) 

# plot
top_n <- 2 # num. largest differences to label
plot_df <- pred_df |>
  group_by(method, spp) |>
  summarize(mspe = mean((height - pred)^2), .groups = "drop") |>
  arrange(spp, mspe) |>
  tidyr::pivot_wider(names_from = method, values_from = mspe) |>
  mutate(
    diff = abs(log10(lmer) - log10(mcpherson))  # difference on log scale
  ) |>
  mutate(label = if_else(rank(-diff) <= top_n, spp, NA_character_))

write_csv(
  plot_df |> select(spp, lmer, mcpherson, diff) |>
    mutate(lmer_rmspe = sqrt(lmer), mcpherson_rmspe = sqrt(mcpherson)),
  paste0(rslt_dir, 'dbh-ht-rmspe-comparison.csv')
)

ggplot(plot_df, aes(x = lmer, y = mcpherson)) +
  geom_point(alpha = 0.75, aes(color = factor(label > 0))) +
  geom_abline(slope = 1, linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = label), 
                           na.rm = TRUE, 
                           size = 3,
                           hjust = 1,
                           direction = 'y',
                           nudge_x = -0.05) +
  scale_y_log10() +
  scale_x_log10() +
  guides(color = guide_none()) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.1, color = "black")
  ) +
  labs(x = "mixed model RMSPE", y = "species model RMSPE")

paste(plt_dir, 'dbh-ht-mcpherson-lmer-rmspe-comparison.png', sep = '') |>
  ggsave(width = 3, height = 3, units = 'in', dpi = 300)

## species-level plots (back-transformed, one png per species) -----------------

spp_plt_dir <- '_figures/dbh-ht/'
# Clear the directory first: this loop writes one PNG per fitted species, and a
# run that fits fewer species than the last one would otherwise leave the old
# files behind, indistinguishable by name from current output.
fs::dir_delete(spp_plt_dir[fs::dir_exists(spp_plt_dir)])
fs::dir_create(spp_plt_dir)

sigma_val <- blues_tbl$estimate[blues_tbl$term == "sigma"]

for (sp in sort(unique(dbh_ht_data$spp))) {
  sp_data <- dbh_ht_data |> filter(spp == sp)
  n_obs   <- nrow(sp_data)
  sp_coef <- sp_blups_se |> filter(spp == sp)
  if (nrow(sp_coef) == 0) next
  
  x_seq   <- seq(min(sp_data$dbh), max(sp_data$dbh), length.out = 200)
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
  
  eq_label <- sprintf("height = %.3f \u00d7 dbh^%.3f",
                      exp(sp_coef$intercept), sp_coef$slope)
  
  p <- ggplot(sp_data, aes(x = dbh, y = height)) +
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
         x = "DBH (cm)", y = "Height (m)") +
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