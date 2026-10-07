library(tidyverse)

data_dir <- 'data/_processed/'

# Hybrid names carry a leading "x" ("x Chitalpa tashkentensis"), so word(spp, 1)
# yields "x" rather than a genus. Resolve that here, at the single point where the
# genus column is created, so every downstream consumer sees a correct value.
fix_genus <- function(df) {
  df |>
    mutate(genus = case_when(
      genus == "x" & str_starts(spp, "x ") ~ word(spp, 2),
      genus == "x"                          ~ word(spp, 1),
      str_starts(genus, "x ")              ~ str_remove(genus, "^x "),
      TRUE                                 ~ genus
    ))
}

# Both scripts that read the inventory CSV take the same override, so the public
# release can be substituted for the private pipeline endpoint in one place:
#   INVENTORY_CSV <- "_zenodo/wca_inventory_public.csv"
if (!exists("INVENTORY_CSV"))
  INVENTORY_CSV <- paste0(data_dir, 'wca_filtered_annotated.csv')

# read in data. any_of() so the lat/lon rename is a no-op on the public release,
# which has no coordinates.
wca <- INVENTORY_CSV |>
  read_csv() |>
  rename_with(~str_remove_all(.x, '_bill')) |>
  rename(any_of(c(dbh = "dbh_cm",
                  height = "height_m",
                  age = "tree_age_yrs",
                  spp = "scientific_name_clean",
                  date = "date_worked",
                  lat = "latitude",
                  lon = "longitude"))) |>
  mutate(genus = word(spp, 1)) |>
  fix_genus()

# check one record per tree?
head(wca)
wca |> count(tree_id)

# retrieve spp with n >= 100
spp_of_interest <- wca |>
  filter(dbh > 2.54) |>
  distinct(spp, tree_id) |>
  count(spp) |>
  filter(n >= 100) |>
  pull(spp)

# filter to species and dbh of interest; retain all records per tree so that
# dbh-growth-rates.R can compute growth rates from repeated measurements
wca_clean <- wca |>
  filter(spp %in% spp_of_interest,
         dbh > 2.54)

## growth-rate quality screen --------------------------------------------------

# Successive measurements of the same tree imply a DBH growth rate. Implausible
# rates indicate a bad record (mis-measurement, mis-linked tree id, bad date), so
# trees failing the screen are flagged here rather than in each model script, and
# the models filter on `growth_qc_pass`. Computing it once keeps every model
# working from the same definition of a clean tree.
#
# A tree passes only if all successive intervals are plausible:
#   growth_flag    DBH growth > 4 cm/yr, or exactly zero
#   delta_flag     gap between measurements > 2 years
#   delta_dbh_flag DBH did not increase
#
# Trees with fewer than two dated measurements cannot be screened and do not pass.
wca_rates <- wca_clean |>
  group_by(spp, tree_id) |>
  filter(n() >= 2, !is.na(age)) |>
  arrange(date, .by_group = TRUE) |>
  mutate(
    delta_days   = as.numeric(difftime(date, lag(date), units = "days")),
    delta_yrs    = delta_days / 365.25,
    delta_dbh_cm = dbh - lag(dbh)
  ) |>
  ungroup() |>
  mutate(
    growth_rate    = delta_dbh_cm / delta_yrs,
    growth_flag    = ifelse(growth_rate > 4 | growth_rate == 0, TRUE, FALSE),
    delta_flag     = ifelse(delta_yrs > 2, FALSE, TRUE),
    delta_dbh_flag = ifelse(delta_dbh_cm > 0, FALSE, TRUE)
  ) |>
  filter(!is.na(delta_yrs), is.finite(delta_yrs), is.finite(delta_dbh_cm))

clean_tree_ids <- wca_rates |>
  filter(growth_flag == FALSE, delta_flag == FALSE, delta_dbh_flag == FALSE) |>
  distinct(tree_id)

wca_clean <- wca_clean |>
  mutate(growth_qc_pass = tree_id %in% clean_tree_ids$tree_id)

message(sprintf("growth QC: %s of %s trees pass",
                format(n_distinct(clean_tree_ids$tree_id), big.mark = ","),
                format(n_distinct(wca_clean$tree_id), big.mark = ",")))

# check fractions and counts of anomalous points by species
wca_clean |>
  group_by(spp) |>
  summarize(across(contains('is_anom'), .fns = list('mean' = mean, 'n' = sum))) |>
  view()

## data partitioning -----------------------------------------------------------

# partition by tree_id (80/20 within species); all records for each tree stay together
set.seed(22626)
pct_train <- 0.8

train_tree_ids <- wca_clean |>
  distinct(spp, tree_id) |>
  group_by(spp) |>
  slice_sample(prop = pct_train) |>
  ungroup() |>
  pull(tree_id)

wca_train <- wca_clean |> filter( tree_id %in% train_tree_ids)
wca_val   <- wca_clean |> filter(!tree_id %in% train_tree_ids)

# export
save(list = c('wca_train', 'wca_val'), file = 'data/_processed/wca-partitions.RData')
