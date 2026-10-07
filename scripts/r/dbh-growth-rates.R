library(tidyverse)
library(terra)
library(sf)


data_dir <- 'data/_processed/'
plt_dir <- '_img/dbh-rates/'
rslt_dir <- 'rslt/dbh-rates/'
fs::dir_create(plt_dir)
fs::dir_create(rslt_dir)



# Both scripts that read the inventory CSV take the same override, so the public
# release can be substituted for the private pipeline endpoint in one place:
#   INVENTORY_CSV <- "_zenodo/wca_inventory_public.csv"
if (!exists("INVENTORY_CSV"))
  INVENTORY_CSV <- paste0(data_dir, 'wca_filtered_annotated.csv')

# read in data. any_of() so the lat/lon rename is a no-op on the public release,
# which has no coordinates; the standardisation below then takes the lookup route.
wca <- INVENTORY_CSV |>
  read_csv() |>
  rename_with(~str_remove_all(.x, '_bill')) |>
  rename(any_of(c(dbh = "dbh_cm",
                  height = "height_m",
                  age = "tree_age_yrs",
                  spp = "scientific_name_clean",
                  date = "date_worked",
                  lat = "latitude",
                  lon = "longitude")))

wca_rates <- wca %>%
  group_by(spp, tree_id) |>
  filter(n() >= 2) |>
  mutate(
    delta_days  = as.numeric(difftime(date, lag(date), units = "days")),
    delta_yrs   = delta_days / 365.25,
    delta_dbh_cm = dbh - lag(dbh),
    dbh_start_cm = lag(dbh),
    dbh_end_cm = dbh
  ) |>
  ungroup() |>
  mutate(growth_rate = delta_dbh_cm/delta_yrs,
         growth_flag = ifelse(growth_rate > 4 | growth_rate == 0, T, F),
         delta_flag = ifelse(delta_yrs > 2, F, T),
         delta_dbh_flag = ifelse(delta_dbh_cm > 0, F, T)
  ) %>%
  filter(!is.na(delta_yrs),
         is.finite(delta_yrs), is.finite(delta_dbh_cm))

spp_of_interest <- unique(wca_rates$spp)

wca_rates_clean <- wca_rates %>% filter(growth_flag == F, delta_flag == F, delta_dbh_flag == F) |>
  filter(spp %in% spp_of_interest,
         dbh > 2.54) %>%
  group_by(spp) %>%
  summarize(median_growth_rate_cm_per_year = median(growth_rate))


# ── Frost-free days standardization ──────────────────────────────────────────
# Standardize growth rates to 153 frost-free days using per-species average
# frost-free days from the USDA Soil frost-free days raster.
#
# Coordinates are used only to reach a per-species mean: the standardization
# consumes one number per species and nothing tree-level survives. So there are
# two equivalent routes, and the public release (which has no coordinates) can
# take the second:
#
#   raster route  extract per tree from data/_raw/frost_free_days.tif, then
#                 average by species. Requires lat/lon.
#   lookup route  read the per-species means directly from
#                 _zenodo/species_frost_free_days.csv, produced by a previous
#                 raster run. Requires no coordinates.
#
# The raster route is used when coordinates are present, so existing behaviour is
# unchanged. NOTE: the raster is 8-bit and saturates at 255; much of urban
# California genuinely exceeds that (roughly 300-365 frost-free days), so the
# divisor is capped and standardized rates are conservative there.

FROST_RASTER <- "data/_raw/frost_free_days.tif"
FROST_LOOKUP <- "_zenodo/species_frost_free_days.csv"

have_coords <- all(c("lat", "lon") %in% names(wca)) &&
               any(is.finite(wca$lat) & is.finite(wca$lon))

if (have_coords && file.exists(FROST_RASTER)) {
  frost_rast <- rast(FROST_RASTER)

  # One location per unique tree with valid coordinates
  tree_pts <- wca %>%
    filter(!is.na(lat), !is.na(lon),
           lat != 0, lon != 0,
           lon < -100,               # drop any swapped coordinates
           is.finite(lat), is.finite(lon)) %>%
    distinct(spp, tree_id, .keep_all = TRUE) %>%
    select(spp, tree_id, lat, lon)

  pts_vect <- vect(tree_pts, geom = c("lon", "lat"), crs = "EPSG:4326")

  # Extract frost-free days at each tree location
  frost_vals <- terra::extract(frost_rast, pts_vect)
  tree_pts$frost_free_days <- frost_vals[, 2]

  # Average frost-free days per species across all valid tree locations
  spp_frost <- tree_pts %>%
    filter(!is.na(frost_free_days), frost_free_days > 0) %>%
    group_by(spp) %>%
    summarize(avg_frost_free_days = mean(frost_free_days))

  message("frost-free days: raster route (", nrow(spp_frost), " species)")

} else if (file.exists(FROST_LOOKUP)) {
  spp_frost <- readr::read_csv(FROST_LOOKUP, show_col_types = FALSE)
  stopifnot(all(c("spp", "avg_frost_free_days") %in% names(spp_frost)))
  message("frost-free days: lookup route, no coordinates needed (",
          nrow(spp_frost), " species)")

} else {
  stop("need either ", FROST_RASTER, " with coordinates, or ", FROST_LOOKUP)
}

# Join and compute standardized rate: rate × (153 / avg_frost_free_days)
wca_rates_clean <- wca_rates_clean %>%
  left_join(spp_frost, by = "spp") %>%
  mutate(standardized_growth_rate_cm_per_year = median_growth_rate_cm_per_year * (153 / avg_frost_free_days))

write.csv(wca_rates_clean, file = paste(rslt_dir, "wca_dbh_growth_rates.csv", sep = ""))

