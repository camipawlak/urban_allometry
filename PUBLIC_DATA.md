# Inventory data — public release

Two files, deposited alongside the paper.

These are the cleaned analysis data behind the published equations. Tree
locations (latitude/longitude) and internal source-system columns have been
removed for data-privacy reasons and under the data-sharing agreement with West
Coast Arborists, Inc.; the raw inventory and all location-bearing intermediate
files are not publicly available. The code that produced these files from the
raw records is at https://github.com/camipawlak/urban_allometry (see `description.txt` for
provenance).

## `wca_inventory_public.csv`

One row per tree measurement. Trees measured more than once appear in multiple
rows, linked by `tree_id`.

| column | type | description |
|---|---|---|
| `tree_id` | integer | Identifier for an individual tree, stable across repeated measurements. |
| `scientific_name_clean` | text | Species name, after resolution against a reference workbook. |
| `date_worked_bill` | date | Date of measurement. |
| `dbh_cm_bill` | numeric | Diameter at breast height, centimetres. |
| `height_m_bill` | numeric | Tree height, metres. |
| `tree_age_yrs` | numeric | Estimated age in years. Present only for trees with an early-measurement record (DBH ≤ 1 inch). |
| `canopy_width_avg_m` | numeric | Mean crown width, metres. Present only where a canopy measurement was available. |
| `dbh_height_is_anom` | boolean | Outlier flag from per-species Isolation Forest on the DBH–height pair. |
| `age_dbh_height_is_anom` | boolean | Outlier flag, age against DBH and height. |
| `dbh_height_canopy_width_is_anom` | boolean | Outlier flag, canopy width against DBH and height. |

Column names carry a `_bill` suffix from the source system.

## `species_frost_free_days.csv`

One row per species.

| column | type | description |
|---|---|---|
| `spp` | text | Species name, matching `scientific_name_clean`. |
| `avg_frost_free_days` | numeric | Mean annual frost-free days across the sampled locations of that species. Values reach a ceiling of 255, the maximum the source raster records. |
