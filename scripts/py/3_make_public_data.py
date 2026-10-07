# -*- coding: utf-8 -*-
"""
scripts/py/3_make_public_data.py

Usage (terminal):
  python scripts/py/3_make_public_data.py

Builds the public release of the inventory data from the endpoint of the private
preprocessing chain, for deposition alongside the paper.

This is a final, additive layer: it reads the output of 2_anomaly.py and writes a
new file. Nothing upstream is modified, so the private pipeline is unaffected.

WHAT IS REMOVED

  1. Location. latitude/longitude are dropped -- the reason this release exists.
  2. Columns no downstream analysis reads. Raw anomaly scores (the boolean flags
     are kept), the alternate canopy geometries (ew/ns, area, imperial units), and
     internal join keys (row_id, poly_id, poly_idx, index_pt, canopy_year).
  3. genus, which make-partitions.R recomputes from the species name anyway.

WHAT IS KEPT, AND WHY

  Every column the R analysis actually reads is retained, so this file can be
  substituted for wca_filtered_annotated.csv without changing any result. The set
  was determined by auditing column references across scripts/r/, not by
  judgement:

    tree_id                          grouping for the per-tree collapse and the
                                     tree-level train/test split
    scientific_name_clean            species; also the source of genus
    date_worked_bill                 required by the growth-rate quality screen
                                     in make-partitions.R (successive measurements)
    dbh_cm_bill, height_m_bill       the response/predictor pairs
    tree_age_yrs                     age-DBH and age-height models
    canopy_width_avg_m               DBH-canopy model
    *_is_anom (3 of them)            per-model anomaly filters

KNOWN LIMITATION

  scripts/r/dbh-growth-rates.R cannot run against this file. It joins trees to a
  frost-free-days raster via terra::vect(geom = c("lon","lat")), so it needs
  coordinates by construction. That analysis is not reproducible from a
  location-free release; everything else in the paper is.
"""

# %% setup
import argparse
from pathlib import Path

import pandas as pd

data_dir = Path("data/_processed")
out_dir = Path("_zenodo")
input_csv = "wca_filtered_annotated.csv"
out_csv = "wca_inventory_public.csv"

# columns retained, in the order they appear in the source
KEEP = [
    "tree_id",
    "scientific_name_clean",
    "date_worked_bill",
    "dbh_cm_bill",
    "height_m_bill",
    "tree_age_yrs",
    "canopy_width_avg_m",
    "dbh_height_is_anom",
    "age_dbh_height_is_anom",
    "dbh_height_canopy_width_is_anom",
]

DROP_LOCATION = ["latitude_bill", "longitude_bill"]


def build(in_path: Path, out_path: Path) -> pd.DataFrame:
    df = pd.read_csv(in_path, low_memory=False)

    missing = [c for c in KEEP if c not in df.columns]
    if missing:
        raise SystemExit(f"input is missing required columns: {missing}")

    # fail loudly if a location column is somehow still present after subsetting
    pub = df[KEEP].copy()
    leaked = [c for c in pub.columns if c in DROP_LOCATION or "lat" in c.lower()
              or "lon" in c.lower()]
    if leaked:
        raise SystemExit(f"location columns leaked into the release: {leaked}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    pub.to_csv(out_path, index=False)
    return pub


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Build the public, location-free inventory release."
    )
    p.add_argument("--data-dir", type=Path, default=data_dir,
                   help="Directory containing the input CSV.")
    p.add_argument("--out-dir", type=Path, default=out_dir,
                   help="Directory to write the public release into.")
    p.add_argument("--input-csv", default=input_csv,
                   help="Input CSV filename (in data-dir).")
    p.add_argument("--out-csv", default=out_csv,
                   help="Output CSV filename (in out-dir).")
    a = p.parse_args()

    src = a.data_dir / a.input_csv
    dst = a.out_dir / a.out_csv
    pub = build(src, dst)

    print(f"read   {src}")
    print(f"wrote  {dst}")
    print(f"  rows    {len(pub):,}")
    print(f"  trees   {pub['tree_id'].nunique():,}")
    print(f"  species {pub['scientific_name_clean'].nunique():,}")
    print(f"  columns {len(pub.columns)} (from {len(pd.read_csv(src, nrows=0).columns)})")
