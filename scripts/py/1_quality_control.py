# -*- coding: utf-8 -*-
"""
scripts/py/1_quality_control.py

Usage (terminal):
  python scripts/py/1_quality_control.py

Usage (terminal, custom paths / filenames):
  python scripts/py/1_quality_control.py \
    --data-dir data/_processed \
    --out-dir data/_processed \
    --input-csv wca_name_resolved.csv \
    --clean-out wca_filtered.csv \
    --removed-out wca_removals.csv \
    --log-dir data/_processing-log \
    --report-dir data/_reports

Interactive (VS Code):
  Run cells top-to-bottom (Shift+Enter on each # %% cell).
  Edit the "config" cell variables (data_dir/out_dir/filenames) as needed.
  The CLI block only runs when executing the file in the terminal.

Logging (terminal runs):
  - All print output is written to:
      <log_dir>/<script_stem>-log.txt

Reports:
  - Removal counts by reason written to:
      <report_dir>/qc_removed_counts.csv
"""

# %% setup
import os
import argparse
from pathlib import Path
import contextlib
from datetime import datetime
import pandas as pd
import numpy as np

# %% config (interactive defaults)
data_dir = Path("data/_processed")
out_dir = Path("data/_processed")
input_csv = "wca_name_resolved.csv"
clean_out = "wca_filtered.csv"
removed_out = "wca_removals.csv"
log_dir = Path("data/processing-log")
report_dir = Path("data/_reports")
chdir_repo_root = True

# %% CLI overrides (runs ONLY in terminal execution)
if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Clean and process tree history data for quality control.")
    p.add_argument("--data-dir", type=Path, default=data_dir, help="Base directory containing input files.")
    p.add_argument("--out-dir", type=Path, default=out_dir, help="Directory for outputs.")
    p.add_argument("--input-csv", default=input_csv, help="Input CSV filename (in data-dir).")
    p.add_argument("--clean-out", default=clean_out, help="Clean output CSV filename (in out-dir).")
    p.add_argument("--removed-out", default=removed_out, help="Removed output CSV filename (in out-dir).")
    p.add_argument("--log-dir", type=Path, default=log_dir, help="Directory for processing logs (default: data/processing-log).")
    p.add_argument("--report-dir", type=Path, default=report_dir, help="Directory for reports (default: data/_reports).")
    p.add_argument("--no-chdir-repo-root", action="store_true", help="Do not auto-chdir to repo root (skip .git search).")
    args = p.parse_args()

    data_dir = args.data_dir
    out_dir = args.out_dir
    input_csv = args.input_csv
    clean_out = args.clean_out
    removed_out = args.removed_out
    log_dir = args.log_dir
    report_dir = args.report_dir
    chdir_repo_root = not args.no_chdir_repo_root

# %% for interactive session: reset directory to repository top level (optional)
if chdir_repo_root:
    # repo root = first folder up the tree with .git or allometry.Rproj (a zip or
    # Zenodo download has no .git); fall back to the current folder
    repo_root = next((p for p in (Path.cwd(), *Path.cwd().parents)
                      if (p / ".git").exists() or (p / "allometry.Rproj").exists()),
                     Path.cwd())
    os.chdir(repo_root)

data_dir = data_dir.expanduser().resolve()
out_dir = out_dir.expanduser().resolve()
log_dir = log_dir.expanduser().resolve()
report_dir = report_dir.expanduser().resolve()

out_dir.mkdir(parents=True, exist_ok=True)
log_dir.mkdir(parents=True, exist_ok=True)
report_dir.mkdir(parents=True, exist_ok=True)

input_path = data_dir / input_csv
assert input_path.exists(), f"Missing input file: {input_path}"
clean_path = out_dir / clean_out
removed_path = out_dir / removed_out

# reports
report_removed_counts = report_dir / "qc_removed_counts.csv"

# log file name
if "__file__" in globals():
    script_stem = Path(__file__).stem
else:
    script_stem = "1_quality_control"
log_path = log_dir / f"{script_stem}-log.txt"

# %% redirect prints to log file when run as a script
redirect_ctx = (
    contextlib.redirect_stdout(open(log_path, "w"))  # noqa: SIM115
    if __name__ == "__main__"
    else contextlib.nullcontext()
)

with redirect_ctx:
    print("Run timestamp:", datetime.now().astimezone().isoformat(timespec="seconds"))
    print("-" * 80)

    print("Using data_dir  :", data_dir)
    print("Using out_dir   :", out_dir)
    print("Using log_dir   :", log_dir)
    print("Using report_dir:", report_dir)
    print("Input file      :", input_path)
    print("Clean output    :", clean_path)
    print("Removed output  :", removed_path)
    print("Log file        :", log_path)
    print("Report (qc_removed_counts):", report_removed_counts)

    # %% read input
    if __name__ == "__main__":
        import sys
        print("Reading input CSV...", file=sys.stderr)
    usecols = [
        "Id", "TreeId",
        "SpeciesBotanical_clean", "CurrentSpeciesBotanical_clean", "genus",
        "TreeHistoryBillDbh", "TreeHistoryBillHeight", "TreeHistoryDateWorked",
        "TreeHistoryLatitude", "TreeHistoryLongitude",
        "CurrentTreeLatitude", "CurrentTreeLongitude",
        "status",
    ]
    df_raw = pd.read_csv(input_path, usecols=lambda c: c in usecols)
    print("Read rows:", len(df_raw))
    print("Columns:", list(df_raw.columns))

    removed_records = []
    n_start = len(df_raw)

    # %% 1) unit conversion + date parse
    if __name__ == "__main__":
        import sys
        print("Converting units and parsing dates...", file=sys.stderr)
    df1 = df_raw.copy()
    df1["TreeHistoryBillDbh_cm"] = df1["TreeHistoryBillDbh"] * 2.54
    df1["TreeHistoryBillHeight_m"] = df1["TreeHistoryBillHeight"] * 0.3048
    df1["TreeHistoryDateWorked"] = pd.to_datetime(df1["TreeHistoryDateWorked"], errors="coerce")
    print("After conversion rows:", len(df1))

    # %% 2) remove missing/invalid values
    mask_invalid = (
        df1["TreeHistoryBillDbh"].isna()
        | df1["TreeHistoryBillHeight"].isna()
        | df1["TreeHistoryDateWorked"].isna()
        | (df1["TreeHistoryBillDbh"] <= 0)
        | (df1["TreeHistoryBillHeight"] <= 0)
    )
    invalid_rows = df1[mask_invalid].copy()
    invalid_rows["reason_removed"] = "missing_or_invalid_value"
    removed_records.append(invalid_rows)
    df2 = df1[~mask_invalid].copy()
    print("Removed invalid:", len(invalid_rows))
    print("Remaining:", len(df2))

    # %% 2b) remove under 5 ft
    mask_under_5 = df2["TreeHistoryBillHeight"] < 5
    under_5_rows = df2[mask_under_5].copy()
    under_5_rows["reason_removed"] = "under_5_ft"
    removed_records.append(under_5_rows)
    df3 = df2[~mask_under_5].copy()
    print("Removed under 5 ft:", len(under_5_rows))
    print("Remaining:", len(df3))

    # %% 2c) remove conflicting DBH on same date (per TreeId)
    dup_groups = (
        df3.groupby(["TreeId", "TreeHistoryDateWorked"])["TreeHistoryBillDbh"]
           .nunique()
    )
    conflict_keys = dup_groups[dup_groups > 1].index
    mask_conflict = (
        df3.set_index(["TreeId", "TreeHistoryDateWorked"]).index.isin(conflict_keys)
    )
    conflict_rows = df3[mask_conflict].copy()
    conflict_rows["reason_removed"] = "conflicting_dbh_same_date"
    removed_records.append(conflict_rows)
    df4 = df3[~mask_conflict].copy()
    print("Removed conflicting DBH groups rows:", len(conflict_rows))
    print("Remaining:", len(df4))

    # %% 3a) remove species changes over time (within TreeId)
    species_nunique = (
        df4.groupby("TreeId")[["SpeciesBotanical_clean", "CurrentSpeciesBotanical_clean"]]
           .nunique()
    )
    consistent_species_ids = species_nunique.query(
        "SpeciesBotanical_clean == 1 and CurrentSpeciesBotanical_clean == 1"
    ).index
    mask_inconsistent_time = ~df4["TreeId"].isin(consistent_species_ids)
    inconsistent_time_rows = df4[mask_inconsistent_time].copy()
    inconsistent_time_rows["reason_removed"] = "species_change_over_time"
    removed_records.append(inconsistent_time_rows)
    df5 = df4[~mask_inconsistent_time].copy()
    print("Removed inconsistent over time:", len(inconsistent_time_rows))
    print("Remaining:", len(df5))

    # %% 3b) remove species mismatch within a record
    mask_species_mismatch = df5["SpeciesBotanical_clean"] != df5["CurrentSpeciesBotanical_clean"]
    mismatch_rows = df5[mask_species_mismatch].copy()
    mismatch_rows["reason_removed"] = "species_changed_within_row"
    removed_records.append(mismatch_rows)
    df6 = df5[~mask_species_mismatch].copy()
    print("Removed mismatch within row:", len(mismatch_rows))
    print("Remaining:", len(df6))

    # %% 4) remove negative growth (per TreeId time series)
    df6_sorted = df6.sort_values(["TreeId", "TreeHistoryDateWorked"]).copy()

    dbh_diff = df6_sorted.groupby("TreeId")["TreeHistoryBillDbh"].diff()
    ht_diff  = df6_sorted.groupby("TreeId")["TreeHistoryBillHeight"].diff()

    bad_step = (dbh_diff < 0) | (ht_diff < -3)
    bad_tree_ids = df6_sorted.loc[bad_step, "TreeId"].unique()
    mask_bad_tree = df6_sorted["TreeId"].isin(bad_tree_ids)

    negative_growth_rows = df6_sorted[mask_bad_tree].copy()
    negative_growth_rows["reason_removed"] = "negative_growth"
    removed_records.append(negative_growth_rows)

    df7 = df6_sorted[~mask_bad_tree].copy()

    print("Removed negative growth rows:", len(negative_growth_rows))
    print("Remaining:", len(df7))
    print("Trees removed (negative growth):", len(bad_tree_ids))

    # %% 5) compute tree age
    df7 = df7.copy()
    df7["date_planted"] = np.where(
        df7["TreeHistoryBillDbh"] <= 1,
        df7["TreeHistoryDateWorked"],
        pd.NA
    )
    df7["date_planted"] = pd.to_datetime(df7["date_planted"], errors="coerce")

    dat_ages = (
        df7.dropna(subset=["date_planted"])
           .loc[:, ["TreeId", "date_planted"]]
           .groupby("TreeId", as_index=False)
           .agg(date_planted_age=("date_planted", "min"))
           .drop_duplicates(subset=["TreeId"])
    )

    df8 = df7.merge(dat_ages, on="TreeId", how="left")
    df8["tree_age_days"] = (df8["TreeHistoryDateWorked"] - df8["date_planted_age"]).dt.days
    df8["tree_age_yrs"] = df8["tree_age_days"] / 365
    print("Rows with computed tree_age_yrs:", df8["tree_age_yrs"].notna().sum())

    # %% 5b) lat/long preference
    df8 = df8.copy()
    df8["latitude_bill"] = df8.get("TreeHistoryLatitude")
    if "CurrentTreeLatitude" in df8.columns:
        df8["latitude_bill"] = df8["latitude_bill"].fillna(df8["CurrentTreeLatitude"])

    df8["longitude_bill"] = df8.get("TreeHistoryLongitude")
    if "CurrentTreeLongitude" in df8.columns:
        df8["longitude_bill"] = df8["longitude_bill"].fillna(df8["CurrentTreeLongitude"])

    # %% 6) keep + rename columns
    keep_cols = [
        "Id", "TreeId", "SpeciesBotanical_clean", "genus", "TreeHistoryDateWorked",
        "TreeHistoryBillDbh_cm", "TreeHistoryBillHeight_m",
        "latitude_bill", "longitude_bill", "status",
        "tree_age_days", "tree_age_yrs",
    ]
    rename_map = {
        "Id": "row_id",
        "TreeId": "tree_id",
        "SpeciesBotanical_clean": "scientific_name_clean",
        "genus": "genus",
        "TreeHistoryDateWorked": "date_worked_bill",
        "TreeHistoryBillDbh_cm": "dbh_cm_bill",
        "TreeHistoryBillHeight_m": "height_m_bill",
        "latitude_bill": "latitude_bill",
        "longitude_bill": "longitude_bill",
        "status": "status",
        "tree_age_days": "tree_age_days",
        "tree_age_yrs": "tree_age_yrs",
    }
    df_clean = df8.loc[:, [c for c in keep_cols if c in df8.columns]].rename(columns=rename_map)
    print("Final clean rows:", len(df_clean))
    print("Final clean columns:", list(df_clean.columns))

    # %% 7) combine removed + summary
    df_removed = (
        pd.concat(removed_records, ignore_index=True)
        if removed_records
        else pd.DataFrame()
    )
    print(f"Total records before cleaning: {n_start:,}")
    print(f"Total clean records: {len(df_clean):,}")
    print(f"Total removed records: {len(df_removed):,}")

    if not df_removed.empty and "reason_removed" in df_removed.columns:
        print("\nRemoved by reason:")
        print(df_removed["reason_removed"].value_counts())

    # %% export primary datasets
    df_clean.to_csv(clean_path, index=False)
    df_removed.to_csv(removed_path, index=False)
    if __name__ == "__main__":
        import sys
        print("Done. Output saved.", file=sys.stderr)
    print(f"\nSaved cleaned data to: {clean_path}")
    print(f"Saved removed records to: {removed_path}")

    # %% report: removal counts by reason
    if not df_removed.empty and "reason_removed" in df_removed.columns:
        removed_counts = (
            df_removed["reason_removed"]
            .value_counts()
            .rename_axis("reason_removed")
            .reset_index(name="n_removed")
        )
        removed_counts["pct_removed"] = removed_counts["n_removed"] / max(len(df_raw), 1)
        removed_counts.to_csv(report_removed_counts, index=False)
        print("\nWrote report:", report_removed_counts)

# %% terminal hint (only prints in interactive mode)
if __name__ != "__main__":
    print(f"(Interactive) If you run this script in terminal, logs will be written to:\n  {log_path}")
    print(f"(Interactive) QC removal counts report will be written to:\n  {report_removed_counts}")