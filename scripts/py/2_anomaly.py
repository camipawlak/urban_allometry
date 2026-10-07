# -*- coding: utf-8 -*-
"""
scripts/py/2_anomaly.py

Usage (terminal):
  python scripts/py/2_anomaly.py

Detects outliers using Isolation Forest per species across multiple model comparisons:
  - dbh_height              : DBH vs Height
  - age_dbh_height          : Tree Age vs (DBH, Height)
  - dbh_height_canopy_width : (DBH, Height) vs Canopy Width (requires canopy data)
  - dbh_height_canopy_area  : (DBH, Height) vs Canopy Area (requires canopy data)

Output: Single CSV with all anomaly scores/flags as columns.
  - <model>_anom_score: Isolation Forest anomaly score (NaN if model skipped)
  - <model>_is_anom:    Boolean flag at threshold (NaN if model skipped)

Logging (terminal runs):
  - All print output is written to:
      <log_dir>/<script_stem>-log.txt
"""

# %% setup
import os
import argparse
from pathlib import Path
import contextlib
from datetime import datetime

import pandas as pd
import numpy as np
from isotree import IsolationForest


# %% Canopy merge utilities
def normalize_row_id(s: pd.Series) -> pd.Series:
    """Normalize row_id for robust merging: handle int/str/float mismatches."""
    s = s.astype(str).str.strip()
    s = s.str.replace(r"\.0$", "", regex=True)
    s = s.replace({"nan": np.nan, "None": np.nan, "": np.nan})
    return s


def canopy_merge_diagnostics(df_before: pd.DataFrame, df_after: pd.DataFrame) -> None:
    """Print diagnostics of canopy merge success/failures."""
    print("\n--- Canopy merge diagnostics ---")
    print("Input rows:", len(df_before), "Merged rows:", len(df_after))

    canopy_cols = [
        "canopy_width_ew_m", "canopy_width_ns_m", "canopy_width_avg_m",
        "canopy_width_ew_ft", "canopy_width_ns_ft", "canopy_width_avg_ft",
        "canopy_area_m2",
    ]

    present = [c for c in canopy_cols if c in df_after.columns]
    if present:
        any_canopy = df_after[present].notna().any(axis=1)
        print("Rows with ANY canopy value:", int(any_canopy.sum()), f"({any_canopy.mean():.1%})")

    print("--- end diagnostics ---\n")


# %% config (interactive defaults)
data_dir = Path("data/_processed")
out_dir = Path("data/_processed")
input_csv = "wca_filtered.csv"  # output of 1_quality_control.py
canopy_csv = "canopy_data.csv"  # in data/_raw

out_csv = "wca_filtered_annotated.csv"

log_dir = Path("data/processing-log")
report_dir = Path("data/_reports")

threshold = 0.85
min_points_per_species = 5
enable_plots = False
chdir_repo_root = True


# %% CLI overrides (runs ONLY in terminal execution)
if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Detect outliers using Isolation Forest per species (multi-model).")
    p.add_argument("--data-dir", type=Path, default=data_dir, help="Directory containing the input CSV.")
    p.add_argument("--out-dir", type=Path, default=out_dir, help="Directory to save annotated dataset.")
    p.add_argument("--input-csv", default=input_csv, help="Input CSV filename (in data-dir).")
    p.add_argument("--canopy-csv", default=canopy_csv, help="Canopy data CSV (in data/_raw).")
    p.add_argument("--out-csv", default=out_csv, help="Annotated output CSV filename (in out-dir).")
    p.add_argument("--log-dir", type=Path, default=log_dir, help="Directory for processing logs.")
    p.add_argument("--report-dir", type=Path, default=report_dir, help="Directory for reports.")
    p.add_argument("--threshold", type=float, default=threshold, help="Anomaly threshold (default: 0.85).")
    p.add_argument("--min-points", type=int, default=min_points_per_species, help="Min points per species (default: 5).")
    p.add_argument("--plots", action="store_true", help="Generate visualization plots per species.")
    p.add_argument("--no-chdir-repo-root", action="store_true")
    args = p.parse_args()

    data_dir = args.data_dir
    out_dir = args.out_dir
    input_csv = args.input_csv
    canopy_csv = args.canopy_csv
    out_csv = args.out_csv
    log_dir = args.log_dir
    report_dir = args.report_dir
    threshold = args.threshold
    min_points_per_species = args.min_points
    enable_plots = args.plots
    chdir_repo_root = not args.no_chdir_repo_root


# %% resolve paths / ensure dirs
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
out_path = out_dir / out_csv
canopy_path = Path("data/_raw") / canopy_csv

assert input_path.exists(), f"Missing input file: {input_path}"

# log file name
script_stem = Path(__file__).stem if "__file__" in globals() else "2_anomaly"
log_path = log_dir / f"{script_stem}-log.txt"

# reports
report_summary = report_dir / "outlier_summary_multimodel.csv"


# %% redirect prints to log file when run as a script
redirect_ctx = (
    contextlib.redirect_stdout(open(log_path, "w"))  # noqa: SIM115
    if __name__ == "__main__"
    else contextlib.nullcontext()
)

with redirect_ctx:
    def log_print(msg: str) -> None:
        print(msg)

    print("Run timestamp:", datetime.now().astimezone().isoformat(timespec="seconds"))
    print("-" * 80)

    print("Using data_dir  :", data_dir)
    print("Using out_dir   :", out_dir)
    print("Using log_dir   :", log_dir)
    print("Using report_dir:", report_dir)
    print("Input path      :", input_path)
    print("Output path     :", out_path)
    print("Canopy path     :", canopy_path)
    print("Log file        :", log_path)
    print("Threshold       :", threshold)
    print("Min points/spec :", min_points_per_species)

    # %% read data
    if __name__ == "__main__":
        import sys
        print("Reading input CSV...", file=sys.stderr)

    df = pd.read_csv(input_path)
    log_print(f"Loaded {len(df):,} rows")
    log_print(f"Columns: {list(df.columns)}")

    # --- PATCH: require genus column and carry it through to outputs ---
    assert "genus" in df.columns, "Missing required column 'genus' in input (expected from 1_quality_control output)."
    # --- END PATCH ---

    # %% try to merge canopy data
    canopy_merged = False
    df_raw = df.copy()
    
    # Attempt canopy merge
    log_print(f"\n{'='*80}")
    log_print("Step: Canopy merge")
    log_print(f"Canopy path: {canopy_path}")
    log_print(f"Canopy exists: {canopy_path.exists()}")
    log_print(f"{'='*80}")

    if canopy_path.exists():
        try:
            if __name__ == "__main__":
                import sys
                print("Merging canopy data...", file=sys.stderr)

            df_canopy = pd.read_csv(canopy_path)
            log_print(f"Loaded canopy file: {len(df_canopy)} rows, {len(df_canopy.columns)} columns")
            log_print(f"Canopy columns: {list(df_canopy.columns)}")

            if "row_id" not in df_canopy.columns:
                log_print("WARNING: Canopy file has no 'row_id' column. Canopy merge skipped.")
            elif "row_id" not in df.columns:
                log_print("WARNING: Main data has no 'row_id' column. Canopy merge skipped.")
            else:
                df["row_id"] = normalize_row_id(df["row_id"])
                df_canopy["row_id"] = normalize_row_id(df_canopy["row_id"])

                dup_ct = df_canopy["row_id"].duplicated().sum()
                if dup_ct > 0:
                    log_print(f"WARNING: canopy file has {dup_ct} duplicate row_id values; keeping first.")
                    df_canopy = df_canopy.drop_duplicates(subset=["row_id"], keep="first")

                canopy_cols = [c for c in df_canopy.columns if c != "row_id"]
                log_print(f"All canopy columns: {canopy_cols}")

                new_canopy_cols = [c for c in canopy_cols if c not in df.columns]
                log_print(f"New canopy columns to add: {new_canopy_cols}")

                if len(new_canopy_cols) == 0:
                    log_print("No new canopy columns to add (all already present).")
                else:
                    log_print(f"Merging {len(new_canopy_cols)} new canopy columns: {new_canopy_cols[:3]}...")
                    df_before_cols = set(df.columns)

                    df_merged = df.merge(df_canopy[["row_id"] + new_canopy_cols], on="row_id", how="left")

                    if all(c in df_merged.columns for c in df.columns):
                        df = df_merged
                        df_after_cols = set(df.columns)
                        added_cols = df_after_cols - df_before_cols
                        log_print(f"Merge successful. Added {len(added_cols)} columns: {list(added_cols)[:5]}...")
                        canopy_merged = True
                    else:
                        missing = [c for c in df.columns if c not in df_merged.columns]
                        log_print(f"ERROR: Merge lost columns: {missing}. Reverting merge.")
        except Exception as e:
            log_print(f"Canopy merge failed: {e}")
    else:
        log_print(f"Canopy file not found: {canopy_path}")

    print("Columns after merge:", list(df.columns))
    canopy_merge_diagnostics(df_raw, df)

    # %% Model configurations
    model_configs = {
        "dbh_height": {
            "feature_cols": ["dbh_cm_bill", "height_m_bill"],
            "labels": ["DBH (cm)", "Height (m)"],
        },
        "age_dbh_height": {
            "feature_cols": ["tree_age_yrs", "dbh_cm_bill", "height_m_bill"],
            "labels": ["Tree Age (years)", "DBH (cm)", "Height (m)"],
        },
        "dbh_height_canopy_width": {
            "feature_cols": ["dbh_cm_bill", "height_m_bill", "canopy_width_avg_m"],
            "labels": ["DBH (cm)", "Height (m)", "Canopy Width (avg, m)"],
        },
        "dbh_height_canopy_area": {
            "feature_cols": ["dbh_cm_bill", "height_m_bill", "canopy_area_m2"],
            "labels": ["DBH (cm)", "Height (m)", "Canopy Area (m²)"],
        },
    }

    # %% Initialize output dataset with all score/flag columns
    log_print(f"\n{'='*80}")
    log_print("Step: Anomaly detection")
    log_print(f"Models: {list(model_configs.keys())}")
    log_print(f"{'='*80}")

    df_ann = df.copy()
    for model_name in model_configs.keys():
        df_ann[f"{model_name}_anom_score"] = np.nan
        df_ann[f"{model_name}_is_anom"] = np.nan

    group_col = "scientific_name_clean"
    assert group_col in df_ann.columns, f"Required column missing: {group_col}"

    # %% Detect anomalies using Isolation Forest per model per species
    if __name__ == "__main__":
        import sys
        print("Detecting anomalies per species...", file=sys.stderr)

    models_fitted = {}

    for model_name, cfg in model_configs.items():
        feature_cols = cfg["feature_cols"]
        log_print(f"\nFitting model: {model_name} with features {feature_cols}")
        score_col = f"{model_name}_anom_score"
        flag_col = f"{model_name}_is_anom"

        # Check if required columns exist
        missing = [c for c in feature_cols if c not in df_ann.columns]
        if missing:
            print(f"Skipping {model_name}: missing {missing}")
            continue

        print(f"\nFitting model: {model_name} (features: {', '.join(feature_cols)})")
        fitted_species = 0
        scored_rows = 0

        for species, group in df_ann.groupby(group_col):
            group_valid = group.dropna(subset=feature_cols)
            if len(group_valid) < min_points_per_species:
                continue

            X = group_valid[feature_cols].to_numpy()
            ndim = len(feature_cols)

            iso = IsolationForest(
                ndim=ndim,
                ntrees=200,
                sample_size="auto",
                nthreads=1,
                random_seed=101825,
            )

            iso.fit(X)
            scores = iso.predict(X, output="score")  # higher = more anomalous
            df_ann.loc[group_valid.index, score_col] = scores

            fitted_species += 1
            scored_rows += len(group_valid)

        print(f"  Fitted species: {fitted_species:,}")
        print(f"  Scored rows   : {scored_rows:,} / {len(df_ann):,}")

        # Apply threshold flag
        df_ann[flag_col] = df_ann[score_col] > float(threshold)
        models_fitted[model_name] = (fitted_species, scored_rows)

    # %% Visualization: Individual plots per species (optional)
    if enable_plots:
        import matplotlib.pyplot as plt

        if __name__ == "__main__":
            import sys
            print("Creating visualizations...", file=sys.stderr)

        for model_name, cfg in model_configs.items():
            score_col = f"{model_name}_anom_score"
            flag_col = f"{model_name}_is_anom"
            feature_cols = cfg["feature_cols"]
            labels = cfg["labels"]
            # Skip if columns missing

            if score_col not in df_ann.columns or not all(c in df_ann.columns for c in feature_cols):
                print(f"Skipping visualization for {model_name}: required columns missing")
                continue
            
            # For 2D models: plot x vs y; for 3D+: plot first two features
            x_col = feature_cols[0]
            y_col = feature_cols[1]
            x_label = labels[0]
            y_label = labels[1]
            
            # Filter to rows with all features and a score
            df_plot = df_ann[feature_cols + [group_col, score_col, flag_col]].dropna(subset=feature_cols + [score_col])
            if len(df_plot) == 0:
                print(f"Skipping visualization for {model_name}: no valid data")
                continue
            # Create one plot per species

            for species in sorted(df_plot[group_col].unique()):
                df_spp = df_plot[df_plot[group_col] == species]
                if len(df_spp) < 2:
                    continue

                fig, ax = plt.subplots(figsize=(8, 6))
                
                # Normal points
                normal = df_spp[flag_col] == False
                ax.scatter(df_spp.loc[normal, x_col], df_spp.loc[normal, y_col], alpha=0.5, s=30, color="blue", label="Normal")

                anom = df_spp[flag_col] == True
                ax.scatter(df_spp.loc[anom, x_col], df_spp.loc[anom, y_col], alpha=0.8, s=50, color="red", marker="X", label="Anomaly")

                ax.set_xlabel(x_label)
                ax.set_ylabel(y_label)
                model_title = model_name.replace("_", " ").title()
                if len(feature_cols) > 2:
                    model_title += f" (showing {x_label} vs {y_label})"
                ax.set_title(f"{species} - {model_title}")
                ax.legend()
                ax.grid(True, alpha=0.3)
                
                # Save figure
                safe_species_name = species.replace(" ", "_").replace("/", "_")
                img_file = report_dir / f"{safe_species_name}_{model_name}.png"
                plt.savefig(img_file, dpi=100, bbox_inches="tight")
                plt.close()
                print(f"Saved: {img_file}")

    # %% Summary report
    summary_data = []
    for model_name, (n_species, n_rows) in models_fitted.items():
        flag_col = f"{model_name}_is_anom"
        n_flagged = (df_ann[flag_col] == True).sum()
        summary_data.append({
            "model": model_name,
            "species_fitted": n_species,
            "rows_scored": n_rows,
            "rows_flagged": n_flagged,
            "pct_flagged": n_flagged / len(df_ann) if len(df_ann) > 0 else 0,
        })

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(report_summary, index=False)
    print("\nModel Summary:")
    print(summary_df)
    print(f"\nWrote report: {report_summary}")

    # %% Save annotated dataset
    log_print(f"\n{'='*80}")
    log_print("Step: Save outputs")
    log_print(f"Output path: {out_path}")

    if __name__ == "__main__":
        import sys
        print("Saving output...", file=sys.stderr)
    df_ann.to_csv(out_path, index=False)
    log_print(f"[OK] Saved annotated dataset: {len(df_ann):,} rows, {len(df_ann.columns)} columns")
    log_print(f"{'='*80}")

# %% terminal hint (only prints in interactive mode)
if __name__ != "__main__":
    print(f"(Interactive) If you run this script in terminal, logs will be written to:\n  {log_path}")
    print(f"(Interactive) Reports will be written to:\n  {report_summary}")