# -*- coding: utf-8 -*-
"""
scripts/py/0_clean_names.py

Usage (terminal):
  python scripts/py/0_clean_names.py

Usage (terminal, custom paths / filenames):
  python scripts/py/0_clean_names.py \
    --data-dir data/_raw \
    --out-dir data/_processed \
    --input-csv wca_All_SpeciesBotanical_accessed012225.csv \
    --resolver-xlsx TreeNameResolver.xlsx \
    --top200-csv cufi_top_200_species_in_california.csv \
    --out-csv wca_name_resolved.csv \
    --log-dir data/_processing-log \
    --report-dir data/_reports

Usage (terminal, do NOT auto-chdir to repo root):
  python scripts/py/0_clean_names.py --no-chdir-repo-root

Interactive (VS Code):
  Run cells top-to-bottom (Shift+Enter on each # %% cell).
  Edit the "config" cell variables (data_dir/out_dir/filenames) as needed.
  The CLI block only runs when executing the file in the terminal.

Logging (terminal runs):
  - All print output is written to:
      <log_dir>/<script_stem>-log.txt

Reports (terminal runs + interactive if you run the cells):
  - Unmatched name counts (actionable backlog) written to:
      <report_dir>/unmatched_speciesbotanical_bycount.csv
      <report_dir>/unmatched_currentspeciesbotanical_bycount.csv
"""

# %% setup
print("Importing packages...")
import os
import argparse
from pathlib import Path
import contextlib
from datetime import datetime
import pandas as pd


# %% config (interactive defaults)
data_dir = Path("data/_raw")
out_dir = Path("data/_processed")
input_csv = "wca_All_SpeciesBotanical_accessed012225.csv"
resolver_xlsx = "TreeNameResolver.xlsx"
top200_csv = "cufi_top_200_species_in_california.csv"

out_csv = "wca_name_resolved.csv"
chdir_repo_root = True
log_dir = Path("data/processing-log")
report_dir = Path("data/_reports")


# %% CLI overrides (runs ONLY in terminal execution)
if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Clean WCA botanical names and subset to top 200 species.")
    p.add_argument("--data-dir", type=Path, default=data_dir)
    p.add_argument("--out-dir", type=Path, default=out_dir)
    p.add_argument("--input-csv", default=input_csv)
    p.add_argument("--resolver-xlsx", default=resolver_xlsx)
    p.add_argument("--top200-csv", default=top200_csv)
    p.add_argument("--out-csv", default=out_csv)
    p.add_argument("--log-dir", type=Path, default=log_dir, help="Directory for processing logs (default: data/processing-log).")
    p.add_argument("--report-dir", type=Path, default=report_dir, help="Directory for reports (default: data/_reports).")
    p.add_argument("--no-chdir-repo-root", action="store_true")
    args = p.parse_args()

    data_dir = args.data_dir
    out_dir = args.out_dir
    input_csv = args.input_csv
    resolver_xlsx = args.resolver_xlsx
    top200_csv = args.top200_csv
    out_csv = args.out_csv
    log_dir = args.log_dir
    report_dir = args.report_dir
    chdir_repo_root = not args.no_chdir_repo_root


# %% reset directory to repository top level (optional)
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
resolver_path = data_dir / resolver_xlsx
top200_path = data_dir / top200_csv
out_path = out_dir / out_csv

assert input_path.exists(), f"Missing input file: {input_path}"
assert resolver_path.exists(), f"Missing resolver file: {resolver_path}"
assert top200_path.exists(), f"Missing top200 file: {top200_path}"

# log file name: <script_stem>-log.txt (only reliable in terminal/script execution)
if "__file__" in globals():
    script_stem = Path(__file__).stem
else:
    script_stem = "0_clean_names"
log_path = log_dir / f"{script_stem}-log.txt"

# report outputs (unmatched counts)
report_unmatched_species = report_dir / "unmatched_speciesbotanical_bycount.csv"
report_unmatched_current = report_dir / "unmatched_currentspeciesbotanical_bycount.csv"

# %% redirect prints to log file when run as a script
redirect_ctx = (
    contextlib.redirect_stdout(open(log_path, "w"))  # noqa: SIM115
    if __name__ == "__main__"
    else contextlib.nullcontext()
)

if __name__ == "__main__":
    import sys
    print("Reading input CSV...", file=sys.stderr)

with redirect_ctx:
    print("Run timestamp:", datetime.now().astimezone().isoformat(timespec="seconds"))
    print("-" * 80)

    print("Using data_dir :", data_dir)
    print("Using out_dir  :", out_dir)
    print("Using log_dir  :", log_dir)
    print("Using report_dir:", report_dir)
    print("Input path     :", input_path)
    print("Resolver path  :", resolver_path)
    print("Top200 path    :", top200_path)
    print("Output path    :", out_path)
    print("Log file       :", log_path)
    print("Report (unmatched SpeciesBotanical) :", report_unmatched_species)
    print("Report (unmatched CurrentSpeciesBotanical) :", report_unmatched_current)

    # %% read in raw data file
    print("\nReading input CSV...")
    df = pd.read_csv(input_path).reset_index(drop=True)
    df["_row_id"] = df.index
    print("Loaded", len(df), "rows")

    # %% load resolver
    names_df = pd.read_excel(resolver_path)
    names_df["Tokenized_name"] = names_df["Name_submitted"].astype(str).str.strip()  # keep case-sensitive
    print("Resolver rows:", len(names_df))

    # %% TNRS helpers
    def run_tnrs(df_in: pd.DataFrame):
        df_in = df_in.copy()
        df_in["Tokenized_name"] = df_in["SCIENTIFIC_ORIG"].astype(str).str.strip()  # keep case-sensitive

        merged_df = (
            df_in.merge(names_df, on="Tokenized_name", how="left")
                 .drop(columns=["Name_submitted", "Tokenized_name"])
        )

        matched_df = merged_df[
            merged_df["Name_corrected"].notna() &
            (merged_df["Name_corrected"] != "")
        ]
        unmatched_df = merged_df.drop(index=matched_df.index)
        return matched_df, unmatched_df


    def build_correction_map(base_df: pd.DataFrame, source_col: str, corrected_col: str):
        temp = base_df.rename(columns={source_col: "SCIENTIFIC_ORIG"})
        matched, unmatched = run_tnrs(temp)

        mapping = (
            matched[["_row_id", "Name_corrected", "Genus_corrected"]]
            .drop_duplicates("_row_id", keep="first")
            .rename(columns={
                "Name_corrected": corrected_col,
                "Genus_corrected": "genus"
            })
        )
        return mapping, matched, unmatched

    def unmatched_counts(unmatched_df: pd.DataFrame, label: str, head_n: int = 60) -> pd.DataFrame:
        counts = (
            unmatched_df["SCIENTIFIC_ORIG"]
            .astype(str).str.strip()
            .value_counts()
            .reset_index()
            .rename(columns={"index": label, "SCIENTIFIC_ORIG": "count"})
        )
        print(counts.head(head_n))
        print(f"\nTotal unmatched unique ({label}): {len(counts):,}")
        return counts

    # %% Apply TNRS to BOTH columns (independent passes)
    print("\nCleaning botanical names...")
    if __name__ == "__main__":
        import sys
        print("Cleaning botanical names...", file=sys.stderr)
    base = df.copy()
    assert base["_row_id"].is_unique, "_row_id must be unique in base df"

    species_map, matched_species, unmatched_species = build_correction_map(
        base, "SpeciesBotanical", "SpeciesBotanical_corrected"
    )
    current_map, matched_current, unmatched_current = build_correction_map(
        base, "CurrentSpeciesBotanical", "CurrentSpeciesBotanical_corrected"
    )

    current_map = current_map.drop(columns=["genus"], errors="ignore")

    df_cleaned = (
        base
        .merge(species_map, on="_row_id", how="left", validate="one_to_one")  # keeps genus
        .merge(current_map, on="_row_id", how="left", validate="one_to_one")  # no genus collision
    )

    # %% combine outputs
    df_cleaned = (
        base
        .merge(species_map, on="_row_id", how="left", validate="one_to_one")
        .merge(current_map, on="_row_id", how="left", validate="one_to_one")
    )

    df_cleaned["SpeciesBotanical_clean"] = (
        df_cleaned["SpeciesBotanical_corrected"].fillna(df_cleaned["SpeciesBotanical"])
    )
    df_cleaned["CurrentSpeciesBotanical_clean"] = (
        df_cleaned["CurrentSpeciesBotanical_corrected"].fillna(df_cleaned["CurrentSpeciesBotanical"])
    )

    df_cleaned = df_cleaned.drop(columns=[
        "SpeciesBotanical_corrected",
        "CurrentSpeciesBotanical_corrected",
    ])

    # %% checks
    assert len(df_cleaned) == len(base), "Row count changed — merge problem!"
    assert df_cleaned["_row_id"].is_unique, "_row_id is not unique!"

    print("SpeciesBotanical cleaned:",
          (df_cleaned["SpeciesBotanical_clean"] != df_cleaned["SpeciesBotanical"]).sum())
    print("CurrentSpeciesBotanical cleaned:",
          (df_cleaned["CurrentSpeciesBotanical_clean"] != df_cleaned["CurrentSpeciesBotanical"]).sum())

    print("Unique SpeciesBotanical_clean:", df_cleaned["SpeciesBotanical_clean"].nunique())
    print("Unique CurrentSpeciesBotanical_clean:", df_cleaned["CurrentSpeciesBotanical_clean"].nunique())

    # %% unmatched reports (console/log + write to reports)
    unmatched_counts_species = unmatched_counts(unmatched_species, "unmatched_speciesbotanical")
    unmatched_counts_current = unmatched_counts(unmatched_current, "unmatched_currentspeciesbotanical")

    unmatched_counts_species.to_csv(report_unmatched_species, index=False)
    unmatched_counts_current.to_csv(report_unmatched_current, index=False)
    print("\nWrote unmatched reports:")
    print(" -", report_unmatched_species)
    print(" -", report_unmatched_current)

    # %% read in top 200 species list
    top_200 = pd.read_csv(top200_path)
    valid_species = set(top_200["species"].astype(str).str.strip())
    print("Number of valid species:", len(valid_species))  # should be ~200

    # %% subset to valid species
    n_rows_before = len(df_cleaned)
    always_keep = {"Lagerstroemia indica"}

    species_clean = df_cleaned["SpeciesBotanical_clean"].astype(str).str.strip()
    df_subset = df_cleaned[
        species_clean.isin(valid_species) | species_clean.isin(always_keep)
    ].copy()

    # %% summarize + export
    print("Rows before subsetting:", n_rows_before)
    print("Rows after subsetting :", len(df_subset))
    print("Unique spp after      :", df_subset["SpeciesBotanical_clean"].nunique())

    df_subset.to_csv(out_path, index=False)
    print("Writing output:", len(df_subset), "rows saved to", out_path)
    if __name__ == "__main__":
        import sys
        print("Done. Output saved to", out_path, file=sys.stderr)

# %% terminal hint (only prints in interactive mode)
if __name__ != "__main__":
    print(f"(Interactive) If you run this script in terminal, logs will be written to:\n  {log_path}")
    print(f"(Interactive) Unmatched reports will be written to:\n  {report_unmatched_species}\n  {report_unmatched_current}")