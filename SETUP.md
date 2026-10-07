# Setup Instructions

Covers dependency installation and Python/R configuration for running the preprocessing and analysis scripts.

---

## Install dependencies

### Option A: pip (quickest)

```bash
python -m pip install -r requirements.txt
```

### Option B: conda (reproducible environment)

```bash
conda env create -f environment.yml
conda activate allometry
```

---

## R dependencies

The R analysis scripts are not covered by `requirements.txt` or `environment.yml`.
Install them from CRAN:

```r
install.packages(c(
  "cowplot", "fs", "ggforce", "ggnewscale", "ggrepel", "here", "lme4", "merTools", "modelr", "patchwork", "readxl", "rlang", "scales", "sf", "stringr", "terra", "tibble", "tidyr", "tidyverse", "writexl"
))
```

The exploratory notebook `notebooks/rmd/age_dbh_filters.Rmd` also uses:

```r
install.packages(c("ggeffects", "plotly"))
```

Versions used to produce the committed results:

| package | version |
|---|---|
| `cowplot` | 1.2.0 | 
| `fs` | 2.1.0 | 
| `ggforce` | 0.5.0 | 
| `ggnewscale` | 0.5.2 | 
| `ggrepel` | 0.9.8 | 
| `lme4` | 2.0.1 | 
| `merTools` | 1.0.0 | 
| `modelr` | 0.1.11 | 
| `patchwork` | 1.3.2 | 
| `readxl` | 1.5.0 | 
| `rlang` | 1.3.0 | 
| `scales` | 1.4.0 | 
| `sf` | 1.1.1 | 
| `stringr` | 1.6.0 | 
| `terra` | 1.9.34 | 
| `tibble` | 3.3.1 | 
| `tidyr` | 1.3.2 | 
| `tidyverse` | 2.0.0 | 
| `writexl` | 2.0.0 | 

R itself: R version 4.6.1 (2026-06-24)

Two of these are easy to miss because they are called with `::` rather than
attached, so a missing install surfaces only mid-run:

- `merTools` — `combine_data.R` uses `merTools::predictInterval()` for the
  tree-level prediction intervals. Without it the script aborts partway.
- `ggnewscale` — `pooled_plots.R` and `regional_plots_by_spp.R` need
  `new_scale_color()`. In `pooled_plots.R` the per-species plots are wrapped in
  `possibly()`, so a missing install reports `skipped N` rather than an error.

---

## Run scripts

From the repo root:

```bash
python scripts/py/0_clean_names.py
python scripts/py/1_quality_control.py
python scripts/py/2_anomaly.py
```

Alternatively, use the provided launchers (conda environment must be active or pip install complete):

```bash
# macOS/Linux
chmod +x run_script.sh
./run_script.sh 0_clean_names.py

# Windows
run_script.bat 0_clean_names.py
```

---

## VS Code integration

After installing dependencies, VS Code should find the Python interpreter automatically via `.vscode/settings.json`. If it doesn't:

1. Open Command Palette: **Ctrl+Shift+P** (or **Cmd+Shift+P**)
2. Select `Python: Select Interpreter`
3. Choose the `allometry` conda environment, or whichever environment has the installed packages

---

## Verify setup

```bash
python scripts/py/0_clean_names.py --help
```

If you see the script's help message with argument descriptions, you're all set.
