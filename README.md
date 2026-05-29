# Code and Data for Zito et al. (2026) - Predicting global biodiversity via Hubbell regression

This repository contains all R code and derived data necessary to reproduce the 
analyses and figures in Zito et al. (2026) - Predicting global biodiversity via Hubbell regression. 
The study applies the Hubbell Regression (`HubbellGLM` R package) to the Global 
Malaise Trap Program (GMTP) arthropod dataset to quantify how energy availability 
(actual evapotranspiration, AET) and 
human pressure (Human Footprint Index, HFP) drive arthropod species richness across 
global climatic zones. 
---

## Repository structure

```
HubbellGLM-paper/
├── preprocessing/
│   ├── 0_GMTP_preprocess_merge.R   # Merges occurrence + climate data (requires raw GMTP files)
│   ├── 1_GMTP_create_dataset.R     # Builds the analysis dataset and Jaccard similarity matrix
│   └── 2_GMTP_create_raster.R      # Assembles the global 2014 raster for out-of-sample prediction
├── analysis/
│   ├── 1_run_models.R              # Fits stepwise models M0–M7; reproduces Tables S2–S4 and Figure S5
│   ├── 2_plot_prediction_results.R # Produces Figures 1b, 1c, 1e, 1f, 4, S1–S4
│   ├── 3_plot_Index_variations.R   # Produces Figure 3 (diversity indices by zone and HFP)
│   └── 4_Simulation_Figure2.R      # 10-fold cross-validation; reproduces Figures 2 and S6
├── data/
│   ├── data_GMTP_clean.rdata                # R workspace: `dataset` + `JaccardSim` matrix
│   ├── data_GMTP_merged_covariates.rds.gzip # Compressed RDS: `data_fields` covariate table
│   ├── abundance_matrix.rds.gzip           # Compressed RDS: species abundance matrix
│   ├── AET_vs_others_dropDeviance.txt       # Pre-computed deviance table (legacy)
│   ├── TableS2_AET_vs_others_dropDeviance.txt  # Pre-computed deviance table (Table S2; read by scripts)
│   ├── CI_indexes_poly.rds                  # Pre-computed bootstrap CIs for Figure 3
│   ├── climate/
│   │   ├── Raster_2014.rds                  # Global raster (before AET merge)
│   │   ├── Raster_2014_merged.rds           # Global raster (with AET; used for Figures 1, 4)
│   │   ├── TerraClimate_aet_2014.nc         # TerraClimate AET NetCDF (downloaded by preprocessing/2)
│   │   └── ecoregions/dataverse_files/      # WWF TNC terrestrial ecoregions shapefile
│   └── simulation_output/
│       ├── folds_random.rdata               # Pre-computed random 10-fold CV splits
│       └── Simulation_10fold_models.tsv     # Pre-computed CV results (Figure 2)
├── figures/                                 # All output figures (PDF and PNG)
├── HubbellGLM_1.0.0.tar.gz                 # Local copy of the HubbellGLM package source
├── session_info.txt                         # R session info (versions of all packages used)
└── HubbellGLM-paper.Rproj                  # RStudio project file
```

---

## Data sources

| Dataset | Source | Access |
|---|---|---|
| GMTP arthropod occurrences | Global Malaise Trap Program (BOLD Systems) | [https://www.boldsystems.org](https://www.boldsystems.org) |
| ERA5 weekly climate | Copernicus Climate Data Store | [https://cds.climate.copernicus.eu](https://cds.climate.copernicus.eu) |
| TerraClimate (AET, PET) | IDAHO EPSCOR / Northwest Knowledge Network | [https://www.climatologylab.org/terraclimate.html](https://www.climatologylab.org/terraclimate.html) |
| GIMMS NDVI | NASA Global Inventory Monitoring and Modeling System | [https://glam1.gsfc.nasa.gov](https://glam1.gsfc.nasa.gov) |
| Human Footprint Index (HFP) | Mu et al. (2022) | [https://www.nature.com/articles/s41597-022-01284-8](https://www.nature.com/articles/s41597-022-01284-8) |
| Wetland cover | Global Surface Water / Copernicus | [https://global-surface-water.appspot.com](https://global-surface-water.appspot.com) |
| WWF TNC Ecoregions | The Nature Conservancy (Harvard Dataverse) | Included in `data/climate/ecoregions/` |
| Koppen-Geiger zones | `kgc` R package (built-in `climatezones` data) | CRAN |

The **analysis-ready dataset** (`data/data_GMTP_clean.rdata`) is provided directly in this repository so that all analysis scripts can be run without access to the raw GMTP occurrence files.

---

## R package dependencies

All scripts were developed in **R ≥ 4.3**. Install all required packages with:

```r
# Core data and plotting
install.packages(c(
  "tidyverse", "lubridate", "splines", "MASS",
  "ggpubr", "patchwork", "scales", "RColorBrewer", "ggsignif"
))

# Spatial and geographic
install.packages(c(
  "sf", "rnaturalearth", "rnaturalearthdata",
  "geosphere", "spdep", "kgc"
))

# Ecology and statistics
install.packages(c(
  "vegan", "lmtest", "multcomp",
  "lmerTest", "stargazer", "conleyreg", "slider", "sandwich"
))

# Parallel computing and cross-validation
install.packages(c("foreach", "doParallel", "caret"))

# BNP vegan - a function to plot acculation curves
remotes::install_github("alessandrozito/BNPvegan")

# HubbellGLM (the core modelling package for this paper)
# Install from the companion repository:
remotes::install_github("alessandrozito/HubbellGLM")

```


---

## How to reproduce the analysis

Scripts must be run in order. All paths are relative to `~/HubbellGLM-paper/` (the RStudio project root).

### Step 0 — Preprocessing (optional if using provided data)

These steps require the raw GMTP occurrence files and climate CSV downloads, which are not included in this repository. **Skip to Step 1 if using the provided `data_GMTP_clean.rdata`.**

```r
# Merge GMTP occurrences with ERA5, TerraClimate, NDVI, HFP, and wetland data
source("preprocessing/0_GMTP_preprocess_merge.R")

# Build the analysis dataset and compute the Jaccard similarity matrix
source("preprocessing/1_GMTP_create_dataset.R")

# Assemble the global 2014 prediction raster
# (automatically downloads TerraClimate_aet_2014.nc if not present)
source("preprocessing/2_GMTP_create_raster.R")
```

### Step 1 — Model fitting (Tables S2–S4, Figure S5)

```r
source("analysis/1_run_models.R")
```

Fits the stepwise model sequence M0–M7 with both the polynomial and canonical links. Reproduces the coefficient tables (Tables S3 and S4) and deviance comparison across environmental predictors (Table S2). Also produces Figure S5 (confidence interval ranges under four sandwich estimators). The `rerun_table` flag inside the script controls whether Table S2 is recomputed from scratch (takes ~30 min); pre-computed results are loaded from `data/TableS2_AET_vs_others_dropDeviance.txt`.

### Step 2 — Figures 1, 4, S1–S4

```r
source("analysis/2_plot_prediction_results.R")
```

Produces all exploratory and model-output figures. Requires `data_GMTP_clean.rdata` and `data/climate/Raster_2014_merged.rds`.

### Step 3 — Figure 3 (diversity index variations)

```r
source("analysis/3_plot_Index_variations.R")
```

Plots how AET and HFP jointly predict six diversity indices across climatic zones. Bootstrap CIs (500 simulations) are cached in `data/CI_indexes_poly.rds`; re-run from scratch by deleting that file.

### Step 4 — Figures 2 and S6 (10-fold cross-validation)

```r
source("analysis/4_Simulation_Figure2.R")
```

Compares HubbellGLM against Poisson, Negative Binomial, Fisher's alpha, and linear richness models via random 10-fold cross-validation. Pre-computed results are loaded by default; set `rerun <- TRUE` to recompute (runs in parallel on 10 cores, takes ~2–4 h). Also produces Figure S6, which overlays fitted accumulation curves from Poisson, Negative Binomial, and HubbellGLM on the empirical rarefaction curve for a selected site.

---

## Figures produced

| Script | Figure | Description |
|---|---|---|
| `2_plot_prediction_results.R` | Figure 1b | Site locations coloured by climatic zone |
| `2_plot_prediction_results.R` | Figure 1c | AET and HFP in the Amazon |
| `2_plot_prediction_results.R` | Figure 1e | Predicted accumulation curves for Amazon sites |
| `2_plot_prediction_results.R` | Figure 1f | Global maps of diversity indices (S_σ, α, Shannon, Simpson, Hill) |
| `4_Simulation_Figure2.R` | Figure 2 | 10-fold CV comparison of regression models |
| `3_plot_Index_variations.R` | Figure 3 | Diversity index predictions by AET, zone, and HFP |
| `2_plot_prediction_results.R` | Figure 4a | Predicted richness change if HFP = 0 |
| `2_plot_prediction_results.R` | Figure 4b | Predicted richness change under +10% HFP |
| `2_plot_prediction_results.R` | Figure 4c | Predicted richness change under +10% AET |
| `2_plot_prediction_results.R` | Figure S1 | Accumulation curve shapes across σ and η |
| `2_plot_prediction_results.R` | Figure S2 | Geographic distance vs. Jaccard similarity |
| `2_plot_prediction_results.R` | Figure S3 | Global AET, HFP, and climatic zone rasters |
| `2_plot_prediction_results.R` | Figure S4 | Study design and diversity–covariate relationships (panels a–d) |
| `1_run_models.R` | Figure S5 | Coefficient CIs under Normal, Quasi, White, and Jaccard sandwich estimators |
| `4_Simulation_Figure2.R` | Figure S6 | Empirical rarefaction vs. fitted accumulation curves (Poisson, NB, HubbellGLM) |

---

