# Climate structures global associations between human footprint and arthropod diversity

Code to reproduce Zito et al. (2026). Models are fitted with the
[`HubbellGLM`](https://github.com/alessandrozito/HubbellGLM) R package on the
Global Malaise Trap Program (GMTP) arthropod data.

## Setup

R 4.3.1.

```r
install.packages(c("tidyverse", "Matrix", "sf", "terra", "rnaturalearth",
                   "rnaturalearthdata", "geosphere", "patchwork", "ggpubr",
                   "lmtest", "sandwich", "multcomp", "scales", "kgc", "blockCV",
                   "foreach", "doParallel", "httr", "httr2", "remotes"))
remotes::install_github("alessandrozito/HubbellGLM")
remotes::install_github("alessandrozito/BNPvegan")   # rarefaction, Fig S6
```

Or install the bundled source: `install.packages("HubbellGLM_1.0.0.tar.gz", repos = NULL)`.

All paths are set in `preprocessing/paths.R` (repository at `~/HubbellGLM-paper`).

## Run

```sh
Rscript analysis/00_run_all.R        # models, tables, figures, from the data in data/
```

Any step can be run on its own; `ONLY=5,6` runs a subset. The cross-validation
(`09_cv_benchmark.R`) is by far the slowest step.

`preprocessing/00_run_all.R` rebuilds `data/` from scratch. It needs the raw GMTP
export from BOLD (`~/gmtp_filtered_bcdm.tsv.zip`) and downloads ~17 GB of covariates.

## Structure

```
preprocessing/   01-15  build the dataset, rasters and co-occurrence matrix
analysis/        01-11  fit models, make every table and figure
figures/                paper figures, named by figure number
data/                   analysis inputs (see below); the rest is generated
output/                 models and tables, generated
```

## Scripts and outputs

| Script | Paper output |
|---|---|
| `analysis/01_fit_models.R` | Tables S2, S4 |
| `analysis/02_table_S3_fixed_sigma.R` | Table S3 |
| `analysis/03_descriptive_figures.R` | Fig 1b, 1c (maps); Figs S2, S3, S4, S8, S10 |
| `analysis/04_accumulation_curves.R` | Fig 1c (curves) |
| `analysis/05_index_curves.R` | Fig 2; Figs S7, S11 |
| `analysis/06_scenario_maps.R` | Fig 3 |
| `analysis/07_ci_comparison.R` | Fig S1 |
| `analysis/08_cv_folds.R` | CV partitions |
| `analysis/09_cv_benchmark.R` | CV fits |
| `analysis/10_cv_figure.R` | Figs S5, S9 |
| `analysis/11_admissibility_curves.R` | Fig S6 |

Figures 1 and S4 are assembled by hand from the panels above. Table S1 has no
data. Tables are written to `output/models/`.

## Data

`data/` ships the inputs the analysis needs:

| File | Content |
|---|---|
| `GMTP_analysis_dataset.rds` | 15,711 collection events with covariates (15,176 up to 2024 are analysed) |
| `JaccardSim_sparse.rds` | shared-BIN similarity between events |
| `rasters_by_year/Raster_2024_full.rds` | global 0.1-degree grid for prediction |
| `FigS6_abundances.rds` | BIN counts for the two events in Fig S6 |

## Data sources

| Data | Source |
|---|---|
| GMTP arthropod records | [BOLD Systems](https://www.boldsystems.org) |
| Human footprint | Mu et al., *Scientific Data* (figshare) |
| AET, PET | TerraClimate |
| Weather, anomalies | ERA5-Land, Copernicus Climate Data Store |
| Wetlands | Copernicus Global Land Cover |
| Realms | WWF / TNC ecoregions, Harvard Dataverse |
| Köppen-Geiger zones | `kgc` R package |
| Landmass and isolation | Natural Earth |
