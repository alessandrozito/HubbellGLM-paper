# Shared paths and constants. Sourced by every script.

REPO   <- path.expand("~/HubbellGLM-paper")

DATA   <- file.path(REPO, "data")                   # inputs and data products
OUT    <- file.path(REPO, "output")                 # analysis outputs
FIG    <- file.path(REPO, "figures")                # paper figures
EXTRA  <- file.path(FIG, "extra")                   # side figures, not in the paper

QC     <- file.path(DATA, "qc")                     # QC tables
RASTER <- file.path(DATA, "rasters_by_year")        # annual grids, Raster_YYYY_full.rds
CACHE  <- file.path(DATA, "cache")                  # downloaded external inputs
MODELS <- file.path(OUT, "models")                  # fitted models, predictions, tables
SIM    <- file.path(OUT, "simulation")              # cross-validation output

for (d in c(DATA, OUT, FIG, EXTRA, QC, RASTER, CACHE, MODELS, SIM))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Climate classification used throughout.
ZONE <- "zone_kgc"
ZLV  <- c("Continental", "Dry", "Polar", "Temperate", "Tropical")

# Analysis window and reference year for the maps.
YEAR_MAX <- as.integer(Sys.getenv("YEAR_MAX", "2024"))
REF_YEAR <- as.integer(Sys.getenv("REF_YEAR", "2024"))

# Global 0.1-degree prediction lattice (coordinates and elevation).
GRID <- file.path(DATA, "grid_coordinates.rds")
