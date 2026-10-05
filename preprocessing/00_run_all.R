# End-to-end rebuild of the analysis dataset, the annual prediction rasters and
# the co-occurrence matrix.
#
# Steps marked [net] download from an external service and are slow; they are
# skipped when their output already exists. Run a single step directly instead
# of the whole driver when iterating.
#
#   Rscript preprocessing/00_run_all.R
#   ONLY=12,14  Rscript preprocessing/00_run_all.R    # a subset

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

STEPS <- c(
  "01_preprocess_sample.R",        # raw BCDM  -> one clean record per fieldid
  "02_build_occurrence.R",         #           -> occurrence + Jaccard matrices
  "03_download_hfp.R",             # [net] human footprint rasters
  "04_download_era5_annual.R",     # [net] annual ERA5-Land
  "05_fix_era5_land_mask.R",       #       repair land-masked sample cells
  "06_download_wetlands.R",        # [net] annual wetland cover
  "07_download_era5_daily.R",      # [net] daily ERA5-Land, 1991-2024
  "08_site_island.R",              #       landmass area and isolation, sites
  "09_merge_site_covariates.R",    # [net] realm, zone, aet, pet, hfp
  "10_fill_zone_gaps.R",           #       fill the Koppen zone where kgc is NA
  "11_merge_annual_covariates.R",  #       -> Raster_YYYY_full.rds
  "12_grid_island.R",              #       landmass area and isolation, grid
  "13_compute_anomalies.R",        #       collection-window weather anomalies
  "14_build_analysis_dataset.R",   #       -> GMTP_analysis_dataset.rds
  "15_grid_mess.R"                 #       extrapolation score into every raster
)

only <- Sys.getenv("ONLY", "")
if (nzchar(only)) {
  pick  <- sprintf("%02d", as.integer(strsplit(only, ",")[[1]]))
  STEPS <- STEPS[substr(STEPS, 1, 2) %in% pick]
}

here <- file.path(REPO, "preprocessing")
for (s in STEPS) {
  message("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70))
  t0 <- Sys.time()
  system2("Rscript", file.path(here, s))
  message("-> ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}
message("\ndone")
