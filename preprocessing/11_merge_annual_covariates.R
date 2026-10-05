# Join the annual ERA5, wetland and PET extracts onto the sample and the grid.
#
# Produces: GMTP_sample_covariates_full.rds, rasters_by_year/Raster_YYYY_full.rds

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

library(tidyverse)

data_dir <- DATA
ras_dir <- RASTER

read_years <- function(pattern) {
  f <- list.files(ras_dir, pattern = pattern, full.names = TRUE)
  if (!length(f)) return(NULL)
  map_dfr(f, readRDS)
}

#---- 1. Sample ---------------------------------------------------------------

message("Sample ...")
samp <- readRDS(file.path(data_dir, "GMTP_sample_covariates.rds"))
e_s  <- read_years("^era5_sample_\\d{4}\\.rds$")
w_s  <- read_years("^wetlands_sample_\\d{4}\\.rds$")

if (is.null(e_s)) message("  no ERA5 files yet - skipping")
if (is.null(w_s)) message("  no wetland files yet - skipping")

full <- samp
if (!is.null(e_s)) {
  message("  ERA5: ", format(nrow(e_s), big.mark = ","), " coord-years, ",
          n_distinct(e_s$year), " years")
  full <- left_join(full, dplyr::select(e_s, -any_of("ID")),
                    by = c("Longitude", "Latitude", "year"))
}
if (!is.null(w_s)) {
  message("  wetlands: ", format(nrow(w_s), big.mark = ","), " coord-years, ",
          n_distinct(w_s$year), " years")
  full <- left_join(full, dplyr::select(w_s, -any_of("ID")),
                    by = c("Longitude", "Latitude", "year"))
}

new_cols <- setdiff(names(full), names(samp))
if (length(new_cols)) {
  message("  added ", length(new_cols), " columns; NA counts:")
  print(colSums(is.na(full[, new_cols, drop = FALSE])))
}

saveRDS(full, file.path(data_dir, "GMTP_sample_covariates_full.rds"))
write_tsv(full, file.path(data_dir, "GMTP_sample_covariates_full.tsv"))
message("  wrote GMTP_sample_covariates_full.rds  (", nrow(full), " rows)")

#---- 2. Grid, one year at a time --------------------------------------------
# 1.5 M rows x ~20 columns per year is fine on its own but not all at once, so
# each year is joined and written before the next is read.

message("\nGrid ...")
grid_files <- list.files(ras_dir, pattern = "^Raster_\\d{4}\\.rds$", full.names = TRUE)
wet_years <- as.integer(str_extract(
  list.files(ras_dir, pattern = "^wetlands_grid_\\d{4}\\.rds$"), "\\d{4}"))
if (length(wet_years)) message("  wetland map years available: ",
                               paste(sort(wet_years), collapse = ", "))
for (f in sort(grid_files)) {
  yr <- as.integer(str_extract(basename(f), "\\d{4}"))
  g  <- readRDS(f)
  fe <- file.path(ras_dir, sprintf("era5_grid_%d.rds", yr))
  add <- character(0)
  # The grid and the sample must carry the SAME classification, or the
  # coefficients are applied to labels they were not estimated under.
  g <- dplyr::mutate(g, zone_kgc = zone, kg_code_kgc = kg_code)
  if (file.exists(fe)) {
    g <- left_join(g, dplyr::select(readRDS(fe), -any_of("ID")),
                   by = c("Longitude", "Latitude", "year")); add <- c(add, "era5")
  }
  # Wetlands falls back to the nearest year that has a land-cover map, matching
  # the rule 14_build_analysis_dataset.R applies to the sample. CCI covers
  # 1992-2022, so years outside that carry the closest one; wetlands_src records
  # which year each cell actually came from, so the sample and the grid can
  # never silently disagree about it.
  if (length(wet_years)) {
    wy <- wet_years[which.min(abs(wet_years - yr))]
    fw <- file.path(ras_dir, sprintf("wetlands_grid_%d.rds", wy))
    # Keep script 07's own label rather than re-deriving it. Years outside the
    # CCI record (2023-2026) get a carry-forward FILE whose `year` equals the
    # target year, so `wy == yr` is true and re-deriving would call them
    # "observed" - which is how the 2024 grid ended up claiming an observed
    # wetland surface that is really the 2022 map.
    wsrc <- readRDS(fw)
    lab  <- if ("wetlands_src" %in% names(wsrc)) wsrc$wetlands_src[1] else "observed"
    g <- left_join(g, wsrc %>%
                     dplyr::select(-any_of(c("ID", "year", "wetlands_src"))),
                   by = c("Longitude", "Latitude")) %>%
      mutate(wetlands_src = if (wy == yr) lab else paste0("held_", wy))
    add <- c(add, if (wy == yr) "wetlands" else paste0("wetlands(", wy, ")"))
  }
  if (!length(add)) { message("  ", yr, ": nothing to add, skipping"); next }
  saveRDS(g, file.path(ras_dir, sprintf("Raster_%d_full.rds", yr)))
  message("  ", yr, ": +", paste(add, collapse = ", "), " -> Raster_", yr, "_full.rds")
  rm(g); gc(FALSE)
}

message("\nDone.")
