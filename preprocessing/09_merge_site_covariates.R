# Merge the time-invariant and annual covariates onto the sample and the grid:
# biogeographic realm, Koppen-Geiger zone, actual and potential
# evapotranspiration, human footprint.
#
# Produces: GMTP_sample_covariates.rds, cache/annual/gmtp_cov_YYYY.rds

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

library(tidyverse)
library(terra)
library(sf)
library(lubridate)
# NOT library(kgc): it attaches plyr, which masks dplyr's mutate(), count(),
# rename(), summarise() and arrange(). Everything kgc provides is reached with
# kgc:: instead.

gmtp_in <- file.path(DATA, "GMTP_sample_clean.rds")
grid_in <- GRID
wwf_shp <- file.path(CACHE, "ecoregions", "dataverse_files",
                     "tnc_terr_ecoregions.shp")
hfp_dir <- file.path(CACHE, "hfp")
clim_dir <- file.path(CACHE, "terraclimate")
data_dir <- DATA
ras_dir <- RASTER
tmp_dir <- file.path(CACHE, "hfp_unzipped")   # scratch, emptied every year

dir.create(ras_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

FORCE     <- nzchar(Sys.getenv("FORCE"))
YEARS_ENV <- Sys.getenv("YEARS")

TC_URL <- paste0("http://thredds.northwestknowledge.net:8080/thredds/",
                 "fileServer/TERRACLIMATE_ALL/data/TerraClimate_%s_%d.nc")

#---- Helpers -----------------------------------------------------------------

kopp_to_zone <- function(koppen_codes) {
  zone_names <- c("A" = "Tropical", "B" = "Dry", "C" = "Temperate",
                  "D" = "Continental", "E" = "Polar")
  ifelse(koppen_codes == "Climate Zone info missing" | is.na(koppen_codes),
         NA_character_, unname(zone_names[substr(koppen_codes, 1, 1)]))
}

# Koppen-Geiger zone by lookup on the kgc 0.5-degree grid.
add_kg_zone <- function(d) {
  e <- new.env(parent = emptyenv())
  utils::data("climatezones", package = "kgc", envir = e)
  cz <- get("climatezones", envir = e)
  d %>%
    dplyr::mutate(rndCoord.lon = kgc::RoundCoordinates(Longitude),
                  rndCoord.lat = kgc::RoundCoordinates(Latitude)) %>%
    dplyr::left_join(cz, by = c("rndCoord.lat" = "Lat", "rndCoord.lon" = "Lon")) %>%
    dplyr::mutate(kg_code = as.character(Cls), zone = kopp_to_zone(kg_code)) %>%
    dplyr::select(-Cls, -rndCoord.lon, -rndCoord.lat)
}

# Biogeographic realm from the WWF/TNC ecoregion polygons. Joined on the
# distinct coordinates only - the grid repeats many of them.
add_realm <- function(d, wwf) {
  pts <- d %>% dplyr::distinct(Longitude, Latitude)
  sfp <- st_as_sf(pts, coords = c("Longitude", "Latitude"), crs = 4326,
                  remove = FALSE)
  j <- st_join(sfp, wwf, join = st_intersects) %>%
    st_drop_geometry() %>%
    dplyr::select(Longitude, Latitude, WWF_REALM2) %>%
    dplyr::distinct(Longitude, Latitude, .keep_all = TRUE) %>%
    # "Indo-Malay" in the shapefile is "Indomalayan" in the source data.
    dplyr::mutate(realm = ifelse(WWF_REALM2 == "Indo-Malay", "Indomalayan",
                                 WWF_REALM2)) %>%
    dplyr::select(-WWF_REALM2)
  dplyr::left_join(d, j, by = c("Longitude", "Latitude"))
}

# Years whose NetCDF is not on the server, so the failed download is attempted
# once rather than once per dataset.
tc_missing <- new.env(parent = emptyenv())

# Annual TerraClimate variable at the given points: sum of the 12 monthly
# layers, x 0.1. aet and pet are the same product and the same operation, so
# they share one extractor.
#
# terra already applies the NetCDF scale factor, so the layers come back in mm
# (a July AET maximum of ~193 mm). The extra x 0.1 is the unit convention the
# paper uses - `aet` is therefore annual AET in units of 10 mm, matching
# dataset$aet exactly, and `pet` is on the same scale.
extract_tc <- function(v, year, pts_ll) {
  key <- paste0(v, year)
  if (!is.null(tc_missing[[key]])) return(rep(NA_real_, nrow(pts_ll)))
  nc <- file.path(clim_dir, sprintf("TerraClimate_%s_%d.nc", v, year))
  if (!file.exists(nc)) {
    message("     downloading ", basename(nc))
    ok <- tryCatch({
      download.file(sprintf(TC_URL, v, year), nc, mode = "wb", quiet = TRUE); TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    # A missing year comes back as a small error page, not a NetCDF file.
    if (!ok || !file.exists(nc) || file.size(nc) < 1e6) {
      unlink(nc); tc_missing[[key]] <- TRUE
      message("     ", toupper(v), " unavailable for ", year, " - filling NA")
      return(rep(NA_real_, nrow(pts_ll)))
    }
  }
  r <- tryCatch(rast(nc), error = function(e) NULL)
  if (is.null(r)) {
    message("     ", toupper(v), " raster unreadable for ", year, " - filling NA")
    return(rep(NA_real_, nrow(pts_ll)))
  }
  crs(r) <- "EPSG:4326"          # TerraClimate is lon/lat; the file says "unknown"
  pv <- vect(pts_ll, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
  e  <- terra::extract(r, pv, method = "bilinear")
  rowSums(e[, -1, drop = FALSE], na.rm = TRUE) * 0.1
}

#---- Load --------------------------------------------------------------------
message("Loading inputs ...")
gmtp <- readRDS(gmtp_in)
gmtp$year <- year(gmtp$collection_start_date)

grid <- readRDS(grid_in) %>%
  distinct(Longitude, Latitude, .keep_all = TRUE)

message("  GMTP sample: ", nrow(gmtp), " fieldids, ",
        nrow(distinct(gmtp, Longitude, Latitude)), " coordinates")
message("  grid:        ", nrow(grid), " coordinates")

years <- if (nzchar(YEARS_ENV)) {
  as.integer(strsplit(YEARS_ENV, ",")[[1]])
} else {
  sort(unique(gmtp$year))
}
message("  years:       ", paste(range(years), collapse = "-"),
        " (", length(years), " of them)")

#---- Time-invariant covariates -----------------------------------------------
message("\nRealm and Koppen-Geiger zone (time-invariant) ...")
wwf <- st_make_valid(st_read(wwf_shp, quiet = TRUE))

grid <- grid %>% add_kg_zone() %>% add_realm(wwf)
gmtp <- gmtp %>% add_kg_zone() %>% add_realm(wwf)

message("  grid: ", sum(is.na(grid$zone)), " points with no KG zone, ",
        sum(is.na(grid$realm)), " with no realm")
message("  GMTP: ", sum(is.na(gmtp$zone)), " fieldids with no KG zone, ",
        sum(is.na(gmtp$realm)), " with no realm")

#---- Project once for HFP ----------------------------------------------------
# The Mollweide projection of the points does not change between years, so it
# is done once and reused.
hfp_any <- list.files(hfp_dir, pattern = "^hfp\\d{4}\\.zip$", full.names = TRUE)
if (!length(hfp_any)) stop("no HFP archives in ", hfp_dir, call. = FALSE)
hfp_crs <- crs(rast(paste0("/vsizip/", hfp_any[1], "/",
                           sub("\\.zip$", ".tif", basename(hfp_any[1])))))

grid_v_moll <- project(vect(grid, geom = c("Longitude", "Latitude"),
                            crs = "EPSG:4326"), hfp_crs)
gmtp_v_moll <- project(vect(gmtp %>% distinct(Longitude, Latitude),
                            geom = c("Longitude", "Latitude"),
                            crs = "EPSG:4326"), hfp_crs)
gmtp_coords <- gmtp %>% distinct(Longitude, Latitude)

#---- Yearly loop -------------------------------------------------------------

for (yy in years) {
  f_ras  <- file.path(ras_dir, sprintf("Raster_%d.rds", yy))
  f_gmtp <- file.path(ras_dir, sprintf("gmtp_cov_%d.rds", yy))

  if (!FORCE && file.exists(f_ras) && file.exists(f_gmtp)) {
    message("\n", yy, ": already done - skipping")
    next
  }
  message("\n", yy, ":")
  t_year <- Sys.time()

  #-- AET and PET (same product, same extractor) -------------------------------
  message("   AET ...")
  aet_grid <- extract_tc("aet", yy, grid[, c("Longitude", "Latitude")])
  aet_gmtp <- extract_tc("aet", yy, gmtp_coords)
  message("   PET ...")
  pet_grid <- extract_tc("pet", yy, grid[, c("Longitude", "Latitude")])
  pet_gmtp <- extract_tc("pet", yy, gmtp_coords)

  #-- HFP ---------------------------------------------------------------------
  zipf <- file.path(hfp_dir, sprintf("hfp%d.zip", yy))
  if (!file.exists(zipf)) {
    message("   HFP unavailable for ", yy, " - filling NA")
    hfp_grid <- rep(NA_real_, nrow(grid))
    hfp_gmtp <- rep(NA_real_, nrow(gmtp_coords))
  } else {
    message("   HFP: unzipping ", basename(zipf), " ...")
    unlink(list.files(tmp_dir, full.names = TRUE))
    unzip(zipf, exdir = tmp_dir, junkpaths = TRUE)
    tif <- list.files(tmp_dir, pattern = "\\.tif$", full.names = TRUE)
    if (length(tif) != 1) {
      stop("expected one .tif in ", basename(zipf), ", found ", length(tif),
           call. = FALSE)
    }
    r <- rast(tif)
    message("   HFP: extracting ", nrow(grid), " grid + ",
            nrow(gmtp_coords), " GMTP points (bilinear) ...")
    hfp_grid <- terra::extract(r, grid_v_moll, method = "bilinear")[, 2]
    hfp_gmtp <- terra::extract(r, gmtp_v_moll, method = "bilinear")[, 2]
    rm(r)
    # Free the 2.4 GB before the next year is unzipped.
    unlink(list.files(tmp_dir, full.names = TRUE))
    message("   HFP: raster deleted")
  }

  #-- write -------------------------------------------------------------------
  saveRDS(grid %>% mutate(year = yy, aet = aet_grid, pet = pet_grid,
                          hfp = hfp_grid), f_ras, compress = "gzip")
  saveRDS(gmtp_coords %>% mutate(year = yy, aet = aet_gmtp, pet = pet_gmtp,
                                 hfp = hfp_gmtp), f_gmtp)

  message("   wrote ", basename(f_ras), " and ", basename(f_gmtp),
          "  [", round(difftime(Sys.time(), t_year, units = "mins"), 1), " min]")
}

#---- Attach the yearly values to the GMTP sample -----------------------------
message("\nAttaching aet/hfp to the GMTP sample ...")
# Every per-year slice on disk, not just the ones this run produced. Reading
# only `years` would blank the covariates for every other year whenever the
# script is run with a YEARS subset.
gmtp_cov <- list.files(ras_dir, pattern = "^gmtp_cov_\\d{4}\\.rds$",
                       full.names = TRUE) %>%
  map_dfr(readRDS)

gmtp_out <- gmtp %>%
  left_join(gmtp_cov, by = c("Longitude", "Latitude", "year"))

stopifnot(nrow(gmtp_out) == nrow(gmtp))

f_out <- file.path(data_dir, "GMTP_sample_covariates.rds")
saveRDS(gmtp_out, f_out)
write_tsv(gmtp_out, file.path(data_dir, "GMTP_sample_covariates.tsv"))

message("\n--- coverage --------------------------------------------------")
message("  fieldids:          ", nrow(gmtp_out))
message("  with realm:        ", sum(!is.na(gmtp_out$realm)))
message("  with KG zone:      ", sum(!is.na(gmtp_out$zone)))
message("  with aet:          ", sum(!is.na(gmtp_out$aet)))
message("  with hfp:          ", sum(!is.na(gmtp_out$hfp)))
message("  complete cases:    ",
        sum(complete.cases(gmtp_out[, c("realm", "zone", "aet", "hfp")])))
message("---------------------------------------------------------------")
message("Wrote ", f_out)
message("Yearly grids in ", ras_dir)

unlink(tmp_dir, recursive = TRUE)
