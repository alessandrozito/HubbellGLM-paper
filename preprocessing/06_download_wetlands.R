# Annual wetland cover for the sample and the global grid.
#
# Produces: cache/annual/wetlands_*_YYYY.rds

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

library(tidyverse)
library(terra)
library(httr2)

gmtp_in <- file.path(DATA, "GMTP_sample_covariates.rds")
grid_in <- GRID
data_dir <- DATA
ras_dir <- RASTER
tmp_dir <- file.path(CACHE, "lc_tmp")

dir.create(ras_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

FORCE     <- nzchar(Sys.getenv("FORCE"))
YEARS_ENV <- Sys.getenv("YEARS")

CDS_BASE   <- "https://cds.climate.copernicus.eu/api"
COLLECTION <- "satellite-land-cover"
WET_CLASSES <- c(160, 170, 180)
FACT        <- 36                    # 1/360 deg -> 0.1 deg

terraOptions(memfrac = 0.6, progress = 0)

# CCI is split across product versions by year. A year outside every range has
# no map at all and is filled by carry-forward at the end.
cci_version <- function(y) {
  if (y >= 1992 && y <= 2015) "v2_0_7cds" else
  if (y >= 2016 && y <= 2022) "v2_1_1"   else NA_character_
}

cds_key <- function() {
  k <- Sys.getenv("CDSAPI_KEY")
  if (nzchar(k)) return(k)
  rc <- path.expand("~/.cdsapirc")
  if (!file.exists(rc)) {
    stop("No CDS credentials. See the header of 04_download_era5_annual.R.", call. = FALSE)
  }
  ln <- readLines(rc, warn = FALSE)
  k  <- sub("^\\s*key\\s*:\\s*", "", grep("^\\s*key\\s*:", ln, value = TRUE))[1]
  if (is.na(k) || !nzchar(k)) stop("~/.cdsapirc has no 'key:' line", call. = FALSE)
  str_squish(k)
}

cds_download <- function(year, dest, key) {
  ver <- cci_version(year)
  if (is.na(ver)) stop("no CCI land cover for ", year, call. = FALSE)
  body <- list(variable = "all", year = list(as.character(year)),
               version = list(ver), download_format = "zip")

  resp <- request(paste0(CDS_BASE, "/retrieve/v1/processes/", COLLECTION, "/execute")) %>%
    req_headers(`PRIVATE-TOKEN` = key, Accept = "application/json") %>%
    req_body_json(list(inputs = body)) %>%
    req_error(is_error = function(r) FALSE) %>%
    req_perform()
  if (resp_status(resp) >= 400) {
    stop("CDS rejected the ", year, " request (HTTP ", resp_status(resp), "): ",
         resp_body_string(resp), call. = FALSE)
  }
  job <- resp_body_json(resp)$jobID
  message("     job ", job, " (", ver, ")")

  wait <- 10
  repeat {
    st <- request(paste0(CDS_BASE, "/retrieve/v1/jobs/", job)) %>%
      req_headers(`PRIVATE-TOKEN` = key, Accept = "application/json") %>%
      req_perform() %>% resp_body_json()
    if (identical(st$status, "successful")) break
    if (st$status %in% c("failed", "dismissed", "rejected")) {
      # "rejected" carries no message field and generally means the account's
      # concurrency cap was hit, so it is worth retrying rather than aborting.
      stop("CDS job for ", year, " ", st$status,
           if (identical(st$status, "rejected"))
             " (concurrency cap? retry when script 06 has finished)" else "",
           ": ", paste(unlist(st$message), collapse = " "), call. = FALSE)
    }
    Sys.sleep(wait); wait <- min(wait * 1.5, 60)
  }

  res <- request(paste0(CDS_BASE, "/retrieve/v1/jobs/", job, "/results")) %>%
    req_headers(`PRIVATE-TOKEN` = key, Accept = "application/json") %>%
    req_perform() %>% resp_body_json()
  href <- res$asset$value$href
  if (is.null(href)) stop("CDS returned no download href for ", year, call. = FALSE)
  message("     downloading ", round(as.numeric(res$asset$value$`file:size`) / 1e6), " MB")
  request(href) %>% req_progress() %>% req_perform(path = dest)
  invisible(dest)
}

#---- Point sets --------------------------------------------------------------

message("Loading the two point sets ...")
samp_pts <- readRDS(gmtp_in) %>% distinct(Longitude, Latitude)
grid_pts <- readRDS(grid_in) %>% distinct(Longitude, Latitude)
message("  sample: ", format(nrow(samp_pts), big.mark = ","), " coordinates")
message("  grid:   ", format(nrow(grid_pts), big.mark = ","), " coordinates")
samp_xy <- as.matrix(samp_pts[, c("Longitude", "Latitude")])
grid_xy <- as.matrix(grid_pts[, c("Longitude", "Latitude")])

years <- if (nzchar(YEARS_ENV)) as.integer(str_split_1(YEARS_ENV, ",")) else {
  yr <- readRDS(gmtp_in)$year
  seq(min(yr, na.rm = TRUE), max(yr, na.rm = TRUE))
}
message("Years: ", paste(range(years), collapse = "-"), " (", length(years), ")")

key <- cds_key()

#---- Main loop ---------------------------------------------------------------

done <- integer(0)
for (yr in years) {
  f_s <- file.path(ras_dir, sprintf("wetlands_sample_%d.rds", yr))
  f_g <- file.path(ras_dir, sprintf("wetlands_grid_%d.rds",   yr))
  if (!FORCE && file.exists(f_s) && file.exists(f_g)) {
    message("\n", yr, ": already done, skipping"); done <- c(done, yr); next
  }
  if (is.na(cci_version(yr))) {
    message("\n", yr, ": outside the CCI record - will carry forward"); next
  }
  message("\n", yr, ":")

  zipf <- file.path(tmp_dir, sprintf("cci_%d.zip", yr))
  ok <- tryCatch({ if (!file.exists(zipf)) cds_download(yr, zipf, key); TRUE },
                 error = function(e) { message("     ", conditionMessage(e)); FALSE })
  if (!ok) next

  ncf <- unzip(zipf, list = TRUE)$Name
  ncf <- ncf[grepl("\\.nc$", ncf)][1]
  unzip(zipf, files = ncf, exdir = tmp_dir, overwrite = TRUE)
  ncp <- file.path(tmp_dir, ncf)

  lc <- rast(ncp, subds = "lccs_class")
  message("     ", paste(dim(lc)[1:2], collapse = " x "), " at ",
          signif(res(lc)[1], 3), " deg")
  t_step <- Sys.time()

  # Binary flooded mask, written to disk so the aggregation streams rather than
  # trying to hold 8.4e9 cells in memory.
  mask_f <- file.path(tmp_dir, sprintf("wet_mask_%d.tif", yr))
  wet <- classify(lc, cbind(0:220, as.integer(0:220 %in% WET_CLASSES)),
                  others = 0, filename = mask_f, overwrite = TRUE,
                  datatype = "INT1U", wopt = list(gdal = "COMPRESS=LZW"))

  message("     reclassified in ",
          round(difftime(Sys.time(), t_step, units = "secs")), " s")
  t_step <- Sys.time()
  message("     aggregating by ", FACT, " to 0.1 degree ...")
  agg_f <- file.path(tmp_dir, sprintf("wet_pct_%d.tif", yr))
  wpct <- aggregate(wet, fact = FACT, fun = "mean", na.rm = TRUE,
                    filename = agg_f, overwrite = TRUE,
                    wopt = list(gdal = "COMPRESS=LZW")) * 100
  names(wpct) <- "wetlands"
  message("     aggregated in ",
          round(difftime(Sys.time(), t_step, units = "mins"), 1), " min")

  # Cell centres coincide with the grid points, so no interpolation is wanted.
  e_s <- terra::extract(wpct, samp_xy, method = "simple")
  e_g <- terra::extract(wpct, grid_xy, method = "simple")

  saveRDS(bind_cols(samp_pts, as_tibble(e_s)) %>%
            mutate(year = yr, wetlands_src = "observed"), f_s)
  saveRDS(bind_cols(grid_pts, as_tibble(e_g)) %>%
            mutate(year = yr, wetlands_src = "observed"), f_g)
  message("     wrote ", basename(f_s), " and ", basename(f_g))
  done <- c(done, yr)

  rm(lc, wet, wpct); gc(FALSE)
  unlink(c(zipf, ncp, mask_f, agg_f))
}

#---- Carry forward years with no map ----------------------------------------
# Years outside 1992-2022 take the nearest year that does have one, labelled so
# the distinction survives into the model frame.

gap <- setdiff(years, done)
if (length(gap) && length(done)) {
  message("\nCarrying forward to ", length(gap), " year(s) with no CCI map: ",
          paste(gap, collapse = ", "))
  for (yr in gap) {
    src_yr <- done[which.min(abs(done - yr))]
    for (tag in c("sample", "grid")) {
      s <- readRDS(file.path(ras_dir, sprintf("wetlands_%s_%d.rds", tag, src_yr)))
      s$year <- yr
      s$wetlands_src <- paste0("carried_from_", src_yr)
      saveRDS(s, file.path(ras_dir, sprintf("wetlands_%s_%d.rds", tag, yr)))
    }
    message("  ", yr, " <- ", src_yr)
  }
}

message("\nDone. Per-year files are in ", ras_dir)
message("Merge them with 11_merge_annual_covariates.R")
