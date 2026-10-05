# Annual ERA5-Land covariates for the sample and the global prediction grid.
#
# Produces: cache/annual/era5_{sample,grid}_YYYY.rds

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

library(tidyverse)
library(terra)
library(httr2)

gmtp_in  <- file.path(DATA, "GMTP_sample_covariates.rds")
grid_in <- GRID
data_dir <- DATA
ras_dir <- RASTER
tmp_dir  <- file.path(CACHE, "era5_tmp")        # scratch, emptied every year

dir.create(ras_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

FORCE     <- nzchar(Sys.getenv("FORCE"))
YEARS_ENV <- Sys.getenv("YEARS")

CDS_BASE   <- "https://cds.climate.copernicus.eu/api"
COLLECTION <- "reanalysis-era5-land-monthly-means"

# CDS name -> the name used downstream. The right-hand side matches the columns
# of the prediction grid so the new values drop into the existing schema.
# Only what the models actually use. skin_reservoir_content and
# surface_solar_radiation_downwards appear in the collaborator's daily extract
# but in no model formula, so they are off by default; EXTRA=1 adds them back
# at ~40% more download.
SHORT <- c("2m_temperature" = "t2m", "2m_dewpoint_temperature" = "d2m",
           "10m_u_component_of_wind" = "u10", "10m_v_component_of_wind" = "v10",
           "total_precipitation" = "tp")
if (nzchar(Sys.getenv("EXTRA"))) {
  SHORT <- c(SHORT, "surface_solar_radiation_downwards" = "ssrd",
             "skin_reservoir_content" = "src")
}
# CDS enforces a per-account cap on simultaneous requests. Six at once were all
# rejected on 2026-08-20; two run fine. "rejected" carries no message field, so
# it cannot be distinguished from a real failure by content - it is treated as
# retryable and the year goes back on the queue.
MAXJOBS <- as.integer(Sys.getenv("MAXJOBS", "2"))
MAXTRY  <- 4L

#---- Credentials -------------------------------------------------------------

cds_key <- function() {
  k <- Sys.getenv("CDSAPI_KEY")
  if (nzchar(k)) return(k)
  rc <- path.expand("~/.cdsapirc")
  if (!file.exists(rc)) {
    stop("No CDS credentials. Register at https://cds.climate.copernicus.eu, ",
         "accept the ERA5-Land licence, then write ~/.cdsapirc:\n",
         "  url: https://cds.climate.copernicus.eu/api\n  key: <token>",
         call. = FALSE)
  }
  ln <- readLines(rc, warn = FALSE)
  k  <- sub("^\\s*key\\s*:\\s*", "", grep("^\\s*key\\s*:", ln, value = TRUE))[1]
  if (is.na(k) || !nzchar(k)) stop("~/.cdsapirc has no 'key:' line", call. = FALSE)
  str_squish(k)
}
key <- cds_key()

#---- CDS, as three separate steps -------------------------------------------
# Submitting and collecting are split so that many years can be in the queue at
# once. Queue latency (~65 s) dominates a single year, so running MAXJOBS of
# them concurrently is most of the speed-up.

cds_submit <- function(year, area) {
  body <- list(
    product_type    = list("monthly_averaged_reanalysis"),
    variable        = as.list(names(SHORT)),
    year            = list(as.character(year)),
    month           = as.list(sprintf("%02d", 1:12)),
    time            = list("00:00"),
    data_format     = "netcdf",
    download_format = "unarchived",
    area            = as.list(area)          # N, W, S, E
  )
  resp <- request(paste0(CDS_BASE, "/retrieve/v1/processes/", COLLECTION, "/execute")) %>%
    req_headers(`PRIVATE-TOKEN` = key, Accept = "application/json") %>%
    req_body_json(list(inputs = body)) %>%
    req_error(is_error = function(r) FALSE) %>%
    req_perform()
  if (resp_status(resp) >= 400) {
    stop("CDS rejected ", year, " (HTTP ", resp_status(resp), "): ",
         resp_body_string(resp), call. = FALSE)
  }
  resp_body_json(resp)$jobID
}

cds_status <- function(job) {
  request(paste0(CDS_BASE, "/retrieve/v1/jobs/", job)) %>%
    req_headers(`PRIVATE-TOKEN` = key, Accept = "application/json") %>%
    req_error(is_error = function(r) FALSE) %>%
    req_perform() %>% resp_body_json()
}

cds_fetch <- function(job, dest) {
  res <- request(paste0(CDS_BASE, "/retrieve/v1/jobs/", job, "/results")) %>%
    req_headers(`PRIVATE-TOKEN` = key, Accept = "application/json") %>%
    req_perform() %>% resp_body_json()
  href <- res$asset$value$href
  if (is.null(href)) stop("no download href", call. = FALSE)
  request(href) %>% req_perform(path = dest)
  as.numeric(res$asset$value$`file:size`)
}

#---- Physics -----------------------------------------------------------------

# Saturation vapour pressure over water, kPa, from temperature in Celsius.
# Tetens/Magnus, the same form used to build the original ERA5 extract.
svp <- function(tc) 0.6108 * exp(17.27 * tc / (tc + 237.3))

#---- Point sets --------------------------------------------------------------

message("Loading the two point sets ...")
samp_pts <- readRDS(gmtp_in) %>% distinct(Longitude, Latitude)
grid_pts <- readRDS(grid_in) %>% distinct(Longitude, Latitude)
message("  sample: ", format(nrow(samp_pts), big.mark = ","), " coordinates")
message("  grid:   ", format(nrow(grid_pts), big.mark = ","), " coordinates")
samp_xy <- as.matrix(samp_pts[, c("Longitude", "Latitude")])
grid_xy <- as.matrix(grid_pts[, c("Longitude", "Latitude")])

# Request only the latitude band the grid actually occupies, padded by one cell
# so bilinear interpolation at the edge still has neighbours on both sides.
AREA <- c(min(90, ceiling(max(grid_pts$Latitude) * 10) / 10 + 0.2), -180,
          max(-90, floor(min(grid_pts$Latitude) * 10) / 10 - 0.2),  180)
message("  area (N,W,S,E): ", paste(AREA, collapse = ", "))

years <- if (nzchar(YEARS_ENV)) as.integer(str_split_1(YEARS_ENV, ",")) else {
  yr <- readRDS(gmtp_in)$year
  seq(min(yr, na.rm = TRUE), max(yr, na.rm = TRUE))
}

todo <- years[FORCE | !(file.exists(file.path(ras_dir, sprintf("era5_sample_%d.rds", years))) &
                        file.exists(file.path(ras_dir, sprintf("era5_grid_%d.rds",   years))))]
message("Years: ", paste(range(years), collapse = "-"), " (", length(years),
        "), ", length(todo), " to fetch")
if (!length(todo)) { message("Nothing to do."); quit(save = "no") }

#---- Per-year processing -----------------------------------------------------

process_year <- function(nc, yr) {
  r <- rast(nc)
  # Layers are named "<short>_valid_time=<epoch>", e.g. "t2m_valid_time=1388534400".
  # Both the variable and the month have to come from that string:
  #   - terra::varnames() returns ONE entry PER VARIABLE (length 5), not per
  #     layer (length 60), so using its index as a layer index silently selects
  #     the wrong field - which is exactly the bug this replaced.
  #   - terra::time() is NA on these files, so the month cannot come from there.
  # Sorting on the epoch is what guarantees layer i is month i for every
  # variable, rather than trusting the order CDS happens to write.
  lyr_var <- sub("_valid_time=.*$", "", names(r))
  ep <- as.numeric(str_match(names(r), "valid_time=(\\d+)")[, 2])
  if (anyNA(ep)) stop("cannot read valid_time from layer names", call. = FALSE)
  mn <- format(as.POSIXct(ep, origin = "1970-01-01", tz = "UTC"), "%Y-%m")

  months_present <- sort(unique(mn))
  n_month <- length(months_present)
  dim_days <- as.integer(days_in_month(as.Date(paste0(months_present, "-01"))))
  wt <- dim_days / sum(dim_days)
  if (n_month < 12) message("     PARTIAL YEAR (", n_month, " months) - sums set to NA")

  grab <- function(short) {
    idx <- which(lyr_var == short)
    if (!length(idx)) stop("variable ", short, " not in the netCDF", call. = FALSE)
    if (length(idx) != n_month) {
      stop(short, ": ", length(idx), " layers for ", n_month, " months", call. = FALSE)
    }
    r[[idx[order(ep[idx])]]]
  }
  st <- lapply(SHORT, grab)
  names(st) <- names(SHORT)

  wmean <- function(s) app(s * wt, sum, na.rm = FALSE)

  t2m <- wmean(st$`2m_temperature`) - 273.15
  d2m <- wmean(st$`2m_dewpoint_temperature`) - 273.15
  u10 <- wmean(st$`10m_u_component_of_wind`)
  v10 <- wmean(st$`10m_v_component_of_wind`)

  # Monthly RH and VPD first, then averaged - see the header note on Jensen.
  tc_m <- st$`2m_temperature` - 273.15
  dc_m <- st$`2m_dewpoint_temperature` - 273.15
  es_m <- svp(tc_m); ea_m <- svp(dc_m)
  # clamp(), not min()/max(): with a scalar first argument those dispatch to
  # the base generic and hand app() an S4 object it cannot use.
  rh   <- app(clamp(100 * ea_m / es_m, upper = 100, values = TRUE) * wt,
              sum, na.rm = FALSE)
  vpd  <- app(clamp(es_m - ea_m, lower = 0, values = TRUE) * wt,
              sum, na.rm = FALSE)

  # Accumulations: the monthly value is a mean daily rate.
  if (n_month == 12) {
    tp <- app(st$total_precipitation * dim_days, sum, na.rm = FALSE) * 1000
  } else {
    tp <- t2m * NA
  }

  ws <- sqrt(u10^2 + v10^2)                        # VECTOR mean, see header
  wd <- (atan2(u10, v10) * 180 / pi) %% 360

  ann <- c(t2m, d2m, rh, vpd, ws, wd, tp)
  nms <- c("temperature_2m", "dewpoint_temperature_2m", "relative_humidity",
           "vpd", "wind_speed", "wind_direction", "total_precipitation_sum")
  if (!is.null(st$surface_solar_radiation_downwards)) {
    ssrd <- if (n_month == 12)
      app(st$surface_solar_radiation_downwards * dim_days, sum, na.rm = FALSE) else t2m * NA
    ann <- c(ann, ssrd); nms <- c(nms, "surface_solar_radiation_downwards_sum")
  }
  if (!is.null(st$skin_reservoir_content)) {
    ann <- c(ann, wmean(st$skin_reservoir_content) * 1000)
    nms <- c(nms, "skin_reservoir_content")
  }
  names(ann) <- nms
  ann
}

#---- Pipeline: keep MAXJOBS in the CDS queue, process as they land -----------

queue   <- todo
pending <- list()                # year -> jobID
tries   <- list()                # year -> submissions so far
t_start <- Sys.time()

repeat {
  while (length(pending) < MAXJOBS && length(queue)) {
    yr <- queue[1]; queue <- queue[-1]
    job <- tryCatch(cds_submit(yr, AREA),
                    error = function(e) { message(yr, ": ", conditionMessage(e)); NA })
    if (!is.na(job)) {
      pending[[as.character(yr)]] <- job
      message("submitted ", yr, "  (", length(pending), " in flight, ",
              length(queue), " queued)")
    } else {
      queue <- c(queue, yr)
    }
    Sys.sleep(2)                     # do not burst-submit into the cap
  }
  if (!length(pending)) break

  ready <- NULL
  for (y in names(pending)) {
    st <- cds_status(pending[[y]])
    if (identical(st$status, "successful")) { ready <- y; break }
    if (st$status %in% c("failed", "dismissed", "rejected")) {
      yi <- as.integer(y)
      tries[[y]] <- (tries[[y]] %||% 0L) + 1L
      msg <- paste(unlist(st$message), collapse = " ")
      pending[[y]] <- NULL
      if (tries[[y]] < MAXTRY) {
        message(y, ": ", st$status, " (attempt ", tries[[y]], ") - requeued", 
                if (nzchar(msg)) paste0(": ", msg) else "")
        queue <- c(queue, yi)
        Sys.sleep(20)                # let the queue drain before retrying
      } else {
        message(y, ": ", st$status, " after ", MAXTRY, " attempts - giving up")
      }
      ready <- NA; break
    }
  }
  if (is.null(ready)) { Sys.sleep(10); next }
  if (is.na(ready))   next

  yr  <- as.integer(ready)
  job <- pending[[ready]]
  pending[[ready]] <- NULL
  nc  <- file.path(tmp_dir, sprintf("era5land_%d.nc", yr))

  message("\n", yr, ":")
  if (file.exists(nc) && file.size(nc) > 1e6) {
    message("     reusing ", basename(nc), " already on disk")
  } else {
    t0 <- Sys.time()
    sz <- cds_fetch(job, nc)
    message("     downloaded ", round(sz / 1e6), " MB in ",
            round(difftime(Sys.time(), t0, units = "secs")), " s")
  }

  t0  <- Sys.time()
  ann <- process_year(nc, yr)
  message("     aggregated in ", round(difftime(Sys.time(), t0, units = "secs")), " s")

  t0  <- Sys.time()
  e_s <- terra::extract(ann, samp_xy, method = "bilinear")
  e_g <- terra::extract(ann, grid_xy, method = "bilinear")
  message("     extracted ", format(nrow(grid_xy), big.mark = ","), " points in ",
          round(difftime(Sys.time(), t0, units = "secs")), " s")

  saveRDS(bind_cols(samp_pts, as_tibble(e_s)) %>% mutate(year = yr),
          file.path(ras_dir, sprintf("era5_sample_%d.rds", yr)))
  saveRDS(bind_cols(grid_pts, as_tibble(e_g)) %>% mutate(year = yr),
          file.path(ras_dir, sprintf("era5_grid_%d.rds", yr)))
  message("     wrote era5_sample_", yr, ".rds and era5_grid_", yr, ".rds")

  rm(ann); gc(FALSE)
  unlink(nc)
}

message("\nDone in ", round(difftime(Sys.time(), t_start, units = "mins"), 1),
        " min. Per-year files are in ", ras_dir)
message("Merge them with 11_merge_annual_covariates.R")
