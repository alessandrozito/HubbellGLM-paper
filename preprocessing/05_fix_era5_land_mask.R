# Repair the sample points that ERA5-Land's land-sea mask leaves empty: coastal
# points snap to the nearest native cell, oceanic islands fall back to ERA5
# single-levels.
#
# Produces: cache/annual/era5_sample_YYYY.rds (in place), qc/qc_era5_*.tsv

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

library(tidyverse)
library(httr2)

data_dir <- DATA
ras_dir <- RASTER
cell_dir<- file.path(CACHE, "era5_daily_cells")
tmp_dir <- file.path(CACHE, "era5_ts_tmp")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

SNAP_MAX_KM <- as.numeric(Sys.getenv("SNAP_MAX_KM", "30"))
ONLY        <- Sys.getenv("ONLY", "")
DATE_RANGE  <- "1991-01-01/2024-12-31"
VARS <- c("2m_temperature", "2m_dewpoint_temperature",
          "10m_u_component_of_wind", "10m_v_component_of_wind",
          "total_precipitation")

#---- CDS ---------------------------------------------------------------------

cds_key <- function() {
  k <- Sys.getenv("CDSAPI_KEY"); if (nzchar(k)) return(k)
  ln <- readLines(path.expand("~/.cdsapirc"), warn = FALSE)
  str_squish(sub("^\\s*key\\s*:\\s*", "", grep("^\\s*key\\s*:", ln, value = TRUE))[1])
}
key      <- cds_key()
CDS_BASE <- "https://cds.climate.copernicus.eu/api"
safe_json <- function(expr) tryCatch(expr, error = function(e) NULL)

cds_submit <- function(collection, lon, lat) {
  resp <- request(paste0(CDS_BASE, "/retrieve/v1/processes/", collection, "/execute")) %>%
    req_headers(`PRIVATE-TOKEN` = key, Accept = "application/json") %>%
    req_body_json(list(inputs = list(
      variable = as.list(VARS),
      location = list(longitude = lon, latitude = lat),
      date = list(DATE_RANGE), data_format = "csv"))) %>%
    req_error(is_error = function(r) FALSE) %>% req_perform()
  if (resp_status(resp) >= 400)
    stop("HTTP ", resp_status(resp), ": ", resp_body_string(resp), call. = FALSE)
  resp_body_json(resp)$jobID
}

cds_status <- function(job) {
  r <- safe_json(request(paste0(CDS_BASE, "/retrieve/v1/jobs/", job)) %>%
    req_headers(`PRIVATE-TOKEN` = key, Accept = "application/json") %>%
    req_error(is_error = function(r) FALSE) %>%
    req_timeout(120) %>% req_perform() %>% resp_body_json())
  if (is.null(r) || is.null(r$status)) list(status = "unknown") else r
}

cds_fetch <- function(job, dest) {
  res <- safe_json(request(paste0(CDS_BASE, "/retrieve/v1/jobs/", job, "/results")) %>%
    req_headers(`PRIVATE-TOKEN` = key, Accept = "application/json") %>%
    req_timeout(120) %>% req_perform() %>% resp_body_json())
  if (is.null(res$asset$value$href)) stop("no download href", call. = FALSE)
  request(res$asset$value$href) %>% req_timeout(900) %>% req_perform(path = dest)
  invisible(dest)
}

svp <- function(tc) 0.6108 * exp(17.27 * tc / (tc + 237.3))   # kPa

reduce_to_daily <- function(zipf, tag) {
  ex <- file.path(tmp_dir, paste0("x_", tag))
  unlink(ex, recursive = TRUE); dir.create(ex, recursive = TRUE)
  parts <- lapply(list.files(ex <- {unzip(zipf, exdir = ex); ex}, full.names = TRUE),
                  function(f) read_csv(f, show_col_types = FALSE, progress = FALSE) %>%
                    dplyr::select(-any_of(c("latitude", "longitude"))))
  h <- reduce(parts, full_join, by = "valid_time")
  unlink(ex, recursive = TRUE)
  h %>%
    mutate(date = as.Date(valid_time), tc = t2m - 273.15, dc = d2m - 273.15,
           rh_h = pmin(100, 100 * svp(dc) / svp(tc))) %>%
    group_by(date) %>%
    summarise(temperature_2m          = mean(tc,   na.rm = TRUE),
              dewpoint_temperature_2m = mean(dc,   na.rm = TRUE),
              relative_humidity       = mean(rh_h, na.rm = TRUE),
              u10 = mean(u10, na.rm = TRUE), v10 = mean(v10, na.rm = TRUE),
              total_precipitation_sum = if (all(is.na(tp))) NA_real_
                                        else sum(tp, na.rm = TRUE) * 1000,
              .groups = "drop") %>%
    mutate(wind_speed = sqrt(u10^2 + v10^2),
           vpd = pmax(0, svp(temperature_2m) - svp(dewpoint_temperature_2m))) %>%
    dplyr::select(-u10, -v10)
}

# Serial: this is a handful of cells, not 818, and the queue is friendlier.
fetch_one <- function(collection, lon, lat, tag) {
  zipf <- file.path(tmp_dir, paste0(tag, ".zip"))
  job  <- cds_submit(collection, lon, lat)
  repeat {
    st <- cds_status(job)$status
    if (identical(st, "successful")) break
    if (st %in% c("failed", "dismissed")) stop("job ", st, call. = FALSE)
    Sys.sleep(10)
  }
  cds_fetch(job, zipf)
  out <- reduce_to_daily(zipf, tag)
  unlink(zipf)
  out
}

#---- Which cells are masked, and how far is valid land? ----------------------

# The mask MUST come from ERA5-Land's native grid, not from era5_grid_*.rds.
# The stored grid samples at 0.05-offset points (-179.95, -179.85, ...) while
# ERA5-Land's cell centres are exact multiples of 0.1. A snap target taken from
# the stored grid therefore sits on the corner of four native cells, and the
# timeseries API resolves it to whichever one it likes - which for a coastal
# site is a coin flip between the land cell and the ocean cell that started the
# problem. Four of the first six snaps came back all-NA that way.
#
# era5land_landmask.tif is one month of 2m_temperature downloaded globally; its
# non-NA footprint is the mask itself, on the native lattice.
mask_f <- file.path(data_dir, "era5land_landmask.tif")
if (!file.exists(mask_f))
  stop("missing ", mask_f, " - see scratchpad/get_mask.R", call. = FALSE)
suppressMessages(library(terra))
mk <- terra::rast(mask_f)
land <- as.data.frame(mk, xy = TRUE) %>% setNames(c("Longitude", "Latitude", "v")) %>%
  filter(v > 0) %>%
  mutate(Longitude = ((Longitude + 180) %% 360) - 180,   # raster is 0..360
         Longitude = round(Longitude, 1), Latitude = round(Latitude, 1)) %>%
  dplyr::select(Longitude, Latitude)
message("native ERA5-Land cells with data: ", format(nrow(land), big.mark = ","))

# Equirectangular distance: exact enough at the tens-of-km scale that decides
# coastal-vs-island, and it keeps this a plain vector search.
nearest_land <- function(lon, lat) {
  map2_dfr(lon, lat, function(x, y) {
    sub <- land %>% filter(abs(Longitude - x) < 6, abs(Latitude - y) < 6)
    if (!nrow(sub)) sub <- land
    dx <- (sub$Longitude - x) * cos(y * pi / 180); dy <- sub$Latitude - y
    k  <- which.min(dx^2 + dy^2)
    tibble(snap_lon = sub$Longitude[k], snap_lat = sub$Latitude[k],
           snap_km  = sqrt(dx[k]^2 + dy[k]^2) * 111.32)
  })
}

#---- A. Daily cells ----------------------------------------------------------

if (!identical(ONLY, "annual")) {
  # ^c...\.rds only - never the .masked.rds backups this script writes.
  files <- list.files(cell_dir, "^c[-+0-9.]+\\.rds$", full.names = TRUE)
  bad <- keep(files, ~ all(is.na(readRDS(.x)$temperature_2m)))
  # A cell whose .rds was renamed to .masked.rds but whose replacement never
  # landed (interrupted run) has no .rds to test, so pick it up by its backup.
  bad <- union(bad, sub("\\.masked\\.rds$", ".rds",
                        list.files(cell_dir, "\\.masked\\.rds$", full.names = TRUE)))
  bad <- keep(bad, ~ !file.exists(.x) || all(is.na(readRDS(.x)$temperature_2m)))
  message("land-masked daily cells: ", length(bad), " of ", length(files))

  if (length(bad)) {
    b <- tibble(path = bad, cell = tools::file_path_sans_ext(basename(bad))) %>%
      mutate(lon = as.numeric(substr(cell, 2, 9)), lat = as.numeric(substr(cell, 10, 16))) %>%
      bind_cols(nearest_land(.$lon, .$lat)) %>%
      mutate(mode = if_else(snap_km <= SNAP_MAX_KM, "snap", "era5"))
    message("  coastal (snap <= ", SNAP_MAX_KM, " km): ", sum(b$mode == "snap"),
            "   oceanic (ERA5 single levels): ", sum(b$mode == "era5"))

    for (i in seq_len(nrow(b))) {
      # Preserve the original as .masked.rds once, so a rerun can start over.
      keepf <- sub("\\.rds$", ".masked.rds", b$path[i])
      if (!file.exists(keepf)) file.rename(b$path[i], keepf)
      if (file.exists(b$path[i]) &&
          !all(is.na(readRDS(b$path[i])$temperature_2m))) next   # already repaired
      unlink(b$path[i])
      # A snap can still come back empty: the monthly-mean mask and the
      # timeseries service do not agree cell-for-cell on the very smallest
      # islands. Rather than leave a hole, fall back to unmasked ERA5 at the
      # true coordinate - the same treatment the oceanic sites get.
      grab <- function(mode) {
        ts <- if (mode == "snap")
          fetch_one("reanalysis-era5-land-timeseries", b$snap_lon[i], b$snap_lat[i], b$cell[i])
        else
          fetch_one("reanalysis-era5-single-levels-timeseries", b$lon[i], b$lat[i],
                    paste0(b$cell[i], "_sl"))
        if (all(is.na(ts$temperature_2m)))
          stop("fill is itself all-NA", call. = FALSE)
        ts
      }
      used <- b$mode[i]
      ts <- tryCatch(grab(b$mode[i]),
                     error = function(e) { message("  . ", b$cell[i], " ", b$mode[i],
                                                   " failed (", conditionMessage(e),
                                                   ") - retrying unmasked ERA5"); NULL })
      if (is.null(ts) && b$mode[i] == "snap") {
        used <- "era5_fallback"
        ts <- tryCatch(grab("era5"), error = function(e) {
          message("  ! ", b$cell[i], ": ", conditionMessage(e)); NULL })
      }
      ok <- !is.null(ts)
      if (ok) {
        attr(ts, "fill") <- list(mode = used, snap_km = b$snap_km[i],
                                 lon = b$snap_lon[i], lat = b$snap_lat[i])
        saveRDS(ts, b$path[i])
      }
      b$mode[i] <- used
      message("  ", i, "/", nrow(b), " ", b$cell[i], " [", used, " ",
              round(b$snap_km[i], 1), " km] ", if (ok) "ok" else "FAILED")
    }
    write_tsv(b %>% dplyr::select(-path), file.path(QC, "qc_era5_landmask_cells.tsv"))
  }
}

#---- B. Annual sample values -------------------------------------------------

if (!identical(ONLY, "daily")) {
  message("\nannual extract")
  yrs <- list.files(ras_dir, "^era5_sample_(\\d{4})\\.rds$") %>%
    str_match("(\\d{4})") %>% {.[, 2]} %>% as.integer() %>% sort()
  EVARS <- c("temperature_2m", "dewpoint_temperature_2m", "relative_humidity",
             "vpd", "wind_speed", "wind_direction", "total_precipitation_sum")

  # The oceanic sites have no valid grid cell to borrow from, so they come from
  # the ERA5 single-levels daily series written in part A, averaged to the year.
  ocean_annual <- function(lon, lat, yr) {
    cid <- sprintf("c%+08.3f%+07.3f", round(lon / 0.1) * 0.1, round(lat / 0.1) * 0.1)
    f <- file.path(cell_dir, paste0(cid, ".rds"))
    if (!file.exists(f)) return(NULL)
    ts <- readRDS(f) %>% filter(lubridate::year(date) == yr)
    if (!nrow(ts) || all(is.na(ts$temperature_2m))) return(NULL)
    tibble(temperature_2m = mean(ts$temperature_2m, na.rm = TRUE),
           dewpoint_temperature_2m = mean(ts$dewpoint_temperature_2m, na.rm = TRUE),
           relative_humidity = mean(ts$relative_humidity, na.rm = TRUE),
           vpd = mean(ts$vpd, na.rm = TRUE),
           wind_speed = mean(ts$wind_speed, na.rm = TRUE),
           wind_direction = NA_real_,
           total_precipitation_sum = sum(ts$total_precipitation_sum, na.rm = TRUE))
  }

  log <- list()
  for (yy in yrs) {
    fs <- file.path(ras_dir, sprintf("era5_sample_%d.rds", yy))
    s  <- readRDS(fs)
    if (!"era5_src" %in% names(s)) s$era5_src <- if_else(is.na(s$temperature_2m), NA_character_, "observed")
    miss <- which(is.na(s$temperature_2m))
    if (!length(miss)) next
    g <- readRDS(file.path(ras_dir, sprintf("era5_grid_%d.rds", yy))) %>%
      filter(!is.na(temperature_2m))
    for (i in miss) {
      x <- s$Longitude[i]; y <- s$Latitude[i]
      dx <- (g$Longitude - x) * cos(y * pi / 180); dy <- g$Latitude - y
      k  <- which.min(dx^2 + dy^2); km <- sqrt(dx[k]^2 + dy[k]^2) * 111.32
      if (km <= SNAP_MAX_KM) {
        s[i, EVARS] <- g[k, EVARS]; s$era5_src[i] <- "nearest_cell"
      } else {
        v <- ocean_annual(x, y, yy)
        if (!is.null(v)) { s[i, EVARS] <- v[, EVARS]; s$era5_src[i] <- "era5_single_levels" }
      }
      log[[length(log) + 1]] <- tibble(year = yy, Longitude = x, Latitude = y,
                                       snap_km = km, src = s$era5_src[i])
    }
    saveRDS(s, fs)
  }
  if (length(log)) {
    L <- bind_rows(log)
    write_tsv(L, file.path(QC, "qc_era5_annual_fills.tsv"))
    nc <- sum(L$src == "nearest_cell", na.rm = TRUE)
    message("  attempted ", nrow(L), " coordinate-years, filled ",
            sum(!is.na(L$src)))
    print(as.data.frame(L %>% count(src, name = "n")), row.names = FALSE)
    if (nc)
      message("  snap distance for nearest_cell fills: median ",
              round(median(L$snap_km[L$src == "nearest_cell"]), 1), " km, max ",
              round(max(L$snap_km[L$src == "nearest_cell"]), 1), " km")
    # Neither route can reach an oceanic site in a year the daily series does
    # not cover: DATE_RANGE stops at the last complete ERA5 year, so anything
    # past it stays NA and is dropped by YEAR_MAX in the analysis.
    if (anyNA(L$src)) {
      u <- L %>% filter(is.na(src)) %>% count(year, name = "n")
      message("  UNFILLED (past the ", sub(".*/", "", DATE_RANGE),
              " end of the daily series): ",
              paste0(u$year, " x", u$n, collapse = ", "))
    }
  }
  left <- map_int(yrs, ~ sum(is.na(readRDS(file.path(ras_dir, sprintf("era5_sample_%d.rds", .x)))$temperature_2m)))
  message("  remaining NA per year: ", paste(left, collapse = " "))
}
