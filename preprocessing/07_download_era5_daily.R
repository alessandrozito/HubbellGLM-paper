# Daily ERA5-Land series at every sample location, 1991-2024, for the
# collection-window weather anomalies.
#
# Produces: cache/era5_daily_cells/

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

library(tidyverse)
library(httr2)

data_dir <- DATA
cell_dir<- file.path(CACHE, "era5_daily_cells")
tmp_dir <- file.path(CACHE, "era5_ts_tmp")
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir,  recursive = TRUE, showWarnings = FALSE)

CDS_BASE   <- "https://cds.climate.copernicus.eu/api"
COLLECTION <- "reanalysis-era5-land-timeseries"
DATE_RANGE <- "1991-01-01/2024-12-31"
VARS <- c("2m_temperature", "2m_dewpoint_temperature",
          "10m_u_component_of_wind", "10m_v_component_of_wind",
          "total_precipitation")
MAXJOBS <- as.integer(Sys.getenv("MAXJOBS", "2"))
MAXTRY  <- 4L

cds_key <- function() {
  k <- Sys.getenv("CDSAPI_KEY"); if (nzchar(k)) return(k)
  ln <- readLines(path.expand("~/.cdsapirc"), warn = FALSE)
  str_squish(sub("^\\s*key\\s*:\\s*", "", grep("^\\s*key\\s*:", ln, value = TRUE))[1])
}
key <- cds_key()

cds_submit <- function(lon, lat) {
  body <- list(variable = as.list(VARS),
               location = list(longitude = lon, latitude = lat),
               date = list(DATE_RANGE), data_format = "csv")
  resp <- request(paste0(CDS_BASE, "/retrieve/v1/processes/", COLLECTION, "/execute")) %>%
    req_headers(`PRIVATE-TOKEN` = key, Accept = "application/json") %>%
    req_body_json(list(inputs = body)) %>%
    req_error(is_error = function(r) FALSE) %>% req_perform()
  if (resp_status(resp) >= 400)
    stop("HTTP ", resp_status(resp), ": ", resp_body_string(resp), call. = FALSE)
  resp_body_json(resp)$jobID
}
# CDS returns an nginx 502 HTML page often enough over a multi-hour run that
# every call has to tolerate it. resp_body_json() throws on unparseable
# content, so an unguarded poll kills the whole job - which is exactly what
# happened on the first attempt, four cells in. A 502 is transient: report it
# as "unknown" and let the caller poll again.
safe_json <- function(expr) tryCatch(expr, error = function(e) NULL)

cds_status <- function(job) {
  r <- safe_json(
    request(paste0(CDS_BASE, "/retrieve/v1/jobs/", job)) %>%
      req_headers(`PRIVATE-TOKEN` = key, Accept = "application/json") %>%
      req_error(is_error = function(r) FALSE) %>%
      req_timeout(120) %>% req_perform() %>% resp_body_json())
  if (is.null(r) || is.null(r$status)) list(status = "unknown") else r
}

cds_fetch <- function(job, dest) {
  res <- safe_json(
    request(paste0(CDS_BASE, "/retrieve/v1/jobs/", job, "/results")) %>%
      req_headers(`PRIVATE-TOKEN` = key, Accept = "application/json") %>%
      req_timeout(120) %>% req_perform() %>% resp_body_json())
  if (is.null(res$asset$value$href)) stop("no download href", call. = FALSE)
  request(res$asset$value$href) %>% req_timeout(900) %>% req_perform(path = dest)
  invisible(dest)
}

svp <- function(tc) 0.6108 * exp(17.27 * tc / (tc + 237.3))   # kPa

# Reduce one downloaded archive to a daily series.
reduce_to_daily <- function(zipf) {
  fs <- unzip(zipf, list = TRUE)$Name
  ex <- file.path(tmp_dir, "x"); unlink(ex, recursive = TRUE); dir.create(ex)
  unzip(zipf, exdir = ex)
  parts <- lapply(list.files(ex, full.names = TRUE), function(f)
    read_csv(f, show_col_types = FALSE, progress = FALSE) %>%
      dplyr::select(-any_of(c("latitude", "longitude"))))
  h <- reduce(parts, full_join, by = "valid_time")
  unlink(ex, recursive = TRUE)

  h %>%
    mutate(date = as.Date(valid_time),
           tc = t2m - 273.15, dc = d2m - 273.15,
           rh_h = pmin(100, 100 * svp(dc) / svp(tc))) %>%
    group_by(date) %>%
    summarise(temperature_2m        = mean(tc,  na.rm = TRUE),
              dewpoint_temperature_2m = mean(dc, na.rm = TRUE),
              relative_humidity     = mean(rh_h, na.rm = TRUE),
              u10 = mean(u10, na.rm = TRUE), v10 = mean(v10, na.rm = TRUE),
              # na.rm on an all-NA day returns 0, which silently turns a
              # land-masked (ocean) cell into a permanently dry one.
              total_precipitation_sum = if (all(is.na(tp))) NA_real_
                                        else sum(tp, na.rm = TRUE) * 1000,
              .groups = "drop") %>%
    mutate(wind_speed = sqrt(u10^2 + v10^2),
           vpd = pmax(0, svp(temperature_2m) - svp(dewpoint_temperature_2m))) %>%
    dplyr::select(-u10, -v10)
}

#---- Cells -------------------------------------------------------------------

d <- readRDS(file.path(data_dir, "GMTP_sample_clean.rds"))
cells <- d %>%
  mutate(clon = round(Longitude / 0.1) * 0.1, clat = round(Latitude / 0.1) * 0.1) %>%
  distinct(clon, clat) %>%
  mutate(id = sprintf("c%+08.3f%+07.3f", clon, clat))
message("sample coordinates: ", nrow(distinct(d, Longitude, Latitude)),
        "  ->  ERA5-Land cells: ", nrow(cells))

todo <- cells %>% filter(!file.exists(file.path(cell_dir, paste0(id, ".rds"))))
message("cells still to fetch: ", nrow(todo), " of ", nrow(cells))
if (!nrow(todo)) { message("Nothing to do."); quit(save = "no") }

#---- Pipeline ----------------------------------------------------------------

queue <- seq_len(nrow(todo)); pending <- list(); tries <- list()
t_start <- Sys.time(); done <- 0L

repeat {
  while (length(pending) < MAXJOBS && length(queue)) {
    i <- queue[1]; queue <- queue[-1]
    job <- tryCatch(cds_submit(todo$clon[i], todo$clat[i]),
                    error = function(e) {
                      message("  submit failed for ", todo$id[i], ": ",
                              substr(gsub("\\s+", " ", conditionMessage(e)), 1, 120))
                      NA })
    if (!is.na(job)) pending[[as.character(i)]] <- job else queue <- c(queue, i)
    Sys.sleep(2)
  }
  if (!length(pending)) break

  ready <- NULL
  for (k in names(pending)) {
    st <- cds_status(pending[[k]])
    if (identical(st$status, "unknown")) next        # transient 502, poll again
    if (identical(st$status, "successful")) { ready <- k; break }
    if (st$status %in% c("failed", "dismissed", "rejected")) {
      tries[[k]] <- (tries[[k]] %||% 0L) + 1L
      pending[[k]] <- NULL
      if (tries[[k]] < MAXTRY) { queue <- c(queue, as.integer(k)); Sys.sleep(20) }
      else message("  cell ", todo$id[as.integer(k)], ": giving up after ", MAXTRY)
      ready <- NA; break
    }
  }
  if (is.null(ready)) { Sys.sleep(10); next }
  if (is.na(ready))   next

  i <- as.integer(ready); job <- pending[[ready]]; pending[[ready]] <- NULL
  zipf <- file.path(tmp_dir, paste0(todo$id[i], ".zip"))
  ok <- tryCatch({ cds_fetch(job, zipf)
                   saveRDS(reduce_to_daily(zipf), file.path(cell_dir, paste0(todo$id[i], ".rds")))
                   TRUE },
                 error = function(e) { message("  ", todo$id[i], ": ", conditionMessage(e)); FALSE })
  unlink(zipf)
  if (ok) {
    done <- done + 1L
    el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    if (done %% 10 == 0 || done == 1)
      message("  ", done, "/", nrow(todo), " cells  (", round(el, 1), " min elapsed, ",
              round(el / done * (nrow(todo) - done) / 60, 1), " h remaining)")
  }
}
message("\nDone in ", round(difftime(Sys.time(), t_start, units = "hours"), 2),
        " h. Daily cell files in ", cell_dir)
message("Next: 13_compute_anomalies.R")
