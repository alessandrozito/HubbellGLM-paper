# Collection-window weather anomalies against each cell's own 1991-2020
# day-of-year normal.
#
# Produces: GMTP_anomalies.rds

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

library(tidyverse)
library(lubridate)

data_dir <- DATA
cell_dir <- file.path(CACHE, "era5_daily_cells")

BASE_FROM    <- 1991L
BASE_TO      <- 2020L
SMOOTH_DAYS  <- as.integer(Sys.getenv("SMOOTH_DAYS", "7"))
VARS <- c("temperature_2m", "total_precipitation_sum", "relative_humidity", "wind_speed")

d <- readRDS(file.path(data_dir, "GMTP_sample_clean.rds")) %>%
  mutate(clon = round(Longitude / 0.1) * 0.1, clat = round(Latitude / 0.1) * 0.1,
         cell = sprintf("c%+08.3f%+07.3f", clon, clat),
         collection_start_date = as.Date(collection_start_date),
         collection_end_date   = as.Date(collection_end_date),
         year = lubridate::year(as.Date(collection_start_date)))

have <- tools::file_path_sans_ext(list.files(cell_dir, "\\.rds$"))
message("cells with daily data: ", length(have), " of ", n_distinct(d$cell))
miss <- setdiff(unique(d$cell), have)
if (length(miss)) message("  ", length(miss), " cells missing - those fieldids get NA")

# Circular +/- SMOOTH_DAYS day-of-year window. Row dd of WIN holds every day of
# year inside dd's window, so the normal is one matrix gather per cell rather
# than 366 filters, and late December borrows from early January.
WIN  <- outer(1:366, -SMOOTH_DAYS:SMOOTH_DAYS,
              function(a, b) ((a - 1 + b) %% 366) + 1L)
ANML <- paste0(VARS, "_anml")
rows_by_cell <- split(seq_len(nrow(d)), d$cell)

anom_one_cell <- function(cid) {
  f <- file.path(cell_dir, paste0(cid, ".rds"))
  if (!file.exists(f)) return(NULL)
  ts   <- readRDS(f)
  doy  <- yday(ts$date)
  keep <- year(ts$date) >= BASE_FROM & year(ts$date) <= BASE_TO
  if (!any(keep)) return(NULL)
  OBS <- as.matrix(ts[, VARS])

  # 1991-2020 normal per day of year, then the smoothing window. rowsum() over
  # the day of year gives the sums and the counts that na.rm implies, and the
  # window is applied by indexing the 366-row table with WIN.
  B  <- OBS[keep, , drop = FALSE]
  M  <- rowsum(ifelse(is.na(B), 0, B), doy[keep], reorder = TRUE)
  Nn <- rowsum(1 * !is.na(B),          doy[keep], reorder = TRUE)
  DM <- matrix(NA_real_, 366, length(VARS), dimnames = list(NULL, VARS))
  DM[as.integer(rownames(M)), ] <- M / Nn
  NORM <- vapply(seq_along(VARS),
                 function(j) rowMeans(matrix(DM[, j][WIN], nrow = 366), na.rm = TRUE),
                 numeric(366))
  colnames(NORM) <- VARS

  idx <- rows_by_cell[[cid]]
  if (is.null(idx)) return(NULL)
  out <- matrix(NA_real_, length(idx), length(VARS), dimnames = list(NULL, ANML))
  for (k in seq_along(idx)) {
    i    <- idx[k]
    days <- seq(d$collection_start_date[i], d$collection_end_date[i], by = "day")
    oi   <- match(days, ts$date)
    oi   <- oi[!is.na(oi)]
    # The normal holds one row per day of year, so a collection window longer
    # than a year still counts each day of year once.
    ni   <- unique(yday(days))
    if (!length(oi) || !length(ni)) next
    out[k, ] <- colMeans(OBS[oi, , drop = FALSE],  na.rm = TRUE) -
                colMeans(NORM[ni, , drop = FALSE], na.rm = TRUE)
  }
  tibble(fieldid = d$fieldid[idx]) %>% bind_cols(as_tibble(out))
}

message("computing anomalies ...")
t0 <- Sys.time()
res <- map_dfr(intersect(unique(d$cell), have), anom_one_cell)
message("  ", nrow(res), " fieldids in ",
        round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

res <- d %>% dplyr::select(fieldid) %>% left_join(res, by = "fieldid")
saveRDS(res, file.path(data_dir, "GMTP_anomalies.rds"))
write_tsv(res, file.path(data_dir, "GMTP_anomalies.tsv"))

#---- Report ------------------------------------------------------------------

message("\n=== anomaly distributions ===")
print(as.data.frame(res %>% dplyr::select(ends_with("_anml")) %>%
  pivot_longer(everything(), names_to = "variable") %>%
  group_by(variable) %>%
  summarise(n = sum(!is.na(value)), NA_n = sum(is.na(value)),
            mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
            q05 = quantile(value, .05, na.rm = TRUE),
            q95 = quantile(value, .95, na.rm = TRUE), .groups = "drop")),
  row.names = FALSE, digits = 4)

message("\n=== does a FIXED baseline make the anomaly encode year? ===")
chk <- d %>% dplyr::select(fieldid, year) %>% left_join(res, by = "fieldid")
print(as.data.frame(map_dfr(paste0(VARS, "_anml"), function(v) {
  ok <- is.finite(chk[[v]])
  m <- lm(chk[[v]][ok] ~ chk$year[ok])
  tibble(variable = v, slope_per_year = coef(m)[2],
         p_value = summary(m)$coefficients[2, 4],
         sd_of_anomaly = sd(chk[[v]], na.rm = TRUE),
         drift_over_18yr_in_sd = 18 * coef(m)[2] / sd(chk[[v]], na.rm = TRUE))
})), row.names = FALSE, digits = 4)
message("A |drift_over_18yr_in_sd| well under ~0.5 means the trend is second-order")
message("relative to week-to-week weather, and the fixed baseline is safe to use.")
message("\nwrote GMTP_anomalies.rds / .tsv")
