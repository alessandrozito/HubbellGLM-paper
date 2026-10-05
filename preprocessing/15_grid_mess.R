# MESS (Elith et al. 2010): how far outside the training envelope each grid cell
# sits. Written into every annual raster, so the analysis can cut on it without
# recomputing or re-joining anything.
#
# Adds to rasters_by_year/Raster_YYYY_full.rds:
#   mess_global         MESS against the whole training set
#   mod_global          the variable that limits it
#   mess_zone           MESS against that zone's training events only
#   mod_zone            the variable that limits it
#   vars_out_of_range   how many of the MESS variables are individually
#                       outside the zone's training range in that cell
#
# MUST RUN AFTER 14_build_analysis_dataset.R, and is the last preprocessing
# step. MESS scores each grid cell against the TRAINING ENVELOPE, and the
# training envelope is GMTP_analysis_dataset.rds itself - which step 14 builds.
# One of the MESS variables, resid_temperature_lat, does not exist anywhere
# before step 14 either: it is a residual against the sample's latitude fit,
# derived at 14_build_analysis_dataset.R:160. So this cannot be moved earlier
# in the sequence, and rerunning step 14 invalidates the MESS columns written
# here - rerun this one after it, with FORCE=1.

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

suppressMessages(library(tidyverse))

data_dir <- DATA
ras_dir  <- RASTER
FORCE    <- nzchar(Sys.getenv("FORCE"))

# The continuous predictors of the fitted model. week_adj and collection_days
# are held fixed when predicting, so they carry no spatial extrapolation risk;
# realm and zone are categorical and handled by the per-zone split below.
VARS <- c("aet", "Latitude", "hfp", "wind_speed", "wetlands", "resid_temperature_lat")
if (nzchar(Sys.getenv("MESS_VARS")))
  VARS <- str_split(Sys.getenv("MESS_VARS"), ",")[[1]]
MESS_COLS <- c("mess_global", "mod_global", "mess_zone", "mod_zone",
               "vars_out_of_range")

#---- Training envelope -------------------------------------------------------

ds_f <- file.path(data_dir, "GMTP_analysis_dataset.rds")
if (!file.exists(ds_f))
  stop("no GMTP_analysis_dataset.rds - this step scores the grid against the\n",
       "  training envelope, so 14_build_analysis_dataset.R must run first.",
       call. = FALSE)

d <- readRDS(ds_f) %>%
  dplyr::filter(year <= YEAR_MAX) %>%
  dplyr::mutate(zone = .data[[ZONE]])
d <- d[complete.cases(d[, c(VARS, "zone")]), ]
message("training envelope: ", format(nrow(d), big.mark = ","), " events, ",
        length(VARS), " variables (", paste(VARS, collapse = ", "), ")")
print(as.data.frame(d %>% count(zone, name = "events")), row.names = FALSE)

# resid_temperature_lat is a residual against the SAMPLE's latitude fit, so the
# grid has to be residualised against that same fit rather than its own.
lat_fit <- if ("resid_temperature_lat" %in% VARS)
  lm(temperature_2m ~ poly(Latitude, 2, raw = TRUE), data = d) else NULL

#---- MESS --------------------------------------------------------------------

# Similarity of p to the training vector v, one predictor. Negative means the
# cell sits outside the training range, scaled by the span of that range.
mess_var <- function(p, v) {
  v <- sort(v[is.finite(v)]); n <- length(v)
  mn <- v[1]; mx <- v[n]; rng <- mx - mn
  if (!is.finite(rng) || rng == 0) return(rep(NA_real_, length(p)))
  # findInterval(left.open = TRUE) counts training values STRICTLY below p.
  f <- 100 * findInterval(p, v, left.open = TRUE) / n
  ifelse(f == 0,  100 * (p - mn) / rng,
  ifelse(f <= 50, 2 * f,
  ifelse(f < 100, 2 * (100 - f),
                  100 * (mx - p) / rng)))
}

# The full per-variable similarity matrix, not only its row minimum: MESS and
# MoD say how bad the worst variable is and which one it is, but not how MANY
# are outside the range at once - which is what separates "one odd covariate"
# from "this place is not like anywhere we sampled".
mess_S <- function(newdata, ref) {
  S <- vapply(VARS, function(v) mess_var(newdata[[v]], ref[[v]]),
              numeric(nrow(newdata)))
  if (is.null(dim(S)))
    S <- matrix(S, nrow = nrow(newdata), dimnames = list(NULL, VARS))
  S
}

#---- One raster at a time ----------------------------------------------------
# 1.5 M rows x ~25 columns is comfortable on its own but not 27 of them at once,
# so each year is read, extended and written before the next is touched.

files <- sort(list.files(ras_dir, "^Raster_[0-9]{4}_full\\.rds$", full.names = TRUE))
if (!length(files))
  stop("no Raster_YYYY_full.rds in ", ras_dir,
       " - run 11_merge_annual_covariates.R first", call. = FALSE)
message("\n", length(files), " annual rasters")

# The MESS columns describe one particular training envelope. If step 14 has
# rebuilt the dataset since they were written, they describe an envelope that no
# longer exists, and skipping on "already has MESS" would keep them.
stale <- files[file.mtime(files) < file.mtime(ds_f)]
if (!FORCE && length(stale))
  message("  NOTE: ", length(stale), " raster(s) predate the analysis dataset. Any ",
          "MESS\n        columns in them are against an older envelope - ",
          "rerun with FORCE=1.")

for (f in files) {
  yr <- as.integer(str_extract(basename(f), "[0-9]{4}"))
  R  <- readRDS(f)
  if (!FORCE && all(MESS_COLS %in% names(R))) {
    message("  ", yr, ": already has MESS - skipping (FORCE=1 to redo)")
    next
  }
  R <- dplyr::select(R, -any_of(MESS_COLS))
  R[MESS_COLS] <- list(NA_real_, NA_character_, NA_real_, NA_character_,
                       NA_integer_)

  needed <- unique(c(setdiff(VARS, "resid_temperature_lat"),
                     if ("resid_temperature_lat" %in% VARS) "temperature_2m"))
  absent <- setdiff(needed, names(R))
  if (length(absent)) {
    message("  ", yr, ": raster has no ", paste(absent, collapse = ", "),
            " - MESS left NA")
    saveRDS(R, f, compress = "gzip")
    next
  }

  if (!is.null(lat_fit) && !"resid_temperature_lat" %in% names(R))
    R$resid_temperature_lat <- R$temperature_2m -
      predict(lat_fit, newdata = R[, "Latitude", drop = FALSE])
  R$zone <- R[[ZONE]]

  ok <- complete.cases(R[, VARS]) & !is.na(R$zone)
  if (!any(ok)) {
    # The columns are present but empty - a year past the END of a product,
    # such as hfp after 2024.
    empty <- VARS[!vapply(VARS, function(v) any(is.finite(R[[v]])), logical(1))]
    message("  ", yr, ": no cell has every predictor",
            if (length(empty)) paste0(" (all NA: ", paste(empty, collapse = ", "), ")"),
            " - MESS left NA")
    saveRDS(R, f, compress = "gzip")
    next
  }

  Rok <- R[ok, ]
  Sg  <- mess_S(Rok, d)
  R$mess_global[ok] <- apply(Sg, 1, min)
  R$mod_global[ok]  <- VARS[apply(Sg, 1, which.min)]

  Sz <- matrix(NA_real_, nrow(Rok), length(VARS), dimnames = list(NULL, VARS))
  zs <- Rok$zone
  for (z in sort(unique(zs))) {
    ref <- d[d$zone == z, ]
    if (nrow(ref) < 2) {
      message("     zone ", z, ": only ", nrow(ref), " training events - skipped")
      next
    }
    i <- which(zs == z)
    Sz[i, ] <- mess_S(Rok[i, ], ref)
  }
  has_z <- rowSums(is.na(Sz)) < length(VARS)
  R$mess_zone[ok][has_z]        <- apply(Sz[has_z, , drop = FALSE], 1, min)
  R$mod_zone[ok][has_z]         <- VARS[apply(Sz[has_z, , drop = FALSE], 1, which.min)]
  R$vars_out_of_range[ok][has_z] <- rowSums(Sz[has_z, , drop = FALSE] < 0, na.rm = TRUE)

  saveRDS(R, f, compress = "gzip")
  message(sprintf("  %d: %s cells scored | extrapolating (MESS < 0) global %.1f%%, zone %.1f%%",
                  yr, format(sum(ok), big.mark = ","),
                  100 * mean(R$mess_global < 0, na.rm = TRUE),
                  100 * mean(R$mess_zone   < 0, na.rm = TRUE)))
}

#---- Report on the reference year --------------------------------------------

f <- file.path(ras_dir, sprintf("Raster_%d_full.rds", REF_YEAR))
if (file.exists(f)) {
  R <- readRDS(f) %>% dplyr::mutate(zone = .data[[ZONE]])
  message("\n=== ", REF_YEAR, ": cumulative % of each zone below each MESS cut ===")
  print(as.data.frame(R %>% dplyr::filter(!is.na(mess_zone)) %>%
    group_by(zone) %>%
    summarise(cells = n(),
              `<0`   = round(100 * mean(mess_zone < 0),   1),
              `<-5`  = round(100 * mean(mess_zone < -5),  1),
              `<-15` = round(100 * mean(mess_zone < -15), 1),
              `<-25` = round(100 * mean(mess_zone < -25), 1),
              `<-50` = round(100 * mean(mess_zone < -50), 1),
              .groups = "drop")), row.names = FALSE)

  message("\n=== ", REF_YEAR, ": which variable binds, where MESS < 0 ===")
  print(as.data.frame(R %>% dplyr::filter(!is.na(mess_zone), mess_zone < 0) %>%
    count(zone, mod_zone, name = "n") %>% group_by(zone) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    slice_max(n, n = 3) %>% ungroup()), row.names = FALSE)

  message("\n=== ", REF_YEAR, ": how many variables are out of range at once ===")
  print(as.data.frame(R %>% dplyr::filter(!is.na(vars_out_of_range)) %>%
    count(zone, vars_out_of_range) %>% group_by(zone) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    dplyr::select(-n) %>%
    pivot_wider(names_from = vars_out_of_range, values_from = pct,
                names_prefix = "off_", values_fill = 0) %>% ungroup()),
    row.names = FALSE)
}

message("\nwrote ", paste(MESS_COLS, collapse = ", "), " into ", length(files),
        " annual rasters")
