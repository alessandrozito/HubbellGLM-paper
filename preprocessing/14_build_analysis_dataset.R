# Assemble the modelling dataset: one row per fieldid, no NA in any covariate
# the models use.
#
# Produces: GMTP_analysis_dataset.rds

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

library(tidyverse)
library(sf)
library(lubridate)

data_dir <- DATA
ras_dir <- RASTER
wwf_shp <- file.path(CACHE, "ecoregions", "dataverse_files",
                     "tnc_terr_ecoregions.shp")

#---- 1. Base sample ----------------------------------------------------------

d <- readRDS(file.path(data_dir, "GMTP_sample_covariates.rds"))
message("sample: ", nrow(d), " rows")
na0 <- colSums(is.na(d[, c("aet", "hfp", "realm", "zone")]))
message("  starting NA: ", paste(names(na0), na0, collapse = "  "))

#---- 2. Zone: take script 05's resolved classification -----------------------

# All THREE classifications are kept, not just the resolved one, so the choice
# stays auditable and a sensitivity check can refit on either source:
#   zone_kgc   the kgc package's 0.5 degree lookup - the published definition
#              needed filling (10_fill_zone_gaps.R)
# zone_kgc reproduces the published classification.
kg <- read_tsv(file.path(data_dir, "kg_comparison_gmtp.tsv"),
               col_types = cols(.default = col_character())) %>%
  dplyr::select(fieldid,
                zone_kgc    = zone,
                kg_code_kgc = kg_code,
                kg_filled_km)
d <- d %>%
  dplyr::select(-any_of(c("zone", "kg_code"))) %>%
  left_join(kg, by = "fieldid") %>%
  mutate(zone         = zone_kgc,
         kg_code      = kg_code_kgc,
         kg_filled_km = as.numeric(kg_filled_km))
message("  zone NA: ", sum(is.na(d$zone)),
        "  (", sum(d$kg_filled_km > 0, na.rm = TRUE),
        " filled from the 1 km raster)")

#---- 3. Realm: intersect, then nearest polygon for the misses ----------------

pts <- d %>% distinct(Longitude, Latitude)
sfp <- st_as_sf(pts, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
wwf <- st_make_valid(st_read(wwf_shp, quiet = TRUE))

hit <- st_join(sfp, wwf["WWF_REALM2"], join = st_intersects) %>%
  st_drop_geometry() %>% distinct(Longitude, Latitude, .keep_all = TRUE)
miss <- which(is.na(hit$WWF_REALM2))
message("  realm: ", length(miss), " of ", nrow(hit), " coordinates outside every polygon")
if (length(miss)) {
  near <- st_join(sfp[miss, ], wwf["WWF_REALM2"], join = st_nearest_feature) %>%
    st_drop_geometry()
  hit$WWF_REALM2[miss] <- near$WWF_REALM2
  # how far the snap moved, for the record
  idx <- st_nearest_feature(sfp[miss, ], wwf)
  hit$realm_snap_km <- 0
  hit$realm_snap_km[miss] <-
    as.numeric(st_distance(sfp[miss, ], wwf[idx, ], by_element = TRUE)) / 1000
  message("    snapped, median ", round(median(hit$realm_snap_km[miss]), 2),
          " km, max ", round(max(hit$realm_snap_km[miss]), 1), " km")
} else {
  hit$realm_snap_km <- 0
}
hit <- hit %>%
  mutate(realm_new = ifelse(WWF_REALM2 == "Indo-Malay", "Indomalayan", WWF_REALM2)) %>%
  dplyr::select(Longitude, Latitude, realm_new, realm_snap_km)

d <- d %>%
  dplyr::select(-realm) %>%
  left_join(hit, by = c("Longitude", "Latitude")) %>%
  rename(realm = realm_new)
message("  realm NA after snapping: ", sum(is.na(d$realm)))

#---- 4. aet and hfp: nearest valid grid cell for the gaps --------------------

d <- d %>% mutate(aet_src = if_else(is.na(aet), NA_character_, "observed"),
                  hfp_src = if_else(is.na(hfp), NA_character_, "observed"))
message("  aet NA ", sum(is.na(d$aet)), ", hfp NA ", sum(is.na(d$hfp)))

# The remainder are coordinates the HFP land mask simply does not cover. Take
# the nearest grid cell that does have a value in the same year - the same
# "closest valid classification" rule used for zone. The prediction grid files
# already hold hfp at 1.54 M points per year, so no raster is re-read.
gap <- d %>% filter(is.na(hfp)) %>% distinct(Longitude, Latitude, year)
if (nrow(gap)) {
  message("  filling ", nrow(gap), " coordinate-years from the nearest grid cell")
  no_layer <- integer(0)
  fills <- map_dfr(sort(unique(gap$year)), function(yy) {
    f <- file.path(ras_dir, sprintf("Raster_%d.rds", yy))
    g <- gap %>% filter(year == yy)
    if (!file.exists(f)) return(mutate(g, hfp_fill = NA_real_, hfp_fill_km = NA_real_))
    r <- readRDS(f) %>% filter(!is.na(hfp)) %>% dplyr::select(Longitude, Latitude, hfp)
    # A year past the last HFP release has the grid but not the layer, so there
    # is no valid cell anywhere to snap to and the gap stays a gap.
    if (!nrow(r)) {
      no_layer <<- c(no_layer, yy)
      return(mutate(g, hfp_fill = NA_real_, hfp_fill_km = NA_real_))
    }
    out <- map_dfr(seq_len(nrow(g)), function(i) {
      # great-circle distance to every valid cell; 1.5 M rows is fast enough
      # for the handful of points that reach here.
      dlat <- (r$Latitude  - g$Latitude[i])
      dlon <- (r$Longitude - g$Longitude[i]) * cos(g$Latitude[i] * pi / 180)
      dd   <- dlat^2 + dlon^2
      k    <- which.min(dd)
      tibble(hfp_fill = r$hfp[k], hfp_fill_km = sqrt(dd[k]) * 111.32)
    })
    bind_cols(g, out)
  })
  d <- d %>%
    left_join(fills, by = c("Longitude", "Latitude", "year")) %>%
    mutate(hfp_src = if_else(is.na(hfp) & !is.na(hfp_fill), "nearest_cell", hfp_src),
           hfp     = coalesce(hfp, hfp_fill)) %>%
    dplyr::select(-hfp_fill)
  if (length(no_layer))
    message("    no HFP layer for ", paste(no_layer, collapse = ", "),
            " - ", sum(gap$year %in% no_layer), " coordinate-years stay NA")
  if (any(!is.na(fills$hfp_fill_km)))
    message("    median snap ", round(median(fills$hfp_fill_km, na.rm = TRUE), 2),
            " km, max ", round(max(fills$hfp_fill_km, na.rm = TRUE), 1), " km")
}
message("  final: aet NA ", sum(is.na(d$aet)), ", hfp NA ", sum(is.na(d$hfp)))

#---- 5. Seasonal term --------------------------------------------------------
# Week of year, shifted by half a year in the southern hemisphere so that the
# Fourier term lines up with the local growing season rather than the calendar.
# Same definition as the published preprocessing/00_newdata.R, lines 124-131.

d <- d %>%
  mutate(collection_start_date = as.Date(collection_start_date),
         week  = isoweek(collection_start_date),
         month = month(collection_start_date),
         week  = case_when(week == 53 & month == 1  ~ 1,
                           week == 53 & month == 12 ~ 52,
                           TRUE ~ week),
         week_adj = case_when(Latitude < 0 ~ (week + 26) %% 52, TRUE ~ week),
         week_adj = if_else(week_adj == 0, 52, week_adj)) %>%
  dplyr::select(-month)

#---- 6. ERA5 and wetlands, if scripts 06/07 have produced them ---------------

era <- list.files(ras_dir, "^era5_sample_\\d{4}\\.rds$", full.names = TRUE)
if (length(era)) {
  e <- map_dfr(era, readRDS) %>% dplyr::select(-any_of("ID"))
  d <- left_join(d, e, by = c("Longitude", "Latitude", "year"))
  message("  joined ERA5 for ", length(era), " years; ",
          sum(!is.na(d$temperature_2m)), " of ", nrow(d), " rows have it")
  # Temperature residualised on latitude, as in 1_GMTP_create_dataset.R:144.
  ok <- !is.na(d$temperature_2m)
  if (any(ok)) {
    m <- lm(temperature_2m ~ poly(Latitude, 2, raw = TRUE), data = d[ok, ])
    d$resid_temperature_lat <- NA_real_
    d$resid_temperature_lat[ok] <- residuals(m)
  }
} else message("  no ERA5 files yet - wind_speed and resid_temperature_lat unavailable")

wet <- list.files(ras_dir, "^wetlands_sample_\\d{4}\\.rds$", full.names = TRUE)
if (length(wet)) {
  # Keep script 07's own label. Years outside the CCI record (2023-2026) get a
  # file written by carry-forward, so the year matches exactly here and the
  # nearest-year rule below would call them "observed" - silently erasing the
  # fact that they are the 2022 surface.
  w <- map_dfr(wet, readRDS) %>%
    dplyr::select(-any_of("ID")) %>%
    dplyr::rename(wet_year = year, wet_src_file = wetlands_src)
  yrs <- sort(unique(w$wet_year))
  d$.wy <- yrs[max.col(-abs(outer(d$year, yrs, "-")), ties.method = "first")]
  d <- d %>%
    left_join(w, by = c("Longitude", "Latitude", ".wy" = "wet_year")) %>%
    mutate(wetlands_src = if_else(year == .wy, wet_src_file,
                                  paste0("held_", .wy))) %>%
    dplyr::select(-.wy, -wet_src_file)
  message("  joined wetlands from ", length(yrs), " map year(s): ",
          paste(yrs, collapse = ", "))
  message("    ", sum(!is.na(d$wetlands)), " of ", nrow(d), " rows have a value")
  print(table(d$wetlands_src, useNA = "ifany"))
} else message("  no wetland files yet - `wetlands` unavailable")

#---- 6c. Weather anomalies, from 13_compute_anomalies.R -----------------------
anom_f <- file.path(data_dir, "GMTP_anomalies.rds")
if (file.exists(anom_f)) {
  d <- left_join(d, readRDS(anom_f), by = "fieldid")
  nm <- grep("_anml$", names(d), value = TRUE)
  message("  joined anomalies: ", paste(nm, collapse = ", "))
  message("    complete on all four: ", sum(complete.cases(d[, nm])), " of ", nrow(d))
} else message("  no anomaly file yet - run 08 then 15")

#---- 6b. Site-level island covariates -----------------------------------------
# Keyed on site_code, from 08_site_island.R: landmass_km2, dist_mainland_km,
# on_continent. Joined here rather than in each analysis script so there is one
# dataset every downstream script reads.
fp <- file.path(data_dir, "GMTP_site_island.rds")
if (file.exists(fp)) {
  add <- readRDS(fp) %>%
    dplyr::select(-any_of(c("Longitude", "Latitude", "events", "country_iso")))
  new_cols <- setdiff(names(add), c("site_code", names(d)))
  d <- dplyr::left_join(d, add[, c("site_code", new_cols)], by = "site_code")
  message("  joined island covariates: ", paste(new_cols, collapse = ", "),
          "  (", sum(!is.na(d[[new_cols[1]]])), "/", nrow(d), " rows)")
} else message("  GMTP_site_island.rds not found - skipped")

#---- 7. Report and write -----------------------------------------------------

core <- c("aet", "hfp", "realm", "zone", "week_adj", "Latitude", "n", "y",
          "collection_days")
message("\nNA in the core model covariates:")
print(colSums(is.na(d[, core])))
message("\ncomplete cases on the core covariates: ",
        sum(complete.cases(d[, core])), " of ", nrow(d))
message("provenance:")
print(table(aet = d$aet_src, useNA = "ifany"))
print(table(hfp = d$hfp_src, useNA = "ifany"))

saveRDS(d, file.path(data_dir, "GMTP_analysis_dataset.rds"))
write_tsv(d, file.path(data_dir, "GMTP_analysis_dataset.tsv"))
message("\nwrote GMTP_analysis_dataset.rds / .tsv  (", nrow(d), " rows, ",
        ncol(d), " columns)")
