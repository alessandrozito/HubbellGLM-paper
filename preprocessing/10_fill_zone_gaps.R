# Fill the Koppen-Geiger zone where the kgc 0.5 degree lookup returns nothing,
# using the Beck et al. 1 km raster.
#
# Produces: kg_comparison_gmtp.tsv

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

library(tidyverse)
library(terra)

data_dir <- DATA
kg_path <- file.path(CACHE, "koppen", "koppen_geiger_0p00833333.tif")
gmtp_in <- file.path(data_dir, "GMTP_sample_covariates.rds")

# The kgc labels this script repairs are time-invariant, so any one of script
# 09's annual grids carries them. GRID itself holds coordinates only.
grid_files <- sort(list.files(RASTER, "^Raster_[0-9]{4}\\.rds$", full.names = TRUE))
if (!length(grid_files))
  stop("no Raster_YYYY.rds in ", RASTER,
       " - run 09_merge_site_covariates.R first", call. = FALSE)
grid_in <- grid_files[1]

kg_legend <- tribble(
  ~KG_ID, ~KG_Code,
  1,"Af", 2,"Am", 3,"Aw", 4,"BWh", 5,"BWk", 6,"BSh", 7,"BSk",
  8,"Csa", 9,"Csb", 10,"Csc", 11,"Cwa", 12,"Cwb", 13,"Cwc",
  14,"Cfa", 15,"Cfb", 16,"Cfc", 17,"Dsa", 18,"Dsb", 19,"Dsc", 20,"Dsd",
  21,"Dwa", 22,"Dwb", 23,"Dwc", 24,"Dwd", 25,"Dfa", 26,"Dfb", 27,"Dfc",
  28,"Dfd", 29,"ET", 30,"EF"
)

kopp_to_zone <- function(codes) {
  z <- c("A" = "Tropical", "B" = "Dry", "C" = "Temperate",
         "D" = "Continental", "E" = "Polar")
  ifelse(is.na(codes) | codes == "Climate Zone info missing",
         NA_character_, unname(z[substr(codes, 1, 1)]))
}

kg <- rast(kg_path)

# MAJORITY (modal) lookup over a square window, then the legend. 0 is ocean.
#
# WHY NOT method = "simple"
#   A single 1 km cell is too fine a target for a trap coordinate. At Los
#   Quetzales, Costa Rica (9.614, -83.819) the point lands on an ET (tundra)
#   sliver, which the first-letter mapping turns into "Polar" - but the site is
#   at ~2,600 m in Talamancan montane forest and 76% of the surrounding cells
#   are Cwb/Cfb. That single sliver put 141 collection events into the Polar
#   zone. The modal rule returns Cwb -> Temperate there, while genuine paramo
#   at Paluguillo, Ecuador (-0.306, -78.231, ~3,500 m) stays ET -> Polar at
#   every window size. So the majority rule separates the artefact from the
#   real high-altitude E, which is the point: E still maps to Polar, it just
#   has to be the majority of the neighbourhood rather than one cell.
#
#   Ties are broken toward the cell the point actually falls in, by giving it
#   half a vote more than any neighbour.
#
#   MODAL_K sets the half-width in cells (default 3, about 2.8 km).
MODAL_K <- as.integer(Sys.getenv("MODAL_K", "3"))

label_1km <- function(d, method = "modal", k = MODAL_K) {
  xy <- as.matrix(d[, c("Longitude", "Latitude")])
  if (!identical(method, "modal")) {
    # xy is a matrix, so extract() returns the layer only - no ID column.
    v  <- terra::extract(kg, xy, method = method)[, 1]
    id <- if (method == "bilinear") round(v) else v
  } else {
    cs  <- terra::res(kg)
    off <- expand.grid(dx = -k:k, dy = -k:k)
    cnt <- matrix(0, nrow(xy), 30)
    for (i in seq_len(nrow(off))) {
      v <- terra::extract(kg, cbind(xy[, 1] + off$dx[i] * cs[1],
                                    xy[, 2] + off$dy[i] * cs[2]))[, 1]
      ok <- which(!is.na(v) & v >= 1 & v <= 30)
      if (!length(ok)) next
      w <- if (off$dx[i] == 0 && off$dy[i] == 0) 1.5 else 1   # centre tie-break
      cnt[cbind(ok, v[ok])] <- cnt[cbind(ok, v[ok])] + w
    }
    id <- max.col(cnt, ties.method = "first")
    id[rowSums(cnt) == 0] <- NA
  }
  id[is.na(id) | id < 1 | id > 30] <- NA
  tibble(KG_ID = id) %>%
    left_join(kg_legend, by = "KG_ID") %>%
    mutate(zone_1km = kopp_to_zone(KG_Code)) %>%
    dplyr::select(KG_ID, KG_Code_1km = KG_Code, zone_1km)
}

#---- Gap filling -------------------------------------------------------------
# The traps are on land, so a point that lands on ocean/no-data is a
# geolocation or rounding artifact rather than a real absence of climate. Those
# points take the class of the nearest valid cell.
#
# The search walks outward one ring at a time in Chebyshev radius. Any cell
# outside a square of half-width k has max(|dx|, |dy|) > k and so lies further
# than k * cellsize; once a hit is found within k * cellsize it is provably the
# nearest and the search stops. The distance moved is recorded for every filled
# point so the substitution stays auditable.

fill_nearest_raster <- function(r, lon, lat, max_k = 25) {
  n <- length(lon)
  val <- rep(NA_real_, n); dist <- rep(NA_real_, n)
  todo <- seq_len(n); cs <- res(r)
  for (k in seq_len(max_k)) {
    if (!length(todo)) break
    g  <- expand.grid(dx = -k:k, dy = -k:k)
    g  <- g[pmax(abs(g$dx), abs(g$dy)) == k, ]
    gd <- sqrt((g$dx * cs[1])^2 + (g$dy * cs[2])^2)
    o  <- order(gd); g <- g[o, ]; gd <- gd[o]
    bv <- rep(NA_real_, length(todo)); bd <- rep(Inf, length(todo))
    for (i in seq_len(nrow(g))) {
      if (all(!is.na(bv) & bd <= gd[i])) break
      v <- terra::extract(r, cbind(lon[todo] + g$dx[i] * cs[1],
                                   lat[todo] + g$dy[i] * cs[2]))[, 1]
      hit <- !is.na(v) & v >= 1 & v <= 30 & gd[i] < bd
      bv[hit] <- v[hit]; bd[hit] <- gd[i]
    }
    got <- !is.na(bv)
    val[todo[got]] <- bv[got]; dist[todo[got]] <- bd[got]
    todo <- todo[!got]
  }
  list(value = val, dist_km = dist * 111, unresolved = todo)
}

# Same idea on the kgc 0.5-degree lookup, which is a regular grid keyed on
# rounded coordinates.
fill_nearest_kgc <- function(cz, lon, lat, step = 0.5, max_k = 12) {
  key <- paste(cz$Lat, cz$Lon)
  map <- setNames(as.character(cz$Cls), key)
  rlat <- kgc::RoundCoordinates(lat); rlon <- kgc::RoundCoordinates(lon)
  n <- length(lon)
  val <- rep(NA_character_, n); dist <- rep(NA_real_, n)
  todo <- seq_len(n)
  for (k in seq_len(max_k)) {
    if (!length(todo)) break
    g  <- expand.grid(dx = -k:k, dy = -k:k)
    g  <- g[pmax(abs(g$dx), abs(g$dy)) == k, ]
    gd <- sqrt((g$dx * step)^2 + (g$dy * step)^2)
    o  <- order(gd); g <- g[o, ]; gd <- gd[o]
    bv <- rep(NA_character_, length(todo)); bd <- rep(Inf, length(todo))
    for (i in seq_len(nrow(g))) {
      if (all(!is.na(bv) & bd <= gd[i])) break
      v <- unname(map[paste(rlat[todo] + g$dy[i] * step,
                            rlon[todo] + g$dx[i] * step)])
      hit <- !is.na(v) & v != "Climate Zone info missing" & gd[i] < bd
      bv[hit] <- v[hit]; bd[hit] <- gd[i]
    }
    got <- !is.na(bv)
    val[todo[got]] <- bv[got]; dist[todo[got]] <- bd[got]
    todo <- todo[!got]
  }
  list(value = val, dist_km = dist * 111, unresolved = todo)
}

report <- function(d, what) {
  message("\n", strrep("=", 64), "\n", what, " (n = ", nrow(d), ")\n", strrep("=", 64))
  both <- d %>% filter(!is.na(zone), !is.na(zone_1km))
  message("  usable on both labels: ", nrow(both),
          "   kgc NA: ", sum(is.na(d$zone)),
          "   1km NA: ", sum(is.na(d$zone_1km)))

  agr <- mean(both$zone == both$zone_1km)
  message("  ZONE agreement (5 classes): ", sprintf("%.1f%%", 100 * agr))

  bc <- both %>% filter(!is.na(kg_code), !is.na(KG_Code_1km))
  message("  FULL CODE agreement (30 classes): ",
          sprintf("%.1f%%", 100 * mean(bc$kg_code == bc$KG_Code_1km)),
          "   (n = ", nrow(bc), ")")

  message("\n  zone confusion (rows = kgc 0.5 deg, cols = 1 km raster):")
  tb <- table(kgc = both$zone, km1 = both$zone_1km)
  print(tb)
  message("\n  per-class recall (share of each kgc zone the 1 km map agrees with):")
  rc <- round(100 * diag(tb) / rowSums(tb), 1)
  print(rc)
  invisible(both)
}

#---- Attach both labels and gap-fill ----------------------------------------
e <- new.env(parent = emptyenv())
utils::data("climatezones", package = "kgc", envir = e)
cz <- get("climatezones", envir = e)

prepare <- function(d, what) {
  d <- bind_cols(d, label_1km(d))
  d$kg_filled_km  <- NA_real_
  d$km1_filled_km <- NA_real_

  # 1 km raster gaps
  i <- which(is.na(d$KG_ID))
  if (length(i)) {
    message("  ", what, ": filling ", length(i), " points with no 1 km class ...")
    f <- fill_nearest_raster(kg, d$Longitude[i], d$Latitude[i])
    d$KG_ID[i]         <- f$value
    d$km1_filled_km[i] <- f$dist_km
    if (length(f$unresolved)) {
      message("    ", length(f$unresolved),
              " still unresolved beyond the search radius (left NA)")
    }
    d <- d %>%
      dplyr::select(-KG_Code_1km, -zone_1km) %>%
      left_join(kg_legend, by = "KG_ID") %>%
      mutate(KG_Code_1km = KG_Code, zone_1km = kopp_to_zone(KG_Code)) %>%
      dplyr::select(-KG_Code)
  }

  # kgc gaps
  j <- which(is.na(d$zone))
  if (length(j)) {
    message("  ", what, ": filling ", length(j), " points with no kgc class ...")
    f <- fill_nearest_kgc(cz, d$Longitude[j], d$Latitude[j])
    d$kg_code[j]      <- f$value
    d$zone[j]         <- kopp_to_zone(f$value)
    d$kg_filled_km[j] <- f$dist_km
    if (length(f$unresolved)) {
      message("    ", length(f$unresolved), " still unresolved (left NA)")
    }
  }

  fk <- d$km1_filled_km[!is.na(d$km1_filled_km)]
  gk <- d$kg_filled_km[!is.na(d$kg_filled_km)]
  if (length(fk)) message("    1 km fill distance (km): median ",
                          round(median(fk), 2), ", max ", round(max(fk), 1))
  if (length(gk)) message("    kgc  fill distance (km): median ",
                          round(median(gk), 2), ", max ", round(max(gk), 1))

  # "Closest label wins": whichever source needed the least extrapolation.
  # A label measured in situ (distance 0) always beats a borrowed one, and the
  # two are only compared when both had to reach. This matters for islands:
  # the 0.5-degree kgc grid does not resolve them, so it borrows a mainland
  # class from as far as 477 km, while the 1 km raster has a class on the
  # island itself. Ties go to the 1 km layer as the finer of the two.
  d <- d %>%
    mutate(
      .d_kgc = ifelse(is.na(zone),     Inf, coalesce(kg_filled_km,  0)),
      .d_1km = ifelse(is.na(zone_1km), Inf, coalesce(km1_filled_km, 0)),
      zone_best_src = case_when(
        is.infinite(.d_kgc) & is.infinite(.d_1km) ~ NA_character_,
        .d_1km <= .d_kgc                          ~ "1km",
        TRUE                                      ~ "kgc"),
      zone_best = case_when(
        is.na(zone_best_src)   ~ NA_character_,
        zone_best_src == "1km" ~ zone_1km,
        TRUE                   ~ zone),
      zone_best_dist_km = pmin(.d_kgc, .d_1km)
    ) %>%
    dplyr::select(-.d_kgc, -.d_1km)

  message("    zone_best source: ",
          paste(names(table(d$zone_best_src)), table(d$zone_best_src),
                sep = "=", collapse = ", "),
          "   NA: ", sum(is.na(d$zone_best)))
  moved <- d %>% filter(zone_best_dist_km > 0)
  message("    labels taken from a neighbouring cell: ", nrow(moved),
          if (nrow(moved)) paste0(" (max ", round(max(moved$zone_best_dist_km), 1),
                                  " km)") else "")
  d
}

#---- Prediction grid ---------------------------------------------------------
message("Loading the prediction grid from ", basename(grid_in), " ...")
grid <- readRDS(grid_in) %>% dplyr::select(Longitude, Latitude, kg_code, zone)
grid_raw <- bind_cols(grid, label_1km(grid))
report(grid_raw, "PREDICTION GRID - before gap filling")
grid <- prepare(grid, "grid")
report(grid, "PREDICTION GRID - after gap filling")

#---- GMTP sample -------------------------------------------------------------
message("\nLoading the GMTP sample ...")
gmtp0 <- readRDS(gmtp_in) %>%
  dplyr::select(fieldid, site_code, Longitude, Latitude, kg_code, zone)
report(bind_cols(gmtp0, label_1km(gmtp0)), "GMTP SAMPLE - before gap filling")
gmtp <- prepare(gmtp0, "sample")
report(gmtp, "GMTP SAMPLE - after gap filling")

#---- What bilinear does to a categorical raster ------------------------------
message("\n", strrep("=", 64))
message("Effect of method = 'bilinear' on these class IDs")
message(strrep("=", 64))
bil <- label_1km(gmtp, method = "bilinear")
cmp <- tibble(simple = gmtp$KG_ID, bilinear = bil$KG_ID,
              zs = gmtp$zone_1km, zb = bil$zone_1km)
ok <- cmp %>% filter(!is.na(simple), !is.na(bilinear))
message("  KG_ID differs from nearest-neighbour at ",
        sprintf("%.1f%%", 100 * mean(ok$simple != ok$bilinear)),
        " of the ", nrow(ok), " sample points")
okz <- ok %>% filter(!is.na(zs), !is.na(zb))
message("  and the 5-class zone differs at ",
        sprintf("%.1f%%", 100 * mean(okz$zs != okz$zb)), " of them")
message("  worst case: a point can be assigned a class neither neighbour has.")

#---- Save --------------------------------------------------------------------
write_tsv(gmtp, file.path(data_dir, "kg_comparison_gmtp.tsv"))
saveRDS(grid, file.path(data_dir, "kg_comparison_grid.rds"), compress = "gzip")
message("\nWrote kg_comparison_gmtp.tsv and kg_comparison_grid.rds to ", data_dir)
