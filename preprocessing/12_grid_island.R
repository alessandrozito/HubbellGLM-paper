# Landmass area and distance to the nearest continent for every grid cell.
#
# Produces: grid_island_covariates.rds, and the same three columns written
#           into every rasters_by_year/Raster_YYYY_full.rds

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

suppressMessages({library(tidyverse); library(sf)})
source(file.path(file.path(REPO, "preprocessing"),
                 "island_covariates.R"))

data_dir <- DATA
MIN_CONT_KM2 <- as.numeric(Sys.getenv("MIN_CONT_KM2", "5e6"))
REF_YEAR     <- as.integer(Sys.getenv("REF_YEAR", "2024"))
CORES        <- as.integer(Sys.getenv("CORES", "23"))

g <- readRDS(file.path(data_dir, sprintf("rasters_by_year/Raster_%d_full.rds", REF_YEAR)))
message("grid: ", format(nrow(g), big.mark = ","), " cells")

ISL   <- c("landmass_km2", "dist_mainland_km", "on_continent")
isl_f <- file.path(data_dir, "grid_island_covariates.rds")
FORCE <- nzchar(Sys.getenv("FORCE"))

# The point-in-polygon pass costs minutes on every core, so a rerun that only
# needs the rasters rewritten reuses the stored result.
reuse <- !FORCE && file.exists(isl_f)
if (reuse) {
  out <- readRDS(isl_f)
  reuse <- nrow(out) == nrow(g) &&
           identical(out$Longitude, g$Longitude) &&
           identical(out$Latitude,  g$Latitude)
  message(if (reuse) "  reusing grid_island_covariates.rds"
          else "  stored grid_island_covariates.rds does not match the grid - recomputing")
}

if (!reuse) {

isl <- make_island_fn(MIN_CONT_KM2 = MIN_CONT_KM2, SNAP_KM = Inf)

# mclapply forks, so `land` and the unioned continents are shared copy-on-write
# rather than rebuilt CORES times. Both are plain R objects - sf geometries are
# lists, not external pointers - so they survive the fork intact.
# CORES * 4 chunks rather than CORES, so one slow chunk cannot leave the rest
# of the cores idle waiting on it.
t0 <- Sys.time()
ix <- split(seq_len(nrow(g)),
            ceiling(seq_len(nrow(g)) / ceiling(nrow(g) / (CORES * 4))))
res <- parallel::mclapply(ix, function(k) isl(g$Longitude[k], g$Latitude[k]),
                          mc.cores = CORES)
bad <- which(!vapply(res, is.data.frame, logical(1)))
if (length(bad))
  stop("chunk(s) failed: ", paste(bad, collapse = ", "), "\n",
       paste(utils::head(unlist(res[bad])), collapse = "\n"))
sn <- unlist(lapply(res, function(z) attr(z, "snap_km")), use.names = FALSE)
if (length(sn))
  message("  ", format(length(sn), big.mark = ","), " cells inside no land ",
          "polygon, snapped to nearest (median ", round(median(sn), 1),
          " km, max ", round(max(sn), 1), " km)")
out <- as_tibble(bind_rows(res))
stopifnot(nrow(out) == nrow(g),
          identical(out$Longitude, g$Longitude),
          identical(out$Latitude,  g$Latitude))
message("  computed in ", round(difftime(Sys.time(), t0, units = "mins"), 2),
        " min on ", CORES, " cores")
stopifnot(!anyNA(out$landmass_km2), !anyNA(out$dist_mainland_km))
saveRDS(out, isl_f)

}

#---- Fold into the annual rasters --------------------------------------------
message("\nwriting the island columns into the annual rasters ...")
for (f in sort(list.files(RASTER, "^Raster_[0-9]{4}_full\\.rds$", full.names = TRUE))) {
  g_y <- readRDS(f)
  stopifnot(identical(g_y$Longitude, out$Longitude),
            identical(g_y$Latitude,  out$Latitude))
  g_y <- dplyr::bind_cols(dplyr::select(g_y, -any_of(ISL)), out[, ISL])
  saveRDS(g_y, f, compress = "gzip")
  message("  ", basename(f), "  (", ncol(g_y), " columns)")
}

message("\n=== grid cells by landmass class ===")
print(as.data.frame(out %>%
  mutate(class = case_when(
    on_continent        ~ "continent",
    landmass_km2 >= 1e5 ~ "large island (>=100k km2)",
    landmass_km2 >= 1e4 ~ "medium island (10k-100k)",
    landmass_km2 >= 1e3 ~ "small island (1k-10k)",
    TRUE                ~ "very small island (<1k km2)")) %>%
  group_by(class) %>%
  summarise(cells = n(), pct = round(100*n()/nrow(out), 2),
            med_dist_km = round(median(dist_mainland_km)), .groups = "drop")),
  row.names = FALSE)

message("\n=== range check against the SAMPLE (this is what MESS uses) ===")
s <- readRDS(file.path(data_dir, "GMTP_site_island.rds"))
for (v in c("landmass_km2", "dist_mainland_km")) {
  rs <- range(s[[v]], na.rm = TRUE); rg <- range(out[[v]], na.rm = TRUE)
  message(sprintf("  %-18s sample [%.4g, %.4g]  grid [%.4g, %.4g]  outside: %.2f%% below, %.2f%% above",
    v, rs[1], rs[2], rg[1], rg[2],
    100*mean(out[[v]] < rs[1]), 100*mean(out[[v]] > rs[2])))
}
message("\nwrote grid_island_covariates.rds and updated the annual rasters")
