# island_covariates(lon, lat) -> landmass_km2, dist_mainland_km, on_continent
#
# Exact point-in-polygon and exact great-circle distance (s2), same definitions
# as 08_site_island.R. Land is loaded once by the factory, so repeat calls
# on the returned function are cheap.
#
#   source("preprocessing/island_covariates.R")
#   isl <- make_island_fn()
#   isl(-87.08, 5.53)                      # one point
#   isl(c(-87.08, -21.54), c(5.53, 64.52)) # or vectors

suppressMessages({library(sf); library(rnaturalearth)})
sf_use_s2(TRUE)

# SNAP_KM: a point falling inside no land polygon is snapped to the nearest one
#   within this many km, else NA. 25 is right for SITES - a trap 100 km out to
#   sea is a data error worth surfacing. For the prediction GRID pass Inf: every
#   cell is already land by the ERA5 mask, so the only question is WHICH
#   landmass it belongs to, and an NA there just breaks the model matrix.
make_island_fn <- function(MIN_CONT_KM2 = 5e6, SNAP_KM = 25, scale = 10) {
  land <- ne_download(scale = scale, type = "land", category = "physical",
                      returnclass = "sf") |>
    st_make_valid() |> st_cast("POLYGON", warn = FALSE)
  land$km2 <- as.numeric(st_area(land)) / 1e6
  land <- land[order(-land$km2), ]
  cont <- st_union(land[land$km2 >= MIN_CONT_KM2, ])
  message("land polygons: ", nrow(land), " | continents: ",
          sum(land$km2 >= MIN_CONT_KM2))

  function(lon, lat) {
    pts <- st_as_sf(data.frame(lon = lon, lat = lat),
                    coords = c("lon", "lat"), crs = 4326)

    # 1. which landmass contains the point -> its area
    i <- st_within(pts, land)
    km2 <- vapply(i, function(k) if (length(k)) land$km2[k[1]] else NA_real_, numeric(1))

    # 2. a coastal point can fall just offshore of the 10m coastline; snap it
    miss <- which(is.na(km2))
    if (length(miss)) {
      nr <- st_nearest_feature(pts[miss, ], land)
      dk <- as.numeric(st_distance(pts[miss, ], land[nr, ], by_element = TRUE)) / 1000
      km2[miss] <- ifelse(dk <= SNAP_KM, land$km2[nr], NA_real_)
      attr(km2, "snap_km") <- dk
    }

    # 3. great-circle distance to the nearest continent; 0 if on one
    dist_km <- as.numeric(st_distance(pts, cont)) / 1000

    out <- data.frame(Longitude = lon, Latitude = lat,
                      landmass_km2 = as.numeric(km2), dist_mainland_km = dist_km,
                      on_continent = !is.na(km2) & km2 >= MIN_CONT_KM2)
    attr(out, "snap_km") <- attr(km2, "snap_km")
    out
  }
}
