# Landmass area and distance to the nearest continent for every site.
#
# Produces: GMTP_site_island.rds

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

suppressMessages({library(tidyverse); library(sf); library(rnaturalearth)})
sf_use_s2(TRUE)

data_dir <- DATA
MIN_CONT_KM2 <- as.numeric(Sys.getenv("MIN_CONT_KM2", "5e6"))
SNAP_KM      <- as.numeric(Sys.getenv("SNAP_KM", "25"))

d <- readRDS(file.path(data_dir, "GMTP_sample_clean.rds"))
sites <- d %>% group_by(site_code) %>%
  summarise(Longitude = mean(Longitude), Latitude = mean(Latitude),
            events = n(),
            country_iso = { t <- sort(table(country_iso), decreasing = TRUE)
                            if (length(t)) names(t)[1] else NA_character_ },
            .groups = "drop")
message("sites: ", nrow(sites))

# 10m scale, because Cocos and the Caribbean sites do not exist at 50m.
land <- ne_download(scale = 10, type = "land", category = "physical", returnclass = "sf") %>%
  st_make_valid() %>% st_cast("POLYGON", warn = FALSE)
land$km2 <- as.numeric(st_area(land)) / 1e6
land <- land[order(-land$km2), ]
land$landmass_id <- seq_len(nrow(land))
message("landmass polygons: ", nrow(land),
        "  largest: ", paste(round(head(land$km2, 5) / 1e6, 1), collapse = ", "), " M km2")

cont <- land %>% dplyr::filter(km2 >= MIN_CONT_KM2)
message("treated as continents: ", nrow(cont), " polygons (>= ",
        format(MIN_CONT_KM2, big.mark = ",", scientific = FALSE), " km2)")

pts <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326)

# Which landmass is each site on? Coastal coordinates can fall just offshore of
# a 10m coastline, so a miss snaps to the nearest polygon within SNAP_KM.
hit <- st_join(pts, land[, c("landmass_id", "km2")], join = st_within) %>% st_drop_geometry()
miss <- which(is.na(hit$landmass_id))
if (length(miss)) {
  nr <- st_nearest_feature(pts[miss, ], land)
  dk <- as.numeric(st_distance(pts[miss, ], land[nr, ], by_element = TRUE)) / 1000
  hit$landmass_id[miss] <- ifelse(dk <= SNAP_KM, land$landmass_id[nr], NA_integer_)
  hit$km2[miss]         <- ifelse(dk <= SNAP_KM, land$km2[nr], NA_real_)
  message("  ", length(miss), " sites not inside a polygon; ",
          sum(dk <= SNAP_KM), " snapped within ", SNAP_KM, " km, ",
          sum(dk > SNAP_KM), " left unresolved (max ", round(max(dk)), " km)")
}

# Distance to the nearest continent; zero if the site is on one.
dmain <- as.numeric(st_distance(pts, st_union(cont))) / 1000

isl <- sites %>%
  mutate(landmass_km2 = hit$km2,
         dist_mainland_km = dmain,
         on_continent = landmass_km2 >= MIN_CONT_KM2 & !is.na(landmass_km2))
write_tsv(isl, file.path(data_dir, "GMTP_site_island.tsv"))
saveRDS(isl, file.path(data_dir, "GMTP_site_island.rds"))

#---- Report ------------------------------------------------------------------
message("\n=== sites and events by landmass class ===")
print(as.data.frame(isl %>%
  mutate(class = case_when(
    on_continent ~ "continent",
    landmass_km2 >= 1e5 ~ "large island (>=100k km2)",
    landmass_km2 >= 1e4 ~ "medium island (10k-100k)",
    landmass_km2 >= 1e3 ~ "small island (1k-10k)",
    TRUE ~ "very small island (<1k km2)")) %>%
  group_by(class) %>%
  summarise(sites = n(), events = sum(events),
            med_dist_km = round(median(dist_mainland_km)), .groups = "drop")),
  row.names = FALSE)

message("\n=== where the places you asked about land ===")
NAMED <- c("IS","GL","CR","GB","JP","NZ","MG","PH","ID","BS","BES","SR","CU","LK","TW","PG")
print(as.data.frame(isl %>% dplyr::filter(country_iso %in% NAMED) %>%
  group_by(country_iso) %>%
  summarise(sites = n(), events = sum(events),
            landmass_km2 = round(median(landmass_km2)),
            dist_km = round(median(dist_mainland_km)),
            on_continent = any(on_continent), .groups = "drop") %>%
  arrange(landmass_km2)), row.names = FALSE)

message("\n=== the Cocos sites specifically ===")
print(as.data.frame(isl %>% dplyr::filter(grepl("ICOCO", site_code)) %>%
  transmute(site_code, events, landmass_km2 = round(landmass_km2, 1),
            dist_km = round(dist_mainland_km))), row.names = FALSE)

message("\n=== what each candidate rule would drop ===")
print(as.data.frame(map_dfr(list(
    list(l = "dist > 100 km",                 k = isl$dist_mainland_km > 100),
    list(l = "dist > 200 km",                 k = isl$dist_mainland_km > 200),
    list(l = "not on a continent",            k = !isl$on_continent),
    list(l = "area < 10,000 km2",             k = isl$landmass_km2 < 1e4),
    list(l = "area < 1,000 km2",              k = isl$landmass_km2 < 1e3),
    list(l = "area < 10,000 AND dist > 100",  k = isl$landmass_km2 < 1e4 & isl$dist_mainland_km > 100)),
  function(r) tibble(rule = r$l, sites = sum(r$k, na.rm = TRUE),
                     events = sum(isl$events[which(r$k)]),
                     pct_events = round(100 * sum(isl$events[which(r$k)]) / sum(isl$events), 2),
                     drops_Iceland = any(isl$country_iso[which(r$k)] == "IS"),
                     drops_Cocos = any(grepl("ICOCO", isl$site_code[which(r$k)]))))),
  row.names = FALSE)
message("\nwrote GMTP_site_island.{tsv,rds}")
