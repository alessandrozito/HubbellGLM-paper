# Create the final raster for out of sample prediction. We rely on the raster
# from 2014, since this is the most frequent year in the collection

#---- Load the packages
library(tidyverse)
library(HubbellGLM)
library(sf)
library(terra)
library("rnaturalearth")
library("rnaturalearthdata")
library(geosphere)
library(kgc)
library(spdep)

#==========================================
# Useful functions
#==========================================

# Function to extract the climatic zone from the KG codes.
kopp_to_zone <- function(koppen_codes) {
  # Mapping vector
  zone_names <- c(
    "A" = "Tropical",
    "B" = "Dry",
    "C" = "Temperate",
    "D" = "Continental",
    "E" = "Polar"
  )
  # Vectorized assignment
  result <- ifelse(
    koppen_codes == "Climate Zone info missing",
    NA,
    zone_names[substr(koppen_codes, 1, 1)]
  )
  return(result)
}


#---- Read-in the raster dataset and the GMTP data
df_Raster <- readRDS("~/HubbellGLM-paper/data/climate/Raster_2014.rds")
load("~/HubbellGLM-paper/data/data_GMTP_clean.rdata")

# Load the climatic zone dataset of the package kgc and merge it
data("climatezones", package = "kgc")

df_Raster <- df_Raster %>%
  dplyr::mutate(rndCoord.lon = kgc::RoundCoordinates(Longitude),
                rndCoord.lat = kgc::RoundCoordinates(Latitude)) %>%
  dplyr::left_join(climatezones, by = c("rndCoord.lat" = "Lat", "rndCoord.lon" = "Lon")) %>%
  dplyr::mutate(kg_code = as.character(Cls),
         zone = kopp_to_zone(kg_code)) %>%
  dplyr::select(-Cls)


#---- Merge with wwf ecoregions
wwf_ecoregions <- st_read("~/HubbellGLM-paper/data/climate/ecoregions/dataverse_files/tnc_terr_ecoregions.shp")
wwf_ecoregions <- st_make_valid(wwf_ecoregions)
my_points_sf <- st_as_sf(df_Raster %>%
                           dplyr::mutate(x = Longitude, y = Latitude) %>%
                           dplyr::select(x, y) %>%
                           distinct(),
                         coords = c("x", "y"), crs = 4326)
joined_data <- st_join(my_points_sf, wwf_ecoregions, join = st_intersects)
coords <- st_coordinates(joined_data)
joined_data$Longitude <- coords[, "X"]
joined_data$Latitude  <- coords[, "Y"]


df_Raster <- df_Raster %>%
  dplyr::left_join(joined_data %>% dplyr::select(Longitude, Latitude, WWF_REALM2),
            by = c("Longitude", "Latitude")) %>%
  dplyr::mutate(realm = case_when(WWF_REALM2 == "Indo-Malay"~ "Indomalayan",
                                  TRUE ~ WWF_REALM2)) %>%
  dplyr::select(-WWF_REALM2, -geometry)


#---- Merge with aet
dest_file <- "~/HubbellGLM-paper/data/climate/TerraClimate_aet_2014.nc"
if (!file.exists(dest_file)) {
  url <- "http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data/TerraClimate_aet_2014.nc"
    download.file(url, dest_file, mode = "wb")
}

aet_2014 <- rast(dest_file)
points_vect <- vect(df_Raster, geom = c("x", "y"), crs = "EPSG:4326")
aet_extracted <- extract(aet_2014, points_vect, method = "bilinear")
df_Raster$aet <- rowSums(aet_extracted[, -1], na.rm = TRUE) * 0.1 # <-- sum over the year and rescale

#---- Merge with residuals of the temperature and latitude
data_MAT <- dataset %>% dplyr::select(Latitude, temperature_2m_year_avg)
fit_MAT <- lm(temperature_2m_year_avg ~ poly(Latitude, 2, raw = TRUE), data = data_MAT)
resid_out <- df_Raster$temperature_2m - predict(fit_MAT, df_Raster)
df_Raster$resid_temperature_lat <- resid_out

saveRDS(df_Raster, file = "~/HubbellGLM-paper/data/climate/Raster_2014_merged.rds",
        compress = "gzip")

