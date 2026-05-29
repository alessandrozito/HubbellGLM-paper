# This file creates the final dataset that will be used for the analysis
# in all subsequent scripts. We also calculate the fraction of shared species
# between each two locations.

#---- Load packages
library(tidyverse)
library(sf)
library(lubridate)
library(kgc) # <-- To extract the Koppen-Geiger Climatic zones
library(vegan)
library("rnaturalearth")
library("rnaturalearthdata")

#============================================================
# Useful functions
#============================================================
# Fast function to calculate the Jaccard Index between two sites
# This function is mathematically equivalent to
#     1 - vegdist(x, method = "jaccard", binary = TRUE)
# in the R package "vegan". However, it is significanly faster for large samples
get_shared_species <- function(Xbins){
  # Find shared species
  Xshared <- tcrossprod(1 * (Xbins > 0))
  # Calculate the fraction
  n_species_samples <- diag(Xshared)
  vec1 <- rep(1, length(n_species_samples))
  Xres <- Xshared/(tcrossprod(n_species_samples, vec1) + tcrossprod(vec1, n_species_samples) - Xshared)
  return(Xres)
}

# Function to exctract the climatic zone from the KG codes.
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


#============================================================
# Merge all the datasets
#============================================================
#---- Load the data
data_fields <- readRDS("~/HubbellGLM-paper/data/data_GMTP_merged_covariates.rds.gzip")
abundance_matrix <- readRDS("~/HubbellGLM-paper/data/abundance_matrix.rds.gzip")

#---- Preprocess data
dataset <- data_fields %>%
  # Filter out locations with at least 20 arthropods, remove traps that are
  # combined, and damaged by a bear
  dplyr::filter(n >= 20,
                combined_traps == FALSE,
                trap_damaged == FALSE) %>%
  # Day of the year fourier transformation
  dplyr::mutate(Latitude = Latitude_clean,
                Longitude = Longitude_clean,
                # We look at annual AET and PET, so multiply by 12
                # We also multiply by 0.1 to change scale according to here: https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_TERRACLIMATE#bands
                aet = Terra_aet_year_avg * 0.1 * 12,
                pet = Terra_pet_year_avg * 0.1 * 12,
                # Day of the year
                yday = yday(collection_start_date_clean) - 1*(yday(collection_start_date_clean) == 366),
                yday_adj = case_when(Latitude_clean < 0 ~ (yday + 365/2) %% 365,
                                     TRUE ~ yday),
                # Month and week of the year
                month = month(collection_start_date_clean),
                week = isoweek(collection_start_date_clean),
                week = case_when(week == 53 & month == 1 ~ 1,
                                 week == 53 & month == 12 ~ 52,
                                 TRUE ~ week),
                week_adj = case_when(Latitude_clean < 0 ~ (week + 26) %% 52,
                                     TRUE ~ week),
                week_adj = case_when(week_adj == 0 ~ 52,
                                     TRUE ~ week_adj),
                # Year
                year = as.factor(year(collection_start_date_clean)),
                # Collection length
                collection_days = as.numeric(collection_end_date_clean - collection_start_date_clean))

#----
# Add Koppen-Geiger zones using the package kgc.
# There are 6 sites where kg_code is missing. We substitute it with "Tropical",
# given the location.
dataset <- dataset %>%
  dplyr::mutate(kg_code = LookupCZ(data.frame(rndCoord.lon = RoundCoordinates(Longitude),
                                              rndCoord.lat = RoundCoordinates(Latitude))),
                zone = kopp_to_zone(kg_code))

# Plot the location of the missing KG zone
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- world[world$name !=  "Antarctica", ]

dataset %>%
  dplyr::filter(is.na(zone)) %>%
  ggplot() +
  theme_bw() +
  geom_point(aes(x = Longitude, y = Latitude)) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.15)

# Substitute it with Tropical. These are only 6 points
dataset <- dataset %>%
  dplyr::mutate(zone = case_when(is.na(zone) ~ "Tropical",
                                 TRUE ~ zone))

# Merge with the wwf ecoregion from the Harvard dataverse. There is a minor
# mismatch between the realm indicated in the ecoregion, and the realm that is
# indicated in the dataset (19 points in Arabia that are labelled as Aftrotropic in dataverse,
# but Palearctic in the source data). We rely on the source data, which are also used
# to fill the NAs. We keep both to verify if the analysis is robust.
wwf_ecoregions <- st_read("~/HubbellGLM-paper/data/climate/ecoregions/dataverse_files/tnc_terr_ecoregions.shp")
wwf_ecoregions <- st_make_valid(wwf_ecoregions)
my_points_sf <- st_as_sf(dataset %>%
                           dplyr::mutate(x = Longitude, y = Latitude) %>%
                           dplyr::select(x, y) %>%
                           distinct(),
                         coords = c("x", "y"), crs = 4326)
joined_data <- st_join(my_points_sf, wwf_ecoregions, join = st_intersects)
coords <- st_coordinates(joined_data)
joined_data$Longitude <- coords[, "X"]
joined_data$Latitude  <- coords[, "Y"]

dataset <- dataset %>%
  dplyr::left_join(joined_data, by = c("Longitude", "Latitude")) %>%
  dplyr::mutate(realm_wwf = case_when(WWF_REALM2 == "Indo-Malay"~ "Indomalayan",
                             TRUE ~ WWF_REALM2),
                realm_wwf = case_when(is.na(realm_wwf) ~ realm,
                                      TRUE ~ realm_wwf))

# Calculate the residuals between the average annual temperature and a
# polynomial regression on the latitude
data_MAT <- dataset %>% dplyr::select(Latitude, temperature_2m_year_avg)
fit_MAT <- lm(temperature_2m_year_avg ~ poly(Latitude, 2, raw = TRUE), data = data_MAT)
dataset$resid_temperature_lat <- residuals(fit_MAT)

#---- Clean the column names, and take all variables needed for the analysis
dataset <- dataset %>%
  dplyr::mutate(aet = Terra_aet_year_avg,
                wind_speed = wind_speed_year_avg,
                site = site_clean,
                country = country_clean,
                collection_start_date = collection_start_date_clean,
                collection_end_date = collection_end_date_clean) %>%
  dplyr::select(fieldid, site, country, Latitude, Longitude, elevation,
                collection_start_date, collection_end_date, n, y, aet, hfp, zone,
                realm, realm_wwf, week_adj, collection_days, wind_speed_year_avg, wetlands,
                resid_temperature_lat, pet, vpd_year_avg, NDVI_year_max,
                temperature_2m_year_avg, total_precipitation_sum_year_avg, relative_humidity_year_avg,
                temperature_2m_anml, total_precipitation_sum_anml, relative_humidity_anml, wind_speed_anml)

#---- Calculate the fraction of shared species between each two collection.

# We first demonstrate that it is equivalent to the vegdist function,
# on the first 10 samples
fieldids <- dataset$fieldid
distVegan <- vegdist(abundance_matrix[fieldids, ][1:10,],
                    method = "jaccard",
                    binary = TRUE)
simOwn <- get_shared_species(abundance_matrix[fieldids, ][1:10,])

plot(simOwn, 1 - as.matrix(distVegan))

# In large datasets, out implementation is much faster as it uses matrix multiplication tricks
JaccardSim <- get_shared_species(abundance_matrix[fieldids, ])

# Save the final dataset and the Jaccard similarity matrix
save(dataset, JaccardSim, file = "~/HubbellGLM-paper/data/data_GMTP_clean.rdata")


