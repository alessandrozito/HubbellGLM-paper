# This file merges all data.
library(tidyverse)
library(lubridate)
library(HubbellGLM)
library(BNPvegan)
library(splines)
library(MASS)
library(sf)
library(pheatmap)
library("rnaturalearth")
library("rnaturalearthdata")
library(geosphere)
library(patchwork)
library(conleyreg)
library(slider)

#----------!!!!!!
# NOTE: the script is illustrative. Files are available upon request.
# To jump to the second preprocessing, run 1_GMTP_create_dataset.R
#----------!!!!!!

#------------------------------------------------------------ Part 1 - load data
# Load the GMPT data
dfGMPT <- readRDS("~/HubbellGLM/tutorials/data/DS-20GMP_01_37_merge_with_Lots2.rds.gzip")
head(dfGMPT)

# Load the environmental data.
dfERA <- read_csv("~/HubbellGLM/tutorials/data/GMTP_climate_data/GMTP_ERA5_climate.csv")
dfTerra <- read_csv("~/HubbellGLM/tutorials/data/GMTP_climate_data/GMTP_Terra_climate.csv")
df_HumanInfluence <- read_csv("~/HubbellGLM/tutorials/data/GMTP_climate_data/GMTP_human_influence.csv")
df_NDVI <- read_csv("~/HubbellGLM/tutorials/data/GMTP_climate_data/GMTP_GIMMS_NDVI.csv")
df_wetlands <- read_csv("~/HubbellGLM/tutorials/data/GMTP_climate_data/GMTP_wetland_cover.csv")

# For each location, we have to calculate
# - MAX NDVI at the annual level for that location
# - Human footprint data (already at the annual level)
# - Annual average temperature/wind speed/humidity/precipitation at each location.
# - Weekly average temperature/wind speed/humidity/precipitation at each that location.
# - Weekly average anomalies temperature/wind speed/humidity/precipitation at each that location.
# - Week of collection and day of the year of collection, and flip it 6 months if the
#   below equator

############################################################
# dfGMPT preprocess to get all locations and collection dates
############################################################
df_collectionDates <- dfGMPT %>%
  dplyr::select(fieldid, Latitude_clean, Longitude_clean,
                collection_start_date_clean, collection_end_date_clean) %>%
  distinct() %>%
  rowwise() %>%
  mutate(date = list(seq(collection_start_date_clean, collection_end_date_clean, by = "day"))) %>%
  ungroup() %>%
  unnest(date)

############################################################
# NDVI data
############################################################
# Preliminary useful dataset containing all dates for each collection period,
# including halves of months (useful to merge later)
dates <- seq(as.Date("2010-01-01"), as.Date("2016-12-31"), by = "day")
df_dates <- data.frame(date = dates,
                       day_of_month = day(dates),
                       year = year(dates),
                       month = month(dates)) %>%
  mutate(half = ifelse(day_of_month <= 15, 1, 2))

# Transform the NDVI data into long format containing all days of the week,
# and calculate annual max, annual min, and annual avg NDVI at each location
df_NDVI_alldates <- df_NDVI %>%
  full_join(df_dates, by = c("year", "month", "half"), relationship = "many-to-many") %>%
  group_by(Latitude_clean, Longitude_clean) %>%
  mutate(site = cur_group_id()) %>%
  arrange(site, date) %>%
  dplyr::select(Latitude_clean, Longitude_clean, year, date, NDVI) %>%
  group_by(Latitude_clean, Longitude_clean, year) %>%
  mutate(NDVI_year_max = max(NDVI),
         NDVI_year_min = min(NDVI),
         NDVI_year_avg = mean(NDVI)) %>%
  ungroup()

############################################################
# Human Influence data
############################################################
# Reshape Human influence to long format
df_HumanInfluence_long <- df_HumanInfluence %>%
  pivot_longer(
    cols = matches("^(hfp\\d{4}|population_density_\\d{4}|night_light_\\d{4})"),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    year = str_extract(variable, "\\d{4}"),
    variable2 = str_replace(variable, "_?\\d{4}$", ""), # remove year from end (and _ if present)
    variable2 = str_replace(variable2, "_$", "")    # clean trailing underscore
  ) %>%
  pivot_wider(
    id_cols = c(Longitude_clean, Latitude_clean, wildareas_human_footprint_index, year),
    names_from = variable2, values_from = value
  ) %>%
  rename(
    hfp = hfp,
    population_density = population_density,
    night_light = night_light,
    wildareas_human_footpr = wildareas_human_footprint_index
  ) %>%
  mutate(
    year = as.integer(year),
    hfp = as.numeric(hfp),
    population_density = as.numeric(population_density),
    night_light = as.numeric(night_light)
  ) %>%
  arrange(year) %>%
  distinct() %>%
  ungroup() %>%
  # Fill the NA values fo hfp. They correspond to a case when night_light is zero.
  mutate(hfp = case_when(is.na(hfp) ~ 0,
                         TRUE ~ hfp))

############################################################
# Climate data (ERA5)
############################################################
# Calculate annual average for each covariate at each location.
# To calculate the average for the wind speed and wind direction, need to
# shift to polar coordinates. The surface_solar_radiation_downwards_sum
# is not available for each location.
dfERA_avg <- dfERA %>%
  dplyr::mutate(year = year(date),
         wind_direction_rad = wind_direction * pi / 180,
         u_wind = wind_speed * cos(wind_direction_rad),
         v_wind = wind_speed * sin(wind_direction_rad),
         u_wind_anml = wind_speed_anml * cos(wind_direction_rad),
         v_wind_anml = wind_speed_anml * sin(wind_direction_rad)) %>%
  dplyr::group_by(Latitude_clean, Longitude_clean, year) %>%
  dplyr::mutate(across(.cols = where(is.numeric),
                ~ mean(.x, na.rm = TRUE),
                .names = "{.col}_year_avg")) %>%
  mutate(wind_speed_year_avg_tr = sqrt(u_wind_year_avg^2 + v_wind_year_avg^2),
         wind_speed_anml_year_avg_tr = sqrt(u_wind_anml_year_avg^2 + v_wind_anml_year_avg^2)) %>%
  dplyr::select(-u_wind_year_avg, -v_wind_year_avg,
                -u_wind_anml_year_avg, -v_wind_anml_year_avg,
                -wind_direction_rad_year_avg, -year) %>%
  ungroup()

############################################################
# Climate data (TERRA)
############################################################
# These are monthly data at 5km resolution. We will probably not use them, but still
# worth including them
dfTerra_avg <- dfTerra %>%
  group_by(Latitude_clean, Longitude_clean, year) %>%
  mutate(across(.cols = where(is.numeric),
                ~ mean(.x, na.rm = TRUE),
                .names = "{.col}_year_avg")) %>%
  full_join(df_dates %>% dplyr::select(-half, -day_of_month),
            by = c("year", "month"), relationship = "many-to-many") %>%
  filter(!is.na(date)) %>%
  ungroup() %>%
  distinct() %>%
  rename_with(~ paste0("Terra_", .),
              .cols = -c(Latitude_clean, Longitude_clean, year, month, date))

############################################################
# Merge everything and calculate average values in the week of collection
############################################################
df_climate <- df_collectionDates %>%
  mutate(month = month(collection_start_date_clean),
         year = year(collection_start_date_clean)) %>%
  # Merge with Wetlands
  left_join(df_wetlands,  by = c("Latitude_clean", "Longitude_clean", "year")) %>%
  # Merge with NDVI
  left_join(df_NDVI_alldates %>%
              dplyr::select(-year), by = c("Latitude_clean", "Longitude_clean", "date")) %>%
  # Merge with Human Influence data
  left_join(df_HumanInfluence_long, by = c("Latitude_clean", "Longitude_clean", "year")) %>%
  # Merge with ERA5
  left_join(dfERA_avg %>%
              dplyr::select(-year), by = c("Latitude_clean", "Longitude_clean", "date")) %>%
  # Merge with Terra
  left_join(dfTerra_avg %>%
              dplyr::select(-year, -month), by = c("Latitude_clean", "Longitude_clean", "date")) %>%
  # Calculate all averages
  group_by(fieldid, Latitude_clean, Longitude_clean, collection_start_date_clean,
           collection_end_date_clean) %>%
  dplyr::select(-year, -month) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(wind_speed = sqrt(u_wind^2 + v_wind^2),
         wind_speed_anml = sqrt(u_wind_anml^2 + v_wind_anml^2))


# Merge with site covariates from ever
df_covariates <- dfGMPT %>%
  dplyr::select(fieldid, site_clean, country_clean, most_frequent_note,
                realm, lot_habitat, elevation) %>%
  distinct() %>%
  left_join(df_climate, by = "fieldid")

# BIN matrix
Xbins <- dfGMPT %>%
  filter(!is.na(bin_uri)) %>%
  group_by(fieldid, bin_uri) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = bin_uri, values_from = n, values_fill = 0) %>%
  ungroup() %>%
  column_to_rownames("fieldid") %>%
  as.matrix()
saveRDS(Xbins, file = "~/HubbellGLM/tutorials/data/BIN_matrix.rds.gzip", compress = "gzip")

df_occurrences <- data.frame(fieldid = rownames(Xbins),
                             n = rowSums(Xbins),
                             y = rowSums(Xbins > 0))

# Final data frame
data <- df_covariates %>%
  left_join(df_occurrences, by = "fieldid") %>%
  mutate(
    realm = case_when(is.na(realm) ~ "Nearctic",
                      TRUE ~ realm),
    combined_traps = case_when(
      most_frequent_note == "2 Malaise samples combined (+GMP#02510)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#02516)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#02520)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#02525)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#02528)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#02532)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#02536)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#02548)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#02552)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#02556)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#02558)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#05618)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#05619)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#05624)" ~ TRUE,
      most_frequent_note == "2 Malaise samples combined (+GMP#05627)" ~ TRUE,
      TRUE ~ FALSE
    ),
    trap_damaged = case_when(
      most_frequent_note == "trap fell, bottle empty" ~ TRUE,
      most_frequent_note == "Grizzly punctured & drained bottle" ~ TRUE,
      most_frequent_note == "Jar dry when collected" ~ TRUE,
      most_frequent_note == "Bear damage-bugs and liquid gone" ~ TRUE,
      most_frequent_note == "Gap in dates due to bear issues" ~ TRUE,
      most_frequent_note == "Gap in dates due to bear issues: trap moved to new location" ~ TRUE,
      most_frequent_note == "Entire month of July missing because of bear" ~ TRUE,
      most_frequent_note == "trap damage from wildlife" ~ TRUE,
      most_frequent_note == "incomplete sample due to animal" ~ TRUE,
      most_frequent_note == "Malaise trap collapsed upon collection" ~ TRUE,
      most_frequent_note == "tent collapsed" ~ TRUE,
      most_frequent_note == "tent collapsed, now pegged in by star pickets" ~ TRUE,
      most_frequent_note == "trap partially fallen over" ~ TRUE,
      most_frequent_note == "windstorm Jul 23-24, trap shook a lot and splashed liquid into upper jar, caused some kind of white crystal precipitate to form" ~ TRUE,
      most_frequent_note == "After trap repair" ~ TRUE,
      TRUE ~ FALSE
    )
  )

data_fields <- data
abundance_matrix <- Xbins
# Save output
saveRDS(abundance_matrix, file = "~/HubbellGLM-paper/data/abundance_matrix.rds.gzip")
saveRDS(data_fields, file = "~/HubbellGLM-paper/data/data_GMTP_merged_covariates.rds.gzip")

