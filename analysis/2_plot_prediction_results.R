# This file plots the prediction results for our model
# The figures produced in this scripts are the followings
#   - Exploratory data analysis
#     -- Figure 1 b, c - location of the points and covariates for Amazon
#     -- Figure S1 - Examples of accumulation curves
#     -- Figure S2 - relationship between the geographic distance and Jaccard similarity
#     -- Figure S3 - covariates on the global scale
#     -- Figure S4 - study design and relationship between covariates and diversity
#   - Model results
#     -- Figure 1 e, f - prediction of accumulation curves and global diversity indices
#     -- Figure 4 - scenario evaluation for hfp and aet

#--- Load the R packages
library(tidyverse)
library(lubridate)
library(ggpubr)
library(HubbellGLM)
library(splines)
library("rnaturalearth")
library("rnaturalearthdata")
library(patchwork)
library(lmtest)
library(multcomp)
library(scales)
library(geosphere)
library(stargazer)

#-------- Load the data
load("~/HubbellGLM-paper/data/data_GMTP_clean.rdata")

#-------- Load the global raster for year 2014
df_Raster <- readRDS("~/HubbellGLM-paper/data/climate/Raster_2014_merged.rds")

#===============================================================================
# Figure 1 - Panel b: location of sites in the world
#===============================================================================

# Load the world, excluding Antarctica
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- world[world$name !=  "Antarctica", ]

# Plot the zones
zone_colors <- c("Tropical" = "red3", "Dry" = "#D6BD8DFF", "Temperate" = "blue3",
                 "Continental" = "lightblue", "Polar" = "gray50")

p_locs <- ggplot() +
  theme_void() +
  #geom_point(data = dataset, aes(x = Longitude, y= Latitude)) +
  geom_sf(data = world, fill = "gray95", color = "gray30", size = 0.15) +
  theme(panel.grid = element_blank(), aspect.ratio = 0.5) +
  coord_sf(expand = FALSE) +
  ylab("Latitude") +
  xlab("Longitude") +
  geom_point(
    data = dataset,
    aes(x = Longitude, y = Latitude, fill = zone), shape = 21, size = 2.5) +
  scale_fill_manual(name = "Climatic zone", values = zone_colors,
                    guide = guide_legend(override.aes = list(size = 4)))
ggsave(p_locs, filename = "~/HubbellGLM-paper/figures/Figure1b_locations.png",
       height = 4.04, width = 8.96)

#===============================================================================
# Figure 1 - Panel c: HFP and AET in the Amazon rainforest
#===============================================================================

#--------------- Palettes
aet_palette <- c("#290AD8FF", "#264DFFFF", "#3FA0FFFF", "#72D9FFFF", "#AAF7FFFF", "#E0FFFFFF",
                 "#FFFFBFFF", "#FFE099FF", "#FFAD72FF", "#F76D5EFF", "#D82632FF", "#A50021FF")
hfp_palette <- c("white", "antiquewhite1", "#F8B58BFF", "#F59E72FF", "#F2855DFF",
                 "#EF6A4CFF", "red2", "darkred", "black")

#--------------- AET in the Amazon
ylims <- c(-35, 20)
xlims <- c(-90, -31)
df_Raster_Amazon <- df_Raster %>%
  dplyr::filter(between(x, xlims[1], xlims[2]),
                between(y, ylims[1], ylims[2]))

pAmazonAET <- ggplot() +
  theme_classic() +
  # Map layer: Use the filtered dataframe
  geom_tile(data = df_Raster_Amazon,
            aes(x = x, y = y, fill = aet)) + # applied scale factor
  # Border layer: Canada outline
  geom_sf(data = world,
          fill = NA, color = "gray30", size = 0.2) +
  scale_fill_gradientn(name = "AET\n(mm/year)", colours = rev(aet_palette),
                       na.value = NA, limits = c(0, 140)) +
  coord_sf(expand = FALSE)+
  ylim(ylims) +
  xlim(xlims) +
  theme(axis.title = element_blank())

ggsave(plot = pAmazonAET, filename = "~/HubbellGLM-paper/figures/Figure1c_Amazon_AET.pdf",
       height = 4.40, width = 7.72)

#--------------- HFP in the Amazon
pAmazonHFP <- ggplot() +
  theme_classic() +
  # Map layer: Use the filtered dataframe
  geom_tile(data = df_Raster_Amazon,
            aes(x = x, y = y, fill = hfp)) + # applied scale factor
  # Border layer: Canada outline
  geom_sf(data = world,
          fill = NA, color = "gray30", size = 0.2) +
  scale_fill_gradientn(name = "HFP", colours = hfp_palette,
                       na.value = NA, limits = c(0, 50)) +
  coord_sf(expand = FALSE)+
  ylim(ylims) +
  xlim(xlims)+
  theme(axis.title = element_blank())

ggsave(plot = pAmazonHFP, filename = "~/HubbellGLM-paper/figures/Figure1c_Amazon_HFP.pdf",
       height = 4.40, width = 7.72)


#===============================================================================
# Figure S1 - Examples of accumulation curves
#===============================================================================
n <- 1:1000
eta <- c(- 1, 0, 1)
alpha <- exp(eta)
sigma <- c(-0.5, -0.25, 0, 0.1, 0.25)
df <- data.frame()

# Build the data frame
for (a in alpha) {
  for (s in sigma) {
    mu_n <- HubbellGLM:::polyseries_mean(size = n,
                                         alpha = rep(a, length(n)),
                                         sigma = s)
    df <- rbind(df, data.frame(n, alpha = a, sigma = factor(s), mu_n))
  }
}

df$eta <- log(df$alpha)
df$eta <- factor(paste("eta ==", df$eta), levels = paste("eta ==", eta))

# Plot
pAccum <- ggplot(df, aes(x = n, y = mu_n, color = sigma)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ eta, labeller = label_parsed) + # Parses the string into Greek letters
  scale_color_manual(values = c("#C5D358", "#318B27", "#FF0000", "#4B71D0", "#000080")) +
  labs(y = "Species richness (mean)",
       color = expression(sigma)) +            # Uses expression() for the legend title
  theme_bw(base_size = 14)

ggsave(pAccum, width = 10.35, height = 3.74,
       filename = "~/HubbellGLM-paper/figures/FigureS1_AccumulationCurves.pdf")


#===============================================================================
# Figure S2 - Relationship between geographic distance and points
#===============================================================================
fieldids <- dataset$fieldid
JaccardSim_upper <- JaccardSim[upper.tri(JaccardSim)]
# Calculate distance between coordinates of each two points
distances_sf <- distm(dataset %>%
                        dplyr::select(Longitude, Latitude) %>%
                        as.matrix())
distances_sf_upper <- distances_sf[upper.tri(distances_sf)]

subset_points <- 1:length(JaccardSim_upper)
df_plot_dist <- data.frame(
  geo_dist_km = distances_sf_upper[subset_points] / 1000,
  similarity = JaccardSim_upper[subset_points])

p_dist <- ggplot(df_plot_dist, aes(x = geo_dist_km, y = similarity)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "Density\n(n. of points)", labels = comma)+
  #geom_smooth(method = "gam", color = "darkred", size = 1) +
  scale_x_continuous(labels = comma) +
  labs(
    x = "Geographic distance (km)",
    y = "Fraction of shared species"
  ) +
  theme_bw(base_size = 14) +
  facet_wrap(~"Relationship between similarity and distance")

ggsave(p_dist, width = 9.34, height = 4.55,
       filename = "~/HubbellGLM-paper/figures/FigureS2_Geographic_Similarity.pdf")


#===============================================================================
# Figure S3 - Covariates on the global scale
#===============================================================================

#--------------- Global AET
p_aet <- ggplot() +
  theme_bw() +
  geom_tile(data = df_Raster,
            aes(x = x, y = y, fill = aet), alpha = 0.65) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.15) +
  theme(axis.title = element_blank()) +
  coord_sf(expand = FALSE) +
  scale_fill_gradientn(name = "AET\n(mm/year)", colours = rev(aet_palette))

ggsave(plot = p_aet,
       filename = "~/HubbellGLM-paper/figures/FigureS3a_AET_raster.png",
       height = 2.98, width = 7.72, dpi = 150)

#--------------- Global HFP
p_hfp <- ggplot() +
  theme_bw() +
  geom_tile(data = df_Raster %>% filter(!is.na(hfp)), na.rm = TRUE,
            width = 0.8, height = 0.8,
            aes(x = x, y = y, fill = hfp), alpha = 0.65) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.15) +
  theme(axis.title = element_blank()) +
  coord_sf(expand = FALSE) +
  scale_fill_gradientn(name = "HFP", colours = hfp_palette)

ggsave(plot = p_hfp,
       filename = "~/HubbellGLM-paper/figures/FigureS3b_HFP_raster.png",
       height = 2.98, width = 7.72, dpi = 150)

#--------------- Global climatic zone
zone_colors <- c("Tropical" = "red3", "Dry" = "#D6BD8DFF", "Temperate" = "blue3",
                 "Continental" = "lightblue", "Polar" = "gray50")

p_zone <- ggplot() +
  theme_bw() +
  geom_tile(data = df_Raster %>% dplyr::filter(!is.na(zone)),
            aes(x = x, y = y, fill = zone), alpha = 0.65) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.15) +
  theme(axis.title = element_blank()) +
  coord_sf(expand = FALSE) +
  scale_fill_manual(name = "Climatic zone", values = zone_colors)

ggsave(plot = p_zone,
       filename = "~/HubbellGLM-paper/figures/FigureS3c_Zone_raster.png",
       height = 2.98, width = 7.72, dpi = 150)

#===============================================================================
# Figure S4 - study design and relationship between covariates and diversity
#===============================================================================

#---- Figure S4a a: Distribution of the sampling sites locations
df <- dataset %>%
  dplyr::select(fieldid, site, zone, collection_start_date,
                collection_end_date, realm,
                NDVI_year_max) %>%
  distinct() %>%
  rowwise() %>%
  dplyr::mutate(date = list(seq(collection_start_date, collection_end_date, by = "day"))) %>%
  ungroup() %>%
  unnest(date) %>%
  dplyr::mutate(min_start = min(collection_start_date, na.rm=TRUE),
                date_diff = as.integer(date - min_start) + 1)

site_levels <- dataset %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(Latitude = mean(Latitude),
                   Longitude = mean(Longitude)) %>%
  arrange(desc(Latitude)) %>%
  pull(site)

year_ticks <- df %>%
  dplyr::group_by(year = year(date)) %>%
  dplyr::summarise(date_diff = min(date_diff)) %>%
  as.data.frame() %>%
  rbind(data.frame(year = 2011,
                   date_diff = 339)) %>%
  arrange(year)

df$site <- factor(df$site, levels = site_levels)

p_sites <- ggplot(df, aes(x = date_diff, y = site)) +
  geom_point(aes(color = zone), size = 1) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_y_discrete(limits = rev(levels(df$site))) +
  scale_x_continuous(
    breaks = year_ticks$date_diff,         # places ticks at start of each year
    labels = year_ticks$year,              # label them with the year
    name = 'Collection year'
  ) +
  ylab("Sampling site (sorted by latitude)") +
  scale_color_manual(name = "Climatic zone", values = zone_colors,
                     guide = guide_legend(override.aes = list(size = 4)))
ggsave(p_sites, filename = "~/HubbellGLM-paper/figures/FigureS4a_sites.pdf",
       height = 10.03, width = 4.83)

#---- Panel b and c
dataset$Ssigma <- HubbellGLM:::inv_polyseries(dataset$y, dataset$n, 0.5692932)/0.5692932

pAlpha_hfp <- ggplot(dataset, aes(x = hfp, y = log(Ssigma), color = zone)) +
  geom_point(alpha = 0.7, size = 2)+
  geom_smooth(method = "lm", colour = "black")+
  facet_grid(~zone) +
  scale_color_manual(values = zone_colors) +
  theme_bw()+
  xlab("HFP")+
  scale_color_manual(name = "Climatic zone", values = zone_colors,
                     guide = guide_legend(override.aes = list(size = 3)))

# Log of S_sigma AET
pAlpha_aet <- ggplot(dataset, aes(x = aet, y = log(Ssigma))) +
  geom_point(alpha = 0.7, aes(color = zone), size = 2)+
  geom_smooth(method = "lm", colour = "black") +
  theme_bw()+
  xlab("AET")+
  facet_grid(~"Annual actual evapotranspiration")+
  scale_color_manual(name = "Climatic zone", values = zone_colors,
                     guide = guide_legend(override.aes = list(size = 4)))+
  stat_cor(
    method = "pearson",
    label.x.npc = "center",
    label.y.npc = "bottom",
    size = 4)

p_joind <- pAlpha_aet/pAlpha_hfp + plot_layout(guides = "collect")
ggsave(p_joind, filename = "~/HubbellGLM-paper/figures/FigureS4bc_Ssigma_hfp_aet.pdf",
       width = 7.23, height = 5.59)

#---- Panel d - relationship between alpha diversity and evapotranspiration
dataset_long <- dataset %>%
  dplyr::select(
    Ssigma,
    zone,
    `Latitude` = Latitude,
    `PET` = pet,
    `Vapor pressure deficit` = vpd_year_avg,
    `Wetlands` = wetlands,
    `Temperature` = temperature_2m_year_avg,
    `NDVI` = NDVI_year_max,
    `Precipitation` = total_precipitation_sum_year_avg,
    `Wind speed` = wind_speed_year_avg,
    `Rel. Humidity` = relative_humidity_year_avg
  ) %>%
  pivot_longer(
    cols = c(Latitude, PET, `Vapor pressure deficit`, Wetlands, Temperature, NDVI, Precipitation, `Wind speed`, `Rel. Humidity`),
    names_to = "Variable",
    values_to = "Value"
  )

p_environ_covariates <- ggplot(dataset_long, aes(x = Value, y = log(Ssigma))) +
  geom_point(alpha = 0.7, aes(color = zone), size = 0.7) +
  #geom_smooth(colour = "black") +
  stat_cor(
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = "bottom",
    size = 3
  ) +
  facet_wrap(~ Variable, scales = "free_x", ncol = 3) +
  scale_color_manual(
    name = "Climatic zone",
    values = zone_colors,
    guide = guide_legend(override.aes = list(size = 4))
  ) +
  theme_bw() +
  theme(legend.position = "none")+
  labs(x = NULL, y = "log(Ssigma)")

ggsave(p_environ_covariates, filename = "~/HubbellGLM-paper/figures/FigureS4d_Ssigma_otherst.pdf",
       width = 7.77, height = 5.92)


#===============================================================================
# Figure 1 e, f - prediction of accumulation curves and global diversity indices
#===============================================================================

# Function to extract the fourier series given a certain week of the year
fourier_week <- function(week, k = 1) {
  X <- NULL
  for(j in 1:k) {
    X <- cbind(X, sin(2 * j *pi * week / 52))
    X <- cbind(X, cos(2 * j * pi * week / 52))
    colnames(X)[(ncol(X) - 1):ncol(X)] <- paste0(c("week_sin", "week_cos"), j)
  }
  return(X)
}

data_for_pred <- dataset %>%
  dplyr::mutate(wind_speed = wind_speed_year_avg) %>%
  dplyr::select(n, y, aet, hfp, zone, realm, week_adj, Latitude, collection_days,
                wind_speed, wetlands, resid_temperature_lat)

# We use the model with no weather anomalies
fitHubbell <- HubbellGLM(cbind(n, y) ~ 1 + aet + realm +
                           fourier_week(week_adj, k = 1) +
                           ns(Latitude, 6) + zone * hfp +
                           wind_speed + wetlands + collection_days +
                           resid_temperature_lat,
                         data = data_for_pred,
                         family = hubbell(sigma = 0.5617625))

#------ Make all the predictions at the global scale
n <- 2000
df_Raster$n <- n
df_Raster$week_adj <- 26
df_Raster$collection_days <- 7
df_Raster$realm[df_Raster$realm %in% c("Antarctic", "Oceania")] <- NA
df_Raster$realm <- factor(df_Raster$realm, levels = unique(df_Raster$realm))

# Current level of diversity
S_sigma_curr = exp(predict(fitHubbell,
                              newdata = df_Raster %>%
                                        dplyr::mutate(hfp = pmin((hfp + 0.01), 50))))/fitHubbell$sigma

#--- Level of diversity given by a 10% increase in HFP
S_sigma_hfp10 = exp(predict(fitHubbell, newdata = df_Raster %>%
                                         dplyr::mutate(hfp = pmin((hfp + 0.01) * 1.1, 50))))/fitHubbell$sigma

#--- Predicted richness at 2000 samples
richness_curr <- predict(fitHubbell, newdata = df_Raster %>%
                                      dplyr::mutate(hfp = pmin((hfp + 0.01), 50)), type = "response")

#--- Predicted alpha
alpha_curr <- rep(NA, length(richness_curr))
alpha_curr[!is.na(richness_curr)] <- HubbellGLM:::inv_mean_dirichlet_process(mu_target = richness_curr[!is.na(richness_curr)],
                                                                size = n)
#--- Predicted Shannon, Simpson and Hill, based on the inferred alpha
Shannon_curr <- digamma(alpha_curr + 1) - digamma(1) # Model-based shannon diversity
Simpson_curr <- 1/(alpha_curr + 1) # Model-based simpson index
Tsallis_curr <- alpha_curr * beta(alpha_curr, 0.5)

# Combine everything in a single raster and plot all indices
df_vars <- data.frame(Latitude = df_Raster$y,
                      Longitude = df_Raster$x,
                      S_sigma_curr = S_sigma_curr,
                      alpha_curr = alpha_curr,
                      Shannon_curr = Shannon_curr,
                      Simpson_curr = Simpson_curr,
                      Tsallis_curr = Tsallis_curr,
                      richness_curr = richness_curr)


#----- Plot of the regression based diversity

pal_emrl <- c("antiquewhite2","#D3F2A3FF", "#97E196FF", "#6CC08BFF",
              "#02734AFF", "#217A79FF",  "#105965FF", "#03334AFF", "black")

p_Sigma <- ggplot() +
  theme_void() +
  geom_tile(data = df_vars %>% drop_na(),
            aes(x = Longitude, y = Latitude, fill = S_sigma_curr),
            width = 1.01, height = 1.01) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 0.5)+
  coord_sf(expand = FALSE) +
  guides(
    color = "none",
    fill = guide_colourbar(
      barheight = unit(6.5, "cm"),
      barwidth  = unit(0.5, "cm"))) +
  scale_fill_gradientn(name = "", colours = c("gray60", pal_emrl))
p_Sigma
ggsave(plot = p_Sigma, filename = "~/HubbellGLM-paper/figures/Figure1f_Sigma.png",
       width = 5.81, height = 2.89, dpi = 200)

#-------- alpha diversity
p_alpha <- ggplot() +
  theme_void() +
  geom_tile(data = df_vars %>% drop_na(),
            aes(x = Longitude, y = Latitude, fill = alpha_curr),
            width = 0.58, height = 0.58) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 0.5)+
  coord_sf(expand = FALSE) +
  guides(
    color = "none",
    fill = guide_colourbar(
      barheight = unit(6.5, "cm"),
      barwidth  = unit(0.5, "cm"))) +
  scale_fill_gradientn(name = "", colours = c("gray60", pal_emrl))
p_alpha
ggsave(plot = p_alpha, filename = "~/HubbellGLM-paper/figures/Figure1f_alpha.png",
       width = 5.81, height = 2.89, dpi = 200)


#----- Shannon
pal_shannon <-c("#1E8E99FF", "#51C3CCFF", "#B2FCFFFF", "#CCFEFFFF", "#E5FFFFFF",
                "#FFE5CCFF", "#FFAD65FF", "#FF8E32FF", "#CC5800FF", "#993F00FF")
p_shannon <- ggplot() +
  theme_void() +
  geom_tile(data = df_vars %>% drop_na(),
            aes(x = Longitude, y = Latitude, fill = Shannon_curr),
            width = 0.58, height = 0.58) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 0.5)+
  coord_sf(expand = FALSE) +
  guides(
    color = "none",
    fill = guide_colourbar(
      barheight = unit(7, "cm"),
      barwidth  = unit(0.5, "cm"))) +
  scale_fill_gradientn(name = "", colours = pal_shannon)
p_shannon
ggsave(plot = p_shannon, filename = "~/HubbellGLM-paper/figures/Figure1f_shannon.png",
       width = 5.81, height = 2.89, dpi = 200)

#----- Simpson
pal_simpson <- c("white", "#E5FFFFFF", "#B2F2FFFF", "#7FD4FFFF", "#65BFFFFF", "#4CA5FFFF", "#003FFFFF")
p_simpson <- ggplot() +
  theme_void() +
  geom_tile(data = df_vars %>% drop_na(),
            aes(x = Longitude, y = Latitude, fill = Simpson_curr),
            width = 0.58, height =  0.58) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 0.5)+
  coord_sf(expand = FALSE) +
  guides(
    color = "none",
    fill = guide_colourbar(
      barheight = unit(6.5, "cm"),
      barwidth  = unit(0.5, "cm"))) +
  scale_fill_gradientn(name = "", colours = pal_simpson)
p_simpson
ggsave(plot = p_simpson, filename = "~/HubbellGLM-paper/figures/Figure1f_simpson.png",
       width = 5.81, height = 2.89, dpi = 200)

#----- Hill
palette_tsallis <- c("#F6C866", "#F2AB67", "#EF8F6B", "#ED7470", "#BF6E97",
                     "#926AC2", "#6667EE", "#4959C7", "#2D4A9F", "#173C78")
p_tsallis <- ggplot() +
  theme_void() +
  geom_tile(data = df_vars %>% drop_na(),
            aes(x = Longitude, y = Latitude, fill = Tsallis_curr),
            width = 0.58, height =  0.58) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 0.5)+
  coord_sf(expand = FALSE) +
  guides(
    color = "none",
    fill = guide_colourbar(
      barheight = unit(6.5, "cm"),
      barwidth  = unit(0.5, "cm"))) +
  scale_fill_gradientn(name = "", colours = palette_tsallis)
p_tsallis
ggsave(plot = p_tsallis, filename = "~/HubbellGLM-paper/figures/Figure1f_hill.png",
       width = 5.81, height = 2.89, dpi = 200)

#===============================================================================
# Figure 4 - simulation scenarios
#===============================================================================

#----------- REPRODUCE FIGURE 4b

#--- Predicted richness at 2000 samples, with a 10% increase in hfp
richness_hfp10 <- predict(fitHubbell,
                          newdata = df_Raster %>%
                            dplyr::mutate(hfp = pmin((hfp + 0.01) * 1.1, 50)),
                          type = "response")
#--- Predicted richness at 2000 samples, with a 5 unit increase in hfp
richness_hfp5 <- predict(fitHubbell,
                         newdata = df_Raster %>%
                           dplyr::mutate(hfp = pmin((hfp + 0.01) + 5, 50)),
                         type = "response")

#--- Calculate variations in % scale
df_vars$variation_richness10 <- 100 * (richness_hfp10 - df_vars$richness_curr)/df_vars$richness_curr
df_vars$variation_richness5 <- 100 * (richness_hfp5 - df_vars$richness_curr)/df_vars$richness_curr

p_variations_hfp <- ggplot() +
  theme_void() +
  geom_tile(data = df_vars %>% drop_na(variation_richness10),
            aes(x = Longitude, y = Latitude, fill = variation_richness10),
            width = 0.58, height = 0.58) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 0.5)+
  coord_sf(expand = FALSE) +
  guides(
    color = "none",
    fill = guide_colourbar(
      barheight = unit(7, "cm"),
      barwidth  = unit(0.5, "cm"))) +
  scale_fill_gradientn(
    name   = "",
    colours = c("darkred", "red2", "#F2855DFF", "antiquewhite1", "lightblue1", "#290AD8FF"),
    limits = c(-5.3, 2.8),
    values = scales::rescale(c(seq(-5.3, 0, length.out = 4)[1:3], 0, 0.2, 2.8), from = c(-5.3, 2.8)),
    oob = scales::squish,
    breaks = seq(-5, 2.5, by = 2.5),
    labels = function(x) paste0(x, " %"))
ggsave(p_variations_hfp, filename = "~/HubbellGLM-paper/figures/Figure4b_change_richness_v1.png",
       height = 4.04, width = 8.96)

p_variations_hfp_2 <- ggplot() +
  theme_void() +
  geom_tile(data = df_vars %>% drop_na(variation_richness5),
            aes(x = Longitude, y = Latitude, fill = variation_richness5),
            width = 0.58, height = 0.58) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 0.5)+
  coord_sf(expand = FALSE) +
  guides(
    color = "none",
    fill = guide_colourbar(
      barheight = unit(7, "cm"),
      barwidth  = unit(0.5, "cm"))) +
  scale_fill_gradientn(
    name   = "",
    colours = c("darkred", "red2", "antiquewhite1", "#3FA0FFFF", "#264DFFFF"),
    limits = c(-8, 8),
    oob = scales::squish,
    breaks = round(seq(-8, 8, by = 2)),
    labels = function(x) paste0(x, " %"))

ggsave(p_variations_hfp_2,
       filename = "~/HubbellGLM-paper/figures/Figure4b_change_richness_v2.png",
       height = 4.04, width = 8.96)


#----------- REPRODUCE FIGURE 4a
richness_curr0 <- predict(fitHubbell,
                          newdata = df_Raster, type = "response")
richness_hfp0 <- predict(fitHubbell,
                         newdata = df_Raster %>%
                                      dplyr::mutate(hfp = 0),
                         type = "response")
df_vars$variation_richness0 <- 100 * (richness_hfp0 - richness_curr0)/richness_curr0
maxabs0 <- max(abs(df_vars$variation_richness0), na.rm = TRUE)

p_hfp_zero <- ggplot() +
  theme_void()+
  geom_tile(data = df_vars %>% drop_na(variation_richness0),
            aes(x = Longitude, y = Latitude, fill = variation_richness0),
            width = 0.58, height = 0.58) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  coord_sf(expand = FALSE) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    aspect.ratio = 0.5
  ) +
  guides(
    color = "none",
    fill = guide_colourbar(
      barheight = unit(7, "cm"),
      barwidth  = unit(0.5, "cm"))) +
  scale_fill_gradientn(
    name   = "",
    colours = c("#290AD8FF", "#264DFFFF", "#3FA0FFFF", "antiquewhite1",
                "#6CC08BFF", "#02734AFF", "#03334AFF"),
    limits = c(-maxabs0, maxabs0),
    oob = scales::squish,
    breaks = round(seq(-60, 60, by = 20)),
    labels = function(x) paste0(x, " %"))
p_hfp_zero

ggsave(plot = p_hfp_zero,
       filename = "~/HubbellGLM-paper/figures/Figure4a_raster2014_hfp0.png",
       height = 4.04, width = 8.96)

#----------- REPRODUCE FIGURE 4c
richness_currAET <- predict(fitHubbell,
                            newdata = df_Raster %>%
                              dplyr::mutate(aet = aet + 0.01),
                            type = "response")
richness_AET10 <- predict(fitHubbell,
                          newdata = df_Raster %>%
                            dplyr::mutate(aet = (aet + 0.01) * 1.1),
                          type = "response")

df_vars$variation_richnessAET <- 100 * (richness_AET10 - richness_currAET)/richness_currAET

pal_emrl <- c("antiquewhite1", "#E5FFB6", "#D3F2A3FF", "#97E196FF", "#6CC08BFF",
              "#02734AFF", "#217A79FF",  "#105965FF", "#03334AFF")

p_variations_aet <- ggplot() +
  theme_void() +
  geom_tile(data = df_vars %>% drop_na(variation_richnessAET),
            aes(x = Longitude, y = Latitude, fill = variation_richnessAET),
            width = 0.58, height = 0.58) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        aspect.ratio = 0.5)+
  coord_sf(expand = FALSE) +
  guides(
    color = "none",
    fill = guide_colourbar(
      barheight = unit(7, "cm"),
      barwidth  = unit(0.5, "cm"))) +
  scale_fill_gradientn(name = "", colours = pal_emrl, labels = function(x) paste0(x, " %"))

ggsave(p_variations_aet, filename = "~/HubbellGLM-paper/figures/Figure4c_change_richness_aet.png",
       height = 4.04, width = 8.96)


#--------------------------------------------------------
# Calculate the difference for large values of hfp
# (median of variations by zone)
#--------------------------------------------------------
df_Raster %>%
  dplyr::mutate(variation_richnessAET = df_vars$variation_richnessAET,
                diff = richness_curr0 - richness_hfp0,
                diff_pch = 100 * (richness_curr0 - richness_hfp0)/richness_hfp0) %>%
  tidyr::drop_na(variation_richnessAET) %>%
  dplyr::filter(hfp > 30) %>%
  dplyr::group_by(zone) %>%
  dplyr::reframe(low_diff = quantile(diff, 0.25),
                 median_diff = quantile(diff, 0.5),
                 high_diff = quantile(diff, 0.75),
                 median_richness = quantile(diff_pch, 0.5))


#===============================================================================
# Figure 1e - accumulation curves at new sites
#===============================================================================

ylims <- c(-35, 20)
xlims <- c(-90, -31)
df_Raster_filter <- df_Raster %>%
  dplyr::filter(between(x, xlims[1], xlims[2]),
                between(y, ylims[1], ylims[2]))

df_Raster_filter$S_sigma <- exp(predict(fitHubbell,
                                        newdata = df_Raster_filter))/fitHubbell$sigma

pal_emrl <- c("gray30", "gray50", "gray80", "#D3F2A3FF", "#97E196FF", "#6CC08BFF",
              "#02734AFF", "#217A79FF",  "#105965FF", "#03334AFF")

pSigmaBrazil <- ggplot() +
  theme_classic() +
  # Map layer: Use the filtered dataframe
  geom_tile(data = df_Raster_filter,
            aes(x = x, y = y, fill = S_sigma)) + # applied scale factor
  # Border layer: Canada outline
  geom_sf(data = world,
          fill = NA, color = "gray30", size = 0.2) +
  scale_fill_gradientn(name = "S_sigma", colours = pal_emrl, na.value = NA) +
  coord_sf(expand = FALSE)+
  geom_point(data = dataset, aes(x = Longitude, y = Latitude), size = 2.5) +
  geom_point(aes(x = -53.05, y = -8.05), fill = "magenta", size = 3, shape = 21,
             color = "antiquewhite", stroke = 1) +
  geom_point(aes(x = -50.05, y = -24.05), fill = "blue", size = 3, shape = 21,
             color = "antiquewhite", stroke = 1) +
  geom_point(aes(x = -70.05, y = 10.00), fill = "darkorange", size = 3, shape = 21,
             color = "antiquewhite", stroke = 1) +
  ylim(ylims) +
  xlim(xlims)

# Useful function to find the closest point in the raster
find_closest_point <- function(lon, lat, df_Raster) {
  coords <- data.frame(lon = df_Raster$x, lat = df_Raster$y)

  # Compute distances (in meters)
  dists <- distHaversine(matrix(c(lon, lat), nrow=1),    # c(lon, lat)
                         matrix(c(coords$lon, coords$lat), ncol=2))

  closest_idx <- which.min(dists)
  closest_point <- df_Raster[closest_idx, ]
  return(closest_point)
}

# Matrix of variance-covariance for sandwich estimator
vcov_Shared <- vcov_shared(fit = fitHubbell, JaccardSim)
# Magenta point
point1 <- find_closest_point(lon = -53.05, lat = -8.05, df_Raster_filter)
curve1 <- predict_curve(fit = fitHubbell, n = 2000, xnew = point1, .vcov = vcov_Shared)
curve1$point <- "1"
# Blue point
point2 <- find_closest_point(lon = -50.05, lat = -24.05, df_Raster_filter)
curve2 <- predict_curve(fit = fitHubbell, n = 2000, xnew = point2, .vcov = vcov_Shared)
curve2$point <- "2"
# Orange point
point3 <- find_closest_point(lon = -70.05, lat = 10.00, df_Raster_filter)
curve3 <- predict_curve(fit = fitHubbell, n = 2000, xnew = point3, .vcov = vcov_Shared)
curve3$point <- "3"

# Merge them all
df_curves <- data.frame(point = as.factor(rep(c(1,2,3), each = length(curve1$n))),
                        n = unlist(c(curve1$n, curve2$n, curve3$n)),
                        pred = unlist(c(curve1$mean, curve2$mean, curve3$mean)),
                        se = unlist(c(curve1$se, curve2$se, curve3$se)))


pCurves <- ggplot(df_curves, aes(x = n, y = pred, color=point, fill = point)) +
  geom_line() +
  geom_ribbon(aes(ymin = pred - 1.96 * se, ymax = pred + 1.96 * se),
              alpha = 0.20, color = NA) +
  theme_classic() +
  ylab("Predcited number of BINs") +
  scale_color_manual(values = c("magenta", "blue", "darkorange"))+
  scale_fill_manual(values = c("magenta", "blue", "darkorange"))

pSigmaBrazil + pCurves
ggsave(plot = pSigmaBrazil + pCurves, filename = "~/HubbellGLM-paper/figures/Figure1e_Accumulation_Amazon.png",
       width = 12.54, height = 5.66)








