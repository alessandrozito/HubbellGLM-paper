# Descriptive figures and maps for REF_YEAR (2024).
#
# Needs:    01_fit_models.R
# Produces: Fig1b, Fig1c (driver maps), FigS2, FigS3, FigS4, FigS8, FigS10;
#           side maps in figures/extra/; grid_predictions_YYYY.rds

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))
source(file.path(REPO, "analysis", "utils_predictions.R"))

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
library(Matrix)

data_dir <- DATA
ras_dir <- RASTER

# ZONE picks the classification, matching 01_fit_models.R;
# REF_YEAR is the year.used for the one-off exploratory figures and for the headline maps;
# YEARS is the set of years mapped individually.
#
# 2007-2024 is the predictable range: 2000-2006 have no ERA5 on the grid and
# 2025-2026 have no HFP raster, so both ends are structurally empty rather than merely sparse.

res_dir <- MODELS
fig_root <- FIG
REF_YEAR <- as.integer(Sys.getenv("REF_YEAR", "2024"))
YEARS    <- if (nzchar(Sys.getenv("YEARS"))) {
  as.integer(str_split_1(Sys.getenv("YEARS"), ","))
} else 2007:2024
fig <- fig_root                                   # paper figures; side figures go to EXTRA
message("zone: ", ZONE, " | reference year: ", REF_YEAR,
        " | mapped years: ", paste(range(YEARS), collapse = "-"))


#-------- Load the data
model      <- load_hubbell_model()
fitHubbell <- model$fit
MODEL_VARS <- model$vars
dataset    <- model$dataset
ZLV        <- model$zlv
JaccardSim <- readRDS(file.path(data_dir, "JaccardSim_sparse.rds"))
JaccardSim <- JaccardSim[dataset$fieldid, dataset$fieldid]
message("dataset: ", nrow(dataset), " rows")

#-------- Load the global raster, prepared by utils_predictions.R
YEAR <- REF_YEAR
df_Raster <- prepare_grid(YEAR, model)
usable <- df_Raster$usable
# The model covariates a cell must have to be scored at all. prepare_grid()
# already flags these as `usable`, but the column is lost wherever the grid is
# subset, so keep the list itself around too.
NEED <- setdiff(model$vars, "y")

# The original addresses grid coordinates as x/y; keep both names.
df_Raster$x <- df_Raster$Longitude
df_Raster$y <- df_Raster$Latitude
TILE   <- 0.11                   # 0.1 degree cells, small overlap to avoid seams
n_pred <- 2000
q      <- 0.5

#===============================================================================
# Some initial EDA
#===============================================================================
sum(dataset$n)
dim(dataset)
dataset %>%
  dplyr::select(Latitude, Longitude) %>%
  unique()%>%
  dim()

#===============================================================================
# Preliminary step - Extrapolation score (MESS)
#===============================================================================

MESS_COLS <- c("mess_global", "mod_global", "mess_zone", "mod_zone",
               "vars_out_of_range")
if (!all(MESS_COLS %in% names(df_Raster))) {
  message("no MESS columns on the ", YEAR, " raster - run 15_grid_mess.R")
} else {
  message("\n=== cumulative % of each zone below each MESS cut ===")
  print(as.data.frame(df_Raster %>% dplyr::filter(!is.na(mess_zone)) %>%
    dplyr::group_by(zone) %>%
    dplyr::summarise(cells = dplyr::n(),
      `<0`   = round(100 * mean(mess_zone < 0),   1),
      `<-15` = round(100 * mean(mess_zone < -15), 1),
      `<-25` = round(100 * mean(mess_zone < -25), 1),
      .groups = "drop")), row.names = FALSE)

  message("\n=== which variable binds, where MESS < 0 ===")
  print(as.data.frame(df_Raster %>%
    dplyr::filter(!is.na(mess_zone), mess_zone < 0) %>%
    dplyr::count(zone, mod_zone, name = "n") %>%
    dplyr::group_by(zone) %>%
    dplyr::mutate(pct = round(100 * n / sum(n), 1)) %>%
    dplyr::slice_max(n, n = 3) %>% dplyr::ungroup()), row.names = FALSE)
}

#===============================================================================
# Fig 1b: sampling sites
#===============================================================================

world <- ne_countries(scale = "medium", returnclass = "sf")
world <- world[world$name != "Antarctica", ]

zone_colors <- c("Tropical" = "red3", "Dry" = "#D6BD8DFF", "Temperate" = "blue3",
                 "Continental" = "lightblue", "Polar" = "gray50")

p_locs <- ggplot() +
  theme_void() +
  geom_sf(data = world, fill = "gray95", color = "gray30", size = 0.15) +
  theme(panel.grid = element_blank(), aspect.ratio = 0.5,
        legend.position = "inside",
        legend.position.inside = c(0.015, 0.30),
        legend.justification = c(0, 0.5),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size = 12, hjust = 0),
        legend.text = element_text(size = 10),
        legend.key.spacing.y = unit(1, "pt"))+
  coord_sf(expand = FALSE) +
  ylab("Latitude") + xlab("Longitude") +
  geom_point(data = dataset %>% distinct(Longitude, Latitude, zone),
             aes(x = Longitude, y = Latitude, fill = zone), shape = 21, size = 2) +
  scale_fill_manual(name = "Climatic zone", values = zone_colors,
                    guide = guide_legend(override.aes = list(size = 4)))
ggsave(p_locs, filename = file.path(fig, "Fig1b_sites.png"),
       height = 4.04, width = 8.96)

#===============================================================================
# Fig 1c: AET, HFP and zone over South America
#===============================================================================

aet_palette <- c("#290AD8FF", "#264DFFFF", "#3FA0FFFF", "#72D9FFFF", "#AAF7FFFF", "#E0FFFFFF",
                 "#FFFFBFFF", "#FFE099FF", "#FFAD72FF", "#F76D5EFF", "#D82632FF", "#A50021FF")
hfp_palette <- c("white", "antiquewhite1", "#F8B58BFF", "#F59E72FF", "#F2855DFF",
                 "#EF6A4CFF", "red2", "darkred", "black")

ylims <- c(-35, 20); xlims <- c(-90, -31)
df_Raster_Amazon <- df_Raster %>%
  dplyr::filter(between(x, xlims[1], xlims[2]), between(y, ylims[1], ylims[2]))

inset_bar <- theme(
  axis.title             = element_blank(),
  legend.position        = "inside",
  legend.position.inside = c(0.02, 0.22),
  legend.justification   = c(0, 0.5),
  legend.background      = element_blank(),
  legend.title.position  = "top",
  legend.title           = element_text(size = 9),
  legend.text            = element_text(size = 8))

inset_guide <- guide_colourbar(theme = theme(
  legend.key.width  = unit(0.30, "cm"),
  legend.key.height = unit(1.40, "cm")))


pAmazonAET <- ggplot() +
  theme_classic() +
  geom_tile(data = df_Raster_Amazon, aes(x = x, y = y, fill = aet)) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  scale_fill_gradientn(name = "AET\n(mm/year)", colours = rev(aet_palette),
                       na.value = NA, limits = c(0, 140),
                       guide = inset_guide) +
  coord_sf(expand = FALSE) + ylim(ylims) + xlim(xlims) +
  inset_bar
ggsave(plot = pAmazonAET, filename = file.path(fig, "Fig1c_amazon_AET.pdf"),
       height = 4.40, width = 7.72)

pAmazonHFP <- ggplot() +
  theme_classic() +
  geom_tile(data = df_Raster_Amazon, aes(x = x, y = y, fill = hfp)) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  scale_fill_gradientn(name = "HFP", colours = hfp_palette,
                       na.value = NA, limits = c(0, 50),
                       guide = inset_guide) +
  coord_sf(expand = FALSE) + ylim(ylims) + xlim(xlims) +
  inset_bar
ggsave(plot = pAmazonHFP, filename = file.path(fig, "Fig1c_amazon_HFP.pdf"),
       height = 4.40, width = 7.72)

zone_keys <- tibble::tibble(x = 150, y = -80, zone = names(zone_colors))

pAmazonZone <- ggplot() +
  theme_classic() +
  geom_tile(data = dplyr::bind_rows(tidyr::drop_na(df_Raster_Amazon, zone), zone_keys),
            aes(x = x, y = y, fill = zone)) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  scale_fill_manual(name = "Climatic zone", values = zone_colors,
                    limits = names(zone_colors), drop = FALSE, na.translate = FALSE,
                    guide = guide_legend(override.aes = list(size = 4))) +
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE)

ggsave(plot = pAmazonZone, filename = file.path(fig, "Fig1c_amazon_zone.pdf"),
       height = 4.40, width = 7.72)


#===============================================================================
# Fig S2: accumulation-curve shapes by sigma
#===============================================================================
n <- 1:1000
eta <- c(-1, 0, 1)
alpha <- exp(eta)
sigma <- c(-0.5, -0.25, 0, 0.1, 0.25)
df <- data.frame()
for (a in alpha) for (s in sigma) {
  mu_n <- HubbellGLM:::polyseries_mean(size = n, alpha = rep(a, length(n)), sigma = s)
  df <- rbind(df, data.frame(n, alpha = a, sigma = factor(s), mu_n))
}
df$eta <- log(df$alpha)
df$eta <- factor(paste("eta ==", df$eta), levels = paste("eta ==", eta))

pAccum <- ggplot(df, aes(x = n, y = mu_n, color = sigma)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ eta, labeller = label_parsed) +
  scale_color_manual(values = c("#C5D358", "#318B27", "#FF0000", "#4B71D0", "#000080")) +
  labs(y = "Species richness (mean)", color = expression(sigma)) +
  theme_bw(base_size = 14)
ggsave(pAccum, width = 10.35, height = 3.74,
       filename = file.path(fig, "FigS2_accumulation_shapes.pdf"))

#===============================================================================
# Fig S3: Jaccard similarity vs geographic distance
#===============================================================================
set.seed(1)
NPAIR <- 3e6
ii <- sample.int(nrow(dataset), NPAIR, replace = TRUE)
jj <- sample.int(nrow(dataset), NPAIR, replace = TRUE)
keep <- ii < jj
ii <- ii[keep]; jj <- jj[keep]
sim <- JaccardSim[cbind(ii, jj)]        # sparse indexing; unstored pairs are 0
xy  <- as.matrix(dataset[, c("Longitude", "Latitude")])
gd  <- distHaversine(xy[ii, ], xy[jj, ]) / 1000

df_plot_dist <- data.frame(geo_dist_km = gd, similarity = as.numeric(sim))
message("Fig S3: ", format(nrow(df_plot_dist), big.mark = ","),
        " sampled pairs, ", round(100 * mean(df_plot_dist$similarity > 0), 1),
        "% with shared species")

p_dist <- ggplot(df_plot_dist, aes(x = geo_dist_km, y = similarity)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "Density\n(n. of points)", labels = comma) +
  scale_x_continuous(labels = comma) +
  labs(x = "Geographic distance (km)", y = "Fraction of shared species") +
  theme_bw(base_size = 14) +
  facet_wrap(~"Relationship between similarity and distance")
ggsave(p_dist, width = 9.34, height = 4.55,
       filename = file.path(fig, "FigS3_geographic_similarity.pdf"))

#===============================================================================
# Fig S10: covariate and MESS maps
#===============================================================================

# One grammar for all four panels, matching the CV partition maps in
# 10_cv_figure.R: grey land base, country outlines on top, no axes, bold
# left-aligned title, cropped to the sampled latitudes.
cov_map <- function(d, mapping, scale, title, ...) {
  ggplot() +
    geom_sf(data = world, fill = "grey96", colour = "grey75", linewidth = 0.1) +
    geom_tile(data = d, mapping = mapping, ...) +
    geom_sf(data = world, fill = NA, colour = "grey75", linewidth = 0.1) +
    scale +
    coord_sf(expand = FALSE, ylim = c(-56, 84)) +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0,
                                    margin = margin(b = 3)),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.key.width = unit(0.30, "cm"),
          plot.background = element_rect(fill = "white", colour = NA))
}
save_panel <- function(p, f) ggsave(p, filename = file.path(EXTRA, paste0(f, ".png")),
                                    height = 3.1, width = 7.72, dpi = 150)

p_aet <- cov_map(df_Raster, aes(x, y, fill = aet),
                 scale_fill_gradientn(name = "mm/year", colours = rev(aet_palette)),
                 "Actual evapotranspiration")
p_hfp <- cov_map(df_Raster %>% dplyr::filter(!is.na(hfp)), aes(x, y, fill = hfp),
                 scale_fill_gradientn(name = "HFP", colours = hfp_palette),
                 "Human footprint", width = TILE, height = TILE, na.rm = TRUE)
p_zone <- cov_map(df_Raster %>% dplyr::filter(!is.na(zone)), aes(x, y, fill = zone),
                  scale_fill_manual(name = NULL, values = zone_colors),
                  "K\u00f6ppen-Geiger climatic zone")
save_panel(p_aet,  "covariate_AET")
save_panel(p_hfp,  "covariate_HFP")
save_panel(p_zone, "covariate_zone")

MESS_BR  <- c(-Inf, -50, -25, -15, -5, 0, 25, 50, Inf)
MESS_LAB <- c("< -50", "-50 to -25", "-25 to -15", "-15 to -5", "-5 to 0",
              "0 to 25", "25 to 50", "> 50")
MESS_COL <- c("< -50" = "#67000d", "-50 to -25" = "#b2182b", "-25 to -15" = "#ef6548",
              "-15 to -5" = "#fdae61", "-5 to 0" = "#fee0b6", "0 to 25" = "#d9e6f2",
              "25 to 50" = "#92c5de", "> 50" = "#2166ac")

if ("mess_zone" %in% names(df_Raster)) {
  df_mess <- df_Raster %>%
    dplyr::mutate(mess_bin = cut(mess_zone, MESS_BR, labels = MESS_LAB,
                                 right = FALSE)) %>%
    tidyr::drop_na(mess_bin)

  p_mess <- cov_map(df_mess, aes(x, y, fill = mess_bin),
                    scale_fill_manual(name = NULL, values = MESS_COL, drop = FALSE,
                                      guide = guide_legend(reverse = TRUE)),
                    "MESS extrapolation")
  save_panel(p_mess, "covariate_MESS")

  p_mess_filter <- cov_map(df_mess %>% dplyr::filter(mess_zone < -25),
                           aes(x, y, fill = mod_zone), NULL,
                           "Most extrapolated covariate (MESS < -25)")
  save_panel(p_mess_filter, "covariate_MESS_variables")

  message("\n=== % of land in each MESS bin ===")
  print(as.data.frame(df_mess %>% dplyr::count(mess_bin, .drop = FALSE) %>%
    dplyr::mutate(pct_of_land = round(100 * n / sum(n), 2))), row.names = FALSE)

  ggsave(plot = (p_aet + p_hfp)/(p_mess + p_zone),
         filename = file.path(fig, "FigS10_covariate_maps.png"),
         height = 6.6, width = 14.5, dpi = 150)

} else {
  message("no mess_zone on the ", YEAR, " raster - run preprocessing/15_grid_mess.R")
}


#===============================================================================
# Figs S8 and S4: collection history; covariates vs diversity
#===============================================================================

#---- Fig S8: collection history by site
df <- dataset %>%
  dplyr::select(fieldid, site = site_code, zone, realm,
                collection_start_date, collection_end_date) %>%
  distinct() %>%
  mutate(collection_start_date = as.Date(collection_start_date),
         collection_end_date   = as.Date(collection_end_date)) %>%
  rowwise() %>%
  dplyr::mutate(date = list(seq(collection_start_date, collection_end_date, by = "day"))) %>%
  ungroup() %>%
  unnest(date) %>%
  dplyr::filter(year(date) <= YEAR) %>%
  dplyr::mutate(min_start = min(collection_start_date, na.rm = TRUE),
                date_diff = as.integer(date - min_start) + 1)

site_lat <- dataset %>%
  dplyr::group_by(site = site_code) %>%
  dplyr::summarise(Latitude = mean(Latitude), .groups = "drop")
site_levels <- site_lat %>% arrange(desc(Latitude)) %>% pull(site)

lat_ticks <- tibble(lat = seq(80, -60, by = -20)) %>%
  dplyr::filter(lat <= max(site_lat$Latitude), lat >= min(site_lat$Latitude)) %>%
  dplyr::mutate(
    site  = site_lat$site[vapply(lat, function(z)
              which.min(abs(site_lat$Latitude - z)), integer(1))],
    label = paste0(abs(lat), "\u00b0",
                   ifelse(lat > 0, "N", ifelse(lat < 0, "S", "")))) %>%
  dplyr::filter(!duplicated(site))

year_ticks <- tibble(year = seq(min(year(df$date)), max(year(df$date)))) %>%
  dplyr::mutate(date_diff = as.integer(as.Date(paste0(year, "-01-01")) -
                                       df$min_start[1]) + 1)

df$site <- factor(df$site, levels = site_levels)

p_sites <- ggplot(df, aes(x = date_diff, y = site)) +
  geom_point(aes(color = zone), size = 0.4) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank()) +
  scale_y_discrete(limits = rev(levels(df$site)),
                   breaks = lat_ticks$site, labels = lat_ticks$label,
                   expand = expansion(mult = 0.02)) +
  scale_x_continuous(breaks = year_ticks$date_diff, labels = year_ticks$year,
                     name = "Collection year") +
  ylab("Sampling site (ordered by latitude)") +
  scale_color_manual(name = "Climatic zone", values = zone_colors,
                     guide = guide_legend(override.aes = list(size = 4)))

ggsave(p_sites, filename = file.path(fig, "FigS8_collection_history.png"),
       height = 9.09, width = 8.00)

#---- Panels b and c
dataset$Ssigma <- HubbellGLM:::inv_polyseries(dataset$y, dataset$n, model$sigma) / model$sigma

YLAB <- expression(log(S[sigma]))

ZONE_SCALE <- scale_color_manual(name = "Climatic zone", values = zone_colors,
                                 guide = guide_legend(override.aes = list(size = 3)))

pAlpha_hfp <- ggplot(dataset, aes(x = hfp, y = log(Ssigma), color = zone)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", colour = "black") +
  facet_grid(~zone) +
  theme_bw() + labs(x = "HFP", y = YLAB) +
  ZONE_SCALE

pAlpha_aet <- ggplot(dataset, aes(x = aet, y = log(Ssigma))) +
  geom_point(alpha = 0.4, aes(color = zone), size = 0.8) +
  geom_smooth(method = "lm", colour = "black") +
  theme_bw() + labs(x = "AET", y = YLAB) +
  facet_grid(~"Annual actual evapotranspiration") +
  ZONE_SCALE +
  stat_cor(method = "pearson", label.x.npc = "center",
           label.y.npc = "bottom", size = 4)


#---- Panel d - other covariates against diversity
avail <- c(Latitude = "Latitude", `Vapor pressure deficit` = "vpd",
           `Pot. evapotranspiration` = "pet",
           Wetlands = "wetlands", Temperature = "temperature_2m",
           Precipitation = "total_precipitation_sum",
           `Wind speed` = "wind_speed", `Rel. Humidity` = "relative_humidity")
avail <- avail[avail %in% names(dataset)]
# Panel d wraps into this many columns; its height below is derived from the
# resulting number of rows, so adding a covariate cannot silently squash it.
FACET_NCOL <- 4

dataset_long <- dataset %>%
  dplyr::select(Ssigma, zone, all_of(unname(avail))) %>%
  rename(!!!setNames(unname(avail), names(avail))) %>%
  pivot_longer(cols = all_of(names(avail)), names_to = "Variable", values_to = "Value")

p_environ_covariates <- ggplot(dataset_long, aes(x = Value, y = log(Ssigma))) +
  geom_point(alpha = 0.4, aes(color = zone), size = 0.4) +
  stat_cor(method = "pearson", label.x.npc = "left",
           label.y.npc = "bottom", size = 3) +
  facet_wrap(~ Variable, scales = "free_x", ncol = FACET_NCOL) +
  ZONE_SCALE +
  theme_bw() +
  labs(x = NULL, y = YLAB)

d_rows <- ceiling(length(avail) / FACET_NCOL)
p_joind <- pAlpha_aet / pAlpha_hfp / p_environ_covariates +
  plot_layout(guides = "collect", heights = c(1, 1, d_rows)) &
  theme(legend.position = "right")

ggsave(p_joind, filename = file.path(fig, "FigS4_covariates_diversity.pdf"),
       width = 8.6, height = 3.1 * (2 + d_rows))

#===============================================================================
# Side maps (figures/extra/): global diversity indices
#===============================================================================

# Indices and counterfactuals, both from utils_predictions.R so this file and
# any other year use identical definitions.
df_vars <- dplyr::bind_cols(
  tibble::tibble(Latitude = df_Raster$y, Longitude = df_Raster$x),
  predict_indices(df_Raster, model, n_pred = n_pred, q = q),
  predict_variations(df_Raster, model))

pal_emrl <- c("antiquewhite2", "#D3F2A3FF", "#97E196FF", "#6CC08BFF",
              "#02734AFF", "#217A79FF", "#105965FF", "#03334AFF", "black")

map_index <- function(var, cols, bh = 6.5) {
  ggplot() +
    theme_void() +
    geom_tile(data = df_vars %>% drop_na(all_of(var)),
              aes(x = Longitude, y = Latitude, fill = .data[[var]]),
              width = TILE, height = TILE) +
    geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
          aspect.ratio = 0.5) +
    coord_sf(expand = FALSE) +
    guides(color = "none",
           fill = guide_colourbar(barheight = unit(bh, "cm"),
                                  barwidth = unit(0.5, "cm"))) +
    scale_fill_gradientn(name = "", colours = cols)
}

pal_shannon <- c("#1E8E99FF", "#51C3CCFF", "#B2FCFFFF", "#CCFEFFFF", "#E5FFFFFF",
                 "#FFE5CCFF", "#FFAD65FF", "#FF8E32FF", "#CC5800FF", "#993F00FF")
pal_simpson <- c("white", "#E5FFFFFF", "#B2F2FFFF", "#7FD4FFFF", "#65BFFFFF",
                 "#4CA5FFFF", "#003FFFFF")
palette_tsallis <- c("#F6C866", "#F2AB67", "#EF8F6B", "#ED7470", "#BF6E97",
                     "#926AC2", "#6667EE", "#4959C7", "#2D4A9F", "#173C78")

ggsave(map_index("S_sigma", c("gray60", pal_emrl)),
       filename = file.path(EXTRA, "map_S_sigma.png"), width = 5.81, height = 2.89, dpi = 200)
ggsave(map_index("alpha", c("gray60", pal_emrl)),
       filename = file.path(EXTRA, "map_alpha.png"), width = 5.81, height = 2.89, dpi = 200)
ggsave(map_index("Shannon", pal_shannon, bh = 7),
       filename = file.path(EXTRA, "map_shannon.png"), width = 5.81, height = 2.89, dpi = 200)
ggsave(map_index("Simpson", pal_simpson),
       filename = file.path(EXTRA, "map_simpson.png"), width = 5.81, height = 2.89, dpi = 200)
ggsave(map_index("Tsallis", palette_tsallis),
       filename = file.path(EXTRA, "map_tsallis.png"), width = 5.81, height = 2.89, dpi = 200)

#===============================================================================
# Side maps (figures/extra/): raw scenario maps; Fig 3 is 06_scenario_maps.R
#===============================================================================

#----------- REPRODUCE FIGURE 4b

p_variations_hfp <- ggplot() +
  theme_void() +
  geom_tile(data = df_vars %>% drop_na(var_hfp10),
            aes(x = Longitude, y = Latitude, fill = var_hfp10),
            width = TILE, height = TILE) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  theme(panel.grid = element_blank(), axis.title = element_blank(), aspect.ratio = 0.5) +
  coord_sf(expand = FALSE) +
  guides(color = "none",
         fill = guide_colourbar(barheight = unit(7, "cm"), barwidth = unit(0.5, "cm"))) +
  scale_fill_gradientn(
    name = "",
    colours = c("darkred", "red2", "#F2855DFF", "antiquewhite1", "lightblue1", "#290AD8FF"),
    limits = c(-5.3, 2.8),
    values = scales::rescale(c(seq(-5.3, 0, length.out = 4)[1:3], 0, 0.2, 2.8), from = c(-5.3, 2.8)),
    oob = scales::squish, breaks = seq(-5, 2.5, by = 2.5),
    labels = function(x) paste0(x, " %"))
ggsave(p_variations_hfp, filename = file.path(EXTRA, "scenario_hfp5_v1.png"),
       height = 4.04, width = 8.96)

p_variations_hfp_2 <- ggplot() +
  theme_void() +
  geom_tile(data = df_vars %>% drop_na(var_hfp5),
            aes(x = Longitude, y = Latitude, fill = var_hfp5),
            width = TILE, height = TILE) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  theme(panel.grid = element_blank(), axis.title = element_blank(), aspect.ratio = 0.5) +
  coord_sf(expand = FALSE) +
  guides(color = "none",
         fill = guide_colourbar(barheight = unit(7, "cm"), barwidth = unit(0.5, "cm"))) +
  scale_fill_gradientn(
    name = "", colours = c("darkred", "red2", "antiquewhite1", "#3FA0FFFF", "#264DFFFF"),
    limits = c(-8, 8), oob = scales::squish, breaks = round(seq(-8, 8, by = 2)),
    labels = function(x) paste0(x, " %"))
ggsave(p_variations_hfp_2, filename = file.path(EXTRA, "scenario_hfp5_v2.png"),
       height = 4.04, width = 8.96)

#----------- REPRODUCE FIGURE 4a
maxabs0 <- max(abs(df_vars$var_hfp0), na.rm = TRUE)

p_hfp_zero <- ggplot() +
  theme_void() +
  geom_tile(data = df_vars %>% drop_na(var_hfp0),
            aes(x = Longitude, y = Latitude, fill = var_hfp0),
            width = TILE, height = TILE) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  coord_sf(expand = FALSE) +
  theme(panel.grid = element_blank(), axis.title = element_blank(), aspect.ratio = 0.5) +
  guides(color = "none",
         fill = guide_colourbar(barheight = unit(7, "cm"), barwidth = unit(0.5, "cm"))) +
  scale_fill_gradientn(
    name = "",
    colours = c("#290AD8FF", "#264DFFFF", "#3FA0FFFF", "antiquewhite1",
                "#6CC08BFF", "#02734AFF", "#03334AFF"),
    limits = c(-maxabs0, maxabs0), oob = scales::squish,
    breaks = round(seq(-60, 60, by = 20)), labels = function(x) paste0(x, " %"))
ggsave(plot = p_hfp_zero, filename = file.path(EXTRA, "scenario_hfp0_raw.png"),
       height = 4.04, width = 8.96)

#----------- REPRODUCE FIGURE 4c
pal_emrl <- c("antiquewhite1", "#E5FFB6", "#D3F2A3FF", "#97E196FF", "#6CC08BFF",
              "#02734AFF", "#217A79FF", "#105965FF", "#03334AFF")

p_variations_aet <- ggplot() +
  theme_void() +
  geom_tile(data = df_vars %>% drop_na(var_aet10),
            aes(x = Longitude, y = Latitude, fill = var_aet10),
            width = TILE, height = TILE) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  theme(panel.grid = element_blank(), axis.title = element_blank(), aspect.ratio = 0.5) +
  coord_sf(expand = FALSE) +
  guides(color = "none",
         fill = guide_colourbar(barheight = unit(7, "cm"), barwidth = unit(0.5, "cm"))) +
  scale_fill_gradientn(name = "", colours = pal_emrl, labels = function(x) paste0(x, " %"))
ggsave(p_variations_aet, filename = file.path(EXTRA, "scenario_aet10_raw.png"),
       height = 4.04, width = 8.96)

#--------------------------------------------------------
# Difference for large values of hfp (median of variations by zone)
#--------------------------------------------------------
# "High HFP" is defined PER ZONE, as a quantile of the zone's own HFP
# distribution, rather than by a single global cut such as hfp > 30. A fixed
# cut is not comparable across zones and, worse, lands outside the data: all 39
# Polar grid cells above 30 exceed the largest Polar HFP ever observed in the
# sample (25.06), so the old table's Polar row was pure extrapolation - the
# model reported a ~350% richness change from a coefficient it had no support
# for out there.
#
# The quantile is taken from the SAMPLE, not the grid. The grid's 80th
# percentile for Polar is 0.3, because the zone is almost entirely wilderness,
# so "the top 20% of Polar cells" would mean "any human presence at all". The
# sample quantile is anchored to where the model actually has observations.
# HIGH_HFP=grid switches to the grid distribution; HIGH_Q sets the quantile.
#
# The quantile alone still does not guarantee support - 94 Polar cells above
# the grid q80 exceed the sample maximum - so cells are additionally capped at
# each zone's observed HFP range. Both bounds are reported in the table.
HIGH_Q   <- as.numeric(Sys.getenv("HIGH_Q", "0.8"))
HIGH_SRC <- Sys.getenv("HIGH_HFP", "sample")

thr <- (if (identical(HIGH_SRC, "grid")) df_Raster else dataset) %>%
  dplyr::filter(!is.na(zone), !is.na(hfp)) %>%
  dplyr::group_by(zone) %>%
  dplyr::summarise(hfp_thr = quantile(hfp, HIGH_Q), .groups = "drop")
sup <- dataset %>%
  dplyr::group_by(zone) %>%
  dplyr::summarise(hfp_sup = max(hfp), .groups = "drop")

tab_hfp <- df_Raster %>%
  dplyr::mutate(var_aet10 = df_vars$var_aet10,
                diff = df_vars$richness_raw - df_vars$richness_hfp0,
                diff_pch = 100 * (df_vars$richness_raw - df_vars$richness_hfp0) /
                           df_vars$richness_hfp0) %>%
  tidyr::drop_na(var_aet10) %>%
  dplyr::left_join(thr, by = "zone") %>%
  dplyr::left_join(sup, by = "zone") %>%
  dplyr::filter(hfp >= hfp_thr, hfp <= hfp_sup) %>%
  dplyr::group_by(zone) %>%
  dplyr::reframe(hfp_threshold = first(hfp_thr),
                 hfp_support_max = first(hfp_sup),
                 cells = dplyr::n(),
                 low_diff = quantile(diff, 0.25),
                 median_diff = quantile(diff, 0.5),
                 high_diff = quantile(diff, 0.75),
                 median_richness = quantile(diff_pch, 0.5))
write_tsv(tab_hfp, file.path(res_dir, "hfp_high_by_zone.tsv"))
message("\n=== high-HFP cells (", HIGH_SRC, " q", HIGH_Q * 100,
        " per zone, capped at observed support), richness lost to human pressure ===")
print(as.data.frame(tab_hfp), row.names = FALSE, digits = 4)

#===============================================================================
# Side figure (figures/extra/): original accumulation curves; Fig 1c is 04_accumulation_curves.R
#===============================================================================

ylims <- c(-35, 20); xlims <- c(-90, -31)
df_Raster_filter <- df_Raster %>%
  dplyr::filter(between(x, xlims[1], xlims[2]), between(y, ylims[1], ylims[2]))
df_Raster_filter$S_sigma <- exp(predict(fitHubbell, newdata = df_Raster_filter)) / fitHubbell$sigma

pal_emrl <- c("gray25", "gray50", "gray80", "#D3F2A3FF", "#97E196FF", "#6CC08BFF",
              "#02734AFF", "#217A79FF", "#105965FF", "#03334AFF")
#pal_emrl <- c("white", "antiquewhite", "#D3F2A3FF", "#97E196FF", "#6CC08BFF",
#                            "#02734AFF", "#217A79FF", "#105965FF", "#03334AFF")

pSigmaBrazil <- ggplot() +
  theme_classic() +
  geom_tile(data = df_Raster_filter, aes(x = x, y = y, fill = S_sigma)) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  scale_fill_gradientn(
    name = expression(S[sigma]), colours = pal_emrl, na.value = NA,
    breaks = c(20, 40, 60, 80),
    guide = guide_colourbar(
      theme = theme(legend.title = element_text(size = 17),
                    legend.title.position = "top",
                    legend.text  = element_text(size = 11),
                    legend.key.width  = unit(0.42, "cm"),
                    legend.key.height = unit(2.7, "cm"),
                    legend.ticks = element_line(colour = "white"),
                    legend.frame = element_blank()))) +
  coord_sf(expand = FALSE) +
  geom_point(data = dataset %>% distinct(Longitude, Latitude),
             aes(x = Longitude, y = Latitude), size = 1.2) +
  geom_point(aes(x = -53.05, y = -8.05), fill = "magenta", size = 3, shape = 21,
             color = "antiquewhite", stroke = 1) +
  geom_point(aes(x = -50.05, y = -24.05), fill = "blue", size = 3, shape = 21,
             color = "antiquewhite", stroke = 1) +
  geom_point(aes(x = -70.05, y = 10.00), fill = "darkorange", size = 3, shape = 21,
             color = "antiquewhite", stroke = 1) +
  ylim(ylims) + xlim(xlims) +
  theme(axis.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.105, 0.34),
        legend.justification = c(0.5, 0.5),
        legend.background = element_blank())

pSigmaBrazil

# Restrict to cells the model can actually score, or predict_curve() builds a
# model matrix with an NA row and returns nothing usable.
find_closest_point <- function(lon, lat, df_Raster) {
  ok <- complete.cases(df_Raster[, NEED])
  d2 <- df_Raster[ok, ]
  dists <- distHaversine(matrix(c(lon, lat), nrow = 1),
                         as.matrix(d2[, c("x", "y")]))
  d2[which.min(dists), ]
}

vcov_Shared <- model$vcov
pts <- list(c(-53.05, -8.05), c(-50.05, -24.05), c(-70.05, 10.00))
curves <- lapply(seq_along(pts), function(k) {
  p  <- find_closest_point(pts[[k]][1], pts[[k]][2], df_Raster_filter)
  cv <- predict_curve(fit = fitHubbell, n = 2000, xnew = p, .vcov = vcov_Shared)
  data.frame(point = as.factor(k),
             n    = as.numeric(cv$n),
             pred = as.numeric(unlist(cv$mean)),
             se   = as.numeric(unlist(cv$se)))
})
df_curves <- bind_rows(curves)

pCurves <- ggplot(df_curves, aes(x = n, y = pred, color = point, fill = point)) +
  geom_line() +
  geom_ribbon(aes(ymin = pred - 1.96 * se, ymax = pred + 1.96 * se),
              alpha = 0.20, color = NA) +
  theme_classic() +
  ylab("Predicted number of BINs") +
  scale_color_manual(values = c("magenta", "blue", "darkorange")) +
  scale_fill_manual(values = c("magenta", "blue", "darkorange"))

ggsave(plot = pSigmaBrazil + pCurves,
       filename = file.path(EXTRA, "accumulation_curves_original.png"),
       width = 12.54, height = 5.66)

saveRDS(df_vars, file.path(res_dir, sprintf("grid_predictions_%d.rds", YEAR)))
message("\nfigures in ", fig)
message("grid predictions in ", res_dir, "/grid_predictions_", YEAR, ".rds")
