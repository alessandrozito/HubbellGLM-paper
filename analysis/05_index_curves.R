# Figure 2: regression diversity and predicted richness vs AET, by zone and HFP.
# Intervals are Monte Carlo draws from the Jaccard covariance.
#
# Needs:    01_fit_models.R
# Produces: Fig2_index_curves.png, FigS7_index_curves_all.png, FigS11_index_maps.png

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))
source(file.path(REPO, "analysis", "utils_predictions.R"))

suppressMessages({library(tidyverse); library(sf); library(rnaturalearth)
                  library(patchwork); library(MASS, exclude = c("select"))})

data_dir <- DATA
res_dir  <- MODELS
fig_dir  <- FIG
YEAR     <- REF_YEAR
RERUN    <- Sys.getenv("RERUN", "1") != "0"
NSIMS    <- as.integer(Sys.getenv("NSIMS", "500"))
NSAMP    <- 2000        # every index is read at this sample size
MASK_CUT <- as.numeric(Sys.getenv("MASK_CUT", "-25"))
TILE     <- 0.11
q        <- 0.5
ZLV      <- c("Continental", "Dry", "Polar", "Temperate", "Tropical")

# Which indices go where. MAIN is the top block; MAPS is the bottom row - the
# four that were dropped from Figure 1, S_sigma having stayed there.
MAIN <- c("S_sigma", "Richness")
MAPS <- c("alpha", "Shannon", "Simpson", "Tsallis")
ALL  <- c("S_sigma", "Richness", "alpha", "Shannon", "Simpson", "Tsallis")

curves_f <- file.path(res_dir, sprintf("figure4_index_curves_%d.rds", YEAR))
grid_f   <- file.path(res_dir, sprintf("grid_indices_%d.rds", YEAR))

model <- load_hubbell_model()
fit   <- model$fit
V     <- model$vcov

#==============================================================================
# a - the response curves
#==============================================================================

# One draw of beta gives one whole curve, so the interval is the pointwise
# quantile of NSIMS curves rather than a band assembled per point.
index_ci <- function(fit, newdata, .vcov, indices, size_n, n_sims, q) {
  Terms <- delete.response(terms(fit))
  Xpred <- model.matrix(Terms, data = newdata, xlev = fit$xlevels)
  beta_sims <- MASS::mvrnorm(n_sims, mu = coef(fit), Sigma = .vcov)
  sim <- setNames(lapply(indices, function(i)
    matrix(NA_real_, nrow(newdata), n_sims)), indices)
  dummy <- fit
  need_rich  <- any(c("Richness", "alpha", "Shannon", "Simpson", "Tsallis") %in% indices)
  need_alpha <- any(c("alpha", "Shannon", "Simpson", "Tsallis") %in% indices)
  pb <- txtProgressBar(min = 0, max = n_sims, style = 3)
  for (i in seq_len(n_sims)) {
    b <- beta_sims[i, ]
    if ("S_sigma" %in% indices)
      sim[["S_sigma"]][, i] <- exp(as.vector(Xpred %*% b)) / fit$sigma
    if (need_rich) {
      dummy$coefficients <- b
      r <- predict(dummy, newdata = newdata, type = "response")
      if ("Richness" %in% indices) sim[["Richness"]][, i] <- r
      if (need_alpha) {
        a <- vapply(r, function(x) tryCatch(
          HubbellGLM:::inv_mean_dirichlet_process(mu_target = x, size = size_n),
          error = function(e) NA_real_), numeric(1))
        if ("alpha"   %in% indices) sim[["alpha"]][, i]   <- a
        if ("Shannon" %in% indices) sim[["Shannon"]][, i] <- digamma(a + 1) - digamma(1)
        if ("Simpson" %in% indices) sim[["Simpson"]][, i] <- 1 / (a + 1)
        if ("Tsallis" %in% indices) sim[["Tsallis"]][, i] <- 1 / (q - 1) * (1 - a * beta(a, q))
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  setNames(lapply(indices, function(i) tibble(
    Lower = apply(sim[[i]], 1, quantile, .025, na.rm = TRUE),
    Upper = apply(sim[[i]], 1, quantile, .975, na.rm = TRUE))), indices)
}

if (!RERUN && file.exists(curves_f)) {
  message("RERUN=0: curves from ", basename(curves_f))
  df_long <- readRDS(curves_f)
} else {
  g <- readRDS(file.path(RASTER, sprintf("Raster_%d_full.rds", YEAR)))
  gz <- if (identical(ZONE, "zone_kgc") && "zone_kgc" %in% names(g)) g$zone_kgc else g$zone
  g$zone <- factor(gz, levels = ZLV)

  # Each zone is drawn only across the AET range it actually occupies: the
  # curves are counterfactual in hfp, not in climate.
  zone_env <- g %>%
    dplyr::filter(!is.na(zone), !is.na(realm), !is.na(aet)) %>%
    group_by(zone) %>%
    summarise(Latitude_mean = mean(Latitude, na.rm = TRUE),
              aet_low  = quantile(aet, 0.05, na.rm = TRUE),
              aet_high = quantile(aet, 0.95, na.rm = TRUE),
              realm    = names(sort(table(realm), decreasing = TRUE))[1],
              .groups  = "drop")
  print(as.data.frame(zone_env), row.names = FALSE, digits = 4)

  d <- model$dataset
  df_out <- expand.grid(zone = ZLV, hfp = c(0, 25, 50), aet = 0:120) %>%
    as_tibble() %>% dplyr::mutate(zone = as.character(zone)) %>%
    dplyr::left_join(zone_env %>% dplyr::mutate(zone = as.character(zone)), by = "zone") %>%
    dplyr::filter(aet >= aet_low, aet <= aet_high) %>%
    dplyr::transmute(zone = factor(zone, levels = ZLV),
                     realm = factor(realm, levels = levels(d$realm)),
                     hfp, aet, Latitude = Latitude_mean,
                     week_adj = 26, collection_days = 7,
                     resid_temperature_lat = 0, wetlands = 0, wind_speed = 0,
                     n = NSAMP, y = 100)
  for (v in grep("_anml$", model$vars, value = TRUE)) df_out[[v]] <- 0
  # Controls added to the model after this figure was first drawn. All enter
  # linearly, so holding them fixed shifts every curve by the same amount and
  # leaves the aet and hfp contrasts untouched. The island terms take the sample
  # median: 0 is not a possible landmass and log10(0) is -Inf.
  if ("year" %in% model$vars) df_out$year <- YEAR
  for (v in c("landmass_km2", "dist_mainland_km"))
    if (v %in% model$vars) df_out[[v]] <- median(d[[v]], na.rm = TRUE)

  # MUST be a base data.frame. predict.HubbellGLM does
  #   family(object)$linkinv(pred, newdata[, object$name_size], sigma)
  # and tbl[, "n"] returns a one-column TIBBLE, so linkinv receives a list:
  #   "Not compatible with requested type: [type=list; target=double]".
  df_out <- as.data.frame(df_out)
  message("design rows: ", nrow(df_out), "   sims: ", NSIMS)

  df_out$eta      <- predict(fit, newdata = df_out)
  df_out$Richness <- predict(fit, newdata = df_out, type = "response")
  df_out$alpha    <- HubbellGLM:::inv_mean_dirichlet_process(mu_target = df_out$Richness,
                                                             size = df_out$n)
  df_out$S_sigma  <- exp(df_out$eta) / fit$sigma
  df_out$Shannon  <- digamma(df_out$alpha + 1) - digamma(1)
  df_out$Simpson  <- 1 / (df_out$alpha + 1)
  df_out$Tsallis  <- 1 / (q - 1) * (1 - df_out$alpha * beta(df_out$alpha, q))

  set.seed(10)
  ci <- index_ci(fit, df_out, V, ALL, size_n = NSAMP, n_sims = NSIMS, q = q)

  df_long <- map_dfr(ALL, function(nm)
    tibble(zone = df_out$zone, hfp = df_out$hfp, aet = df_out$aet,
           Index = nm, Value = df_out[[nm]],
           Lower = ci[[nm]]$Lower, Upper = ci[[nm]]$Upper)) %>%
    # Simpson is a probability; a simulated upper limit above 1 is a Monte Carlo
    # artefact, not a claim.
    dplyr::mutate(Upper = ifelse(Index == "Simpson" & Upper > 1, 1, Upper),
                  Index = factor(Index, levels = ALL),
                  hfp_factor = factor(hfp, levels = c(0, 25, 50),
                                      labels = c("HFP 0", "HFP 25", "HFP 50")))
  saveRDS(df_long, curves_f)
}

#---- The index strips, as maths ----------------------------------------------
# label_parsed wants plotmath, so each label is one expression: atop() stacks
# the name over the formula.

IDX_MATH <- c(
  S_sigma  = paste0('atop(displaystyle(atop(bold("Regression"), bold("diversity"))), ',
                    'displaystyle(S[sigma](eta) == e^eta/sigma))'),
  Richness = paste0('atop(displaystyle(atop(bold("Predicted"), bold("richness"))), ',
                    'displaystyle(hat(y)[n]))'),
  alpha    = 'atop(bold("Alpha"), hat(alpha)[n])',
  Shannon  = 'atop(bold("Shannon"), psi(hat(alpha)[n]+1) - psi(1))',
  Simpson  = 'atop(bold("Simpson"), 1/(hat(alpha)[n]+1))',
  Tsallis  = 'atop(bold("Tsallis"), (1 - hat(alpha)[n]*B(hat(alpha)[n], q))/(q-1))')

HFP_COL <- c("HFP 0" = "antiquewhite3", "HFP 25" = "#CC0000", "HFP 50" = "#330000")

variation_panel <- function(which_idx, base_size = 13, strip_size = base_size + 3) {
  dd <- df_long %>% dplyr::filter(Index %in% which_idx) %>%
    dplyr::mutate(Index = factor(IDX_MATH[as.character(Index)],
                                 levels = IDX_MATH[which_idx]))
  ggplot(dd, aes(aet, Value, colour = hfp_factor, fill = hfp_factor)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, linewidth = 0.2,
                linetype = "dashed") +
    geom_line(linewidth = 1) +
    facet_grid(Index ~ zone, scales = "free",
               labeller = labeller(Index = label_parsed)) +
    scale_x_continuous(n.breaks = 4) +
    scale_colour_manual(values = HFP_COL, name = "Human Footprint\n(HFP)") +
    scale_fill_manual(values = HFP_COL, name = "Human Footprint\n(HFP)") +
    labs(x = "Actual evapotranspiration (AET)", y = "Diversity index value") +
    theme_classic(base_size = base_size) +
    theme(strip.background.x = element_blank(),
          strip.background.y = element_rect(fill = "grey95", colour = "grey60"),
          # "Continental" at base_size + 3 is wider than a tight panel, and
          # ggplot clips strip text to the strip without warning
          strip.text.x = element_text(size = strip_size, face = "bold"),
          strip.clip = "off",
          strip.text.y = element_text(size = base_size + 3, angle = 0),
          axis.text  = element_text(colour = "black", size = base_size),
          axis.title = element_text(size = base_size + 2),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          panel.spacing = unit(7, "pt"),
          # rotated: at this panel width "60" of one zone met "0" of the next
          axis.text.x = element_text(angle = 45, hjust = 1, size = base_size - 1),
          legend.position = "right",
          legend.title = element_text(size = base_size + 1),
          legend.text  = element_text(size = base_size),
          legend.key.height = unit(0.75, "cm"),
          plot.tag = element_text(size = 21, face = "bold"),
          # 0.971, not 0.99: it sets a level with the maps block's b in the
          # side-by-side version, whose plots carry a title above the panel
          plot.tag.position = c(0.004, 0.971),
          plot.margin = margin(2, 4, 8, 2))
}

#==============================================================================
# b-e - the global rasters
#==============================================================================

df_map <- if (!RERUN && file.exists(grid_f)) readRDS(grid_f) else NULL
# a cache written before the mask existed has no mess_zone to mask with
if (!is.null(df_map) && "mess_zone" %in% names(df_map)) {
  message("RERUN=0: grid indices from ", basename(grid_f))
} else {
  gg <- prepare_grid(YEAR, model, n_pred = NSAMP)
  gg <- gg[gg$usable, ]
  t0 <- Sys.time()
  df_map <- dplyr::bind_cols(
    tibble(Longitude = gg$Longitude, Latitude = gg$Latitude,
           mess_zone = if ("mess_zone" %in% names(gg)) gg$mess_zone else NA_real_),
    predict_indices(gg, model, n_pred = NSAMP, q = q))
  message("indices on ", format(nrow(df_map), big.mark = ","), " cells in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  saveRDS(df_map, grid_f)
}

world <- tryCatch(ne_countries(scale = "medium", returnclass = "sf") %>%
                    dplyr::filter(name != "Antarctica"), error = function(e) NULL)

# Same mask as the counterfactual figure: cells whose environment is outside the
# sampled envelope are drawn grey, not coloured. An index map is a prediction
# like any other, and these are the cells the model was never asked about.
GREY <- "grey86"
df_map$masked <- if (all(is.na(df_map$mess_zone))) FALSE else
  is.na(df_map$mess_zone) | df_map$mess_zone < MASK_CUT
if (all(is.na(df_map$mess_zone)))
  message("NO MESS on the raster - run 15_grid_mess.R; nothing will be masked")
message(sprintf("mask at MESS < %g: %s of %s cells (%.1f%%)", MASK_CUT,
                format(sum(df_map$masked), big.mark = ","),
                format(nrow(df_map), big.mark = ","), 100 * mean(df_map$masked)))

# The palettes are the published ones, one per index.
PAL <- list(
  alpha   = c("gray60", "antiquewhite2", "#D3F2A3FF", "#97E196FF", "#6CC08BFF",
              "#02734AFF", "#217A79FF", "#105965FF", "#03334AFF", "black"),
  Shannon = c("#1E8E99FF", "#51C3CCFF", "#B2FCFFFF", "#CCFEFFFF", "#E5FFFFFF",
              "#FFE5CCFF", "#FFAD65FF", "#FF8E32FF", "#CC5800FF", "#993F00FF"),
  # the published ramp starts at pure white, which on a white page reads as
  # missing data rather than as a low value
  Simpson = c("#F2FBFF", "#E5FFFFFF", "#B2F2FFFF", "#7FD4FFFF", "#65BFFFFF",
              "#4CA5FFFF", "#003FFFFF"),
  Tsallis = c("#F6C866", "#F2AB67", "#EF8F6B", "#ED7470", "#BF6E97",
              "#926AC2", "#6667EE", "#4959C7", "#2D4A9F", "#173C78"))
MAP_TITLE <- c(
  alpha   = 'atop(bold("Alpha"), hat(alpha)[n])',
  Shannon = 'atop(bold("Shannon"), psi(hat(alpha)[n]+1) - psi(1))',
  Simpson = 'atop(bold("Simpson"), 1/(hat(alpha)[n]+1))',
  Tsallis = 'atop(bold("Tsallis")~(q == 0.5), (1 - hat(alpha)[n]*B(hat(alpha)[n], q))/(q-1))')
TRANS <- c(alpha = "sqrt", Shannon = "identity", Simpson = "log10",
           Tsallis = "identity")
# under sqrt and log the automatic breaks bunch up at one end and the labels
# collide, so the transformed scales get their own
BREAKS <- list(alpha = c(50, 250, 750, 1500), Shannon = waiver(),
               Simpson = c(0.001, 0.01, 0.1), Tsallis = waiver())
# No panel tag: the maps are their own figure now, not panel b of a composite.

map_index <- function(var) {
  ggplot() +
    geom_tile(data = df_map %>% dplyr::filter(masked),
              aes(Longitude, Latitude), fill = GREY,
              width = TILE, height = TILE) +
    geom_tile(data = df_map %>% dplyr::filter(!masked) %>% tidyr::drop_na(all_of(var)),
              aes(Longitude, Latitude, fill = .data[[var]]),
              width = TILE, height = TILE) +
    {if (!is.null(world)) geom_sf(data = world, fill = NA, colour = "grey40",
                                  linewidth = 0.12, inherit.aes = FALSE)} +
    scale_fill_gradientn(
      name = NULL, colours = PAL[[var]], na.value = NA,
      transform = TRANS[[var]], breaks = BREAKS[[var]], n.breaks = 3,
      labels = function(x) formatC(scales::label_number(drop0trailing = TRUE)(x),
                                   width = 5),
      # sizing a colourbar has to happen on the guide: a plot-level
      # legend.key.height does not reach it in ggplot2 4.0
      guide = guide_colourbar(theme = theme(
        legend.key.width  = unit(0.42, "cm"),
        legend.key.height = unit(3.4, "cm"),
        legend.text = element_text(size = 11)))) +
    coord_sf(expand = FALSE, ylim = c(-56, 84)) +
    theme_void(base_size = 10) +
    labs(title = parse(text = MAP_TITLE[[var]])[[1]]) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 3)),
          plot.tag = element_text(size = 21, face = "bold"),
          plot.tag.position = c(0.004, 0.99),
          plot.margin = margin(2, 2, 2, 2),
          legend.position = "right",
          legend.box.spacing = unit(3, "pt"),
          legend.margin = margin(0, 0, 0, 0))
}

#==============================================================================
# Composition
#==============================================================================

FIG_W  <- as.numeric(Sys.getenv("FIG_W", "15"))
FIG_H  <- as.numeric(Sys.getenv("FIG_H", "5.6"))
MAP_W  <- as.numeric(Sys.getenv("MAP_W", "14"))
MAP_H  <- as.numeric(Sys.getenv("MAP_H", "6.4"))

# The curves and the maps are two separate figures. They used to share a canvas,
# and the maps' fixed aspect ratio then drove the whole layout: it set the canvas
# height, and the curve panels had to live in whatever width was left. Apart they
# can each be sized for what they show.

# Fig 2: the two headline indices against AET and human pressure.
fig <- variation_panel(MAIN)
ggsave(file.path(fig_dir, "Fig2_index_curves.png"), fig,
       width = FIG_W, height = FIG_H, dpi = 300, limitsize = FALSE)
message("wrote Fig2_index_curves.png")

# Fig S11: the same fitted eta read through four indices, globally.
# wrap_elements keeps each map out of patchwork's panel alignment; without it a
# fixed-aspect map is pinned to the row height of its neighbour and shrinks.
maps_grid <- wrap_plots(
  lapply(MAPS, function(v) wrap_elements(full = map_index(v))), nrow = 2)
ggsave(file.path(fig_dir, "FigS11_index_maps.png"), maps_grid,
       width = MAP_W, height = MAP_H, dpi = 300, limitsize = FALSE)
message("wrote FigS11_index_maps.png")

#---- SI: the same curves, all six indices ------------------------------------

si <- variation_panel(ALL, base_size = 10)
ggsave(file.path(fig_dir, "FigS7_index_curves_all.png"), si,
       width = 11.5, height = 11, dpi = 200)
message("wrote FigS7_index_curves_all.png")

#---- What the curves say, in numbers -----------------------------------------

message("\n=== index at the top of each zone's AET range, HFP 0 vs HFP 50 ===")
print(as.data.frame(
  df_long %>% dplyr::filter(Index %in% MAIN, hfp %in% c(0, 50)) %>%
    group_by(zone, Index) %>% dplyr::filter(aet == max(aet)) %>%
    dplyr::select(zone, Index, hfp, Value) %>%
    tidyr::pivot_wider(names_from = hfp, values_from = Value,
                       names_prefix = "hfp") %>%
    dplyr::mutate(pct_lost = round(100 * (hfp50 - hfp0) / hfp0, 1),
                  across(c(hfp0, hfp50), ~ signif(.x, 4))) %>%
    dplyr::arrange(Index, zone)), row.names = FALSE)
