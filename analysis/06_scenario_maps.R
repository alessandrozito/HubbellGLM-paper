# Figure 3: modelled % change in richness under HFP = 0, HFP + 5 and AET + 10%.
# Cells with MESS < MASK_CUT are grey; non-significant cells are stippled.
#
# Needs:    01_fit_models.R
# Produces: Fig3_scenario_maps.png, grid_change_YYYY.rds
#           (RERUN=0 replots from the cached grid)

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))
source(file.path(REPO, "analysis", "utils_predictions.R"))

suppressMessages({library(tidyverse); library(sf); library(rnaturalearth)
                  library(patchwork); library(scales)})

data_dir <- DATA
res_dir  <- MODELS
fig_dir  <- FIG
YEAR     <- REF_YEAR
RERUN    <- Sys.getenv("RERUN", "1") != "0"
MASK_CUT <- as.numeric(Sys.getenv("MASK_CUT", "-25"))
CHUNK    <- as.integer(Sys.getenv("CHUNK", "150000"))
TILE     <- 0.11        # grid is 0.1 deg; a hair over, so tiles never gap
NSAMP    <- 2000        # the sample size every map is drawn at
NMAX     <- 2000        # curves run out to the same n
FIG_W    <- as.numeric(Sys.getenv("FIG_W", "17"))
FIG_H    <- as.numeric(Sys.getenv("FIG_H", "12.6"))
MAP_SH   <- as.numeric(Sys.getenv("MAP_SH", "1.78"))  # map : panels width share
PANEL_H  <- as.numeric(Sys.getenv("PANEL_H", "0.74")) # right panels, as a
                                                      # fraction of the row
STIP     <- 2.0         # stipple block, degrees
STIP_F   <- 0.5         # block is stippled when this fraction of its cells
                        # have an interval covering zero

D_f  <- file.path(res_dir, sprintf("grid_change_%d.rds",  YEAR))
CV_f <- file.path(res_dir, sprintf("figure4_curves_%d.rds", YEAR))

#---- The scenarios -----------------------------------------------------------
# The hfp + 0.01 nudge matches utils_predictions.R: hfp = 0 is a boundary of the
# fitted range, so "current" is read just inside it.

SCEN <- list(
  hfp0 = list(
    title = "Variation in richness when HFP is zero everywhere",
    sub   = "modelled richness with HFP set to zero everywhere",
    alt_lab = "HFP = 0",
    base  = function(d) d,
    alt   = function(d) dplyr::mutate(d, hfp = 0),
    pal   = c("#290AD8FF", "#264DFFFF", "#3FA0FFFF", "antiquewhite1",
              "#6CC08BFF", "#02734AFF", "#03334AFF")),
  hfp5 = list(
    title = "Variation in richness from 5 units increase in HFP",
    sub   = "modelled richness with 5 HFP units added everywhere",
    alt_lab = "HFP +5 units",
    base  = function(d) dplyr::mutate(d, hfp = pmin(hfp + 0.01, 50)),
    alt   = function(d) dplyr::mutate(d, hfp = pmin(hfp + 0.01 + 5, 50)),
    pal   = c("darkred", "red2", "antiquewhite1", "#3FA0FFFF", "#264DFFFF"),
    lim   = c(-8, 8), breaks = seq(-8, 8, by = 2)),
  aet10 = list(
    title = "Variation in richness from 10% increase in AET",
    sub   = "modelled richness with AET raised by 10% of its local value",
    alt_lab = "AET +10%",
    base  = function(d) dplyr::mutate(d, aet = aet + 0.01),
    alt   = function(d) dplyr::mutate(d, aet = (aet + 0.01) * 1.1),
    seq   = TRUE,
    pal   = c("antiquewhite1", "#E5FFB6", "#D3F2A3FF", "#97E196FF", "#6CC08BFF",
              "#02734AFF", "#217A79FF", "#105965FF", "#03334AFF")))

#---- The four locations the curves are drawn at ------------------------------
# Chosen to span the contrast the maps are making, not to be representative:
# three human-dominated cells in different zones and one boreal cell where the
# model has nothing to say. The tropical cell is on the Rondonia deforestation
# arc rather than in the untouched central basin - at HFP 0.2 the counterfactual
# moves central Amazonia by 0.1%, which is a true number and an empty panel.

# Row order IS the panel order, and it runs west to east, so the four panels
# read left to right in the same order as the four points on the map.
SITES <- tibble::tribble(
  ~site,      ~label,           ~lon,   ~lat,
  "boreal",   "Boreal Canada",-100.0,   55.0,
  "rondonia", "Rondonia",      -63.0,  -10.0,
  "europe",   "Central Europe", 10.5,   49.0,
  "mekong",   "Indochina",     105.0,   15.0)

# Four hues the maps do not use: the palettes own red, orange, green, blue and
# cream, which leaves violet, magenta, turquoise and ink.
SITE_COL <- c(boreal   = "#7C4DFF",   # electric violet
              rondonia = "#FF2D95",   # hot pink
              europe   = "#111111",   # ink
              mekong   = "#00D0C0")   # turquoise

reuse <- !RERUN && file.exists(D_f) && file.exists(CV_f)
if (reuse) {
  D <- readRDS(D_f); CV <- readRDS(CV_f)
  # a cache written for a different scenario set is not a cache
  reuse <- all(paste0(names(SCEN), "_aest") %in% names(D)) &&
           all(names(SCEN) %in% unique(CV$curves$scenario)) &&
           "abs" %in% names(CV$curves)
  if (reuse) { SITES <- CV$sites; CURVES <- CV$curves } else
    message("cached scenarios do not match SCEN - recomputing")
}

#==============================================================================
# Computation
#==============================================================================

if (!reuse) {

model <- load_hubbell_model()
g     <- prepare_grid(YEAR, model, n_pred = NSAMP)
g     <- g[g$usable, ]
message("cells with every model covariate: ", format(nrow(g), big.mark = ","))

# The delta method, not Monte Carlo. Two linkinv calls per draw at ~3.5 s each on
# this grid is an hour per scenario; the change is a smooth function of beta, so
#
#   R  = 100 (mu1/mu0 - 1)
#   dR/dbeta = 100 r [ (g'(e1)/mu1) X1 - (g'(e0)/mu0) X0 ],   r = mu1/mu0
#   Var(R_i) = grad_i' V grad_i
#
# V is the Jaccard variance-covariance, so these intervals carry the same
# shared-species inflation as every other number in the paper.

fam   <- hubbell(sigma = model$sigma)
Terms <- delete.response(terms(model$fit))
beta  <- coef(model$fit)
V     <- model$vcov
stopifnot(identical(names(beta), colnames(V)))

mm <- function(d) model.matrix(Terms, data = d,
                               xlev = model$fit$xlevels)[, names(beta), drop = FALSE]

# Both scales come out of the same four evaluations. The absolute change is
#   A = mu1 - mu0,  dA/dbeta = g'(e1) X1 - g'(e0) X0
# which is the gradient the percentage version divides through by mu0.
change_se <- function(g, f_base, f_alt) {
  n  <- nrow(g)
  est <- se <- aest <- ase <- rep(NA_real_, n)
  for (from in seq(1, n, by = CHUNK)) {
    i  <- from:min(from + CHUNK - 1, n)
    gi <- g[i, ]
    X0 <- mm(f_base(gi)); X1 <- mm(f_alt(gi))
    e0 <- as.numeric(X0 %*% beta);  e1 <- as.numeric(X1 %*% beta)
    sz <- gi$n
    m0 <- fam$linkinv(e0, sz, model$sigma); m1 <- fam$linkinv(e1, sz, model$sigma)
    d0 <- fam$mu.eta(e0, sz, model$sigma);  d1 <- fam$mu.eta(e1, sz, model$sigma)
    r  <- m1 / m0
    est[i]  <- 100 * (r - 1)
    aest[i] <- m1 - m0
    G  <- 100 * r * ((d1 / m1) * X1 - (d0 / m0) * X0)
    se[i] <- sqrt(pmax(rowSums((G %*% V) * G), 0))
    A  <- d1 * X1 - d0 * X0
    ase[i] <- sqrt(pmax(rowSums((A %*% V) * A), 0))
    rm(X0, X1, G, A); gc(verbose = FALSE)
  }
  tibble(est = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se,
         aest = aest, ase = ase, alo = aest - 1.96 * ase, ahi = aest + 1.96 * ase)
}

D <- tibble(Longitude = g$Longitude, Latitude = g$Latitude, zone = g$zone,
            hfp = g$hfp, aet = g$aet,
            mess_zone = if ("mess_zone" %in% names(g)) g$mess_zone else NA_real_,
            mod_zone  = if ("mod_zone"  %in% names(g)) g$mod_zone  else NA_character_)

for (s in names(SCEN)) {
  message("\n", s, " ...")
  t0 <- Sys.time()
  r  <- change_se(g, SCEN[[s]]$base, SCEN[[s]]$alt)
  D[[paste0(s, "_est")]] <- r$est
  D[[paste0(s, "_se")]]  <- r$se
  D[[paste0(s, "_lo")]]  <- r$lo
  D[[paste0(s, "_hi")]]  <- r$hi
  # A cell whose interval covers zero is a cell where the model cannot say which
  # way richness moves. Those are stippled on the map, not hidden.
  D[[paste0(s, "_sig")]] <- !is.na(r$lo) & (r$lo > 0 | r$hi < 0)
  D[[paste0(s, "_aest")]] <- r$aest
  D[[paste0(s, "_ase")]]  <- r$ase
  D[[paste0(s, "_asig")]] <- !is.na(r$alo) & (r$alo > 0 | r$ahi < 0)
  message(sprintf("  median %+.2f%%  median SE %.2f  |  %.1f%% of cells have a CI clear of zero  [%.0f s]",
                  median(r$est, na.rm = TRUE), median(r$se, na.rm = TRUE),
                  100 * mean(!is.na(r$lo) & (r$lo > 0 | r$hi < 0)),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
saveRDS(D, D_f)

#---- Accumulation curves at the four sites -----------------------------------
# Same delta method, evaluated along n instead of across cells. The curve band
# is g'(eta) * se(eta), which is what HubbellGLM::predict_curve draws; the band
# on the % change needs the full gradient because the two curves share beta.

D$masked <- if (all(is.na(D$mess_zone))) FALSE else
  is.na(D$mess_zone) | D$mess_zone < MASK_CUT

# Snap each target to the nearest cell the model is allowed to speak about. A
# site that landed on a masked cell would be a curve drawn from extrapolation.
snap <- function(lon, lat) {
  d2 <- ((D$Longitude - lon) * cos(lat * pi / 180))^2 + (D$Latitude - lat)^2
  which.min(ifelse(D$masked, Inf, d2))
}
SITES$row <- vapply(seq_len(nrow(SITES)),
                    function(i) snap(SITES$lon[i], SITES$lat[i]), integer(1))
SITES$Longitude <- D$Longitude[SITES$row]
SITES$Latitude  <- D$Latitude [SITES$row]
SITES$zone      <- as.character(D$zone[SITES$row])
SITES$hfp       <- D$hfp[SITES$row]
SITES$panel_lab <- sprintf("%s (%s)", SITES$label, SITES$zone)

ns_grid <- unique(round(seq(1, NMAX, length.out = 140)))

curve_one <- function(row, f_base, f_alt) {
  X0 <- mm(f_base(row)); X1 <- mm(f_alt(row))
  e0 <- as.numeric(X0 %*% beta); e1 <- as.numeric(X1 %*% beta)
  k  <- length(ns_grid)
  m0 <- fam$linkinv(rep(e0, k), ns_grid, model$sigma)
  m1 <- fam$linkinv(rep(e1, k), ns_grid, model$sigma)
  d0 <- fam$mu.eta(rep(e0, k), ns_grid, model$sigma)
  d1 <- fam$mu.eta(rep(e1, k), ns_grid, model$sigma)
  s0 <- sqrt(drop(X0 %*% V %*% t(X0)))
  s1 <- sqrt(drop(X1 %*% V %*% t(X1)))
  r  <- m1 / m0
  G  <- 100 * r * (matrix(d1 / m1, ncol = 1) %*% X1 -
                   matrix(d0 / m0, ncol = 1) %*% X0)
  A  <- matrix(d1, ncol = 1) %*% X1 - matrix(d0, ncol = 1) %*% X0
  tibble(n = ns_grid, mu0 = m0, mu1 = m1, se0 = d0 * s0, se1 = d1 * s1,
         rel = 100 * (r - 1), rel_se = sqrt(pmax(rowSums((G %*% V) * G), 0)),
         abs = m1 - m0,       abs_se = sqrt(pmax(rowSums((A %*% V) * A), 0)))
}

CURVES <- map_dfr(names(SCEN), function(s) {
  map_dfr(seq_len(nrow(SITES)), function(i) {
    row <- g[SITES$row[i], ]
    curve_one(row, SCEN[[s]]$base, SCEN[[s]]$alt) %>%
      mutate(scenario = s, site = SITES$site[i], .before = 1)
  })
})
saveRDS(list(sites = SITES, curves = CURVES, mask_cut = MASK_CUT), CV_f)

} else {
  message("RERUN=0: replotting from ", basename(D_f), " and ", basename(CV_f))
}

#---- The extrapolation mask --------------------------------------------------

if (all(is.na(D$mess_zone))) {
  message("\nNO MESS on the raster - run 15_grid_mess.R. Nothing will be masked, ",
          "and this figure would then show extrapolation as if it were supported.")
  D$masked <- FALSE
} else {
  D$masked <- is.na(D$mess_zone) | D$mess_zone < MASK_CUT
}
PCT_MASK <- 100 * mean(D$masked)
message("\nmask at MESS < ", MASK_CUT, ": ", format(sum(D$masked), big.mark = ","),
        " of ", format(nrow(D), big.mark = ","), " cells (",
        round(PCT_MASK, 2), "%)")
# pct BEFORE the count: summarise() evaluates in order, so computing
# `masked = sum(masked)` first makes the later mean(masked) read the sum.
print(as.data.frame(D %>% group_by(zone) %>%
  summarise(cells = n(), pct_masked = round(100 * mean(masked), 1),
            masked = sum(masked), .groups = "drop") %>%
  dplyr::select(zone, cells, masked, pct_masked)), row.names = FALSE)

message("\n=== the four focal cells ===")
print(as.data.frame(SITES %>% dplyr::select(label, zone, Longitude, Latitude, hfp)),
      row.names = FALSE, digits = 4)

#---- Zone summaries, on the retained cells only ------------------------------

zone_tab <- map_dfr(names(SCEN), function(s) {
  D %>% dplyr::filter(!masked) %>% group_by(zone) %>%
    summarise(scenario = s, cells = n(),
              median_pct = median(.data[[paste0(s, "_est")]], na.rm = TRUE),
              q25 = quantile(.data[[paste0(s, "_est")]], .25, na.rm = TRUE),
              q75 = quantile(.data[[paste0(s, "_est")]], .75, na.rm = TRUE),
              median_se = median(.data[[paste0(s, "_se")]], na.rm = TRUE),
              median_bins = median(.data[[paste0(s, "_aest")]], na.rm = TRUE),
              pct_cells_CI_clear_of_zero =
                round(100 * mean(.data[[paste0(s, "_sig")]], na.rm = TRUE), 1),
              .groups = "drop")
})
write_tsv(zone_tab, file.path(res_dir, "figure4_zone_summary.tsv"))
message("\n=== by zone, masked cells excluded ===")
print(as.data.frame(zone_tab), row.names = FALSE, digits = 3)

#==============================================================================
# Left column: the maps
#==============================================================================

world <- tryCatch(ne_countries(scale = "medium", returnclass = "sf") %>%
                    dplyr::filter(name != "Antarctica"), error = function(e) NULL)
GREY <- "grey86"

# The panels take a mode: "pct" is the relative change the figure reports and the
# published figures used, "abs" the change in expected BINs at n = NSAMP. The two
# differ only in which column they read and how the axis is labelled, so one set
# of functions serves both.
COL <- function(s, what, mode)
  paste0(s, "_", if (mode == "abs") "a" else "", what)

# Limits: fixed where the published figure fixed them, otherwise symmetric about
# zero at the 99th percentile of |change| - the tails are squished into the end
# colour rather than dropped, and one extreme cell cannot flatten the map.
fill_limits <- function(s, v, mode) {
  if (mode == "pct" && !is.null(SCEN[[s]]$lim)) return(SCEN[[s]]$lim)
  if (isTRUE(SCEN[[s]]$seq)) c(0, quantile(v, 0.99, na.rm = TRUE))
  else { q <- quantile(abs(v), 0.99, na.rm = TRUE); c(-q, q) }
}

# Tags by hand: plot_annotation(tag_levels) also tags the facet strips.
TAG <- list(hfp0 = "a", hfp5 = "b", aet10 = "c")

# Three layers, because a plain dot disappears: on the maps these sit on dark
# green, on red and on stipple, and no single colour survives all three. The
# black rim reads against a light cell, the white ring against a dark one.
SITE_PTS <- list(
  geom_point(data = SITES, aes(Longitude, Latitude),
             colour = "black", size = 5.4),
  geom_point(data = SITES, aes(Longitude, Latitude),
             fill = SITE_COL[SITES$site], shape = 21, size = 4.4,
             colour = "white", stroke = 1.5))

# IPCC-style stippling: a 2-degree block is marked when most of its cells have an
# interval covering zero. One dot per cell would be 1.5M dots and unreadable.
stipple_pts <- function(s, mode) {
  D %>% dplyr::filter(!masked) %>%
    dplyr::mutate(gx = floor(Longitude / STIP) * STIP + STIP / 2,
                  gy = floor(Latitude  / STIP) * STIP + STIP / 2) %>%
    group_by(gx, gy) %>%
    summarise(f = mean(!.data[[COL(s, "sig", mode)]]), cells = n(),
              .groups = "drop") %>%
    dplyr::filter(cells >= 20, f > STIP_F)
}

# No titles, no subtitles, tight margins: the map is the panel. What the grey and
# the stipple mean goes in the caption, which the script prints when it finishes.
map_theme <- list(
  coord_sf(expand = FALSE, ylim = c(-56, 84)),
  theme_void(base_size = 10),
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(2, 4, 10, 2),
        # the title is indented past the tag, so the two share the line
        plot.title = element_text(size = 15, hjust = 0, face = "bold",
                                  margin = margin(l = 34, b = 4)),
        plot.tag.position = c(0.004, 0.975),
        plot.tag = element_text(size = 21, face = "bold"),
        legend.position = "right",
        # the bar sits against the map: the default box spacing is 11pt of air
        legend.box.spacing = unit(2, "pt"),
        legend.margin = margin(0, 0, 0, 0)))

map_change <- function(s, mode) {
  keep <- D %>% dplyr::filter(!masked)
  v    <- keep[[COL(s, "est", mode)]]
  lim  <- fill_limits(s, v, mode)
  stp  <- stipple_pts(s, mode)
  ggplot() +
    # Masked cells are drawn, in grey, so the reader sees WHERE the model was not
    # asked rather than an empty ocean.
    geom_tile(data = D %>% dplyr::filter(masked),
              aes(Longitude, Latitude), fill = GREY, width = TILE, height = TILE) +
    geom_tile(data = keep, aes(Longitude, Latitude, fill = .data[[COL(s, "est", mode)]]),
              width = TILE, height = TILE) +
    geom_point(data = stp, aes(gx, gy), size = 0.12, colour = "grey20",
               alpha = 0.65, shape = 16) +
    {if (!is.null(world)) geom_sf(data = world, fill = NA, colour = "grey40",
                                  linewidth = 0.13, inherit.aes = FALSE)} +
    SITE_PTS +
    scale_fill_gradientn(
      name = if (mode == "abs") "BINs" else NULL,
      colours = SCEN[[s]]$pal, limits = lim, oob = scales::squish,
      breaks = if (mode == "pct" && !is.null(SCEN[[s]]$breaks))
                 SCEN[[s]]$breaks else waiver(),
      labels = if (mode == "pct") function(x) paste0(x, " %") else waiver(),
      # the bar's size lives on the guide, not the plot theme: a plot-level
      # legend.key.height does not constrain a colourbar in ggplot2 4.0
      guide = guide_colourbar(theme = theme(
        legend.key.width  = unit(0.55, "cm"),
        legend.key.height = unit(4.8, "cm"),
        legend.text = element_text(size = 12)))) +
    map_theme + labs(tag = TAG[[s]], title = SCEN[[s]]$title)
}

#==============================================================================
# Right column: the same change, along n, at the four marked cells
#==============================================================================

# Facet order follows SITES, not the alphabet, and the label breaks over two
# lines: a 1.3-inch strip truncates "Central Europe (Temperate)" silently.
CURVES$site <- factor(CURVES$site, levels = SITES$site)
FACET_LAB   <- setNames(sprintf("%s\n(%s)", SITES$label, SITES$zone), SITES$site)

curve_theme <- theme_bw(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey96", colour = "grey70"),
        strip.text = element_text(size = 8, lineheight = 1.05,
                                  margin = margin(2, 2, 2, 2)),
        panel.spacing.x = unit(7, "pt"),
        # rotated, so "2000" and the next panel's "0" stop colliding
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.margin = margin(2, 8, 10, 12),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold",
                                  margin = margin(b = 5)),
        axis.title = element_text(size = 9))

change_panel <- function(s, mode) {
  cv <- CURVES %>% dplyr::filter(scenario == s)
  y  <- if (mode == "abs") "abs" else "rel"
  se <- if (mode == "abs") "abs_se" else "rel_se"
  ggplot(cv, aes(n, .data[[y]], colour = site, fill = site)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3, linetype = 2) +
    geom_ribbon(aes(ymin = .data[[y]] - 1.96 * .data[[se]],
                    ymax = .data[[y]] + 1.96 * .data[[se]]),
                colour = NA, alpha = 0.22) +
    geom_line(linewidth = 0.7) +
    # One y axis across the four: the panels are meant to be compared, and the
    # flat purple line in the last one is the comparison.
    facet_wrap(~ site, nrow = 1, labeller = labeller(site = FACET_LAB)) +
    scale_colour_manual(values = SITE_COL, guide = "none") +
    scale_fill_manual(values = SITE_COL, guide = "none") +
    scale_x_continuous(breaks = c(0, 1000, 2000),
                       expand = expansion(mult = c(0.03, 0.05))) +
    labs(x = "Number of individuals sampled (n)",
         title = "Variation at four locations",
         y = if (mode == "abs") "Change in expected BINs"
             else "% change in expected BINs") +
    curve_theme
}

#==============================================================================
# Composition
#==============================================================================

# wrap_elements() is what makes the maps large. patchwork aligns panels across a
# row, so without it the map's panel row is the facet plot's panel row - and
# coord_sf then pins the map's WIDTH to that height times 2.5, leaving an inch of
# white on either side. Freed from the alignment, the map fills its own cell:
# 6.8 x 2.7 in becomes 9.9 x 4.0 in on the same canvas.
#
# FIG_W / FIG_H / MAP_SH then have to agree with the map's aspect, or the slack
# comes straight back: map width = MAP_SH/(MAP_SH+1) * FIG_W - legend, and
# FIG_H/3 has to match that width / 2.5.
compose <- function(mode) {
  rows <- lapply(names(SCEN), function(s)
    (wrap_elements(full = map_change(s, mode)) |
       # spacers above and below shrink the panels without shrinking the row,
       # which is what gives the map the visual weight
       wrap_plots(plot_spacer(), change_panel(s, mode), plot_spacer(), ncol = 1,
                  heights = c((1 - PANEL_H) / 2, PANEL_H, (1 - PANEL_H) / 2))) +
      plot_layout(widths = c(MAP_SH, 1)))
  wrap_plots(rows, ncol = 1)
}

# The figure is the relative change. compose("abs") draws the same figure in
# BINs per NSAMP individuals - the absolute change is still computed and carried
# through D, CURVES and median_bins in the zone table, it is just not the figure.
#
# PNG only: 4.5M tiles across three panels is a vector file no journal system
# will open. The SI maps are written the same way.
ggsave(file.path(fig_dir, "Fig3_scenario_maps.png"), compose("pct"),
       width = FIG_W, height = FIG_H, dpi = 300, limitsize = FALSE)
message("wrote Fig3_scenario_maps.png")

#---- SI: the full standard-error maps ----------------------------------------

map_uncert <- function(s) {
  keep <- D %>% dplyr::filter(!masked)
  hw   <- 1.96 * keep[[paste0(s, "_se")]]
  ggplot() +
    geom_tile(data = D %>% dplyr::filter(masked),
              aes(Longitude, Latitude), fill = GREY, width = TILE, height = TILE) +
    geom_tile(data = keep, aes(Longitude, Latitude, fill = 1.96 * .data[[paste0(s, "_se")]]),
              width = TILE, height = TILE) +
    SITE_PTS +
    scale_fill_viridis_c(
      name = "95% CI\nhalf-width\n(% points)", option = "magma", direction = -1,
      trans = "sqrt", limits = c(0, quantile(hw, 0.99, na.rm = TRUE)),
      oob = scales::squish,
      guide = guide_colourbar(theme = theme(legend.key.width  = unit(0.30, "cm"),
                                            legend.key.height = unit(2.2, "cm")))) +
    {if (!is.null(world)) geom_sf(data = world, fill = NA, colour = "grey40",
                                  linewidth = 0.13, inherit.aes = FALSE)} +
    map_theme + labs(tag = NULL)
}

ggsave(file.path(EXTRA, "Fig3_scenario_uncertainty.png"),
       wrap_plots(lapply(names(SCEN), map_uncert), ncol = 1) +
         plot_annotation(tag_levels = "a"),
       width = 9.5, height = 12, dpi = 200, limitsize = FALSE)
message("wrote Fig3_scenario_uncertainty.png (extra) and figure4_zone_summary.tsv")

#---- The caption, which now carries what the panels no longer say ------------

message("\n--- caption ---\n",
  sprintf(paste0(
    "Figure 3 | Modelled response of arthropod richness to human pressure and to\n",
    "productivity, %d. a,c,e, modelled change in richness under three counterfactuals:\n",
    "%s (a), %s (c) and %s (e), evaluated at n = %d individuals. Grey: cells outside\n",
    "the sampled environmental envelope (MESS < %g; %.1f%% of land), where the model\n",
    "would be extrapolating. Stippling: 2-degree blocks in which most cells have a 95%%\n",
    "interval covering zero. Coloured points mark the four cells of b,d,f. b,d,f, the\n",
    "same change along the number of individuals sampled, at those four cells; bands\n",
    "are 95%% intervals from the Jaccard variance-covariance. These are model\n",
    "projections under a counterfactual, not observed change over time."),
    YEAR, SCEN$hfp0$alt_lab, SCEN$hfp5$alt_lab, SCEN$aet10$alt_lab,
    NSAMP, MASK_CUT, PCT_MASK))
