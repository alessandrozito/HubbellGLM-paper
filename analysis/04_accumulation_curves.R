# Figure 1c, accumulation curves: three tropical cells at low/medium/high HFP,
# matched on latitude, AET, wind and realm.
#
# Needs:    01_fit_models.R
# Produces: Fig1c_accumulation_curves.png (counterfactual version in figures/extra/)

source(file.path(path.expand("~/HubbellGLM-paper"), "preprocessing", "paths.R"))
suppressMessages({library(tidyverse); library(HubbellGLM); library(splines)
                  library(rnaturalearth); library(rnaturalearthdata); library(sf)
                  library(geosphere); library(patchwork)})
source(file.path(REPO, "analysis", "utils_predictions.R"))

YEAR  <- 2024
xlims <- c(-90, -31); ylims <- c(-35, 20)
HFP_TARGETS <- c(0, 22, 45)          # low / medium / high
AET_REF     <- 97                    # ~ the tropical median in this window
MIN_SEP_KM  <- 900                   # far apart, so the Keynote funnels don't cross
LAT_BAND    <- c(-14, 6)             # candidate band
LAT_SPREAD  <- 4                     # max latitude difference WITHIN the chosen triple.
                                     # Longitude is not a model covariate, so the points
                                     # can be spread east-west for free; latitude enters
                                     # through ns(Latitude, 6) and must be held close.
AET_TOL     <- 5                     # |aet - AET_REF|
WIND_MAX    <- 1.5                   # wind enters at -0.095, so match it too
PAL         <- c("antiquewhite3", "#CC0000", "#330000")   # as Fig S7

model      <- load_hubbell_model()
fitHubbell <- model$fit
dataset    <- model$dataset
NEED       <- setdiff(model$vars, "y")

g <- prepare_grid(YEAR, model)
g$x <- g$Longitude; g$y <- g$Latitude

df_win <- g %>% filter(between(x, xlims[1], xlims[2]), between(y, ylims[1], ylims[2]))
df_win$S_sigma <- exp(predict(fitHubbell, newdata = df_win)) / fitHubbell$sigma

# --- candidate pool: tropical, scoreable, inside the training envelope ---------
pool <- df_win %>%
  filter(zone == "Tropical", complete.cases(across(all_of(NEED))),
         !is.na(mess_zone), mess_zone >= -25,
         between(y, LAT_BAND[1], LAT_BAND[2]),
         abs(aet - AET_REF) <= AET_TOL, wind_speed <= WIND_MAX)
message("candidate cells (tropical, in envelope, matched on lat/AET/wind): ", nrow(pool))
message("  HFP range ", paste(round(range(pool$hfp),1), collapse = " - "),
        " | cells with HFP >= 25: ", sum(pool$hfp >= 25))

# --- pick a WELL-SEPARATED triple that is also monotone in predicted richness --
# Spreading the points out (for the funnels) lets them differ in latitude, and the
# latitude spline can then outweigh HFP and make the curves cross. So the search
# maximises the minimum pairwise separation SUBJECT TO predicted richness falling
# as HFP rises - the triple is real cells and real predictions, just not a
# cherry-picked pair that happens to invert.
pool$pred2000 <- predict(fitHubbell, newdata = pool, type = "response")

band <- function(lo, hi, k = 120) pool %>% filter(between(hfp, lo, hi)) %>%
  mutate(d = abs(aet - AET_REF)) %>% arrange(d) %>% slice(1:min(k, n()))
B_lo <- band(0, 1.5); B_md <- band(18, 28); B_hi <- band(38, 50)
message("band sizes: low ", nrow(B_lo), "  med ", nrow(B_md), "  high ", nrow(B_hi))

km <- function(a, b) as.numeric(distHaversine(c(a$x, a$y), c(b$x, b$y))) / 1000
best <- NULL
for (i in seq_len(nrow(B_hi))) for (j in seq_len(nrow(B_md))) {
  h <- B_hi[i, ]; m <- B_md[j, ]
  if (m$pred2000 <= h$pred2000) next
  if (abs(h$y - m$y) > LAT_SPREAD) next
  d_hm <- km(h, m); if (d_hm < MIN_SEP_KM) next
  ok <- B_lo %>% filter(pred2000 > m$pred2000,
                        abs(y - h$y) <= LAT_SPREAD, abs(y - m$y) <= LAT_SPREAD)
  if (!nrow(ok)) next
  d1 <- vapply(seq_len(nrow(ok)), function(r) km(ok[r, ], h), numeric(1))
  d2 <- vapply(seq_len(nrow(ok)), function(r) km(ok[r, ], m), numeric(1))
  sc <- pmin(d1, d2, d_hm); w <- which.max(sc)
  if (sc[w] >= MIN_SEP_KM && (is.null(best) || sc[w] > best$score))
    best <- list(score = sc[w], rows = bind_rows(ok[w, ], m, h))
}
if (is.null(best)) stop("no separated, monotone triple found - relax MIN_SEP_KM or LAT_SPREAD")
message("min pairwise separation: ", round(best$score), " km")
mklab <- function(h, a) sprintf("HFP %.0f  \u00b7  AET %.0f mm/yr", h, a * 10)
pts <- best$rows %>% arrange(hfp) %>%
  mutate(lab = mklab(hfp, aet), lab = factor(lab, levels = mklab(hfp, aet)))

message("\nchosen points (for the caption - everything but HFP should be close):")
print(as.data.frame(pts %>% transmute(lon = round(x,2), lat = round(y,2),
        HFP = round(hfp,1), AET_mm = round(aet*10), wind = round(wind_speed,2),
        wetlands = round(wetlands,2), tempres = round(resid_temperature_lat,2),
        realm = as.character(realm))), row.names = FALSE)
message("latitude spread: ", round(diff(range(pts$y)),2), " deg   AET spread: ",
        round(diff(range(pts$aet))*10), " mm/yr")

# --- accumulation curves ------------------------------------------------------
df_curves <- bind_rows(lapply(seq_len(nrow(pts)), function(k) {
  cv <- predict_curve(fit = fitHubbell, n = 2000, xnew = pts[k, ], .vcov = model$vcov)
  data.frame(lab  = pts$lab[k],
             n    = as.numeric(cv$n),
             pred = as.numeric(unlist(cv$mean)),
             se   = as.numeric(unlist(cv$se)))
}))
message("\npredicted richness at n = 2000:")
print(as.data.frame(df_curves %>% group_by(lab) %>% slice_max(n, n = 1) %>%
                    transmute(lab, n, pred = round(pred), lo = round(pred-1.96*se),
                              hi = round(pred+1.96*se))), row.names = FALSE)

# --- panels -------------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- world[world$name != "Antarctica", ]
pal_emrl <- c("gray25","gray50","gray80","#D3F2A3FF","#97E196FF","#6CC08BFF",
              "#02734AFF","#217A79FF","#105965FF","#03334AFF")

pMap <- ggplot() +
  theme_classic() +
  geom_tile(data = df_win, aes(x = x, y = y, fill = S_sigma)) +
  geom_sf(data = world, fill = NA, color = "gray30", size = 0.2) +
  scale_fill_gradientn(
    name = expression(S[sigma]), colours = pal_emrl, na.value = NA,
    breaks = c(20, 40, 60, 80),
    guide = guide_colourbar(theme = theme(
      legend.title = element_text(size = 17), legend.title.position = "top",
      legend.text = element_text(size = 11),
      legend.key.width = unit(0.42,"cm"), legend.key.height = unit(2.7,"cm"),
      legend.ticks = element_line(colour = "white"), legend.frame = element_blank()))) +
  coord_sf(expand = FALSE) +
  geom_point(data = dataset %>% distinct(Longitude, Latitude),
             aes(x = Longitude, y = Latitude), size = 1.2) +
  geom_point(data = pts, aes(x = x, y = y, fill = NULL, colour = lab),
             size = 4.2, shape = 21, fill = PAL, colour = "black", stroke = 1.1) +
  ylim(ylims) + xlim(xlims) +
  theme(axis.title = element_blank(), legend.position = "inside",
        legend.position.inside = c(0.105, 0.34),
        legend.justification = c(0.5, 0.5), legend.background = element_blank())

pCurves <- ggplot(df_curves, aes(x = n, y = pred, color = lab, fill = lab)) +
  geom_line(linewidth = 0.9) +
  geom_ribbon(aes(ymin = pred - 1.96*se, ymax = pred + 1.96*se), alpha = 0.20, color = NA) +
  theme_classic() +
  xlab("Number of individuals sampled (n)") + ylab("Predicted number of BINs") +
  scale_color_manual(name = NULL, values = setNames(PAL, levels(pts$lab))) +
  scale_fill_manual(name = NULL, values = setNames(PAL, levels(pts$lab))) +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.background = element_blank(),
        legend.text = element_text(size = 12))

out <- file.path(FIG, "Fig1c_accumulation_curves.png")
ggsave(plot = pMap + pCurves, filename = out, width = 12.54, height = 5.66, dpi = 300, bg = "white")
message("\nwrote ", out)

#===============================================================================
# Version B: ONE site, HFP set counterfactually. Everything else is held at that
# cell's own values, so the three curves differ only through human footprint.
#===============================================================================
base_cell <- pool %>% mutate(d = abs(hfp - 0)) %>% arrange(d, abs(aet - AET_REF)) %>% slice(1)
message("\ncounterfactual base cell: lon ", round(base_cell$x,2), " lat ", round(base_cell$y,2),
        "  AET ", round(base_cell$aet,1))
cf <- bind_rows(lapply(HFP_TARGETS, function(h) { z <- base_cell; z$hfp <- h; z }))
cf$lab <- factor(mklab(HFP_TARGETS, base_cell$aet),
                 levels = mklab(HFP_TARGETS, base_cell$aet))

df_cf <- bind_rows(lapply(seq_len(nrow(cf)), function(k) {
  cv <- predict_curve(fit = fitHubbell, n = 2000, xnew = cf[k, ], .vcov = model$vcov)
  data.frame(lab = cf$lab[k], n = as.numeric(cv$n),
             pred = as.numeric(unlist(cv$mean)), se = as.numeric(unlist(cv$se)))
}))
message("counterfactual richness at n = 2000:")
print(as.data.frame(df_cf %>% group_by(lab) %>% slice_max(n, n = 1) %>%
                    transmute(lab, pred = round(pred))), row.names = FALSE)

pMapB <- pMap
pMapB$layers[[length(pMapB$layers)]] <- NULL
pMapB <- pMapB + geom_point(data = base_cell, aes(x = x, y = y), size = 3.4,
                            shape = 21, fill = "white", colour = "black", stroke = 1.3)
# rebuild rather than %+%: the counterfactual levels ("HFP 30") differ from the
# matched-cell ones ("HFP 31"), so reusing pCurves' scale drops the third colour.
pCurvesB <- ggplot(df_cf, aes(x = n, y = pred, color = lab, fill = lab)) +
  geom_line(linewidth = 0.9) +
  geom_ribbon(aes(ymin = pred - 1.96*se, ymax = pred + 1.96*se), alpha = 0.20, color = NA) +
  theme_classic() +
  xlab("Number of individuals sampled (n)") + ylab("Predicted number of BINs") +
  scale_color_manual(name = NULL, values = setNames(PAL, levels(cf$lab))) +
  scale_fill_manual(name = NULL, values = setNames(PAL, levels(cf$lab))) +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.background = element_blank(),
        legend.text = element_text(size = 12))
ggsave(plot = pMapB + pCurvesB,
       filename = file.path(EXTRA, "Fig1c_accumulation_curves_counterfactual.png"),
       width = 12.54, height = 5.66, dpi = 300, bg = "white")
message("wrote Fig1c_accumulation_curves_counterfactual.png")
