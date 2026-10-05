# Figures S5 and S9: paired % RMSE difference against Hubbell regression,
# per fold; and the four partitions on the map.
#
# Needs:    09_cv_benchmark.R
# Produces: FigS5_cv_benchmark.png, FigS9_cv_partitions.png, cv_deltas.tsv, cv_win_rates.tsv

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

suppressMessages({library(tidyverse); library(scales); library(sf)
                  library(rnaturalearth); library(patchwork)})

res_dir <- SIM
fig_dir <- FIG
REF     <- "HubbellPoly"

B <- read_tsv(file.path(res_dir, "cv_benchmark.tsv"), show_col_types = FALSE)

# Where the scheme drops `realm`, the "+realm" specification is the "+aet" one
# under another name - identical to the last decimal. Reporting both would show
# a null realm effect for a term that was never in the model.
dup <- B %>% dplyr::filter(spec %in% c("+aet", "+realm")) %>%
  dplyr::select(scheme, rep, fold, model, spec, rmse) %>%
  tidyr::pivot_wider(names_from = spec, values_from = rmse) %>%
  dplyr::group_by(scheme) %>%
  dplyr::summarise(same = isTRUE(all.equal(`+aet`, `+realm`)), .groups = "drop") %>%
  dplyr::filter(same) %>% dplyr::pull(scheme)
if (length(dup)) {
  message("dropping the duplicated +realm rows for: ", paste(dup, collapse = ", "))
  B <- B %>% dplyr::filter(!(scheme %in% dup & spec == "+realm"))
}
SPEC <- Sys.getenv("SPEC", tail(unique(B$spec), 1))
B <- B %>% dplyr::filter(spec == SPEC)
message("specification: ", SPEC, "   rows: ", nrow(B))

PRETTY <- c(Hubbell = "Hubbell (canonical, sigma = 0)",
            Poisson = "Poisson + log(n)", NB = "Neg. binomial + log(n)",
            Richness = "log richness + ns(n, 10)", Alpha = "Fisher's alpha",
            Poisson_off = "Poisson (offset)", NB_off = "Neg. binomial (offset)")
# Shorter forms for the axis: seven of these have to fit across four panels.
SHORT  <- c(Hubbell = "Hubbell (σ = 0)",
            Poisson = "Poisson + log(n)", NB = "NB + log(n)",
            Richness = "log y + ns(n, 10)", Alpha = "Fisher's α",
            Poisson_off = "Poisson (offset)", NB_off = "NB (offset)")
SCHEME_LAB <- c(random = "Random folds", block500 = "blockCV 500 km blocks",
                kmeans = "blockCV k-means clusters", realm = "Leave-one-realm-out")
# Order here drives everything downstream: the dodged pair, the legend and the
# summary tables all take their order from these names.
SET_LAB <- c(`in` = "In sample", out = "Held out")

#---- Held-out and in-sample, stacked ------------------------------------------
# One row per fold, estimator AND evaluation set. 09_cv_benchmark.R scores each
# fit twice and stores the two side by side (rmse / rmse_in); this turns that
# into the long form the pairing and the facets both want. Guarded, so the
# script still runs against a benchmark table written before *_in existed.

HAS_IN <- all(c("rmse_in", "sse_in", "n_in") %in% names(B)) && any(is.finite(B$rmse_in))
if (!HAS_IN)
  message("note: no in-sample columns in cv_benchmark.tsv - held-out row only")
L <- dplyr::bind_rows(
  B %>% dplyr::transmute(scheme, rep, fold, spec, model, err,
                         set = "out", n, sse, rmse),
  if (HAS_IN) B %>% dplyr::transmute(scheme, rep, fold, spec, model, err,
                                     set = "in", n = n_in, sse = sse_in,
                                     rmse = rmse_in))

#---- Pairing ------------------------------------------------------------------
# rep MUST be in the key: without it a competitor's fold 3 pairs against the
# reference's fold 3 from every repeat, which multiplies the rows and compares
# models fitted on different training sets. `set` too, or an in-sample RMSE
# would be scored against a held-out reference.
KEY <- c("scheme", "rep", "fold", "spec", "set")

ref <- L %>% dplyr::filter(model == REF) %>%
  dplyr::select(all_of(KEY), n, ref_rmse = rmse, ref_sse = sse)
# a fold where either side failed to fit is not a paired comparison
n_drop <- sum(is.na(L$rmse))
if (n_drop) message("dropping ", n_drop, " fold-model rows that failed to fit")
delta <- L %>% dplyr::filter(model != REF, !is.na(rmse)) %>%
  dplyr::inner_join(ref, by = KEY, suffix = c("", "_ref")) %>%
  dplyr::mutate(d_rmse = rmse - ref_rmse,
                pct    = 100 * (rmse - ref_rmse) / ref_rmse,
                model  = factor(model, levels = names(PRETTY)),
                scheme = factor(scheme, levels = names(SCHEME_LAB)),
                set    = factor(set, levels = intersect(names(SET_LAB), set)))
write_tsv(delta, file.path(res_dir, "cv_deltas.tsv"))

#---- Win rates and effect sizes, within repeat then across --------------------

per_rep <- delta %>% group_by(set, scheme, model, rep) %>%
  summarise(folds = dplyr::n(),
            wins = sum(d_rmse > 0, na.rm = TRUE),      # reference wins the fold
            med_pct = median(pct, na.rm = TRUE),
            wins_na = sum(is.na(d_rmse)),
            # n-weighted: one RMSE over all rows in the repeat
            pooled_pct = 100 * (sqrt(sum(sse, na.rm = TRUE) / sum(n, na.rm = TRUE)) /
                                sqrt(sum(ref_sse, na.rm = TRUE) / sum(n, na.rm = TRUE)) - 1),
            .groups = "drop")
win <- per_rep %>% group_by(set, scheme, model) %>%
  summarise(reps = dplyr::n(), folds = sum(folds),
            win_rate = round(100 * sum(wins) / sum(folds), 1),
            median_pct = round(median(med_pct, na.rm = TRUE), 2),
            iqr_lo = round(quantile(med_pct, .25, na.rm = TRUE), 2),
            iqr_hi = round(quantile(med_pct, .75, na.rm = TRUE), 2),
            pooled_pct = round(median(pooled_pct, na.rm = TRUE), 2), .groups = "drop") %>%
  dplyr::arrange(set, scheme, median_pct)
write_tsv(win, file.path(res_dir, "cv_win_rates.tsv"))
for (s in levels(delta$set)) {
  message("\n=== ", SET_LAB[[s]], ": ", REF, " vs each estimator (positive = ",
          REF, " is better) ===")
  print(as.data.frame(win %>% dplyr::filter(set == s) %>% dplyr::select(-set) %>%
                        dplyr::mutate(model = PRETTY[as.character(model)])),
        row.names = FALSE)
}

#---- The figure ---------------------------------------------------------------
# Estimator order on the x axis. Roughly cheapest to most expensive, but fixed
# rather than sorted on the fold medians: the two effort-as-covariate GLMs
# (Poisson, NB) belong next to each other, and sorting split them around the
# spline. The three near-ties come first, then the sigma = 0 contrast, then the
# estimators that are not competitive.
ORDER <- c("Poisson", "NB", "Richness", "Hubbell", "Alpha", "Poisson_off", "NB_off")
stopifnot(setdiff(as.character(delta$model), ORDER) == character(0))
delta <- delta %>% dplyr::mutate(model = factor(as.character(model), levels = ORDER))

# Fold counts are not in the strip: they belong in the caption and the
# win-rate table, and repeating them four times crowds the panel titles.

# One hue per evaluation set, light for the box and saturated for the folds
# drawn on top of it. Colour carries meaning here, so this legend stays.
FILL <- c(out = "#E5B79F", `in`  = "#9DC3C2")
LINE <- c(out = "#A8401E", `in`  = "#215457")
DODGE <- 0.78

p <- ggplot(delta, aes(model, pct, fill = set, colour = set)) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.4) +
  geom_boxplot(position = position_dodge(width = DODGE), width = 0.66,
               outlier.shape = NA, linewidth = 0.32) +
  # every fold drawn. jitter.height = 0 so a point keeps its exact difference
  # and only moves sideways, inside its own dodged box.
  geom_point(position = position_jitterdodge(jitter.width = 0.3,
                                             jitter.height = 0,
                                             dodge.width = DODGE, seed = 1),
             size = 0.9, alpha = 0.5, stroke = 0) +
  facet_wrap(~ scheme, nrow = 1, labeller = labeller(scheme = SCHEME_LAB)) +
  scale_fill_manual(values = FILL, labels = SET_LAB, name = NULL) +
  scale_colour_manual(values = LINE, labels = SET_LAB, name = NULL) +
  scale_x_discrete(labels = SHORT) +
  scale_y_continuous(transform = scales::pseudo_log_trans(sigma = 2),
                     breaks = c(-20, -5, 0, 5, 20, 50, 100, 200),
                     minor_breaks = NULL,
                     labels = function(x) paste0(x, "%")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, linewidth = 0.4)),
         colour = "none") +
  labs(x = element_blank(),
       y = sprintf("Relative difference in RMSE\nwith polynomial Hubbell")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey96", colour = "grey70"),
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "right",
        #legend.margin = margin(b = -4),
        axis.text.x = element_text(size = 8.5, angle = 45, hjust = 1,
                                   colour = "black"),
        axis.text.y = element_text(size = 9),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(fig_dir, "FigS5_cv_benchmark.png"), p,
       width = 11.51, height = 3.73)

#---- The partitions themselves ------------------------------------------------
# The boxplots are only interpretable next to the geography they average over: a
# 500 km block and a k-means cluster are both "spatial", and they separate the
# sample very differently. First repeat of each scheme, one colour per fold.

FD <- readRDS(file.path(res_dir, "cv_folds.rds"))
M  <- readRDS(file.path(MODELS, "models.rds"))
d  <- readRDS(file.path(DATA, "GMTP_analysis_dataset.rds")) %>%
  dplyr::filter(fieldid %in% M$ids) %>% dplyr::distinct(fieldid, .keep_all = TRUE)
d  <- d[match(M$ids, d$fieldid), ]
stopifnot(identical(d$fieldid, FD$ids))

sites <- d %>% group_by(site_code) %>%
  summarise(Longitude = mean(Longitude), Latitude = mean(Latitude),
            events = dplyr::n(), .groups = "drop")
site_of_event <- match(d$site_code, sites$site_code)

world <- tryCatch(ne_countries(scale = "medium", returnclass = "sf") %>%
                    dplyr::filter(name != "Antarctica"), error = function(e) NULL)
MAP_TITLE <- c(random = "Random, 10 folds over events",
               block500 = sprintf("blockCV spatial, %g km blocks",
                                  FD$block_m / 1000),
               kmeans = "blockCV k-means clusters",
               realm = "Leave-one-realm-out")

map_panel <- function(nm) {
  g <- FD$folds %>% dplyr::filter(scheme == nm, rep == 1)
  # one fold per SITE: under every scheme but `random` a site's events share a
  # fold, and plotting events would overplot 15,711 points onto 1,196 places.
  pd <- sites %>% dplyr::mutate(fold = factor(tapply(g$fold, site_of_event,
                                                     function(z) z[1])))
  ggplot() +
    {if (!is.null(world)) geom_sf(data = world, fill = "grey96",
                                  colour = "grey75", linewidth = 0.1)} +
    geom_point(data = pd, aes(Longitude, Latitude, colour = fold),
               size = 0.7, alpha = 0.9) +
    scale_colour_viridis_d(option = "turbo", guide = "none") +
    coord_sf(expand = FALSE, ylim = c(-56, 84)) +
    labs(title = MAP_TITLE[[nm]]) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0,
                                    margin = margin(b = 3)),
          plot.background = element_rect(fill = "white", colour = NA))
}
schemes_present <- intersect(names(MAP_TITLE), unique(FD$folds$scheme))
ggsave(file.path(fig_dir, "FigS9_cv_partitions.png"),
       wrap_plots(lapply(schemes_present, map_panel), ncol = 2),
       width = 13, height = 5.6, dpi = 200)

message("\nwrote FigS5_cv_benchmark.png, FigS9_cv_partitions.png, ",
        "cv_deltas.tsv, cv_win_rates.tsv")
