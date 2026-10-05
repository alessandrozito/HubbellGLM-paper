# Figure S6: observed rarefaction vs predicted accumulation curves at two events.
# Poisson and NB (as in 09_cv_benchmark.R) can exceed y = n; Hubbell cannot.
#
# Needs:    01_fit_models.R
# Produces: FigS6_admissibility_curves.pdf; data/FigS6_abundances.rds (BIN counts)

source(file.path(path.expand("~/HubbellGLM-paper"), "preprocessing", "paths.R"))
suppressMessages({library(tidyverse); library(HubbellGLM); library(splines)
                  library(MASS); library(patchwork)})

EVENTS <- c("GMP#41359", "GMP#37397")   # Manu transect (PE), Khao Yai (TH)
STEP   <- 5                             # plot every 5th n

fourier_week <- function(week, k = 1) {
  X <- NULL
  for (j in 1:k) { X <- cbind(X, sin(2*j*pi*week/52), cos(2*j*pi*week/52))
    colnames(X)[(ncol(X)-1):ncol(X)] <- paste0(c("week_sin","week_cos"), j) }
  X
}

M   <- readRDS(file.path(MODELS, "models.rds")); K <- length(M$covariates)
hub <- M$poly[[K]]$fit
d   <- as.data.frame(hub$data)
str_cov <- M$covariates[K]

fp  <- glm(as.formula(paste("y - 1 ~ log(n) +", str_cov)), family = poisson("log"), data = d)
fnb <- MASS::glm.nb(as.formula(paste("y - 1 ~ log(n) +", str_cov)), data = d)

#---- BIN counts per event, from the raw BOLD export (filters as in preprocessing/02)
ab_f <- file.path(DATA, "FigS6_abundances.rds")
if (!file.exists(ab_f) || !setequal(names(readRDS(ab_f)), EVENTS)) {
  chr <- cols(.default = col_character())
  orphan <- read_tsv(file.path(QC, "qc_records_without_fieldid.tsv"), col_types = chr) %>%
    dplyr::select(site_code, coord, collection_date_start, collection_date_end, fieldid_new)
  far <- read_tsv(file.path(QC, "qc_records_dropped_far.tsv"), col_types = chr) %>%
    dplyr::select(fieldid, coord) %>% distinct()
  occ <- read_tsv(path.expand("~/gmtp_filtered_bcdm.tsv.zip"), col_types = chr,
                  col_select = c(fieldid, bin_uri, sampling_protocol, site_code, coord,
                                 collection_date_start, collection_date_end)) %>%
    filter(sampling_protocol %in% c("Malaise Trap", "Malaise trap"),
           !is.na(bin_uri), bin_uri != "None") %>%
    mutate(fieldid = na_if(fieldid, "None")) %>%
    left_join(orphan, by = c("site_code", "coord", "collection_date_start",
                             "collection_date_end")) %>%
    mutate(fieldid = coalesce(fieldid, fieldid_new)) %>%
    filter(fieldid %in% EVENTS) %>%
    anti_join(far, by = c("fieldid", "coord")) %>%
    count(fieldid, bin_uri)
  saveRDS(split(occ$n, occ$fieldid), ab_f)
}
ab <- readRDS(ab_f)

#---- Curves

curves_at <- function(ev, nn) {
  x  <- d[d$fieldid == ev, ]
  nd <- x[rep(1, length(nn)), ]; nd$n <- nn
  ci <- predict_curve(hub, xnew = x, n = x$n, npoints = min(x$n, 100), .vcov = M$poly[[K]]$vcov)
  bind_rows(
    tibble(n = nn, y = predict(fp,  nd, type = "response") + 1, curve = "Poisson", side = "mid"),
    tibble(n = nn, y = predict(fnb, nd, type = "response") + 1, curve = "Neg. Binomial", side = "mid"),
    tibble(n = nn, y = predict(hub, newdata = nd, type = "response"), curve = "Hubbell", side = "mid"),
    tibble(n = as.numeric(unlist(ci$n)), y = as.numeric(unlist(ci$mean)) + 1.96 * as.numeric(unlist(ci$se)), curve = "Hubbell 95% CI", side = "hi"),
    tibble(n = as.numeric(unlist(ci$n)), y = as.numeric(unlist(ci$mean)) - 1.96 * as.numeric(unlist(ci$se)), curve = "Hubbell 95% CI", side = "lo")
  ) %>% mutate(event = ev)
}
curves <- map_dfr(EVENTS, function(ev) curves_at(ev, seq(1, d$n[d$fieldid == ev], by = STEP)))
rar <- map_dfr(EVENTS, function(ev) {
  r <- BNPvegan::rarefaction(ab[[ev]])
  tibble(event = ev, n = seq_along(r), y = as.numeric(r))
})

LEV <- c("Poisson", "Neg. Binomial", "Hubbell", "Hubbell 95% CI")
COL <- c("Poisson" = "blue", "Neg. Binomial" = "forestgreen",
         "Hubbell" = "red", "Hubbell 95% CI" = "red")
LTY <- c("Poisson" = "solid", "Neg. Binomial" = "solid",
         "Hubbell" = "solid", "Hubbell 95% CI" = "dashed")

layers <- function(cc, rr) list(
  geom_abline(slope = 1, intercept = 0, linetype = "dotted"),
  geom_point(data = rr, aes(n, y, shape = "Rarefaction"), colour = "gray50", size = 1.6),
  geom_line(data = mutate(cc, curve = factor(curve, LEV)),
            aes(n, y, colour = curve, linetype = curve, group = interaction(curve, side))),
  scale_colour_manual(values = COL, name = NULL, breaks = LEV),
  scale_linetype_manual(values = LTY, name = NULL, breaks = LEV),
  scale_shape_manual(values = c(Rarefaction = 1), name = NULL),
  theme_classic(base_size = 13))

panel <- function(ev) {
  lim  <- d$n[d$fieldid == ev]
  ytop <- max(lim, filter(curves, event == ev)$y)          # show the full excess over y = n
  main <- ggplot() + layers(filter(curves, event == ev),
                            filter(rar, event == ev, n %% STEP == 1)) +
    coord_cartesian(xlim = c(0, lim), ylim = c(0, ytop)) +
    labs(title = ev, x = "Number of individuals sampled (n)", y = "Number of BINs") +
    theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02),
          legend.justification = c(1, 0), legend.background = element_blank(),
          legend.spacing.y = unit(0, "pt"), plot.title = element_text(face = "bold"))
  main
}

ggsave(file.path(FIG, "FigS6_admissibility_curves.pdf"),
       panel(EVENTS[1]) + panel(EVENTS[2]), width = 12, height = 5.4)

#---- Check: cached counts match the regression data
for (ev in EVENTS) stopifnot(length(ab[[ev]]) == d$y[d$fieldid == ev],
                             sum(ab[[ev]])    == d$n[d$fieldid == ev])
message("wrote FigS6_admissibility_curves.pdf")
