# Figure S1: 95% CIs for AET and zone-specific HFP effects under four
# standard-error estimators (Normal, Quasi, White, Jaccard), both links.
#
# Needs:    01_fit_models.R
# Produces: FigS1_ci_comparison.png/.pdf, ci_comparison.tsv

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

suppressMessages({library(tidyverse); library(HubbellGLM); library(splines)
                  library(sandwich)})

res_dir <- MODELS
fig_dir <- FIG

M <- readRDS(file.path(res_dir, "models.rds"))
K <- length(M$covariates)                      # the final specification, M7

# The zone-specific HFP slope is hfp + zoneX:hfp, so each effect is a contrast
# vector over the coefficient names. AET is a single coefficient.
EFFECTS <- list(
  "AET"                = "aet",
  "HFP (Continental)"  = "hfp",
  "HFP (Dry)"          = c("hfp", "zoneDry:hfp"),
  "HFP (Polar)"        = c("hfp", "zonePolar:hfp"),
  "HFP (Temperate)"    = c("hfp", "zoneTemperate:hfp"),
  "HFP (Tropical)"     = c("hfp", "zoneTropical:hfp"))

# One row per effect: the estimate and the standard error implied by V.
contrasts_of <- function(b, V) {
  purrr::imap_dfr(EFFECTS, function(k, lab) {
    stopifnot(all(k %in% names(b)))
    tibble(term = lab, estimate = sum(b[k]), std_error = sqrt(sum(V[k, k, drop = FALSE])))
  })
}

# fourier_week lives in the formula of the stored fits, so the quasi refit needs
# it in scope here as well.
fourier_week <- function(week, k = 1) {
  X <- NULL
  for (j in 1:k) { X <- cbind(X, sin(2*j*pi*week/52), cos(2*j*pi*week/52))
    colnames(X)[(ncol(X)-1):ncol(X)] <- paste0(c("week_sin","week_cos"), j) }
  X
}

LINKS <- list(Polynomial = M$poly[[K]], Canonical = M$log[[K]])

ci <- purrr::imap_dfr(LINKS, function(m, link) {
  b <- coef(m$fit)
  f <- stats::as.formula(paste0("cbind(n, y) ~ ", M$covariates[K]), env = environment())
  fq <- HubbellGLM(formula = f, data = m$fit$data,
                   family = quasihubbell(sigma = m$sigma))
  V <- list(`1. Normal`  = vcov(m$fit),
            `2. Quasi`   = vcov(fq),
            `3. White`   = sandwich::vcovHC(m$fit, type = "HC0"),
            `4. Jaccard` = m$vcov)
  purrr::imap_dfr(V, ~ contrasts_of(b, .x) %>% mutate(estimator = .y)) %>%
    mutate(link = sprintf("%s link (sigma = %.3f)", link, m$sigma))
}) %>%
  mutate(lo = estimate - 1.96 * std_error,
         hi = estimate + 1.96 * std_error,
         width = hi - lo,
         term = factor(term, levels = rev(names(EFFECTS))))
write_tsv(ci, file.path(res_dir, "ci_comparison.tsv"))

# Ratio taken WITHIN link and term, then summarised: the Jaccard width has to be
# the one for the same effect, not the median over all of them.
message("\n=== interval width relative to the Jaccard interval ===")
print(as.data.frame(ci %>% group_by(link, term) %>%
  mutate(ratio = width / width[estimator == "4. Jaccard"]) %>%
  group_by(link, estimator) %>%
  summarise(median_ratio = round(median(ratio), 2),
            min_ratio = round(min(ratio), 2),
            max_ratio = round(max(ratio), 2), .groups = "drop")),
  row.names = FALSE)

p <- ggplot(ci, aes(estimate, term, colour = estimator, shape = estimator)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35", linewidth = 0.4) +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.75),
                  linewidth = 0.5, size = 0.35) +
  facet_wrap(~ link, scales = "free_x") +
  scale_colour_manual(values = c("1. Normal"  = "#8C8C8C", "2. Quasi" = "#2E6F73",
                                 "3. White"   = "#C7522B", "4. Jaccard" = "#111111"),
                      name = "SE estimator") +
  scale_shape_manual(values = c(0, 1, 2, 16), name = "SE estimator") +
  labs(x = "Marginal effect on log-diversity (95% CI)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "right",
        strip.background = element_rect(fill = "grey96", colour = "grey70"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(fig_dir, "FigS1_ci_comparison.png"), p, width = 11, height = 3.9, dpi = 300)
ggsave(file.path(fig_dir, "FigS1_ci_comparison.pdf"), p, width = 11, height = 3.9)
message("\nwrote FigS1_ci_comparison.png/.pdf, ci_comparison.tsv")
