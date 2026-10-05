# Fit the nested models M0-M7 (sigma per model), with Jaccard standard errors.
#
# Needs:    preprocessing/ (GMTP_analysis_dataset.rds, JaccardSim_sparse.rds)
# Produces: models.rds; Table S2 (TableS2_deviance.tex), Table S4 (TableS4_models_log.tex)

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

library(tidyverse)
library(Matrix)
library(HubbellGLM)
library(splines)
library(lmtest)
library(sandwich)
library(multcomp)

# Set to FALSE to reuse the models already in models.rds instead of refitting
# M0-M7. RERUN=0 in the environment does the same thing.
rerun <- Sys.getenv("RERUN", "0") != "0"   # RERUN=1 refits; default reuses models.rds

data_dir <- DATA
res_dir <- MODELS
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
message("zone column: ", ZONE, "   results -> ", normalizePath(res_dir, mustWork = FALSE))

#---- Fourier term, as in analysis/1_run_models.R -----------------------------

fourier_week <- function(week, k = 1) {
  X <- NULL
  for (j in 1:k) {
    X <- cbind(X, sin(2 * j * pi * week / 52), cos(2 * j * pi * week / 52))
    colnames(X)[(ncol(X) - 1):ncol(X)] <- paste0(c("week_sin", "week_cos"), j)
  }
  X
}

#---- Data --------------------------------------------------------------------

d <- readRDS(file.path(data_dir, "GMTP_analysis_dataset.rds"))
J <- readRDS(file.path(data_dir, "JaccardSim_sparse.rds"))
message("dataset ", nrow(d), " rows; Jaccard ", nrow(J), " x ", ncol(J),
        " (", round(100 * length(J@x) / prod(dim(J)), 1), "% filled, ",
        round(as.numeric(object.size(J)) / 1e6), " MB)")

#   Set YEAR_MAX=2026 to include them as a sensitivity check.
YEAR_MAX <- as.integer(Sys.getenv("YEAR_MAX", "2024"))
n_pre <- nrow(d)
d <- d %>% filter(year <= YEAR_MAX)
if (nrow(d) < n_pre) {
  message("YEAR_MAX=", YEAR_MAX, ": dropped ", n_pre - nrow(d),
          " rows after ", YEAR_MAX, " (", nrow(d), " remain)")
}

OBSERVED_ONLY <- Sys.getenv("OBSERVED_ONLY", "0") != "0"
if (OBSERVED_ONLY) {
  n0 <- nrow(d)
  d <- d %>% filter(hfp_src == "observed", aet_src == "observed")
  message("OBSERVED_ONLY: kept ", nrow(d), " of ", n0, " rows (dropped ",
          n0 - nrow(d), " with modelled aet/hfp)")
}

d <- d %>%
  filter(fieldid %in% rownames(J)) %>%
  distinct(fieldid, .keep_all = TRUE) %>%
  mutate(zone  = factor(.data[[ZONE]], levels = c("Continental", "Dry", "Polar",
                                                  "Temperate", "Tropical")),
         realm = factor(realm))
message("  modelling rows: ", nrow(d))

ids <- d$fieldid
Js  <- J[ids, ids]
stopifnot(identical(rownames(Js), ids))
message("  similarity subset ", nrow(Js), " x ", ncol(Js), ", ",
        format(length(Js@x), big.mark = ","), " stored")

#---- Model sequence ----------------------------------------------------------

has <- function(v) v %in% names(d) && sum(!is.na(d[[v]])) > 0.5 * nrow(d)
CONTROLS <- Sys.getenv("CONTROLS", "1") != "0"

covariates <- c(
  "1",
  "1 + aet",
  "1 + aet + realm",
  "1 + aet + realm + fourier_week(week_adj, k = 1)",
  "1 + aet + realm + fourier_week(week_adj, k = 1) + ns(Latitude, 6)",
  "1 + aet + realm + fourier_week(week_adj, k = 1) + ns(Latitude, 6) + zone * hfp"
)
col_names <- c("Base", "+aet", "+realm", "+week", "+nsLat", "+zone*hfp")

# Terms may be expressions - I(year - 2015), log10(landmass_km2) - so resolve
# them to the underlying column names before any availability/completeness check.
vars_in <- function(x) if (!length(x)) character(0) else {
  unique(all.vars(stats::as.formula(paste("~", paste(x, collapse = " + ")))))
}

ctrl <- c("wind_speed", "wetlands", "collection_days", "resid_temperature_lat",
          "I(year - 2015)", "log10(landmass_km2)", "log10(dist_mainland_km + 1)")
ctrl <- ctrl[vapply(ctrl, function(t) all(vapply(vars_in(t), has, logical(1))), logical(1))]
if (CONTROLS && length(ctrl)) {
  covariates <- c(covariates, paste(covariates[6], "+", paste(ctrl, collapse = " + ")))
  col_names  <- c(col_names, "+controls")
  message("  controls model includes: ", paste(ctrl, collapse = ", "))
} else {
  message("  no controls model (available: ",
          if (length(ctrl)) paste(ctrl, collapse = ", ") else "none", ")")
}

# M7 adds the weather anomalies: departures of the COLLECTION WINDOW's weather
# from that cell's own 1991-2020 day-of-year normal (13_compute_anomalies.R).
anml <- c("temperature_2m_anml", "total_precipitation_sum_anml",
          "relative_humidity_anml", "wind_speed_anml")
anml <- anml[vapply(anml, has, logical(1))]
if (Sys.getenv("ANOMALIES", "1") != "0" && length(anml) == 4) {
  covariates <- c(covariates, paste(covariates[length(covariates)], "+",
                                    paste(anml, collapse = " + ")))
  col_names  <- c(col_names, "+anomalies")
  message("  anomaly model M", length(covariates) - 1, " includes: ",
          paste(anml, collapse = ", "))
} else if (length(anml) && length(anml) < 4) {
  message("  only ", length(anml), " of 4 anomaly columns usable - M7 skipped")
}

keep <- complete.cases(d[, unique(c("n", "y", vars_in(covariates[length(covariates)])))])
if (!all(keep)) {
  message("  dropping ", sum(!keep), " of ", nrow(d),
          " rows incomplete on the full covariate set")
  d <- d[keep, ]; ids <- d$fieldid; Js <- J[ids, ids]
}

#---- sigma -------------------------------------------------------------------
# SIGMA_PER_MODEL=1 (the default) estimates sigma SEPARATELY for every
# specification M0--M7
#
# SIGMA_PER_MODEL=0 restores a single sigma; SIGMA_MODEL then picks whether it
# comes from the full specification (default) or from "~ 1"
SIGMA_PER_MODEL <- Sys.getenv("SIGMA_PER_MODEL", "1") != "0"
SIGMA_MODEL     <- Sys.getenv("SIGMA_MODEL", "full")
sigma_formula   <- if (identical(SIGMA_MODEL, "intercept")) "cbind(n, y) ~ 1" else
  paste0("cbind(n, y) ~ ", covariates[length(covariates)])

if (nzchar(Sys.getenv("SIGMA"))) {
  best_sigma <- as.numeric(Sys.getenv("SIGMA"))
  SIGMA_PER_MODEL <- FALSE
  message("sigma (given): ", best_sigma)
} else if (SIGMA_PER_MODEL) {
  message("sigma will be estimated separately for each specification")
  best_sigma <- NA_real_          # filled per model below
} else {
  message("estimating one sigma on the ", SIGMA_MODEL, " model:")
  message("  ", sigma_formula)
  t0 <- Sys.time()
  best_sigma <- estimate_sigma(formula = sigma_formula, data = d, verbose = FALSE)
  message("  sigma = ", signif(best_sigma, 7), "  (",
          round(difftime(Sys.time(), t0, units = "secs")), " s)")
}

#---- Fit ---------------------------------------------------------------------

# `sigma` is either a single value used for every model, or NA meaning
# "estimate it for this specification". The canonical log link always passes 0.
estimate_models <- function(covariates, dataset, sigma, similarity, per_model = FALSE) {
  lapply(seq_along(covariates), function(i) {
    f <- paste0("cbind(n, y) ~ ", covariates[i])
    s_i <- if (per_model) {
      t0 <- Sys.time()
      v <- estimate_sigma(formula = f, data = dataset, verbose = TRUE)
      message("  M", i - 1, ": sigma = ", signif(v, 7), "  (",
              round(difftime(Sys.time(), t0, units = "secs")), " s)")
      v
    } else {
      message("  M", i - 1, ": ", covariates[i]); sigma
    }
    fit <- HubbellGLM(formula = f, data = dataset, family = hubbell(sigma = s_i))
    list(fit = fit, vcov = vcov_shared(fit, similarity), sigma = s_i)
  })
}

# Reuse is refused if the stored models were fitted on a different
# specification or a different set of events - otherwise every table below would
# describe a model this script no longer defines.
mods_f <- file.path(res_dir, "models.rds")
reuse  <- !rerun && file.exists(mods_f)
if (reuse) {
  M <- readRDS(mods_f)
  ok_spec <- identical(M$covariates, covariates)
  ok_rows <- identical(sort(M$ids), sort(ids))
  reuse <- ok_spec && ok_rows
  if (!reuse)
    message("\ncannot reuse models.rds - refitting:",
            if (!ok_spec) "\n  the specification has changed" else "",
            if (!ok_rows) paste0("\n  fitted on ", length(M$ids),
                                 " events, this run has ", length(ids)) else "")
} else if (!rerun) {
  message("\nrerun = FALSE but ", basename(mods_f), " does not exist - fitting")
}

if (reuse) {
  message("\nreusing models.rds: ", length(M$poly), " models, sigma = ",
          signif(M$sigma, 6), "   [rerun = FALSE]")
  Models_poly    <- M$poly
  Models_log     <- M$log
  sigma_by_model <- M$sigma_by_model
  best_sigma     <- M$sigma
  ids <- M$ids
  d   <- d[match(ids, d$fieldid), ]
  Js  <- J[ids, ids]
  stopifnot(identical(d$fieldid, ids))
} else {
  message("\nPolynomial link (sigma ",
          if (SIGMA_PER_MODEL) "per specification" else paste("=", signif(best_sigma, 6)), ")")
  t0 <- Sys.time()
  Models_poly <- estimate_models(covariates, d, best_sigma, Js, per_model = SIGMA_PER_MODEL)
  message("  ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

  # The final specification's sigma is the one script 11 predicts with.
  sigma_by_model <- vapply(Models_poly, function(m) m$sigma, numeric(1))
  best_sigma <- sigma_by_model[length(sigma_by_model)]

  message("\nCanonical log link (sigma = 0)")
  t0 <- Sys.time()
  Models_log <- estimate_models(covariates, d, 0, Js)
  message("  ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}

#---- Tidy coefficient tables -------------------------------------------------

tidy_models <- function(Models, label) {
  map_dfr(seq_along(Models), function(i) {
    m  <- Models[[i]]
    se <- sqrt(diag(m$vcov))
    b  <- coef(m$fit)
    tibble(link = label, model = col_names[i], term = names(b),
           estimate = as.numeric(b), std_error = as.numeric(se),
           z = as.numeric(b) / as.numeric(se),
           p_value = 2 * pnorm(-abs(as.numeric(b) / as.numeric(se))),
           sigma = m$fit$sigma,
           deviance = m$fit$deviance, df_residual = m$fit$df.residual,
           aic = tryCatch(AIC(m$fit), error = function(e) NA_real_))
  })
}
coefs <- bind_rows(tidy_models(Models_poly, "polyseries"),
                   tidy_models(Models_log,  "canonical_log"))
write_tsv(coefs, file.path(res_dir, "model_coefficients.tsv"))

#---- ... and the same thing as the manuscript's stepwise LaTeX tables --------
# tab:models_poly and tab:models_log. The published versions were stargazer
# output edited by hand, so the labels and the grouping could drift from the
# models they described; here they are generated from the fits themselves.
# Standard errors are the Jaccard ones, as everywhere else in this script.

TERM_GROUPS <- list(
  "Water--energy"                    = c(aet = "AET"),
  "Biogeographic realm"              = c(realmAustralasia = "Australasia",
                                         realmIndomalayan = "Indomalaya",
                                         realmNearctic    = "Nearctic",
                                         realmNeotropic   = "Neotropic",
                                         realmPalearctic  = "Palearctic"),
  "Seasonality"                      = setNames(
    c("Week (sine)", "Week (cosine)"),
    paste0("fourier_week(week_adj, k = 1)week_", c("sin1", "cos1"))),
  "Latitude splines"                 = setNames(paste("Spline", 1:6),
                                                paste0("ns(Latitude, 6)", 1:6)),
  "K\\\"oppen zone"                    = c(zoneDry = "Dry", zonePolar = "Polar",
                                         zoneTemperate = "Temperate",
                                         zoneTropical = "Tropical"),
  "Human footprint"                  = c(hfp = "HFP",
                                         `zoneDry:hfp`       = "Dry $\\times$ HFP",
                                         `zonePolar:hfp`     = "Polar $\\times$ HFP",
                                         `zoneTemperate:hfp` = "Temperate $\\times$ HFP",
                                         `zoneTropical:hfp`  = "Tropical $\\times$ HFP"),
  "Landscape and sampling controls"  = c(wind_speed = "Wind speed",
                                         wetlands = "Wetlands",
                                         collection_days = "Collection days",
                                         resid_temperature_lat = "Temp residuals",
                                         `I(year - 2015)` = "Year ($-$2015)",
                                         `log10(landmass_km2)` = "Landmass area (log$_{10}$)",
                                         `log10(dist_mainland_km + 1)` = "Isolation (log$_{10}$)"),
  "Weather anomalies"                = c(temperature_2m_anml = "Temp anomaly",
                                         total_precipitation_sum_anml = "Precip anomaly",
                                         relative_humidity_anml = "Rel hum anomaly",
                                         wind_speed_anml = "Wind anomaly"))

COL_TEX <- c(Base = "Base", `+aet` = "{+}\\,AET", `+realm` = "{+}\\,Realm",
             `+week` = "{+}\\,Week", `+nsLat` = "{+}\\,Lat.",
             `+zone*hfp` = "{+}\\,Zone$\\times$HFP", `+controls` = "{+}\\,Controls",
             `+anomalies` = "{+}\\,Weather")

stepwise_tex <- function(Models, label, title, file) {
  est <- lapply(Models, function(m) {
    b <- coef(m$fit); se <- sqrt(diag(m$vcov))
    list(b = b, se = se, p = 2 * pnorm(-abs(b / se)))
  })
  K <- length(Models)
  cell <- function(k, term) {
    e <- est[[k]]
    if (!term %in% names(e$b)) return("")
    st <- if (e$p[term] < 0.01) "^{***}" else if (e$p[term] < 0.05) "^{**}" else
          if (e$p[term] < 0.10) "^{*}" else ""
    sprintf("$%.3f%s$\\,(%.3f)", e$b[term], st, e$se[term])
  }
  row_for <- function(term, lab) {
    v <- vapply(seq_len(K), cell, character(1), term = term)
    if (all(v == "")) return(NULL)                 # term in no model - skip it
    paste0("\t\t", lab, " & ", paste(v, collapse = " & "), "\\\\")
  }
  hdr <- COL_TEX[names(COL_TEX) %in% col_names]
  L <- c("\\begin{table}", "\t\\centering",
         paste0("\t\\caption{\\textbf{", title, "} Coefficients are added in nested ",
                "blocks from the base intercept-only model (M0) to the full ",
                "weather-anomaly model (M", K - 1, "). Realm and climate-zone ",
                "coefficients are relative to the omitted reference category.}"),
         paste0("\t\\label{", label, "}"),
         "\t\\resizebox{\\textwidth}{!}{\\linespread{1}\\selectfont",   # single-space: fits one page
         paste0("\t\\begin{tabular}{l ", strrep("c", K), "}"),
         "\t\t\\\\[-0.5em]", "\t\t\\hline",
         paste0("\t\t& ", paste(paste0("M", seq_len(K) - 1), collapse = " & "), "\\\\"),
         paste0("\t\t& ", paste(hdr, collapse = " & "), "\\\\[0.5em]"),
         "\t\t\\hline\\\\[-0.5em]")
  for (g in names(TERM_GROUPS)) {
    rows <- compact(imap(TERM_GROUPS[[g]], function(lab, term) row_for(term, lab)))
    if (!length(rows)) next
    rows[[length(rows)]] <- paste0(rows[[length(rows)]], "[0.5em]")
    L <- c(L, paste0("\t\t\\multicolumn{", K + 1, "}{l}{\\textit{", g, "}}\\\\"),
           unlist(rows, use.names = FALSE))
  }
  # Anything the model estimated that no group claims, so nothing is dropped
  # silently when the specification changes.
  known <- c("(Intercept)", unlist(lapply(TERM_GROUPS, names), use.names = FALSE))
  extra <- setdiff(unique(unlist(lapply(est, function(e) names(e$b)))), known)
  if (length(extra)) {
    message("  NOTE: ", length(extra), " term(s) not in TERM_GROUPS, listed as Other: ",
            paste(extra, collapse = ", "))
    rows <- compact(lapply(extra, function(t) row_for(t, gsub("_", "\\\\_", t))))
    rows[[length(rows)]] <- paste0(rows[[length(rows)]], "[0.5em]")
    L <- c(L, paste0("\t\t\\multicolumn{", K + 1, "}{l}{\\textit{Other}}\\\\"),
           unlist(rows, use.names = FALSE))
  }
  L <- c(L, "\t\t\\hline\\\\[-0.5em]",
         row_for("(Intercept)", "Constant"),
         paste0("\t\tDeviance & ",
                paste(formatC(vapply(Models, function(m) m$fit$deviance, numeric(1)),
                              format = "d", big.mark = ","), collapse = " & "), "\\\\"),
         paste0("\t\t$\\sigma$ & ",
                paste(sprintf("%.4f", vapply(Models, function(m) m$sigma, numeric(1))),
                      collapse = " & "), "\\\\[0.5em]"),
         "\t\t\\hline", "\t\\end{tabular}", "\t}", "\t\\vspace{1ex}", "\t\\raggedright",
         paste0("\t\\footnotesize Response: \\texttt{cbind(n, y)}. Standard errors ",
                "(Jaccard-adjusted for shared species) in parentheses. ",
                "$^{*}p<0.1$;\\, $^{**}p<0.05$;\\, $^{***}p<0.01$."),
         "\\end{table}")
  writeLines(L, file)
  message("  wrote ", file)
}

stepwise_tex(Models_poly, "tab:models_poly",
             "Stepwise Hubbell regression estimates (polynomial link).",
             file.path(res_dir, "models_poly_sigma_per_model.tex"))   # per-model sigma; not in the paper
stepwise_tex(Models_log, "tab:models_log",
             "Stepwise Hubbell regression estimates (canonical link).",
             file.path(res_dir, "TableS4_models_log.tex"))

#---- Total HFP effect within each zone ---------------------------------------
# hfp alone is the effect in the reference zone (Continental); every other zone
# needs its interaction added, which is what these contrasts do.

zone_lv <- levels(d$zone)
hyp <- c("hfp = 0", paste0("hfp + zone", zone_lv[-1], ":hfp = 0"))
final_poly <- Models_poly[[length(Models_poly)]]
gl <- summary(glht(final_poly$fit, linfct = hyp, vcov = final_poly$vcov))
hfp_zone <- tibble(
  zone     = c(zone_lv[1], zone_lv[-1]),
  estimate = as.numeric(gl$test$coefficients),
  std_error= as.numeric(gl$test$sigma),
  z        = as.numeric(gl$test$tstat),
  p_value  = as.numeric(gl$test$pvalues)) %>%
  mutate(pct_per_5_units = 100 * (exp(5 * estimate) - 1))
write_tsv(hfp_zone, file.path(res_dir, "hfp_effect_by_zone.tsv"))

#---- Variance estimator comparison on the final model ------------------------

fq <- HubbellGLM(formula = paste0("cbind(n, y) ~ ", covariates[length(covariates)]),
                 data = d, family = quasihubbell(sigma = best_sigma))
V <- list(`1. Normal` = vcov(final_poly$fit),
          `2. Quasi`  = vcov(fq),
          `3. White`  = vcovHC(final_poly$fit, type = "HC0"),
          `4. Jaccard`= final_poly$vcov)
se_cmp <- map_dfr(names(V), ~ tibble(estimator = .x, term = names(coef(final_poly$fit)),
                                     std_error = sqrt(diag(V[[.x]]))))
write_tsv(se_cmp, file.path(res_dir, "se_comparison.tsv"))

if (!reuse) saveRDS(list(sigma = best_sigma, zone_column = ZONE, year_max = YEAR_MAX,
             observed_only = OBSERVED_ONLY,
             sigma_model = SIGMA_MODEL,
             sigma_per_model = SIGMA_PER_MODEL, sigma_by_model = sigma_by_model,
             sigma_formula = sigma_formula,
             covariates = covariates, col_names = col_names,
             poly = Models_poly, log = Models_log, ids = ids), mods_f)

#---- Report ------------------------------------------------------------------

message("\n=== HFP effect by zone (final model, Jaccard SEs) ===")
print(as.data.frame(hfp_zone), row.names = FALSE, digits = 4)

message("\n=== SE inflation from shared species (final model) ===")
infl <- se_cmp %>% pivot_wider(names_from = estimator, values_from = std_error) %>%
  mutate(inflation = `4. Jaccard` / `1. Normal`)
print(as.data.frame(infl %>% arrange(desc(inflation)) %>% head(8)),
      row.names = FALSE, digits = 4)
message("median inflation: ", round(median(infl$inflation), 3))

message("\n=== sigma and deviance by model (polyseries) ===")
sig_tbl <- coefs %>% filter(link == "polyseries") %>%
  distinct(model, sigma, deviance, df_residual, aic) %>%
  mutate(model = factor(model, levels = col_names)) %>% arrange(model)
print(as.data.frame(sig_tbl), row.names = FALSE, digits = 7)
write_tsv(sig_tbl, file.path(res_dir, "sigma_by_model.tsv"))
if (SIGMA_PER_MODEL) {
  message("NOTE: each model sits at its own sigma, so these deviances are a ",
          "profile comparison, not a nested likelihood-ratio chain.")
}

message("\nwrote ", res_dir, "/{model_coefficients,hfp_effect_by_zone,se_comparison}.tsv")
message("      ", res_dir, if (reuse) "/models.rds (unchanged, reused)" else "/models.rds")

#===============================================================================

DO <- trimws(strsplit(Sys.getenv("DO", "fit,strata,deviance"), ",")[[1]])

#---- Does the model fit? -----------------------------------------------------
# Calibration of observed against predicted, and the observed/predicted ratio
# across the range of sampling effort.

if ("fit" %in% DO) {
  suppressMessages(library(patchwork))
  m_fin <- final_poly$fit
  resp  <- m_fin$model[[1]]                              # the cbind(n, y) matrix
  ff <- tibble(y = as.numeric(resp[, "y"]), n = as.numeric(resp[, "n"]),
               zone = d$zone, mu = as.numeric(fitted(m_fin)))
  stopifnot(all(ff$y == d$y), all(ff$n == d$n))

  pdev <- function(y, mu) 2 * (ifelse(y > 0, y * log(y / pmax(mu, 1e-8)), 0) - (y - mu))
  message("\n=== model fit: ", col_names[length(Models_poly)],
          "  sigma = ", signif(best_sigma, 4), "  n = ", nrow(ff), " ===")
  message(sprintf("  cor(obs, pred)     = %.3f", cor(ff$y, ff$mu)))
  message(sprintf("  cor on log scale   = %.3f", cor(log(ff$y), log(ff$mu))))
  message(sprintf("  deviance explained = %.1f%%",
                  100 * (1 - sum(pdev(ff$y, ff$mu)) / sum(pdev(ff$y, mean(ff$y))))))

  binned <- function(df, v, k = 20) df %>%
    mutate(b = ntile(.data[[v]], k)) %>% group_by(b) %>%
    summarise(x = mean(.data[[v]]), obs = mean(y), pred = mean(mu),
              lo = quantile(y, .25), hi = quantile(y, .75), n_bin = n(),
              .groups = "drop")
  b_mu <- binned(ff, "mu")
  b_n  <- binned(ff, "n")

  fit_p1 <- ggplot(b_mu, aes(pred, obs)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = 2) +
    geom_linerange(aes(ymin = lo, ymax = hi), colour = "grey75", linewidth = .5) +
    geom_point(size = 2.6, colour = "#2166AC") +
    scale_x_log10() + scale_y_log10() +
    labs(x = "predicted richness (bin mean)", y = "observed richness (bin mean)",
         title = "Calibration: observed vs predicted",
         subtitle = "20 equal-count bins of the fitted value; bars = IQR of observed") +
    theme_bw(base_size = 10)

  fit_p2 <- ggplot(b_n, aes(x)) +
    geom_hline(yintercept = 1, colour = "grey40", linetype = 2) +
    geom_point(aes(y = obs / pred), size = 2.6, colour = "#B2182B") +
    geom_line(aes(y = obs / pred), colour = "#B2182B", linewidth = .4) +
    scale_x_log10() + coord_cartesian(ylim = c(0.7, 1.3)) +
    labs(x = "individuals sampled, n (bin mean)", y = "observed / predicted",
         title = "Does the effort curve hold?",
         subtitle = "ratio should sit on 1 across the whole range of sampling effort") +
    theme_bw(base_size = 10)

  fit_p3 <- ggplot(ff, aes(mu, y)) +
    geom_point(alpha = .08, size = .5, stroke = 0) +
    geom_abline(slope = 1, intercept = 0, colour = "#B2182B", linetype = 2) +
    scale_x_log10() + scale_y_log10() + facet_wrap(~ zone, nrow = 1) +
    labs(x = "predicted richness", y = "observed richness",
         title = "All events, by climate zone") +
    theme_bw(base_size = 9)

  message("\n  by zone:")
  print(as.data.frame(ff %>% group_by(zone) %>%
    summarise(obs = mean(y), pred = mean(mu), ratio = mean(y) / mean(mu),
              n = n(), .groups = "drop")), row.names = FALSE, digits = 4)
  message("\n  by n decile (observed/predicted):")
  print(as.data.frame(b_n %>% transmute(n_mean = round(x), n_bin,
    obs = round(obs, 1), pred = round(pred, 1), ratio = round(obs / pred, 3))),
    row.names = FALSE)

  # Print the plot
  (fit_p1 | fit_p2) / fit_p3 + plot_layout(heights = c(1.15, 1))
  # Save the plot
  ggsave(file.path(EXTRA, "model_fit.png"),
         (fit_p1 | fit_p2) / fit_p3 + plot_layout(heights = c(1.15, 1)) +
           plot_annotation(theme = theme(plot.background =
             element_rect(fill = "white", colour = NA))),
         width = 11, height = 8.5, dpi = 200)
  message("\n  -> ", file.path(EXTRA, "model_fit.png"))
}

#---- Does sigma stabilise across strata once covariates are controlled for? ---
# If the covariates absorb the between-stratum differences, the M6/M7 spread
# should be markedly smaller than the M0 spread.

if ("strata" %in% DO) {
  MIN_N  <- as.integer(Sys.getenv("MIN_N", "150"))
  LAT_DF <- as.integer(Sys.getenv("LAT_DF", "6"))

  # M0/M6/M7 for one stratum, omitting whatever is constant inside it.
  strat_formula <- function(dd, strat_var, spec) {
    if (spec == "M0") return("cbind(n, y) ~ 1")
    t <- c("1", "aet")
    if (strat_var != "realm" && nlevels(droplevels(dd$realm)) > 1) t <- c(t, "realm")
    t <- c(t, "fourier_week(week_adj, k = 1)")
    # A latitude spline needs enough distinct values; small strata get fewer df.
    ndf <- min(LAT_DF, max(1, floor(length(unique(dd$Latitude)) / 20)))
    if (ndf >= 2) t <- c(t, sprintf("ns(Latitude, %d)", ndf))
    t <- c(t, if (strat_var != "zone" && nlevels(droplevels(dd$zone)) > 1)
                "zone * hfp" else "hfp")
    # Same rule as realm and zone above: a control with no variation inside this
    # stratum cannot be estimated, so drop it rather than let it fail the fit.
    varies <- function(term) {
      v <- tryCatch(eval(str2lang(term), dd), error = function(e) NULL)
      !is.null(v) && length(unique(v[!is.na(v)])) > 1
    }
    t <- c(t, Filter(varies, ctrl))
    if (spec == "M7") t <- c(t, anml)
    paste("cbind(n, y) ~", paste(t, collapse = " + "))
  }
  strat_sigma <- function(dd, f)
    tryCatch(suppressMessages(estimate_sigma(f, data = droplevels(dd),
                                             startpoint = 0.55)),
             error = function(e) NA_real_)

  message("\n=== sigma by stratum ===")
  strat <- list()
  for (sv in c("zone", "realm")) {
    for (L in sort(unique(as.character(d[[sv]])))) {
      dd <- d[as.character(d[[sv]]) == L, ]
      if (nrow(dd) < MIN_N) {
        message("  skip ", sv, "=", L, " (", nrow(dd), " events)"); next
      }
      for (spec in c("M0", "M6", "M7")) {
        s <- strat_sigma(dd, strat_formula(dd, sv, spec))
        strat[[length(strat) + 1]] <- tibble(stratify_by = sv, stratum = L,
                                             spec = spec, events = nrow(dd), sigma = s)
        message(sprintf("  %-6s %-12s %-3s n=%-5d sigma=%s", sv, L, spec, nrow(dd),
                        ifelse(is.na(s), "FAILED", sprintf("%.4f", s))))
      }
    }
  }
  # The pooled fits use the full specification, zone * hfp included.
  strat_res <- bind_rows(bind_rows(strat), map_dfr(c("M0", "M6", "M7"), function(spec)
    tibble(stratify_by = "pooled", stratum = "all", spec = spec, events = nrow(d),
           sigma = strat_sigma(d, strat_formula(d, "none", spec)))))
  write_tsv(strat_res, file.path(res_dir, "sigma_by_stratum.tsv"))

  message("\n=== sigma by stratum and specification ===")
  print(as.data.frame(strat_res %>% pivot_wider(names_from = spec, values_from = sigma) %>%
          arrange(stratify_by, M0)), row.names = FALSE, digits = 4)
  message("\n=== spread (max - min) within each stratification ===")
  print(as.data.frame(strat_res %>% filter(stratify_by != "pooled") %>%
    group_by(stratify_by, spec) %>%
    summarise(k = n(), spread = max(sigma, na.rm = TRUE) - min(sigma, na.rm = TRUE),
              sd = sd(sigma, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = spec, values_from = c(spread, sd))),
    row.names = FALSE, digits = 3)
}

#---- AET against the competing climate drivers (Table S2) ---------------------
# Marginal explanatory power of each candidate driver, added to progressively
# richer baselines. sigma is held at the main model's value for every cell, so
# the deviances are comparable across the table.

if ("deviance" %in% DO) {
  competitors <- c(AET = "aet", PET = "pet", VPD = "vpd",
                   Temp = "temperature_2m", Precip = "total_precipitation_sum",
                   Humid = "relative_humidity")
  absent <- competitors[!competitors %in% names(d)]
  if (length(absent))
    message("NOT AVAILABLE, skipped: ", paste(names(absent), collapse = ", "))
  competitors <- competitors[competitors %in% names(d)]

  dev_cov <- c("1", "1 + realm",
               "1 + realm + fourier_week(week_adj, k = 1)",
               "1 + realm + fourier_week(week_adj, k = 1) + ns(Latitude, 6)",
               "1 + realm + fourier_week(week_adj, k = 1) + ns(Latitude, 6) + zone * hfp")
  dev_spec <- c("1", "+ Realm", "+ Week", "+ ns(Lat, 6)", "+ Zone x HFP")
  if (length(ctrl)) {
    dev_cov  <- c(dev_cov, paste(dev_cov[5], "+", paste(ctrl, collapse = " + ")))
    dev_spec <- c(dev_spec, "+ Controls")
  }
  if (length(anml) == 4) {
    dev_cov  <- c(dev_cov, paste(dev_cov[length(dev_cov)], "+",
                                 paste(anml, collapse = " + ")))
    dev_spec <- c(dev_spec, "+ Anomalies")
  }

  # Every driver must be complete on the same rows, or the columns of the table
  # are fitted on different samples and their deviances are not comparable.
  ok <- complete.cases(d[, unname(competitors)])
  dd <- d[ok, ]; Jd <- Js[dd$fieldid, dd$fieldid]
  message("\n=== driver deviance table ===")
  message("  rows: ", nrow(dd), " of ", nrow(d), " complete on all ",
          length(competitors), " drivers | specs: ", length(dev_cov))

  # sigma is held fixed across every cell so the deviances share a scale, but it
  # is estimated on the DRIVER-FREE baseline, not on the final model.
  dev_sigma_src <- Sys.getenv("DEV_SIGMA", "baseline")
  if (!is.na(suppressWarnings(as.numeric(dev_sigma_src)))) {
    dev_sigma <- as.numeric(dev_sigma_src)
    message("  sigma given: ", signif(dev_sigma, 7))
  } else if (identical(dev_sigma_src, "main")) {
    dev_sigma <- best_sigma
    message("  sigma from the final model: ", signif(dev_sigma, 7),
            "  (NOT neutral - it was fitted with aet in the model)")
  } else {
    t0 <- Sys.time()
    dev_sigma <- estimate_sigma(paste("cbind(n, y) ~", dev_cov[length(dev_cov)]),
                                data = dd, verbose = FALSE)
    message("  sigma on the driver-free baseline: ", signif(dev_sigma, 7),
            "   (the final model's, with aet, is ", signif(best_sigma, 7), ")   [",
            round(difftime(Sys.time(), t0, units = "secs")), " s]")
  }

  fam      <- hubbell(sigma = dev_sigma)
  dev_null <- deviance(HubbellGLM(cbind(n, y) ~ 1, data = dd, family = fam))
  message("  null deviance: ", format(round(dev_null), big.mark = ","))

  dev_rows <- list()
  for (i in seq_along(dev_cov)) {
    fit0 <- HubbellGLM(as.formula(paste("cbind(n, y) ~", dev_cov[i])),
                       data = dd, family = fam)
    dev0 <- deviance(fit0); r2_base <- 1 - dev0 / dev_null
    message("\n  ", dev_spec[i], "  baseline deviance ",
            format(round(dev0), big.mark = ","), "  pseudo-R2 ", round(r2_base, 4))
    for (nm in names(competitors)) {
      for (j in 1:2) {
        term  <- if (j == 1) competitors[[nm]] else sprintf("ns(%s, 3)", competitors[[nm]])
        shape <- if (j == 1) "Linear" else "Spline"
        fit1 <- tryCatch(HubbellGLM(as.formula(paste("cbind(n, y) ~", dev_cov[i],
                                                     "+", term)),
                                    data = dd, family = fam), error = function(e) NULL)
        if (is.null(fit1)) next
        dev1 <- deviance(fit1)
        p <- tryCatch(waldtest(fit0, fit1, vcov = vcov_shared(fit = fit1, Jd),
                               test = "Chisq")$`Pr(>Chisq)`[2],
                      error = function(e) NA_real_)
        dev_rows[[length(dev_rows) + 1]] <- tibble(
          spec = dev_spec[i], spec_i = i, driver = nm, form = shape,
          dev_base = dev0, dev_model = dev1, dev_drop = dev0 - dev1,
          pseudo_r2 = 1 - dev1 / dev_null,
          added_r2 = (1 - dev1 / dev_null) - r2_base, p_value = p)
        message(sprintf("     %-7s %-7s  drop %9.0f   R2 %.4f   P %.3g",
                        nm, shape, dev0 - dev1, 1 - dev1 / dev_null, p))
      }
    }
  }
  dev_res <- bind_rows(dev_rows) %>%
    mutate(sigma = dev_sigma, sigma_src = dev_sigma_src, n = nrow(dd), zone_col = ZONE)
  write_tsv(dev_res, file.path(res_dir, "TableS2_driver_deviance.tsv"))

  # Cells with P > 0.05 are blanked, as in the published table.
  dev_tab <- dev_res %>%
    mutate(cell = ifelse(!is.na(p_value) & p_value <= 0.05,
                         sprintf("%.3f (%s)", pseudo_r2,
                                 format(round(dev_drop), big.mark = ",")), "")) %>%
    dplyr::select(driver, form, spec, cell) %>%
    pivot_wider(names_from = spec, values_from = cell)
  write_tsv(dev_tab, file.path(res_dir, "TableS2_formatted.tsv"))
  message("\n=== pseudo-R2 (deviance drop), blank where P > 0.05 ===")
  print(as.data.frame(dev_tab), row.names = FALSE)

  #---- ... rendered as the supplement's LaTeX table --------------------------

  spec_levels   <- dev_spec
  driver_levels <- c("AET", "PET", "VPD", "Temp", "Precip", "Humid")
  driver_levels <- driver_levels[driver_levels %in% unique(dev_res$driver)]
  form_levels   <- c("Linear", "Spline")

  r <- dev_res %>% mutate(spec = factor(spec, levels = spec_levels),
                          driver = factor(driver, levels = driver_levels),
                          form = factor(form, levels = form_levels),
                          keep = !is.na(p_value) & p_value <= 0.05)
  # Best (significant) cell per specification - marks both the R2 and the drop.
  r <- r %>%
    left_join(r %>% filter(keep) %>% group_by(spec) %>%
                slice_max(dev_drop, n = 1, with_ties = FALSE) %>% ungroup() %>%
                transmute(spec, driver, form, is_best = TRUE),
              by = c("spec", "driver", "form")) %>%
    mutate(is_best = replace_na(is_best, FALSE))

  cm  <- function(x) formatC(round(x), format = "d", big.mark = ",")
  bf  <- function(s, on) if (on) paste0("\\textbf{", s, "}") else s
  esc <- function(s) str_replace_all(s, fixed("+ Zone x HFP"), "+ Zone$\\times$HFP")
  cell_r2  <- function(x) if (!nrow(x) || !x$keep) "" else
    bf(sprintf("%.3f", x$pseudo_r2), x$is_best)
  cell_dev <- function(x) if (!nrow(x) || !x$keep) "" else
    bf(sprintf("(%s)", cm(x$dev_drop)), x$is_best)

  L <- c(
    "\\begin{table}",
    "\t\\centering",
    paste0("\t\\caption{\\footnotesize\\textbf{Marginal explanatory power (Pseudo-$R^2$) across ",
           "model specifications, with $\\sigma = ", sprintf("%.3f", dev_sigma), "$.} ",
           "As in Table~\\ref{tab:models_poly}, $\\sigma$ is held fixed across all ",
           "specifications at its maximum likelihood value under the most complex model. ",
           "Here, however, this reference model excludes AET and the other candidate drivers ",
           "compared in the table, so that the choice of $\\sigma$ does not favor any of them; ",
           "this is why $\\sigma$ differs slightly from the value $",
           sprintf("%.4f", Models_poly[[length(Models_poly)]]$sigma), "$ used elsewhere. ",
           "The baseline deviance for each specification is reported below its name. ",
           "The marginal drop in deviance for each predictor is reported below the $R^2$ ",
           "in parentheses. Empty cells indicate variables that were not statistically ",
           "significant ($P > 0.05$, Chi-square test with Jaccard-adjusted $P$ values). ",
           "The highest $R^2$ and largest deviance drop per specification are highlighted ",
           "in bold. The column Form indicates the test function form for the covariate ",
           "(linear or natural cubic spline). ",
           "\\textbf{AET}: annual actual evapotranspiration; ",
           "\\textbf{PET}: annual potential evapotranspiration; ",
           "\\textbf{VPD}: annual mean vapor pressure deficit; ",
           "\\textbf{Temp}: annual mean temperature; ",
           "\\textbf{Precip}: annual total precipitation; ",
           "\\textbf{Humid}: annual mean relative humidity. ",
           "$n = ", sub(",", "{,}", cm(nrow(dd)), fixed = TRUE), "$ collection events.}"),
    "\t\\label{tab:deviance_reduction}",
    "\t\\resizebox{0.86\\textwidth}{!}{",
    "\t\\footnotesize",
    paste0("\t\\begin{tabular}{ll *{", length(spec_levels), "}{p{2.00cm}}}"),
    "\t\t\\\\[-0.5em]",
    "\t\t\\hline",
    paste0("\t\t\\textbf{Driver} & \\textbf{Form} & ",
           paste(sprintf("\\textbf{%s}", esc(spec_levels)), collapse = " & "), " \\\\"),
    paste0("\t\t& & ",
           paste(sprintf("(%s)", cm(r %>% distinct(spec, dev_base) %>%
                                      arrange(spec) %>% pull(dev_base))), collapse = " & "),
           " \\\\[0.3em]"),
    "\t\t\\hline\\\\[-0.5em]")

  for (dv in driver_levels) {
    forms <- form_levels[form_levels %in% r$form[r$driver == dv]]
    for (fi in seq_along(forms)) {
      fm  <- forms[fi]
      lab <- if (fi == 1) sprintf("\\textbf{%s}", dv) else ""
      r2  <- vapply(spec_levels, function(sp)
        cell_r2(r %>% filter(driver == dv, form == fm, spec == sp)), character(1))
      ddp <- vapply(spec_levels, function(sp)
        cell_dev(r %>% filter(driver == dv, form == fm, spec == sp)), character(1))
      L <- c(L,
             paste0("\t\t", lab, " & ", fm, " & ", paste(r2, collapse = " & "), " \\\\"),
             paste0("\t\t& & ", paste(ddp, collapse = " & "), " \\\\[0.3em]"))
    }
    L <- c(L, if (dv == tail(driver_levels, 1)) "\t\t\\hline" else "\t\t\\hline\\\\[-0.5em]")
  }
  L <- c(L, "\t\\end{tabular}", "\t}", "\\end{table}")
  writeLines(L, file.path(res_dir, "TableS2_deviance.tex"))
  cat("\n", paste(L, collapse = "\n"), "\n", sep = "")
  message("\nwrote ", file.path(res_dir, "TableS2_deviance.tex"))
}
