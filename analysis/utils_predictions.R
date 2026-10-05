# Shared helpers for predicting on the global grid. Sourced, not run.

suppressMessages({library(tidyverse); library(HubbellGLM); library(splines)})

ZLV_DEFAULT <- c("Continental", "Dry", "Polar", "Temperate", "Tropical")

fourier_week <- function(week, k = 1) {
  X <- NULL
  for (j in 1:k) {
    X <- cbind(X, sin(2 * j * pi * week / 52), cos(2 * j * pi * week / 52))
    colnames(X)[(ncol(X) - 1):ncol(X)] <- paste0(c("week_sin", "week_cos"), j)
  }
  X
}

#---- The fitted model and the sample it was fitted on -------------------------

load_hubbell_model <- function(models_file = file.path(MODELS, "models.rds"),
                               dataset_file = file.path(DATA, "GMTP_analysis_dataset.rds"),
                               zlv = ZLV_DEFAULT, quiet = FALSE) {
  if (!file.exists(models_file))
    stop("no ", basename(models_file), " - run 01_fit_models.R first", call. = FALSE)
  M   <- readRDS(models_file)
  fit <- M$poly[[length(M$poly)]]$fit

  d <- readRDS(dataset_file) %>%
    dplyr::filter(fieldid %in% M$ids) %>%
    dplyr::distinct(fieldid, .keep_all = TRUE) %>%
    dplyr::mutate(zone  = factor(.data[[ZONE]], levels = zlv),
                  realm = factor(realm))
  d <- d[match(M$ids, d$fieldid), ]
  stopifnot(identical(d$fieldid, M$ids))

  # resid_temperature_lat is a residual against the SAMPLE's latitude fit, so the
  # grid must be residualised against this same fit rather than against its own.
  ok   <- !is.na(d$temperature_2m)
  mlat <- lm(temperature_2m ~ poly(Latitude, 2, raw = TRUE), data = d[ok, ])

  if (!quiet)
    message("model: ", M$col_names[length(M$poly)], "   sigma = ",
            signif(M$sigma, 7), "   fitted on ", nrow(d), " events")

  list(fit = fit, vcov = M$poly[[length(M$poly)]]$vcov,
       vars = all.vars(formula(fit)), sigma = M$sigma,
       label = M$col_names[length(M$poly)], dataset = d, mlat = mlat,
       zlv = zlv, ids = M$ids)
}

#---- One year's grid, ready to predict on -------------------------------------

prepare_grid <- function(year, model, n_pred = 2000, week = 26, days = 7,
                         ras_dir = RASTER, drop_incomplete = FALSE, quiet = FALSE) {
  gf <- file.path(ras_dir, sprintf("Raster_%d_full.rds", year))
  if (!file.exists(gf))
    stop("no ", basename(gf), " - run 11_merge_annual_covariates.R", call. = FALSE)
  g <- readRDS(gf)
  V <- model$vars

  # Scenario constants: a typical sample, mid-season, of a fixed size.
  g$n <- n_pred; g$week_adj <- week; g$collection_days <- days
  # Anomalies are departures of the collection window's weather from that cell's
  # own normal. On the grid there is no collection window, so they are zero -
  # the map is drawn under typical conditions.
  for (v in grep("_anml$", V, value = TRUE)) g[[v]] <- 0

  # Realms the sample never covers cannot be predicted into.
  g$realm[g$realm %in% c("Antarctic", "Oceania")] <- NA
  g$realm <- factor(g$realm, levels = levels(model$dataset$realm))

  # The grid carries both classifications. Predicting under a different one than
  # the coefficients were fitted with is a silent fit/predict mismatch, so the
  # grid follows ZONE exactly as the sample does.
  gz <- if (identical(ZONE, "zone_kgc") && "zone_kgc" %in% names(g)) g$zone_kgc else g$zone
  g$zone <- factor(gz, levels = model$zlv)

  if ("resid_temperature_lat" %in% V)
    g$resid_temperature_lat <- g$temperature_2m -
      predict(model$mlat, newdata = data.frame(Latitude = g$Latitude))

  absent <- setdiff(V, c(names(g), "n", "y"))
  if (length(absent))
    stop("the ", year, " grid is missing model covariates: ",
         paste(absent, collapse = ", "),
         "\n  the maps cannot be produced until these are on the grid", call. = FALSE)

  # `year` rides on the raster itself, so I(year - 2015) evaluates with no help.
  # Assert it: a silently wrong year would shift every prediction.
  stopifnot(all(g$year == year))

  usable <- complete.cases(g[, setdiff(V, "y")])
  if (!quiet)
    message("  ", year, ": ", format(nrow(g), big.mark = ","), " cells, ",
            format(sum(usable), big.mark = ","), " usable (",
            round(100 * mean(usable), 1), "%)")
  if (drop_incomplete) g <- g[usable, ] else g$usable <- usable
  g
}

#---- Predicted indices --------------------------------------------------------
# The hfp + 0.01 nudge is the original scenario definition and is kept: hfp = 0
# is a boundary of the fitted range, and the indices are reported just inside it.

predict_indices <- function(g, model, n_pred = 2000, q = 0.5) {
  fit <- model$fit
  cur <- g %>% dplyr::mutate(hfp = pmin(hfp + 0.01, 50))
  S_sigma  <- exp(predict(fit, newdata = cur)) / fit$sigma
  richness <- predict(fit, newdata = cur, type = "response")
  alpha <- rep(NA_real_, length(richness))
  ok <- is.finite(richness)
  alpha[ok] <- HubbellGLM:::inv_mean_dirichlet_process(mu_target = richness[ok],
                                                       size = n_pred)
  tibble::tibble(S_sigma = S_sigma, richness = richness, alpha = alpha,
                 Shannon = digamma(alpha + 1) - digamma(1),
                 Simpson = 1 / (alpha + 1),
                 Tsallis = 1 / (q - 1) * (1 - alpha * beta(alpha, q)))
}

#---- Counterfactual scenarios -------------------------------------------------
# var_hfp0 is measured against the UN-nudged grid, the others against the nudged
# one, because "if HFP were zero" is a comparison with the world as it is while
# "a 10% increase" is a comparison with the same nudged baseline it perturbs.
# Both conventions come from the published script; keeping them apart matters.

predict_variations <- function(g, model) {
  fit <- model$fit
  p <- function(dd) predict(fit, newdata = dd, type = "response")
  base_raw   <- p(g)
  base_hfp   <- p(g %>% dplyr::mutate(hfp = pmin(hfp + 0.01, 50)))
  base_aet   <- p(g %>% dplyr::mutate(aet = aet + 0.01))
  rich_hfp0  <- p(g %>% dplyr::mutate(hfp = 0))
  tibble::tibble(
    # The two hfp = 0 baselines are returned as well: the richness LOST to human
    # pressure is reported both as a count and as a percentage of the pristine
    # value, and neither is recoverable from var_hfp0 alone.
    richness_raw = base_raw, richness_hfp0 = rich_hfp0,
    var_hfp0  = 100 * (rich_hfp0 - base_raw) / base_raw,
    var_hfp5  = 100 * (p(g %>% dplyr::mutate(hfp = pmin(hfp + 0.01 + 5, 50))) - base_hfp) / base_hfp,
    var_hfp10 = 100 * (p(g %>% dplyr::mutate(hfp = pmin((hfp + 0.01) * 1.1, 50))) - base_hfp) / base_hfp,
    var_aet10 = 100 * (p(g %>% dplyr::mutate(aet = (aet + 0.01) * 1.1)) - base_aet) / base_aet)
}

#---- Everything, for one year -------------------------------------------------

grid_predictions <- function(year, model, n_pred = 2000, week = 26, days = 7,
                             q = 0.5, ras_dir = RASTER, quiet = FALSE,
                             keep = c("zone", "realm", "hfp", "aet", "elevation",
                                      "mess_global", "mod_global", "mess_zone",
                                      "mod_zone", "vars_out_of_range")) {
  g <- prepare_grid(year, model, n_pred = n_pred, week = week, days = days,
                    ras_dir = ras_dir, quiet = quiet)
  dplyr::bind_cols(
    tibble::tibble(Longitude = g$Longitude, Latitude = g$Latitude, year = year),
    g[, intersect(keep, names(g)), drop = FALSE],
    predict_indices(g, model, n_pred = n_pred, q = q),
    predict_variations(g, model))
}
