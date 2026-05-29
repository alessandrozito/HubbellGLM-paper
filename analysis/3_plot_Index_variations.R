# This file plots the variations in biodiversity indices in our model
#   - Model results
#     -- Figure 4 - simulation scenarios

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
library(MASS)

#-------- Load the data
load("~/HubbellGLM-paper/data/data_GMTP_clean.rdata")

#-------- Load the global raster for year 2014
df_Raster <- readRDS("~/HubbellGLM-paper/data/climate/Raster_2014_merged.rds")

#=============================================
# Useful functions
#=============================================
# Function for the fourier week
fourier_week <- function(week, k = 1) {
  X <- NULL
  for(j in 1:k) {
    X <- cbind(X, sin(2 * j *pi * week / 52))
    X <- cbind(X, cos(2 * j * pi * week / 52))
    colnames(X)[(ncol(X) - 1):ncol(X)] <- paste0(c("week_sin", "week_cos"), j)
  }
  return(X)
}

# Function to extract the CIs for all indices via bootstrap
get_multi_index_ci <- function(fit, newdata, .vcov = NULL,
                               indices = c("reg_div", "richness", "alpha", "Shannon", "Simpson", "Tsallis"),
                               size_n = 2000,
                               n_sims = 500,
                               q = 0.5) {

  # Allow multiple indices to be passed
  allowed_indices <- c("reg_div", "richness", "alpha", "Shannon", "Simpson", "Tsallis")
  indices <- match.arg(indices, allowed_indices, several.ok = TRUE)

  # 1. Setup model matrix for prediction
  Terms <- delete.response(terms(fit))
  Xpred <- model.matrix(Terms, data = newdata, xlev = fit$xlevels)

  # 2. Extract coefficients and sample from Multivariate Normal
  beta_hat <- coef(fit)
  if (is.null(.vcov)) .vcov <- vcov(fit)

  beta_sims <- MASS::mvrnorm(n_sims, mu = beta_hat, Sigma = .vcov)

  # 3. Prepare output storage
  # Create a named list of matrices, one for each requested index
  n_obs <- nrow(newdata)
  sim_mats <- list()
  for (idx in indices) {
    sim_mats[[idx]] <- matrix(NA, nrow = n_obs, ncol = n_sims)
  }

  dummy_fit <- fit

  # Logic checks to avoid unnecessary computations
  needs_richness <- any(c("richness", "alpha", "Shannon", "Simpson", "Tsallis") %in% indices)
  needs_alpha    <- any(c("alpha", "Shannon", "Simpson", "Tsallis") %in% indices)

  cat(sprintf("\nRunning single-pass bootstrap for: %s \n(%d simulations)...\n",
              paste(indices, collapse = ", "), n_sims))
  pb <- txtProgressBar(min = 0, max = n_sims, style = 3)

  # 4. Bootstrap loop (Single Pass)
  for (i in 1:n_sims) {
    beta_i <- beta_sims[i, ]
    eta_i <- as.vector(Xpred %*% beta_i) # Linear Predictor

    # 4a. Calculate S_sigma (reg_div)
    if ("reg_div" %in% indices) {
      sim_mats[["reg_div"]][, i] <- exp(eta_i) / fit$sigma
    }

    # 4b. Calculate Richness (if needed)
    if (needs_richness) {
      dummy_fit$coefficients <- beta_i
      rich_i <- predict(dummy_fit, newdata = newdata, type = "response")

      if ("richness" %in% indices) {
        sim_mats[["richness"]][, i] <- rich_i
      }

      # 4c. Calculate Alpha & downstream indices (if needed)
      if (needs_alpha) {
        alpha_i <- sapply(rich_i, function(r) {
          tryCatch(HubbellGLM:::inv_mean_dirichlet_process(mu_target = r, size = size_n),
                   error = function(e) NA)
        })

        # Populate all requested non-linear indices instantly
        if ("alpha" %in% indices)   sim_mats[["alpha"]][, i]   <- alpha_i
        if ("Shannon" %in% indices) sim_mats[["Shannon"]][, i] <- digamma(alpha_i + 1) - digamma(1)
        if ("Simpson" %in% indices) sim_mats[["Simpson"]][, i] <- 1 / (alpha_i + 1)
        if ("Tsallis" %in% indices) sim_mats[["Tsallis"]][, i] <- alpha_i * beta(alpha_i, q)
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)

  # 5. Extract the 95% CI bounds for ALL requested indices
  results_list <- list()
  for (idx in indices) {
    results_list[[idx]] <- data.frame(
      Lower = apply(sim_mats[[idx]], 1, quantile, probs = 0.025, na.rm = TRUE),
      Upper = apply(sim_mats[[idx]], 1, quantile, probs = 0.975, na.rm = TRUE)
    )
  }

  return(results_list)
}


#=============================================
# Estimate the variations
#=============================================
# Sequence of covariates
covariates <- rep(NA, 8)
covariates[1] <- "1"
covariates[2] <- paste0(covariates[1], " + aet")
covariates[3] <- paste0(covariates[2], " + realm")
covariates[4] <- paste0(covariates[3], " + fourier_week(week_adj, k = 1)")
covariates[5] <- paste0(covariates[4], " + ns(Latitude, 6)")
covariates[6] <- paste0(covariates[5], " + zone * hfp")
covariates[7] <- paste0(covariates[6], " + wind_speed_year_avg + wetlands + collection_days + resid_temperature_lat")
covariates[8] <- paste0(covariates[7], " + temperature_2m_anml + total_precipitation_sum_anml + relative_humidity_anml + wind_speed_anml")

# Estimate Hubbell regression for the last models and the variance-covariance
fitHubbell <- HubbellGLM(formula = paste0("cbind(n, y) ~ ", covariates[8]),
                       data = dataset,
                       family = hubbell(sigma = 0.5692932))
vcov_fit <- vcov_shared(fitHubbell, JaccardSim)


# ---- 1. Find plausible AET ranges per zone and pick a reference realm
df_plausible_zone <- df_Raster %>%
  dplyr::filter(!is.na(zone), !is.na(realm)) %>%
  dplyr::group_by(zone) %>%
  dplyr::summarise(
    Latitude_mean = mean(y, na.rm = TRUE),
    wind_speed_mean = mean(wind_speed, na.rm = TRUE),
    aet_low       = quantile(aet, 0.05, na.rm = TRUE),
    aet_high      = quantile(aet, 0.95, na.rm = TRUE),
    # Pick the most frequently occurring realm in this zone as the reference
    realm         = names(sort(table(realm), decreasing = TRUE))[1],
    .groups       = "drop")

# ---- 2. Expand grid only over Zones
df_out <- expand.grid(
  zone = unique(df_plausible_zone$zone),
  hfp  = c(0, 10, 25, 50),
  aet  = 0:120
) %>%
  as_tibble() %>%
  # Join to get the zone-specific baseline realm, latitude, and aet limits
  left_join(df_plausible_zone, by = "zone") %>%
  filter(aet >= aet_low & aet <= aet_high) %>%
  transmute(
    zone,
    realm,
    hfp,
    aet,
    Latitude                     = Latitude_mean,
    week_adj                     = 26,
    resid_temperature_lat        = 0,
    wetlands                     = 0,
    collection_days              = 7,
    wind_speed_year_avg          = 0,
    temperature_2m_anml          = 0,
    total_precipitation_sum_anml = 0,
    relative_humidity_anml       = 0,
    wind_speed_anml              = 0,
    n = 2000,
    y = 100
  )

df_out <- as.data.frame(df_out)


# ---- 3. Predict the fit and calculate indices
df_out$pred0 <- predict(fitHubbell, newdata = df_out)
df_out$Richness <- predict(fitHubbell, newdata = df_out, type = "response")
df_out$alpha <- HubbellGLM:::inv_mean_dirichlet_process(mu_target = df_out$Richness,
                                                        size = df_out$n)

df_out$Shannon <- digamma(df_out$alpha + 1) - digamma(1)
df_out$Simpson <- 1/(df_out$alpha + 1)
df_out$Tsallis <- df_out$alpha * beta(df_out$alpha, 0.5)
df_out$S_sigma <- exp(df_out$pred0) / fitHubbell$sigma

# Create custom labels for the legend
df_out$hfp_factor <- factor(df_out$hfp,
                            levels = c(0, 25, 50),
                            labels = c("HFP 0", "HFP 25", "HFP 50"))

# Calculate CIs using bootstrap
n_sim_test <- 500
set.seed(10)
ci_all <- get_multi_index_ci(fitHubbell, newdata = df_out, .vcov = vcov_fit,
                             indices = c("reg_div", "alpha", "richness", "Shannon", "Simpson", "Tsallis"),
                             size_n = 2000,
                             n_sims = n_sim_test,
                             q = 0.5)
saveRDS(ci_all, file = "~/HubbellGLM-paper/data/CI_indexes_poly.rds")

# Safely Stack everything into Long Format
build_long_df <- function(index_name, values, ci_df) {
  df_out %>%
    dplyr::select(zone, hfp_factor, aet) %>%
    dplyr::mutate(
      Index = index_name,
      Value = values,
      Lower = ci_df$Lower,
      Upper = ci_df$Upper
    )
}

df_long <- bind_rows(
  build_long_df("S_sigma",  df_out$S_sigma,  ci_all$reg_div),
  build_long_df("alpha",   df_out$alpha,  ci_all$alpha),
  build_long_df("Richness", df_out$Richness, ci_all$richness),
  build_long_df("Shannon",  df_out$Shannon,  ci_all$Shannon),
  build_long_df("Simpson",  df_out$Simpson,  ci_all$Simpson),
  build_long_df("Tsallis",  df_out$Tsallis,  ci_all$Tsallis)
)

# Cap Simpson upper bound to strictly 1.0
df_long$Upper <- ifelse(df_long$Index == "Simpson" & df_long$Upper > 1, 1, df_long$Upper)

# Lock the row order for the plot
df_long$Index <- factor(df_long$Index, levels = c("S_sigma", "Richness", "alpha", "Shannon", "Simpson", "Tsallis"))


pVariations <- ggplot(df_long %>% filter(hfp_factor!="HFP 10"), aes(x = aet, y = Value,
                                                          color = hfp_factor, fill = hfp_factor)) +

  # Add the Confidence Intervals underneath the lines
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, linewidth = 0.2, linetype = "dashed") +

  geom_line(linewidth = 1) +

  facet_grid(Index ~ zone, scales = "free") +

  scale_color_manual(values = c("antiquewhite3", "#CC0000", "#330000"), name = "Human Footprint\n(HFP)") +
  scale_fill_manual(values =  c("antiquewhite3", "#CC0000", "#330000"), name = "Human Footprint\n(HFP)") +

  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  labs(
    x = "Actual Evapotranspiration (AET)",
    y = "Diversity Index Value"
  )

print(pVariations)
ggsave(plot = pVariations,
       filename = "~/HubbellGLM-paper/figures/Figure3_Variations_Index.pdf",
       width = 8.19, height = 5.34)

