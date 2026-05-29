# This scripts runs the sequence of models M0 to M7 reported in the main paper
# It replicates
# ---- Table S3 in the supplementary material
# ---- Table S4 in the supplementary material
# ---- Table S2 in the supplementary material
# ---- Figure S5 in the supplementary material

#--- Load the R packages
library(tidyverse)
library(lubridate)
library(HubbellGLM)
library(splines)
library("rnaturalearth")
library("rnaturalearthdata")
library(patchwork)
library(lmtest)
library(sandwich)
library(multcomp)
library(stargazer)

#-------- Load the data
load("~/HubbellGLM-paper/data/data_GMTP_clean.rdata")

#===============================================================================
# Useful functions
#===============================================================================
# Function to extract the fourier series given a certain week of the year
fourier_week <- function(week, k = 1) {
  X <- NULL
  for(j in 1:k) {
    X <- cbind(X, sin(2 * j *pi * week / 52))
    X <- cbind(X, cos(2 * j * pi * week / 52))
    colnames(X)[(ncol(X) - 1):ncol(X)] <- paste0(c("week_sin", "week_cos"), j)
  }
  return(X)
}

# Function to estimate all models
estimate_models <- function(covariates, dataset, sigma = 0, est_sigma = FALSE, similarity = NULL){
  if(est_sigma){
    sigma_used <- estimate_sigma(formula = "cbind(n, y) ~ 1", data = dataset)
  } else {
    sigma_used <- sigma
  }

  if(is.null(similarity)){
    similarity <- diag(nrow(dataset))
  }

  # Sequence of models
  Models <- lapply(1:length(covariates), function(i) {
    message(paste0("Estimate model ", i))
    formula <- paste0("cbind(n, y) ~ ", covariates[i])
    fit <- HubbellGLM(formula = formula, data = dataset, family = hubbell(sigma = sigma_used))
    return(list(fit = fit, vcov = vcov_shared(fit, similarity)))
  })
  return(Models)
}

#===============================================================================
# Run regression models
#===============================================================================
#--- Sequence of models from M0 to M7 as described in the paper
covariates <- rep(NA, 8)
covariates[1] <- "1"
covariates[2] <- paste0(covariates[1], " + aet")
covariates[3] <- paste0(covariates[2], " + realm")
covariates[4] <- paste0(covariates[3], " + fourier_week(week_adj, k = 1)")
covariates[5] <- paste0(covariates[4], " + ns(Latitude, 6)")
covariates[6] <- paste0(covariates[5], " + zone * hfp")
covariates[7] <- paste0(covariates[6], " + wind_speed_year_avg + wetlands + collection_days + resid_temperature_lat")
covariates[8] <- paste0(covariates[7], " + temperature_2m_anml + total_precipitation_sum_anml + relative_humidity_anml + wind_speed_anml")

# Names of the models
col_names <- c("Base", "+aet", "+realm", "+week", "+nsLat", "+zone*hfp", "+controls", "+weather")

#--------------------------------------------
# Part 1 - Polynomial link (the actual model)
#--------------------------------------------
# Estimate sigma in the intercept-only model
run_sigma <- FALSE
if(run_sigma){
  best_sigma <- estimate_sigma("cbind(n, y) ~ 1", data = dataset, startpoint = 0.54, verbose = TRUE)
}
best_sigma <- 0.5692932

# Run models
Models_poly <- estimate_models(covariates, dataset, sigma = best_sigma,
                               est_sigma = FALSE, similarity = JaccardSim)
# Extract the coefficients and the deviance
coefs_poly <- lapply(Models_poly, function(m) coeftest(m$fit, vcov. = m$vcov))
dev_poly <- sapply(Models_poly, function(m) formatC(m$fit$deviance, format="f", digits=1))

# Display the coefficients in all models
stargazer(coefs_poly,
          type = "text", # Change to "latex" for the table
          title = "Stepwise Model Estimates (Polynomial Link)",
          label = "tab:models_poly",
          #keep = c("aet", "zone", "hfp", "Constant"),
          column.labels = col_names,
          add.lines = list(c("Deviance", dev_poly)),
          omit.stat = "all",
          no.space = TRUE,
          font.size = "small",
          dep.var.labels.include = FALSE,
          model.numbers = TRUE)


# Interpretation of the coefficients in the last model
fitHubbell <- HubbellGLM(formula = paste0("cbind(n, y) ~ ", covariates[8]),
                         data = dataset,
                         family = hubbell(sigma = best_sigma))
glht(fitHubbell, linfct = c("aet = 0",
                            "hfp = 0",
                            "hfp + zonePolar:hfp = 0",
                            "hfp + zoneTropical:hfp = 0",
                            "hfp + zoneDry:hfp = 0",
                            "hfp + zoneTemperate:hfp = 0"),
     vcov = vcov_shared(fit = fitHubbell, JaccardSim)) %>%
  summary()


#--------------------------------------------
# Part 2 - Canonical link (further test)
#--------------------------------------------
# Run models
Models_canonical <- estimate_models(covariates, dataset, sigma = 0,
                                    similarity = JaccardSim)
# Extract the coefficients and the deviance
coefs_canonical  <- lapply(Models_canonical,  function(m) coeftest(m$fit, vcov. = m$vcov))
dev_canonical  <- sapply(Models_canonical,  function(m) formatC(m$fit$deviance, format="f", digits=1))

# Display the coefficients in all models
stargazer(coefs_canonical,
          type = "text", # Change to "latex" for the table
          title = "Stepwise Model Estimates (Canonical Link)",
          label = "tab:models_log",
          column.labels = col_names,
          add.lines = list(c("Deviance", dev_canonical)), # Adds deviance neatly at the bottom
          omit.stat = "all",         # Removes clutter
          no.space = TRUE,           # Compacts the rows
          font.size = "small",       # Helps fit 8 columns
          dep.var.labels.include = FALSE,
          model.numbers = TRUE)

# Interpretation of the coefficients in the last model
fitCanonical <- HubbellGLM(formula = paste0("cbind(n, y) ~ ", covariates[8]),
                         data = dataset,
                         family = hubbell(sigma = 0))
glht(fitCanonical, linfct = c("aet = 0",
                            "hfp = 0",
                            "hfp + zonePolar:hfp = 0",
                            "hfp + zoneTropical:hfp = 0",
                            "hfp + zoneDry:hfp = 0",
                            "hfp + zoneTemperate:hfp = 0"),
     vcov = vcov_shared(fit = fitCanonical, JaccardSim)) %>%
  summary()


#===============================================================================
# Figure S5 - Differences between the CIs under different sandwich estimators
#===============================================================================
covariates <- rep(NA, 8)
covariates[1] <- "1"
covariates[2] <- paste0(covariates[1], " + aet")
covariates[3] <- paste0(covariates[2], " + realm")
covariates[4] <- paste0(covariates[3], " + fourier_week(week_adj, k = 1)")
covariates[5] <- paste0(covariates[4], " + ns(Latitude, 6)")
covariates[6] <- paste0(covariates[5], " + zone * hfp")
covariates[7] <- paste0(covariates[6], " + wind_speed_year_avg + wetlands + collection_days + resid_temperature_lat")
covariates[8] <- paste0(covariates[7], " + temperature_2m_anml + total_precipitation_sum_anml + relative_humidity_anml + wind_speed_anml")

# Rescale the two main covariates for visualization
df_scale <- dataset %>%
  mutate(aet = scale(aet),
         hfp = scale(hfp))

# Fit all four models
fitHubbell_log <- HubbellGLM(formula = paste0("cbind(n, y) ~ ", covariates[8]),
                             data = df_scale,
                             family = hubbell(sigma = 0))

fitQuasiHubbell_log <- HubbellGLM(formula = paste0("cbind(n, y) ~ ", covariates[8]),
                                  data = df_scale,
                                  family = quasihubbell(sigma = 0))

fitHubbell <- HubbellGLM(formula = paste0("cbind(n, y) ~ ", covariates[8]),
                         data = df_scale,
                         family = hubbell(sigma = 0.5692932))

fitQuasiHubbell <- HubbellGLM(formula = paste0("cbind(n, y) ~ ", covariates[8]),
                              data = df_scale,
                              family = quasihubbell(sigma = 0.5692932))

# Name the hypothesis tested
hypotheses <- c(
  "AET" = "aet = 0",
  "HFP (Continental)" = "hfp = 0",
  "HFP (Polar)" = "hfp + zonePolar:hfp = 0",
  "HFP (Tropical)" = "hfp + zoneTropical:hfp = 0",
  "HFP (Dry)" = "hfp + zoneDry:hfp = 0",
  "HFP (Temperate)" = "hfp + zoneTemperate:hfp = 0"
)

# 2. Extract the variance-covariance matrices for the Log (Canonical) models
vcov_normal_log <- vcov(fitHubbell_log)
vcov_quasi_log  <- vcov(fitQuasiHubbell_log)
vcov_white_log  <- vcovHC(fitHubbell_log, type = "HC0")
vcov_shared_log <- vcov_shared(fitHubbell_log, JaccardSim)

# 3. Extract the variance-covariance matrices for the Poly models
vcov_normal_poly <- vcov(fitHubbell)
vcov_quasi_poly  <- vcov(fitQuasiHubbell)
vcov_white_poly  <- vcovHC(fitHubbell, type = "HC0")
vcov_shared_poly <- vcov_shared(fitHubbell, JaccardSim)

# 4. Corrected helper function
get_glht_df <- function(model_fit, vcov_mat, model_name, se_name) {
  # Run glht
  test <- glht(model_fit, linfct = hypotheses, vcov = vcov_mat)
  sm <- summary(test)

  # Build dataframe using names(hypotheses) to prevent NAs
  data.frame(
    Term = names(hypotheses), # <--- THIS FIXES THE NA ISSUE
    Estimate = as.numeric(sm$test$coefficients),
    SE = as.numeric(sm$test$sigma),
    Model = model_name,
    SE_Type = se_name,
    stringsAsFactors = FALSE
  )
}

# 5. Extract results for all 8 combinations
results_list <- list(
  get_glht_df(fitHubbell_log, vcov_normal_log, "Canonical (Log)", "1. Normal"),
  get_glht_df(fitQuasiHubbell_log, vcov_quasi_log, "Canonical (Log)", "2. Quasi"),
  get_glht_df(fitHubbell_log, vcov_white_log, "Canonical (Log)", "3. White"),
  get_glht_df(fitHubbell_log, vcov_shared_log, "Canonical (Log)", "4. Jaccard"),

  get_glht_df(fitHubbell, vcov_normal_poly, "Polynomial", "1. Normal"),
  get_glht_df(fitQuasiHubbell, vcov_quasi_poly, "Polynomial", "2. Quasi"),
  get_glht_df(fitHubbell, vcov_white_poly, "Polynomial", "3. White"),
  get_glht_df(fitHubbell, vcov_shared_poly, "Polynomial", "4. Jaccard")
)

# 6. Combine everything, calculate CIs, and order factors
df_plot <- bind_rows(results_list) %>%
  mutate(
    Conf_Low = Estimate - 1.96 * SE,
    Conf_High = Estimate + 1.96 * SE,
    # Reverse factor levels so AET is at the top of the plot
    Term = factor(Term, levels = rev(names(hypotheses)))
  )

# 7. Plot side-by-side
p_rangesCI <- ggplot(df_plot, aes(y = Term, x = Estimate, color = SE_Type, shape = SE_Type)) +
  geom_pointrange(aes(xmin = Conf_Low, xmax = Conf_High),
                  position = position_dodge(width = 0.9),
                  size = 0.3, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  facet_wrap(~ Model) +
  #scale_color_brewer(palette = "Set1") +
  labs(
    x = "Standardized Coefficient Estimate (± 95% CI)",
    y = "Marginal Effect",
    color = "SE Estimator",
    shape = "SE Estimator"
  ) +
  theme_bw(base_size = 14) +
  scale_shape_manual(values = c(0, 1,2, 4))+
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text()
  )
ggsave(plot = p_rangesCI,
       filename = "~/HubbellGLM-paper/figures/FigureS5_Ranges_CI.pdf",
       width = 10.50, height = 3.68)


#===============================================================================
# Table S2 in Supplementary material
#===============================================================================
# Testing the effect that evapotranspiration has the largest drop
# in deviance for every single specification

# 1. Sequence of covariates
covariates <- rep(NA, 7)
covariates[1] <-"1"
covariates[2] <- paste0(covariates[1], " + realm")
covariates[3] <- paste0(covariates[2], " + fourier_week(week_adj, k = 1)")
covariates[4] <- paste0(covariates[3], " + ns(Latitude, 6)")
covariates[5] <- paste0(covariates[4], " + zone * hfp")
covariates[6] <- paste0(covariates[5], " + wind_speed_year_avg + wetlands + collection_days + resid_temperature_lat")
covariates[7] <- paste0(covariates[6], " + temperature_2m_anml + total_precipitation_sum_anml + relative_humidity_anml + wind_speed_anml")

# 2. Explanatory environmental factors for evapotranspiration
competitors <- c("aet", "pet", "vpd_year_avg", "NDVI_year_max",
                 "temperature_2m_year_avg", "total_precipitation_sum_year_avg",
                 "relative_humidity_year_avg")

# 3. Calculate Absolute Null Deviance (Needed for Pseudo R-squared)
fit_null <- HubbellGLM(cbind(n, y) ~ 1, data = dataset,
                       family = hubbell(sigma = best_sigma))
dev_null <- deviance(fit_null)

rerun_table <- TRUE # <---- Set to TRUE if want to rerun
if(rerun_table) {
  # 4. Create empty list to store results
  results_list <- list()

  # 5. Loop through each baseline specification
  for (i in 1:length(covariates)) {

    cat("\nRunning Baseline Spec", i, "...\n")

    # Fit the baseline model (fit0)
    form0 <- as.formula(paste("cbind(n, y) ~", covariates[i]))
    fit0 <- HubbellGLM(form0, data = dataset, family = hubbell(sigma = best_sigma))
    dev0 <- deviance(fit0)
    pseudo_r2_base <- 1 - (dev0 / dev_null)

    # Loop through each competitor variable
    for (var in competitors) {
      cat("\nTesting variable", var, "...\n")
      # Define both the linear and nonlinear (spline) terms to test
      terms_to_test <- c(var, paste0("ns(", var, ", 3)"))
      form_labels <- c("Linear", "Spline (df=3)")

      # Test both shapes for the current variable
      for (j in 1:2) {
        term <- terms_to_test[j]
        shape <- form_labels[j]

        # Fit the updated model (fit1)
        form1 <- as.formula(paste("cbind(n, y) ~", covariates[i], "+", term))
        fit1 <- HubbellGLM(form1, data = dataset, family = hubbell(sigma = best_sigma))
        dev1 <- deviance(fit1)

        # Calculate Deviance reductions
        dev_drop <- dev0 - dev1
        pct_reduction <- (dev_drop / dev0) * 100

        # Calculate Pseudo-R2
        pseudo_r2_model <- 1 - (dev1 / dev_null)
        delta_pseudo_r2 <- pseudo_r2_model - pseudo_r2_base # The isolated R2 added by this variable

        # Run the custom Wald test (Make sure JaccardSim is in your environment)
        wt <- waldtest(fit0, fit1, vcov = vcov_shared(fit = fit1, JaccardSim), test = "Chisq")
        p_val <- wt$`Pr(>Chisq)`[2] # Extract p-value

        # Store the results
        results_list[[length(results_list) + 1]] <- data.frame(
          Baseline_Spec = paste("Spec", i),
          Variable = var,
          Relationship = shape,
          Dev_Drop = round(dev_drop, 2),
          Pct_Dev_Reduction = round(pct_reduction, 4),
          Total_Pseudo_R2 = round(pseudo_r2_model, 4),
          Added_Pseudo_R2 = round(delta_pseudo_r2, 4),
          P_Value = signif(p_val, 4)
        )
      }
    }
  }
  # 6. Combine all results into one dataframe
  final_results <- bind_rows(results_list)
  # Save the results.
  write_tsv(final_results, "~/HubbellGLM-paper/data/TableS2_AET_vs_others_dropDeviance.txt")

}

# Load the dataset with the results
final_results <- read_tsv("~/HubbellGLM-paper/data/TableS2_AET_vs_others_dropDeviance.txt")

summary_table <- final_results %>%
  group_by(Baseline_Spec) %>%
  arrange(Baseline_Spec, desc(Pct_Dev_Reduction)) %>%
  dplyr::mutate(Is_Best_Predictor = ifelse(row_number() == 1, "★ YES ★", "")) %>%
  ungroup() %>%
  dplyr::filter(!(Variable == "aet" & Relationship == "Spline (df=3)")) %>%
  as.data.frame()

# View the final formatted table
data.frame(summary_table)

# Add the deviance of the models
deviance_all <- sapply(1:length(covariates), function(i) {
  form0 <- as.formula(paste("cbind(n, y) ~", covariates[i]))
  fit0 <- HubbellGLM(form0, data = dataset, family = hubbell(sigma = best_sigma))
  deviance(fit0)
})


