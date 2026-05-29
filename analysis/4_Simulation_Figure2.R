# This file runs the 10-fold cross validation analysis on the dataset.
# Running this file reproduces:
# -- Figure 2 of the main paper
# -- Figure S6 of, differences in accumulation curves
library(tidyverse)
library(lubridate)
library(HubbellGLM)
library(BNPvegan)
library(splines)
library(MASS)
library(patchwork)
library(lmtest)
library(ggsignif)
library(RColorBrewer)
library(lmerTest)
library(sf)
library(multcomp)

# This code runs in parallel
library(foreach)
library(doParallel)
registerDoParallel(10)

#================================================
# Useful functions
#================================================
# Function to extract Fisher's alpha
fisher_alpha <- function(y, n, lower = 1e-8, upper = 1e5) {
  # y = richness (S), n = number of individuals (N)
  if (y <= 0 || n <= 0) stop("y and n must be positive.")

  # Equation: y = alpha * log(1 + n/alpha)
  f <- function(alpha) {
    alpha * log(1 + n / alpha) - y
  }

  # Find a root for f(alpha) = 0
  root <- uniroot(f, interval = c(lower, upper))
  return(root$root)
}
fisher_alpha <- Vectorize(fisher_alpha, vectorize.args = c("y", "n"))

richness_alpha <- function(alpha, n){
  alpha * log(1 + n / alpha)
}

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

#-------- Main function to estimate the models
estimate_models <- function(train_data, test_data, str_covariates){
  ###### Hubbell regression with canonical link
  print("Hubbell")
  fitHubbell <- HubbellGLM(formula = paste0("cbind(n, y) ~ ", str_covariates),
                           family = hubbell(sigma = 0),
                           data = train_data)
  # Predict
  #pred_train_Hubbell <- predict.HubbellGLM(fitHubbell, type = "response")
  #pred_test_Hubbell <- predict.HubbellGLM(fitHubbell, type = "response", newdata = test_data)
  pred_train_Hubbell <- predict(fitHubbell, type = "response")
  pred_test_Hubbell <- predict(fitHubbell, type = "response", newdata = test_data)

  ###### Hubbell with polynomial link.
  print("Hubbell Poly")
  best_sigma <- estimate_sigma(paste0("cbind(n, y) ~ ", str_covariates), data = train_data, verbose = TRUE, startpoint = 0.551)
  fitHubbellPoly <- HubbellGLM(formula = paste0("cbind(n, y) ~ ", str_covariates),
                               family = hubbell(sigma = best_sigma),
                               data = train_data)
  # Predict
  #pred_train_HubbellPoly <- predict.HubbellGLM(fitHubbellPoly, type = "response")
  #pred_test_HubbellPoly <- predict.HubbellGLM(fitHubbellPoly, type = "response", newdata = test_data)
  pred_train_HubbellPoly <- predict(fitHubbellPoly, type = "response")
  pred_test_HubbellPoly <- predict(fitHubbellPoly, type = "response", newdata = test_data)

  ###### Poisson with log(n)
  print("Poisson")
  fitPois <- glm(formula = paste0("y - 1 ~ log(n) + ", str_covariates),
                 family = poisson(link = "log"),
                 data = train_data)
  # Predict
  pred_train_Poisson <- predict(fitPois, type = "response") + 1
  pred_test_Poisson <- predict(fitPois, type = "response", newdata = test_data) + 1

  ###### Poisson with offset
  print("Poisson offset")
  fitPois_off <- glm(formula = paste0("y - 1 ~ offset(log(-digamma(1) + log(n))) + ", str_covariates),
                     family = poisson(link = "log"),
                     data = train_data)
  # Predict
  pred_train_Poisson_off <- predict(fitPois_off, type = "response") + 1
  pred_test_Poisson_off <- predict(fitPois_off, type = "response", newdata = test_data) + 1


  ###### Negative binomial with log(n)
  print("NegBin")
  fitNB <- glm.nb(formula = paste0("y - 1 ~ log(n) + ", str_covariates),
                  link = "log",
                  data = train_data)
  # Predict
  pred_train_NB <- predict(fitNB, type = "response") + 1
  pred_test_NB <- predict(fitNB, type = "response", newdata = test_data) + 1

  ###### Negative binomial with offset
  print("NegBin offset")
  fitNB_off <- glm.nb(formula = paste0("y - 1 ~ offset(log(-digamma(1) + log(n))) + ", str_covariates),
                      link = "log",
                      data = train_data)
  # Predict
  pred_train_NB_off <- predict(fitNB_off, type = "response") + 1
  pred_test_NB_off <- predict(fitNB_off, type = "response", newdata = test_data) + 1

  ###### Linear model for alpha diversity
  print("Alpha")
  fitAlpha <- lm(formula = paste0("log(alpha) ~ ", str_covariates),
                 data = train_data %>% mutate(alpha = fisher_alpha(y = y, n = n)))
  # Predict
  pred_train_Alpha <- richness_alpha(alpha = exp(predict(fitAlpha)), n = train_data$n)
  pred_test_Alpha <- richness_alpha(alpha = exp(predict(fitAlpha, newdata = test_data)), n = test_data$n)

  ###### Linear model for richness, with n as natural spline
  print("Richness")
  fitRichness <- lm(formula = paste0("log(y) ~ ns(n, 10) + ", str_covariates),
                    data = train_data)
  # Predict
  pred_train_rich <- exp(predict(fitRichness))
  pred_test_rich <- exp(predict(fitRichness, newdata = test_data))

  # Calculate all RMSE in training at test
  # Calculate all RMSE in training and test
  results_fold <- data.frame(
    model = c("Hubbell", "HubbellPoly", "Poisson", "Poisson_off",
              "NB", "NB_off", "Alpha", "Richness"),
    rmse_train = c(
      sqrt(mean((pred_train_Hubbell      - train_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_train_HubbellPoly  - train_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_train_Poisson      - train_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_train_Poisson_off  - train_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_train_NB           - train_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_train_NB_off       - train_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_train_Alpha        - train_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_train_rich         - train_data$y)^2, na.rm = TRUE))
    ),
    rmse_test = c(
      sqrt(mean((pred_test_Hubbell       - test_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_test_HubbellPoly   - test_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_test_Poisson       - test_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_test_Poisson_off   - test_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_test_NB            - test_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_test_NB_off        - test_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_test_Alpha         - test_data$y)^2, na.rm = TRUE)),
      sqrt(mean((pred_test_rich          - test_data$y)^2, na.rm = TRUE))
    )
  )

  return(results_fold)
}

#--------  Load the data
load("~/HubbellGLM-paper/data/data_GMTP_clean.rdata")
abundance_matrix <- readRDS("~/HubbellGLM-paper/data/abundance_matrix.rds.gzip")

# Create wind_speed column expected by covariates[6]
dataset$wind_speed <- dataset$wind_speed_year_avg

# Covariates for models M1 to M7
covariates <- rep(NA, 7)
covariates[1] <- "aet"
covariates[2] <- paste0(covariates[1], " + realm")
covariates[3] <- paste0(covariates[2], " + fourier_week(week_adj, k = 1)")
covariates[4] <- paste0(covariates[3], " + ns(Latitude, 6)")
covariates[5] <- paste0(covariates[4], " + zone * hfp")
covariates[6] <- paste0(covariates[5], " + wind_speed + wetlands + collection_days + resid_temperature_lat")
covariates[7] <- paste0(covariates[6], " + temperature_2m_anml + total_precipitation_sum_anml + relative_humidity_anml + wind_speed_anml")

#====================================
# Run the simulation
#====================================
resample_folds <- FALSE # <----- Set to TRUE if rerun is needed
if(resample_folds) {
  set.seed(10)
  folds <- caret::createFolds(dataset$y, k = 10, list = TRUE, returnTrain = FALSE)
  #save(folds, file = "~/HubbellGLM-paper/data/simulation_output/folds_random.rdata")
} else {
  load("~/HubbellGLM-paper/data/simulation_output/folds_random.rdata")
}


rerun <- FALSE # <----- Set to TRUE if rerun is needed
if(rerun){
  df_results <- data.frame()
  for(m in 1:length(covariates)){
    print(m)
    str_covariates <- covariates[m]
    df_model <- foreach(i = 1:length(folds), .combine = "rbind") %do% {

      # Find train and test
      test_idx  <- folds[[i]]
      train_idx <- setdiff(seq_len(nrow(dataset)), test_idx)
      train_data <- as.data.frame(dataset[train_idx, ])
      test_data  <- as.data.frame(dataset[test_idx, ])
      # Estimate model
      df_tmp <- estimate_models(train_data, test_data, str_covariates)
      df_tmp$fold <- i
      df_tmp
    }
    df_model$setting <- paste0("M", m)
    df_results <- rbind(df_results, df_model)
  }

  # Save the results
  write_csv(df_results, "~/HubbellGLM-paper/data/simulation_output/Simulation_10fold_models.tsv")
}

#====================================
# Plot the results
#====================================
df_results <- read_csv("~/HubbellGLM-paper/data/simulation_output/Simulation_10fold_models.tsv")

#==================================== Detect the significance
my_comparisons <- list(
  c("HubbellPoly", "Poisson"),
  c("HubbellPoly", "NB"),
  c("HubbellPoly", "Richness")
)
# Example: Poisson=***, Richness=*, NB=ns
my_stars <- c("***", "***", "***")

# Extract the significance from the plots
#--- Hubbell vs Poisson
# Training
fitTrain_HubbPois <-lmerTest::lmer(log2(rmse_train) ~ model + (1|setting) + (1|fold:setting), df_results %>%
                                     dplyr::filter(model %in% c("HubbellPoly", "Poisson")))
summary(fitTrain_HubbPois)
# Test
fitTest_HubbPois <-lmerTest::lmer(log2(rmse_test) ~ model + (1|setting) + (1|fold:setting), df_results %>%
                                    dplyr::filter(model %in% c("HubbellPoly", "Poisson")))
summary(fitTest_HubbPois)

#--- Hubbell vs NB
# Training
fitTrain_HubbNB <-lmerTest::lmer(log2(rmse_train) ~ model + (1|setting) + (1|fold:setting), df_results %>%
                                   dplyr::filter(model %in% c("HubbellPoly", "NB")))
summary(fitTrain_HubbNB)
# Test
fitTest_HubbNB <-lmerTest::lmer(log2(rmse_test) ~ model + (1|setting) + (1|fold:setting), df_results %>%
                                  dplyr::filter(model %in% c("HubbellPoly", "NB")))
summary(fitTest_HubbNB)

#--- Hubbell vs Richness
# Training
fitTrain_HubbRich <-lmerTest::lmer(log2(rmse_train) ~ model + (1|setting) + (1|fold:setting),
                                   df_results %>%
                                     dplyr::filter(model %in% c("HubbellPoly", "Richness")))
summary(fitTrain_HubbRich)
# Test
fitTest_HubbRich <-lmerTest::lmer(log2(rmse_test) ~ model + (1|setting) + (1|fold:setting), df_results %>%
                                    dplyr::filter(model %in% c("HubbellPoly", "Richness")))
summary(fitTest_HubbRich)

#==========================================
# Replicate Figure 2
#==========================================
df_plot <- df_results %>%
  mutate(setting = as.factor(setting),
         model = fct_reorder(model, rmse_train, .fun = median))
my_fills <- setNames(ifelse(levels(df_plot$model) %in% c("HubbellPoly", "Hubbell"), "grey", "white"),
                     levels(df_plot$model))

#----- Training set
star_y_train <- max(log2(df_plot$rmse_train), na.rm = TRUE) - 0.5

p_train <- ggplot(df_plot, aes(x = model, y = log2(rmse_train), fill = model)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, color = "black") +
  geom_line(aes(group = interaction(fold, setting), color = setting),
            alpha = 0.4, linewidth = 0.5) +
  geom_point(aes(color = setting), size = 1.5, alpha = 0.8) +
  geom_signif(
    comparisons = my_comparisons,
    annotations = my_stars,
    y_position = star_y_train,
    step_increase = 0.1,    # <--- This prevents brackets from overlapping!
    tip_length = 0.02,
    vjust = 0.5, color = "gray35", textsize = 5
  ) +
  scale_fill_manual(values = my_fills, guide = "none") +
  theme_classic() +
  labs(y = expression(log[2]("Train RMSE")), color = "Setting") +
  theme(
    axis.title.x = element_blank(),     # Hide x-axis title (Test plot will have it)
    axis.text.x = element_blank(),      # Hide x-axis text to cleanly stack
    axis.ticks.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.5)
  )

#----- Test set
star_y_test <- max(log2(df_plot$rmse_test), na.rm = TRUE) - 0.5

p_test <- ggplot(df_plot,
                 aes(x = model, y = log2(rmse_test), fill = model)) +
  geom_boxplot(width = 0.4, outlier.shape = NA, color = "black") +
  geom_line(aes(group = interaction(fold, setting), color = setting), alpha = 0.4, linewidth = 0.5) +
  geom_point(aes(color = setting), size = 1.5, alpha = 0.8) +
  geom_signif(
    comparisons = my_comparisons,
    annotations = my_stars,
    y_position = star_y_test,
    step_increase = 0.1,    # <--- Stacks the brackets cleanly
    tip_length = 0.02,
    vjust = 0.5, color = "gray35", textsize = 5
  ) +
  scale_fill_manual(values = my_fills, guide = "none") +
  theme_classic() +
  labs(x = "Model (Ordered by Median Test RMSE)",
       y = expression(log[2]("Test RMSE")),
       color = "Setting") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.5)
  )

# ---- Combined plots
combined_plot <- p_train / p_test + plot_layout(guides = "collect")

combined_plot

ggsave(filename = "~/HubbellGLM-paper/figures/Figure2_Hubbell_vs_others_10fold.pdf",
       plot = combined_plot, width = 7.32, height = 4.82)

# ==========================================
# Figure S6
# ==========================================

dataset$wind_speed <- dataset$wind_speed_year_avg
str_covariates <- covariates[7]
fitPois <- glm(formula = paste0("y - 1 ~ log(n) + ", str_covariates),
               family = poisson(link = "log"),
               data = dataset)

fitNB <- glm.nb(formula = paste0("y - 1 ~ log(n) + ", str_covariates),
                link = "log",
                data = dataset)

#0.5586668
#best_sigma <- estimate_sigma(paste0("cbind(n, y) ~ ", str_covariates), data = dataset)
fitHubbellPoly <- HubbellGLM(formula = paste0("cbind(n, y) ~ ", str_covariates),
                             family = hubbell(sigma = 0.5586668),
                             data = as.data.frame(dataset))
vcov_poly <- vcov_shared(fit = fitHubbellPoly, JaccardSim)

# Make the plot
#i <- 1793 # Choose indexes
i <- 2175
nobs <- dataset[i, ]$n
yobs <- dataset[i, ]$y
abundances <- abundance_matrix[dataset$fieldid[i], ]
abundances <- abundances[abundances > 0]
sum(abundances)
length(abundances)
subset_abs <- seq(1, sum(abundances), by = 5)
rar <- BNPvegan::rarefaction(abundances)[subset_abs]

df_pred <- cbind(dataset[i, ] %>% dplyr::select(-n) , n = seq(1, nobs, by = 5))
ypred <- predict(fitPois, newdata = df_pred, type = "response") + 1
ypredNB <- predict(fitNB,  newdata = df_pred, type = "response") + 1
ypredHub <- predict(fitHubbellPoly, newdata = df_pred, type = "response")

res <- predict_curve(fitHubbellPoly, xnew = as.data.frame(dataset[i, ]), n = nobs, npoints = pmin(nobs, 100), .vcov = vcov_poly)
plot(subset_abs, rar, col = "gray50", ylim = c(1, nobs), xlim = c(1, nobs), xlab = "n", ylab = "richness")
abline(a = 0, b = 1, lty = "dotted")
lines(df_pred$n, ypred, type = "l", col = "blue")
lines(df_pred$n, ypredNB, type = "l", col ="forestgreen")
lines(df_pred$n, ypredHub, col = "red")
lines(res$n, (res$mean + 1.96 * res$se)[, 1], col = "red", lty ="dashed")
lines(res$n, (res$mean - 1.96 * res$se)[, 1], col = "red", lty ="dashed")
legend("bottomright",
       legend = c("Rarefaction", "Poisson", "Neg. Binomial", "Hubbell", "Hubbell 95% CI"),
       col    = c("gray50", "blue", "forestgreen", "red", "red"),
       lty    = c(NA, 1, 1, 1, 2),
       pch    = c(1, NA, NA, NA, NA),
       bty    = "n")

