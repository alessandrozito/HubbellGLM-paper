# Table S3: stepwise M0-M7 with sigma fixed at the M7 estimate.
#
# Needs:    01_fit_models.R (models.rds)
# Produces: output/models/TableS3_models_poly.tex

source(file.path(path.expand("~/HubbellGLM-paper"), "preprocessing", "paths.R"))
suppressMessages({library(tidyverse); library(Matrix); library(HubbellGLM)
                  library(splines); library(lmtest); library(sandwich)})

fourier_week <- function(week, k = 1) {
  X <- NULL
  for (j in 1:k) { X <- cbind(X, sin(2*j*pi*week/52), cos(2*j*pi*week/52))
    colnames(X)[(ncol(X)-1):ncol(X)] <- paste0(c("week_sin","week_cos"), j) }
  X
}

M  <- readRDS(file.path(MODELS, "models.rds"))
d  <- readRDS(file.path(DATA, "GMTP_analysis_dataset.rds"))
d  <- d[match(M$ids, d$fieldid), ] %>%             # same rows, same order as M7
  mutate(zone = factor(.data[[ZONE]], levels = ZLV), realm = droplevels(factor(realm)))
J  <- readRDS(file.path(DATA, "JaccardSim_sparse.rds"))[M$ids, M$ids]
stopifnot(identical(d$fieldid, M$ids))

covariates <- M$covariates
col_names  <- c("Base","+aet","+realm","+week","+nsLat","+zone*hfp","+controls","+anomalies")
SIGMA      <- M$sigma_by_model[length(M$sigma_by_model)]
message("sigma fixed at M7: ", signif(SIGMA, 7))

Models_poly <- lapply(covariates, function(cv) {
  fit <- HubbellGLM(stats::as.formula(paste("cbind(n, y) ~", cv)), data = d,
                    family = hubbell(sigma = SIGMA))
  list(fit = fit, vcov = vcov_shared(fit, J), sigma = SIGMA)
})

# Table generator shared with 01_fit_models.R.
src <- readLines(file.path(REPO, "analysis", "01_fit_models.R"))
a <- grep("^TERM_GROUPS <- list\\(", src); b <- grep("^stepwise_tex\\(Models_poly", src) - 1
eval(parse(text = src[a:b]))

out <- file.path(MODELS, "TableS3_models_poly.tex")
stepwise_tex(Models_poly, "tab:models_poly",
             "Stepwise Hubbell regression estimates (polynomial link).", out)

tex <- readLines(out)
i <- grep("omitted reference category\\.\\}$", tex)
tex[i] <- sub("\\}$", paste0(" The accumulation-rate parameter $\\\\sigma$ is held fixed at its ",
  "maximum likelihood value in the full specification (M7) for every model, so that ",
  "coefficients are directly comparable across the nested sequence.}"), tex[i])
writeLines(tex, out)
message("wrote ", out)
