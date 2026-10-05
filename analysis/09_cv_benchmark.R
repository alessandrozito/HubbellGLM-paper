# Fit every estimator on every fold; score in- and out-of-sample RMSE.
# Sigma is re-estimated inside each fold. SPECS=all runs M0-M7 (8x slower).
#
# Needs:    08_cv_folds.R
# Produces: cv_benchmark.tsv

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

# One BLAS thread per worker. The BLAS is OpenBLAS-pthread, and a HubbellGLM
# fit does not benefit from it (3.43 s single-threaded vs 3.62 s multi), so
# leaving it free only lets 20 workers fight over 24 cores. Set before the
# cluster starts: PSOCK workers inherit the environment.
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressMessages({library(tidyverse); library(HubbellGLM); library(splines)
                  library(MASS, exclude = c("select")); library(foreach)
                  library(doParallel)})
select <- dplyr::select

data_dir <- DATA
res_dir  <- SIM
CORES     <- as.integer(Sys.getenv("CORES", "10"))
SIGMA_TOL <- as.numeric(Sys.getenv("SIGMA_TOL", "1e-3"))
SIGMA_INT <- as.numeric(strsplit(Sys.getenv("SIGMA_INT", "0.2,0.85"), ",")[[1]])
RERUN    <- Sys.getenv("RERUN", "1") != "0"
SCHEMES  <- trimws(strsplit(Sys.getenv("SCHEMES", "random,block500,kmeans,realm"),
                            ",")[[1]])
# Schemes whose folds can contain a realm the training data never sees, so the
# term has to go. Obvious for leave-one-realm-out; k-means needs it too, because
# a cluster repeatedly swallowed the whole of Australasia (817 events) and every
# estimator then returned non-finite predictions - 18 of 100 folds, silently.
DROP_REALM <- trimws(strsplit(Sys.getenv("DROP_REALM", "realm,kmeans"), ",")[[1]])
out_f    <- file.path(res_dir, "cv_benchmark.tsv")

#---- Specifications, taken from the fitted models so they cannot drift --------

M <- readRDS(file.path(MODELS, "models.rds"))
SPECS <- if (identical(Sys.getenv("SPECS", ""), "all")) seq_along(M$covariates) else
  length(M$covariates)
message("specifications: ", paste(M$col_names[SPECS], collapse = ", "))

d <- readRDS(file.path(data_dir, "GMTP_analysis_dataset.rds")) %>%
  dplyr::filter(fieldid %in% M$ids) %>% dplyr::distinct(fieldid, .keep_all = TRUE)
d <- d[match(M$ids, d$fieldid), ]
stopifnot(identical(d$fieldid, M$ids))
d <- d %>% dplyr::mutate(zone = factor(.data[[ZONE]],
                                       levels = c("Continental", "Dry", "Polar",
                                                  "Temperate", "Tropical")),
                         realm = factor(realm))
d <- as.data.frame(d)   # predict.HubbellGLM indexes newdata[, size] and a
                        # tibble returns a one-column tibble, not a vector

FD <- readRDS(file.path(res_dir, "cv_folds.rds"))
folds <- FD$folds %>% dplyr::filter(scheme %in% SCHEMES)
if (!nrow(folds)) stop("no folds for ", paste(SCHEMES, collapse = ", "),
                       " - run 08_cv_folds.R", call. = FALSE)
stopifnot(identical(FD$ids, M$ids))

#---- Helpers -----------------------------------------------------------------

fisher_alpha <- Vectorize(function(y, n, lower = 1e-8, upper = 1e5) {
  if (y <= 0 || n <= 0) return(NA_real_)
  tryCatch(uniroot(function(a) a * log(1 + n / a) - y,
                   interval = c(lower, upper))$root, error = function(e) NA_real_)
}, vectorize.args = c("y", "n"))

fourier_week <- function(week, k = 1) {
  X <- NULL
  for (j in 1:k) { X <- cbind(X, sin(2*j*pi*week/52), cos(2*j*pi*week/52))
    colnames(X)[(ncol(X)-1):ncol(X)] <- paste0(c("week_sin","week_cos"), j) }
  X
}

# Failures are recorded, not swallowed. The first version of this script lost
# four estimators in every fold to a silent NA, and the only symptom was an
# empty column in the summary table.
ERRS <- new.env(parent = emptyenv())
safe <- function(expr, n, tag = NULL) tryCatch(expr, error = function(e) {
  if (!is.null(tag)) assign(tag, conditionMessage(e), envir = ERRS)
  rep(NA_real_, n)
})
err_of <- function(tag) if (exists(tag, envir = ERRS, inherits = FALSE))
  get(tag, envir = ERRS) else NA_character_

# Keep the test factors on the TRAINING level set, or model.matrix builds a
# different design and predict() fails on a level the fit never saw.
align_levels <- function(test, train) {
  for (fv in c("realm", "zone"))
    if (is.factor(train[[fv]]))
      test[[fv]] <- factor(as.character(test[[fv]]), levels = levels(train[[fv]]))
  test
}

score <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  if (!any(ok)) return(list(n = 0L, sse = NA_real_, rmse = NA_real_,
                            mdare = NA_real_, mpe = NA_real_))
  list(n     = sum(ok),
       sse   = sum((pred[ok] - obs[ok])^2),
       rmse  = sqrt(mean((pred[ok] - obs[ok])^2)),
       mdare = median(abs(pred[ok] - obs[ok]) / obs[ok]),
       mpe   = mean((pred[ok] - obs[ok]) / obs[ok]))
}

# Held-out and in-sample side by side, one row per fold and estimator. In-sample
# is the same eight estimators scored on their OWN training rows, so the gap
# between the two columns is the optimism of each estimator - the quantity that
# says whether a competitor only looks close out-of-sample because it is
# overfitting, and the reason the *_in columns exist at all.
metrics <- function(model, pred, obs, pred_in = NULL, obs_in = NULL,
                    sigma_fold = NA_real_, tag = model) {
  o <- score(pred, obs)
  # An all-NA prediction vector raises no error - predict() simply returns NA
  # for rows whose factor level was never trained - so it has to be named here
  # or the fold disappears from the table with nothing to explain it.
  e <- err_of(tag)
  if (is.na(e) && o$n == 0L) e <- "no finite predictions (unseen factor level?)"
  i <- if (is.null(pred_in)) score(numeric(0), numeric(0)) else score(pred_in, obs_in)
  tibble(model = model, n = o$n, sigma_fold = sigma_fold, err = e,
         sse = o$sse, rmse = o$rmse, mdare = o$mdare, mpe = o$mpe,
         n_in = i$n, sse_in = i$sse, rmse_in = i$rmse,
         mdare_in = i$mdare, mpe_in = i$mpe)
}

fit_fold <- function(train, test, str_cov) {
  rm(list = ls(envir = ERRS), envir = ERRS)
  # as.formula() here, not a character string: a fitter handed a string converts
  # it in its own environment, which inside a parallel worker does not contain
  # fourier_week(). Built here, the formula carries this frame with it.
  ff <- function(txt) stats::as.formula(txt, env = environment())
  f  <- ff(paste0("cbind(n, y) ~ ", str_cov))

  fitH <- safe(HubbellGLM(formula = f, family = hubbell(sigma = 0), data = train), 1, "Hubbell")

  # sigma by golden section rather than HubbellGLM::estimate_sigma's nlminb.
  # Same answer - 0.581194 both ways on the full sample, to six decimals - in 8
  # objective evaluations instead of 10, and each evaluation is a whole IRLS
  # fit, which is 84% of this function's runtime. Note that warm-starting
  # nlminb AT the optimum is a trap: it then takes 40 evaluations, not 10.
  bs   <- tryCatch(
    optimize(function(x) BIC(HubbellGLM(f, family = hubbell(sigma = x), data = train)),
             interval = SIGMA_INT, tol = SIGMA_TOL)$minimum,
    error = function(e) {
      assign("HubbellPoly", conditionMessage(e), envir = ERRS); NA_real_ })
  # A sigma sitting on an endpoint is not an estimate, it is a search that ran
  # out of interval - record it rather than reporting the boundary as a fit.
  if (!is.na(bs) && min(abs(bs - SIGMA_INT)) < 10 * SIGMA_TOL)
    assign("HubbellPoly", sprintf("sigma at search boundary (%.3f)", bs), envir = ERRS)
  fitP <- if (is.na(bs)) NULL else
    safe(HubbellGLM(formula = f, family = hubbell(sigma = bs), data = train), 1, "HubbellPoly")

  off <- "offset(log(-digamma(1) + log(n)))"
  fp  <- safe(glm(ff(paste0("y - 1 ~ log(n) + ", str_cov)), family = poisson("log"),
                  data = train), 1, "Poisson")
  fpo <- safe(glm(ff(paste0("y - 1 ~ ", off, " + ", str_cov)), family = poisson("log"),
                  data = train), 1, "Poisson_off")
  fnb <- safe(MASS::glm.nb(ff(paste0("y - 1 ~ log(n) + ", str_cov)), data = train), 1, "NB")
  fnbo<- safe(MASS::glm.nb(ff(paste0("y - 1 ~ ", off, " + ", str_cov)), data = train), 1, "NB_off")
  fa  <- safe(lm(ff(paste0("log(alpha) ~ ", str_cov)),
                 data = train %>% dplyr::mutate(alpha = fisher_alpha(y = y, n = n)) %>%
                   dplyr::filter(is.finite(alpha), alpha > 0)), 1, "Alpha")
  fr  <- safe(lm(ff(paste0("log(y) ~ ns(n, 10) + ", str_cov)), data = train), 1, "Richness")

  # One prediction routine, called twice: once on the held-out rows and once on
  # the training rows the fits were built from. Nothing is refitted for the
  # in-sample pass, so it costs eight predict() calls, not eight fits.
  preds <- function(nd) {
    ns_ <- nrow(nd)
    gh <- function(m) if (inherits(m, "HubbellGLM"))
      safe(predict(m, type = "response", newdata = nd), ns_) else rep(NA_real_, ns_)
    gp <- function(m, add1 = TRUE) {
      if (is.null(m) || inherits(m, "try-error") || (length(m) == 1 && all(is.na(m))))
        return(rep(NA_real_, ns_))
      v <- safe(predict(m, newdata = nd, type = "response"), ns_)
      if (add1) v + 1 else v
    }
    # Alpha and Richness are fitted on a log scale, so their predictions come
    # back through exp(); Fisher's alpha then has to be turned into an expected
    # count.
    pA <- if (inherits(fa, "lm")) {
      a <- exp(safe(predict(fa, newdata = nd), ns_)); a * log(1 + nd$n / a)
    } else rep(NA_real_, ns_)
    pR <- if (inherits(fr, "lm")) exp(safe(predict(fr, newdata = nd), ns_)) else
      rep(NA_real_, ns_)
    list(HubbellPoly = gh(fitP), Hubbell = gh(fitH),
         Poisson = gp(fp), Poisson_off = gp(fpo),
         NB = gp(fnb), NB_off = gp(fnbo), Alpha = pA, Richness = pR)
  }
  po <- preds(test)
  pi <- preds(train)

  SIG <- c(HubbellPoly = bs, Hubbell = 0)
  dplyr::bind_rows(lapply(names(po), function(m)
    metrics(m, po[[m]], test$y, pi[[m]], train$y,
            sigma_fold = if (m %in% names(SIG)) SIG[[m]] else NA_real_)))
}

#---- Run ---------------------------------------------------------------------

if (!RERUN && file.exists(out_f)) {
  message("RERUN=0 and ", basename(out_f), " exists - nothing to do")
} else {

jobs <- folds %>% count(scheme, rep, fold, name = "test_n") %>%
  dplyr::mutate(train_n = nrow(d) - test_n)
message("jobs: ", nrow(jobs), " folds x ", length(SPECS), " specification(s) x 8 estimators")

cl <- makeCluster(CORES); registerDoParallel(cl)
on.exit({stopCluster(cl); registerDoSEQ()}, add = TRUE)

# Dropping `realm` collapses M2 into M1: "1 + aet + realm" stripped IS
# "1 + aet", bit for bit. Running both wastes a fold set and puts a duplicate
# column in the table that reads as "realm changed nothing" when realm was never
# fitted. Schedule each distinct stripped specification once, per scheme.
spec_plan <- function(scheme) {
  cov <- M$covariates[SPECS]
  if (scheme %in% DROP_REALM) cov <- gsub("\\s*\\+\\s*realm", "", cov, perl = TRUE)
  keep <- SPECS[!duplicated(cov)]
  if (length(keep) < length(SPECS))
    message("  ", scheme, ": ", paste(M$col_names[setdiff(SPECS, keep)], collapse = ", "),
            " collapses onto an earlier specification once realm is dropped - skipped")
  keep
}

all_res <- list()
for (si in SPECS) {
  str_cov <- M$covariates[si]
  # A held-out realm's level cannot be estimated from training data that never
  # contains it, so leave-one-realm-out drops the term - for every estimator
  # equally, which keeps the comparison fair but not comparable across schemes.
  t0 <- Sys.time()
  # Longest jobs first + dynamic scheduling: fold training sets differ by an
  # order of magnitude under the spatial schemes, and static blocks would leave
  # cores idle behind one slow fold.
  jobs <- jobs[order(-jobs$train_n), ]
  res <- foreach(j = seq_len(nrow(jobs)), .combine = dplyr::bind_rows,
                 .options.snow = list(preschedule = FALSE),
                 # foreach auto-exports only what it can SEE in this
                 # expression - fit_fold, align_levels, d, folds, jobs. The
                 # helpers those call, and fourier_week (which appears only
                 # inside the covariate STRING), have to be named or the workers
                 # fail with "could not find function fourier_week".
                 .export = c("metrics", "safe", "err_of", "ERRS",
                             "fisher_alpha", "fourier_week"),
                 .packages = c("dplyr", "tibble", "splines", "HubbellGLM", "MASS")) %dopar% {
    jb <- jobs[j, c("scheme", "rep", "fold")]
    if (!si %in% spec_plan(jb$scheme)) return(NULL)
    key <- folds$scheme == jb$scheme & folds$rep == jb$rep
    test_rows  <- folds$row_id[key & folds$fold == jb$fold]
    train_rows <- folds$row_id[key & folds$fold != jb$fold]
    train <- d[train_rows, ]; test <- align_levels(d[test_rows, ], train)
    sc <- if (jb$scheme %in% DROP_REALM)
      gsub("\\s*\\+\\s*realm", "", str_cov, perl = TRUE) else str_cov
    out <- tryCatch(fit_fold(train, test, sc),
                    error = function(e) tibble(model = NA_character_, n = NA_integer_,
                                               sigma_fold = NA_real_, err = conditionMessage(e),
                                               sse = NA_real_, rmse = NA_real_,
                                               mdare = NA_real_, mpe = NA_real_,
                                               n_in = NA_integer_, sse_in = NA_real_,
                                               rmse_in = NA_real_, mdare_in = NA_real_,
                                               mpe_in = NA_real_))
    dplyr::bind_cols(jb[rep(1, nrow(out)), ], out)
  }
  res$spec <- M$col_names[si]
  all_res[[length(all_res) + 1]] <- res
  message(sprintf("  %-12s done in %.1f min", M$col_names[si],
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  # Written after every specification, so a crash three hours in costs one spec.
  # Rows for the scheme/spec combinations being recomputed are replaced; every
  # other combination already in the file is kept, so one scheme can be re-run
  # on its own. Via a temp file: write_tsv in place would leave a truncated
  # table if the process died mid-write.
  new <- dplyr::bind_rows(all_res)
  if (file.exists(out_f)) {
    old <- readr::read_tsv(out_f, show_col_types = FALSE)
    old <- dplyr::anti_join(old, dplyr::distinct(new, scheme, spec),
                            by = c("scheme", "spec"))
    new <- dplyr::bind_rows(old, new)
  }
  tmp <- paste0(out_f, ".tmp")
  write_tsv(new, tmp); file.rename(tmp, out_f)
}

B <- dplyr::bind_rows(all_res)
message("\nwrote ", out_f, "  (", nrow(B), " rows)")

message("\n=== median RMSE by scheme and estimator (", M$col_names[max(SPECS)],
        ", held-out / in-sample) ===")
print(as.data.frame(B %>% dplyr::filter(spec == M$col_names[max(SPECS)]) %>%
  group_by(scheme, model) %>%
  summarise(folds = dplyr::n(),
            rmse = sprintf("%.1f / %.1f", median(rmse, na.rm = TRUE),
                           median(rmse_in, na.rm = TRUE)),
            failed = sum(is.na(rmse)),
            why = paste(unique(na.omit(err)), collapse = " | "), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = scheme, values_from = c(rmse, failed))),
  row.names = FALSE)
}
