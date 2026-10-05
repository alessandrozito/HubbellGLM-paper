# Cross-validation partitions: random, 500 km blocks, k-means, leave-one-realm-out.
#
# Needs:    01_fit_models.R
# Produces: cv_folds.rds, cv_fold_diagnostics.tsv

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

suppressMessages({library(tidyverse); library(sf); library(blockCV)})

data_dir <- DATA
res_dir  <- SIM
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

K        <- as.integer(Sys.getenv("K", "10"))
REPEATS  <- as.integer(Sys.getenv("REPEATS", "10"))
SEED     <- as.integer(Sys.getenv("SEED", "42"))
BLOCK_M  <- as.numeric(Sys.getenv("BLOCK_M", "500000"))   # metres
# Project before any distance-based partitioning: on raw lon/lat a degree of
# longitude is 111 km at the equator and 19 km at 80 N, which distorts both the
# block grid and k-means exactly where the Polar sites are. Mollweide is
# equal-area and its units are metres, so `size` needs no conversion.
CRS_EQ   <- "ESRI:54009"

#---- The sample: exactly the rows the model was fitted on --------------------

M <- readRDS(file.path(MODELS, "models.rds"))
d <- readRDS(file.path(data_dir, "GMTP_analysis_dataset.rds")) %>%
  dplyr::filter(fieldid %in% M$ids) %>%
  dplyr::distinct(fieldid, .keep_all = TRUE)
d <- d[match(M$ids, d$fieldid), ]
stopifnot(identical(d$fieldid, M$ids))
d$row_id <- seq_len(nrow(d))
message("events: ", format(nrow(d), big.mark = ","),
        "   sites: ", dplyr::n_distinct(d$site_code),
        "   realms: ", dplyr::n_distinct(d$realm))

sites <- d %>% group_by(site_code) %>%
  summarise(Longitude = mean(Longitude), Latitude = mean(Latitude),
            events = dplyr::n(), .groups = "drop")
pts_p <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(CRS_EQ)
site_of_event <- match(d$site_code, sites$site_code)

# site-level fold ids -> one fold id per event
spread_to_events <- function(site_ids) site_ids[site_of_event]

#---- The four schemes --------------------------------------------------------

build <- list(
  random = function(seed) {
    set.seed(seed)
    # 10 folds over events, sizes as equal as the count allows
    sample(rep(seq_len(K), length.out = nrow(d)))
  },
  block500 = function(seed) {
    cv <- cv_spatial(x = pts_p, size = BLOCK_M, k = K, selection = "random",
                     iteration = 50, seed = seed, progress = FALSE,
                     plot = FALSE, report = FALSE)
    spread_to_events(cv$folds_ids)
  },
  kmeans = function(seed) {
    cv <- cv_cluster(x = pts_p, k = K, seed = seed, report = FALSE,
                     progress = FALSE)
    spread_to_events(cv$folds_ids)
  },
  # Deterministic: the realms ARE the folds. Note the benchmark drops `realm`
  # from every specification here - a held-out realm's level cannot be estimated
  # from training data that never contains it.
  realm = function(seed) as.integer(factor(d$realm))
)
N_REP <- c(random = REPEATS, block500 = REPEATS, kmeans = REPEATS, realm = 1L)

folds <- imap_dfr(build, function(f, nm) {
  map_dfr(seq_len(N_REP[[nm]]), function(r) {
    t0 <- Sys.time()
    ids <- f(SEED + r - 1L)
    message(sprintf("  %-9s rep %2d: %d folds, %.0f s", nm, r,
                    dplyr::n_distinct(ids),
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    tibble(scheme = nm, rep = r, row_id = d$row_id, fold = as.integer(ids))
  })
})
saveRDS(list(folds = folds, ids = M$ids, seed = SEED, k = K,
             block_m = BLOCK_M, repeats = N_REP),
        file.path(res_dir, "cv_folds.rds"))

#---- Diagnostics -------------------------------------------------------------
# Great-circle distance from each site to the nearest site in a DIFFERENT fold.
# It is the number that says whether a "spatial" partition actually separates
# anything: if the median is a few km, the test fold has a training point next
# door and the split is spatial in name only.

rad <- pi / 180
dx <- outer(sites$Longitude, sites$Longitude, "-") *
      cos(outer(sites$Latitude, sites$Latitude, "+") / 2 * rad)
dy <- outer(sites$Latitude, sites$Latitude, "-")
D  <- sqrt(dx^2 + dy^2) * 111.32
diag(D) <- Inf

diag_tbl <- folds %>% group_by(scheme, rep) %>% group_modify(function(g, key) {
  ev <- tapply(rep(1L, nrow(g)), g$fold, sum)
  sf_ids <- tapply(g$fold, site_of_event, function(z) z[1])   # site -> fold
  nd <- vapply(seq_len(nrow(sites)), function(i) {
    j <- which(sf_ids != sf_ids[i])
    if (!length(j)) NA_real_ else min(D[i, j])
  }, numeric(1))
  tibble(folds = length(ev), min_events = min(ev), max_events = max(ev),
         imbalance = round(max(ev) / min(ev), 1),
         med_nn_km = round(median(nd, na.rm = TRUE), 1),
         pct_nn_under_50km = round(100 * mean(nd < 50, na.rm = TRUE)))
}) %>% ungroup()
write_tsv(diag_tbl, file.path(res_dir, "cv_fold_diagnostics.tsv"))

message("\n=== partition diagnostics (median over repeats) ===")
print(as.data.frame(diag_tbl %>% group_by(scheme) %>%
  summarise(reps = dplyr::n(), folds = first(folds),
            min_events = median(min_events), max_events = median(max_events),
            imbalance = median(imbalance), med_nn_km = median(med_nn_km),
            pct_nn_under_50km = median(pct_nn_under_50km), .groups = "drop")),
  row.names = FALSE)

#---- The map of the partitions lives in 10_cv_figure.R -----------------------

message("\nwrote cv_folds.rds, cv_fold_diagnostics.tsv",
        "  (partition map: 10_cv_figure.R)")
