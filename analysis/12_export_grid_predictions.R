# Global grid of the mapped predictions (Fig 3, Fig S11), at n = 2000.
# Joins to Raster_YYYY_full.rds by Longitude, Latitude.
#
# Needs:    05_index_curves.R, 06_scenario_maps.R
# Produces: data/rasters_by_year/Raster_YYYY_predictions.rds

source(file.path(path.expand("~/HubbellGLM-paper"), "preprocessing", "paths.R"))
suppressMessages(library(tidyverse))

YEAR     <- REF_YEAR
MASK_CUT <- as.numeric(Sys.getenv("MASK_CUT", "-25"))   # as in 05 and 06
DIGITS   <- 4                                           # significant digits; keeps the file < 50 MB
SCEN     <- c("hfp0", "hfp5", "aet10")

f_idx <- file.path(MODELS, sprintf("grid_indices_%d.rds", YEAR))
f_chg <- file.path(MODELS, sprintf("grid_change_%d.rds",  YEAR))
for (f in c(f_idx, f_chg))
  if (!file.exists(f)) stop("no ", basename(f), " - run 05 and 06 first", call. = FALSE)
I <- readRDS(f_idx); C <- readRDS(f_chg)
stopifnot(identical(I$Longitude, C$Longitude), identical(I$Latitude, C$Latitude))

# Scenario columns: % change in richness (95% CI) and change in number of species
chg <- map(SCEN, function(s)
  setNames(C[paste0(s, c("_est", "_lo", "_hi", "_aest"))],
           paste0(s, c("_pct", "_pct_lo", "_pct_hi", "_abs")))) %>% bind_cols()

out <- bind_cols(
  I %>% transmute(Longitude, Latitude, mess_zone,
                  masked = is.na(mess_zone) | mess_zone < MASK_CUT),   # grey in the maps
  I %>% select(S_sigma, richness, alpha, Shannon, Simpson, Tsallis),
  chg) %>%
  mutate(across(where(is.double) & !c(Longitude, Latitude), ~ signif(unname(.x), DIGITS)))

f_out <- file.path(RASTER, sprintf("Raster_%d_predictions.rds", YEAR))
saveRDS(out, f_out, compress = "xz")
message("wrote ", basename(f_out), ": ", format(nrow(out), big.mark = ","), " cells, ",
        round(file.size(f_out) / 2^20, 1), " MiB")
