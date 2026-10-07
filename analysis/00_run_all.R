# Run the analysis in order, from the data shipped in data/.
#
#   Rscript analysis/00_run_all.R
#   ONLY=5,6 Rscript analysis/00_run_all.R      # a subset
#
# Figures 1 and S4 are assembled by hand from the panels these scripts write.

source(file.path(path.expand("~/HubbellGLM-paper"), "preprocessing", "paths.R"))

STEPS <- c(
  "01_fit_models.R",            # models.rds; Tables S2, S4
  "02_table_S3_fixed_sigma.R",  # Table S3
  "03_descriptive_figures.R",   # Fig 1b, 1c maps; Figs S2, S3, S4, S8, S10
  "04_accumulation_curves.R",   # Fig 1c curves
  "05_index_curves.R",          # Fig 2; Figs S7, S11
  "06_scenario_maps.R",         # Fig 3
  "07_ci_comparison.R",         # Fig S1
  "08_cv_folds.R",              # CV partitions
  "09_cv_benchmark.R",          # CV fits (slow)
  "10_cv_figure.R",             # Figs S5, S9
  "11_admissibility_curves.R",  # Fig S6
  "12_export_grid_predictions.R" # Fig 3 and S11 values on the grid
)

only <- Sys.getenv("ONLY", "")
if (nzchar(only)) {
  pick  <- sprintf("%02d", as.integer(strsplit(only, ",")[[1]]))
  STEPS <- STEPS[substr(STEPS, 1, 2) %in% pick]
}

here <- file.path(REPO, "analysis")
for (s in STEPS) {
  message("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70))
  t0 <- Sys.time()
  if (system2("Rscript", file.path(here, s)) != 0) stop(s, " failed")
  message("-> ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}
message("\ndone")
