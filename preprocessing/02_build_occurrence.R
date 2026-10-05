# Sparse fieldid x BIN occurrence matrix and the Jaccard similarity matrix that
# vcov_shared() uses for spatially-correlated standard errors.
#
# Produces: occurrence_sparse.rds, JaccardSim_sparse.rds

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

library(tidyverse)
library(Matrix)

raw_zip <- path.expand("~/gmtp_filtered_bcdm.tsv.zip")
data_dir <- DATA
sample_rds <- file.path(data_dir, "GMTP_sample_clean.rds")

MALAISE <- c("Malaise Trap", "Malaise trap")

#---- Jaccard on a sparse incidence matrix ------------------------------------

#' Fraction of shared species between every pair of samples
#'
#' @param B sparse binary incidence matrix, samples in rows.
#' @return symmetric sparse matrix (dsCMatrix), 1 on the diagonal.
get_shared_species_sparse <- function(B) {
  C  <- tcrossprod(B)                 # dsCMatrix: shared-species counts
  dg <- Matrix::diag(C)               # species per sample
  Ct <- methods::as(C, "TsparseMatrix")
  # Jaccard_ij = C_ij / (d_i + d_j - C_ij), evaluated only where C_ij > 0.
  Ct@x <- Ct@x / (dg[Ct@i + 1L] + dg[Ct@j + 1L] - Ct@x)
  methods::as(Ct, "CsparseMatrix")
}

#---- 1. Which fieldids ------------------------------------------------------
message("Reading the analysis sample ...")
samp <- readRDS(sample_rds)
keep <- samp$fieldid
message("  ", length(keep), " fieldids")

#---- 2. Stream the occurrences ----------------------------------------------
orphan_map <- read_tsv(file.path(QC, "qc_records_without_fieldid.tsv"),
                       col_types = cols(.default = col_character())) %>%
  dplyr::select(site_code, coord, collection_date_start, collection_date_end,
                fieldid_new)

message("Reading occurrences from ", basename(raw_zip), " ...")
occ <- read_tsv(
  raw_zip,
  col_select = c(fieldid, bin_uri, sampling_protocol, site_code, coord,
                 collection_date_start, collection_date_end),
  col_types  = cols(.default = col_character())
) %>%
  filter(sampling_protocol %in% MALAISE, !is.na(bin_uri), bin_uri != "None") %>%
  mutate(fieldid = na_if(fieldid, "None")) %>%
  left_join(orphan_map, by = c("site_code", "coord",
                               "collection_date_start", "collection_date_end")) %>%
  # Only a blank fieldid is replaced; a record that already has one keeps it.
  mutate(fieldid = coalesce(fieldid, fieldid_new)) %>%
  filter(fieldid %in% keep) %>%
  # Script 01 drops records further than COORD_RADIUS_KM from their fieldid's
  # modal coordinate, because a handful of fieldids are reused at unrelated
  # sites (Aberdeen + Costa Rica, Louisiana + Israel). Those records must be
  # excluded HERE too, or the foreign BINs stay in the occurrence matrix and
  # the Jaccard similarity is computed from species the sample never caught -
  # while n and y no longer count them. Same exclusion, same key.
  anti_join(read_tsv(file.path(QC, "qc_records_dropped_far.tsv"),
                     col_types = cols(.default = col_character())) %>%
              dplyr::select(fieldid, coord) %>% distinct(),
            by = c("fieldid", "coord")) %>%
  dplyr::select(fieldid, bin_uri) %>%
  distinct()                                   # incidence, not abundance

message("  ", nrow(occ), " distinct fieldid-BIN pairs")

fids <- sort(unique(occ$fieldid))
bins <- sort(unique(occ$bin_uri))
message("  ", length(fids), " fieldids x ", length(bins), " BINs")
if (length(fids) < length(keep)) {
  message("  note: ", length(keep) - length(fids),
          " sampled fieldids have no usable BIN and are absent from the matrix")
}

B <- sparseMatrix(
  i = match(occ$fieldid, fids),
  j = match(occ$bin_uri, bins),
  x = 1,
  dims = c(length(fids), length(bins)),
  dimnames = list(fids, bins)
)
message("  occurrence matrix: ", format(length(B@x), big.mark = ","),
        " nonzeros, ", round(as.numeric(object.size(B)) / 1e6), " MB in memory ",
        "(dense would be ", round(prod(dim(B)) * 8 / 1e9, 1), " GB)")

saveRDS(B, file.path(data_dir, "occurrence_sparse.rds"), compress = "gzip")

#---- 3. Jaccard --------------------------------------------------------------
message("\nComputing the Jaccard similarity ...")
t0 <- Sys.time()
J  <- get_shared_species_sparse(B)
message("  ", round(difftime(Sys.time(), t0, units = "secs")), "s")
message("  nnz ", format(length(J@x), big.mark = ","), " of ",
        format(prod(dim(J)), big.mark = ","),
        " (", round(100 * length(J@x) / prod(dim(J)), 1), "% filled)")
message("  ", round(as.numeric(object.size(J)) / 1e6), " MB in memory ",
        "(dense would be ", round(prod(dim(J)) * 8 / 1e9, 2), " GB)")

# The checks vcov_shared() will make.
stopifnot(all(abs(Matrix::diag(J) - 1) < 1e-12))
Jt  <- methods::as(J, "TsparseMatrix")   # symmetric: one triangle
off <- Jt@i != Jt@j
stopifnot(all(Jt@x[off] >= 0), all(Jt@x[off] < 1))
message("  diagonal is 1, off-diagonal in [0, 1): OK")
message("  max off-diagonal similarity: ", round(max(Jt@x[off]), 4))

saveRDS(J, file.path(data_dir, "JaccardSim_sparse.rds"), compress = "gzip")

message("\nWrote:")
message("  ", file.path(data_dir, "occurrence_sparse.rds"))
message("  ", file.path(data_dir, "JaccardSim_sparse.rds"))
message("\nUse with:  vcov_shared(fit, readRDS('JaccardSim_sparse.rds'))")
message("The rows of the similarity must line up with the rows of the model ",
        "frame, so subset it by fieldid before fitting:")
message("  J <- J[fit_ids, fit_ids]")
