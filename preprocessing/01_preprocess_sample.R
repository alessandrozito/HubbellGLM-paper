# Build the analysis sample from the raw BCDM export: one clean record per
# fieldid, with one latitude/longitude, one collection window and one specimen
# count.
#
# Produces: GMTP_sample_clean.rds, GMTP_allsites.tsv, qc/*.tsv

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

library(tidyverse)
library(geosphere)

#---- Paths -------------------------------------------------------------------
# Read the ZIP, not a loose .tsv. An extracted ~/gmtp_filtered_bcdm.tsv on this
# machine was a truncated 2.3 MB fragment of the real 4.2 GB file; readr reads
# the archive directly and warns that it is skipping the __MACOSX sibling entry.
raw_zip <- path.expand("~/gmtp_filtered_bcdm.tsv.zip")
data_dir <- DATA
qc_dir  <- QC
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

MALAISE <- c("Malaise Trap", "Malaise trap")

#---- Helpers -----------------------------------------------------------------

# Most frequent value of `col` within each fieldid. Ties break on the value
# itself, so the answer does not depend on the order rows happen to arrive in.
resolve_modal <- function(data, col, out_name) {
  data %>%
    count(fieldid, .data[[col]], name = "n_rec") %>%
    arrange(fieldid, desc(n_rec), .data[[col]]) %>%
    group_by(fieldid) %>%
    summarise(
      "{out_name}"              := first(.data[[col]]),
      "{out_name}_n_distinct"   := n(),
      "{out_name}_frac_modal"   := first(n_rec) / sum(n_rec),
      .groups = "drop"
    )
}

# Largest distance (km) between any two coordinates recorded for one fieldid.
max_spread_km <- function(lon, lat) {
  if (length(lon) < 2) return(0)
  max(distm(cbind(lon, lat), fun = distHaversine)) / 1000
}

#---- 1. Load -----------------------------------------------------------------
# Only the columns needed here: the full table is 12.2 M rows x 26 columns.
# Everything is read as character and parsed explicitly, so a stray "None"
# becomes an honest NA instead of silently changing a column's type.
message("Reading ", raw_zip, " ...")
dfGMTP <- read_tsv(
  raw_zip,
  col_select = c(fieldid, site_code, coord, ecoregion, country_iso, bin_uri,
                 collection_date_start, collection_date_end, collection_notes,
                 sampling_protocol),
  col_types  = cols(.default = col_character())
)

df <- dfGMTP %>%
  filter(sampling_protocol %in% MALAISE) %>%
  mutate(across(everything(), ~ na_if(.x, "None")))

message("  ", nrow(dfGMTP), " records, ", nrow(df), " from Malaise traps, ",
        n_distinct(df$fieldid), " fieldids")

#---- 2. Records with no fieldid ----------------------------------------------
# Some records carry no fieldid (blank, or the literal string "None"). Left
# alone they all collapse into one phantom NA "fieldid" spanning many sites.
#
# A trap deployment is identified by site, position and collection window, so
# where an orphan record shares all four with exactly one named fieldid it
# belongs to that deployment and is merged into it. Where it matches none - or
# matches more than one, which makes the choice arbitrary - it is a deployment
# whose label was lost, and gets a synthetic id GMP#NA_1, GMP#NA_2, ...
# Ids are assigned in a fixed order so re-runs are reproducible.
dep_key <- function(d) {
  paste(d$site_code, d$coord, d$collection_date_start, d$collection_date_end,
        sep = "\r")
}
df$.key <- dep_key(df)

named_key <- df %>%
  filter(!is.na(fieldid)) %>%
  distinct(.key, fieldid) %>%
  group_by(.key) %>%
  summarise(n_fieldid = n(), fieldid_match = first(fieldid), .groups = "drop")

orphan <- df %>%
  filter(is.na(fieldid)) %>%
  distinct(.key, site_code, coord, collection_date_start, collection_date_end) %>%
  left_join(named_key, by = ".key") %>%
  mutate(n_fieldid = replace_na(n_fieldid, 0L),
         status = case_when(n_fieldid == 1L ~ "merged",
                            n_fieldid > 1L  ~ "ambiguous",
                            TRUE            ~ "unmatched")) %>%
  arrange(site_code, collection_date_start, collection_date_end)

new_id <- orphan$status != "merged"
orphan$fieldid_new <- ifelse(
  new_id,
  paste0("GMP#NA_", cumsum(new_id)),
  orphan$fieldid_match
)

if (nrow(orphan)) {
  n_orphan_rec <- sum(is.na(df$fieldid))
  message("  ", n_orphan_rec, " records have no fieldid, in ", nrow(orphan),
          " deployments: ", sum(orphan$status == "merged"),
          " merged into an existing fieldid, ", sum(new_id),
          " given a synthetic id (", sum(orphan$status == "unmatched"),
          " matched nothing, ", sum(orphan$status == "ambiguous"),
          " matched more than one fieldid)")
  write_tsv(orphan, file.path(qc_dir, "qc_records_without_fieldid.tsv"))

  lookup <- setNames(orphan$fieldid_new, orphan$.key)
  df <- df %>%
    mutate(fieldid = ifelse(is.na(fieldid), unname(lookup[.key]), fieldid))
}
stopifnot(!any(is.na(df$fieldid)))
df$.key <- NULL

#---- 3. Parse the coordinates ------------------------------------------------
# `coord` arrives as the string "[lat, lon]". Parse the distinct values once and
# join back - str_match over 12 M rows is needlessly slow.
coord_lookup <- tibble(coord = unique(df$coord))
m <- str_match(
  coord_lookup$coord,
  "^\\s*\\[\\s*(-?[0-9.eE+-]+)\\s*,\\s*(-?[0-9.eE+-]+)\\s*\\]\\s*$"
)
coord_lookup$Latitude  <- as.numeric(m[, 2])
coord_lookup$Longitude <- as.numeric(m[, 3])

bad_coord <- coord_lookup %>% filter(!is.na(coord), is.na(Latitude) | is.na(Longitude))
if (nrow(bad_coord)) {
  warning(nrow(bad_coord), " coordinate strings did not parse; see qc_bad_coord.tsv")
  write_tsv(bad_coord, file.path(qc_dir, "qc_bad_coord.tsv"))
}

df <- df %>% left_join(coord_lookup, by = "coord")

# Coordinates must be on the globe; a swapped lat/lon shows up here.
stopifnot(all(is.na(df$Latitude)  | abs(df$Latitude)  <= 90),
          all(is.na(df$Longitude) | abs(df$Longitude) <= 180))

#---- 4. Resolve coordinates and dates per fieldid ----------------------------
# Where a fieldid carries more than one coordinate the extra values are almost
# always a single stray specimen against hundreds of consistent ones, so the
# modal value is the trustworthy one. Section 7 flags every case so the choice
# stays visible.
coord_counts <- df %>% count(fieldid, coord, Latitude, Longitude, name = "n_rec")

# The spread is computed in its own pass.
coord_spread <- coord_counts %>%
  group_by(fieldid) %>%
  summarise(coord_spread_km = max_spread_km(Longitude, Latitude), .groups = "drop")

# Records further than COORD_RADIUS_KM from the fieldid's modal coordinate are DROPPED
COORD_RADIUS_KM <- as.numeric(Sys.getenv("COORD_RADIUS_KM", "5"))

modal_coord <- coord_counts %>%
  arrange(fieldid, desc(n_rec), coord) %>%
  group_by(fieldid) %>%
  summarise(.mlon = first(Longitude), .mlat = first(Latitude), .groups = "drop")

df <- df %>% left_join(modal_coord, by = "fieldid")
df$.dist_km <- NA_real_
ok <- stats::complete.cases(df[, c("Longitude", "Latitude", ".mlon", ".mlat")])
df$.dist_km[ok] <- distHaversine(as.matrix(df[ok, c("Longitude", "Latitude")]),
                                 as.matrix(df[ok, c(".mlon", ".mlat")])) / 1000

far <- which(!is.na(df$.dist_km) & df$.dist_km > COORD_RADIUS_KM)
message("  dropping ", length(far), " records further than ", COORD_RADIUS_KM,
        " km from their fieldid's modal coordinate (",
        dplyr::n_distinct(df$fieldid[far]), " fieldids affected)")
df %>%
  dplyr::slice(far) %>%
  dplyr::count(fieldid, site_code, coord, dist_km = round(.dist_km, 2), name = "records") %>%
  dplyr::arrange(dplyr::desc(dist_km)) %>%
  write_tsv(file.path(qc_dir, "qc_records_dropped_far.tsv"))

if (length(far)) df <- df[-far, ]
df <- df %>% dplyr::select(-.mlon, -.mlat, -.dist_km)

# Recount on the cleaned records so n_coord_distinct and coord_spread_km
# describe what actually went into the sample; the pre-filter spread is kept
# as coord_spread_km_raw so the collisions stay auditable.
coord_spread_raw <- coord_spread %>% dplyr::rename(coord_spread_km_raw = coord_spread_km)
coord_counts <- df %>% count(fieldid, coord, Latitude, Longitude, name = "n_rec")
coord_spread <- coord_counts %>%
  group_by(fieldid) %>%
  summarise(coord_spread_km = max_spread_km(Longitude, Latitude), .groups = "drop")

coord_resolved <- coord_counts %>%
  arrange(fieldid, desc(n_rec), coord) %>%
  group_by(fieldid) %>%
  summarise(
    Latitude          = first(Latitude),
    Longitude         = first(Longitude),
    n_coord_distinct  = n(),
    coord_frac_modal  = first(n_rec) / sum(n_rec),
    .groups = "drop"
  ) %>%
  left_join(coord_spread, by = "fieldid") %>%
  left_join(coord_spread_raw, by = "fieldid")

# Dates: same rule, resolved on the (start, end) pair so a window is never
# assembled from two different collection events.
date_resolved <- df %>%
  count(fieldid, collection_date_start, collection_date_end, name = "n_rec") %>%
  arrange(fieldid, desc(n_rec), collection_date_start, collection_date_end) %>%
  group_by(fieldid) %>%
  summarise(
    collection_start_date = as.Date(first(collection_date_start)),
    collection_end_date   = as.Date(first(collection_date_end)),
    n_dates_distinct      = n(),
    dates_frac_modal      = first(n_rec) / sum(n_rec),
    .groups = "drop"
  )

#---- 5. Resolve the remaining site-level fields ------------------------------
site_resolved  <- resolve_modal(df, "site_code",        "site_code")
eco_resolved   <- resolve_modal(df, "ecoregion",        "ecoregion")
iso_resolved   <- resolve_modal(df, "country_iso",      "country_iso")
# Named to match `most_frequent_note` in 0_GMTP_preprocess_merge.R, which the
# combined_traps / trap_damaged flags are built from.
notes_resolved <- resolve_modal(df, "collection_notes", "most_frequent_note")

effort <- df %>%
  group_by(fieldid) %>%
  summarise(n_records = n(),
            n_bins    = n_distinct(bin_uri[!is.na(bin_uri)]),
            .groups   = "drop")

fieldid_clean <- coord_resolved %>%
  left_join(date_resolved,  by = "fieldid") %>%
  left_join(site_resolved,  by = "fieldid") %>%
  left_join(eco_resolved,   by = "fieldid") %>%
  left_join(iso_resolved,   by = "fieldid") %>%
  left_join(notes_resolved, by = "fieldid") %>%
  left_join(effort,         by = "fieldid") %>%
  mutate(collection_days = as.numeric(collection_end_date - collection_start_date)) %>%
  relocate(fieldid, site_code, Latitude, Longitude,
           collection_start_date, collection_end_date, collection_days,
           n_records, n_bins, ecoregion, country_iso, most_frequent_note)

stopifnot(!any(duplicated(fieldid_clean$fieldid)))

#---- 6. Trap condition from collection_notes ---------------------------------
# WARNING: the note vocabulary in this export does not match the strings that
# 0_GMTP_preprocess_merge.R hard-codes. None of its 15 `combined_traps` strings
# and none of its 15 `trap_damaged` strings occur anywhere in this file -
# "2 Malaise samples combined (+GMP#02510)" is now written
# "2 malaise traps combined (+GMP#01513)". Copying those case_when() lists
# across would flag nothing at all and quietly let combined and damaged traps
# into the analysis.
#
# The field is pipe-delimited, roughly trap|weather|temperature|habitat, and the
# habitat segment describes the forest - phrases like "wind damage" and "fallen
# logs" there refer to the stand, not the trap. Only the first segment is read.
note_trap_segment <- function(x) str_squish(sub("\\|.*$", "", replace_na(x, "")))

DAMAGE_RX <- regex(paste(
  "bear|grizzly|wildlife|vandal|stolen|damag|collaps|fell down|fallen over",
  "blown over|knocked over|punctur|drain|rip in|ripped|torn|bottle empty",
  "jar dry|poles bent|trap fell", sep = "|"), ignore_case = TRUE)

fieldid_clean <- fieldid_clean %>%
  mutate(
    note_trap      = note_trap_segment(most_frequent_note),
    combined_traps = str_detect(note_trap, regex("combin", ignore_case = TRUE)),
    trap_damaged   = str_detect(note_trap, DAMAGE_RX)
  )

message("trap flags from notes: combined = ", sum(fieldid_clean$combined_traps),
        ", damaged = ", sum(fieldid_clean$trap_damaged))

# Audit trail: every distinct trap-note and how it was classified, so the
# regexes above can be checked by eye rather than trusted.
fieldid_clean %>%
  count(note_trap, combined_traps, trap_damaged, name = "n_fieldid") %>%
  arrange(desc(combined_traps | trap_damaged), desc(n_fieldid)) %>%
  write_tsv(file.path(qc_dir, "qc_trap_notes_classified.tsv"))

#---- 7. QC -------------------------------------------------------------------
# Nothing here is dropped automatically. These are the rows where the raw data
# disagreed with itself and the modal value was taken.
qc <- fieldid_clean %>%
  mutate(
    flag_coord_mismatch  = n_coord_distinct > 1,
    # A few hundred metres is GPS/rounding noise; a jump of kilometres is a
    # genuine data-entry error in one of the specimen records.
    flag_coord_far       = coord_spread_km > 1,
    flag_date_mismatch   = n_dates_distinct > 1,
    flag_site_mismatch   = site_code_n_distinct > 1,
    flag_eco_mismatch    = ecoregion_n_distinct > 1,
    flag_iso_mismatch    = country_iso_n_distinct > 1,
    flag_note_mismatch   = most_frequent_note_n_distinct > 1,
    flag_zero_length     = !is.na(collection_days) & collection_days == 0,
    flag_negative_length = !is.na(collection_days) & collection_days < 0,
    flag_missing_coord   = is.na(Latitude) | is.na(Longitude),
    flag_missing_date    = is.na(collection_start_date) | is.na(collection_end_date)
  ) %>%
  filter(if_any(starts_with("flag_")))

message("\n--- QC summary ------------------------------------------------")
message("fieldids resolved:            ", nrow(fieldid_clean))
message("  >1 coordinate:              ", sum(fieldid_clean$n_coord_distinct > 1))
message("    of which spread > 1 km:   ", sum(fieldid_clean$coord_spread_km > 1))
message("    of which spread > 100 km: ", sum(fieldid_clean$coord_spread_km > 100))
message("  >1 collection window:       ", sum(fieldid_clean$n_dates_distinct > 1))
message("  >1 site_code:               ", sum(fieldid_clean$site_code_n_distinct > 1))
message("  >1 ecoregion:               ", sum(fieldid_clean$ecoregion_n_distinct > 1))
message("  >1 country_iso:             ", sum(fieldid_clean$country_iso_n_distinct > 1))
message("  zero-length collections:    ",
        sum(!is.na(fieldid_clean$collection_days) & fieldid_clean$collection_days == 0))
message("  missing coordinate:         ", sum(is.na(fieldid_clean$Latitude)))
message("  missing date:               ", sum(is.na(fieldid_clean$collection_start_date)))
message("unique coordinates:           ",
        nrow(distinct(fieldid_clean, Latitude, Longitude)))
message("unique site_code:             ", n_distinct(fieldid_clean$site_code))
message("---------------------------------------------------------------\n")

write_tsv(qc, file.path(qc_dir, "qc_fieldid_mismatches.tsv"))

# The full per-fieldid list of competing coordinates, for the flagged cases.
coord_counts %>%
  semi_join(filter(fieldid_clean, n_coord_distinct > 1), by = "fieldid") %>%
  arrange(fieldid, desc(n_rec)) %>%
  write_tsv(file.path(qc_dir, "qc_coord_candidates.tsv"))

#---- 8. Site-level check -----------------------------------------------------
# Coordinates are really a property of the site, not the fieldid: in the
# published dataset 2415 fieldids share only 141 distinct coordinates. A site
# whose fieldids disagree spatially is worth a look.
site_check <- fieldid_clean %>%
  filter(!is.na(Latitude)) %>%
  group_by(site_code) %>%
  summarise(n_fieldid       = n(),
            n_coord         = n_distinct(paste(Latitude, Longitude)),
            site_spread_km  = max_spread_km(Longitude, Latitude),
            .groups = "drop") %>%
  arrange(desc(site_spread_km))

write_tsv(site_check, file.path(qc_dir, "qc_site_coord_spread.tsv"))
message("sites with fieldids > 1 km apart: ", sum(site_check$site_spread_km > 1))

#---- 9. Richness counts, filtering and output --------------------------------
# n = individuals sequenced in the sample, y = distinct BINs (richness).
# These were already computed as n_records / n_bins; rename to the names the
# rest of the pipeline and the published dataset use.
fieldid_clean <- fieldid_clean %>%
  rename(n = n_records, y = n_bins)

# The full resolved table, before any exclusion, so the filtering can be
# revisited without re-reading the 4.2 GB export.
saveRDS(fieldid_clean, file.path(data_dir, "GMTP_fieldid_locations_dates.rds"))

# Exclusions, applied one at a time so each one's cost is visible.
#   - combined traps pool two Malaise samples into one record, so n and y are
#     not a single trap's effort
#   - damaged traps lost part of their catch
#   - n > 20 follows the published preprocessing/00_newdata.R (which dataset used n >= 20)
#   - y < n drops samples where every individual is its own BIN; there the
#     accumulation curve is saturated and alpha is not identified
steps <- tibble(
  step = c("resolved fieldids", "drop combined traps", "drop damaged traps",
           "n > 20", "y < n"),
  kept = NA_integer_
)
d <- fieldid_clean;                        steps$kept[1] <- nrow(d)
d <- d %>% filter(!combined_traps);        steps$kept[2] <- nrow(d)
d <- d %>% filter(!trap_damaged);          steps$kept[3] <- nrow(d)
d <- d %>% filter(n > 20);                 steps$kept[4] <- nrow(d)
d <- d %>% filter(y < n);                  steps$kept[5] <- nrow(d)
steps$dropped <- c(NA_integer_, -diff(steps$kept))

dataset_new <- d %>%
  dplyr::select(fieldid, site_code, Latitude, Longitude,
                collection_start_date, collection_end_date, collection_days,
                n, y, ecoregion, country_iso, most_frequent_note,
                combined_traps, trap_damaged,
                n_coord_distinct, coord_spread_km, n_dates_distinct)

message("\n--- filtering ------------------------------------------------")
for (i in seq_len(nrow(steps))) {
  message(sprintf("  %-22s %6d%s", steps$step[i], steps$kept[i],
                  if (is.na(steps$dropped[i])) "" else
                    sprintf("   (-%d)", steps$dropped[i])))
}
message("--------------------------------------------------------------")
message("final sample: ", nrow(dataset_new), " fieldids, ",
        n_distinct(dataset_new$site_code), " sites, ",
        nrow(distinct(dataset_new, Latitude, Longitude)), " coordinates")
message("  n range: ", min(dataset_new$n), " - ", max(dataset_new$n),
        " | y range: ", min(dataset_new$y), " - ", max(dataset_new$y))
message("  dates:   ", format(min(dataset_new$collection_start_date)), " to ",
        format(max(dataset_new$collection_end_date)))

write_tsv(dataset_new, file.path(data_dir, "GMTP_allsites.tsv"))
saveRDS(dataset_new, file.path(data_dir, "GMTP_sample_clean.rds"))
write_tsv(steps, file.path(qc_dir, "qc_filter_steps.tsv"))

message("\nWrote ", file.path(data_dir, "GMTP_fieldid_locations_dates.rds"),
        " (all ", nrow(fieldid_clean), " fieldids)")
message("Wrote ", file.path(data_dir, "GMTP_sample_clean.rds"),
        " (filtered sample)")
message("Wrote ", file.path(data_dir, "GMTP_allsites.tsv"))
message("QC files in ", qc_dir)

