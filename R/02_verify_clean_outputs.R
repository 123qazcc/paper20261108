#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 120)

if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
library(data.table)

stage_root <- Sys.getenv("GA_STAGE_ROOT", "private/stage")
output_root <- Sys.getenv("GA_OUTPUT_ROOT", file.path(stage_root, "output"))
data_root <- file.path(output_root, "data_clean")
qc_root <- file.path(output_root, "qc")

expected <- c(
  "hrs_clean", "elsa_clean", "charls_clean", "hrs_main_analysis", "elsa_main_analysis",
  "charls_mechanism_fiml", "charls_mechanism_complete3",
  "elsa_mechanism_fiml", "elsa_mechanism_complete3"
)

results <- data.table(check = character(), status = character(), detail = character())
record <- function(check, ok, detail = "") {
  results <<- rbind(results, data.table(
    check = check,
    status = if (isTRUE(ok)) "PASS" else "FAIL",
    detail = as.character(detail)
  ))
  if (!isTRUE(ok)) stop(check, ": ", detail)
}

for (stem in expected) {
  rds_path <- file.path(data_root, paste0(stem, ".rds"))
  csv_path <- file.path(data_root, paste0(stem, ".csv.gz"))
  record(paste0(stem, "_files_exist"), file.exists(rds_path) && file.exists(csv_path))
  x <- readRDS(rds_path)
  record(paste0(stem, "_rds_readable"), is.data.frame(x) && nrow(x) > 0, paste(nrow(x), "rows"))
  record(paste0(stem, "_required_keys"), all(c("cohort", "id", "wave", "event_time") %in% names(x)))
}

counts <- fread(file.path(qc_root, "dataset_counts.csv"))
record("dataset_counts_complete", setequal(counts$dataset, expected))
for (stem in expected) {
  x <- readRDS(file.path(data_root, paste0(stem, ".rds")))
  target <- counts[dataset == stem]
  record(
    paste0(stem, "_counts_match"),
    nrow(x) == target$rows && uniqueN(x$id) == target$ids,
    paste(nrow(x), "rows /", uniqueN(x$id), "people")
  )
}

qc <- fread(file.path(qc_root, "quality_checks.csv"))
record("all_quality_checks_zero", all(qc$value == 0), paste(sum(qc$value), "total violations"))

for (stem in c("hrs_clean", "elsa_clean", "charls_clean")) {
  x <- as.data.table(readRDS(file.path(data_root, paste0(stem, ".rds"))))
  visits <- x[event_type == "visit"]
  record(paste0(stem, "_unique_id_wave"), !anyDuplicated(visits[, .(id, wave)]))
  record(paste0(stem, "_state_range"), all(is.na(x$state) | x$state %in% 0:4))
  setorder(x, id, event_time, wave)
  reversals <- x[event_type == "visit" & !is.na(state), any(diff(state) < 0), by = id][V1 == TRUE, .N]
  record(paste0(stem, "_no_living_state_reversal"), reversals == 0)
  nonlast_deaths <- x[, any(dead == 1L & seq_len(.N) != .N), by = id][V1 == TRUE, .N]
  record(paste0(stem, "_death_absorbing"), nonlast_deaths == 0)
}

for (stem in c("hrs_main_analysis", "elsa_main_analysis")) {
  x <- as.data.table(readRDS(file.path(data_root, paste0(stem, ".rds"))))
  eligibility <- x[, .(
    first_age = age[which(event_type == "visit")[which.min(event_time[event_type == "visit"])]],
    n_states = sum(!is.na(state)),
    n_living = sum(event_type == "visit" & state %in% 0:3),
    n_insomnia = sum(!is.na(insomnia_score)),
    min_time = min(time_since_entry)
  ), by = id]
  record(paste0(stem, "_eligibility"), all(
    eligibility$first_age >= 50 & eligibility$n_states >= 2 & eligibility$n_living >= 1 &
      eligibility$n_insomnia >= 1 & abs(eligibility$min_time) < 1e-10
  ))
}

for (prefix in c("charls", "elsa")) {
  complete <- as.data.table(readRDS(file.path(data_root, paste0(prefix, "_mechanism_complete3.rds"))))
  expected_waves <- if (prefix == "charls") c(1L, 2L, 3L) else c(4L, 6L, 8L)
  complete_check <- complete[, .(
    waves_ok = setequal(wave, expected_waves),
    n_rows = .N,
    complete_sleep = all(!is.na(sleep_hours_observed)),
    complete_sbp = all(!is.na(systolic_bp))
  ), by = id]
  record(paste0(prefix, "_mechanism_complete3_integrity"), all(
    complete_check$waves_ok & complete_check$n_rows == 3L &
      complete_check$complete_sleep & complete_check$complete_sbp
  ))
}

elsa <- as.data.table(readRDS(file.path(data_root, "elsa_clean.rds")))
locf <- elsa[event_type == "visit" & !is.na(insomnia_score)]
record("elsa_insomnia_source_wave_valid", all(locf$insomnia_source_wave %in% c(4L, 6L, 8L)))
record("elsa_observed_sleep_only_w468", all(elsa[!is.na(insomnia_score_observed), wave] %in% c(4L, 6L, 8L)))

fwrite(results, file.path(qc_root, "verification_results.csv"), na = "")
message("Verification passed: ", nrow(results), " checks")
