suppressPackageStartupMessages(library(data.table))

stage_root <- Sys.getenv("GA_STAGE_ROOT", unset = "private/stage")
output_root <- Sys.getenv(
  "GA_ANALYSIS_OUTPUT_ROOT",
  unset = file.path(stage_root, "analysis_output")
)
derived_root <- file.path(output_root, "data_derived")
qc_root <- file.path(output_root, "qc")
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

results <- list()
add_check <- function(cohort, check, passed, detail) {
  results[[length(results) + 1L]] <<- data.table(
    cohort = cohort,
    check = check,
    status = if (isTRUE(passed)) "PASS" else "FAIL",
    detail = as.character(detail)
  )
}

for (cohort in c("hrs", "elsa")) {
  panel <- as.data.table(readRDS(file.path(derived_root, paste0(cohort, "_exposure_panel.rds"))))
  baseline <- as.data.table(readRDS(file.path(derived_root, paste0(cohort, "_baseline.rds"))))
  cohort_label <- toupper(cohort)

  bad_primary <- panel[!is.na(insomnia_score_primary) & !insomnia_score_primary %in% 0:3, .N]
  add_check(cohort_label, "Primary score is in 0-3", bad_primary == 0L, paste("violations", bad_primary))

  bad_sum <- panel[
    !is.na(insomnia_score_observed) &
      insomnia_score_observed != insomnia_onset_primary + insomnia_maintenance_primary + insomnia_nonrestorative_primary,
    .N
  ]
  add_check(cohort_label, "Observed primary score equals three-domain sum", bad_sum == 0L, paste("violations", bad_sum))

  bad_cumavg <- panel[!is.na(insomnia_score_cumavg) & (insomnia_score_cumavg < 0 | insomnia_score_cumavg > 3), .N]
  add_check(cohort_label, "Cumulative mean is in 0-3", bad_cumavg == 0L, paste("violations", bad_cumavg))

  valid_patterns <- c("no prior assessment", "persistent low", "incident high", "remission", "persistent high")
  bad_pattern <- panel[!is.na(insomnia_dynamic_pattern) & !insomnia_dynamic_pattern %chin% valid_patterns, .N]
  add_check(cohort_label, "Dynamic pattern uses frozen categories", bad_pattern == 0L, paste("violations", bad_pattern))

  bad_code <- panel[
    (!is.na(insomnia_dynamic_code) & !insomnia_dynamic_code %in% 0:3) |
      (insomnia_dynamic_pattern == "persistent low" & insomnia_dynamic_code != 0L) |
      (insomnia_dynamic_pattern == "incident high" & insomnia_dynamic_code != 1L) |
      (insomnia_dynamic_pattern == "remission" & insomnia_dynamic_code != 2L) |
      (insomnia_dynamic_pattern == "persistent high" & insomnia_dynamic_code != 3L),
    .N
  ]
  add_check(cohort_label, "Dynamic codes match labels", bad_code == 0L, paste("violations", bad_code))

  bad_future <- panel[event_type == "visit" & !is.na(insomnia_source_wave) & insomnia_source_wave > wave, .N]
  add_check(cohort_label, "No future insomnia assessment carried backward", bad_future == 0L, paste("violations", bad_future))

  eligible_ids <- panel[analysis_sample == 1L, uniqueN(id)]
  baseline_ids <- uniqueN(baseline$id)
  add_check(cohort_label, "Exactly one baseline per final participant", eligible_ids == baseline_ids && nrow(baseline) == baseline_ids,
            paste("eligible", eligible_ids, "baseline", baseline_ids))

  bad_baseline <- baseline[is.na(insomnia_score_primary) | !state %in% 0:3 | n_state_records_from_exposure_entry < 2L, .N]
  add_check(cohort_label, "Baseline has exposure, state, and follow-up", bad_baseline == 0L, paste("violations", bad_baseline))

  burden_total <- baseline[, sum(insomnia_burden_code %in% 0:2)]
  add_check(cohort_label, "Burden categories cover all final participants", burden_total == nrow(baseline),
            paste("categorized", burden_total, "of", nrow(baseline)))

  state_total <- baseline[, sum(state %in% 0:3)]
  add_check(cohort_label, "Baseline states cover all final participants", state_total == nrow(baseline),
            paste("classified", state_total, "of", nrow(baseline)))

  if (cohort == "hrs") {
    bad_extended <- panel[!is.na(insomnia_score_extended4) & !insomnia_score_extended4 %in% 0:4, .N]
    add_check(cohort_label, "HRS extended score is in 0-4", bad_extended == 0L, paste("violations", bad_extended))
  } else {
    unavailable_ok <- all(is.na(panel$insomnia_score_extended4)) && all(is.na(panel$sleep_medication_current))
    add_check(cohort_label, "ELSA-only unavailable secondary fields remain missing", unavailable_ok,
              "extended score and sleep medication")
  }
}

burden_table <- fread(file.path(derived_root, "table1_by_insomnia_burden.csv"))
state_table <- fread(file.path(derived_root, "table1_by_baseline_state.csv"))
add_check("BOTH", "Burden table contains both cohort panels", setequal(unique(burden_table$cohort), c("HRS", "ELSA")),
          paste(unique(burden_table$cohort), collapse = ", "))
add_check("BOTH", "State table contains both cohort panels", setequal(unique(state_table$cohort), c("HRS", "ELSA")),
          paste(unique(state_table$cohort), collapse = ", "))

verification <- rbindlist(results, use.names = TRUE, fill = TRUE)
fwrite(verification, file.path(qc_root, "exposure_baseline_verification.csv"), na = "")
print(verification)

if (verification[status == "FAIL", .N] > 0L) {
  stop("Exposure/baseline verification failed; inspect exposure_baseline_verification.csv")
}

cat("All exposure/baseline verification checks passed.\n")
