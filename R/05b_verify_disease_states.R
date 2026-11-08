#!/usr/bin/env Rscript

# Independent disease-state audit for the HRS/ELSA multi-state analysis.
# This script deliberately recomputes cumulative disease histories and states
# from the three component diagnoses instead of trusting the stored state field.

suppressPackageStartupMessages({
  library(data.table)
  library(GlobalAgingR)
})

options(stringsAsFactors = FALSE, width = 160)

stage_root <- Sys.getenv("GA_STAGE_ROOT", unset = "private/stage")
clean_root <- Sys.getenv(
  "GA_CLEAN_INPUT_ROOT",
  unset = file.path(stage_root, "output", "data_clean")
)
exposure_root <- Sys.getenv(
  "GA_EXPOSURE_INPUT_ROOT",
  unset = file.path(stage_root, "analysis_output", "data_derived")
)
output_root <- Sys.getenv(
  "GA_MULTISTATE_OUTPUT_ROOT",
  unset = file.path(stage_root, "multistate_output")
)
data_root <- file.path(output_root, "data_multistate")
result_root <- file.path(output_root, "results")
qc_root <- file.path(output_root, "qc")
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

state_labels <- c(
  `0` = "S0: no hypertension, heart disease, or stroke",
  `1` = "S1: hypertension only",
  `2` = "S2: exactly one of heart disease or stroke, without hypertension",
  `3` = "S3: at least two of hypertension, heart disease, and stroke",
  `4` = "S4: death (absorbing)"
)

expected_q_edges <- c("0->1", "0->2", "0->4", "1->3", "1->4", "2->3", "2->4", "3->4")
expected_panel_edges <- c(expected_q_edges, "0->3")

same_with_na <- function(x, y) {
  (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
}

cum_diagnosis_independent <- function(x) {
  x <- as.integer(x)
  out <- rep(NA_integer_, length(x))
  observed <- FALSE
  ever <- FALSE
  for (i in seq_along(x)) {
    if (!is.na(x[i])) {
      observed <- TRUE
      if (x[i] == 1L) ever <- TRUE
    }
    if (ever) {
      out[i] <- 1L
    } else if (observed) {
      out[i] <- 0L
    }
  }
  out
}

derive_state_independent <- function(hypertension, heart_disease, stroke, dead) {
  h <- as.integer(hypertension)
  c <- as.integer(heart_disease)
  s <- as.integer(stroke)
  d <- as.integer(dead)
  components <- cbind(h, c, s)
  yes_count <- rowSums(components == 1L, na.rm = TRUE)
  known_count <- rowSums(!is.na(components))
  out <- rep(NA_integer_, length(h))
  out[d == 1L] <- 4L
  alive <- is.na(d) | d != 1L
  complete <- known_count == 3L
  out[alive & complete & yes_count == 0L] <- 0L
  out[alive & complete & h == 1L & c == 0L & s == 0L] <- 1L
  out[alive & complete & h == 0L & (c + s) == 1L] <- 2L
  # Two observed positive diagnoses are sufficient to establish multimorbidity,
  # even if the third diagnosis is missing.
  out[alive & yes_count >= 2L] <- 3L
  out
}

yes1_no2_independent <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_integer_, length(x))
  out[x == 1] <- 1L
  out[x == 2] <- 0L
  out
}

binary01_independent <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_integer_, length(x))
  out[x == 0] <- 0L
  out[x == 1] <- 1L
  out
}

combine_binary_any_independent <- function(...) {
  fields <- cbind(...)
  out <- rep(NA_integer_, nrow(fields))
  any_yes <- rowSums(fields == 1L, na.rm = TRUE) > 0L
  all_known_no <- rowSums(!is.na(fields)) == ncol(fields) &
    rowSums(fields == 0L, na.rm = TRUE) == ncol(fields)
  out[any_yes] <- 1L
  out[!any_yes & all_known_no] <- 0L
  out
}

is_explicit_dispute_independent <- function(x) {
  suppressWarnings(as.numeric(x)) %in% c(2, 3, 4)
}

ever_with_feedforward_independent <- function(ever, prior_confirmation) {
  ever <- suppressWarnings(as.numeric(ever))
  prior_confirmation <- suppressWarnings(as.numeric(prior_confirmation))
  out <- rep(NA_integer_, length(ever))
  out[ever == 1 | prior_confirmation == 1] <- 1L
  rejected_or_not_routed <- prior_confirmation %in% c(-2, -1, 2, 3, 4)
  out[ever == 2 & rejected_or_not_routed] <- 0L
  out
}

derive_elsa_w10_heart_components_independent <- function(x) {
  mi_route <- cbind(
    suppressWarnings(as.numeric(x$HEAgeMI)),
    suppressWarnings(as.numeric(x$HEOthHA)),
    suppressWarnings(as.numeric(x$HEMedHA))
  )
  mi_route_observed <- apply(mi_route, 1L, function(value) {
    any(!is.na(value) & !(value %in% c(-1, -2)))
  })
  mi_route_not_applicable <- apply(mi_route, 1L, function(value) {
    all(!is.na(value) & value %in% c(-1, -2))
  })
  data.table(
    expected_heart_angina_raw = ever_with_feedforward_independent(x$HEEverAN, x$HENvrAN),
    expected_heart_mi_raw = fcase(
      mi_route_observed, 1L,
      mi_route_not_applicable, 0L,
      default = NA_integer_
    ),
    expected_heart_failure_raw = ever_with_feedforward_independent(x$HEEverHF, x$HENvrHF),
    expected_heart_murmur_raw = ever_with_feedforward_independent(x$HEEverHM, x$HENvrHM),
    expected_heart_rhythm_raw = ever_with_feedforward_independent(x$HEEverAH, x$HENvrAH),
    expected_heart_other_raw = ever_with_feedforward_independent(x$HEEverOH, x$HENvrOH)
  )
}

derive_elsa_w10_heart_independent <- function(x) {
  components <- derive_elsa_w10_heart_components_independent(x)
  combine_binary_any_independent(
    components$expected_heart_angina_raw,
    components$expected_heart_mi_raw,
    components$expected_heart_failure_raw,
    components$expected_heart_murmur_raw,
    components$expected_heart_rhythm_raw,
    components$expected_heart_other_raw
  )
}

checks <- list()
add_check <- function(cohort, source, check_id, criterion, violations, detail = "") {
  violations <- as.integer(violations)
  checks[[length(checks) + 1L]] <<- data.table(
    cohort = cohort,
    source = source,
    check_id = check_id,
    criterion = criterion,
    violations = violations,
    status = ifelse(violations == 0L, "PASS", "FAIL"),
    detail = as.character(detail)
  )
}

combination_tables <- list()
component_tables <- list()
w10_heart_audit <- data.table()
elsa_backcast_audit <- data.table()
w10_source_rows <- list()

audit_dataset <- function(cohort, source_name, path) {
  if (!file.exists(path)) {
    add_check(cohort, source_name, "INPUT_EXISTS", "Input dataset exists", 1L, path)
    return(invisible(NULL))
  }
  x <- as.data.table(readRDS(path))
  required <- c(
    "id", "wave", "event_time", "event_type", "dead", "state",
    "hypertension_raw", "heart_disease_raw", "stroke_raw",
    "hypertension", "heart_disease", "stroke"
  )
  missing_fields <- setdiff(required, names(x))
  add_check(
    cohort, source_name, "REQUIRED_FIELDS", "All disease-state audit fields are present",
    length(missing_fields), paste(missing_fields, collapse = "; ")
  )
  if (length(missing_fields)) return(invisible(NULL))

  x[, audit_row_id := .I]
  setorder(x, id, event_time, wave, na.last = TRUE)
  duplicate_visit <- x[event_type == "visit", sum(duplicated(paste(id, wave)))]
  add_check(cohort, source_name, "UNIQUE_ID_WAVE", "Living visits are unique by participant and wave", duplicate_visit)

  for (component in c("hypertension", "heart_disease", "stroke")) {
    raw_field <- paste0(component, "_raw")
    raw_bad <- x[event_type == "visit" & !is.na(get(raw_field)) & !get(raw_field) %in% 0:1, .N]
    cumulative_bad <- x[event_type == "visit" & !is.na(get(component)) & !get(component) %in% 0:1, .N]
    add_check(cohort, source_name, paste0(toupper(component), "_RAW_DOMAIN"),
              paste0(raw_field, " is coded 0/1/NA"), raw_bad)
    add_check(cohort, source_name, paste0(toupper(component), "_CUMULATIVE_DOMAIN"),
              paste0(component, " is coded 0/1/NA"), cumulative_bad)

    living <- copy(x[event_type == "visit"])
    if (identical(source_name, "exposure_panel")) {
      # The exposure panel is left-truncated by analysis eligibility.  A
      # participant can therefore enter it after a positive diagnosis that is
      # visible only in the full cleaned history.  Recompute the cumulative
      # component from that full raw history, then match the expected value at
      # each retained id-wave; recomputing only inside the restricted panel
      # would incorrectly erase pre-entry diagnoses.
      history_path <- file.path(clean_root, paste0(tolower(cohort), "_clean.rds"))
      history <- as.data.table(readRDS(history_path))
      setorder(history, id, event_time, wave, na.last = TRUE)
      history <- history[event_type == "visit", .(id, wave, raw_value = get(raw_field))]
      history[, expected_component := cum_diagnosis_independent(raw_value), by = id]
      expected_lookup <- unique(history[, .(id, wave, expected_component)], by = c("id", "wave"))
      living <- merge(living, expected_lookup, by = c("id", "wave"), all.x = TRUE, sort = FALSE)
    } else {
      living[, expected_component := cum_diagnosis_independent(get(raw_field)), by = id]
    }
    component_mismatch <- living[!same_with_na(get(component), expected_component), .N]
    reversals <- living[!is.na(get(component)), .(
      reversal = any(diff(as.integer(get(component))) < 0L)
    ), by = id][reversal == TRUE, .N]
    add_check(
      cohort, source_name, paste0(toupper(component), "_CUMULATIVE_RECOMPUTE"),
      paste0(component, " matches an independent cumulative-ever recomputation"),
      component_mismatch
    )
    add_check(
      cohort, source_name, paste0(toupper(component), "_NO_REVERSAL"),
      paste0(component, " never changes from 1 back to 0 within participant"),
      reversals
    )
    component_tables[[length(component_tables) + 1L]] <<- data.table(
      cohort = cohort,
      source = source_name,
      component = component,
      rows = nrow(living),
      participants = uniqueN(living$id),
      raw_nonmissing = sum(!is.na(living[[raw_field]])),
      cumulative_nonmissing = sum(!is.na(living[[component]])),
      positive_rows = sum(living[[component]] == 1L, na.rm = TRUE),
      recomputation_mismatches = component_mismatch,
      participants_with_reversal = reversals
    )
  }

  x[, expected_state := derive_state_independent(hypertension, heart_disease, stroke, dead)]
  state_mismatch <- x[!same_with_na(as.integer(state), expected_state), .N]
  add_check(
    cohort, source_name, "STATE_INDEPENDENT_RECOMPUTE",
    "Stored state exactly matches a fresh derivation from the three cumulative diagnoses and death",
    state_mismatch
  )

  live <- x[event_type == "visit"]
  live[, `:=`(
    known_components = rowSums(!is.na(cbind(hypertension, heart_disease, stroke))),
    positive_components = rowSums(cbind(hypertension, heart_disease, stroke) == 1L, na.rm = TRUE)
  )]
  ambiguous_classification <- live[
    (!is.na(state) & state %in% 0:2 & known_components < 3L) |
      (!is.na(state) & state == 3L & positive_components < 2L),
    .N
  ]
  add_check(
    cohort, source_name, "NO_AMBIGUOUS_CLASSIFICATION",
    "S0-S2 require all three components; S3 requires at least two observed positive diagnoses",
    ambiguous_classification
  )

  state_reversals <- live[!is.na(state), .(
    reversal = any(diff(as.integer(state)) < 0L)
  ), by = id][reversal == TRUE, .N]
  add_check(
    cohort, source_name, "NO_LIVING_STATE_REVERSAL",
    "Living disease state never moves to a lower state within participant",
    state_reversals
  )

  death_wrong_state <- x[dead == 1L & state != 4L, .N]
  death_not_final <- x[, .(
    violation = any(dead == 1L & seq_len(.N) != .N)
  ), by = id][violation == TRUE, .N]
  add_check(cohort, source_name, "DEATH_IS_S4", "Every death record is coded S4", death_wrong_state)
  add_check(cohort, source_name, "DEATH_FINAL", "S4/death is the final record for each participant", death_not_final)

  combo <- live[, .(
    rows = .N,
    participants = uniqueN(id)
  ), by = .(hypertension, heart_disease, stroke, state, expected_state, known_components, positive_components)]
  combo[, `:=`(
    cohort = cohort,
    source = source_name,
    state_label = unname(state_labels[as.character(state)])
  )]
  setcolorder(combo, c(
    "cohort", "source", "hypertension", "heart_disease", "stroke", "known_components",
    "positive_components", "state", "expected_state", "state_label", "rows", "participants"
  ))
  combination_tables[[length(combination_tables) + 1L]] <<- combo

  add_check(
    cohort, source_name, "EIGHT_COMPLETE_COMBINATIONS",
    "All eight fully observed binary disease combinations occur and map deterministically",
    as.integer(uniqueN(combo[known_components == 3L, paste(hypertension, heart_disease, stroke)]) != 8L),
    paste0("observed_combinations=", uniqueN(combo[known_components == 3L, paste(hypertension, heart_disease, stroke)]))
  )
  invisible(NULL)
}

for (cohort in c("HRS", "ELSA")) {
  stem <- tolower(cohort)
  audit_dataset(cohort, "clean_longitudinal", file.path(clean_root, paste0(stem, "_clean.rds")))
  audit_dataset(cohort, "exposure_panel", file.path(exposure_root, paste0(stem, "_exposure_panel.rds")))
}

# Wave 10 is appended from the original ELSA interviewer file because the
# installed GlobalAgingR 2.0.0 spine ends at wave 9.  Reconstruct the broad
# Harmonized ELSA heart-problem phenotype independently from the raw questions
# so a routed myocardial-infarction response cannot be mistaken for absence.
w10_path <- file.path(stage_root, "elsa", "w10.dta")
w10_required <- c(
  "idauniq", "HEEverBP", "HENvrBP", "HEEverAN", "HENvrAN",
  "HEEverHF", "HENvrHF", "HEEverHM", "HENvrHM", "HEEverAH", "HENvrAH",
  "HEEverOH", "HENvrOH", "HEEverST", "HENvrST",
  "HEAgeMI", "HEOthHA", "HEMedHA"
)
if (!file.exists(w10_path)) {
  add_check("ELSA", "w10_raw_heart_mapping", "INPUT_EXISTS", "ELSA Wave 10 raw interviewer file exists", 1L, w10_path)
} else {
  w10_raw <- haven::read_dta(w10_path, col_select = dplyr::any_of(w10_required))
  missing_w10 <- setdiff(w10_required, names(w10_raw))
  add_check(
    "ELSA", "w10_raw_heart_mapping", "REQUIRED_FIELDS",
    "Wave 10 broad heart-problem and myocardial-infarction routing fields are present",
    length(missing_w10), paste(missing_w10, collapse = "; ")
  )
  if (!length(missing_w10)) {
    w10_dt <- as.data.table(w10_raw)
    w10_dt[, expected_hypertension_raw := ever_with_feedforward_independent(HEEverBP, HENvrBP)]
    w10_components <- derive_elsa_w10_heart_components_independent(w10_dt)
    w10_component_names <- names(w10_components)
    w10_dt[, (w10_component_names) := w10_components]
    w10_dt[, expected_heart_disease_raw := combine_binary_any_independent(
      expected_heart_angina_raw, expected_heart_mi_raw, expected_heart_failure_raw,
      expected_heart_murmur_raw, expected_heart_rhythm_raw, expected_heart_other_raw
    )]
    w10_dt[, expected_stroke_raw := ever_with_feedforward_independent(HEEverST, HENvrST)]
    w10_dt[, id := suppressWarnings(as.numeric(idauniq))]
    w10_dt[, mi_route_evidence := as.integer(
      (!is.na(as.numeric(HEAgeMI)) & !(as.numeric(HEAgeMI) %in% c(-1, -2))) |
        (!is.na(as.numeric(HEOthHA)) & !(as.numeric(HEOthHA) %in% c(-1, -2))) |
        (!is.na(as.numeric(HEMedHA)) & !(as.numeric(HEMedHA) %in% c(-1, -2)))
    )]
    add_check(
      "ELSA", "w10_raw_heart_mapping", "MI_ROUTE_IMPLIES_HEART_HISTORY",
      "Every substantive HEAgeMI, HEOthHA or HEMedHA response establishes heart-attack history",
      w10_dt[mi_route_evidence == 1L & expected_heart_disease_raw != 1L, .N],
      paste0("routed_records=", w10_dt[mi_route_evidence == 1L, .N])
    )
    feedforward_pairs <- list(
      hypertension = c("HEEverBP", "HENvrBP"),
      angina = c("HEEverAN", "HENvrAN"),
      heart_failure = c("HEEverHF", "HENvrHF"),
      heart_murmur = c("HEEverHM", "HENvrHM"),
      abnormal_rhythm = c("HEEverAH", "HENvrAH"),
      other_heart = c("HEEverOH", "HENvrOH"),
      stroke = c("HEEverST", "HENvrST")
    )
    feedforward_audit <- rbindlist(lapply(names(feedforward_pairs), function(component) {
      fields <- feedforward_pairs[[component]]
      reconstructed <- ever_with_feedforward_independent(
        w10_dt[[fields[1L]]], w10_dt[[fields[2L]]]
      )
      confirmed <- suppressWarnings(as.numeric(w10_dt[[fields[2L]]])) == 1
      data.table(
        component = component,
        confirmed_prior_diagnoses = sum(confirmed, na.rm = TRUE),
        violations = sum(confirmed & reconstructed != 1L, na.rm = TRUE)
      )
    }))
    add_check(
      "ELSA", "w10_raw_heart_mapping", "FEEDFORWARD_CONFIRMATION_IMPLIES_HISTORY",
      "Every HENvrXX=1 confirmation is retained as an ever-diagnosis",
      sum(feedforward_audit$violations),
      paste0(
        feedforward_audit$component, "=", feedforward_audit$confirmed_prior_diagnoses,
        collapse = "; "
      )
    )

    clean_w10_path <- file.path(clean_root, "elsa_clean.rds")
    exposure_w10_path <- file.path(exposure_root, "elsa_exposure_panel.rds")
    for (mapping_source in c("clean_longitudinal", "exposure_panel")) {
      mapping_path <- if (mapping_source == "clean_longitudinal") clean_w10_path else exposure_w10_path
      if (!file.exists(mapping_path)) {
        add_check("ELSA", "w10_raw_heart_mapping", paste0(toupper(mapping_source), "_EXISTS"),
                  paste0(mapping_source, " exists for Wave 10 comparison"), 1L, mapping_path)
        next
      }
      mapped <- as.data.table(readRDS(mapping_path))
      required_mapped <- c(
        "id", "wave", "hypertension_raw", "heart_disease_raw", "stroke_raw",
        "heart_angina_raw", "heart_mi_raw", "heart_failure_raw",
        "heart_murmur_raw", "heart_rhythm_raw", "heart_other_raw"
      )
      mapped_missing <- setdiff(required_mapped, names(mapped))
      if (length(mapped_missing)) {
        add_check("ELSA", "w10_raw_heart_mapping", paste0(toupper(mapping_source), "_SCHEMA"),
                  "Mapped dataset contains all three Wave 10 raw disease fields", length(mapped_missing),
                  paste(mapped_missing, collapse = "; "))
        next
      }
      mapped <- mapped[wave == 10L, .(
        id = suppressWarnings(as.numeric(id)),
        stored_hypertension_raw = as.integer(hypertension_raw),
        stored_heart_disease_raw = as.integer(heart_disease_raw),
        stored_stroke_raw = as.integer(stroke_raw),
        stored_heart_angina_raw = as.integer(heart_angina_raw),
        stored_heart_mi_raw = as.integer(heart_mi_raw),
        stored_heart_failure_raw = as.integer(heart_failure_raw),
        stored_heart_murmur_raw = as.integer(heart_murmur_raw),
        stored_heart_rhythm_raw = as.integer(heart_rhythm_raw),
        stored_heart_other_raw = as.integer(heart_other_raw),
        mapping_present = 1L
      )]
      extra_mapped_ids <- setdiff(mapped$id, w10_dt$id)
      add_check(
        "ELSA", "w10_raw_heart_mapping", paste0(toupper(mapping_source), "_IDS_VALID"),
        paste0(mapping_source, " Wave 10 rows are a subset of original interviewer records"),
        length(extra_mapped_ids), paste0("mapped_records=", nrow(mapped))
      )
      compared <- merge(
        w10_dt[, .(
          id, expected_hypertension_raw, expected_heart_disease_raw,
          expected_stroke_raw,
          expected_heart_angina_raw, expected_heart_mi_raw,
          expected_heart_failure_raw, expected_heart_murmur_raw,
          expected_heart_rhythm_raw, expected_heart_other_raw
        )], mapped,
        by = "id", all.x = TRUE, sort = FALSE
      )
      # clean_longitudinal retains every Wave 10 interviewer record, whereas
      # exposure_panel is intentionally restricted to the eligible analytic
      # cohort.  Therefore compare values only for rows present in the mapped
      # dataset and audit the subset relation separately above.
      for (component in c("hypertension", "heart_disease", "stroke")) {
        expected_column <- paste0("expected_", component, "_raw")
        stored_column <- paste0("stored_", component, "_raw")
        mismatch <- compared[
          mapping_present == 1L &
            !same_with_na(get(expected_column), get(stored_column)),
          .N
        ]
        add_check(
          "ELSA", "w10_raw_heart_mapping",
          paste0(toupper(mapping_source), "_", toupper(component), "_MATCH"),
          paste0(
            mapping_source, " Wave 10 ", component,
            "_raw matches an independent raw-question reconstruction"
          ),
          mismatch,
          paste0(
            "compared_records=", compared[mapping_present == 1L, .N],
            "; raw_records_not_in_analytic_subset=", compared[is.na(mapping_present), .N]
          )
        )
      }
      for (component in c("heart_angina", "heart_mi", "heart_failure", "heart_murmur", "heart_rhythm", "heart_other")) {
        expected_column <- paste0("expected_", component, "_raw")
        stored_column <- paste0("stored_", component, "_raw")
        mismatch <- compared[
          mapping_present == 1L & !same_with_na(get(expected_column), get(stored_column)),
          .N
        ]
        add_check(
          "ELSA", "w10_raw_heart_mapping",
          paste0(toupper(mapping_source), "_", toupper(component), "_MATCH"),
          paste0(mapping_source, " Wave 10 ", component, " matches its independent component reconstruction"),
          mismatch,
          paste0("compared_records=", compared[mapping_present == 1L, .N])
        )
      }
      if (mapping_source == "clean_longitudinal") {
        setnames(
          mapped,
          c("stored_hypertension_raw", "stored_heart_disease_raw", "stored_stroke_raw"),
          c("clean_hypertension_raw", "clean_heart_disease_raw", "clean_stroke_raw")
        )
        w10_heart_audit <- merge(w10_dt, mapped, by = "id", all.x = TRUE, sort = FALSE)
        w10_heart_audit[, mapping_match_all_three := as.integer(
          same_with_na(expected_hypertension_raw, clean_hypertension_raw) &
            same_with_na(expected_heart_disease_raw, clean_heart_disease_raw) &
            same_with_na(expected_stroke_raw, clean_stroke_raw)
        )]
      }
    }

    # Independently reconstruct Waves 4-9 from the wide Harmonized ELSA G3
    # source and then apply Wave-10 component disputes.  This checks the two
    # details that an internal cumulative-ever audit cannot detect: the empty
    # packaged Wave-5 angina column and the official retrospective dispute
    # rule.
    harmonized_path <- file.path(stage_root, "elsa_harmonized", "h_elsa_g3.dta")
    clean_history_path <- file.path(clean_root, "elsa_clean.rds")
    if (!file.exists(harmonized_path)) {
      add_check(
        "ELSA", "harmonized_component_backcast", "SOURCE_EXISTS",
        "Staged Harmonized ELSA G3 source exists", 1L, harmonized_path
      )
    } else if (!file.exists(clean_history_path)) {
      add_check(
        "ELSA", "harmonized_component_backcast", "CLEAN_EXISTS",
        "Clean ELSA longitudinal dataset exists", 1L, clean_history_path
      )
    } else {
      disease_bases <- c(
        "hibpe", "hearte", "stroke", "angine", "hrtatte", "conhrtfe", "hrtmre", "hrtrhme"
      )
      source_variables <- c(
        "idauniqc",
        unlist(lapply(4:9, function(w) paste0("r", w, disease_bases)), use.names = FALSE)
      )
      source_wide <- haven::read_dta(
        harmonized_path,
        col_select = dplyr::all_of(source_variables)
      )
      source_wide <- as.data.table(source_wide)
      source_panel <- rbindlist(lapply(4:9, function(w) {
        z <- data.table(
          id = suppressWarnings(as.numeric(source_wide[["idauniqc"]])),
          wave = as.integer(w),
          source_key_present = 1L
        )
        for (variable in disease_bases) {
          z[[paste0("source_", variable)]] <- binary01_independent(
            source_wide[[paste0("r", w, variable)]]
          )
        }
        z
      }))
      clean_history <- as.data.table(readRDS(clean_history_path))
      component_fields <- c(
        "hypertension_raw", "stroke_raw", "heart_disease_raw",
        "heart_angina_raw", "heart_mi_raw", "heart_failure_raw",
        "heart_murmur_raw", "heart_rhythm_raw", "heart_other_raw"
      )
      missing_component_fields <- setdiff(component_fields, names(clean_history))
      add_check(
        "ELSA", "harmonized_component_backcast", "CLEAN_COMPONENT_SCHEMA",
        "Clean ELSA retains every component needed for an independent broad-heart reconstruction",
        length(missing_component_fields), paste(missing_component_fields, collapse = "; ")
      )
      if (!length(missing_component_fields)) {
        clean_history <- clean_history[
          event_type == "visit" & wave %in% 4:9,
          c(list(id = suppressWarnings(as.numeric(id)), wave = as.integer(wave)), mget(component_fields))
        ]
        expected <- merge(clean_history, source_panel, by = c("id", "wave"), all.x = TRUE, sort = FALSE)
        add_check(
          "ELSA", "harmonized_component_backcast", "SOURCE_KEY_COVERAGE",
          "Every packaged ELSA Wave 4-9 visit has a matching Harmonized ELSA G3 id/wave key",
          expected[is.na(source_key_present), .N], paste0("clean_rows=", nrow(expected))
        )

        named_source <- cbind(
          expected$source_angine, expected$source_hrtatte, expected$source_conhrtfe,
          expected$source_hrtmre, expected$source_hrtrhme
        )
        named_source_yes <- rowSums(named_source == 1L, na.rm = TRUE) > 0L
        expected[, source_other := fcase(
          source_hearte == 1L & !named_source_yes, 1L,
          source_hearte == 0L, 0L,
          default = NA_integer_
        )]

        backcast_specs <- list(
          hypertension = list(expected = "source_hibpe", stored = "hypertension_raw", henvr = "HENvrBP"),
          angina = list(expected = "source_angine", stored = "heart_angina_raw", henvr = "HENvrAN"),
          myocardial_infarction = list(expected = "source_hrtatte", stored = "heart_mi_raw", henvr = NA_character_),
          heart_failure = list(expected = "source_conhrtfe", stored = "heart_failure_raw", henvr = "HENvrHF"),
          heart_murmur = list(expected = "source_hrtmre", stored = "heart_murmur_raw", henvr = "HENvrHM"),
          abnormal_rhythm = list(expected = "source_hrtrhme", stored = "heart_rhythm_raw", henvr = "HENvrAH"),
          other_heart = list(expected = "source_other", stored = "heart_other_raw", henvr = "HENvrOH"),
          stroke = list(expected = "source_stroke", stored = "stroke_raw", henvr = "HENvrST")
        )
        elsa_backcast_audit <- rbindlist(lapply(names(backcast_specs), function(component) {
          spec <- backcast_specs[[component]]
          disputed_ids <- numeric()
          source_positive_before <- sum(expected[[spec$expected]] == 1L, na.rm = TRUE)
          rows_backcast <- 0L
          ambiguous_other_resolved <- 0L
          if (!is.na(spec$henvr)) {
            disputed_ids <- unique(w10_dt[
              is_explicit_dispute_independent(get(spec$henvr)),
              id
            ])
            recode <- expected$id %in% disputed_ids & expected[[spec$expected]] == 1L
            rows_backcast <- sum(recode, na.rm = TRUE)
            if (any(recode, na.rm = TRUE)) {
              set(expected, i = which(recode %in% TRUE), j = spec$expected, value = 0L)
            }
            if (component == "other_heart") {
              resolve_other <- expected$id %in% disputed_ids &
                is.na(expected[[spec$expected]]) & expected$source_hearte == 1L
              ambiguous_other_resolved <- sum(resolve_other, na.rm = TRUE)
              if (any(resolve_other, na.rm = TRUE)) {
                set(expected, i = which(resolve_other %in% TRUE), j = spec$expected, value = 0L)
              }
            }
          }
          mismatch <- !same_with_na(expected[[spec$expected]], expected[[spec$stored]])
          out <- data.table(
            component = component,
            source_field = spec$expected,
            stored_field = spec$stored,
            henvr_field = spec$henvr,
            disputed_ids = length(disputed_ids),
            source_positive_rows_before_backcast = source_positive_before,
            positive_rows_backcast_to_no = rows_backcast,
            ambiguous_other_rows_resolved_to_no = ambiguous_other_resolved,
            expected_positive_rows_after_backcast = sum(expected[[spec$expected]] == 1L, na.rm = TRUE),
            mismatches = sum(mismatch)
          )
          add_check(
            "ELSA", "harmonized_component_backcast",
            paste0(toupper(component), "_MATCH"),
            paste0(
              spec$stored,
              " matches the Harmonized source after independent Wave-10 dispute backcasting"
            ),
            out$mismatches,
            paste0("disputed_ids=", out$disputed_ids)
          )
          out
        }))

        expected[, expected_broad_heart := combine_binary_any_independent(
          source_angine, source_hrtatte, source_conhrtfe,
          source_hrtmre, source_hrtrhme, source_other
        )]
        broad_mismatch <- !same_with_na(expected$expected_broad_heart, expected$heart_disease_raw)
        add_check(
          "ELSA", "harmonized_component_backcast", "BROAD_HEART_MATCH",
          "Historical broad heart disease is the OR of five named components plus residual other heart trouble after backcasting",
          sum(broad_mismatch),
          paste0(
            "source_hearte_zero_but_named_component_positive=",
            expected[source_hearte == 0L & rowSums(named_source == 1L, na.rm = TRUE) > 0L, .N]
          )
        )
        add_check(
          "ELSA", "harmonized_component_backcast", "WAVE5_ANGINA_RESTORED",
          "Wave-5 angina in the clean data matches non-empty Harmonized ELSA source values",
          expected[wave == 5L & !same_with_na(source_angine, heart_angina_raw), .N],
          paste0("wave5_source_nonmissing=", expected[wave == 5L & !is.na(source_angine), .N])
        )
      }
    }

    for (variable in w10_required[-1L]) {
      value_labels <- attr(w10_raw[[variable]], "labels")
      w10_source_rows[[length(w10_source_rows) + 1L]] <- data.table(
        cohort = "ELSA",
        package = "ELSA Wave 10 original interviewer data",
        package_version = NA_character_,
        source_object = basename(w10_path),
        source_variable = variable,
        source_label = enc2utf8(as.character(attr(w10_raw[[variable]], "label"))),
        coding = paste0(names(value_labels), "=", unname(value_labels), collapse = "; "),
        analytic_component = fcase(
          grepl("BP$", variable), "hypertension_wave10",
          grepl("ST$", variable), "stroke_wave10",
          default = "heart_disease_wave10"
        )
      )
    }
  }
}

# Audit the formal continuous-time Markov transition graph.
catalog_path <- file.path(result_root, "transition_catalog.csv")
if (!file.exists(catalog_path)) {
  add_check("BOTH", "transition_graph", "CATALOG_EXISTS", "Transition catalog exists", 1L, catalog_path)
} else {
  catalog <- fread(catalog_path)
  required <- c("from_state", "to_state")
  missing_fields <- setdiff(required, names(catalog))
  add_check("BOTH", "transition_graph", "CATALOG_SCHEMA", "Catalog has from_state and to_state", length(missing_fields), paste(missing_fields, collapse = "; "))
  if (!length(missing_fields)) {
    catalog[, edge := paste0(as.integer(from_state), "->", as.integer(to_state))]
    edge_violations <- length(setdiff(expected_q_edges, catalog$edge)) +
      length(setdiff(catalog$edge, expected_q_edges)) +
      sum(duplicated(catalog$edge))
    add_check(
      "BOTH", "transition_graph", "EXACT_EIGHT_Q_EDGES",
      "Q contains exactly 0->1, 0->2, 0->4, 1->3, 1->4, 2->3, 2->4, and 3->4",
      edge_violations, paste(sort(catalog$edge), collapse = "; ")
    )
    add_check(
      "BOTH", "transition_graph", "NO_DIRECT_S0_S3_INTENSITY",
      "The CTMC Q matrix does not estimate a direct S0->S3 intensity",
      catalog[edge == "0->3", .N]
    )
    add_check(
      "BOTH", "transition_graph", "DEATH_ABSORBING",
      "S4/death has no outgoing transition intensity",
      catalog[from_state == 4L, .N]
    )
  }
}

panel_edge_tables <- list()
for (cohort in c("HRS", "ELSA")) {
  stem <- tolower(cohort)
  interval_path <- file.path(data_root, paste0(stem, "_multistate_intervals.rds"))
  if (!file.exists(interval_path)) {
    add_check(cohort, "analysis_intervals", "INTERVALS_EXIST", "Prepared interval dataset exists", 1L, interval_path)
    next
  }
  z <- as.data.table(readRDS(interval_path))
  z[, edge := paste0(as.integer(from_state), "->", as.integer(to_state))]
  changed <- z[as.integer(from_state) != as.integer(to_state)]
  invalid <- changed[!edge %chin% expected_panel_edges, .N]
  add_check(
    cohort, "analysis_intervals", "OBSERVED_EDGE_SET",
    "Changed panel states use the eight Q edges or interval-censored S0->S3 skip",
    invalid, paste(sort(unique(changed$edge)), collapse = "; ")
  )
  s0_s3_n <- changed[edge == "0->3", .N]
  add_check(
    cohort, "analysis_intervals", "S0_S3_PANEL_SKIP_DOCUMENTED",
    "Observed S0->S3 records are retained as panel skips but excluded from the Q edge catalog",
    0L, paste0("observed_S0_to_S3=", s0_s3_n)
  )
  edge_counts <- changed[, .(
    observed_intervals = .N,
    participants = uniqueN(id)
  ), by = .(from_state, to_state, edge)]
  edge_counts[, `:=`(
    cohort = cohort,
    q_intensity_edge = as.integer(edge %chin% expected_q_edges),
    panel_skip_only = as.integer(edge == "0->3")
  )]
  setcolorder(edge_counts, c("cohort", "from_state", "to_state", "edge", "observed_intervals", "participants", "q_intensity_edge", "panel_skip_only"))
  panel_edge_tables[[length(panel_edge_tables) + 1L]] <- edge_counts
}

# Record source labels directly from GlobalAgingR 2.0.0, making provenance
# inspectable without relying on the cleaning script's comments.
source_rows <- list()
for (cohort in c("HRS", "ELSA")) {
  object_name <- tolower(cohort)
  env <- new.env(parent = baseenv())
  utils::data(list = object_name, package = "GlobalAgingR", envir = env)
  package_data <- get(object_name, envir = env, inherits = FALSE)
  for (variable in c("hibpe", "hearte", "stroke")) {
    value_labels <- attr(package_data[[variable]], "labels")
    source_rows[[length(source_rows) + 1L]] <- data.table(
      cohort = cohort,
      package = "GlobalAgingR",
      package_version = as.character(utils::packageVersion("GlobalAgingR")),
      source_object = object_name,
      source_variable = variable,
      source_label = enc2utf8(as.character(attr(package_data[[variable]], "label"))),
      coding = paste0(names(value_labels), "=", unname(value_labels), collapse = "; "),
      analytic_component = c(hibpe = "hypertension", hearte = "heart_disease", stroke = "stroke")[[variable]]
    )
  }
}

check_table <- rbindlist(checks, use.names = TRUE, fill = TRUE)
combination_table <- rbindlist(combination_tables, use.names = TRUE, fill = TRUE)
component_table <- rbindlist(component_tables, use.names = TRUE, fill = TRUE)
panel_edge_table <- rbindlist(panel_edge_tables, use.names = TRUE, fill = TRUE)
source_table <- rbindlist(c(source_rows, w10_source_rows), use.names = TRUE, fill = TRUE)

setorder(check_table, cohort, source, check_id)
setorder(combination_table, cohort, source, state, hypertension, heart_disease, stroke)
setorder(component_table, cohort, source, component)
setorder(panel_edge_table, cohort, from_state, to_state)

fwrite(check_table, file.path(qc_root, "disease_state_verification.csv"), na = "")
fwrite(combination_table, file.path(qc_root, "disease_state_combination_counts.csv"), na = "")
fwrite(component_table, file.path(qc_root, "disease_component_cumulative_audit.csv"), na = "")
fwrite(panel_edge_table, file.path(qc_root, "disease_state_observed_edge_audit.csv"), na = "")
fwrite(source_table, file.path(qc_root, "disease_component_source_provenance.csv"), na = "")
if (nrow(w10_heart_audit)) {
  fwrite(w10_heart_audit, file.path(qc_root, "elsa_w10_heart_mapping_audit.csv"), na = "")
}
if (nrow(elsa_backcast_audit)) {
  fwrite(
    elsa_backcast_audit,
    file.path(qc_root, "elsa_harmonized_component_backcast_verification.csv"),
    na = ""
  )
}

w10_summary_lines <- if (exists("w10_dt", inherits = FALSE) && nrow(w10_dt)) {
  c(
    "",
    "ELSA Wave 10 raw disease reconstruction:",
    paste0(
      "- Hypertension: positive=", w10_dt[expected_hypertension_raw == 1L, .N],
      "; negative=", w10_dt[expected_hypertension_raw == 0L, .N],
      "; unknown=", w10_dt[is.na(expected_hypertension_raw), .N], "."
    ),
    paste0(
      "- Broad heart problem: positive=", w10_dt[expected_heart_disease_raw == 1L, .N],
      "; negative=", w10_dt[expected_heart_disease_raw == 0L, .N],
      "; unknown=", w10_dt[is.na(expected_heart_disease_raw), .N], "."
    ),
    paste0(
      "- Stroke: positive=", w10_dt[expected_stroke_raw == 1L, .N],
      "; negative=", w10_dt[expected_stroke_raw == 0L, .N],
      "; unknown=", w10_dt[is.na(expected_stroke_raw), .N], "."
    ),
    paste0("- Myocardial-infarction route evidence: ", w10_dt[mi_route_evidence == 1L, .N], " records."),
    "- Direct heart components: angina, heart failure, heart murmur, abnormal heart rhythm, and other heart trouble.",
    "- HENvrXX=1 confirmations retain diagnoses carried forward from prior waves.",
    "- Labelled HENvrXX=2/3/4 disputes clear the matching prior component-positive reports; undocumented codes 85/86 remain unknown.",
    "- The archived Wave-10 file has no HEEverMI/HENvrMI fields, so a Wave-10 historical-MI dispute cannot be backcast."
  )
} else {
  c("", "ELSA Wave 10 raw disease reconstruction was unavailable.")
}

summary_lines <- c(
  "# Independent HRS/ELSA disease-state verification",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("GlobalAgingR version: ", as.character(utils::packageVersion("GlobalAgingR"))),
  "",
  "State definitions:",
  paste0("- ", names(state_labels), ": ", unname(state_labels)),
  "",
  paste0("Formal Q edges: ", paste(expected_q_edges, collapse = ", ")),
  "Observed S0->S3 jumps are interval-censored skips over unobserved intermediate disease onset and are not a ninth Q intensity.",
  w10_summary_lines,
  "",
  paste0("Checks passed: ", sum(check_table$status == "PASS"), "/", nrow(check_table)),
  paste0("Checks failed: ", sum(check_table$status == "FAIL"))
)
writeLines(summary_lines, file.path(qc_root, "disease_state_verification.md"), useBytes = TRUE)

print(check_table)
if (any(check_table$status == "FAIL")) {
  stop(
    "Disease-state verification failed: ",
    paste(check_table[status == "FAIL", paste(cohort, source, check_id, sep = "/")], collapse = "; ")
  )
}

cat("Independent disease-state verification passed.\n")
cat("QC directory: ", normalizePath(qc_root, winslash = "/", mustWork = TRUE), "\n", sep = "")
