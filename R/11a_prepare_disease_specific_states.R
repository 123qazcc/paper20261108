#!/usr/bin/env Rscript

# Disease-specific extension: data preparation and audit only.
#
# This script deliberately does not modify or refit the frozen S0-S4 analysis.
# It reads the final delivered HRS/ELSA panels, expands the three cumulative
# disease components into the monotone D0-D7 combination lattice, appends D8
# death for HRS, and writes model-ready panels plus reproducibility/QC artefacts.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

options(stringsAsFactors = FALSE, width = 160)

# This Windows host can start R with an unavailable C.UTF-8 locale.  Switch the
# character type locale only when a non-ASCII working path would otherwise be
# unreadable; LC_NUMERIC remains unchanged.
if (.Platform$OS.type == "windows" &&
    identical(Sys.getlocale("LC_CTYPE"), "C") &&
    grepl("[^ -~]", getwd())) {
  suppressWarnings(try(Sys.setlocale("LC_CTYPE", "Chinese"), silent = TRUE))
}

script_version <- "2026-09-09-v1"
run_timestamp_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

script_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", script_args, value = TRUE)
project_root_env <- Sys.getenv("GA_PROJECT_ROOT", unset = "")
if (nzchar(project_root_env)) {
  # On Windows, dirname() can fail when R starts under a broken C.UTF-8 locale
  # and the absolute project path contains non-ASCII characters.  An explicit
  # project-root override avoids deriving the root from --file in that case.
  project_root <- project_root_env
  if (!dir.exists(project_root)) {
    stop("GA_PROJECT_ROOT does not exist: ", project_root, call. = FALSE)
  }
  script_file <- file.path(
    project_root, "R", "11a_prepare_disease_specific_states.R"
  )
} else if (file.exists(file.path("R", "11a_prepare_disease_specific_states.R"))) {
  # The normal invocation is from the project root.  Keeping this path relative
  # also avoids unnecessary native-encoding round trips on Windows.
  project_root <- "."
  script_file <- file.path("R", "11a_prepare_disease_specific_states.R")
} else {
  script_file <- if (length(script_arg)) {
    sub("^--file=", "", script_arg[1L])
  } else {
    file.path(getwd(), "R", "11a_prepare_disease_specific_states.R")
  }
  script_file <- normalizePath(script_file, winslash = "/", mustWork = TRUE)
  project_root <- dirname(dirname(script_file))
}

final_delivery_root <- file.path(
  project_root, "outputs", "hrs_elsa_multistate_final_20260909"
)
input_root <- file.path(final_delivery_root, "data_multistate")
nextsteps_root <- Sys.getenv(
  "GA_NEXTSTEPS_OUTPUT_ROOT",
  unset = file.path(project_root, "outputs", "hrs_elsa_multistate_nextsteps_20260909")
)
data_output_root <- file.path(nextsteps_root, "data_disease_specific")
qc_output_root <- file.path(nextsteps_root, "qc", "disease_specific")

dir.create(data_output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_output_root, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  HRS = file.path(input_root, "hrs_msm_panel.rds"),
  ELSA = file.path(input_root, "elsa_msm_panel.rds")
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop(
    "Authoritative final panel input(s) missing: ",
    paste(unname(missing_inputs), collapse = "; "),
    call. = FALSE
  )
}

sha256_file <- function(path) {
  if (!file.exists(path)) stop("Cannot fingerprint missing file: ", path, call. = FALSE)
  digest::digest(file = path, algo = "sha256")
}

atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  suffix <- if (grepl("\\.gz$", path, ignore.case = TRUE)) ".csv.gz" else ".csv"
  temporary <- tempfile(pattern = ".tmp_", tmpdir = dirname(path), fileext = suffix)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  data.table::fwrite(x, temporary, na = "")
  if (!file.copy(temporary, path, overwrite = TRUE)) {
    stop("Atomic CSV publication failed: ", path, call. = FALSE)
  }
  invisible(path)
}

atomic_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".tmp_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(x, temporary, compress = "xz")
  if (!file.copy(temporary, path, overwrite = TRUE)) {
    stop("Atomic RDS publication failed: ", path, call. = FALSE)
  }
  invisible(path)
}

state_dictionary <- data.table(
  state = 0:8,
  state_code = paste0("D", 0:8),
  hypertension = c(0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, NA_integer_),
  heart_disease = c(0L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, NA_integer_),
  stroke = c(0L, 0L, 0L, 1L, 0L, 1L, 1L, 1L, NA_integer_),
  death = c(rep(0L, 8L), 1L),
  state_label = c(
    "No target disease",
    "Hypertension only",
    "Heart disease only",
    "Stroke only",
    "Hypertension and heart disease",
    "Hypertension and stroke",
    "Heart disease and stroke without hypertension",
    "Hypertension, heart disease and stroke",
    "Death"
  ),
  hrs_primary = 1L,
  elsa_w4_w10_nonfatal_primary = c(rep(1L, 8L), 0L)
)
state_dictionary[, component_key := fifelse(
  death == 1L,
  NA_character_,
  paste0(hypertension, heart_disease, stroke)
)]

nonfatal_edges <- data.table(
  from_state = c(0L, 0L, 0L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 5L, 6L),
  to_state = c(1L, 2L, 3L, 4L, 5L, 4L, 6L, 5L, 6L, 7L, 7L, 7L)
)
target_keys <- c("0->2", "0->3", "1->4", "1->5")

make_allowed_edges <- function(cohort_name) {
  x <- copy(nonfatal_edges)
  if (cohort_name == "HRS") {
    x <- rbind(
      x,
      data.table(from_state = 0:7, to_state = 8L),
      use.names = TRUE
    )
  }
  x[, `:=`(
    cohort = cohort_name,
    analysis_id = if (cohort_name == "HRS") {
      "hrs_9state_full"
    } else {
      "elsa_8state_w4_w10_nonfatal"
    },
    edge_key = paste0(from_state, "->", to_state),
    transition = paste0("D", from_state, "->D", to_state),
    target_edge = as.integer(paste0(from_state, "->", to_state) %chin% target_keys),
    death_edge = as.integer(to_state == 8L),
    edge_type = fifelse(
      to_state == 8L,
      "death",
      fifelse(
        paste0(from_state, "->", to_state) %chin% target_keys,
        "confirmatory_disease_specific",
        "non_target_process"
      )
    )
  )]
  x <- merge(
    x,
    state_dictionary[, .(from_state = state, from_label = state_label)],
    by = "from_state", all.x = TRUE, sort = FALSE
  )
  x <- merge(
    x,
    state_dictionary[, .(to_state = state, to_label = state_label)],
    by = "to_state", all.x = TRUE, sort = FALSE
  )
  setorder(x, from_state, to_state)
  x[, edge_id := sprintf("DS_%s_%02d", cohort_name, seq_len(.N))]
  setcolorder(x, c(
    "edge_id", "cohort", "analysis_id", "from_state", "to_state",
    "transition", "from_label", "to_label", "edge_type",
    "target_edge", "death_edge", "edge_key"
  ))
  x[]
}

allowed_edges <- rbindlist(
  lapply(c("HRS", "ELSA"), make_allowed_edges),
  use.names = TRUE
)

derive_combination_state <- function(hypertension, heart_disease, stroke, dead, include_death) {
  h <- suppressWarnings(as.integer(hypertension))
  c <- suppressWarnings(as.integer(heart_disease))
  s <- suppressWarnings(as.integer(stroke))
  d <- suppressWarnings(as.integer(dead))
  out <- rep(NA_integer_, length(h))
  living_known <- d != 1L & !is.na(h) & !is.na(c) & !is.na(s)
  keys <- paste0(h[living_known], c[living_known], s[living_known])
  mapping <- c(
    `000` = 0L, `100` = 1L, `010` = 2L, `001` = 3L,
    `110` = 4L, `101` = 5L, `011` = 6L, `111` = 7L
  )
  out[living_known] <- unname(mapping[keys])
  if (include_death) out[d == 1L] <- 8L
  out
}

state_label_from_code <- function(x) {
  out <- rep(NA_character_, length(x))
  known <- !is.na(x)
  out[known] <- state_dictionary$state_label[match(x[known], state_dictionary$state)]
  out
}

valid_component_code_count <- function(x, dead) {
  sum(!is.na(x) & !(as.integer(x) %in% 0:1) & as.integer(dead) != 1L)
}

required_common <- c(
  "id", "cohort", "wave", "event_type", "event_time", "exposure_entry_time",
  "dead", "hypertension", "heart_disease", "stroke",
  "insomnia_score_primary", "age_bl", "female_bl", "race_bl",
  "education_level_bl", "wealth_quintile_bl", "married_bl",
  "current_smoker_bl", "current_drinker_bl", "inactive_bl", "bmi_bl",
  "diabetes_bl", "lung_disease_bl", "frequent_pain_bl"
)

m4_covariates <- c(
  "age_bl", "female_bl", "race_bl", "education_level_bl",
  "wealth_quintile_bl", "married_bl", "current_smoker_bl",
  "current_drinker_bl", "inactive_bl", "bmi_bl", "diabetes_bl",
  "lung_disease_bl", "frequent_pain_bl"
)

exposure_columns <- c(
  "insomnia_score_primary", "insomnia_burden_code", "insomnia_high_primary",
  "insomnia_onset_primary", "insomnia_maintenance_primary",
  "insomnia_nonrestorative_primary", "insomnia_score_cumavg",
  "insomnia_score_prev_valid", "insomnia_score_lag1_row",
  "insomnia_dynamic_pattern", "insomnia_dynamic_code",
  "insomnia_assessment_complete", "insomnia_assessment_number",
  "insomnia_score_extended4", "insomnia_early_waking_current",
  "sleep_medication_current", "insomnia_score_bl", "insomnia_burden_bl"
)

trace_columns <- c(
  "id", "cohort", "wave", "event_type", "event_time", "event_time_original",
  "radyear", "interview_year",
  "exposure_entry_time", "analysis_sample", "analysis_baseline", "dead",
  "death_year", "death_source", "death_time_midyear", "death_time_adjusted",
  "death_time_method_analysis", "last_alive_time", "death_time_from_span",
  "state", "state_msm", "time_msm", "entry_time", "time_since_entry",
  "analysis_time_since_entry", "entry_wave", "entry_year", "weight_bl"
)

component_columns <- c(
  "hibpe", "hearte", "rxhibp", "source_hibpe", "source_hearte",
  "source_stroke", "heart_aggregate_source_raw",
  "hypertension_raw", "heart_disease_raw", "stroke_raw",
  "hypertension", "heart_disease", "stroke",
  "heart_disease_lower", "heart_disease_upper",
  "heart_angina_raw", "heart_mi_raw", "heart_failure_raw",
  "heart_murmur_raw", "heart_rhythm_raw", "heart_other_raw",
  "heart_other_ambiguous_raw", "heart_other_lower_raw",
  "heart_other_upper_raw", "heart_disease_raw_lower",
  "heart_disease_raw_upper", "state_heart_lower", "state_heart_upper"
)

input_fingerprint <- function(cohort_name, path, phase, n_rows = NA_integer_, n_ids = NA_integer_) {
  info <- file.info(path)
  data.table(
    cohort = cohort_name,
    phase = phase,
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(info$size),
    modified_utc = format(info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    sha256 = sha256_file(path),
    n_rows = as.integer(n_rows),
    n_participants = as.integer(n_ids)
  )
}

pre_fingerprints <- rbindlist(lapply(names(input_files), function(cohort_name) {
  input_fingerprint(cohort_name, input_files[[cohort_name]], "before_read")
}))

prepare_cohort <- function(cohort_name, input_path) {
  message("Reading authoritative ", cohort_name, " panel: ", input_path)
  source <- as.data.table(readRDS(input_path))
  source_rows <- nrow(source)
  source_participants <- uniqueN(source$id)

  missing_columns <- setdiff(required_common, names(source))
  if (length(missing_columns)) {
    stop(
      cohort_name, " panel is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  source[, source_row_number := .I]
  source[, `:=`(
    id = as.numeric(id),
    wave = as.integer(wave),
    event_time = as.numeric(event_time),
    exposure_entry_time = as.numeric(exposure_entry_time),
    dead = as.integer(dead),
    hypertension = as.integer(hypertension),
    heart_disease = as.integer(heart_disease),
    stroke = as.integer(stroke)
  )]

  source_death_audit <- source[dead == 1L, .(
    source_death_rows = .N,
    source_death_participants = uniqueN(id),
    minimum_death_year = suppressWarnings(min(as.numeric(death_year), na.rm = TRUE)),
    maximum_death_year = suppressWarnings(max(as.numeric(death_year), na.rm = TRUE))
  )]
  for (column in c("minimum_death_year", "maximum_death_year")) {
    if (!is.finite(source_death_audit[[column]])) source_death_audit[[column]] <- NA_real_
  }

  selector <- !is.na(source$id) & is.finite(source$event_time) &
    is.finite(source$exposure_entry_time) & source$event_time >= source$exposure_entry_time
  if ("analysis_sample" %in% names(source)) {
    selector <- selector & as.integer(source$analysis_sample) == 1L
  }
  if (cohort_name == "HRS") {
    selector <- selector & source$event_type %chin% c("visit", "death")
  } else {
    selector <- selector & source$event_type == "visit" &
      source$wave %in% 4:10 & source$dead != 1L
  }
  panel <- copy(source[selector])
  setorder(panel, id, event_time, dead, wave, na.last = TRUE)
  panel[, is_first_analysis_record := as.integer(seq_len(.N) == 1L), by = id]

  panel[, disease_specific_state := derive_combination_state(
    hypertension, heart_disease, stroke, dead,
    include_death = cohort_name == "HRS"
  )]
  panel[, `:=`(
    disease_specific_state_code = fifelse(
      is.na(disease_specific_state), NA_character_, paste0("D", disease_specific_state)
    ),
    disease_specific_state_label = state_label_from_code(disease_specific_state),
    disease_specific_state_msm = disease_specific_state + 1L,
    disease_specific_time = event_time - exposure_entry_time,
    disease_specific_state_known = as.integer(!is.na(disease_specific_state)),
    m4_complete_case = as.integer(complete.cases(.SD))
  ), .SDcols = m4_covariates]

  if (cohort_name == "ELSA" &&
      all(c("heart_disease_lower", "heart_disease_upper") %in% names(panel))) {
    panel[, `:=`(
      disease_specific_state_lower = derive_combination_state(
        hypertension, heart_disease_lower, stroke, dead, include_death = FALSE
      ),
      disease_specific_state_upper = derive_combination_state(
        hypertension, heart_disease_upper, stroke, dead, include_death = FALSE
      )
    )]
  } else {
    panel[, `:=`(
      disease_specific_state_lower = disease_specific_state,
      disease_specific_state_upper = disease_specific_state
    )]
  }

  panel[, valid_state_records := sum(!is.na(disease_specific_state)), by = id]
  panel[, disease_specific_subject_eligible := as.integer(valid_state_records >= 2L)]

  # Retain a deliberately explicit analysis schema rather than serialising every
  # source field.  ELSA component-level heart variables are kept when available.
  cohort_exposure_columns <- grep("^insomnia_", names(panel), value = TRUE)
  keep_columns <- unique(c(
    trace_columns, component_columns, exposure_columns,
    cohort_exposure_columns, m4_covariates,
    "source_row_number", "disease_specific_state",
    "disease_specific_state_code", "disease_specific_state_label",
    "disease_specific_state_msm", "disease_specific_time",
    "disease_specific_state_known", "disease_specific_state_lower",
    "disease_specific_state_upper", "valid_state_records",
    "disease_specific_subject_eligible", "m4_complete_case",
    "is_first_analysis_record"
  ))
  keep_columns <- intersect(keep_columns, names(panel))
  panel <- panel[, ..keep_columns]
  setorder(panel, id, event_time, dead, wave, na.last = TRUE)

  list(
    cohort = cohort_name,
    source_rows = source_rows,
    source_participants = source_participants,
    source_exposure_columns = grep("^insomnia_", names(source), value = TRUE),
    source_death_audit = source_death_audit,
    panel = panel,
    invalid_component_codes = data.table(
      cohort = cohort_name,
      component = c("hypertension", "heart_disease", "stroke"),
      invalid_rows = c(
        valid_component_code_count(source$hypertension, source$dead),
        valid_component_code_count(source$heart_disease, source$dead),
        valid_component_code_count(source$stroke, source$dead)
      )
    )
  )
}

prepared <- lapply(names(input_files), function(cohort_name) {
  prepare_cohort(cohort_name, input_files[[cohort_name]])
})
names(prepared) <- names(input_files)

post_fingerprints <- rbindlist(lapply(names(input_files), function(cohort_name) {
  item <- prepared[[cohort_name]]
  input_fingerprint(
    cohort_name, input_files[[cohort_name]], "after_read",
    item$source_rows, item$source_participants
  )
}))
input_fingerprints <- rbind(pre_fingerprints, post_fingerprints, use.names = TRUE)

make_intervals <- function(panel, cohort_name, edges) {
  x <- copy(panel)
  setorder(x, id, event_time, dead, wave, na.last = TRUE)
  next_columns <- c(
    "disease_specific_state", "disease_specific_state_code",
    "hypertension", "heart_disease", "stroke", "event_time", "wave",
    "event_type", "source_row_number"
  )
  next_names <- paste0("next_", next_columns)
  x[, (next_names) := shift(.SD, type = "lead"), by = id, .SDcols = next_columns]
  pairs <- x[!is.na(next_event_time)]
  pairs[, `:=`(
    interval_years = next_event_time - event_time,
    from_state = disease_specific_state,
    to_state = next_disease_specific_state,
    from_state_code = disease_specific_state_code,
    to_state_code = next_disease_specific_state_code
  )]
  pairs[, transition := fifelse(
    !is.na(from_state) & !is.na(to_state),
    paste0("D", from_state, "->D", to_state),
    NA_character_
  )]
  pairs[, `:=`(
    state_changed = as.integer(!is.na(from_state) & !is.na(to_state) & from_state != to_state),
    destination_is_death = as.integer(to_state == 8L),
    n_components_added = fifelse(
      !is.na(from_state) & !is.na(to_state) & to_state != 8L,
      (next_hypertension - hypertension) +
        (next_heart_disease - heart_disease) +
        (next_stroke - stroke),
      NA_integer_
    ),
    components_monotone = as.integer(
      to_state == 8L |
        (!is.na(hypertension) & !is.na(next_hypertension) &
           !is.na(heart_disease) & !is.na(next_heart_disease) &
           !is.na(stroke) & !is.na(next_stroke) &
           next_hypertension >= hypertension &
           next_heart_disease >= heart_disease &
           next_stroke >= stroke)
    )
  )]
  edge_keys <- edges$edge_key
  pairs[, `:=`(
    direct_q_edge = as.integer(
      !is.na(from_state) & !is.na(to_state) &
        paste0(from_state, "->", to_state) %chin% edge_keys
    ),
    target_direct_edge = as.integer(
      !is.na(from_state) & !is.na(to_state) &
        paste0(from_state, "->", to_state) %chin% target_keys
    ),
    interval_censored_multicomponent_jump = as.integer(
      !is.na(n_components_added) & n_components_added >= 2L
    ),
    transition_state_known = as.integer(!is.na(from_state) & !is.na(to_state)),
    admissible_monotone_observation = as.integer(
      !is.na(from_state) & !is.na(to_state) &
        (from_state == to_state | to_state == 8L | components_monotone == 1L)
    )
  )]
  pairs[]
}

intervals <- lapply(names(prepared), function(cohort_name) {
  make_intervals(
    prepared[[cohort_name]]$panel,
    cohort_name,
    allowed_edges[cohort == cohort_name]
  )
})
names(intervals) <- names(prepared)

wave_state_counts <- rbindlist(lapply(names(prepared), function(cohort_name) {
  x <- copy(prepared[[cohort_name]]$panel)
  x[, reporting_state := fifelse(
    is.na(disease_specific_state_code), "D?", disease_specific_state_code
  )]
  out <- x[, .(
    rows = .N,
    participants = uniqueN(id),
    model_eligible_rows = sum(disease_specific_subject_eligible == 1L &
      disease_specific_state_known == 1L),
    m4_complete_case_rows = sum(m4_complete_case == 1L)
  ), by = .(wave, event_type, reporting_state,
            disease_specific_state, disease_specific_state_label)]
  out[, cohort := cohort_name]
  setcolorder(out, "cohort")
  out[]
}), use.names = TRUE, fill = TRUE)
setorder(wave_state_counts, cohort, wave, disease_specific_state)

combination_counts <- rbindlist(lapply(names(prepared), function(cohort_name) {
  x <- prepared[[cohort_name]]$panel
  out <- x[, .(
    rows = .N,
    participants = uniqueN(id),
    first_record_participants = sum(is_first_analysis_record == 1L)
  ), by = .(
    hypertension, heart_disease, stroke, dead,
    disease_specific_state, disease_specific_state_code,
    disease_specific_state_label
  )]
  out[, cohort := cohort_name]
  setcolorder(out, "cohort")
  out[]
}), use.names = TRUE, fill = TRUE)
setorder(combination_counts, cohort, disease_specific_state)

observed_transitions <- rbindlist(lapply(names(intervals), function(cohort_name) {
  x <- intervals[[cohort_name]][transition_state_known == 1L]
  out <- x[, .(
    intervals = .N,
    participants = uniqueN(id),
    person_years = sum(interval_years, na.rm = TRUE)
  ), by = .(
    from_state, to_state, transition,
    state_changed, direct_q_edge, target_direct_edge,
    destination_is_death, n_components_added,
    interval_censored_multicomponent_jump,
    components_monotone, admissible_monotone_observation
  )]
  out[, cohort := cohort_name]
  setcolorder(out, "cohort")
  out[]
}), use.names = TRUE, fill = TRUE)
setorder(observed_transitions, cohort, from_state, to_state)

target_edge_catalog <- data.table(
  target_id = sprintf("DST%02d", 1:4),
  from_state = c(0L, 0L, 1L, 1L),
  to_state = c(2L, 3L, 4L, 5L),
  target_disease = c("heart_disease", "stroke", "heart_disease", "stroke"),
  target_label = c(
    "No disease to heart disease only",
    "No disease to stroke only",
    "Hypertension only to hypertension and heart disease",
    "Hypertension only to hypertension and stroke"
  )
)

target_class_selector <- function(x, from_state, to_state) {
  if (from_state == 0L && to_state == 2L) return(x$to_state %in% c(2L, 4L))
  if (from_state == 0L && to_state == 3L) return(x$to_state %in% c(3L, 5L))
  x$to_state == to_state
}

joint_sequence_unknown_selector <- function(x, from_state) {
  if (from_state == 0L) return(x$to_state %in% c(6L, 7L))
  if (from_state == 1L) return(x$to_state == 7L)
  rep(FALSE, nrow(x))
}

target_edge_support <- rbindlist(lapply(names(intervals), function(cohort_name) {
  x <- intervals[[cohort_name]][
    transition_state_known == 1L & interval_years > 0 &
      admissible_monotone_observation == 1L
  ]
  rbindlist(lapply(c("MI20_eligible_full_panel", "M4_complete_case"), function(sample_name) {
    sample_x <- if (sample_name == "M4_complete_case") {
      x[m4_complete_case == 1L]
    } else x
    rbindlist(lapply(seq_len(nrow(target_edge_catalog)), function(i) {
      target <- target_edge_catalog[i]
      risk <- sample_x[from_state == target$from_state]
      direct_event <- risk$to_state == target$to_state
      target_class_event <- target_class_selector(
        risk, target$from_state, target$to_state
      )
      joint_unknown <- joint_sequence_unknown_selector(risk, target$from_state)
      data.table(
        cohort = cohort_name,
        sample = sample_name,
        target_id = target$target_id,
        target_label = target$target_label,
        target_disease = target$target_disease,
        from_state = target$from_state,
        to_state = target$to_state,
        transition = paste0("D", target$from_state, "->D", target$to_state),
        risk_intervals = nrow(risk),
        risk_participants = uniqueN(risk$id),
        person_years = sum(risk$interval_years, na.rm = TRUE),
        direct_state_events = sum(direct_event, na.rm = TRUE),
        direct_state_event_participants = uniqueN(risk$id[direct_event]),
        target_class_events_including_concurrent_hypertension = sum(
          target_class_event, na.rm = TRUE
        ),
        target_class_event_participants = uniqueN(risk$id[target_class_event]),
        joint_heart_stroke_sequence_unknown_events = sum(joint_unknown, na.rm = TRUE),
        joint_heart_stroke_sequence_unknown_participants = uniqueN(
          risk$id[joint_unknown]
        )
      )
    }))
  }))
}), use.names = TRUE, fill = TRUE)
setorder(target_edge_support, cohort, sample, target_id)

multijump_audit <- rbindlist(lapply(names(intervals), function(cohort_name) {
  x <- intervals[[cohort_name]][
    transition_state_known == 1L &
      interval_censored_multicomponent_jump == 1L
  ]
  if (!nrow(x)) {
    return(data.table(
      cohort = character(), from_state = integer(), to_state = integer(),
      transition = character(), destination_wave = integer(),
      n_components_added = integer(), intervals = integer(),
      participants = integer(), person_years = numeric()
    ))
  }
  out <- x[, .(
    intervals = .N,
    participants = uniqueN(id),
    person_years = sum(interval_years, na.rm = TRUE)
  ), by = .(
    from_state, to_state, transition,
    destination_wave = next_wave, n_components_added
  )]
  out[, cohort := cohort_name]
  setcolorder(out, "cohort")
  out[]
}), use.names = TRUE, fill = TRUE)
setorder(multijump_audit, cohort, from_state, to_state, destination_wave)

unknown_rows <- rbindlist(lapply(names(prepared), function(cohort_name) {
  x <- prepared[[cohort_name]]$panel[is.na(disease_specific_state)]
  wanted <- intersect(c(
    "id", "cohort", "wave", "event_type", "event_time", "dead",
    "hypertension", "heart_disease", "stroke", "heart_disease_lower",
    "heart_disease_upper", "disease_specific_state_lower",
    "disease_specific_state_upper", "source_row_number"
  ), names(x))
  out <- x[, ..wanted]
  out[, unknown_reason := fifelse(
    dead == 1L,
    "death row not mapped",
    paste0(
      "missing:",
      paste0(
        c("hypertension", "heart_disease", "stroke")[
          c(any(is.na(x$hypertension)), any(is.na(x$heart_disease)), any(is.na(x$stroke)))
        ],
        collapse = "|"
      )
    )
  )]
  out
}), use.names = TRUE, fill = TRUE)

make_reverse_rows <- function(panel, cohort_name) {
  x <- copy(panel[event_type != "death"])
  setorder(x, id, event_time, wave)
  out <- rbindlist(lapply(c("hypertension", "heart_disease", "stroke"), function(component) {
    previous <- x[, shift(get(component)), by = id]$V1
    bad <- !is.na(previous) & !is.na(x[[component]]) & previous == 1L & x[[component]] == 0L
    if (!any(bad)) {
      return(data.table(
        cohort = character(), id = numeric(), wave = integer(), event_time = numeric(),
        component = character(), previous_value = integer(), current_value = integer()
      ))
    }
    data.table(
      cohort = cohort_name,
      id = x$id[bad],
      wave = x$wave[bad],
      event_time = x$event_time[bad],
      component = component,
      previous_value = previous[bad],
      current_value = x[[component]][bad]
    )
  }), use.names = TRUE, fill = TRUE)
  out
}

reverse_rows <- rbindlist(lapply(names(prepared), function(cohort_name) {
  make_reverse_rows(prepared[[cohort_name]]$panel, cohort_name)
}), use.names = TRUE, fill = TRUE)

unknown_reverse_summary <- rbindlist(lapply(names(prepared), function(cohort_name) {
  panel <- prepared[[cohort_name]]$panel
  ints <- intervals[[cohort_name]]
  source_death <- prepared[[cohort_name]]$source_death_audit
  invalid <- prepared[[cohort_name]]$invalid_component_codes
  unknown_by_id <- panel[, {
    known_positions <- which(!is.na(disease_specific_state))
    internal_unknown <- if (length(known_positions) >= 2L) {
      sum(is.na(disease_specific_state) &
            seq_len(.N) > min(known_positions) & seq_len(.N) < max(known_positions))
    } else 0L
    list(internal_unknown_rows = internal_unknown)
  }, by = id]
  data.table(
    cohort = cohort_name,
    analysis_rows = nrow(panel),
    analysis_participants = uniqueN(panel$id),
    unknown_state_rows = panel[is.na(disease_specific_state), .N],
    unknown_state_participants = panel[is.na(disease_specific_state), uniqueN(id)],
    internal_unknown_rows = sum(unknown_by_id$internal_unknown_rows),
    invalid_component_code_rows = sum(invalid$invalid_rows),
    component_reversal_rows = reverse_rows[cohort == cohort_name, .N],
    nonpositive_interval_rows = ints[!is.finite(interval_years) | interval_years <= 0, .N],
    nonmonotone_known_transition_rows = ints[
      transition_state_known == 1L & admissible_monotone_observation != 1L, .N
    ],
    unknown_transition_endpoint_rows = ints[transition_state_known != 1L, .N],
    source_death_rows = source_death$source_death_rows,
    source_death_participants = source_death$source_death_participants,
    minimum_source_death_year = source_death$minimum_death_year,
    maximum_source_death_year = source_death$maximum_death_year
  )
}), use.names = TRUE, fill = TRUE)

qc_rows <- list()
add_check <- function(cohort, check_id, description, condition, observed, expected, critical = TRUE) {
  qc_rows[[length(qc_rows) + 1L]] <<- data.table(
    cohort = cohort,
    check_id = check_id,
    description = description,
    status = if (isTRUE(condition)) "PASS" else "FAIL",
    critical = as.integer(critical),
    observed = as.character(observed),
    expected = as.character(expected)
  )
}

for (cohort_name in names(prepared)) {
  panel <- prepared[[cohort_name]]$panel
  ints <- intervals[[cohort_name]]
  edges <- allowed_edges[cohort == cohort_name]
  invalid <- prepared[[cohort_name]]$invalid_component_codes

  duplicate_keys <- panel[, .N, by = .(id, event_time, wave, event_type)][N > 1L, .N]
  bad_order <- panel[, any(diff(event_time) <= 0), by = id][V1 == TRUE, .N]
  mapping_mismatch <- panel[
    !is.na(disease_specific_state) & disease_specific_state !=
      derive_combination_state(
        hypertension, heart_disease, stroke, dead,
        include_death = cohort_name == "HRS"
      ),
    .N
  ]
  unexplained_unknown <- panel[
    is.na(disease_specific_state) & dead != 1L &
      !is.na(hypertension) & !is.na(heart_disease) & !is.na(stroke),
    .N
  ]
  death_not_final <- if (cohort_name == "HRS") {
    panel[, any(dead == 1L & seq_len(.N) != .N), by = id][V1 == TRUE, .N]
  } else 0L
  target_support <- target_edge_support[
    cohort == cohort_name & sample == "MI20_eligible_full_panel"
  ]
  required_output_columns <- unique(c(
    required_common, "event_time_original", "hypertension_raw",
    "heart_disease_raw", "stroke_raw", "disease_specific_state",
    "disease_specific_state_code", "disease_specific_state_msm",
    "disease_specific_time", "disease_specific_state_known",
    "disease_specific_state_lower", "disease_specific_state_upper",
    "m4_complete_case"
  ))
  missing_output_columns <- setdiff(required_output_columns, names(panel))
  missing_exposure_columns <- setdiff(
    prepared[[cohort_name]]$source_exposure_columns, names(panel)
  )

  add_check(
    cohort_name, "DS_INPUT_DIRECTORY",
    "Input is the frozen final-delivery data_multistate directory",
    startsWith(normalizePath(input_files[[cohort_name]], winslash = "/"),
               normalizePath(input_root, winslash = "/")),
    dirname(normalizePath(input_files[[cohort_name]], winslash = "/")),
    normalizePath(input_root, winslash = "/")
  )
  add_check(cohort_name, "DS_COMPONENT_CODES", "Disease components use only 0, 1 or missing",
            sum(invalid$invalid_rows) == 0L, sum(invalid$invalid_rows), 0L)
  add_check(cohort_name, "DS_REQUIRED_OUTPUT_SCHEMA",
            "Event time, raw/cumulative components, primary exposure and M4 covariates are retained",
            length(missing_output_columns) == 0L,
            if (length(missing_output_columns)) paste(missing_output_columns, collapse = "|") else "none",
            "none")
  add_check(cohort_name, "DS_ALL_EXPOSURE_COLUMNS_PRESERVED",
            "Every insomnia exposure/provenance column in the authoritative panel is retained",
            length(missing_exposure_columns) == 0L,
            if (length(missing_exposure_columns)) paste(missing_exposure_columns, collapse = "|") else "none",
            "none")
  add_check(cohort_name, "DS_FIRST_RECORD_COVERAGE",
            "Exactly one first analysis record is identified for every participant",
            panel[is_first_analysis_record == 1L, .N] == uniqueN(panel$id),
            panel[is_first_analysis_record == 1L, .N], uniqueN(panel$id))
  add_check(cohort_name, "DS_UNIQUE_RECORD_KEYS", "No duplicated id-time-wave-event records",
            duplicate_keys == 0L, duplicate_keys, 0L)
  add_check(cohort_name, "DS_STRICT_TIME_ORDER", "Event time is strictly increasing within participant",
            bad_order == 0L, bad_order, 0L)
  add_check(cohort_name, "DS_STATE_MAPPING", "Known D states reproduce the three-component mapping",
            mapping_mismatch == 0L, mapping_mismatch, 0L)
  add_check(cohort_name, "DS_UNKNOWN_EXPLAINED", "Every unknown living D state has a missing component",
            unexplained_unknown == 0L, unexplained_unknown, 0L)
  add_check(cohort_name, "DS_COMPONENT_REVERSAL", "No cumulative disease component reverses from 1 to 0",
            reverse_rows[cohort == cohort_name, .N] == 0L,
            reverse_rows[cohort == cohort_name, .N], 0L)
  add_check(cohort_name, "DS_POSITIVE_INTERVALS", "All consecutive record intervals are positive",
            ints[!is.finite(interval_years) | interval_years <= 0, .N] == 0L,
            ints[!is.finite(interval_years) | interval_years <= 0, .N], 0L)
  add_check(cohort_name, "DS_MONOTONE_TRANSITIONS", "Every fully known observation is monotone or death",
            ints[transition_state_known == 1L & admissible_monotone_observation != 1L, .N] == 0L,
            ints[transition_state_known == 1L & admissible_monotone_observation != 1L, .N], 0L)
  add_check(cohort_name, "DS_TARGET_EDGE_CATALOG", "Exactly four confirmatory disease-specific edges are present",
            edges[target_edge == 1L, .N] == 4L, edges[target_edge == 1L, .N], 4L)
  add_check(cohort_name, "DS_TARGET_EVENT_SUPPORT", "Every target edge has at least one directly observed event",
            nrow(target_support) == 4L && all(target_support$direct_state_events > 0L),
            paste(target_support$direct_state_events, collapse = "|"), ">0 for each of four edges")
  add_check(cohort_name, "DS_NO_MULTICOMPONENT_Q_EDGE",
            "Allowed Q contains only one-component disease additions or death",
            edges[death_edge == 0L, {
              from <- state_dictionary[match(from_state, state), .(hypertension, heart_disease, stroke)]
              to <- state_dictionary[match(to_state, state), .(hypertension, heart_disease, stroke)]
              sum((to$hypertension - from$hypertension) +
                    (to$heart_disease - from$heart_disease) +
                    (to$stroke - from$stroke) != 1L)
            }] == 0L,
            edges[death_edge == 0L, .N], "all non-death edges add exactly one disease")
  add_check(cohort_name, "DS_Q_EDGE_COUNT", "Allowed Q edge count matches the prespecified lattice",
            nrow(edges) == if (cohort_name == "HRS") 20L else 12L,
            nrow(edges), if (cohort_name == "HRS") 20L else 12L)
  add_check(cohort_name, "DS_DEATH_SCOPE", "Cohort-specific death scope is enforced",
            if (cohort_name == "HRS") {
              panel[dead == 1L & disease_specific_state != 8L, .N] == 0L && death_not_final == 0L
            } else {
              panel[dead == 1L | event_type == "death" | !wave %in% 4:10, .N] == 0L
            },
            if (cohort_name == "HRS") paste0("death_not_final=", death_not_final) else
              paste0("out_of_scope_rows=", panel[dead == 1L | event_type == "death" | !wave %in% 4:10, .N]),
            if (cohort_name == "HRS") "all deaths map to D8 and are final" else
              "W4-W10 visits only; no death rows"
  )
}

for (cohort_name in names(input_files)) {
  before <- pre_fingerprints[cohort == cohort_name, sha256]
  after <- post_fingerprints[cohort == cohort_name, sha256]
  add_check(
    cohort_name, "DS_INPUT_IMMUTABLE",
    "Authoritative input fingerprint is unchanged by preparation",
    identical(before, after), after, before
  )
}

qc_checks <- rbindlist(qc_rows, use.names = TRUE, fill = TRUE)
critical_failures <- qc_checks[critical == 1L & status == "FAIL"]
if (nrow(critical_failures)) {
  print(critical_failures)
  stop(
    "Disease-specific data preparation fail-fast gate failed before publication (",
    nrow(critical_failures), " critical checks).",
    call. = FALSE
  )
}

published_files <- character()
publish_csv <- function(x, filename, root = qc_output_root) {
  path <- file.path(root, filename)
  atomic_fwrite(x, path)
  published_files <<- c(published_files, path)
  path
}
publish_rds <- function(x, filename, root = data_output_root) {
  path <- file.path(root, filename)
  atomic_save_rds(x, path)
  published_files <<- c(published_files, path)
  path
}

publish_csv(state_dictionary, "disease_specific_state_dictionary.csv", data_output_root)
publish_csv(allowed_edges, "disease_specific_allowed_q_edges.csv", data_output_root)
publish_csv(target_edge_catalog, "disease_specific_target_edge_catalog.csv", data_output_root)

for (cohort_name in names(prepared)) {
  stem <- tolower(cohort_name)
  panel <- prepared[[cohort_name]]$panel
  ints <- intervals[[cohort_name]]
  publish_rds(panel, paste0(stem, "_disease_specific_panel.rds"))
  publish_csv(panel, paste0(stem, "_disease_specific_panel.csv.gz"), data_output_root)
  publish_rds(ints, paste0(stem, "_disease_specific_intervals.rds"))
  publish_csv(ints, paste0(stem, "_disease_specific_intervals.csv.gz"), data_output_root)
}

publish_csv(input_fingerprints, "input_file_fingerprints.csv")
publish_csv(wave_state_counts, "disease_specific_wave_state_counts.csv")
publish_csv(combination_counts, "disease_specific_combination_counts.csv")
publish_csv(observed_transitions, "disease_specific_observed_transitions.csv")
publish_csv(target_edge_support, "disease_specific_target_edge_support.csv")
publish_csv(multijump_audit, "disease_specific_multijump_audit.csv")
publish_csv(unknown_rows, "disease_specific_unknown_rows.csv")
publish_csv(reverse_rows, "disease_specific_reverse_rows.csv")
publish_csv(unknown_reverse_summary, "disease_specific_unknown_reverse_summary.csv")

run_metadata <- data.table(
  field = c(
    "script_version", "run_timestamp_utc", "script_file", "script_sha256",
    "project_root", "authoritative_input_root", "output_root", "R_version",
    "data.table_version", "digest_version", "analysis_scope",
    "state_definition", "elsa_mortality_scope"
  ),
  value = c(
    script_version, run_timestamp_utc, script_file, sha256_file(script_file),
    project_root, normalizePath(input_root, winslash = "/"),
    normalizePath(nextsteps_root, winslash = "/", mustWork = TRUE),
    R.version.string, as.character(packageVersion("data.table")),
    as.character(packageVersion("digest")),
    "Data preparation and audit only; no model fitted",
    "D0-D7 are all combinations of cumulative hypertension, broad heart disease and stroke; HRS D8 is death",
    "Primary ELSA disease-specific panel is restricted to living W4-W10 visits; source deaths are audited but excluded"
  )
)
publish_csv(run_metadata, "run_metadata.csv")

# Verify serialised row counts before publishing the final QC gate.
for (cohort_name in names(prepared)) {
  stem <- tolower(cohort_name)
  panel_path <- file.path(data_output_root, paste0(stem, "_disease_specific_panel.rds"))
  interval_path <- file.path(data_output_root, paste0(stem, "_disease_specific_intervals.rds"))
  panel_rows_readback <- nrow(readRDS(panel_path))
  interval_rows_readback <- nrow(readRDS(interval_path))
  add_check(
    cohort_name, "DS_OUTPUT_PANEL_READBACK",
    "Published panel RDS row count reproduces the in-memory object",
    panel_rows_readback == nrow(prepared[[cohort_name]]$panel),
    panel_rows_readback, nrow(prepared[[cohort_name]]$panel)
  )
  add_check(
    cohort_name, "DS_OUTPUT_INTERVAL_READBACK",
    "Published interval RDS row count reproduces the in-memory object",
    interval_rows_readback == nrow(intervals[[cohort_name]]),
    interval_rows_readback, nrow(intervals[[cohort_name]])
  )
  add_check(
    cohort_name, "DS_INPUT_IMMUTABLE_AFTER_WRITE",
    "Authoritative input fingerprint remains unchanged after output publication",
    sha256_file(input_files[[cohort_name]]) == pre_fingerprints[cohort == cohort_name, sha256],
    sha256_file(input_files[[cohort_name]]),
    pre_fingerprints[cohort == cohort_name, sha256]
  )
}

qc_checks <- rbindlist(qc_rows, use.names = TRUE, fill = TRUE)
critical_failures <- qc_checks[critical == 1L & status == "FAIL"]
if (nrow(critical_failures)) {
  print(critical_failures)
  stop(
    "Disease-specific output verification fail-fast gate failed (",
    nrow(critical_failures), " critical checks).",
    call. = FALSE
  )
}
publish_csv(qc_checks, "disease_specific_qc_checks.csv")

output_manifest <- rbindlist(lapply(sort(unique(published_files)), function(path) {
  info <- file.info(path)
  data.table(
    relative_path = substring(
      normalizePath(path, winslash = "/", mustWork = TRUE),
      nchar(normalizePath(nextsteps_root, winslash = "/", mustWork = TRUE)) + 2L
    ),
    bytes = as.numeric(info$size),
    sha256 = sha256_file(path),
    modified_utc = format(info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
}))
atomic_fwrite(output_manifest, file.path(qc_output_root, "output_manifest.csv"))

cat("Disease-specific state preparation completed.\n")
cat("Authoritative inputs:", normalizePath(input_root, winslash = "/"), "\n")
cat("Output data:", normalizePath(data_output_root, winslash = "/"), "\n")
cat("QC output:", normalizePath(qc_output_root, winslash = "/"), "\n")
cat("Critical QC failures: 0\n")
print(target_edge_support[sample == "MI20_eligible_full_panel", .(
  cohort, target_id, transition, direct_state_events,
  target_class_events_including_concurrent_hypertension,
  joint_heart_stroke_sequence_unknown_events, person_years
)])
