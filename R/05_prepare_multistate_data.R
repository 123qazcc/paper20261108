#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GlobalAgingR)
})

options(stringsAsFactors = FALSE, width = 140)

stage_root <- Sys.getenv("GA_STAGE_ROOT", unset = "private/stage")
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
dir.create(data_root, recursive = TRUE, showWarnings = FALSE)
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

raw_numeric <- function(x) {
  if (is.null(x)) return(numeric())
  y <- unclass(x)
  if (is.factor(y)) y <- as.character(y)
  suppressWarnings(as.numeric(y))
}

valid_number <- function(x) {
  x <- raw_numeric(x)
  x[!is.finite(x)] <- NA_real_
  x
}

locf_down <- function(x) {
  if (!length(x)) return(x)
  last <- NA
  for (i in seq_along(x)) {
    if (!is.na(x[i])) last <- x[i] else if (!is.na(last)) x[i] <- last
  }
  x
}

read_package_data <- function(item) {
  env <- new.env(parent = baseenv())
  utils::data(list = item, package = "GlobalAgingR", envir = env)
  get(item, envir = env, inherits = FALSE)
}

make_package_supplement <- function(cohort) {
  x <- read_package_data(tolower(cohort))
  if (cohort == "HRS") {
    out <- data.table(
      id = raw_numeric(x$hhidpn),
      wave = as.integer(raw_numeric(x$wave)),
      household_wealth = valid_number(x$atotw),
      household_income = valid_number(x$itot)
    )
  } else {
    out <- data.table(
      id = raw_numeric(x$idauniqc),
      wave = as.integer(raw_numeric(x$wave)),
      household_wealth = valid_number(x$hatotb),
      household_income = valid_number(x$hitot)
    )
  }
  unique(out[!is.na(id) & !is.na(wave)], by = c("id", "wave"))
}

rank_quintile <- function(x) {
  out <- rep(NA_integer_, length(x))
  ok <- !is.na(x) & is.finite(x)
  if (!any(ok)) return(out)
  r <- rank(x[ok], ties.method = "average")
  out[ok] <- pmin(5L, pmax(1L, ceiling(5 * r / length(r))))
  out
}

state_labels <- c(
  `0` = "S0: No target disease",
  `1` = "S1: Hypertension only",
  `2` = "S2: One CVD without hypertension",
  `3` = "S3: Multimorbidity",
  `4` = "S4: Death"
)

transition_catalog <- data.table(
  transition_id = sprintf("T%02d", 1:8),
  from_state = c(0L, 0L, 0L, 1L, 1L, 2L, 2L, 3L),
  to_state = c(1L, 2L, 4L, 3L, 4L, 3L, 4L, 4L)
)
transition_catalog[, `:=`(
  from_label = unname(state_labels[as.character(from_state)]),
  to_label = unname(state_labels[as.character(to_state)]),
  transition = paste0("S", from_state, "→S", to_state),
  death_transition = as.integer(to_state == 4L)
)]

prepare_cohort <- function(cohort_name) {
  message("Preparing multistate data for ", cohort_name, "...")
  panel_file <- file.path(exposure_root, paste0(tolower(cohort_name), "_exposure_panel.rds"))
  panel <- as.data.table(readRDS(panel_file))
  supplement <- make_package_supplement(cohort_name)
  panel <- merge(panel, supplement, by = c("id", "wave"), all.x = TRUE, sort = FALSE)
  panel[, event_time_original := as.numeric(event_time)]
  panel[, death_time_method_analysis := fifelse(
    dead == 1L, "Death-year midpoint from harmonized data", NA_character_
  )]
  if (cohort_name == "HRS") {
    death_timing <- panel[dead == 0L, .(
      last_alive_time = max(as.numeric(event_time), na.rm = TRUE),
      death_span_months = {
        z <- as.numeric(radtimtdth[!is.na(radtimtdth)])
        if (length(z)) z[1L] else NA_real_
      }
    ), by = id]
    death_timing[, death_time_from_span := last_alive_time + death_span_months / 12]
    panel <- merge(panel, death_timing, by = "id", all.x = TRUE, sort = FALSE)
    panel[
      dead == 1L & is.finite(death_time_from_span),
      `:=`(
        event_time = pmax(death_time_from_span, last_alive_time + 1 / 365.25),
        death_time_method_analysis = "Last interview plus RADTIMTDTH months"
      )
    ]
  } else {
    panel[, `:=`(
      last_alive_time = NA_real_, death_span_months = NA_real_, death_time_from_span = NA_real_
    )]
  }
  setorder(panel, id, event_time, dead, wave, na.last = TRUE)
  panel[, `:=`(
    household_wealth_current = locf_down(household_wealth),
    household_income_current = locf_down(household_income)
  ), by = id]

  baseline <- panel[analysis_baseline == 1L]
  setorder(baseline, id, event_time)
  baseline <- baseline[, .SD[1L], by = id]
  baseline[, female_bl := fifelse(!is.na(sex_male), as.integer(sex_male == 0L), NA_integer_)]
  baseline[, race_bl := if (cohort_name == "HRS") {
    fcase(
      hispanic == 1L, "Hispanic",
      hispanic == 0L & race_code == 1L, "Non-Hispanic White",
      hispanic == 0L & race_code == 2L, "Non-Hispanic Black",
      hispanic == 0L & race_code == 3L, "Other",
      default = NA_character_
    )
  } else {
    fifelse(is.na(race_white), NA_character_, fifelse(race_white == 1L, "White", "Non-White"))
  }]
  baseline[, education_level_bl := fifelse(raeducl %in% 1:3, as.integer(raeducl), NA_integer_)]
  baseline[, married_bl := fifelse(
    is.na(marital_status), NA_integer_, as.integer(marital_status %in% c(1L, 2L, 3L))
  )]
  baseline[, inactive_bl := if (cohort_name == "HRS") {
    fifelse(
      !is.na(vigorous_activity) & !is.na(moderate_activity),
      as.integer(vigorous_activity == 0L & moderate_activity == 0L),
      NA_integer_
    )
  } else {
    fifelse(
      !is.na(vigorous_activity) & !is.na(moderate_activity),
      as.integer(vigorous_activity == 5L & moderate_activity == 5L),
      NA_integer_
    )
  }]
  baseline[, weight_bl := as.numeric(sampling_weight)]
  baseline[!is.finite(weight_bl) | weight_bl <= 0, weight_bl := NA_real_]
  baseline[, wealth_quintile_bl := rank_quintile(household_wealth_current)]
  baseline[, income_quintile_bl := rank_quintile(household_income_current)]
  baseline[, `:=`(
    age_bl = as.numeric(age),
    bmi_bl = as.numeric(bmi),
    current_smoker_bl = as.integer(current_smoker),
    current_drinker_bl = as.integer(current_drinker),
    diabetes_bl = as.integer(diabetes),
    lung_disease_bl = as.integer(lung_disease),
    frequent_pain_bl = as.integer(frequent_pain),
    depression_score_bl = as.numeric(depression_score),
    insomnia_score_bl = as.integer(insomnia_score_primary),
    insomnia_burden_bl = as.integer(insomnia_burden_code),
    state_bl = as.integer(state),
    entry_wave = as.integer(wave),
    entry_year = as.numeric(event_time)
  )]

  baseline_vars <- c(
    "id", "cohort", "entry_wave", "entry_year", "state_bl", "weight_bl",
    "insomnia_score_bl", "insomnia_burden_bl", "age_bl", "female_bl", "race_bl",
    "education_level_bl", "wealth_quintile_bl", "income_quintile_bl", "married_bl",
    "current_smoker_bl", "current_drinker_bl", "inactive_bl", "bmi_bl", "diabetes_bl",
    "lung_disease_bl", "frequent_pain_bl", "depression_score_bl", "household_wealth_current",
    "household_income_current"
  )
  baseline <- baseline[, ..baseline_vars]

  state_rows <- panel[
    analysis_sample == 1L & !is.na(exposure_entry_time) &
      event_time >= exposure_entry_time & !is.na(state)
  ]
  setorder(state_rows, id, event_time, dead, wave, na.last = TRUE)
  state_rows[, `:=`(
    next_state = shift(as.integer(state), type = "lead"),
    next_time = shift(as.numeric(event_time), type = "lead"),
    next_wave = shift(as.integer(wave), type = "lead"),
    next_event_type = shift(as.character(event_type), type = "lead")
  ), by = id]

  exposure_vars <- c(
    "insomnia_score_primary", "insomnia_burden_code", "insomnia_high_primary",
    "insomnia_onset_primary", "insomnia_maintenance_primary", "insomnia_nonrestorative_primary",
    "insomnia_score_cumavg", "insomnia_score_prev_valid", "insomnia_score_lag1_row",
    "insomnia_dynamic_pattern", "insomnia_dynamic_code", "insomnia_assessment_complete",
    "insomnia_assessment_number", "insomnia_dynamic_source_wave", "insomnia_assessment_gap_waves",
    "insomnia_source_wave", "insomnia_locf_reason", "insomnia_score_extended4",
    "insomnia_early_waking_current", "sleep_medication_current", "proxy_interview"
  )
  missing_exposure_vars <- setdiff(exposure_vars, names(state_rows))
  if (length(missing_exposure_vars)) {
    state_rows[, (missing_exposure_vars) := NA]
  }
  interval_cols <- c(
    "id", "cohort", "wave", "next_wave", "event_type", "next_event_type",
    "event_time", "next_time", "exposure_entry_time", "state", "next_state", exposure_vars
  )
  intervals <- state_rows[!is.na(next_state) & !is.na(next_time), ..interval_cols]
  intervals[, `:=`(
    tstart = as.numeric(event_time - exposure_entry_time),
    tstop = as.numeric(next_time - exposure_entry_time),
    interval_years = as.numeric(next_time - event_time),
    from_state = as.integer(state),
    to_state = as.integer(next_state),
    state_changed = as.integer(next_state != state)
  )]
  intervals[, observed_transition := fifelse(
    state_changed == 1L, paste0("S", from_state, "→S", to_state), "No transition"
  )]
  intervals <- merge(intervals, baseline, by = c("id", "cohort"), all.x = TRUE, sort = FALSE)
  setorder(intervals, id, tstart, tstop)

  intensity_keys <- transition_catalog$transition
  admissible_observation_keys <- c(intensity_keys, "S0→S3")
  intervals[, allowed_intensity_edge := as.integer(observed_transition %chin% intensity_keys)]
  intervals[, interval_censored_multistep_jump := as.integer(observed_transition == "S0→S3")]
  intervals[, admissible_observation := as.integer(
    state_changed == 0L | observed_transition %chin% admissible_observation_keys
  )]
  intervals[, early_period := as.integer(tstop <= 2)]
  intervals[, interval_starts_at_observed_sleep := as.integer(insomnia_assessment_complete == 1L)]

  transition_counts <- rbindlist(lapply(seq_len(nrow(transition_catalog)), function(i) {
    tr <- transition_catalog[i]
    risk <- intervals[from_state == tr$from_state]
    events <- risk[to_state == tr$to_state & state_changed == 1L]
    data.table(
      cohort = cohort_name,
      transition_id = tr$transition_id,
      transition = tr$transition,
      from_state = tr$from_state,
      to_state = tr$to_state,
      risk_intervals = nrow(risk),
      risk_participants = uniqueN(risk$id),
      events = nrow(events),
      event_participants = uniqueN(events$id),
      person_years = sum(risk$interval_years, na.rm = TRUE),
      rate_per_1000_py = 1000 * nrow(events) / sum(risk$interval_years, na.rm = TRUE),
      weighted_events = sum(events$weight_bl, na.rm = TRUE),
      weighted_person_years = sum(risk$weight_bl * risk$interval_years, na.rm = TRUE),
      weighted_rate_per_1000_py = 1000 * sum(events$weight_bl, na.rm = TRUE) /
        sum(risk$weight_bl * risk$interval_years, na.rm = TRUE)
    )
  }))

  burden_counts <- rbindlist(lapply(seq_len(nrow(transition_catalog)), function(i) {
    tr <- transition_catalog[i]
    risk <- intervals[from_state == tr$from_state]
    risk[, .(
      risk_intervals = .N,
      risk_participants = uniqueN(id),
      events = sum(to_state == tr$to_state & state_changed == 1L),
      person_years = sum(interval_years),
      rate_per_1000_py = 1000 * sum(to_state == tr$to_state & state_changed == 1L) / sum(interval_years)
    ), by = .(cohort, insomnia_burden_code)][, `:=`(
      transition_id = tr$transition_id,
      transition = tr$transition
    )]
  }))
  setcolorder(burden_counts, c(
    "cohort", "transition_id", "transition", "insomnia_burden_code",
    "risk_intervals", "risk_participants", "events", "person_years", "rate_per_1000_py"
  ))

  observed_counts <- intervals[, .(
    observed_intervals = .N,
    participants = uniqueN(id),
    person_years = sum(interval_years, na.rm = TRUE)
  ), by = .(cohort, from_state, to_state, observed_transition)]
  observed_counts[, `:=`(
    state_changed = as.integer(from_state != to_state),
    model_intensity_edge = as.integer(observed_transition %chin% intensity_keys),
    interval_censored_multistep_jump = as.integer(observed_transition == "S0→S3")
  )]

  panel_msm <- copy(state_rows)
  panel_msm[, c("next_state", "next_time", "next_wave", "next_event_type") := NULL]
  fill_exposure_vars <- intersect(exposure_vars, names(panel_msm))
  panel_msm[, (fill_exposure_vars) := lapply(.SD, locf_down), by = id, .SDcols = fill_exposure_vars]
  panel_msm[, `:=`(
    time_msm = as.numeric(event_time - exposure_entry_time),
    state_msm = as.integer(state) + 1L,
    age_tv = as.numeric(age),
    bmi_tv = as.numeric(bmi),
    current_smoker_tv = as.integer(current_smoker),
    current_drinker_tv = as.integer(current_drinker),
    diabetes_tv = as.integer(diabetes),
    lung_disease_tv = as.integer(lung_disease),
    frequent_pain_tv = as.integer(frequent_pain),
    depression_score_tv = as.numeric(depression_score)
  )]
  panel_msm[, inactive_tv := if (cohort_name == "HRS") {
    fifelse(
      !is.na(vigorous_activity) & !is.na(moderate_activity),
      as.integer(vigorous_activity == 0L & moderate_activity == 0L), NA_integer_
    )
  } else {
    fifelse(
      !is.na(vigorous_activity) & !is.na(moderate_activity),
      as.integer(vigorous_activity == 5L & moderate_activity == 5L), NA_integer_
    )
  }]
  tv_fill_vars <- c(
    "age_tv", "bmi_tv", "current_smoker_tv", "current_drinker_tv", "inactive_tv",
    "diabetes_tv", "lung_disease_tv", "frequent_pain_tv", "depression_score_tv"
  )
  panel_msm[, (tv_fill_vars) := lapply(.SD, locf_down), by = id, .SDcols = tv_fill_vars]
  panel_msm <- merge(panel_msm, baseline, by = c("id", "cohort"), all.x = TRUE, sort = FALSE)
  setorder(panel_msm, id, time_msm, state_msm)

  required_primary <- c(
    "age_bl", "female_bl", "race_bl", "education_level_bl", "wealth_quintile_bl",
    "married_bl", "current_smoker_bl", "current_drinker_bl", "inactive_bl", "bmi_bl",
    "diabetes_bl", "lung_disease_bl", "frequent_pain_bl"
  )
  flow <- data.table(
    cohort = cohort_name,
    step = c(
      "Exposure panel participants", "Exposure-valid analysis participants",
      "Positive baseline sampling weight", "Complete primary covariates",
      "At least two complete sleep assessments", "At least one observed transition",
      "Observed death transitions"
    ),
    n = c(
      uniqueN(panel$id), nrow(baseline), sum(!is.na(baseline$weight_bl)),
      sum(complete.cases(baseline[, ..required_primary])),
      panel[, uniqueN(id[insomnia_assessment_number >= 2L])],
      intervals[state_changed == 1L, uniqueN(id)],
      intervals[to_state == 4L & state_changed == 1L, uniqueN(id)]
    )
  )

  qc <- data.table(
    cohort = cohort_name,
    check = c(
      "One baseline row per final participant", "Positive interval duration",
      "No forbidden or reverse panel observation", "Primary exposure present at interval start",
      "Interval times ordered", "State code in S0-S4", "All transition events catalogued",
      "Fixed baseline weight does not vary within participant",
      "MSM time is nonnegative and ordered", "MSM death rows retain interval-start exposure",
      "Death is the final record", "S0-S3 jump is not parameterized as a direct intensity"
    ),
    violations = c(
      nrow(baseline) - uniqueN(baseline$id),
      intervals[!is.finite(interval_years) | interval_years <= 0, .N],
      intervals[admissible_observation != 1L, .N],
      intervals[is.na(insomnia_score_primary), .N],
      intervals[tstop <= tstart, .N],
      intervals[!from_state %in% 0:3 | !to_state %in% 0:4, .N],
      intervals[state_changed == 1L & !observed_transition %chin% admissible_observation_keys, .N],
      intervals[, .(n_weight_values = uniqueN(weight_bl, na.rm = FALSE)), by = id][
        n_weight_values > 1L, .N
      ],
      panel_msm[, .(bad = any(diff(time_msm) <= 0) | any(time_msm < 0)), by = id][bad == TRUE, .N],
      panel_msm[state_msm == 5L & is.na(insomnia_score_primary), .N],
      panel_msm[, .(death_not_final = any(state_msm == 5L & seq_len(.N) != .N)), by = id][death_not_final == TRUE, .N],
      nrow(transition_catalog[from_state == 0L & to_state == 3L])
    )
  )
  qc[, status := fifelse(violations == 0L, "PASS", "FAIL")]

  saveRDS(panel, file.path(data_root, paste0(tolower(cohort_name), "_exposure_panel_augmented.rds")), compress = "xz")
  saveRDS(baseline, file.path(data_root, paste0(tolower(cohort_name), "_baseline_covariates.rds")), compress = "xz")
  saveRDS(intervals, file.path(data_root, paste0(tolower(cohort_name), "_multistate_intervals.rds")), compress = "xz")
  saveRDS(panel_msm, file.path(data_root, paste0(tolower(cohort_name), "_msm_panel.rds")), compress = "xz")
  fwrite(baseline, file.path(data_root, paste0(tolower(cohort_name), "_baseline_covariates.csv.gz")), na = "")
  fwrite(intervals, file.path(data_root, paste0(tolower(cohort_name), "_multistate_intervals.csv.gz")), na = "")
  fwrite(panel_msm, file.path(data_root, paste0(tolower(cohort_name), "_msm_panel.csv.gz")), na = "")
  fwrite(transition_counts, file.path(result_root, paste0(tolower(cohort_name), "_transition_counts_rates.csv")), na = "")
  fwrite(burden_counts, file.path(result_root, paste0(tolower(cohort_name), "_transition_rates_by_burden.csv")), na = "")
  fwrite(observed_counts, file.path(result_root, paste0(tolower(cohort_name), "_observed_panel_transition_counts.csv")), na = "")
  fwrite(flow, file.path(result_root, paste0(tolower(cohort_name), "_analysis_flow.csv")), na = "")
  fwrite(qc, file.path(qc_root, paste0(tolower(cohort_name), "_multistate_data_qc.csv")), na = "")

  list(
    cohort = cohort_name,
    baseline = baseline,
    intervals = intervals,
    panel_msm = panel_msm,
    transition_counts = transition_counts,
    burden_counts = burden_counts,
    observed_counts = observed_counts,
    flow = flow,
    qc = qc
  )
}

hrs <- prepare_cohort("HRS")
elsa <- prepare_cohort("ELSA")

combined_counts <- rbindlist(list(hrs$transition_counts, elsa$transition_counts), use.names = TRUE)
combined_burden <- rbindlist(list(hrs$burden_counts, elsa$burden_counts), use.names = TRUE)
combined_observed <- rbindlist(list(hrs$observed_counts, elsa$observed_counts), use.names = TRUE)
combined_flow <- rbindlist(list(hrs$flow, elsa$flow), use.names = TRUE)
combined_qc <- rbindlist(list(hrs$qc, elsa$qc), use.names = TRUE)

fwrite(transition_catalog, file.path(result_root, "transition_catalog.csv"), na = "")
fwrite(combined_counts, file.path(result_root, "transition_counts_rates.csv"), na = "")
fwrite(combined_burden, file.path(result_root, "transition_rates_by_burden.csv"), na = "")
fwrite(combined_observed, file.path(result_root, "observed_panel_transition_counts.csv"), na = "")
fwrite(combined_flow, file.path(result_root, "analysis_flow.csv"), na = "")
fwrite(combined_qc, file.path(qc_root, "multistate_data_qc.csv"), na = "")

globalagingr_usage <- data.table(
  component = c(
    "Longitudinal cohort spine", "Filtering/selecting harmonized variables",
    "Markov Q-matrix initialization", "Continuous-time Markov fitting",
    "Markov diagnostics/prediction", "Transition-specific hazard models"
  ),
  GlobalAgingR_function = c(
    "hrs / elsa package datasets", "CH_filter / CH_select_columns",
    "estimate_Q_matrix", "msm_transition_model",
    "msm_diagnostics / predict_msm_transitions", "No complete package wrapper"
  ),
  implementation = c(
    "Package data are the harmonized longitudinal spine; original ELSA core items remain the sleep source.",
    "Stage-1 cleaning reused package selection/filtering logic.",
    "Called on recoded S0-S4 data, then constrained to the protocol transition graph.",
    "Called for cohort-specific continuous-time Markov models with the constrained Q matrix.",
    "Called where model convergence permits; outputs supplemented with tidy tables and publication plots.",
    "Supplemented with survival::coxph to support sampling weights, clustered robust SEs and transition-specific effects."
  )
)
fwrite(globalagingr_usage, file.path(result_root, "globalagingr_functions_used.csv"), na = "")

print(combined_flow)
print(combined_counts)
print(combined_qc)
if (combined_qc[status == "FAIL", .N] > 0L) {
  stop("Multistate data QC failed. Inspect multistate_data_qc.csv")
}
cat("Multistate interval data prepared successfully.\n")
