#!/usr/bin/env Rscript

# Supplementary protocol sensitivity analyses for the HRS-ELSA insomnia/CVD
# multistate study.
#
# Design principles:
#   * participant id is the independent unit; repeated risk intervals use
#     cluster-robust covariance;
#   * the frozen eight-edge graph is never expanded (S0->S3 remains a censored
#     panel-observed multistep jump, never a target intensity);
#   * effect estimates, 95% confidence intervals and exact two-sided p values
#     are retained together;
#   * unsupported outcome redefinitions and emission/misclassification models
#     are reported as not estimable rather than populated with assumed values.

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

options(stringsAsFactors = FALSE, width = 180)

script_version <- "2026-09-09-v3"
stage_root <- Sys.getenv("GA_STAGE_ROOT", unset = "private/stage")
output_root <- Sys.getenv(
  "GA_MULTISTATE_OUTPUT_ROOT",
  unset = file.path(stage_root, "multistate_output")
)
data_root <- file.path(output_root, "data_multistate")
result_root <- file.path(output_root, "results", "protocol_sensitivities")
report_root <- file.path(output_root, "reports")
qc_root <- file.path(output_root, "qc", "protocol_sensitivities")
model_root <- file.path(output_root, "models", "protocol_sensitivities")

for (path in c(result_root, report_root, qc_root, model_root)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

run_mode <- tolower(trimws(Sys.getenv("GA_PROTOCOL_SENS_RUN_MODE", unset = "full")))
if (!run_mode %in% c("audit", "light", "full")) {
  stop("GA_PROTOCOL_SENS_RUN_MODE must be audit, light, or full")
}

split_env <- function(name, default) {
  value <- trimws(Sys.getenv(name, unset = default))
  value <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  value[nzchar(value)]
}

requested_cohorts <- toupper(split_env("GA_PROTOCOL_SENS_COHORTS", "HRS,ELSA"))
requested_scenarios <- tolower(split_env(
  "GA_PROTOCOL_SENS_SCENARIOS",
  paste(
    c(
      "rcs", "alternate_death", "attained_age_period", "ipcw",
      "proxy_exclusion", "expanded_hypertension", "heart_lower",
      "heart_upper", "semi_markov"
    ),
    collapse = ","
  )
))
max_models <- suppressWarnings(as.integer(Sys.getenv(
  "GA_PROTOCOL_SENS_MAX_MODELS",
  unset = if (run_mode == "light") "2" else "0"
)))
if (is.na(max_models) || max_models < 0L) max_models <- 0L
min_events <- suppressWarnings(as.integer(Sys.getenv("GA_PROTOCOL_SENS_MIN_EVENTS", unset = "10")))
if (is.na(min_events) || min_events < 1L) min_events <- 10L
min_epp <- suppressWarnings(as.numeric(Sys.getenv("GA_PROTOCOL_SENS_MIN_EPP", unset = "5")))
if (!is.finite(min_epp) || min_epp <= 0) min_epp <- 5
elsa_death_cutoff <- suppressWarnings(as.numeric(Sys.getenv(
  "GA_ELSA_DEATH_CUTOFF", unset = "2012.999"
)))
if (!is.finite(elsa_death_cutoff)) stop("GA_ELSA_DEATH_CUTOFF must be a finite numeric year")
save_models <- tolower(Sys.getenv("GA_PROTOCOL_SENS_SAVE_MODELS", unset = "1")) %in%
  c("1", "true", "yes", "y", "on")

expected_edges <- data.table(
  transition_id = sprintf("T%02d", 1:8),
  from_state = c(0L, 0L, 0L, 1L, 1L, 2L, 2L, 3L),
  to_state = c(1L, 2L, 4L, 3L, 4L, 3L, 4L, 4L)
)
expected_edges[, `:=`(
  transition = paste0("S", from_state, "->S", to_state),
  death_transition = as.integer(to_state == 4L),
  confirmatory_nonfatal = as.integer(
    paste0(from_state, "->", to_state) %chin% c("0->1", "0->2", "1->3", "2->3")
  )
)]

# Protocol item 7 is an ELSA-only bounded-definition analysis.  The lower and
# upper state columns are constructed upstream from the same source-backed
# component diagnoses as the primary state; this script deliberately does not
# reinterpret the component data a second time.
heart_definition_scenarios <- data.table(
  scenario = c("heart_lower", "heart_upper"),
  sensitivity_id = 7L,
  state_column = c("state_heart_lower", "state_heart_upper"),
  label = c(
    "ELSA lower-bound heart phenotype for unresolved historical other-heart status",
    "ELSA upper-bound heart phenotype for unresolved historical other-heart status"
  ),
  bound_assumption = c(
    "Unresolved historical other-heart component coded absent",
    "Unresolved historical other-heart component coded present"
  )
)

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

normalise_status <- function(x) {
  z <- gsub("[^a-z0-9]+", "_", tolower(clean_text(x)))
  fcase(
    grepl("success|fitted|completed|rubin_pooled_complete", z), "completed",
    grepl("not_estim|skip|disabled|not_run|planned", z), "not_estimable",
    grepl("fail|error|timeout|interrupt", z), "failed",
    default = z
  )
}

atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = basename(path), tmpdir = dirname(path), fileext = ".tmp")
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  fwrite(as.data.table(x), tmp, na = "")
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not write ", path)
  invisible(path)
}

safe_fread <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0L) return(data.table())
  tryCatch(suppressWarnings(fread(path, showProgress = FALSE)), error = function(e) data.table())
}

read_required_rds <- function(path) {
  if (!file.exists(path)) stop("Missing required staged input: ", path)
  as.data.table(readRDS(path))
}

cum_diagnosis <- function(x) {
  x <- as.integer(x)
  out <- rep(NA_integer_, length(x))
  ever <- FALSE
  known <- FALSE
  for (i in seq_along(x)) {
    if (!is.na(x[i])) {
      known <- TRUE
      if (x[i] == 1L) ever <- TRUE
    }
    if (ever) out[i] <- 1L else if (known) out[i] <- 0L
  }
  out
}

derive_state <- function(hypertension, heart, stroke, dead) {
  h <- as.integer(hypertension)
  c <- as.integer(heart)
  s <- as.integer(stroke)
  d <- as.integer(dead)
  values <- cbind(h, c, s)
  positive <- rowSums(values == 1L, na.rm = TRUE)
  known <- rowSums(!is.na(values))
  state <- rep(NA_integer_, length(h))
  state[d == 1L] <- 4L
  alive <- d != 1L
  state[alive & known == 3L & positive == 0L] <- 0L
  state[alive & h == 1L & c == 0L & s == 0L] <- 1L
  state[alive & known == 3L & h == 0L & positive == 1L & (c == 1L | s == 1L)] <- 2L
  state[alive & positive >= 2L] <- 3L
  state
}

wald_quadratic <- function(beta, variance) {
  if (!length(beta) || !all(is.finite(beta)) || !all(is.finite(variance))) {
    return(list(statistic = NA_real_, df = NA_integer_, p_value = NA_real_))
  }
  variance <- (variance + t(variance)) / 2
  eig <- eigen(variance, symmetric = TRUE)
  tolerance <- max(abs(eig$values), 1) * max(dim(variance)) * sqrt(.Machine$double.eps)
  keep <- eig$values > tolerance
  rank <- sum(keep)
  if (rank < 1L) return(list(statistic = NA_real_, df = 0L, p_value = NA_real_))
  projected <- as.numeric(crossprod(eig$vectors[, keep, drop = FALSE], beta))
  statistic <- sum(projected^2 / eig$values[keep])
  list(statistic = statistic, df = rank, p_value = pchisq(statistic, rank, lower.tail = FALSE))
}

rcs_basis <- function(x, knots = c(0, 1, 2, 3)) {
  x <- as.numeric(x)
  knots <- sort(unique(as.numeric(knots)))
  if (length(knots) != 4L || any(!is.finite(knots)) || diff(range(knots)) <= 0) {
    stop("The prespecified insomnia RCS requires four distinct knots at scores 0, 1, 2 and 3")
  }
  positive_cube <- function(z) pmax(z, 0)^3
  k <- length(knots)
  scale <- (knots[k] - knots[1L])^2
  nonlinear <- lapply(seq_len(k - 2L), function(j) {
    (positive_cube(x - knots[j]) -
       ((knots[k] - knots[j]) / (knots[k] - knots[k - 1L])) * positive_cube(x - knots[k - 1L]) +
       ((knots[k - 1L] - knots[j]) / (knots[k] - knots[k - 1L])) * positive_cube(x - knots[k])) / scale
  })
  data.table(
    rcs_linear = x,
    rcs_nl1 = nonlinear[[1L]],
    rcs_nl2 = nonlinear[[2L]]
  )
}

prepare_model_columns <- function(x) {
  x <- copy(as.data.table(x))
  numeric_columns <- intersect(c(
    "tstart", "tstop", "event_time", "age_bl", "bmi_bl", "female_bl",
    "married_bl", "current_smoker_bl", "current_drinker_bl", "inactive_bl",
    "diabetes_bl", "lung_disease_bl", "frequent_pain_bl",
    "insomnia_score_primary", "proxy_interview", "ipcw"
  ), names(x))
  x[, (numeric_columns) := lapply(.SD, as.numeric), .SDcols = numeric_columns]
  x[, `:=`(
    age_bl_5 = (age_bl - 65) / 5,
    bmi_bl_5 = (bmi_bl - 25) / 5,
    female_factor = factor(female_bl, levels = c(0, 1), labels = c("Male", "Female")),
    race_factor = factor(race_bl),
    education_factor = factor(education_level_bl),
    wealth_factor = factor(wealth_quintile_bl),
    married_factor = factor(married_bl, levels = c(0, 1)),
    current_smoker_factor = factor(current_smoker_bl, levels = c(0, 1)),
    current_drinker_factor = factor(current_drinker_bl, levels = c(0, 1)),
    inactive_factor = factor(inactive_bl, levels = c(0, 1)),
    diabetes_factor = factor(diabetes_bl, levels = c(0, 1)),
    lung_factor = factor(lung_disease_bl, levels = c(0, 1)),
    pain_factor = factor(frequent_pain_bl, levels = c(0, 1)),
    calendar_period = cut(
      event_time,
      breaks = c(-Inf, 2012, 2016, Inf), right = FALSE,
      labels = c("before_2012", "2012_2015", "2016_plus")
    )
  )]
  rcs <- rcs_basis(x$insomnia_score_primary)
  x[, `:=`(
    rcs_linear = rcs$rcs_linear,
    rcs_nl1 = rcs$rcs_nl1,
    rcs_nl2 = rcs$rcs_nl2
  )]
  x[]
}

m4_terms <- c(
  "age_bl_5", "female_factor", "race_factor", "education_factor", "wealth_factor",
  "married_factor", "current_smoker_factor", "current_drinker_factor",
  "inactive_factor", "bmi_bl_5", "diabetes_factor", "lung_factor", "pain_factor"
)
m4_terms_no_age <- setdiff(m4_terms, "age_bl_5")

build_intervals_from_panel <- function(panel, state_column) {
  x <- copy(panel)
  if (!state_column %in% names(x)) stop("State column missing: ", state_column)
  x[, state_sensitivity := as.integer(get(state_column))]
  x <- x[
    analysis_sample == 1L & !is.na(exposure_entry_time) &
      event_time >= exposure_entry_time & !is.na(state_sensitivity)
  ]
  setorder(x, id, event_time, dead, wave, na.last = TRUE)
  x[, `:=`(
    next_state_sensitivity = shift(state_sensitivity, type = "lead"),
    next_time_sensitivity = shift(as.numeric(event_time), type = "lead")
  ), by = id]
  keep_columns <- intersect(c(
    "id", "cohort", "wave", "event_time", "exposure_entry_time",
    "insomnia_score_primary", "proxy_interview"
  ), names(x))
  out <- x[!is.na(next_state_sensitivity) & !is.na(next_time_sensitivity), c(
    keep_columns, "state_sensitivity", "next_state_sensitivity", "next_time_sensitivity"
  ), with = FALSE]
  setnames(out, c("state_sensitivity", "next_state_sensitivity"), c("from_state", "to_state"))
  out[, `:=`(
    tstart = as.numeric(event_time - exposure_entry_time),
    tstop = as.numeric(next_time_sensitivity - exposure_entry_time),
    interval_years = as.numeric(next_time_sensitivity - event_time),
    state_changed = as.integer(from_state != to_state)
  )]
  baseline <- read_required_rds(file.path(
    data_root, paste0(tolower(out$cohort[1L]), "_baseline_covariates.rds")
  ))
  out <- merge(out, baseline, by = c("id", "cohort"), all.x = TRUE, sort = FALSE)
  setorder(out, id, tstart, tstop)
  out[]
}

state_value_label <- function(x) {
  fifelse(is.na(x), "Missing", paste0("S", as.integer(x)))
}

make_heart_state_audit <- function(panel, scenario, state_column, bound_assumption) {
  required <- c(
    "id", "state", state_column, "analysis_sample", "exposure_entry_time",
    "event_time", "dead"
  )
  missing_columns <- setdiff(required, names(panel))
  if (length(missing_columns)) {
    stop("Heart-definition state audit inputs missing: ", paste(missing_columns, collapse = ", "))
  }

  x <- copy(panel)
  x[, `:=`(
    primary_state = as.integer(state),
    sensitivity_state = as.integer(get(state_column))
  )]
  analysis_followup <- x$analysis_sample %in% 1L &
    is.finite(x$exposure_entry_time) & is.finite(x$event_time) &
    x$event_time >= x$exposure_entry_time
  scoped <- rbindlist(list(
    copy(x)[, analysis_scope := "all_panel_rows"],
    copy(x[analysis_followup])[, analysis_scope := "analysis_followup_all_states"],
    copy(x[analysis_followup & dead %in% 0L])[, analysis_scope := "analysis_followup_living_states"]
  ), use.names = TRUE, fill = TRUE)
  scoped[, `:=`(
    scenario = scenario,
    state_column = state_column,
    bound_assumption = bound_assumption,
    both_observed = as.integer(!is.na(primary_state) & !is.na(sensitivity_state)),
    observed_state_changed = as.integer(
      !is.na(primary_state) & !is.na(sensitivity_state) &
        primary_state != sensitivity_state
    ),
    identical_including_joint_missing = as.integer(
      (is.na(primary_state) & is.na(sensitivity_state)) |
        (!is.na(primary_state) & !is.na(sensitivity_state) &
           primary_state == sensitivity_state)
    )
  )]

  summary <- scoped[, .(
    rows = .N,
    participants = uniqueN(id),
    primary_state_observed_rows = sum(!is.na(primary_state)),
    sensitivity_state_observed_rows = sum(!is.na(sensitivity_state)),
    both_state_observed_rows = sum(both_observed),
    identical_rows_including_joint_missing = sum(identical_including_joint_missing),
    changed_observed_state_rows = sum(observed_state_changed),
    participants_with_changed_observed_state = uniqueN(id[observed_state_changed == 1L]),
    primary_missing_sensitivity_observed_rows = sum(is.na(primary_state) & !is.na(sensitivity_state)),
    primary_observed_sensitivity_missing_rows = sum(!is.na(primary_state) & is.na(sensitivity_state)),
    joint_missing_rows = sum(is.na(primary_state) & is.na(sensitivity_state))
  ), by = .(scenario, state_column, bound_assumption, analysis_scope)]

  crosswalk <- scoped[, .(
    rows = .N,
    participants = uniqueN(id)
  ), by = .(
    scenario, state_column, bound_assumption, analysis_scope,
    primary_state, sensitivity_state
  )]
  crosswalk[, `:=`(
    primary_state_label = state_value_label(primary_state),
    sensitivity_state_label = state_value_label(sensitivity_state),
    observed_state_changed = as.integer(
      !is.na(primary_state) & !is.na(sensitivity_state) &
        primary_state != sensitivity_state
    )
  )]
  setcolorder(crosswalk, c(
    "scenario", "state_column", "bound_assumption", "analysis_scope",
    "primary_state", "primary_state_label", "sensitivity_state",
    "sensitivity_state_label", "observed_state_changed", "rows", "participants"
  ))
  setorder(crosswalk, scenario, analysis_scope, primary_state, sensitivity_state, na.last = TRUE)
  list(summary = summary, crosswalk = crosswalk)
}

make_heart_transition_audit <- function(primary_intervals, sensitivity_intervals,
                                        scenario, state_column, bound_assumption) {
  audit_edges <- rbindlist(list(
    copy(expected_edges)[, edge_role := "frozen_model_intensity_edge"],
    data.table(
      transition_id = "PANEL_S0_S3", from_state = 0L, to_state = 3L,
      transition = "S0->S3", death_transition = 0L,
      confirmatory_nonfatal = 0L,
      edge_role = "panel_observed_interval_censored_multistep_jump"
    )
  ), use.names = TRUE, fill = TRUE)
  audit_edges[, followup_scope := fifelse(
    death_transition == 1L,
    paste0("ELSA mortality risk truncated at ", elsa_death_cutoff),
    "ELSA nonfatal follow-up through W10"
  )]

  summarise_edges <- function(intervals, prefix) {
    out <- rbindlist(lapply(seq_len(nrow(audit_edges)), function(i) {
      edge <- audit_edges[i]
      risk_all <- copy(intervals[from_state == edge$from_state])
      risk_all[, calendar_interval_end := exposure_entry_time + tstop]
      post_cutoff_excluded <- 0L
      if (edge$death_transition == 1L) {
        post_cutoff_excluded <- risk_all[event_time >= elsa_death_cutoff, .N]
        risk <- risk_all[event_time < elsa_death_cutoff]
        risk[, audit_interval_years := pmax(
          0, pmin(calendar_interval_end, elsa_death_cutoff) - event_time
        )]
        risk <- risk[audit_interval_years > 0]
        events <- risk[
          state_changed == 1L & to_state == edge$to_state &
            calendar_interval_end <= elsa_death_cutoff
        ]
      } else {
        risk <- risk_all
        risk[, audit_interval_years := interval_years]
        events <- risk[state_changed == 1L & to_state == edge$to_state]
      }
      cbind(edge, data.table(
        risk_intervals = nrow(risk),
        risk_participants = uniqueN(risk$id),
        target_events = nrow(events),
        event_participants = uniqueN(events$id),
        risk_person_years = sum(risk$audit_interval_years, na.rm = TRUE),
        post_cutoff_risk_intervals_excluded = post_cutoff_excluded
      ))
    }), use.names = TRUE, fill = TRUE)
    metric_columns <- c(
      "risk_intervals", "risk_participants", "target_events",
      "event_participants", "risk_person_years", "post_cutoff_risk_intervals_excluded"
    )
    setnames(out, metric_columns, paste0(prefix, metric_columns))
    out
  }

  primary <- summarise_edges(primary_intervals, "primary_")
  sensitivity <- summarise_edges(sensitivity_intervals, "sensitivity_")
  key_columns <- intersect(names(audit_edges), names(primary))
  out <- merge(primary, sensitivity, by = key_columns, all = TRUE, sort = FALSE)
  out[, `:=`(
    scenario = scenario,
    state_column = state_column,
    bound_assumption = bound_assumption,
    risk_interval_delta = sensitivity_risk_intervals - primary_risk_intervals,
    risk_participant_delta = sensitivity_risk_participants - primary_risk_participants,
    event_delta = sensitivity_target_events - primary_target_events,
    event_participant_delta = sensitivity_event_participants - primary_event_participants,
    risk_person_year_delta = sensitivity_risk_person_years - primary_risk_person_years,
    post_cutoff_exclusion_delta =
      sensitivity_post_cutoff_risk_intervals_excluded -
        primary_post_cutoff_risk_intervals_excluded
  )]
  setcolorder(out, c(
    "scenario", "state_column", "bound_assumption", "transition_id", "transition",
    "from_state", "to_state", "edge_role", "death_transition", "confirmatory_nonfatal",
    setdiff(names(out), c(
      "scenario", "state_column", "bound_assumption", "transition_id", "transition",
      "from_state", "to_state", "edge_role", "death_transition", "confirmatory_nonfatal"
    ))
  ))
  setorder(out, scenario, from_state, to_state)
  out[]
}

make_heart_interval_crosswalk <- function(primary_intervals, sensitivity_intervals,
                                          scenario, state_column, bound_assumption) {
  make_signature <- function(intervals, prefix) {
    out <- intervals[, .(
      id, tstart = as.numeric(tstart), tstop = as.numeric(tstop),
      from_state = as.integer(from_state), to_state = as.integer(to_state),
      present = 1L
    )]
    if (out[, anyDuplicated(paste(id, sprintf("%.10f", tstart), sprintf("%.10f", tstop)))]) {
      stop("Interval signatures are not unique by participant, tstart and tstop")
    }
    setnames(
      out,
      c("from_state", "to_state", "present"),
      paste0(prefix, c("from_state", "to_state", "present"))
    )
    out
  }

  primary <- make_signature(primary_intervals, "primary_")
  sensitivity <- make_signature(sensitivity_intervals, "sensitivity_")
  out <- merge(primary, sensitivity, by = c("id", "tstart", "tstop"), all = TRUE, sort = FALSE)
  out[, primary_transition := fcase(
    primary_present != 1L | is.na(primary_present), "<interval absent>",
    primary_from_state == primary_to_state, "No transition",
    default = paste0("S", primary_from_state, "->S", primary_to_state)
  )]
  out[, sensitivity_transition := fcase(
    sensitivity_present != 1L | is.na(sensitivity_present), "<interval absent>",
    sensitivity_from_state == sensitivity_to_state, "No transition",
    default = paste0("S", sensitivity_from_state, "->S", sensitivity_to_state)
  )]
  out[, transition_class_changed := as.integer(primary_transition != sensitivity_transition)]
  out <- out[, .(
    intervals = .N,
    participants = uniqueN(id),
    interval_person_years = sum(tstop - tstart, na.rm = TRUE)
  ), by = .(primary_transition, sensitivity_transition, transition_class_changed)]
  out[, `:=`(
    scenario = scenario,
    state_column = state_column,
    bound_assumption = bound_assumption
  )]
  setcolorder(out, c(
    "scenario", "state_column", "bound_assumption", "primary_transition",
    "sensitivity_transition", "transition_class_changed", "intervals",
    "participants", "interval_person_years"
  ))
  setorder(out, scenario, -transition_class_changed, primary_transition, sensitivity_transition)
  out[]
}

add_state_duration <- function(x) {
  x <- copy(x)
  setorder(x, id, tstart, tstop)
  x[, spell_number := cumsum(seq_len(.N) == 1L | from_state != shift(from_state)), by = id]
  x[, state_entry_tstart := min(tstart), by = .(id, spell_number)]
  x[, state_duration_years := pmax(0, tstart - state_entry_tstart)]
  x[, log1p_state_duration := log1p(state_duration_years)]
  x[]
}

make_alternate_death_intervals <- function(intervals, panel) {
  x <- copy(intervals)
  x[, primary_tstop := as.numeric(tstop)]
  deaths <- panel[dead == 1L & is.finite(event_time_original) & is.finite(exposure_entry_time), .(
    id, alternate_death_tstop = as.numeric(event_time_original - exposure_entry_time),
    original_analysis_death_time = as.numeric(event_time),
    alternate_calendar_death_time = as.numeric(event_time_original)
  )]
  deaths <- unique(deaths, by = "id")
  x[deaths, on = "id", `:=`(
    alternate_death_tstop = i.alternate_death_tstop,
    original_analysis_death_time = i.original_analysis_death_time,
    alternate_calendar_death_time = i.alternate_calendar_death_time
  )]
  x[to_state == 4L & is.finite(alternate_death_tstop), tstop := alternate_death_tstop]
  x[, interval_years := tstop - tstart]
  x[, alternate_death_time_valid := as.integer(
    to_state != 4L | (is.finite(alternate_death_tstop) & tstop > tstart)
  )]
  x[]
}

make_expanded_hypertension_intervals <- function(panel) {
  required <- c(
    "hypertension_raw", "antihypertensive_medication", "systolic_bp",
    "diastolic_bp", "heart_disease", "stroke", "dead"
  )
  if (!all(required %in% names(panel))) {
    stop("Expanded hypertension inputs missing: ", paste(setdiff(required, names(panel)), collapse = ", "))
  }
  x <- copy(panel)
  x[, measured_high_bp := fifelse(
    !is.na(systolic_bp) | !is.na(diastolic_bp),
    as.integer((!is.na(systolic_bp) & systolic_bp >= 140) |
      (!is.na(diastolic_bp) & diastolic_bp >= 90)),
    NA_integer_
  )]
  x[, expanded_hypertension_visit := fcase(
    hypertension_raw == 1L | antihypertensive_medication == 1L | measured_high_bp == 1L, 1L,
    hypertension_raw == 0L & antihypertensive_medication == 0L & measured_high_bp == 0L, 0L,
    default = NA_integer_
  )]
  setorder(x, id, event_time, dead, wave, na.last = TRUE)
  x[, expanded_hypertension := cum_diagnosis(expanded_hypertension_visit), by = id]
  x[, state_expanded_hypertension := derive_state(
    expanded_hypertension, heart_disease, stroke, dead
  )]
  list(panel = x, intervals = build_intervals_from_panel(x, "state_expanded_hypertension"))
}

make_denial_audit <- function(panel, cohort) {
  rows <- list()
  for (disease in c("hypertension_raw", "heart_disease_raw", "stroke_raw")) {
    if (!disease %in% names(panel)) next
    x <- panel[event_type == "visit", .(id, wave, event_time, response = as.integer(get(disease)))]
    setorder(x, id, event_time, wave)
    x[, prior_positive := shift(cummax(fifelse(response == 1L, 1L, 0L)), fill = 0L), by = id]
    denial <- x[prior_positive == 1L & response == 0L]
    rows[[length(rows) + 1L]] <- data.table(
      cohort = cohort,
      disease = sub("_raw$", "", disease),
      observed_visit_rows = nrow(x),
      positive_visit_rows = x[response == 1L, .N],
      later_negative_after_positive_rows = nrow(denial),
      participants_with_later_negative_after_positive = uniqueN(denial$id),
      interpretation = "Descriptive denial/inconsistency audit only; cumulative diagnosis in the primary state prevents reverse disease transitions"
    )
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

fit_ipcw <- function(panel, intervals, cohort) {
  p <- copy(panel)
  setorder(p, id, time_msm, state_msm)
  maximum_wave <- max(p$wave[p$wave < 900L], na.rm = TRUE)
  p[, `:=`(
    next_time = shift(time_msm, type = "lead"),
    retained_next = as.integer(!is.na(shift(time_msm, type = "lead")))
  ), by = id]
  p[, retention_risk := as.integer(state_msm != 5L & wave < maximum_wave)]
  p[, `:=`(
    age_retention_5 = (as.numeric(age_tv) - 65) / 5,
    state_retention = factor(state_msm),
    wave_retention = factor(wave),
    female_retention = factor(female_bl),
    race_retention = factor(race_bl),
    education_retention = factor(education_level_bl),
    wealth_retention = factor(wealth_quintile_bl),
    smoking_retention = factor(current_smoker_bl),
    diabetes_retention = factor(diabetes_bl)
  )]
  denominator_terms <- c(
    "wave_retention", "age_retention_5", "female_retention", "state_retention",
    "insomnia_score_primary", "race_retention", "education_retention",
    "wealth_retention", "smoking_retention", "diabetes_retention"
  )
  if (cohort == "HRS" && "proxy_interview" %in% names(p) && any(!is.na(p$proxy_interview))) {
    p[, proxy_retention := factor(proxy_interview)]
    denominator_terms <- c(denominator_terms, "proxy_retention")
  }
  required <- unique(c("retained_next", denominator_terms))
  risk <- p[retention_risk == 1L & complete.cases(p[, ..required])]
  if (nrow(risk) < 100L || uniqueN(risk$retained_next) < 2L) {
    return(list(
      intervals = copy(intervals),
      diagnostics = data.table(
        cohort = cohort, status = "not_estimable", reason = "Retention model lacks at least 100 complete risk rows or both retention outcomes",
        retention_risk_rows = nrow(risk), retained = sum(risk$retained_next == 1L),
        censored = sum(risk$retained_next == 0L)
      )
    ))
  }
  denominator_formula <- reformulate(denominator_terms, response = "retained_next")
  numerator_formula <- retained_next ~ wave_retention
  denominator_fit <- tryCatch(glm(denominator_formula, data = risk, family = binomial()), error = function(e) e)
  numerator_fit <- tryCatch(glm(numerator_formula, data = risk, family = binomial()), error = function(e) e)
  if (inherits(denominator_fit, "error") || inherits(numerator_fit, "error")) {
    reason <- paste(
      if (inherits(denominator_fit, "error")) conditionMessage(denominator_fit) else "denominator fit available",
      if (inherits(numerator_fit, "error")) conditionMessage(numerator_fit) else "numerator fit available",
      sep = " | "
    )
    return(list(
      intervals = copy(intervals),
      diagnostics = data.table(
        cohort = cohort, status = "failed", reason = reason,
        retention_risk_rows = nrow(risk), retained = sum(risk$retained_next == 1L),
        censored = sum(risk$retained_next == 0L)
      )
    ))
  }
  risk[, `:=`(
    probability_denominator = pmin(0.99, pmax(0.01, predict(denominator_fit, type = "response"))),
    probability_numerator = pmin(0.99, pmax(0.01, predict(numerator_fit, type = "response")))
  )]
  p[risk, on = c("id", "time_msm"), `:=`(
    probability_denominator = i.probability_denominator,
    probability_numerator = i.probability_numerator
  )]
  p[, retention_ratio := fifelse(
    retention_risk == 1L & retained_next == 1L &
      is.finite(probability_denominator) & is.finite(probability_numerator),
    probability_numerator / probability_denominator,
    1
  )]
  p[, ipcw_untrimmed := shift(cumprod(retention_ratio), fill = 1), by = id]
  finite_weights <- p[state_msm != 5L & is.finite(ipcw_untrimmed), ipcw_untrimmed]
  limits <- quantile(finite_weights, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  p[, ipcw := pmin(limits[2L], pmax(limits[1L], ipcw_untrimmed))]
  out <- copy(intervals)
  weight_map <- unique(p[, .(id, tstart = time_msm, ipcw, ipcw_untrimmed)], by = c("id", "tstart"))
  out[weight_map, on = c("id", "tstart"), `:=`(
    ipcw = i.ipcw, ipcw_untrimmed = i.ipcw_untrimmed
  )]
  diagnostics <- data.table(
    cohort = cohort,
    status = "completed",
    reason = "",
    retention_risk_rows = nrow(risk),
    retained = sum(risk$retained_next == 1L),
    censored = sum(risk$retained_next == 0L),
    denominator_parameters = length(coef(denominator_fit)),
    denominator_converged = isTRUE(denominator_fit$converged),
    numerator_converged = isTRUE(numerator_fit$converged),
    intervals_with_ipcw = sum(is.finite(out$ipcw)),
    ipcw_min = min(out$ipcw, na.rm = TRUE),
    ipcw_p01 = quantile(out$ipcw, 0.01, na.rm = TRUE),
    ipcw_median = median(out$ipcw, na.rm = TRUE),
    ipcw_mean = mean(out$ipcw, na.rm = TRUE),
    ipcw_p99 = quantile(out$ipcw, 0.99, na.rm = TRUE),
    ipcw_max = max(out$ipcw, na.rm = TRUE),
    trimming = "Stabilized observation weights truncated at cohort-specific 1st and 99th percentiles",
    independent_unit = "Participant id; retention model uses repeated person-wave rows only to estimate observation probabilities"
  )
  list(intervals = out, diagnostics = diagnostics)
}

# -----------------------------------------------------------------------------
# Load and audit staged cohort data
# -----------------------------------------------------------------------------

base_intervals <- list()
panels <- list()
panel_msm <- list()
for (cohort in intersect(c("HRS", "ELSA"), requested_cohorts)) {
  prefix <- tolower(cohort)
  base_intervals[[cohort]] <- read_required_rds(file.path(data_root, paste0(prefix, "_multistate_intervals.rds")))
  panels[[cohort]] <- read_required_rds(file.path(data_root, paste0(prefix, "_exposure_panel_augmented.rds")))
  panel_msm[[cohort]] <- read_required_rds(file.path(data_root, paste0(prefix, "_msm_panel.rds")))
}

data_variants <- list()
data_audit <- list()
denial_audit <- list()
ipcw_diagnostics <- list()

for (cohort in names(base_intervals)) {
  base <- add_state_duration(base_intervals[[cohort]])
  base <- prepare_model_columns(base)
  data_variants[[paste(cohort, "base", sep = "::")]] <- base
  denial_audit[[cohort]] <- make_denial_audit(panels[[cohort]], cohort)

  score_support <- sort(unique(base$insomnia_score_primary[is.finite(base$insomnia_score_primary)]))
  proxy_available <- "proxy_interview" %in% names(base) && any(!is.na(base$proxy_interview))
  data_audit[[length(data_audit) + 1L]] <- data.table(
    cohort = cohort,
    audit_item = "base_interval_and_RCS_support",
    rows = nrow(base),
    participants = uniqueN(base$id),
    events_on_eight_edges = base[state_changed == 1L &
      paste0(from_state, "->", to_state) %chin% paste0(expected_edges$from_state, "->", expected_edges$to_state), .N],
    score_support = paste(score_support, collapse = ","),
    rcs_rank_supported = as.integer(setequal(as.numeric(score_support), 0:3)),
    proxy_flag_available = as.integer(proxy_available),
    detail = "Frozen eight-edge outcome; S0->S3 is not a fitted target"
  )

  ipcw <- fit_ipcw(panel_msm[[cohort]], base, cohort)
  ipcw$intervals <- prepare_model_columns(ipcw$intervals)
  data_variants[[paste(cohort, "ipcw", sep = "::")]] <- ipcw$intervals
  ipcw_diagnostics[[cohort]] <- ipcw$diagnostics

  if (cohort == "HRS") {
    alternate <- make_alternate_death_intervals(base, panels[[cohort]])
    alternate <- prepare_model_columns(alternate)
    data_variants[[paste(cohort, "alternate_death", sep = "::")]] <- alternate
    data_audit[[length(data_audit) + 1L]] <- data.table(
      cohort = cohort,
      audit_item = "alternate_death_timing",
      rows = alternate[to_state == 4L, .N],
      participants = alternate[to_state == 4L, uniqueN(id)],
      events_on_eight_edges = alternate[to_state == 4L & state_changed == 1L, .N],
      score_support = NA_character_,
      rcs_rank_supported = NA_integer_,
      proxy_flag_available = NA_integer_,
      detail = paste0(
        "death intervals changed from R/05 timing to harmonized death-year midpoint=",
        alternate[to_state == 4L & is.finite(alternate_death_tstop) &
          abs(alternate_death_tstop - primary_tstop) > 1e-10, .N],
        "; invalid alternate intervals=", alternate[to_state == 4L & alternate_death_time_valid != 1L, .N]
      )
    )
  }
}

expanded_state_counts <- data.table(
  cohort = character(),
  rows = integer(),
  primary_state_missing = integer(),
  expanded_state_missing = integer(),
  changed_state_rows = integer(),
  expanded_hypertension_positive = integer()
)
if ("ELSA" %in% names(panels)) {
  expanded <- tryCatch(make_expanded_hypertension_intervals(panels$ELSA), error = function(e) e)
  if (!inherits(expanded, "error")) {
    expanded_intervals <- add_state_duration(expanded$intervals)
    expanded_intervals <- prepare_model_columns(expanded_intervals)
    data_variants[["ELSA::expanded_hypertension"]] <- expanded_intervals
    expanded_state_counts <- expanded$panel[dead != 1L, .(
      rows = .N,
      primary_state_missing = sum(is.na(state)),
      expanded_state_missing = sum(is.na(state_expanded_hypertension)),
      changed_state_rows = sum(state != state_expanded_hypertension, na.rm = TRUE),
      expanded_hypertension_positive = sum(expanded_hypertension == 1L, na.rm = TRUE)
    ), by = cohort]
    data_audit[[length(data_audit) + 1L]] <- data.table(
      cohort = "ELSA",
      audit_item = "expanded_hypertension_state",
      rows = nrow(expanded_intervals),
      participants = uniqueN(expanded_intervals$id),
      events_on_eight_edges = expanded_intervals[state_changed == 1L &
        paste0(from_state, "->", to_state) %chin% paste0(expected_edges$from_state, "->", expected_edges$to_state), .N],
      score_support = NA_character_,
      rcs_rank_supported = NA_integer_,
      proxy_flag_available = NA_integer_,
      detail = paste0(
        "Alternative living-state rows changed=", expanded_state_counts$changed_state_rows,
        "; definition=self-report OR antihypertensive medication OR SBP>=140 OR DBP>=90, then cumulative"
      )
    )
  } else {
    data_audit[[length(data_audit) + 1L]] <- data.table(
      cohort = "ELSA", audit_item = "expanded_hypertension_state",
      rows = NA_integer_, participants = NA_integer_, events_on_eight_edges = NA_integer_,
      score_support = NA_character_, rcs_rank_supported = 0L, proxy_flag_available = NA_integer_,
      detail = paste("Not constructible:", conditionMessage(expanded))
    )
  }
}

heart_state_summary_rows <- list()
heart_state_crosswalk_rows <- list()
heart_transition_audit_rows <- list()
heart_interval_crosswalk_rows <- list()
if ("ELSA" %in% names(panels)) {
  for (i in seq_len(nrow(heart_definition_scenarios))) {
    definition <- heart_definition_scenarios[i]
    scenario_value <- definition$scenario[[1L]]
    state_column_value <- definition$state_column[[1L]]
    assumption_value <- definition$bound_assumption[[1L]]

    state_audit <- tryCatch(
      make_heart_state_audit(
        panels$ELSA, scenario_value, state_column_value, assumption_value
      ),
      error = function(e) e
    )
    if (!inherits(state_audit, "error")) {
      heart_state_summary_rows[[length(heart_state_summary_rows) + 1L]] <- state_audit$summary
      heart_state_crosswalk_rows[[length(heart_state_crosswalk_rows) + 1L]] <- state_audit$crosswalk
    }

    variant <- tryCatch({
      intervals <- build_intervals_from_panel(panels$ELSA, state_column_value)
      intervals <- add_state_duration(intervals)
      prepare_model_columns(intervals)
    }, error = function(e) e)

    if (inherits(state_audit, "error") || inherits(variant, "error")) {
      reason <- paste(
        if (inherits(state_audit, "error")) conditionMessage(state_audit) else "State audit available",
        if (inherits(variant, "error")) conditionMessage(variant) else "Interval variant available",
        sep = " | "
      )
      data_audit[[length(data_audit) + 1L]] <- data.table(
        cohort = "ELSA",
        audit_item = paste0(scenario_value, "_heart_definition"),
        rows = NA_integer_, participants = NA_integer_, events_on_eight_edges = NA_integer_,
        score_support = NA_character_, rcs_rank_supported = NA_integer_,
        proxy_flag_available = NA_integer_,
        detail = paste("Not constructible:", reason)
      )
      next
    }

    data_variants[[paste("ELSA", scenario_value, sep = "::")]] <- variant
    transition_audit <- make_heart_transition_audit(
      data_variants[["ELSA::base"]], variant, scenario_value,
      state_column_value, assumption_value
    )
    interval_crosswalk <- make_heart_interval_crosswalk(
      data_variants[["ELSA::base"]], variant, scenario_value,
      state_column_value, assumption_value
    )
    heart_transition_audit_rows[[length(heart_transition_audit_rows) + 1L]] <- transition_audit
    heart_interval_crosswalk_rows[[length(heart_interval_crosswalk_rows) + 1L]] <- interval_crosswalk

    followup_summary <- state_audit$summary[
      analysis_scope == "analysis_followup_all_states"
    ]
    data_audit[[length(data_audit) + 1L]] <- data.table(
      cohort = "ELSA",
      audit_item = paste0(scenario_value, "_heart_definition"),
      rows = nrow(variant),
      participants = uniqueN(variant$id),
      events_on_eight_edges = variant[
        state_changed == 1L &
          paste0(from_state, "->", to_state) %chin%
            paste0(expected_edges$from_state, "->", expected_edges$to_state),
        .N
      ],
      score_support = paste(
        sort(unique(variant$insomnia_score_primary[is.finite(variant$insomnia_score_primary)])),
        collapse = ","
      ),
      rcs_rank_supported = NA_integer_,
      proxy_flag_available = NA_integer_,
      detail = paste0(
        assumption_value,
        "; follow-up state rows changed versus primary=",
        followup_summary$changed_observed_state_rows,
        "; participants changed=",
        followup_summary$participants_with_changed_observed_state,
        "; absolute event-count delta across eight edges=",
        transition_audit[edge_role == "frozen_model_intensity_edge", sum(abs(event_delta))]
      )
    )
  }
}

heart_state_summary <- if (length(heart_state_summary_rows)) {
  rbindlist(heart_state_summary_rows, use.names = TRUE, fill = TRUE)
} else {
  data.table(
    scenario = character(), state_column = character(), bound_assumption = character(),
    analysis_scope = character(), rows = integer(), participants = integer(),
    primary_state_observed_rows = integer(), sensitivity_state_observed_rows = integer(),
    both_state_observed_rows = integer(), identical_rows_including_joint_missing = integer(),
    changed_observed_state_rows = integer(), participants_with_changed_observed_state = integer(),
    primary_missing_sensitivity_observed_rows = integer(),
    primary_observed_sensitivity_missing_rows = integer(), joint_missing_rows = integer()
  )
}
heart_state_crosswalk <- if (length(heart_state_crosswalk_rows)) {
  rbindlist(heart_state_crosswalk_rows, use.names = TRUE, fill = TRUE)
} else {
  data.table(
    scenario = character(), state_column = character(), bound_assumption = character(),
    analysis_scope = character(), primary_state = integer(), primary_state_label = character(),
    sensitivity_state = integer(), sensitivity_state_label = character(),
    observed_state_changed = integer(), rows = integer(), participants = integer()
  )
}
heart_transition_audit <- if (length(heart_transition_audit_rows)) {
  rbindlist(heart_transition_audit_rows, use.names = TRUE, fill = TRUE)
} else {
  data.table(
    scenario = character(), state_column = character(), bound_assumption = character(),
    transition_id = character(), transition = character(), from_state = integer(),
    to_state = integer(), edge_role = character(), death_transition = integer(),
    confirmatory_nonfatal = integer(), followup_scope = character(),
    primary_risk_intervals = integer(),
    primary_risk_participants = integer(), primary_target_events = integer(),
    primary_event_participants = integer(), primary_risk_person_years = numeric(),
    primary_post_cutoff_risk_intervals_excluded = integer(),
    sensitivity_risk_intervals = integer(), sensitivity_risk_participants = integer(),
    sensitivity_target_events = integer(), sensitivity_event_participants = integer(),
    sensitivity_risk_person_years = numeric(),
    sensitivity_post_cutoff_risk_intervals_excluded = integer(),
    risk_interval_delta = integer(),
    risk_participant_delta = integer(), event_delta = integer(),
    event_participant_delta = integer(), risk_person_year_delta = numeric(),
    post_cutoff_exclusion_delta = integer()
  )
}
heart_interval_crosswalk <- if (length(heart_interval_crosswalk_rows)) {
  rbindlist(heart_interval_crosswalk_rows, use.names = TRUE, fill = TRUE)
} else {
  data.table(
    scenario = character(), state_column = character(), bound_assumption = character(),
    primary_transition = character(), sensitivity_transition = character(),
    transition_class_changed = integer(), intervals = integer(), participants = integer(),
    interval_person_years = numeric()
  )
}

data_audit <- rbindlist(data_audit, use.names = TRUE, fill = TRUE)
denial_audit <- rbindlist(denial_audit, use.names = TRUE, fill = TRUE)
ipcw_diagnostics <- rbindlist(ipcw_diagnostics, use.names = TRUE, fill = TRUE)

# -----------------------------------------------------------------------------
# Prespecified supplementary model catalog
# -----------------------------------------------------------------------------

scenario_definitions <- rbindlist(list(
  data.table(
    scenario = c(
      "rcs", "alternate_death", "attained_age_period", "ipcw",
      "proxy_exclusion", "expanded_hypertension", "semi_markov"
    ),
    sensitivity_id = c(1L, 10L, 11L, 13L, 14L, 6L, 12L),
    label = c(
      "Restricted cubic spline/nonlinearity for insomnia score",
      "HRS harmonized death-year midpoint instead of last interview plus months-to-death",
      "Attained-age time scale with calendar-period-stratified baseline hazard",
      "Stabilized inverse-probability-of-censoring weighting",
      "Exclusion of recorded proxy interviews",
      "Expanded hypertension using diagnosis, medication or measured blood pressure",
      "Semi-Markov Cox approximation with time since current-state entry"
    ),
    data_variant = c(
      "base", "alternate_death", "base", "ipcw", "base",
      "expanded_hypertension", "base"
    )
  ),
  heart_definition_scenarios[, .(
    scenario, sensitivity_id, label, data_variant = scenario
  )]
), use.names = TRUE, fill = TRUE)

make_specs <- function() {
  rows <- list()
  for (i in seq_len(nrow(scenario_definitions))) {
    scenario <- scenario_definitions[i]
    cohorts <- switch(
      scenario$scenario,
      alternate_death = "HRS",
      proxy_exclusion = "HRS",
      expanded_hypertension = "ELSA",
      heart_lower = "ELSA",
      heart_upper = "ELSA",
      c("HRS", "ELSA")
    )
    cohorts <- intersect(cohorts, names(base_intervals))
    edges <- if (scenario$scenario == "alternate_death") expected_edges[to_state == 4L] else expected_edges
    for (cohort in cohorts) {
      value <- copy(edges)
      value[, `:=`(
        cohort = cohort,
        scenario = scenario$scenario,
        sensitivity_id = scenario$sensitivity_id,
        scenario_label = scenario$label,
        data_variant = scenario$data_variant
      )]
      rows[[length(rows) + 1L]] <- value
    }
  }
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  out <- out[scenario %chin% requested_scenarios]
  out[, spec_id := sprintf("PS%02d_%s_%s_%s", sensitivity_id, cohort, scenario, transition_id)]
  out[, `:=`(
    outcome_state_column = fcase(
      scenario == "heart_lower", "state_heart_lower",
      scenario == "heart_upper", "state_heart_upper",
      default = "state"
    ),
    outcome_definition = fcase(
      scenario == "heart_lower", heart_definition_scenarios[
        scenario == "heart_lower", bound_assumption
      ],
      scenario == "heart_upper", heart_definition_scenarios[
        scenario == "heart_upper", bound_assumption
      ],
      default = "Primary disease-state definition"
    ),
    analysis_role = "Protocol sensitivity; primary inference remains the panel-observed continuous-time Markov model",
    model_level = "M4",
    model_adjustment = paste(
      "Age, sex, race/ethnicity, education, wealth quintile, marital status,",
      "smoking, alcohol, physical inactivity, BMI, diabetes, lung disease and pain"
    ),
    independent_unit = "Participant id with cluster-robust covariance across repeated risk intervals",
    mortality_scope = fifelse(
      cohort == "ELSA" & death_transition == 1L,
      paste0(
        "ELSA limited mortality validation; death-risk follow-up truncated at ",
        elsa_death_cutoff
      ),
      fifelse(cohort == "ELSA", "ELSA W4-W10 nonfatal follow-up", "HRS follow-up through 2020")
    )
  )]
  out[]
}

specifications <- make_specs()

scenario_model_terms <- function(scenario) {
  exposure <- if (scenario == "rcs") c("rcs_linear", "rcs_nl1", "rcs_nl2") else "insomnia_score_primary"
  covariates <- if (scenario == "attained_age_period") m4_terms_no_age else m4_terms
  if (scenario == "semi_markov") covariates <- c(covariates, "log1p_state_duration")
  list(exposure = exposure, covariates = covariates)
}

prepare_spec_data <- function(spec) {
  key <- paste(spec$cohort, spec$data_variant, sep = "::")
  x <- data_variants[[key]]
  if (is.null(x)) {
    return(list(ok = FALSE, reason = paste("Data variant unavailable:", key), data = data.table()))
  }
  d <- copy(x[from_state == spec$from_state])
  if (spec$scenario == "proxy_exclusion") {
    if (!"proxy_interview" %in% names(d) || !any(!is.na(d$proxy_interview))) {
      return(list(ok = FALSE, reason = "No populated proxy-interview flag is available", data = data.table()))
    }
    d <- d[proxy_interview == 0L]
  }
  if (spec$scenario == "alternate_death") {
    d <- d[alternate_death_time_valid == 1L]
  }
  if (spec$scenario == "ipcw") {
    if (!"ipcw" %in% names(d) || !any(is.finite(d$ipcw) & d$ipcw > 0)) {
      return(list(ok = FALSE, reason = "Stabilized IPCW was not estimable for this cohort", data = data.table()))
    }
    d[, model_weight := as.numeric(ipcw)]
  } else {
    d[, model_weight := 1]
  }
  d[, calendar_interval_end := exposure_entry_time + tstop]
  mortality_cutoff_applied <- NA_real_
  post_cutoff_intervals_excluded <- 0L
  intervals_truncated_at_cutoff <- 0L
  if (spec$cohort == "ELSA" && spec$death_transition == 1L) {
    mortality_cutoff_applied <- elsa_death_cutoff
    post_cutoff_intervals_excluded <- d[event_time >= elsa_death_cutoff, .N]
    d <- d[event_time < elsa_death_cutoff]
    intervals_truncated_at_cutoff <- d[calendar_interval_end > elsa_death_cutoff, .N]
    d[, outcome_tstop := pmin(calendar_interval_end, elsa_death_cutoff) - exposure_entry_time]
    d[, event_target := as.integer(
      state_changed == 1L & to_state == spec$to_state &
        calendar_interval_end <= elsa_death_cutoff
    )]
  } else {
    d[, outcome_tstop := tstop]
    d[, event_target := as.integer(state_changed == 1L & to_state == spec$to_state)]
  }
  if (spec$scenario == "attained_age_period") {
    d[, `:=`(
      analysis_tstart = age_bl + tstart,
      analysis_tstop = age_bl + outcome_tstop
    )]
  } else {
    d[, `:=`(analysis_tstart = tstart, analysis_tstop = outcome_tstop)]
  }
  terms <- scenario_model_terms(spec$scenario)
  required <- unique(c(
    "id", "analysis_tstart", "analysis_tstop", "event_target", "model_weight",
    terms$exposure, terms$covariates,
    if (spec$scenario == "attained_age_period") "calendar_period" else character()
  ))
  missing_columns <- setdiff(required, names(d))
  if (length(missing_columns)) {
    return(list(ok = FALSE, reason = paste("Required columns absent:", paste(missing_columns, collapse = ", ")), data = data.table()))
  }
  d <- d[
    complete.cases(d[, ..required]) & is.finite(analysis_tstart) &
      is.finite(analysis_tstop) & analysis_tstop > analysis_tstart &
      is.finite(model_weight) & model_weight > 0
  ]
  factor_columns <- names(d)[vapply(d, is.factor, logical(1L))]
  if (length(factor_columns)) d[, (factor_columns) := lapply(.SD, droplevels), .SDcols = factor_columns]
  constant_terms <- c(terms$covariates[vapply(terms$covariates, function(name) {
    uniqueN(d[[name]], na.rm = TRUE) < 2L
  }, logical(1L))])
  terms$covariates <- setdiff(terms$covariates, constant_terms)
  if (!nrow(d) || uniqueN(d$id) < 2L) {
    return(list(ok = FALSE, reason = "Fewer than two participants remain after complete-case filtering", data = d))
  }
  if (spec$scenario == "rcs" && uniqueN(d$insomnia_score_primary) < 4L) {
    return(list(ok = FALSE, reason = "All four score values 0-3 are required for the prespecified RCS", data = d))
  }
  regular_terms <- c(terms$exposure, terms$covariates)
  matrix_result <- tryCatch(model.matrix(reformulate(regular_terms), data = d), error = function(e) e)
  if (inherits(matrix_result, "error")) {
    return(list(ok = FALSE, reason = paste("Design matrix failed:", conditionMessage(matrix_result)), data = d))
  }
  matrix_result <- matrix_result[, colnames(matrix_result) != "(Intercept)", drop = FALSE]
  n_parameters <- qr(matrix_result)$rank
  n_events_value <- sum(d$event_target == 1L)
  epp <- if (n_parameters > 0L) n_events_value / n_parameters else NA_real_
  reasons <- character()
  if (n_events_value < min_events) reasons <- c(reasons, paste0("events ", n_events_value, " < ", min_events))
  if (!is.finite(epp) || epp < min_epp) reasons <- c(reasons, paste0("events per parameter ", format(epp, digits = 3), " < ", min_epp))
  if (!n_parameters) reasons <- c(reasons, "zero-rank design matrix")
  list(
    ok = !length(reasons),
    reason = paste(reasons, collapse = "; "),
    data = d,
    exposure_terms = terms$exposure,
    covariate_terms = terms$covariates,
    omitted_constant_terms = constant_terms,
    n_records = nrow(d),
    n_subjects = uniqueN(d$id),
    n_events = n_events_value,
    n_parameters = n_parameters,
    events_per_parameter = epp,
    mortality_cutoff_applied = mortality_cutoff_applied,
    post_cutoff_intervals_excluded = post_cutoff_intervals_excluded,
    intervals_truncated_at_cutoff = intervals_truncated_at_cutoff
  )
}

fit_one_spec <- function(spec, prepared) {
  d <- prepared$data
  regular_terms <- c(prepared$exposure_terms, prepared$covariate_terms)
  fit_terms <- c(
    regular_terms,
    if (spec$scenario == "attained_age_period") "strata(calendar_period)" else character(),
    "cluster(id)"
  )
  formula <- as.formula(paste(
    "Surv(analysis_tstart, analysis_tstop, event_target) ~",
    paste(fit_terms, collapse = " + ")
  ))
  started <- proc.time()[[3L]]
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      coxph(
        formula, data = d, weights = model_weight, ties = "efron",
        robust = TRUE, singular.ok = FALSE,
        control = coxph.control(iter.max = 40, eps = 1e-8)
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[[3L]] - started
  if (inherits(fit, "error")) {
    return(list(
      fit = NULL, estimates = data.table(), joint = data.table(),
      status = "failed", reason = conditionMessage(fit), warnings = warnings,
      elapsed_seconds = elapsed, ph_global_p_value = NA_real_
    ))
  }
  beta <- coef(fit)
  variance <- fit$var
  if (!length(beta) || any(!is.finite(beta)) || !all(is.finite(variance))) {
    return(list(
      fit = fit, estimates = data.table(), joint = data.table(),
      status = "failed", reason = "Cox fit returned non-finite coefficients or robust covariance",
      warnings = warnings, elapsed_seconds = elapsed, ph_global_p_value = NA_real_
    ))
  }
  ph_global <- tryCatch({
    z <- cox.zph(fit, transform = "km")$table
    if ("GLOBAL" %in% rownames(z)) as.numeric(z["GLOBAL", "p"]) else NA_real_
  }, error = function(e) NA_real_)
  base_metadata <- data.table(
    spec_id = spec$spec_id,
    sensitivity_id = spec$sensitivity_id,
    cohort = spec$cohort,
    scenario = spec$scenario,
    scenario_label = spec$scenario_label,
    outcome_state_column = spec$outcome_state_column,
    outcome_definition = spec$outcome_definition,
    transition_id = spec$transition_id,
    transition = spec$transition,
    from_state = spec$from_state,
    to_state = spec$to_state,
    death_transition = spec$death_transition,
    confirmatory_nonfatal = spec$confirmatory_nonfatal,
    mortality_scope = spec$mortality_scope,
    mortality_cutoff_applied = prepared$mortality_cutoff_applied,
    post_cutoff_intervals_excluded = prepared$post_cutoff_intervals_excluded,
    intervals_truncated_at_cutoff = prepared$intervals_truncated_at_cutoff,
    model_level = spec$model_level,
    model_adjustment = spec$model_adjustment,
    n_records = prepared$n_records,
    n_subjects = prepared$n_subjects,
    events = prepared$n_events,
    n_parameters = prepared$n_parameters,
    events_per_parameter = prepared$events_per_parameter,
    independent_unit = spec$independent_unit,
    status = "completed",
    reason = ""
  )
  estimate_rows <- list()
  joint_rows <- list()
  if (spec$scenario == "rcs") {
    exposure_names <- c("rcs_linear", "rcs_nl1", "rcs_nl2")
    index <- match(exposure_names, names(beta))
    if (anyNA(index)) stop("RCS coefficients were not retained in the fitted model")
    basis_values <- cbind(
      score = 0:3,
      as.data.frame(rcs_basis(0:3))
    )
    basis_matrix <- as.matrix(basis_values[, exposure_names])
    for (score in 1:3) {
      contrast <- basis_matrix[score + 1L, ] - basis_matrix[1L, ]
      log_hr <- sum(contrast * beta[index])
      se <- sqrt(as.numeric(t(contrast) %*% variance[index, index, drop = FALSE] %*% contrast))
      z_value <- log_hr / se
      estimate_rows[[length(estimate_rows) + 1L]] <- cbind(base_metadata, data.table(
        term = paste0("insomnia_score_", score, "_vs_0"),
        estimate_type = "RCS contrast",
        comparison = paste0(score, " symptoms versus 0 symptoms"),
        log_hazard_ratio = log_hr,
        robust_standard_error = se,
        hazard_ratio = exp(log_hr),
        conf_low = exp(log_hr - qnorm(0.975) * se),
        conf_high = exp(log_hr + qnorm(0.975) * se),
        z_value = z_value,
        p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE),
        rcs_knots = "0,1,2,3"
      ))
    }
    nonlinear_index <- match(c("rcs_nl1", "rcs_nl2"), names(beta))
    nonlinear <- wald_quadratic(beta[nonlinear_index], variance[nonlinear_index, nonlinear_index, drop = FALSE])
    overall <- wald_quadratic(beta[index], variance[index, index, drop = FALSE])
    joint_rows[[1L]] <- cbind(base_metadata, data.table(
      test = c("RCS nonlinearity", "RCS overall association"),
      statistic = c(nonlinear$statistic, overall$statistic),
      df = c(nonlinear$df, overall$df),
      p_value = c(nonlinear$p_value, overall$p_value),
      method = "Robust Wald chi-square test",
      tested_terms = c("rcs_nl1 | rcs_nl2", "rcs_linear | rcs_nl1 | rcs_nl2")
    ))
  } else {
    index <- match("insomnia_score_primary", names(beta))
    if (is.na(index)) stop("Insomnia score coefficient was not retained in the fitted model")
    log_hr <- beta[index]
    se <- sqrt(variance[index, index])
    z_value <- log_hr / se
    estimate_rows[[1L]] <- cbind(base_metadata, data.table(
      term = "insomnia_score_primary",
      estimate_type = "linear transition hazard ratio",
      comparison = "per one additional insomnia symptom",
      log_hazard_ratio = log_hr,
      robust_standard_error = se,
      hazard_ratio = exp(log_hr),
      conf_low = exp(log_hr - qnorm(0.975) * se),
      conf_high = exp(log_hr + qnorm(0.975) * se),
      z_value = z_value,
      p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE),
      rcs_knots = NA_character_
    ))
  }
  list(
    fit = fit,
    estimates = rbindlist(estimate_rows, use.names = TRUE, fill = TRUE),
    joint = rbindlist(joint_rows, use.names = TRUE, fill = TRUE),
    status = "completed",
    reason = "",
    warnings = warnings,
    elapsed_seconds = elapsed,
    ph_global_p_value = ph_global
  )
}

# -----------------------------------------------------------------------------
# Fit models (or complete a configuration/data audit only)
# -----------------------------------------------------------------------------

model_registry_rows <- list()
estimate_rows <- list()
joint_rows <- list()
models_started <- 0L

for (i in seq_len(nrow(specifications))) {
  spec <- specifications[i]
  prepared <- prepare_spec_data(spec)
  base_registry <- data.table(
    script_version = script_version,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    spec_id = spec$spec_id,
    sensitivity_id = spec$sensitivity_id,
    cohort = spec$cohort,
    scenario = spec$scenario,
    scenario_label = spec$scenario_label,
    outcome_state_column = spec$outcome_state_column,
    outcome_definition = spec$outcome_definition,
    transition_id = spec$transition_id,
    transition = spec$transition,
    from_state = spec$from_state,
    to_state = spec$to_state,
    mortality_scope = spec$mortality_scope,
    mortality_cutoff_applied = if (!is.null(prepared$mortality_cutoff_applied)) prepared$mortality_cutoff_applied else NA_real_,
    post_cutoff_intervals_excluded = if (!is.null(prepared$post_cutoff_intervals_excluded)) prepared$post_cutoff_intervals_excluded else NA_integer_,
    intervals_truncated_at_cutoff = if (!is.null(prepared$intervals_truncated_at_cutoff)) prepared$intervals_truncated_at_cutoff else NA_integer_,
    model_level = spec$model_level,
    model_adjustment = spec$model_adjustment,
    n_records = if (!is.null(prepared$n_records)) prepared$n_records else nrow(prepared$data),
    n_subjects = if (!is.null(prepared$n_subjects)) prepared$n_subjects else uniqueN(prepared$data$id),
    events = if (!is.null(prepared$n_events)) prepared$n_events else if (nrow(prepared$data)) sum(prepared$data$event_target) else NA_integer_,
    n_parameters = if (!is.null(prepared$n_parameters)) prepared$n_parameters else NA_integer_,
    events_per_parameter = if (!is.null(prepared$events_per_parameter)) prepared$events_per_parameter else NA_real_,
    independent_unit = spec$independent_unit
  )
  if (!prepared$ok) {
    model_registry_rows[[length(model_registry_rows) + 1L]] <- cbind(base_registry, data.table(
      status = "not_estimable", reason = prepared$reason,
      warnings = "", elapsed_seconds = 0, ph_global_p_value = NA_real_,
      omitted_constant_terms = paste(prepared$omitted_constant_terms, collapse = " | ")
    ))
    next
  }
  if (run_mode == "audit" || (max_models > 0L && models_started >= max_models)) {
    reason <- if (run_mode == "audit") {
      "Configuration/data audit mode: model fitting intentionally not started"
    } else {
      paste0("Light-run model limit reached (GA_PROTOCOL_SENS_MAX_MODELS=", max_models, ")")
    }
    model_registry_rows[[length(model_registry_rows) + 1L]] <- cbind(base_registry, data.table(
      status = "not_run", reason = reason,
      warnings = "", elapsed_seconds = 0, ph_global_p_value = NA_real_,
      omitted_constant_terms = paste(prepared$omitted_constant_terms, collapse = " | ")
    ))
    next
  }
  models_started <- models_started + 1L
  message(sprintf("[%d/%d] %s %s %s", i, nrow(specifications), spec$cohort, spec$scenario, spec$transition))
  fitted <- tryCatch(fit_one_spec(spec, prepared), error = function(e) list(
    fit = NULL, estimates = data.table(), joint = data.table(), status = "failed",
    reason = conditionMessage(e), warnings = character(), elapsed_seconds = 0,
    ph_global_p_value = NA_real_
  ))
  model_registry_rows[[length(model_registry_rows) + 1L]] <- cbind(base_registry, data.table(
    status = fitted$status,
    reason = fitted$reason,
    warnings = paste(unique(fitted$warnings), collapse = " | "),
    elapsed_seconds = fitted$elapsed_seconds,
    ph_global_p_value = fitted$ph_global_p_value,
    omitted_constant_terms = paste(prepared$omitted_constant_terms, collapse = " | ")
  ))
  if (nrow(fitted$estimates)) estimate_rows[[length(estimate_rows) + 1L]] <- fitted$estimates
  if (nrow(fitted$joint)) joint_rows[[length(joint_rows) + 1L]] <- fitted$joint
  if (save_models && !is.null(fitted$fit)) {
    saveRDS(
      list(
        specification = spec,
        fit = fitted$fit,
        data_summary = base_registry,
        script_version = script_version
      ),
      file.path(model_root, paste0(spec$spec_id, ".rds")),
      compress = "xz"
    )
  }
  rm(fitted)
  gc(verbose = FALSE)
}

model_registry <- rbindlist(model_registry_rows, use.names = TRUE, fill = TRUE)
estimates <- if (length(estimate_rows)) rbindlist(estimate_rows, use.names = TRUE, fill = TRUE) else data.table(
  spec_id = character(), sensitivity_id = integer(), cohort = character(), scenario = character(),
  outcome_state_column = character(), outcome_definition = character(), model_level = character(),
  transition = character(), from_state = integer(), to_state = integer(), term = character(),
  hazard_ratio = numeric(), conf_low = numeric(), conf_high = numeric(), p_value = numeric(),
  status = character(), reason = character()
)
joint_tests <- if (length(joint_rows)) rbindlist(joint_rows, use.names = TRUE, fill = TRUE) else data.table(
  spec_id = character(), sensitivity_id = integer(), cohort = character(), scenario = character(),
  transition = character(), test = character(), statistic = numeric(), df = numeric(),
  p_value = numeric(), status = character(), reason = character()
)

if (nrow(estimates)) {
  estimates[, `:=`(p_fdr_bh_sensitivity = NA_real_, fdr_family_size = NA_integer_)]
  estimates[is.finite(p_value), `:=`(
    p_fdr_bh_sensitivity = p.adjust(p_value, method = "BH"),
    fdr_family_size = .N
  ), by = .(cohort, scenario)]
}

if (nrow(model_registry)) {
  model_registry[, event_parameter_risk := fcase(
    is.finite(events_per_parameter) & events_per_parameter >= 10, "adequate: EPP >=10",
    is.finite(events_per_parameter) & events_per_parameter >= 5, "caution: EPP 5-<10",
    is.finite(events_per_parameter), "sparse: EPP <5",
    default = "not quantifiable"
  )]
}

# -----------------------------------------------------------------------------
# Exact coverage audit for the protocol's 16 sensitivity families
# -----------------------------------------------------------------------------

protocol_items <- data.table(
  sensitivity_id = 1:16,
  protocol_sensitivity = c(
    "Categorical dose-response and restricted cubic spline/nonlinearity",
    "Lagged exposure",
    "Exclude first two years",
    "Cumulative exposure and dynamic pattern",
    "HRS 0-4 extension and cohort-specific symptom domains",
    "Expanded hypertension definition",
    "Strict and broad heart-disease definitions",
    "Stroke excluding transient ischaemic attack",
    "Denial/outcome-misclassification sensitivity",
    "Alternative death timing",
    "Attained-age time scale and piecewise/calendar-period non-homogeneity",
    "Semi-Markov sensitivity",
    "Complete case, multiple imputation and IPCW",
    "Proxy-interview exclusion",
    "Sex, age, obesity, diabetes and wealth interactions",
    "Emission misclassification matrix if identifiable"
  )
)

coverage <- CJ(sensitivity_id = 1:16, cohort = intersect(c("HRS", "ELSA"), requested_cohorts))
coverage <- merge(coverage, protocol_items, by = "sensitivity_id", all.x = TRUE, sort = FALSE)
coverage[, `:=`(
  data_support = "",
  existing_implementation = "",
  supplementary_implementation = "",
  status = "not_estimable",
  reason = "",
  evidence_files = "",
  attempted_models = 0L,
  completed_models = 0L,
  not_estimable_models = 0L,
  failed_models = 0L,
  mortality_scope = fifelse(
    cohort == "ELSA",
    "Nonfatal follow-up through W10; any death result is limited mortality validation through 2012",
    "HRS follow-up through 2020"
  ),
  independent_unit = "Participant id; repeated intervals clustered within participant",
  multiplicity = "Exact two-sided p values with BH FDR within cohort-scenario for supplementary exposure estimates"
)]

cox_catalog_path <- file.path(output_root, "results", "cox_sensitivity", "transition_cox_model_specifications.csv")
cox_catalog <- safe_fread(cox_catalog_path)
upstream_patterns <- list(
  `1` = c("primary_burden_tv"),
  `2` = c("lag1_row", "previous_valid_assessment"),
  `3` = c("landmark2"),
  `4` = c("cumulative_average", "dynamic_pattern_delayed_entry"),
  `5` = c("onset_symptom", "maintenance_symptom", "nonrestorative_symptom", "hrs_extended_score4"),
  `13` = c("primary_score_tv", "primary_burden_tv"),
  `15` = c("interaction_sex", "interaction_age", "interaction_obesity", "interaction_diabetes", "interaction_wealth")
)

for (id in as.integer(names(upstream_patterns))) {
  patterns <- upstream_patterns[[as.character(id)]]
  for (cohort_value in coverage[sensitivity_id == id, cohort]) {
    available <- if (nrow(cox_catalog)) cox_catalog[
      toupper(get("cohort")) == toupper(cohort_value) &
        Reduce(`|`, lapply(patterns, function(pattern) grepl(pattern, analysis_id, fixed = TRUE)))
    ] else data.table()
    coverage[sensitivity_id == id & cohort == cohort_value, `:=`(
      data_support = paste0("Prespecified Cox catalog rows=", nrow(available)),
      existing_implementation = "R/06_fit_transition_cox_models.R",
      status = if (nrow(available)) "covered_upstream" else "not_estimable",
      reason = if (nrow(available)) {
        "Specification is present in the upstream Cox sensitivity catalog; execution status remains in its run-labelled outputs"
      } else {
        "No matching specification row was found in the staged upstream catalog"
      },
      evidence_files = if (nrow(available)) "results/cox_sensitivity/transition_cox_model_specifications.csv" else ""
    )]
  }
}

# New executable scenarios mapped to protocol items.
new_scenario_map <- list(
  `1` = "rcs",
  `6` = "expanded_hypertension",
  `7` = c("heart_lower", "heart_upper"),
  `10` = "alternate_death",
  `11` = "attained_age_period",
  `12` = "semi_markov",
  `13` = "ipcw",
  `14` = "proxy_exclusion"
)
for (id in as.integer(names(new_scenario_map))) {
  scenario_values <- new_scenario_map[[as.character(id)]]
  for (cohort_value in coverage[sensitivity_id == id, cohort]) {
    specs <- model_registry[
      scenario %chin% scenario_values & toupper(cohort) == toupper(cohort_value)
    ]
    if (!nrow(specs)) next
    attempted <- specs[status %chin% c("completed", "not_estimable", "failed"), .N]
    complete <- specs[status == "completed", .N]
    inestimable <- specs[status == "not_estimable", .N]
    failed <- specs[status == "failed", .N]
    not_run_count <- specs[status == "not_run", .N]
    coverage[sensitivity_id == id & cohort == cohort_value, `:=`(
      supplementary_implementation = "R/07b_protocol_sensitivities.R",
      attempted_models = attempted,
      completed_models = complete,
      not_estimable_models = inestimable,
      failed_models = failed,
      status = fcase(
        complete > 0L & complete == nrow(specs), "completed",
        complete > 0L, "partially_completed",
        not_run_count > 0L & failed == 0L, "implemented_pending_execution",
        failed > 0L, "failed",
        default = "not_estimable"
      ),
      reason = paste0(
        "Supplementary scenarios=", paste(scenario_values, collapse = " | "),
        "; specifications=", nrow(specs), "; completed=", complete,
        ", not_estimable=", inestimable, ", failed=", failed,
        ", not run=", not_run_count
      ),
      evidence_files = "results/protocol_sensitivities/protocol_sensitivity_model_registry.csv"
    )]
  }
}

# Item-specific support and non-identifiability boundaries.
coverage[sensitivity_id == 1L, data_support := "Insomnia score has observed support at all four integer values 0,1,2,3; categorical models are upstream and a four-knot restricted cubic spline is implemented here"]

coverage[sensitivity_id == 6L & cohort == "HRS", `:=`(
  data_support = "HRS staged panel contains physician-diagnosed hypertension only; no harmonized antihypertensive medication plus measured SBP/DBP fields",
  supplementary_implementation = "Not implementable from the staged HRS fields",
  status = "not_estimable",
  reason = "An expanded HRS hypertension state cannot be distinguished from the primary diagnosis-only state without medication or measured blood-pressure inputs",
  evidence_files = "data_multistate/hrs_exposure_panel_augmented.rds"
)]
coverage[sensitivity_id == 6L & cohort == "ELSA", data_support := "ELSA supplies hypertension_raw, antihypertensive_medication, systolic_bp and diastolic_bp; expanded status is cumulative after an observed positive component"]

coverage[sensitivity_id == 7L & cohort == "HRS", `:=`(
  data_support = "Only the harmonized composite hearte/heart_disease_raw indicator is retained longitudinally",
  status = "not_estimable",
  reason = "No longitudinal myocardial-infarction, angina and heart-failure components are available to construct distinct strict and broad HRS definitions",
  evidence_files = "data_multistate/hrs_exposure_panel_augmented.rds"
)]
heart_variant_keys <- paste("ELSA", heart_definition_scenarios$scenario, sep = "::")
if (all(heart_variant_keys %chin% names(data_variants))) {
  coverage[sensitivity_id == 7L & cohort == "ELSA", `:=`(
    data_support = paste(
      "The upstream ELSA panel retains primary, lower-bound and upper-bound heart-state trajectories.",
      "Bounds differ only where the source-backed named components cannot resolve historical",
      "coexisting other-heart trouble after component-specific Wave-10 dispute backcasting."
    ),
    supplementary_implementation = "R/07b_protocol_sensitivities.R; ELSA-only transition-specific M4 Cox models",
    evidence_files = paste(
      "data_multistate/elsa_exposure_panel_augmented.rds;",
      "qc/protocol_sensitivities/heart_definition_state_summary.csv;",
      "qc/protocol_sensitivities/heart_definition_state_crosswalk.csv;",
      "qc/protocol_sensitivities/heart_definition_transition_audit.csv;",
      "qc/protocol_sensitivities/heart_definition_interval_crosswalk.csv;",
      "results/protocol_sensitivities/protocol_sensitivity_model_registry.csv;",
      "results/protocol_sensitivities/protocol_sensitivity_estimates.csv"
    )
  )]
} else {
  coverage[sensitivity_id == 7L & cohort == "ELSA", `:=`(
    data_support = paste0(
      "Required upstream state columns/interval variants unavailable: ",
      paste(setdiff(heart_variant_keys, names(data_variants)), collapse = " | ")
    ),
    supplementary_implementation = "Implemented but input state bounds unavailable",
    status = "not_estimable",
    reason = "R/01 and R/03-R/05 must be rerun so state_heart_lower and state_heart_upper reach the staged augmented panel",
    evidence_files = "data_multistate/elsa_exposure_panel_augmented.rds"
  )]
}

coverage[sensitivity_id == 8L & cohort == "HRS", `:=`(
  data_support = "GlobalAgingR HRS stroke is labelled physician-diagnosed stroke or transient ischaemic attack as one inseparable binary field",
  status = "not_estimable",
  reason = "TIA cannot be removed because no separate longitudinal TIA indicator is present",
  evidence_files = "GlobalAgingR::hrs$stroke; data_multistate/hrs_exposure_panel_augmented.rds"
)]
coverage[sensitivity_id == 8L & cohort == "ELSA", `:=`(
  data_support = "The harmonized ELSA stroke field is a stroke diagnosis; no separate TIA field is supplied",
  status = "covered_by_primary_definition",
  reason = "The available ELSA primary stroke indicator is already stroke-specific, but no distinct TIA-inclusive definition exists for a within-cohort contrast",
  evidence_files = "GlobalAgingR::elsa$stroke; data_multistate/elsa_exposure_panel_augmented.rds"
)]

coverage[sensitivity_id == 9L, `:=`(
  data_support = "Repeated raw diagnosis responses permit a descriptive count of later negative reports after a prior positive, but there is no adjudicated disease gold standard or validation subset",
  supplementary_implementation = "Descriptive denial-pattern audit only",
  status = "not_estimable",
  reason = "Outcome sensitivity/specificity and false-positive/false-negative mechanisms are not identified; arbitrary relabelling would fabricate assumptions",
  evidence_files = "results/protocol_sensitivities/diagnosis_denial_pattern_audit.csv"
)]

coverage[sensitivity_id == 10L & cohort == "ELSA", `:=`(
  data_support = "Available ELSA linkage supplies year of death through 2012 and the analysis uses a year midpoint; no month/day or interval-to-death variable is retained",
  supplementary_implementation = "Not implementable as a distinct ELSA timing model",
  status = "not_estimable",
  reason = "There is no second empirically supported ELSA death-time definition to compare with the year-midpoint analysis",
  evidence_files = "data_multistate/elsa_exposure_panel_augmented.rds"
)]
coverage[sensitivity_id == 10L & cohort == "HRS", data_support := "HRS supplies both the R/05 last-interview plus RADTIMTDTH analysis time and the retained harmonized death-year midpoint"]

coverage[sensitivity_id == 11L, data_support := "Attained age is reconstructed as baseline age plus elapsed interval time; calendar period is stratified as before 2012, 2012-2015 and 2016 onward. This is a Cox non-homogeneity sensitivity, not a fitted msm pci model."]
coverage[sensitivity_id == 12L, data_support := "Time since entry into the current observed state is reconstructed for every interval and entered as log1p duration in a transition-specific Cox model"]
coverage[sensitivity_id == 13L, data_support := "Complete-case and MI specifications are upstream; stabilized observation IPCW is estimated here from participant-wave retention and trimmed at the 1st/99th percentiles"]

coverage[sensitivity_id == 14L & cohort == "ELSA", `:=`(
  data_support = "The staged ELSA panel has no populated proxy-interview flag",
  supplementary_implementation = "Not implementable from current ELSA fields",
  status = "not_estimable",
  reason = "Proxy interviews cannot be distinguished from self interviews; complete sleep-assessment restrictions are available upstream but are not a proxy flag",
  evidence_files = "data_multistate/elsa_multistate_intervals.rds"
)]
coverage[sensitivity_id == 14L & cohort == "HRS", data_support := "HRS has a populated interval-start proxy_interview indicator; recorded proxy rows are excluded"]

coverage[sensitivity_id == 16L, `:=`(
  data_support = "Observed states are deterministically constructed from self-reported cumulative diagnoses; no validation replicate, adjudicated state, or prespecified emission probabilities are available",
  supplementary_implementation = "No ematrix fitted",
  status = "not_estimable",
  reason = "A free emission matrix is not separately identifiable from transition intensities in these monotone panel states, while fixing an ematrix would require unsupported sensitivity/specificity assumptions",
  evidence_files = "data_multistate/*_msm_panel.rds"
)]

setorder(coverage, sensitivity_id, cohort)

# -----------------------------------------------------------------------------
# Write auditable outputs
# -----------------------------------------------------------------------------

atomic_fwrite(specifications, file.path(result_root, "protocol_sensitivity_specifications.csv"))
atomic_fwrite(model_registry, file.path(result_root, "protocol_sensitivity_model_registry.csv"))
atomic_fwrite(estimates, file.path(result_root, "protocol_sensitivity_estimates.csv"))
atomic_fwrite(joint_tests, file.path(result_root, "protocol_sensitivity_joint_tests.csv"))
atomic_fwrite(data_audit, file.path(qc_root, "protocol_sensitivity_data_audit.csv"))
atomic_fwrite(ipcw_diagnostics, file.path(result_root, "ipcw_weight_diagnostics.csv"))
atomic_fwrite(denial_audit, file.path(result_root, "diagnosis_denial_pattern_audit.csv"))
atomic_fwrite(expanded_state_counts, file.path(result_root, "expanded_hypertension_state_audit.csv"))
atomic_fwrite(heart_state_summary, file.path(qc_root, "heart_definition_state_summary.csv"))
atomic_fwrite(heart_state_crosswalk, file.path(qc_root, "heart_definition_state_crosswalk.csv"))
atomic_fwrite(heart_transition_audit, file.path(qc_root, "heart_definition_transition_audit.csv"))
atomic_fwrite(heart_interval_crosswalk, file.path(qc_root, "heart_definition_interval_crosswalk.csv"))
atomic_fwrite(coverage, file.path(report_root, "protocol_sensitivity_coverage.csv"))

session_lines <- c(
  paste("script_version:", script_version),
  paste("generated_at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste("run_mode:", run_mode),
  paste("cohorts:", paste(requested_cohorts, collapse = ",")),
  paste("scenarios:", paste(requested_scenarios, collapse = ",")),
  paste("minimum_events:", min_events),
  paste("minimum_events_per_parameter:", min_epp),
  paste("ELSA_death_cutoff:", elsa_death_cutoff),
  paste("R:", R.version.string),
  paste("survival:", as.character(packageVersion("survival"))),
  "Independent unit: participant id; cluster-robust covariance for repeated intervals",
  "Inference: hazard ratios with 95% confidence intervals and exact two-sided p values",
  "Multiplicity: BH FDR within cohort-scenario for supplementary exposure estimates",
  "Heart-definition bounds: ELSA-only primary versus lower/upper state trajectories; transition-specific M4 Cox models",
  paste0("ELSA death boundary: risk intervals truncated at ", elsa_death_cutoff)
)
writeLines(session_lines, file.path(qc_root, "protocol_sensitivity_session_info.txt"), useBytes = TRUE)

cat("Protocol sensitivity workflow completed.\n")
cat("Run mode:", run_mode, "\n")
cat("Specifications:", nrow(specifications), "\n")
if (nrow(model_registry)) print(model_registry[, .N, by = status][order(status)])
cat("Coverage report:", file.path(report_root, "protocol_sensitivity_coverage.csv"), "\n")
