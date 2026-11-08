suppressPackageStartupMessages({
  library(data.table)
})

output_root <- normalizePath(
  Sys.getenv(
    "GA_MULTISTATE_OUTPUT_ROOT",
    unset = "private/stage/multistate_output"
  ),
  winslash = "/",
  mustWork = TRUE
)

result_root <- file.path(output_root, "results", "protocol_sensitivities")
qc_root <- file.path(output_root, "qc", "protocol_sensitivities")
data_root <- file.path(output_root, "data_multistate")
model_root <- file.path(output_root, "models", "protocol_sensitivities")
coverage_path <- file.path(output_root, "reports", "protocol_sensitivity_coverage.csv")
death_cutoff <- 2012.999

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  fread(path, na.strings = c("", "NA"), showProgress = FALSE)
}

checks <- list()
add_check <- function(name, pass, detail) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name,
    status = if (isTRUE(pass)) "PASS" else "FAIL",
    detail = as.character(detail)
  )
}

table_differences <- function(observed, expected, keys, tolerance = 1e-8) {
  observed <- copy(as.data.table(observed))
  expected <- copy(as.data.table(expected))
  if (!setequal(names(observed), names(expected))) {
    return(paste0(
      "column mismatch; observed-only=", paste(setdiff(names(observed), names(expected)), collapse = "|"),
      "; expected-only=", paste(setdiff(names(expected), names(observed)), collapse = "|")
    ))
  }
  setcolorder(expected, names(observed))
  setorderv(observed, keys, na.last = TRUE)
  setorderv(expected, keys, na.last = TRUE)
  if (nrow(observed) != nrow(expected)) {
    return(sprintf("row-count mismatch: observed=%d expected=%d", nrow(observed), nrow(expected)))
  }
  failures <- character()
  for (column in names(observed)) {
    a <- observed[[column]]
    b <- expected[[column]]
    if (is.numeric(a) || is.integer(a) || is.numeric(b) || is.integer(b)) {
      a <- suppressWarnings(as.numeric(a))
      b <- suppressWarnings(as.numeric(b))
      same <- (is.na(a) & is.na(b)) |
        (!is.na(a) & !is.na(b) & abs(a - b) <= tolerance * pmax(1, abs(a), abs(b)))
    } else {
      a <- as.character(a)
      b <- as.character(b)
      same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
    }
    if (any(!same)) failures <- c(failures, paste0(column, "[", sum(!same), "]"))
  }
  paste(failures, collapse = "; ")
}

registry <- read_csv(file.path(result_root, "protocol_sensitivity_model_registry.csv"))
specifications <- read_csv(file.path(result_root, "protocol_sensitivity_specifications.csv"))
estimates <- read_csv(file.path(result_root, "protocol_sensitivity_estimates.csv"))
joint_tests <- read_csv(file.path(result_root, "protocol_sensitivity_joint_tests.csv"))
coverage <- read_csv(coverage_path)

edges <- data.table(
  transition_id = sprintf("T%02d", 1:8),
  from_state = c(0L, 0L, 0L, 1L, 1L, 2L, 2L, 3L),
  to_state = c(1L, 2L, 4L, 3L, 4L, 3L, 4L, 4L)
)
edges[, transition := paste0("S", from_state, "->S", to_state)]

scenario_cohorts <- list(
  rcs = c("HRS", "ELSA"),
  alternate_death = "HRS",
  attained_age_period = c("HRS", "ELSA"),
  ipcw = c("HRS", "ELSA"),
  proxy_exclusion = "HRS",
  expanded_hypertension = "ELSA",
  semi_markov = c("HRS", "ELSA"),
  heart_lower = "ELSA",
  heart_upper = "ELSA"
)
expected_rows <- list()
for (scenario in names(scenario_cohorts)) {
  scenario_edges <- if (scenario == "alternate_death") edges[to_state == 4L] else edges
  for (cohort in scenario_cohorts[[scenario]]) {
    expected_rows[[length(expected_rows) + 1L]] <- scenario_edges[, .(
      cohort = cohort,
      scenario = scenario,
      transition_id,
      transition,
      from_state,
      to_state
    )]
  }
}
expected_specs <- rbindlist(expected_rows)
observed_specs <- specifications[, .(cohort, scenario, transition_id, transition, from_state, to_state)]
spec_diff <- table_differences(
  observed_specs,
  expected_specs,
  c("cohort", "scenario", "transition_id")
)
add_check(
  "Expected cohort-scenario-transition catalog",
  !nzchar(spec_diff) && nrow(specifications) == 100L && uniqueN(specifications$spec_id) == 100L,
  if (nzchar(spec_diff)) spec_diff else "100 unique specifications exactly match the prespecified catalog"
)

add_check(
  "Registry one-to-one with specifications",
  nrow(registry) == 100L && uniqueN(registry$spec_id) == 100L &&
    setequal(registry$spec_id, specifications$spec_id),
  sprintf("registry rows=%d; unique IDs=%d", nrow(registry), uniqueN(registry$spec_id))
)

terminal_statuses <- c("completed", "not_estimable", "failed")
nonterminal <- registry[!status %chin% terminal_statuses]
add_check(
  "All specifications have a terminal status",
  !nrow(nonterminal),
  if (nrow(nonterminal)) paste(nonterminal$spec_id, collapse = "|") else
    paste(registry[, .N, by = status][order(status)], collapse = "; ")
)

failures <- registry[status == "failed"]
add_check(
  "No model failures",
  !nrow(failures),
  if (nrow(failures)) paste(failures$spec_id, collapse = "|") else "0 failed models"
)

not_estimable <- registry[status == "not_estimable"]
expected_inestimable <- not_estimable[
  cohort == "ELSA" & to_state == 4L &
    grepl("events per parameter .* < 5", reason)
]
add_check(
  "Not-estimable models are the expected sparse ELSA death specifications",
  nrow(not_estimable) == 28L && nrow(expected_inestimable) == 28L,
  sprintf("not estimable=%d; expected sparse ELSA death specifications=%d", nrow(not_estimable), nrow(expected_inestimable))
)

elsa_death <- registry[cohort == "ELSA" & to_state == 4L]
death_fields_ok <- nrow(elsa_death) == 28L &&
  all(abs(elsa_death$mortality_cutoff_applied - death_cutoff) < 1e-12) &&
  all(grepl("2012\\.999", elsa_death$mortality_scope)) &&
  all(is.finite(elsa_death$post_cutoff_intervals_excluded) & elsa_death$post_cutoff_intervals_excluded > 0) &&
  all(is.finite(elsa_death$intervals_truncated_at_cutoff) & elsa_death$intervals_truncated_at_cutoff > 0)
add_check(
  "ELSA death cutoff and audit fields",
  death_fields_ok,
  sprintf(
    "rows=%d; cutoff range=%s; excluded range=%d-%d; truncated range=%d-%d",
    nrow(elsa_death),
    paste(range(elsa_death$mortality_cutoff_applied), collapse = "-"),
    min(elsa_death$post_cutoff_intervals_excluded), max(elsa_death$post_cutoff_intervals_excluded),
    min(elsa_death$intervals_truncated_at_cutoff), max(elsa_death$intervals_truncated_at_cutoff)
  )
)

elsa_nonfatal <- registry[cohort == "ELSA" & to_state != 4L]
add_check(
  "ELSA nonfatal specifications remain through W10",
  nrow(elsa_nonfatal) == 28L && all(is.na(elsa_nonfatal$mortality_cutoff_applied)) &&
    all(elsa_nonfatal$post_cutoff_intervals_excluded == 0L) &&
    all(elsa_nonfatal$intervals_truncated_at_cutoff == 0L) &&
    all(elsa_nonfatal$mortality_scope == "ELSA W4-W10 nonfatal follow-up"),
  sprintf("rows=%d; cutoff fields are blank and exclusion/truncation counts are zero", nrow(elsa_nonfatal))
)

completed <- registry[status == "completed"]
expected_estimate_rows <- completed[, .(
  expected_terms = if (scenario == "rcs") 3L else 1L
), by = .(spec_id, cohort, scenario)]
observed_estimate_rows <- estimates[, .(observed_terms = .N), by = .(spec_id, cohort, scenario)]
estimate_term_audit <- merge(expected_estimate_rows, observed_estimate_rows,
  by = c("spec_id", "cohort", "scenario"), all = TRUE
)
add_check(
  "Each completed model has all expected exposure estimates",
  nrow(estimate_term_audit) == nrow(completed) &&
    all(estimate_term_audit$expected_terms == estimate_term_audit$observed_terms) &&
    all(estimates$spec_id %chin% completed$spec_id) &&
    !anyDuplicated(estimates[, paste(spec_id, term)]),
  sprintf("completed models=%d; estimate rows=%d", nrow(completed), nrow(estimates))
)

rcs_terms_ok <- estimates[
  scenario == "rcs",
  all(sort(unique(term)) == sort(c(
    "insomnia_score_1_vs_0", "insomnia_score_2_vs_0", "insomnia_score_3_vs_0"
  ))),
  by = spec_id
][, all(V1)]
linear_terms_ok <- estimates[scenario != "rcs", all(term == "insomnia_score_primary")]
numerics_ok <- estimates[, all(
  is.finite(hazard_ratio) & hazard_ratio > 0 &
    is.finite(conf_low) & conf_low > 0 &
    is.finite(conf_high) & conf_high > 0 &
    conf_low <= hazard_ratio & hazard_ratio <= conf_high &
    is.finite(p_value) & p_value >= 0 & p_value <= 1
)]
add_check(
  "Estimate term labels and numerical bounds",
  rcs_terms_ok && linear_terms_ok && numerics_ok && all(estimates$status == "completed"),
  "RCS has three prespecified score-versus-zero contrasts; other scenarios have one score term; HR/CI/p values are finite and bounded"
)

fdr_audit <- estimates[, {
  recomputed <- p.adjust(p_value, method = "BH")
  list(
    rows = .N,
    registered_family_size_min = min(fdr_family_size),
    registered_family_size_max = max(fdr_family_size),
    max_absolute_difference = max(abs(recomputed - p_fdr_bh_sensitivity))
  )
}, by = .(cohort, scenario)]
fdr_expected <- completed[, .(
  expected_family_size = .N * if (first(scenario) == "rcs") 3L else 1L
), by = .(cohort, scenario)]
fdr_audit <- merge(fdr_audit, fdr_expected, by = c("cohort", "scenario"), all = TRUE)
fdr_ok <- fdr_audit[, all(
  rows == expected_family_size &
    registered_family_size_min == rows & registered_family_size_max == rows &
    max_absolute_difference < 1e-12
)]
add_check(
  "BH families are structurally complete and reproducible",
  fdr_ok,
  paste(
    fdr_audit[order(cohort, scenario), paste0(cohort, "/", scenario, "=", rows)],
    collapse = "; "
  )
)

completed_rcs <- completed[scenario == "rcs"]
joint_counts <- joint_tests[, .N, by = spec_id]
add_check(
  "RCS joint tests are complete",
  nrow(joint_counts) == nrow(completed_rcs) && all(joint_counts$N == 2L) &&
    setequal(joint_counts$spec_id, completed_rcs$spec_id) &&
    all(joint_tests$status == "completed") &&
    all(joint_tests$test %chin% c("RCS nonlinearity", "RCS overall association")),
  sprintf("completed RCS models=%d; joint-test rows=%d", nrow(completed_rcs), nrow(joint_tests))
)

model_files <- list.files(model_root, pattern = "\\.rds$", full.names = FALSE)
expected_model_files <- paste0(completed$spec_id, ".rds")
add_check(
  "Saved model objects match completed registry",
  setequal(model_files, expected_model_files) && length(model_files) == nrow(completed),
  sprintf("saved model objects=%d; completed registry rows=%d", length(model_files), nrow(completed))
)

add_check(
  "Coverage matrix is complete and failure-free",
  nrow(coverage) == 32L && uniqueN(coverage[, paste(sensitivity_id, cohort)]) == 32L &&
    setequal(coverage$sensitivity_id, 1:16) && setequal(coverage$cohort, c("HRS", "ELSA")) &&
    !any(coverage$status == "failed") && !any(coverage$failed_models > 0, na.rm = TRUE),
  paste(coverage[, .N, by = status][order(status), paste0(status, "=", N)], collapse = "; ")
)

state_label <- function(x) fifelse(is.na(x), "Missing", paste0("S", as.integer(x)))

rebuild_state_audit <- function(panel, scenario, state_column, assumption) {
  x <- copy(panel)
  x[, `:=`(primary_state = as.integer(state), sensitivity_state = as.integer(get(state_column)))]
  in_followup <- x$analysis_sample %in% 1L & is.finite(x$exposure_entry_time) &
    is.finite(x$event_time) & x$event_time >= x$exposure_entry_time
  scoped <- rbindlist(list(
    copy(x)[, analysis_scope := "all_panel_rows"],
    copy(x[in_followup])[, analysis_scope := "analysis_followup_all_states"],
    copy(x[in_followup & dead %in% 0L])[, analysis_scope := "analysis_followup_living_states"]
  ))
  scoped[, `:=`(
    scenario = scenario,
    state_column = state_column,
    bound_assumption = assumption,
    both_observed = as.integer(!is.na(primary_state) & !is.na(sensitivity_state)),
    observed_state_changed = as.integer(!is.na(primary_state) & !is.na(sensitivity_state) & primary_state != sensitivity_state),
    identical_including_joint_missing = as.integer(
      (is.na(primary_state) & is.na(sensitivity_state)) |
        (!is.na(primary_state) & !is.na(sensitivity_state) & primary_state == sensitivity_state)
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
  crosswalk <- scoped[, .(rows = .N, participants = uniqueN(id)), by = .(
    scenario, state_column, bound_assumption, analysis_scope, primary_state, sensitivity_state
  )]
  crosswalk[, `:=`(
    primary_state_label = state_label(primary_state),
    sensitivity_state_label = state_label(sensitivity_state),
    observed_state_changed = as.integer(!is.na(primary_state) & !is.na(sensitivity_state) & primary_state != sensitivity_state)
  )]
  setcolorder(crosswalk, c(
    "scenario", "state_column", "bound_assumption", "analysis_scope",
    "primary_state", "primary_state_label", "sensitivity_state",
    "sensitivity_state_label", "observed_state_changed", "rows", "participants"
  ))
  list(summary = summary, crosswalk = crosswalk)
}

build_intervals <- function(panel, state_column) {
  x <- copy(panel)
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
  out <- x[!is.na(next_state_sensitivity) & !is.na(next_time_sensitivity), .(
    id, cohort, wave, event_time, exposure_entry_time,
    from_state = state_sensitivity,
    to_state = next_state_sensitivity,
    next_time_sensitivity
  )]
  out[, `:=`(
    tstart = as.numeric(event_time - exposure_entry_time),
    tstop = as.numeric(next_time_sensitivity - exposure_entry_time),
    interval_years = as.numeric(next_time_sensitivity - event_time),
    state_changed = as.integer(from_state != to_state)
  )]
  setorder(out, id, tstart, tstop)
  out[]
}

audit_edges <- rbindlist(list(
  edges[, .(
    transition_id, from_state, to_state, transition,
    death_transition = as.integer(to_state == 4L),
    confirmatory_nonfatal = as.integer(paste0(from_state, "->", to_state) %chin% c("0->1", "0->2", "1->3", "2->3")),
    edge_role = "frozen_model_intensity_edge"
  )],
  data.table(
    transition_id = "PANEL_S0_S3", from_state = 0L, to_state = 3L,
    transition = "S0->S3", death_transition = 0L, confirmatory_nonfatal = 0L,
    edge_role = "panel_observed_interval_censored_multistep_jump"
  )
), use.names = TRUE)
audit_edges[, followup_scope := fifelse(
  death_transition == 1L,
  paste0("ELSA mortality risk truncated at ", death_cutoff),
  "ELSA nonfatal follow-up through W10"
)]

summarise_edges <- function(intervals, prefix) {
  out <- rbindlist(lapply(seq_len(nrow(audit_edges)), function(i) {
    edge <- audit_edges[i]
    risk_all <- copy(intervals[from_state == edge$from_state])
    risk_all[, calendar_interval_end := exposure_entry_time + tstop]
    excluded <- 0L
    if (edge$death_transition == 1L) {
      excluded <- risk_all[event_time >= death_cutoff, .N]
      risk <- risk_all[event_time < death_cutoff]
      risk[, audit_interval_years := pmax(0, pmin(calendar_interval_end, death_cutoff) - event_time)]
      risk <- risk[audit_interval_years > 0]
      events <- risk[state_changed == 1L & to_state == edge$to_state & calendar_interval_end <= death_cutoff]
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
      post_cutoff_risk_intervals_excluded = excluded
    ))
  }))
  metrics <- c(
    "risk_intervals", "risk_participants", "target_events", "event_participants",
    "risk_person_years", "post_cutoff_risk_intervals_excluded"
  )
  setnames(out, metrics, paste0(prefix, metrics))
  out
}

rebuild_transition_audit <- function(primary, sensitivity, scenario, state_column, assumption) {
  key <- names(audit_edges)
  out <- merge(summarise_edges(primary, "primary_"), summarise_edges(sensitivity, "sensitivity_"),
    by = key, all = TRUE, sort = FALSE
  )
  out[, `:=`(
    scenario = scenario,
    state_column = state_column,
    bound_assumption = assumption,
    risk_interval_delta = sensitivity_risk_intervals - primary_risk_intervals,
    risk_participant_delta = sensitivity_risk_participants - primary_risk_participants,
    event_delta = sensitivity_target_events - primary_target_events,
    event_participant_delta = sensitivity_event_participants - primary_event_participants,
    risk_person_year_delta = sensitivity_risk_person_years - primary_risk_person_years,
    post_cutoff_exclusion_delta = sensitivity_post_cutoff_risk_intervals_excluded - primary_post_cutoff_risk_intervals_excluded
  )]
  setcolorder(out, c(
    "scenario", "state_column", "bound_assumption", "transition_id", "transition",
    "from_state", "to_state", "edge_role", "death_transition", "confirmatory_nonfatal",
    setdiff(names(out), c(
      "scenario", "state_column", "bound_assumption", "transition_id", "transition",
      "from_state", "to_state", "edge_role", "death_transition", "confirmatory_nonfatal"
    ))
  ))
  out[]
}

rebuild_interval_crosswalk <- function(primary, sensitivity, scenario, state_column, assumption) {
  signature <- function(x, prefix) {
    z <- x[, .(id, tstart = as.numeric(tstart), tstop = as.numeric(tstop),
      from_state = as.integer(from_state), to_state = as.integer(to_state), present = 1L
    )]
    setnames(z, c("from_state", "to_state", "present"), paste0(prefix, c("from_state", "to_state", "present")))
    z
  }
  out <- merge(signature(primary, "primary_"), signature(sensitivity, "sensitivity_"),
    by = c("id", "tstart", "tstop"), all = TRUE, sort = FALSE
  )
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
  out[, `:=`(scenario = scenario, state_column = state_column, bound_assumption = assumption)]
  setcolorder(out, c(
    "scenario", "state_column", "bound_assumption", "primary_transition",
    "sensitivity_transition", "transition_class_changed", "intervals", "participants",
    "interval_person_years"
  ))
  out[]
}

panel <- as.data.table(readRDS(file.path(data_root, "elsa_exposure_panel_augmented.rds")))
primary_intervals <- as.data.table(readRDS(file.path(data_root, "elsa_multistate_intervals.rds")))
definitions <- data.table(
  scenario = c("heart_lower", "heart_upper"),
  state_column = c("state_heart_lower", "state_heart_upper"),
  bound_assumption = c(
    "Unresolved historical other-heart component coded absent",
    "Unresolved historical other-heart component coded present"
  )
)

state_summaries <- list()
state_crosswalks <- list()
transition_audits <- list()
interval_crosswalks <- list()
for (i in seq_len(nrow(definitions))) {
  d <- definitions[i]
  state_audit <- rebuild_state_audit(panel, d$scenario, d$state_column, d$bound_assumption)
  variant <- build_intervals(panel, d$state_column)
  state_summaries[[i]] <- state_audit$summary
  state_crosswalks[[i]] <- state_audit$crosswalk
  transition_audits[[i]] <- rebuild_transition_audit(
    primary_intervals, variant, d$scenario, d$state_column, d$bound_assumption
  )
  interval_crosswalks[[i]] <- rebuild_interval_crosswalk(
    primary_intervals, variant, d$scenario, d$state_column, d$bound_assumption
  )
}

heart_tables <- list(
  state_summary = list(
    observed = read_csv(file.path(qc_root, "heart_definition_state_summary.csv")),
    expected = rbindlist(state_summaries),
    keys = c("scenario", "analysis_scope")
  ),
  state_crosswalk = list(
    observed = read_csv(file.path(qc_root, "heart_definition_state_crosswalk.csv")),
    expected = rbindlist(state_crosswalks),
    keys = c("scenario", "analysis_scope", "primary_state", "sensitivity_state")
  ),
  transition_audit = list(
    observed = read_csv(file.path(qc_root, "heart_definition_transition_audit.csv")),
    expected = rbindlist(transition_audits),
    keys = c("scenario", "transition_id")
  ),
  interval_crosswalk = list(
    observed = read_csv(file.path(qc_root, "heart_definition_interval_crosswalk.csv")),
    expected = rbindlist(interval_crosswalks),
    keys = c("scenario", "primary_transition", "sensitivity_transition")
  )
)

heart_differences <- vapply(names(heart_tables), function(name) {
  item <- heart_tables[[name]]
  table_differences(item$observed, item$expected, item$keys)
}, character(1L))
add_check(
  "Heart-bound QC tables independently reproduce from the current final panel",
  all(!nzchar(heart_differences)),
  if (all(!nzchar(heart_differences))) {
    paste(names(heart_tables), vapply(heart_tables, function(x) nrow(x$observed), integer(1L)), sep = "=", collapse = "; ")
  } else {
    paste(names(heart_differences)[nzchar(heart_differences)], heart_differences[nzchar(heart_differences)], sep = ": ", collapse = "; ")
  }
)

living <- panel[dead != 1L]
ordering_ok <- living[, all(
  is.na(state_heart_lower) | is.na(state) | state_heart_lower <= state
)] && living[, all(
  is.na(state) | is.na(state_heart_upper) | state <= state_heart_upper
)]
state_output <- heart_tables$state_summary$observed
changed_lower <- state_output[scenario == "heart_lower" & analysis_scope == "analysis_followup_all_states", changed_observed_state_rows]
changed_upper <- state_output[scenario == "heart_upper" & analysis_scope == "analysis_followup_all_states", changed_observed_state_rows]
changed_upper_people <- state_output[scenario == "heart_upper" & analysis_scope == "analysis_followup_all_states", participants_with_changed_observed_state]
transition_output <- heart_tables$transition_audit$observed
lower_abs_event_delta <- transition_output[
  scenario == "heart_lower" & edge_role == "frozen_model_intensity_edge",
  sum(abs(event_delta))
]
upper_abs_event_delta <- transition_output[
  scenario == "heart_upper" & edge_role == "frozen_model_intensity_edge",
  sum(abs(event_delta))
]
add_check(
  "Heart-bound ordering and materiality",
  ordering_ok && changed_lower == 0L && changed_upper == 134L &&
    changed_upper_people == 41L && lower_abs_event_delta == 0L && upper_abs_event_delta == 40L,
  sprintf(
    "lower changed rows=%d and absolute eight-edge event delta=%d; upper changed rows=%d in %d participants and absolute event delta=%d",
    changed_lower, lower_abs_event_delta, changed_upper, changed_upper_people, upper_abs_event_delta
  )
)

check_table <- rbindlist(checks)
print(check_table)
cat("\nSummary:\n")
print(check_table[, .N, by = status][order(status)])
if (check_table[status == "FAIL", .N]) quit(status = 1L)
