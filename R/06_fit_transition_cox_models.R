#!/usr/bin/env Rscript

# Transition-specific Cox models for the HRS/ELSA multistate study.
#
# IMPORTANT: these models treat the destination state observed at a visit as if
# the transition occurred at that visit.  They are therefore a semiparametric
# sensitivity analysis for panel-observed disease states, not the primary
# continuous-time panel Markov analysis.  For each of the eight protocol edges,
# another exit from the same origin state is treated as a competing censoring
# event.  S0->S3 is never a Cox target event.

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(mice)
})

options(stringsAsFactors = FALSE, width = 160, warn = 1)

`%||%` <- function(x, y) if (is.null(x) || !length(x) || all(is.na(x))) y else x

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_argument)) {
  sub("^--file=", "", script_argument[1L])
} else {
  file.path(getwd(), "R", "06_fit_transition_cox_models.R")
}
project_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = FALSE)

as_flag <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) return(default)
  tolower(trimws(x)) %in% c("1", "true", "t", "yes", "y", "on")
}

as_integer_setting <- function(name, default, lower = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || value < lower) default else value
}

as_numeric_setting <- function(name, default, lower = 0) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value) || value < lower) default else value
}

split_setting <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  unique(trimws(strsplit(value, ",", fixed = TRUE)[[1L]]))
}

safe_label <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

sha256_object <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)

canonical_scalar <- function(name, value) {
  if (length(value) != 1L) stop("Canonical hash fields must be scalar: ", name)
  type <- typeof(value)
  rendered <- if (is.na(value)) "<NA>" else enc2utf8(as.character(value))
  paste0(nchar(enc2utf8(name), type = "bytes"), ":", enc2utf8(name), "|",
         nchar(type, type = "bytes"), ":", type, "|",
         nchar(rendered, type = "bytes"), ":", rendered)
}

canonical_named_payload <- function(values) {
  if (is.null(names(values)) || any(!nzchar(names(values)))) stop("Canonical payload must be fully named")
  fields <- sort(names(values))
  paste(vapply(fields, function(name) canonical_scalar(name, values[[name]]), character(1L)), collapse = "||")
}

sha256_named_fields <- function(values) {
  digest::digest(canonical_named_payload(values), algo = "sha256", serialize = FALSE)
}

sha256_data_table <- function(x, order_by = names(x)) {
  x <- as.data.table(copy(x))
  order_by <- intersect(order_by, names(x))
  if (length(order_by)) setorderv(x, order_by)
  fields <- sort(names(x))
  rows <- vapply(seq_len(nrow(x)), function(i) {
    canonical_named_payload(as.list(x[i, ..fields]))
  }, character(1L))
  digest::digest(paste(rows, collapse = "\n<ROW>\n"), algo = "sha256", serialize = FALSE)
}

same_with_na <- function(lhs, rhs) {
  lhs <- as.character(lhs)
  rhs <- as.character(rhs)
  if (length(lhs) != length(rhs)) return(rep(FALSE, max(length(lhs), length(rhs))))
  both_na <- is.na(lhs) & is.na(rhs)
  neither_na <- !is.na(lhs) & !is.na(rhs)
  both_na | (neither_na & lhs == rhs)
}

stage_root <- Sys.getenv("GA_STAGE_ROOT", unset = "private/stage")
configured_multistate_root <- Sys.getenv("GA_MULTISTATE_OUTPUT_ROOT", unset = "")
multistate_root <- if (nzchar(configured_multistate_root)) {
  configured_multistate_root
} else if (file.exists(file.path(project_root, "data_multistate", "hrs_multistate_intervals.rds"))) {
  project_root
} else {
  file.path(stage_root, "multistate_output")
}
input_root <- Sys.getenv(
  "GA_MULTISTATE_DATA_ROOT",
  unset = file.path(multistate_root, "data_multistate")
)
output_root <- Sys.getenv("GA_COX_OUTPUT_ROOT", unset = multistate_root)
mice_cache_root <- Sys.getenv(
  "GA_COX_MICE_CACHE_ROOT",
  unset = file.path(multistate_root, "data_multistate", "cox_mice_cache")
)

analysis_tier <- tolower(Sys.getenv("GA_COX_TIER", unset = "core"))
if (!analysis_tier %in% c("spec", "core", "full")) {
  stop("GA_COX_TIER must be one of: spec, core, full")
}
analysis_modules <- tolower(split_setting(
  "GA_COX_MODULES",
  if (analysis_tier == "full") "all" else "main"
))
valid_modules <- c("all", "main", "design", "alternative", "dynamic", "hrs_specific", "interaction")
if (any(!analysis_modules %in% valid_modules)) {
  stop("Unknown GA_COX_MODULES value(s): ", paste(setdiff(analysis_modules, valid_modules), collapse = ", "))
}

selected_cohorts <- toupper(split_setting("GA_COX_COHORTS", "HRS,ELSA"))
selected_transitions <- toupper(split_setting("GA_COX_TRANSITIONS", "ALL"))

config <- list(
  tier = analysis_tier,
  modules = analysis_modules,
  cohorts = selected_cohorts,
  transitions = selected_transitions,
  run_mice = as_flag(Sys.getenv("GA_RUN_MICE", unset = "true"), TRUE),
  run_weighted = as_flag(Sys.getenv("GA_RUN_WEIGHTED", unset = "false"), FALSE),
  weighted_only = as_flag(Sys.getenv("GA_COX_WEIGHTED_ONLY", unset = "false"), FALSE),
  run_zph = as_flag(Sys.getenv("GA_RUN_ZPH", unset = "true"), TRUE),
  save_model_objects = as_flag(Sys.getenv("GA_SAVE_COX_OBJECTS", unset = "false"), FALSE),
  save_imputation_estimates = as_flag(Sys.getenv("GA_SAVE_IMPUTATION_ESTIMATES", unset = "false"), FALSE),
  mice_m = as_integer_setting("GA_MICE_M", 20L, 2L),
  mice_maxit = as_integer_setting("GA_MICE_MAXIT", 10L, 1L),
  seed = as_integer_setting("GA_COX_SEED", 20260908L, 1L),
  min_events = as_integer_setting("GA_MIN_EVENTS", 10L, 1L),
  min_events_per_parameter = as_numeric_setting("GA_MIN_EVENTS_PER_PARAMETER", 10, 0),
  min_interaction_cell_events = as_integer_setting("GA_MIN_INTERACTION_CELL_EVENTS", 10L, 1L),
  elsa_death_cutoff = as.numeric(Sys.getenv("GA_ELSA_DEATH_CUTOFF", unset = "2012.999")),
  resume_checkpoints = as_flag(Sys.getenv("GA_COX_RESUME_CHECKPOINTS", unset = "true"), TRUE),
  publish_canonical = as_flag(Sys.getenv("GA_COX_PUBLISH_CANONICAL", unset = "false"), FALSE),
  force_publish_canonical = as_flag(Sys.getenv("GA_COX_FORCE_PUBLISH_CANONICAL", unset = "false"), FALSE),
  spec_ids = split_setting("GA_COX_SPEC_IDS", "ALL")
)

if (!is.finite(config$elsa_death_cutoff)) stop("GA_ELSA_DEATH_CUTOFF must be a finite numeric year")

run_label_default <- paste(
  config$tier,
  paste(config$modules, collapse = "-"),
  paste(config$cohorts, collapse = "-"),
  sep = "_"
)
run_label_requested <- safe_label(Sys.getenv("GA_COX_RUN_LABEL", unset = run_label_default))

build_input_manifest <- function(cohorts) {
  rbindlist(lapply(sort(unique(cohorts)), function(cohort_name) {
    path <- file.path(input_root, paste0(tolower(cohort_name), "_multistate_intervals.rds"))
    info <- suppressWarnings(file.info(path))
    data.table(
      cohort = cohort_name,
      path = normalizePath(path, winslash = "/", mustWork = FALSE),
      exists = file.exists(path),
      bytes = if (file.exists(path)) as.numeric(info$size) else NA_real_,
      sha256 = sha256_file(path)
    )
  }), use.names = TRUE, fill = TRUE)
}

input_manifest <- build_input_manifest(config$cohorts)

script_sha256 <- sha256_file(script_file)
# A fit-only hotfix does not change the imputation frame, methods, predictor
# matrix, or protected-field rules.  This explicit, audited override permits a
# cache created by the verified pre-hotfix script to be reused; all substantive
# frame and plan hashes below must still match exactly.
mice_cache_script_sha256 <- tolower(trimws(Sys.getenv(
  "GA_COX_MICE_CACHE_SCRIPT_SHA256", unset = script_sha256
)))
if (!grepl("^[0-9a-f]{64}$", mice_cache_script_sha256)) {
  stop("GA_COX_MICE_CACHE_SCRIPT_SHA256 must be a 64-character SHA-256 digest")
}
fingerprint_config <- config
fingerprint_config$cohorts <- sort(unique(fingerprint_config$cohorts))
fingerprint_config$modules <- if ("all" %chin% fingerprint_config$modules) {
  sort(setdiff(valid_modules, "all"))
} else {
  sort(unique(fingerprint_config$modules))
}
fingerprint_config$transitions <- if ("ALL" %chin% fingerprint_config$transitions) {
  sprintf("T%02d", 1:8)
} else {
  sort(unique(fingerprint_config$transitions))
}
fingerprint_config$resume_checkpoints <- NULL
fingerprint_config$publish_canonical <- NULL
fingerprint_config$force_publish_canonical <- NULL
run_signature <- list(
  signature_version = "2026-09-08-v2",
  script_sha256 = script_sha256,
  input_manifest = input_manifest,
  analytic_configuration = fingerprint_config,
  R_version = as.character(getRversion()),
  survival_version = as.character(packageVersion("survival")),
  mice_version = as.character(packageVersion("mice")),
  data_table_version = as.character(packageVersion("data.table"))
)
base_run_fingerprint <- sha256_object(run_signature)
run_fingerprint <- base_run_fingerprint
run_label <- paste0(run_label_requested, "__", substr(run_fingerprint, 1L, 16L))
run_instance_default <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_pid", Sys.getpid())
run_instance_id <- safe_label(Sys.getenv("GA_COX_RUN_INSTANCE", unset = run_instance_default))
catalog_root <- file.path(output_root, "results", "cox_sensitivity")
result_root <- file.path(catalog_root, run_label, run_instance_id)
qc_root <- file.path(output_root, "qc", "cox_sensitivity", run_label, run_instance_id)
checkpoint_root <- file.path(result_root, "checkpoints")

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
  transition = paste0("S", from_state, "->S", to_state),
  from_label = unname(state_labels[as.character(from_state)]),
  to_label = unname(state_labels[as.character(to_state)]),
  death_transition = as.integer(to_state == 4L),
  confirmatory_nonfatal = as.integer(to_state != 4L)
)]

if (nrow(transition_catalog) != 8L ||
    transition_catalog[from_state == 0L & to_state == 3L, .N] != 0L) {
  stop("The protocol transition catalog must contain exactly eight edges and no direct S0->S3 intensity edge")
}

baseline_impute_vars <- c(
  "age_bl", "female_bl", "race_bl", "education_level_bl", "wealth_quintile_bl",
  "married_bl", "current_smoker_bl", "current_drinker_bl", "inactive_bl", "bmi_bl",
  "diabetes_bl", "lung_disease_bl", "frequent_pain_bl", "depression_score_bl"
)

protected_baseline_vars <- c(
  "id", "cohort", "entry_wave", "entry_year", "state_bl", "weight_bl",
  "insomnia_score_bl", "insomnia_burden_bl"
)

optional_interval_vars <- list(
  insomnia_score_cumavg = NA_real_,
  insomnia_score_prev_valid = NA_real_,
  insomnia_score_lag1_row = NA_real_,
  insomnia_dynamic_pattern = NA_character_,
  insomnia_dynamic_code = NA_integer_,
  insomnia_assessment_complete = NA_integer_,
  insomnia_assessment_number = NA_integer_,
  insomnia_source_wave = NA_integer_,
  insomnia_locf_reason = NA_character_,
  insomnia_score_extended4 = NA_real_,
  insomnia_early_waking_current = NA_integer_,
  sleep_medication_current = NA_integer_,
  proxy_interview = NA_integer_,
  interval_starts_at_observed_sleep = NA_integer_,
  interval_censored_multistep_jump = NA_integer_,
  allowed_intensity_edge = NA_integer_
)

atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  fwrite(x, tmp, na = "")
  if (file.exists(path) && !file.remove(path)) stop("Could not replace output: ", path)
  if (!file.rename(tmp, path)) stop("Could not atomically write output: ", path)
  invisible(path)
}

atomic_saveRDS <- function(x, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, compress = compress)
  if (file.exists(path) && !file.remove(path)) stop("Could not replace output: ", path)
  if (!file.rename(tmp, path)) stop("Could not atomically write output: ", path)
  invisible(path)
}

collapse_messages <- function(x) {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) paste(x, collapse = " | ") else ""
}

finite_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (any(is.finite(x))) min(x[is.finite(x)]) else NA_real_
}

finite_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
}

empty_table <- function(columns) {
  as.data.table(setNames(replicate(length(columns), character(), simplify = FALSE), columns))
}

bind_or_empty <- function(x, columns) {
  if (!length(x)) return(empty_table(columns))
  rbindlist(x, use.names = TRUE, fill = TRUE)
}

normalize_intervals <- function(x, cohort_name) {
  x <- as.data.table(copy(x))
  required <- c(
    "id", "cohort", "wave", "next_wave", "event_time", "next_time",
    "exposure_entry_time", "tstart", "tstop", "interval_years", "from_state", "to_state",
    "state_changed", "insomnia_score_primary", "insomnia_burden_code",
    "insomnia_onset_primary", "insomnia_maintenance_primary",
    "insomnia_nonrestorative_primary", protected_baseline_vars, baseline_impute_vars
  )
  missing_required <- setdiff(required, names(x))
  if (length(missing_required)) {
    stop(cohort_name, " interval input is missing required field(s): ", paste(missing_required, collapse = ", "))
  }
  for (nm in names(optional_interval_vars)) {
    if (!nm %in% names(x)) x[, (nm) := optional_interval_vars[[nm]]]
  }
  x[, cohort := toupper(as.character(cohort))]
  x[, `:=`(
    id = as.numeric(id),
    from_state = as.integer(from_state),
    to_state = as.integer(to_state),
    state_changed = as.integer(state_changed),
    tstart = as.numeric(tstart),
    tstop = as.numeric(tstop),
    interval_years = as.numeric(interval_years),
    event_time = as.numeric(event_time),
    next_time = as.numeric(next_time),
    exposure_entry_time = as.numeric(exposure_entry_time)
  )]
  setorder(x, id, tstart, tstop)
  x[]
}

model_terms <- list(
  M1 = c("age5", "female_factor"),
  M2 = c(
    "age5", "female_factor", "race_factor", "education_factor",
    "wealth_factor", "married_bl"
  ),
  M3 = c(
    "age5", "female_factor", "race_factor", "education_factor",
    "wealth_factor", "married_bl", "current_smoker_bl", "current_drinker_bl",
    "inactive_bl", "bmi5"
  ),
  M4 = c(
    "age5", "female_factor", "race_factor", "education_factor",
    "wealth_factor", "married_bl", "current_smoker_bl", "current_drinker_bl",
    "inactive_bl", "bmi5", "diabetes_bl", "lung_disease_bl", "frequent_pain_bl"
  )
)

make_spec_catalog <- function() {
  specs <- list()
  add_spec <- function(module, analysis_family, analysis_id, analysis_label, scenario,
                       exposure_var, exposure_type, model_level, basis, weighted = FALSE,
                       interaction_var = "", cohorts = c("HRS", "ELSA"), tier_required = "full") {
    for (cohort_name in cohorts) {
      for (i in seq_len(nrow(transition_catalog))) {
        tr <- transition_catalog[i]
        specs[[length(specs) + 1L]] <<- data.table(
          module = module,
          tier_required = tier_required,
          cohort = cohort_name,
          transition_id = tr$transition_id,
          transition = tr$transition,
          from_state = tr$from_state,
          to_state = tr$to_state,
          death_transition = tr$death_transition,
          confirmatory_nonfatal = tr$confirmatory_nonfatal,
          analysis_family = analysis_family,
          analysis_id = analysis_id,
          analysis_label = analysis_label,
          scenario = scenario,
          exposure_var = exposure_var,
          exposure_type = exposure_type,
          interaction_var = interaction_var,
          model_level = model_level,
          basis = basis,
          weighted = as.integer(weighted)
        )
      }
    }
  }

  main_exposures <- data.table(
    analysis_id = c("primary_score_tv", "primary_burden_tv"),
    analysis_label = c("Time-varying insomnia score, per 1 symptom", "Time-varying burden: 0 / 1 / >=2 symptoms"),
    exposure_var = c("insomnia_score_primary", "insomnia_burden_code"),
    exposure_type = c("continuous", "categorical")
  )
  for (j in seq_len(nrow(main_exposures))) {
    ex <- main_exposures[j]
    for (m in names(model_terms)) {
      bases <- c("mice_pooled", "complete_case")
      for (basis_name in bases) {
        add_spec(
          "main", "main", ex$analysis_id, ex$analysis_label, "standard",
          ex$exposure_var, ex$exposure_type, m, basis_name, FALSE,
          tier_required = "core"
        )
        add_spec(
          "main", "weight_sensitivity", paste0(ex$analysis_id, "_weighted"),
          paste0(ex$analysis_label, "; baseline-weight sensitivity"), "standard",
          ex$exposure_var, ex$exposure_type, m, basis_name, TRUE,
          tier_required = "core"
        )
      }
    }
  }

  design_sensitivities <- data.table(
    analysis_id = c(
      "landmark2_score", "landmark2_burden", "observed_sleep_score",
      "observed_sleep_burden", "interval_le4_score", "interval_le4_burden"
    ),
    analysis_label = c(
      "Two-year landmark; time-varying score", "Two-year landmark; time-varying burden",
      "Intervals starting at a complete sleep assessment; recorded proxies excluded; score",
      "Intervals starting at a complete sleep assessment; recorded proxies excluded; burden",
      "Intervals no longer than four years; score", "Intervals no longer than four years; burden"
    ),
    scenario = c("landmark2", "landmark2", "observed_sleep", "observed_sleep", "interval_le4", "interval_le4"),
    exposure_var = rep(c("insomnia_score_primary", "insomnia_burden_code"), 3L),
    exposure_type = rep(c("continuous", "categorical"), 3L)
  )
  for (j in seq_len(nrow(design_sensitivities))) {
    ex <- design_sensitivities[j]
    for (basis_name in c("mice_pooled", "complete_case")) {
      add_spec(
        "design", "design_sensitivity", ex$analysis_id, ex$analysis_label,
        ex$scenario, ex$exposure_var, ex$exposure_type, "M4", basis_name
      )
    }
  }

  alternative_exposures <- data.table(
    analysis_id = c(
      "baseline_score", "baseline_burden", "cumulative_average", "lag1_row",
      "previous_valid_assessment", "onset_symptom", "maintenance_symptom",
      "nonrestorative_symptom", "depression_adjusted_score", "depression_adjusted_burden"
    ),
    analysis_label = c(
      "Baseline-fixed insomnia score", "Baseline-fixed insomnia burden",
      "Cumulative mean of complete insomnia scores", "One ordered-record lag of insomnia score",
      "Immediately preceding complete insomnia assessment score", "Difficulty initiating sleep",
      "Night waking", "Nonrestorative sleep",
      "M4 plus harmonized CES-D summary; sleep-item overlap sensitivity",
      "M4 categorical burden plus harmonized CES-D summary; sleep-item overlap sensitivity"
    ),
    scenario = c(
      "standard", "standard", "standard", "lag_available", "previous_valid_available",
      "standard", "standard", "standard", "standard", "standard"
    ),
    exposure_var = c(
      "insomnia_score_bl", "insomnia_burden_bl", "insomnia_score_cumavg",
      "insomnia_score_lag1_row", "insomnia_score_prev_valid", "insomnia_onset_primary",
      "insomnia_maintenance_primary", "insomnia_nonrestorative_primary",
      "insomnia_score_primary", "insomnia_burden_code"
    ),
    exposure_type = c(
      "continuous", "categorical", "continuous", "continuous", "continuous",
      "binary", "binary", "binary", "continuous", "categorical"
    )
  )
  for (j in seq_len(nrow(alternative_exposures))) {
    ex <- alternative_exposures[j]
    family_name <- if (grepl("^depression_adjusted_", ex$analysis_id)) {
      "depression_sensitivity"
    } else {
      "alternative_exposure"
    }
    for (basis_name in c("mice_pooled", "complete_case")) {
      add_spec(
        "alternative", family_name, ex$analysis_id, ex$analysis_label,
        ex$scenario, ex$exposure_var, ex$exposure_type, "M4", basis_name
      )
    }
  }

  for (basis_name in c("mice_pooled", "complete_case")) {
    add_spec(
      "dynamic", "dynamic_pattern", "dynamic_pattern_delayed_entry",
      "Adjacent complete-assessment dynamic pattern; delayed entry at second complete assessment",
      "dynamic_delayed", "insomnia_dynamic_pattern", "dynamic", "M4", basis_name
    )
    add_spec(
      "hrs_specific", "hrs_specific", "hrs_extended_score4",
      "HRS extended insomnia score (0-4)", "standard", "insomnia_score_extended4",
      "continuous", "M4", basis_name, cohorts = "HRS"
    )
    add_spec(
      "hrs_specific", "hrs_specific", "hrs_sleep_medication",
      "HRS sleep-medication use", "standard", "sleep_medication_current",
      "binary", "M4", basis_name, cohorts = "HRS"
    )
  }

  interactions <- data.table(
    analysis_id = c("interaction_sex", "interaction_age", "interaction_obesity", "interaction_diabetes", "interaction_wealth"),
    analysis_label = c(
      "Insomnia score by sex", "Insomnia score by age", "Insomnia score by obesity",
      "Insomnia score by diabetes", "Insomnia score by wealth quintile"
    ),
    interaction_var = c("female_factor", "age5", "obesity_factor", "diabetes_factor", "wealth_factor")
  )
  for (j in seq_len(nrow(interactions))) {
    ex <- interactions[j]
    for (basis_name in c("mice_pooled", "complete_case")) {
      add_spec(
        "interaction", "effect_modification", ex$analysis_id, ex$analysis_label,
        "standard", "insomnia_score_primary", "interaction", "M4", basis_name,
        interaction_var = ex$interaction_var
      )
    }
  }

  out <- rbindlist(specs, use.names = TRUE, fill = TRUE)
  out[, covariate_set := vapply(seq_len(.N), function(i) {
    values <- covariate_rhs(out[i])
    if (length(values)) paste(values, collapse = " + ") else "none"
  }, character(1))]
  out[, formula_rhs := vapply(seq_len(.N), function(i) {
    paste(c(exposure_rhs(out[i]), covariate_rhs(out[i])), collapse = " + ")
  }, character(1))]
  out[, `:=`(
    estimation_role = "Visit-detection-time semiparametric sensitivity; not primary panel Markov model",
    event_definition = paste0("Target event only when from_state=", from_state, " and observed to_state=", to_state),
    competing_exit_rule = "Any other observed exit from the same origin state is censored at that visit",
    mortality_scope = fifelse(
      cohort == "ELSA" & to_state == 4L,
      paste0("limited mortality coverage; risk truncated at ", config$elsa_death_cutoff),
      fifelse(cohort == "ELSA", "nonfatal model may include W10", "HRS follow-up")
    )
  )]
  identity_fields <- c(
    "cohort", "transition_id", "from_state", "to_state", "analysis_family", "analysis_id",
    "scenario", "exposure_var", "exposure_type", "interaction_var", "model_level", "basis",
    "weighted", "formula_rhs"
  )
  identity <- lapply(seq_len(nrow(out)), function(i) {
    payload <- as.list(out[i, ..identity_fields])
    payload$specification_schema <- "cox-transition-v2-cluster-id-efron"
    payload$elsa_death_cutoff <- if (out$cohort[i] == "ELSA" && out$to_state[i] == 4L) {
      config$elsa_death_cutoff
    } else {
      NA_real_
    }
    canonical <- canonical_named_payload(payload)
    list(canonical = canonical, hash = digest::digest(canonical, algo = "sha256", serialize = FALSE))
  })
  out[, `:=`(
    spec_hash_payload = vapply(identity, `[[`, character(1L), "canonical"),
    spec_hash_full = vapply(identity, `[[`, character(1L), "hash")
  )]
  out[, spec_id := paste0("CXS_", substr(spec_hash_full, 1L, 20L))]
  if (anyDuplicated(out$spec_id)) {
    duplicate_ids <- unique(out$spec_id[duplicated(out$spec_id)])
    stop("Stable specification hash collision/duplicate content: ", paste(duplicate_ids, collapse = ", "))
  }
  setcolorder(out, c("spec_id", setdiff(names(out), "spec_id")))
  out[]
}

validate_spec_catalog <- function(catalog) {
  expected_counts <- data.table(
    module = c("main", "design", "alternative", "dynamic", "hrs_specific", "interaction"),
    expected_rows = c(512L, 192L, 320L, 32L, 32L, 160L)
  )
  observed_counts <- catalog[, .(observed_rows = .N), by = module]
  module_checks <- merge(expected_counts, observed_counts, by = "module", all = TRUE)
  module_checks[is.na(expected_rows), expected_rows := 0L]
  module_checks[is.na(observed_rows), observed_rows := 0L]
  module_checks[, `:=`(
    check = paste0("module_row_count_", module),
    status = fifelse(observed_rows == expected_rows, "pass", "fail"),
    detail = paste0("observed=", observed_rows, "; expected=", expected_rows)
  )]

  main_unweighted <- catalog[module == "main" & analysis_family == "main" & weighted == 0L]
  main_expected <- CJ(
    cohort = c("ELSA", "HRS"),
    transition_id = transition_catalog$transition_id,
    analysis_id = c("primary_burden_tv", "primary_score_tv"),
    model_level = names(model_terms),
    basis = c("complete_case", "mice_pooled"),
    unique = TRUE
  )
  main_observed <- unique(main_unweighted[, .(cohort, transition_id, analysis_id, model_level, basis)])
  setkeyv(main_expected, names(main_expected))
  setkeyv(main_observed, names(main_observed))
  missing_main <- fsetdiff(main_expected, main_observed)

  global_checks <- data.table(
    check = c(
      "catalog_total_rows", "unique_spec_ids", "eight_target_edges", "no_direct_S0_S3_target",
      "main_M1_M4_MI_and_CC_complete", "all_cohort_edge_pairs_present"
    ),
    status = c(
      if (nrow(catalog) == sum(expected_counts$expected_rows)) "pass" else "fail",
      if (uniqueN(catalog$spec_id) == nrow(catalog)) "pass" else "fail",
      if (uniqueN(catalog$transition_id) == 8L) "pass" else "fail",
      if (catalog[from_state == 0L & to_state == 3L, .N] == 0L) "pass" else "fail",
      if (!nrow(missing_main)) "pass" else "fail",
      if (uniqueN(catalog, by = c("cohort", "transition_id")) == 16L) "pass" else "fail"
    ),
    detail = c(
      paste0("observed=", nrow(catalog), "; expected=", sum(expected_counts$expected_rows)),
      paste0("unique=", uniqueN(catalog$spec_id), "; rows=", nrow(catalog)),
      paste0("observed=", uniqueN(catalog$transition_id), "; expected=8"),
      paste0("observed=", catalog[from_state == 0L & to_state == 3L, .N], "; expected=0"),
      paste0("missing combinations=", nrow(missing_main)),
      paste0("observed=", uniqueN(catalog, by = c("cohort", "transition_id")), "; expected=16")
    )
  )
  checks <- rbindlist(list(
    global_checks,
    module_checks[, .(check, status, detail)]
  ), use.names = TRUE)
  if (checks[status == "fail", .N]) {
    stop(
      "Cox specification catalog validation failed: ",
      paste(checks[status == "fail", paste0(check, " (", detail, ")")], collapse = "; ")
    )
  }
  checks[]
}

build_person_baseline <- function(intervals) {
  baseline_columns <- unique(c(protected_baseline_vars, baseline_impute_vars))
  baseline <- intervals[order(tstart), .SD[1L], by = id, .SDcols = setdiff(baseline_columns, "id")]
  auxiliary <- intervals[, .(
    followup_years = max(tstop, na.rm = TRUE),
    any_death = as.integer(any(state_changed == 1L & to_state == 4L)),
    any_nonfatal_transition = as.integer(any(state_changed == 1L & to_state != 4L)),
    n_observed_transitions = sum(state_changed == 1L),
    n_intervals = .N
  ), by = id]
  transition_auxiliary <- intervals[, {
    values <- lapply(seq_len(nrow(transition_catalog)), function(i) {
      tr <- transition_catalog[i]
      as.integer(any(state_changed == 1L & from_state == tr$from_state & to_state == tr$to_state))
    })
    names(values) <- paste0("event_", transition_catalog$transition_id)
    values
  }, by = id]
  out <- merge(baseline, auxiliary, by = "id", all.x = TRUE, sort = FALSE)
  merge(out, transition_auxiliary, by = "id", all.x = TRUE, sort = FALSE)
}

prepare_mice_frame <- function(baseline) {
  d <- as.data.frame(copy(baseline))
  d$cohort <- factor(d$cohort)
  d$state_bl <- factor(d$state_bl, levels = 0:3)
  d$insomnia_burden_bl <- factor(d$insomnia_burden_bl, levels = 0:2)
  d$female_bl <- factor(d$female_bl, levels = 0:1)
  d$race_bl <- factor(d$race_bl)
  d$education_level_bl <- factor(d$education_level_bl, levels = 1:3)
  d$wealth_quintile_bl <- factor(d$wealth_quintile_bl, levels = 1:5)
  binary_vars <- c(
    "married_bl", "current_smoker_bl", "current_drinker_bl", "inactive_bl",
    "diabetes_bl", "lung_disease_bl", "frequent_pain_bl", "any_death",
    "any_nonfatal_transition"
  )
  for (v in intersect(binary_vars, names(d))) d[[v]] <- factor(d[[v]], levels = 0:1)
  d
}

factor_code_numeric <- function(x) {
  if (is.factor(x) || is.character(x)) return(suppressWarnings(as.numeric(as.character(x))))
  suppressWarnings(as.numeric(x))
}

canonicalize_completed_baseline <- function(d) {
  d <- as.data.table(copy(d))
  numeric_coded <- intersect(c(
    "id", "entry_wave", "entry_year", "state_bl", "weight_bl", "insomnia_score_bl",
    "insomnia_burden_bl", "age_bl", "female_bl", "education_level_bl",
    "wealth_quintile_bl", "married_bl", "current_smoker_bl", "current_drinker_bl",
    "inactive_bl", "bmi_bl", "diabetes_bl", "lung_disease_bl", "frequent_pain_bl",
    "depression_score_bl", "followup_years", "any_death", "any_nonfatal_transition",
    "n_observed_transitions", "n_intervals"
  ), names(d))
  for (v in numeric_coded) set(d, j = v, value = factor_code_numeric(d[[v]]))
  if ("cohort" %in% names(d)) d[, cohort := as.character(cohort)]
  if ("race_bl" %in% names(d)) d[, race_bl := as.character(race_bl)]
  d[]
}

count_protected_changes <- function(completed, baseline, protected_fields) {
  protected_fields <- intersect(protected_fields, intersect(names(baseline), names(completed)))
  if (!length(protected_fields)) return(0L)
  reference <- copy(baseline)
  candidate <- copy(completed)
  setorder(reference, id)
  setorder(candidate, id)
  if (!identical(as.character(reference$id), as.character(candidate$id))) return(.Machine$integer.max)
  sum(vapply(protected_fields, function(v) {
    sum(!same_with_na(candidate[[v]], reference[[v]]))
  }, integer(1L)))
}

run_baseline_mice <- function(intervals, cohort_name) {
  baseline <- build_person_baseline(intervals)
  setorder(baseline, id)
  frame <- prepare_mice_frame(baseline)
  frame <- frame[order(frame$id), , drop = FALSE]
  method <- rep("", ncol(frame))
  names(method) <- names(frame)
  continuous <- c("age_bl", "bmi_bl", "depression_score_bl")
  binary <- c(
    "female_bl", "married_bl", "current_smoker_bl", "current_drinker_bl",
    "inactive_bl", "diabetes_bl", "lung_disease_bl", "frequent_pain_bl"
  )
  multicategory <- c("race_bl", "education_level_bl", "wealth_quintile_bl")
  for (v in intersect(continuous, names(frame))) if (anyNA(frame[[v]])) method[v] <- "pmm"
  for (v in intersect(binary, names(frame))) if (anyNA(frame[[v]])) method[v] <- "logreg"
  for (v in intersect(multicategory, names(frame))) if (anyNA(frame[[v]])) method[v] <- "polyreg"

  predictor <- mice::make.predictorMatrix(frame)
  predictor[,] <- 1L
  diag(predictor) <- 0L
  protected <- intersect(c(
    protected_baseline_vars, "followup_years", "any_death", "any_nonfatal_transition",
    "n_observed_transitions", "n_intervals", grep("^event_T", names(frame), value = TRUE)
  ), names(frame))
  predictor[protected, ] <- 0L
  predictor[, intersect(c("id", "cohort", "weight_bl"), names(frame))] <- 0L
  for (v in names(frame)) {
    observed <- frame[[v]][!is.na(frame[[v]])]
    if (length(unique(observed)) < 2L) {
      predictor[, v] <- 0L
      method[v] <- ""
    }
  }
  where <- is.na(frame)
  where[, protected] <- FALSE

  baseline_digest_frame <- baseline[, sort(names(baseline)), with = FALSE]
  imputation_digest_frame <- frame[, sort(names(frame)), drop = FALSE]
  cache_signature <- list(
    signature_version = "2026-09-08-v2-full-frame",
    cohort = cohort_name,
    participants = nrow(baseline),
    baseline_sha256 = sha256_object(baseline_digest_frame),
    imputation_frame_sha256 = sha256_object(imputation_digest_frame),
    imputation_plan_sha256 = sha256_object(list(method = method, predictor = predictor, where = where)),
    protected_fields = protected,
    m = config$mice_m,
    maxit = config$mice_maxit,
    seed = config$seed,
    mice_version = as.character(packageVersion("mice")),
    script_sha256 = mice_cache_script_sha256
  )
  cache_key <- sha256_object(cache_signature)
  cache_file <- file.path(
    mice_cache_root,
    paste0(tolower(cohort_name), "__m", config$mice_m, "__", cache_key, ".rds")
  )
  if (file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    cache_valid <- !is.null(cached) && identical(cached$signature, cache_signature) &&
      identical(cached$result$status, "success") &&
      identical(as.integer(cached$result$protection_violations), 0L) &&
      length(cached$result$completed) == config$mice_m
    if (cache_valid) {
      cache_protection <- vapply(cached$result$completed, function(d) {
        count_protected_changes(as.data.table(copy(d)), baseline, protected)
      }, integer(1L))
      completed_digest <- sha256_object(cached$result$completed)
      cache_valid <- all(cache_protection == 0L) &&
        identical(cached$result$completed_sha256, completed_digest)
    }
    if (cache_valid) {
      cached$result$cache_status <- "reused_exact_signature"
      cached$result$cache_file <- cache_file
      return(cached$result)
    }
  }

  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      mice::mice(
        frame,
        m = config$mice_m,
        maxit = config$mice_maxit,
        method = method,
        predictorMatrix = predictor,
        where = where,
        seed = config$seed + match(cohort_name, c("HRS", "ELSA")),
        printFlag = FALSE,
        remove.constant = TRUE,
        remove.collinear = TRUE
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(
      status = "failed", message = conditionMessage(fit), completed = list(),
      baseline = baseline, log = data.table(), warnings = collapse_messages(warnings),
      protection_violations = NA_integer_, cache_status = "not_created_mice_failed",
      cache_file = cache_file, cache_signature = cache_signature
    ))
  }

  completed <- mice::complete(fit, action = "all", include = FALSE)
  completed <- lapply(completed, canonicalize_completed_baseline)
  protected_compare <- intersect(protected, names(baseline))
  protection_violations <- sum(vapply(completed, function(d) {
    count_protected_changes(d, baseline, protected_compare)
  }, integer(1L)))
  logged <- if (is.null(fit$loggedEvents)) data.table() else as.data.table(fit$loggedEvents)
  status <- if (protection_violations > 0L) "failed_protected_field_changed" else "success"
  message <- if (protection_violations > 0L) {
    paste0(protection_violations, " protected baseline value(s) changed across completed data sets")
  } else {
    ""
  }
  result <- list(
    status = status, message = message, completed = completed, baseline = baseline,
    log = logged, warnings = collapse_messages(warnings),
    protection_violations = protection_violations,
    baseline_sha256 = cache_signature$baseline_sha256,
    imputation_frame_sha256 = cache_signature$imputation_frame_sha256,
    imputation_plan_sha256 = cache_signature$imputation_plan_sha256,
    completed_sha256 = sha256_object(completed),
    cache_status = "created",
    cache_file = cache_file,
    cache_signature = cache_signature
  )
  dir.create(mice_cache_root, recursive = TRUE, showWarnings = FALSE)
  cache_lock <- paste0(cache_file, ".lock")
  if (dir.create(cache_lock, showWarnings = FALSE)) {
    on.exit(unlink(cache_lock, recursive = TRUE, force = TRUE), add = TRUE)
    atomic_saveRDS(
      list(signature = cache_signature, result = result),
      cache_file,
      compress = "gzip"
    )
  } else {
    result$cache_status <- "created_in_memory_cache_lock_unavailable"
  }
  result
}

replace_imputed_covariates <- function(intervals, completed_baseline) {
  d <- copy(intervals)
  replacement <- intersect(baseline_impute_vars, names(completed_baseline))
  idx <- match(d$id, completed_baseline$id)
  if (anyNA(idx)) stop("Completed baseline could not be matched to all interval IDs")
  for (v in replacement) set(d, j = v, value = completed_baseline[[v]][idx])
  d
}

derive_model_variables <- function(d, cohort_name) {
  d[, `:=`(
    age5 = as.numeric(age_bl) / 5,
    bmi5 = as.numeric(bmi_bl) / 5,
    obesity_bl = fifelse(is.na(bmi_bl), NA_integer_, as.integer(as.numeric(bmi_bl) >= 30))
  )]
  d[, `:=`(
    female_factor = factor(as.integer(female_bl), levels = 0:1, labels = c("Male", "Female")),
    education_factor = factor(as.integer(education_level_bl), levels = 1:3,
                              labels = c("Below upper secondary", "Upper secondary/vocational", "Tertiary")),
    wealth_factor = factor(as.integer(wealth_quintile_bl), levels = 1:5,
                           labels = paste0("Q", 1:5)),
    obesity_factor = factor(obesity_bl, levels = 0:1, labels = c("BMI<30", "BMI>=30")),
    diabetes_factor = factor(as.integer(diabetes_bl), levels = 0:1, labels = c("No diabetes", "Diabetes"))
  )]
  race_levels <- if (cohort_name == "HRS") {
    c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other")
  } else {
    c("White", "Non-White")
  }
  d[, race_factor := factor(as.character(race_bl), levels = race_levels)]
  d[]
}

apply_scenario <- function(intervals, spec) {
  d <- intervals[from_state == spec$from_state]
  d[, `:=`(tstart_fit = as.numeric(tstart), tstop_fit = as.numeric(tstop))]

  mortality_coverage <- "not_applicable"
  if (spec$cohort == "ELSA" && spec$to_state == 4L) {
    cutoff <- config$elsa_death_cutoff
    d <- d[event_time < cutoff]
    d[, tstop_fit := pmin(next_time, cutoff) - exposure_entry_time]
    mortality_coverage <- paste0("limited_mortality_coverage_through_", floor(cutoff))
  } else if (spec$cohort == "ELSA") {
    mortality_coverage <- "nonfatal_followup_includes_W10"
  } else if (spec$to_state == 4L) {
    mortality_coverage <- "HRS_death_timing_from_R05"
  }

  if (spec$scenario == "landmark2") {
    d <- d[tstop_fit > 2]
    d[, tstart_fit := pmax(tstart_fit, 2)]
  } else if (spec$scenario == "observed_sleep") {
    d <- d[interval_starts_at_observed_sleep == 1L &
             (is.na(proxy_interview) | proxy_interview == 0L)]
  } else if (spec$scenario == "interval_le4") {
    d <- d[interval_years <= 4]
  } else if (spec$scenario == "lag_available") {
    d <- d[!is.na(insomnia_score_lag1_row)]
  } else if (spec$scenario == "previous_valid_available") {
    d <- d[!is.na(insomnia_score_prev_valid)]
  } else if (spec$scenario == "dynamic_delayed") {
    dynamic_levels <- c("persistent low", "incident high", "remission", "persistent high")
    d <- d[insomnia_assessment_number >= 2L & insomnia_dynamic_pattern %chin% dynamic_levels]
  }

  if (spec$cohort == "ELSA" && spec$to_state == 4L) {
    d[, event_target := as.integer(
      state_changed == 1L & to_state == spec$to_state & next_time <= config$elsa_death_cutoff
    )]
    d[, competing_exit_observed := as.integer(
      state_changed == 1L & to_state != spec$to_state & next_time <= config$elsa_death_cutoff
    )]
  } else {
    d[, event_target := as.integer(state_changed == 1L & to_state == spec$to_state)]
    d[, competing_exit_observed := as.integer(state_changed == 1L & event_target == 0L)]
  }
  d <- d[is.finite(tstart_fit) & is.finite(tstop_fit)]
  # survival::Surv/aeqSurv compares times with a square-root machine-epsilon
  # tolerance.  After landmark truncation, floating-point subtraction can leave
  # nominally positive intervals (~1e-13 years) that Surv correctly treats as
  # zero length.  Apply the same scale-aware tolerance before model fitting.
  unique_times <- unique(c(d$tstart_fit, d$tstop_fit))
  duration_tolerance <- sqrt(.Machine$double.eps) * max(1, mean(abs(unique_times)))
  d <- d[(tstop_fit - tstart_fit) > duration_tolerance]
  attr(d, "mortality_coverage") <- mortality_coverage
  d[]
}

exposure_rhs <- function(spec) {
  switch(
    spec$exposure_type,
    continuous = "exposure_value",
    binary = "exposure_value",
    categorical = "exposure_factor",
    dynamic = "exposure_dynamic",
    interaction = paste0("exposure_value * ", spec$interaction_var),
    stop("Unknown exposure type: ", spec$exposure_type)
  )
}

covariate_rhs <- function(spec) {
  terms <- model_terms[[spec$model_level]]
  if (is.null(terms)) stop("Unknown model level: ", spec$model_level)
  if (spec$analysis_id == "interaction_obesity") {
    terms <- c(setdiff(terms, "bmi5"), "obesity_factor")
  }
  if (spec$analysis_id == "interaction_diabetes") {
    # The interaction expansion already supplies the diabetes main effect.
    # Keeping diabetes_bl as a second coding of the same binary variable makes
    # the design matrix structurally rank deficient.
    terms <- c(setdiff(terms, "diabetes_bl"), "diabetes_factor")
  }
  if (spec$analysis_id %chin% c("depression_adjusted_score", "depression_adjusted_burden")) {
    terms <- c(terms, "depression_score_bl")
  }
  unique(terms)
}

prepare_model_data <- function(intervals, spec, completed_baseline = NULL) {
  d <- if (is.null(completed_baseline)) copy(intervals) else replace_imputed_covariates(intervals, completed_baseline)
  d <- derive_model_variables(d, spec$cohort)
  d <- apply_scenario(d, spec)
  mortality_coverage <- attr(d, "mortality_coverage")

  if (!spec$exposure_var %in% names(d)) {
    stop("Exposure field unavailable: ", spec$exposure_var)
  }
  value <- d[[spec$exposure_var]]
  if (spec$exposure_type %in% c("continuous", "binary", "interaction")) {
    d[, exposure_value := as.numeric(value)]
  } else if (spec$exposure_type == "categorical") {
    d[, exposure_factor := factor(as.integer(value), levels = 0:2,
                                  labels = c("0 symptoms", "1 symptom", ">=2 symptoms"))]
  } else if (spec$exposure_type == "dynamic") {
    d[, exposure_dynamic := factor(
      as.character(value),
      levels = c("persistent low", "incident high", "remission", "persistent high")
    )]
  }
  attr(d, "mortality_coverage") <- mortality_coverage
  d[]
}

build_formula_parts <- function(spec) {
  exposure <- exposure_rhs(spec)
  covariates <- covariate_rhs(spec)
  rhs <- paste(c(exposure, covariates), collapse = " + ")
  list(
    rhs = rhs,
    formula = as.formula(paste0(
      "Surv(tstart_fit, tstop_fit, event_target) ~ ", rhs, " + cluster(id)"
    )),
    matrix_formula = as.formula(paste0("~ ", rhs))
  )
}

complete_model_frame <- function(d, spec, formula_parts) {
  required <- unique(c(
    "id", "tstart_fit", "tstop_fit", "event_target",
    all.vars(formula_parts$matrix_formula),
    if (spec$weighted == 1L) "weight_bl" else character()
  ))
  required <- intersect(required, names(d))
  keep <- complete.cases(d[, ..required])
  if (spec$weighted == 1L) keep <- keep & is.finite(d$weight_bl) & d$weight_bl > 0
  d <- d[keep]
  if (spec$weighted == 1L && nrow(d)) {
    d[, model_weight := as.numeric(weight_bl) / mean(as.numeric(weight_bl))]
  } else {
    d[, model_weight := 1]
  }
  d[]
}

interaction_cell_check <- function(d, spec) {
  if (spec$exposure_type != "interaction") return(list(ok = TRUE, detail = "not an interaction"))
  modifier <- spec$interaction_var
  if (!modifier %in% names(d)) return(list(ok = FALSE, detail = paste("modifier unavailable:", modifier)))
  event_values <- d[event_target == 1L, get(modifier)]
  if (is.factor(d[[modifier]]) || is.character(d[[modifier]])) {
    counts <- table(event_values, useNA = "no")
    if (length(counts) < 2L || min(counts) < config$min_interaction_cell_events) {
      return(list(
        ok = FALSE,
        detail = paste0("target-event counts by modifier level: ", paste(names(counts), counts, sep = "=", collapse = ", "))
      ))
    }
  } else if (length(event_values) < config$min_events || !is.finite(sd(as.numeric(event_values))) || sd(as.numeric(event_values)) == 0) {
    return(list(ok = FALSE, detail = "modifier has inadequate variation among target events"))
  }
  list(ok = TRUE, detail = "interaction event-cell check passed")
}

preflight_model <- function(d, spec, formula_parts) {
  n_risk_rows_raw <- nrow(d)
  n_risk_ids_raw <- uniqueN(d$id)
  n_events_raw <- sum(d$event_target == 1L, na.rm = TRUE)
  d <- complete_model_frame(d, spec, formula_parts)
  n_events <- sum(d$event_target == 1L)
  n_competing <- sum(d$competing_exit_observed == 1L, na.rm = TRUE)
  matrix_result <- tryCatch(model.matrix(formula_parts$matrix_formula, data = d), error = function(e) e)
  if (inherits(matrix_result, "error")) {
    return(list(
      ok = FALSE, status = "not_estimable_design_matrix", reason = conditionMessage(matrix_result),
      data = d, n_risk_rows_raw = n_risk_rows_raw, n_risk_ids_raw = n_risk_ids_raw,
      n_events_raw = n_events_raw, n_risk_rows = nrow(d), n_risk_ids = uniqueN(d$id),
      n_events = n_events, n_competing = n_competing, n_parameters = NA_integer_,
      events_per_parameter = NA_real_, interaction_check = "not run: design matrix failed"
    ))
  }
  mm <- matrix_result
  if (ncol(mm) && colnames(mm)[1L] == "(Intercept)") mm <- mm[, -1L, drop = FALSE]
  n_parameters <- if (ncol(mm)) qr(mm)$rank else 0L
  epp <- if (n_parameters > 0L) n_events / n_parameters else NA_real_
  interaction_check <- interaction_cell_check(d, spec)
  reasons <- character()
  status <- "estimable"
  if (!nrow(d) || uniqueN(d$id) < 2L) {
    status <- "not_estimable_no_risk_set"
    reasons <- c(reasons, "fewer than two participants after required complete-case filtering")
  }
  if (n_events < config$min_events) {
    status <- "not_estimable_low_events"
    reasons <- c(reasons, paste0("target events ", n_events, " < minimum ", config$min_events))
  }
  if (n_parameters < 1L) {
    status <- "not_estimable_no_exposure_variation"
    reasons <- c(reasons, "design matrix has zero rank")
  }
  if (n_parameters > 0L && is.finite(config$min_events_per_parameter) &&
      epp < config$min_events_per_parameter) {
    status <- "not_estimable_low_events_per_parameter"
    reasons <- c(reasons, sprintf(
      "events per estimable parameter %.2f < %.2f", epp, config$min_events_per_parameter
    ))
  }
  if (!interaction_check$ok) {
    status <- "not_estimable_interaction_cells"
    reasons <- c(reasons, interaction_check$detail)
  }
  list(
    ok = identical(status, "estimable"), status = status,
    reason = collapse_messages(reasons), data = d,
    n_risk_rows_raw = n_risk_rows_raw, n_risk_ids_raw = n_risk_ids_raw,
    n_events_raw = n_events_raw, n_risk_rows = nrow(d), n_risk_ids = uniqueN(d$id),
    n_events = n_events, n_competing = n_competing,
    n_parameters = n_parameters, events_per_parameter = epp,
    interaction_check = interaction_check$detail
  )
}

classify_term <- function(term) {
  if (grepl(":", term, fixed = TRUE)) return("interaction")
  if (grepl("^exposure_(value|factor|dynamic)", term)) return("exposure")
  "covariate"
}

extract_cox_estimates <- function(fit) {
  beta <- coef(fit)
  if (!length(beta)) return(data.table())
  robust_var <- fit$var
  se <- sqrt(diag(robust_var))
  z <- beta / se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  data.table(
    term = names(beta),
    term_role = vapply(names(beta), classify_term, character(1)),
    estimate = as.numeric(beta),
    robust_se = as.numeric(se),
    within_variance = as.numeric(se^2),
    statistic = as.numeric(z),
    df = Inf,
    p_value = as.numeric(p),
    hazard_ratio = exp(as.numeric(beta)),
    conf_low = exp(as.numeric(beta) - qnorm(0.975) * as.numeric(se)),
    conf_high = exp(as.numeric(beta) + qnorm(0.975) * as.numeric(se))
  )
}

generalized_inverse <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  x <- (x + t(x)) / 2
  eig <- eigen(x, symmetric = TRUE)
  keep <- eig$values > max(abs(eig$values), 1) * tolerance
  if (!any(keep)) return(list(inverse = matrix(NA_real_, nrow(x), ncol(x)), rank = 0L))
  inverse <- eig$vectors[, keep, drop = FALSE] %*%
    diag(1 / eig$values[keep], nrow = sum(keep)) %*%
    t(eig$vectors[, keep, drop = FALSE])
  list(inverse = inverse, rank = sum(keep))
}

pool_joint_wald <- function(fits, coefficient_names, block_name) {
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (!length(fits)) return(data.table())
  available <- lapply(fits, function(f) intersect(coefficient_names, names(coef(f))))
  common <- Reduce(intersect, available)
  if (!length(common)) return(data.table())
  q <- do.call(rbind, lapply(fits, function(f) as.numeric(coef(f)[common])))
  if (is.null(dim(q))) q <- matrix(q, ncol = length(common))
  u <- lapply(fits, function(f) {
    idx <- match(common, names(coef(f)))
    f$var[idx, idx, drop = FALSE]
  })
  valid <- vapply(seq_along(u), function(i) {
    all(is.finite(q[i, ])) && all(is.finite(u[[i]]))
  }, logical(1))
  q <- q[valid, , drop = FALSE]
  u <- u[valid]
  m <- nrow(q)
  if (!m) return(data.table())
  qbar <- colMeans(q)
  ubar <- Reduce(`+`, u) / m
  between <- if (m > 1L) stats::cov(q) else matrix(0, nrow = length(common), ncol = length(common))
  if (length(common) == 1L) between <- matrix(between, 1L, 1L)
  total <- ubar + (1 + 1 / m) * between
  inverse <- generalized_inverse(total)
  if (inverse$rank < 1L || any(!is.finite(inverse$inverse))) return(data.table())
  statistic <- as.numeric(t(qbar) %*% inverse$inverse %*% qbar)
  data.table(
    term_block = block_name,
    terms = paste(common, collapse = " | "),
    statistic = statistic,
    df = inverse$rank,
    p_value = pchisq(statistic, df = inverse$rank, lower.tail = FALSE),
    m_used = m,
    method = if (m > 1L) {
      "Large-sample joint Wald test using Rubin total covariance"
    } else {
      "Cluster-robust joint Wald test"
    }
  )
}

make_joint_term_tests <- function(fits) {
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (!length(fits)) return(data.table())
  common_names <- Reduce(intersect, lapply(fits, function(f) names(coef(f))))
  exposure_main <- common_names[grepl("^exposure_(value|factor|dynamic)", common_names) & !grepl(":", common_names, fixed = TRUE)]
  interaction <- common_names[grepl(":", common_names, fixed = TRUE)]
  rbindlist(list(
    pool_joint_wald(fits, exposure_main, "exposure_overall"),
    pool_joint_wald(fits, interaction, "interaction_overall")
  ), use.names = TRUE, fill = TRUE)
}

extract_zph <- function(fit, imputation_index) {
  if (!config$run_zph) return(data.table())
  z <- tryCatch(survival::cox.zph(fit, transform = "km", global = TRUE), error = function(e) e)
  if (inherits(z, "error")) {
    return(data.table(
      imputation = imputation_index, ph_term = "<cox.zph failed>",
      chisq = NA_real_, df = NA_real_, p_value = NA_real_, message = conditionMessage(z)
    ))
  }
  out <- as.data.table(z$table, keep.rownames = "ph_term")
  setnames(out, old = intersect(c("chisq", "df", "p"), names(out)),
           new = c("chisq", "df", "p_value")[match(intersect(c("chisq", "df", "p"), names(out)), c("chisq", "df", "p"))])
  out[, `:=`(imputation = imputation_index, message = "")]
  setcolorder(out, c("imputation", "ph_term", "chisq", "df", "p_value", "message"))
  out[]
}

fit_one_cox <- function(intervals, spec, completed_baseline = NULL, imputation_index = 0L) {
  warnings <- character()
  prepared <- tryCatch(
    prepare_model_data(intervals, spec, completed_baseline),
    error = function(e) e
  )
  if (inherits(prepared, "error")) {
    return(list(
      status = "failed_data_preparation", message = conditionMessage(prepared),
      estimates = data.table(), ph = data.table(), fit = NULL,
      qc = data.table(
        imputation = imputation_index, status = "failed_data_preparation", reason = conditionMessage(prepared),
        n_risk_rows_raw = NA_integer_, n_risk_ids_raw = NA_integer_, n_events_raw = NA_integer_,
        n_risk_rows = NA_integer_, n_risk_ids = NA_integer_, n_events = NA_integer_,
        n_competing_exits = NA_integer_, n_parameters = NA_integer_, events_per_parameter = NA_real_,
        interaction_check = NA_character_, mortality_coverage = NA_character_
      )
    ))
  }
  mortality_coverage <- attr(prepared, "mortality_coverage")
  parts <- build_formula_parts(spec)
  preflight <- preflight_model(prepared, spec, parts)
  qc <- data.table(
    imputation = imputation_index,
    status = preflight$status,
    reason = preflight$reason,
    n_risk_rows_raw = preflight$n_risk_rows_raw,
    n_risk_ids_raw = preflight$n_risk_ids_raw,
    n_events_raw = preflight$n_events_raw,
    n_risk_rows = preflight$n_risk_rows,
    n_risk_ids = preflight$n_risk_ids,
    n_events = preflight$n_events,
    n_competing_exits = preflight$n_competing,
    n_parameters = preflight$n_parameters,
    events_per_parameter = preflight$events_per_parameter,
    interaction_check = preflight$interaction_check,
    mortality_coverage = mortality_coverage
  )
  if (!preflight$ok) {
    return(list(
      status = preflight$status, message = preflight$reason,
      estimates = data.table(), ph = data.table(), fit = NULL, qc = qc
    ))
  }

  cox_arguments <- list(
    formula = parts$formula,
    data = preflight$data,
    ties = "efron",
    singular.ok = TRUE,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )
  # Passing `preflight$data$model_weight` unevaluated through coxph() leaves a
  # call that model.frame later evaluates outside this local scope.  Supply an
  # already evaluated vector through do.call() so weighted fits are restartable
  # and do not fail with `object 'preflight' not found`.
  if (spec$weighted == 1L) cox_arguments$weights <- preflight$data$model_weight
  fit <- tryCatch(
    withCallingHandlers(
      do.call(survival::coxph, cox_arguments),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    qc[, `:=`(status = "failed_coxph", reason = conditionMessage(fit))]
    return(list(
      status = "failed_coxph", message = conditionMessage(fit),
      estimates = data.table(), ph = data.table(), fit = NULL, qc = qc
    ))
  }

  estimates <- extract_cox_estimates(fit)
  if (!nrow(estimates) || any(!is.finite(estimates$estimate)) || any(!is.finite(estimates$robust_se))) {
    status <- "failed_nonfinite_coefficients"
    message <- collapse_messages(c(warnings, "one or more coefficients or robust standard errors are non-finite"))
  } else {
    status <- if (length(warnings)) "success_with_warning" else "success"
    message <- collapse_messages(warnings)
  }
  status_value <- status
  message_value <- message
  qc[, `:=`(
    status = status_value,
    reason = message_value,
    model_n = fit$n,
    model_events = fit$nevent,
    concordance = unname(summary(fit)$concordance[1L])
  )]
  list(
    status = status, message = message, estimates = estimates,
    ph = extract_zph(fit, imputation_index),
    fit = fit,
    qc = qc
  )
}

pool_scalar_rubin <- function(q, u) {
  ok <- is.finite(q) & is.finite(u) & u >= 0
  q <- q[ok]
  u <- u[ok]
  m <- length(q)
  if (!m) return(list(m = 0L, estimate = NA_real_, se = NA_real_, df = NA_real_, p = NA_real_, low = NA_real_, high = NA_real_, within = NA_real_, between = NA_real_, fmi = NA_real_))
  qbar <- mean(q)
  ubar <- mean(u)
  between <- if (m > 1L) stats::var(q) else 0
  total <- ubar + (1 + 1 / m) * between
  se <- sqrt(total)
  if (!is.finite(between) || between <= .Machine$double.eps) {
    df <- Inf
    crit <- qnorm(0.975)
    fmi <- 0
  } else {
    relative <- (1 + 1 / m) * between / max(ubar, .Machine$double.eps)
    df <- (m - 1) * (1 + 1 / relative)^2
    crit <- qt(0.975, df = df)
    fmi <- (1 + 1 / m) * between / total
  }
  statistic <- qbar / se
  p <- if (is.finite(df)) 2 * pt(abs(statistic), df = df, lower.tail = FALSE) else 2 * pnorm(abs(statistic), lower.tail = FALSE)
  list(
    m = m, estimate = qbar, se = se, df = df, p = p,
    low = qbar - crit * se, high = qbar + crit * se,
    within = ubar, between = between, fmi = fmi
  )
}

pool_mi_estimates <- function(estimates) {
  if (!nrow(estimates)) return(data.table())
  estimates[, {
    pooled <- pool_scalar_rubin(estimate, within_variance)
    list(
      m_used = pooled$m,
      estimate = pooled$estimate,
      robust_se = pooled$se,
      within_variance = pooled$within,
      between_variance = pooled$between,
      fraction_missing_information = pooled$fmi,
      statistic = pooled$estimate / pooled$se,
      df = pooled$df,
      p_value = pooled$p,
      hazard_ratio = exp(pooled$estimate),
      conf_low = exp(pooled$low),
      conf_high = exp(pooled$high)
    )
  }, by = .(term, term_role)]
}

summarize_mi_ph <- function(ph) {
  if (!nrow(ph)) return(data.table())
  ph[ph_term != "<cox.zph failed>", .(
    imputations_with_test = sum(is.finite(p_value)),
    p_min = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
    p_median = if (any(is.finite(p_value))) median(p_value, na.rm = TRUE) else NA_real_,
    proportion_p_lt_0_05 = if (any(is.finite(p_value))) mean(p_value < 0.05, na.rm = TRUE) else NA_real_,
    interpretation = "Descriptive across imputations; Schoenfeld tests are not Rubin-pooled"
  ), by = ph_term]
}

spec_metadata <- function(spec) {
  data.table(
    spec_id = spec$spec_id,
    spec_hash_full = spec$spec_hash_full,
    module = spec$module,
    cohort = spec$cohort,
    transition_id = spec$transition_id,
    transition = spec$transition,
    from_state = spec$from_state,
    to_state = spec$to_state,
    death_transition = spec$death_transition,
    confirmatory_nonfatal = spec$confirmatory_nonfatal,
    analysis_family = spec$analysis_family,
    analysis_id = spec$analysis_id,
    analysis_label = spec$analysis_label,
    scenario = spec$scenario,
    exposure_var = spec$exposure_var,
    exposure_type = spec$exposure_type,
    interaction_var = spec$interaction_var,
    model_level = spec$model_level,
    basis = spec$basis,
    weighted = spec$weighted,
    covariate_set = spec$covariate_set,
    formula_rhs = spec$formula_rhs,
    estimation_role = spec$estimation_role,
    mortality_scope = spec$mortality_scope
  )
}

is_claim_bearing_primary <- function(x) {
  x$module == "main" &
    x$analysis_family == "main" &
    x$analysis_id == "primary_score_tv" &
    x$scenario == "standard" &
    x$exposure_type == "continuous" &
    x$model_level == "M4" &
    x$basis == "mice_pooled" &
    as.integer(x$weighted) == 0L &
    as.integer(x$confirmatory_nonfatal) == 1L
}

run_spec <- function(intervals, spec, mice_result) {
  metadata <- spec_metadata(spec)
  fit_results <- list()
  if (spec$basis == "mice_pooled") {
    if (!config$run_mice) {
      return(list(
        run = cbind(metadata, data.table(status = "not_run_mice_disabled", message = "GA_RUN_MICE is false")),
        estimates = data.table(), term_tests = data.table(), ph = data.table(), ph_summary = data.table(), qc = data.table(),
        individual_estimates = data.table(), objects = list()
      ))
    }
    if (is.null(mice_result) || mice_result$status != "success") {
      message <- if (is.null(mice_result)) "MICE result unavailable" else mice_result$message
      return(list(
        run = cbind(metadata, data.table(status = "failed_mice_unavailable", message = message)),
        estimates = data.table(), term_tests = data.table(), ph = data.table(), ph_summary = data.table(), qc = data.table(),
        individual_estimates = data.table(), objects = list()
      ))
    }
    fit_results <- lapply(seq_along(mice_result$completed), function(i) {
      fit_one_cox(intervals, spec, mice_result$completed[[i]], i)
    })
  } else {
    fit_results <- list(fit_one_cox(intervals, spec, NULL, 0L))
  }

  qc <- rbindlist(lapply(fit_results, `[[`, "qc"), use.names = TRUE, fill = TRUE)
  if (nrow(qc)) qc <- cbind(metadata[rep(1L, nrow(qc))], qc)
  individual <- rbindlist(lapply(seq_along(fit_results), function(i) {
    z <- fit_results[[i]]$estimates
    if (!nrow(z)) return(NULL)
    z[, imputation := if (spec$basis == "mice_pooled") i else 0L]
    z
  }), use.names = TRUE, fill = TRUE)
  ph <- rbindlist(lapply(fit_results, `[[`, "ph"), use.names = TRUE, fill = TRUE)

  success <- vapply(fit_results, function(z) z$status %in% c("success", "success_with_warning"), logical(1))
  successful_fits <- lapply(fit_results[success], `[[`, "fit")
  term_tests <- make_joint_term_tests(successful_fits)
  messages <- collapse_messages(vapply(fit_results, `[[`, character(1), "message"))
  if (spec$basis == "mice_pooled") {
    minimum_success <- config$mice_m
    if (sum(success) < minimum_success) {
      status <- "failed_insufficient_successful_imputations"
      estimates <- data.table()
    } else {
      estimates <- pool_mi_estimates(individual)
      status <- if (any(vapply(fit_results, `[[`, character(1), "status") == "success_with_warning")) {
        "success_with_warning"
      } else {
        "success"
      }
    }
    ph_summary <- summarize_mi_ph(ph)
  } else {
    status <- if (length(fit_results)) fit_results[[1L]]$status else "failed_no_fit"
    estimates <- copy(individual)
    if ("imputation" %in% names(estimates)) estimates[, imputation := NULL]
    if (nrow(estimates)) estimates[, `:=`(
      m_used = 1L,
      between_variance = 0,
      fraction_missing_information = 0
    )]
    ph_summary <- if (nrow(ph)) ph[, .(
      imputations_with_test = sum(is.finite(p_value)),
      p_min = if (any(is.finite(p_value))) min(p_value, na.rm = TRUE) else NA_real_,
      p_median = if (any(is.finite(p_value))) median(p_value, na.rm = TRUE) else NA_real_,
      proportion_p_lt_0_05 = if (any(is.finite(p_value))) mean(p_value < 0.05, na.rm = TRUE) else NA_real_,
      interpretation = "Complete-case Schoenfeld residual test"
    ), by = ph_term] else data.table()
  }

  if (nrow(estimates)) estimates <- cbind(metadata[rep(1L, nrow(estimates))], estimates)
  if (nrow(term_tests)) term_tests <- cbind(metadata[rep(1L, nrow(term_tests))], term_tests)
  if (nrow(ph)) ph <- cbind(metadata[rep(1L, nrow(ph))], ph)
  if (nrow(ph_summary)) ph_summary <- cbind(metadata[rep(1L, nrow(ph_summary))], ph_summary)
  claim_bearing <- isTRUE(is_claim_bearing_primary(spec))
  primary_eligible <- claim_bearing && identical(status, "success") &&
    spec$basis == "mice_pooled" && sum(success) == config$mice_m
  primary_reason <- if (!claim_bearing) {
    "not_claim_bearing_primary"
  } else if (identical(status, "success_with_warning")) {
    "withheld_from_primary_inference_due_to_fit_warning"
  } else if (!identical(status, "success")) {
    paste0("withheld_from_primary_inference_status_", status)
  } else if (sum(success) != config$mice_m) {
    paste0("withheld_from_primary_inference_only_", sum(success), "_of_", config$mice_m, "_imputations_valid")
  } else {
    "eligible_before_four_edge_omnibus_and_FDR_gate"
  }
  run_summary <- cbind(metadata, data.table(
    status = status,
    message = messages,
    imputations_requested = if (spec$basis == "mice_pooled") config$mice_m else 1L,
    imputations_successful = sum(success),
    minimum_events = finite_min(qc$n_events),
    maximum_parameters = finite_max(qc$n_parameters),
    minimum_events_per_parameter = finite_min(qc$events_per_parameter),
    claim_bearing_primary = as.integer(claim_bearing),
    primary_inference_eligible_before_family_gate = as.integer(primary_eligible),
    primary_inference_reason = primary_reason,
    mortality_coverage = if (spec$cohort == "ELSA" && spec$to_state == 4L) {
      paste0("limited_mortality_coverage_through_", floor(config$elsa_death_cutoff))
    } else if (spec$cohort == "ELSA") {
      "nonfatal_followup_includes_W10"
    } else {
      "HRS"
    }
  ))
  objects <- if (config$save_model_objects) lapply(fit_results, `[[`, "fit") else list()
  list(
    run = run_summary, estimates = estimates, term_tests = term_tests, ph = ph, ph_summary = ph_summary,
    qc = qc,
    individual_estimates = if (config$save_imputation_estimates) individual else data.table(),
    objects = objects
  )
}

select_specs <- function(catalog) {
  x <- copy(catalog)
  x[, selected := TRUE]
  x[!cohort %chin% config$cohorts, selected := FALSE]
  if (!"ALL" %chin% config$transitions) x[!transition_id %chin% config$transitions, selected := FALSE]
  if (!"all" %chin% config$modules) x[!module %chin% config$modules, selected := FALSE]
  if (config$tier == "spec") x[, selected := FALSE]
  if (config$tier == "core") x[tier_required != "core", selected := FALSE]
  x[weighted == 1L & !config$run_weighted, selected := FALSE]
  if (config$weighted_only) x[weighted == 0L, selected := FALSE]
  x[basis == "mice_pooled" & !config$run_mice, selected := FALSE]
  if (!"ALL" %chin% config$spec_ids) x[!spec_id %chin% config$spec_ids, selected := FALSE]
  x[, selection_reason := fifelse(selected, "selected", "not selected")]
  if (config$tier == "spec") {
    x[selected == FALSE, selection_reason := "specification-only tier"]
  } else {
    x[selected == FALSE & selection_reason == "not selected" & !cohort %chin% config$cohorts,
      selection_reason := "cohort excluded by configuration"]
    x[selected == FALSE & selection_reason == "not selected" &
        !"ALL" %chin% config$transitions & !transition_id %chin% config$transitions,
      selection_reason := "transition excluded by configuration"]
    x[selected == FALSE & selection_reason == "not selected" &
        !"all" %chin% config$modules & !module %chin% config$modules,
      selection_reason := "module excluded by configuration"]
    x[selected == FALSE & selection_reason == "not selected" &
        config$tier == "core" & tier_required != "core",
      selection_reason := "requires full tier"]
    x[selected == FALSE & selection_reason == "not selected" & weighted == 1L & !config$run_weighted,
      selection_reason := "baseline-weight sensitivity disabled"]
    if (config$weighted_only) {
      x[selected == FALSE & selection_reason == "not selected" & weighted == 0L,
        selection_reason := "weighted-only run"]
    }
    x[selected == FALSE & selection_reason == "not selected" & basis == "mice_pooled" & !config$run_mice,
      selection_reason := "MICE disabled"]
    if (!"ALL" %chin% config$spec_ids) {
      x[selected == FALSE & selection_reason == "not selected" & !spec_id %chin% config$spec_ids,
        selection_reason := "specification ID excluded by configuration"]
    }
  }
  x[]
}

add_bh_fdr <- function(estimates, runs, term_tests, alpha = 0.05) {
  out <- copy(estimates)
  joint_table <- copy(term_tests)
  joint_required <- c(
    "spec_id", "module", "cohort", "transition_id", "analysis_family", "analysis_id",
    "scenario", "exposure_type", "model_level", "basis", "weighted",
    "confirmatory_nonfatal", "term_block", "terms", "df", "m_used", "p_value"
  )
  for (name in setdiff(joint_required, names(joint_table))) joint_table[, (name) := NA]
  out[, `:=`(
    claim_bearing_primary = as.integer(
      is_claim_bearing_primary(out) & term_role == "exposure" & term == "exposure_value"
    ),
    confirmatory_family_id = NA_character_,
    fdr_family_expected = NA_integer_,
    fdr_family_observed = NA_integer_,
    fdr_family_complete = NA_integer_,
    joint_overall_gate_method = NA_character_,
    joint_overall_p_value = NA_real_,
    joint_overall_gate_pass = NA_integer_,
    primary_inference_eligible = 0L,
    primary_inference_reason = "not_claim_bearing_primary",
    p_fdr_bh_confirmatory = NA_real_,
    fdr_family_size = NA_integer_
  )]

  expected_edges <- transition_catalog[confirmatory_nonfatal == 1L, transition_id]
  expected_n <- length(expected_edges)
  candidate_cohorts <- sort(unique(c(
    as.character(runs[claim_bearing_primary == 1L, cohort]),
    as.character(out[claim_bearing_primary == 1L, cohort])
  )))
  candidate_cohorts <- candidate_cohorts[!is.na(candidate_cohorts) & nzchar(candidate_cohorts)]
  audits <- list()

  for (cohort_name in candidate_cohorts) {
    family_id <- paste0(cohort_name, "__primary_M4_mice_score__four_nonfatal_edges")
    family_runs <- runs[
      cohort == cohort_name & claim_bearing_primary == 1L & transition_id %chin% expected_edges
    ]
    family_estimates <- out[
      cohort == cohort_name & claim_bearing_primary == 1L & transition_id %chin% expected_edges &
        term_role == "exposure" & term == "exposure_value"
    ]
    family_joint <- joint_table[
      cohort == cohort_name & transition_id %chin% expected_edges &
        module == "main" & analysis_family == "main" & analysis_id == "primary_score_tv" &
        scenario == "standard" & exposure_type == "continuous" &
        model_level == "M4" & basis == "mice_pooled" & weighted == 0L &
        confirmatory_nonfatal == 1L & term_block == "exposure_overall"
    ]

    runs_complete <- nrow(family_runs) == expected_n && uniqueN(family_runs$transition_id) == expected_n
    estimates_complete <- nrow(family_estimates) == expected_n &&
      uniqueN(family_estimates$transition_id) == expected_n && all(is.finite(family_estimates$p_value)) &&
      all(family_estimates$m_used == config$mice_m)
    joint_complete <- nrow(family_joint) == expected_n &&
      uniqueN(family_joint$transition_id) == expected_n && all(is.finite(family_joint$p_value)) &&
      all(family_joint$df == 1) && all(family_joint$terms == "exposure_value") &&
      all(family_joint$m_used == config$mice_m)
    matching_spec_sets <- runs_complete && estimates_complete && joint_complete &&
      setequal(family_runs$spec_id, family_estimates$spec_id) &&
      setequal(family_runs$spec_id, family_joint$spec_id)
    clean_runs <- runs_complete && all(family_runs$status == "success") &&
      all(family_runs$imputations_successful == config$mice_m) &&
      all(family_runs$primary_inference_eligible_before_family_gate == 1L)
    family_complete <- runs_complete && estimates_complete && joint_complete && matching_spec_sets && clean_runs
    joint_p <- if (estimates_complete) {
      min(1, expected_n * min(family_estimates$p_value))
    } else {
      NA_real_
    }
    gate_pass <- family_complete && is.finite(joint_p) && joint_p <= alpha
    reason <- if (!runs_complete) {
      paste0("withheld_family_runs_observed_", uniqueN(family_runs$transition_id), "_expected_", expected_n)
    } else if (!clean_runs) {
      "withheld_not_all_four_models_completed_cleanly_with_all_imputations"
    } else if (!estimates_complete) {
      paste0("withheld_family_estimates_observed_", uniqueN(family_estimates$transition_id), "_expected_", expected_n)
    } else if (!joint_complete) {
      paste0("withheld_joint_tests_observed_", uniqueN(family_joint$transition_id), "_expected_", expected_n)
    } else if (!matching_spec_sets) {
      "withheld_run_estimate_joint_spec_sets_do_not_match"
    } else if (!gate_pass) {
      "withheld_four_edge_bonferroni_omnibus_gate_not_significant"
    } else {
      "eligible_four_edge_family_complete_and_omnibus_gate_passed"
    }

    target_idx <- which(
      out$cohort == cohort_name & out$claim_bearing_primary == 1L &
        out$transition_id %chin% expected_edges & out$term_role == "exposure" &
        out$term == "exposure_value"
    )
    if (length(target_idx)) {
      out[target_idx, `:=`(
        confirmatory_family_id = family_id,
        fdr_family_expected = expected_n,
        fdr_family_observed = uniqueN(family_estimates$transition_id),
        fdr_family_complete = as.integer(family_complete),
        joint_overall_gate_method = "Bonferroni minimum-p omnibus across four edge-specific scalar Rubin tests; joint-Wald rows required as structural audit",
        joint_overall_p_value = joint_p,
        joint_overall_gate_pass = as.integer(gate_pass),
        primary_inference_eligible = as.integer(gate_pass),
        primary_inference_reason = reason,
        fdr_family_size = if (gate_pass) expected_n else NA_integer_
      )]
      if (gate_pass) {
        adjusted <- p.adjust(family_estimates$p_value, method = "BH")
        adjustment <- family_estimates[, .(spec_id, term)]
        adjustment[, p_fdr_bh_confirmatory := adjusted]
        out[adjustment, on = .(spec_id, term),
            p_fdr_bh_confirmatory := i.p_fdr_bh_confirmatory]
      }
    }

    audits[[length(audits) + 1L]] <- data.table(
      cohort = cohort_name,
      confirmatory_family_id = family_id,
      expected_edges = paste(expected_edges, collapse = ","),
      expected_family_size = expected_n,
      run_edges_observed = uniqueN(family_runs$transition_id),
      estimate_edges_observed = uniqueN(family_estimates$transition_id),
      joint_test_edges_observed = uniqueN(family_joint$transition_id),
      clean_models = as.integer(clean_runs),
      family_complete = as.integer(family_complete),
      joint_overall_gate_method = "Bonferroni minimum-p omnibus across four edge-specific scalar Rubin tests; joint-Wald rows required as structural audit",
      joint_overall_p_value = joint_p,
      alpha = alpha,
      joint_overall_gate_pass = as.integer(gate_pass),
      bh_applied = as.integer(gate_pass),
      inference_reason = reason
    )
  }
  list(
    estimates = out[],
    audit = bind_or_empty(audits, c(
      "cohort", "confirmatory_family_id", "expected_edges", "expected_family_size",
      "run_edges_observed", "estimate_edges_observed", "joint_test_edges_observed",
      "clean_models", "family_complete", "joint_overall_gate_method", "joint_overall_p_value",
      "alpha", "joint_overall_gate_pass", "bh_applied", "inference_reason"
    ))
  )
}

make_run_metadata <- function(catalog_fingerprint = NA_character_, protocol_catalog_fingerprint = NA_character_) {
  metadata <- data.table(
    field = c(
      "analysis_role", "transition_time_assumption", "competing_exit_rule", "S0_to_S3_rule",
      "independent_unit", "variance", "imputation", "multiplicity", "ELSA_mortality",
      "direct_sleep_sensitivity", "dynamic_pattern_entry",
      "tier", "modules", "mice_m", "mice_maxit", "weighted_sensitivity", "weighted_only", "R_version",
      "survival_version", "mice_version", "data.table_version"
    ),
    value = c(
      "Semiparametric visit-detection-time sensitivity analysis; not the primary panel Markov model",
      "Disease transition is approximated at the destination visit date",
      "Another exit from the same origin state is censored at its observed time",
      "Observed S0->S3 jumps are competing censoring exits and never target Cox events",
      "Participant; repeated risk intervals use cluster(id)",
      "Cluster-robust sandwich standard errors",
      paste0(config$mice_m, " baseline-covariate imputations with scalar Rubin pooling; all ", config$mice_m, " fits required; exposure/state/time/cohort/weight protected"),
      "Claim-bearing family only: unweighted M4 mice-pooled continuous score; four nonfatal edges per cohort; BH only after a complete clean family and conservative four-edge omnibus gate",
      paste0("Death-risk intervals truncated at ", config$elsa_death_cutoff, "; estimates labelled limited mortality coverage; W10 retained only for nonfatal models"),
      "Require a complete sleep assessment and exclude recorded proxy interviews; ELSA has no populated proxy flag, so completeness is the available directness criterion",
      "Dynamic-pattern risk time begins only at assessment number 2 or later with a valid adjacent-assessment pattern",
      config$tier,
      paste(config$modules, collapse = ","),
      as.character(config$mice_m),
      as.character(config$mice_maxit),
      as.character(config$run_weighted),
      as.character(config$weighted_only),
      as.character(getRversion()),
      as.character(packageVersion("survival")),
      as.character(packageVersion("mice")),
      as.character(packageVersion("data.table"))
    )
  )
  extra <- data.table(
    field = c(
      "run_label_requested", "run_label_fingerprinted", "run_instance_id", "run_fingerprint",
      "script_sha256", "selected_catalog_fingerprint", "protocol_catalog_fingerprint", "cohorts", "transitions", "seed",
      "min_events", "min_events_per_parameter", "ELSA_death_cutoff_numeric",
      "resume_checkpoints", "publish_canonical", "force_publish_canonical", "spec_ids",
      "mice_cache_script_sha256"
    ),
    value = c(
      run_label_requested, run_label, run_instance_id, run_fingerprint,
      script_sha256, catalog_fingerprint, protocol_catalog_fingerprint, paste(sort(config$cohorts), collapse = ","),
      paste(sort(config$transitions), collapse = ","), as.character(config$seed),
      as.character(config$min_events), as.character(config$min_events_per_parameter),
      as.character(config$elsa_death_cutoff), as.character(config$resume_checkpoints),
      as.character(config$publish_canonical), as.character(config$force_publish_canonical),
      paste(sort(config$spec_ids), collapse = ","), mice_cache_script_sha256
    )
  )
  input_rows <- input_manifest[, .(
    field = paste0("input_", cohort, "_sha256"),
    value = paste0(sha256, " | bytes=", bytes, " | path=", path)
  )]
  rbindlist(list(metadata, extra, input_rows), use.names = TRUE, fill = TRUE)
}

make_model_adjustment_dictionary <- function() {
  data.table(
    model_level = names(model_terms),
    description = c(
      "M1: age and sex",
      "M2: M1 plus race/ethnicity, comparable education, wealth quintile and marital status",
      "M3: M2 plus smoking, drinking, physical inactivity and BMI",
      "M4 primary: M3 plus diabetes, lung disease and frequent pain"
    ),
    covariate_terms = vapply(model_terms, function(x) paste(x, collapse = " + "), character(1)),
    reporting_note = c(
      "Age hazard ratio is per 5 years",
      "Categorical covariates retain fixed reference levels within cohort",
      "BMI hazard ratio is per 5 kg/m2",
      "Harmonized CES-D is omitted because exclusion of a sleep item cannot be guaranteed; separately specified depression-adjusted sensitivity models add it"
    )
  )
}

build_edge_registry <- function(intervals_by_cohort, catalog) {
  rows <- list()
  for (cohort_name in names(intervals_by_cohort)) {
    reference_specs <- catalog[
      cohort == cohort_name & analysis_id == "primary_score_tv" & model_level == "M4" &
        basis == "complete_case" & weighted == 0L
    ]
    for (i in seq_len(nrow(transition_catalog))) {
      tr <- transition_catalog[i]
      spec <- reference_specs[transition_id == tr$transition_id]
      if (!nrow(spec)) next
      spec <- spec[1L]
      prepared <- tryCatch(
        prepare_model_data(intervals_by_cohort[[cohort_name]], spec, NULL),
        error = function(e) e
      )
      if (inherits(prepared, "error")) {
        rows[[length(rows) + 1L]] <- data.table(
          cohort = cohort_name, transition_id = tr$transition_id, transition = tr$transition,
          from_state = tr$from_state, to_state = tr$to_state,
          reference_spec_id = spec$spec_id,
          risk = NA_integer_, risk_intervals = NA_integer_, risk_participants = NA_integer_,
          target_events = NA_integer_, competing_exits = NA_integer_, n_parameters = NA_integer_,
          events_per_parameter = NA_real_, status = "failed_data_preparation",
          reason = conditionMessage(prepared), mortality_coverage = NA_character_
        )
        next
      }
      mortality_coverage <- attr(prepared, "mortality_coverage")
      preflight <- preflight_model(prepared, spec, build_formula_parts(spec))
      rows[[length(rows) + 1L]] <- data.table(
        cohort = cohort_name, transition_id = tr$transition_id, transition = tr$transition,
        from_state = tr$from_state, to_state = tr$to_state,
        reference_spec_id = spec$spec_id,
        risk = preflight$n_risk_rows,
        risk_intervals = preflight$n_risk_rows,
        risk_participants = preflight$n_risk_ids,
        target_events = preflight$n_events,
        competing_exits = preflight$n_competing,
        n_parameters = preflight$n_parameters,
        events_per_parameter = preflight$events_per_parameter,
        status = preflight$status,
        reason = preflight$reason,
        mortality_coverage = mortality_coverage
      )
    }
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

make_canonical_cox_diagnostics <- function(edge_registry, runs, ph) {
  out <- copy(edge_registry)
  if (!nrow(out)) return(out)
  run_status <- if (nrow(runs) && all(c("spec_id", "status", "message") %in% names(runs))) {
    unique(runs[, .(
      reference_spec_id = spec_id,
      run_status = as.character(status),
      run_reason = as.character(message)
    )], by = "reference_spec_id")
  } else {
    data.table(reference_spec_id = character(), run_status = character(), run_reason = character())
  }
  out <- merge(out, run_status, by = "reference_spec_id", all.x = TRUE, sort = FALSE)

  global_ph <- if (nrow(ph) && all(c("spec_id", "ph_term", "p_value") %in% names(ph))) {
    ph[ph_term == "GLOBAL", .(
      ph_global_p_value = finite_min(p_value)
    ), by = .(reference_spec_id = spec_id)]
  } else {
    data.table(reference_spec_id = character(), ph_global_p_value = numeric())
  }
  out <- merge(out, global_ph, by = "reference_spec_id", all.x = TRUE, sort = FALSE)
  out[, preflight_status := as.character(status)]
  out[, preflight_reason := as.character(reason)]
  out[, status := fcase(
    run_status %chin% c("success", "success_with_warning"), "completed",
    grepl("^not_estimable", run_status), "not_estimable",
    grepl("^failed", run_status), "failed",
    grepl("^not_estimable", preflight_status), "not_estimable",
    grepl("^failed", preflight_status), "failed",
    default = "not_estimable"
  )]
  out[, reason := fcase(
    status == "completed", "",
    !is.na(run_reason) & nzchar(trimws(run_reason)), run_reason,
    !is.na(preflight_reason) & nzchar(trimws(preflight_reason)), preflight_reason,
    default = "Reference M4 complete-case diagnostic model was not selected or not run"
  )]
  out[, `:=`(
    n_events = as.numeric(target_events),
    event_parameter_risk = fcase(
      is.finite(events_per_parameter) & events_per_parameter >= config$min_events_per_parameter,
      "adequate",
      is.finite(events_per_parameter), "sparse_events",
      default = "not_quantifiable"
    ),
    mortality_note = mortality_coverage
  )]
  setcolorder(out, c(
    "cohort", "transition_id", "transition", "from_state", "to_state",
    "reference_spec_id", "status", "reason", "risk", "risk_intervals",
    "risk_participants", "n_events", "target_events", "competing_exits",
    "n_parameters", "events_per_parameter", "event_parameter_risk",
    "ph_global_p_value", "mortality_note", "preflight_status", "preflight_reason",
    "run_status", "run_reason"
  ))
  out[]
}

checkpoint_result_digest <- function(result) {
  keep <- intersect(
    c("run", "estimates", "term_tests", "ph", "ph_summary", "qc", "individual_estimates"),
    names(result)
  )
  sha256_object(result[keep])
}

checkpoint_file_for <- function(spec) {
  file.path(checkpoint_root, paste0(spec$spec_id, ".rds"))
}

load_spec_checkpoint <- function(spec) {
  if (!config$resume_checkpoints) return(NULL)
  path <- checkpoint_file_for(spec)
  if (!file.exists(path)) return(NULL)
  checkpoint <- tryCatch(readRDS(path), error = function(e) NULL)
  valid <- !is.null(checkpoint) &&
    identical(checkpoint$checkpoint_schema, "cox-spec-checkpoint-v1") &&
    identical(checkpoint$run_fingerprint, run_fingerprint) &&
    identical(checkpoint$spec_id, spec$spec_id) &&
    identical(checkpoint$spec_hash_full, spec$spec_hash_full) &&
    is.list(checkpoint$result) &&
    identical(checkpoint$result_sha256, checkpoint_result_digest(checkpoint$result))
  if (!valid) return(NULL)
  checkpoint$result
}

save_spec_checkpoint <- function(spec, result) {
  dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)
  checkpoint <- list(
    checkpoint_schema = "cox-spec-checkpoint-v1",
    run_fingerprint = run_fingerprint,
    spec_id = spec$spec_id,
    spec_hash_full = spec$spec_hash_full,
    result_sha256 = checkpoint_result_digest(result),
    written_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    result = result
  )
  atomic_saveRDS(checkpoint, checkpoint_file_for(spec), compress = "gzip")
  invisible(checkpoint)
}

write_run_state <- function(status, expected_specs = NA_integer_, checkpoint_specs = NA_integer_,
                            detail = "", catalog_fingerprint = NA_character_) {
  atomic_fwrite(data.table(
    status = status,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    run_fingerprint = run_fingerprint,
    catalog_fingerprint = catalog_fingerprint,
    run_label = run_label,
    run_instance_id = run_instance_id,
    expected_specifications = expected_specs,
    checkpoint_specifications = checkpoint_specs,
    detail = detail
  ), file.path(result_root, "RUN_STATE.csv"))
}

canonical_publish_eligibility <- function(catalog, runs, fdr_audit) {
  expected_edges <- transition_catalog$transition_id
  configured_edges <- if ("ALL" %chin% config$transitions) expected_edges else config$transitions
  scope_ok <- setequal(config$cohorts, c("HRS", "ELSA")) &&
    setequal(configured_edges, expected_edges) &&
    ("all" %chin% config$modules || "main" %chin% config$modules) &&
    config$run_mice && !config$weighted_only
  selected_primary <- catalog[selected == TRUE & is_claim_bearing_primary(catalog)]
  run_primary <- runs[claim_bearing_primary == 1L]
  family_ok <- nrow(fdr_audit) == 2L && setequal(fdr_audit$cohort, c("HRS", "ELSA")) &&
    all(fdr_audit$family_complete == 1L)
  primary_ok <- nrow(selected_primary) == 8L && uniqueN(selected_primary$spec_id) == 8L &&
    nrow(run_primary) == 8L && uniqueN(run_primary$spec_id) == 8L &&
    all(run_primary$status == "success") && all(run_primary$imputations_successful == config$mice_m)
  eligible <- scope_ok && family_ok && primary_ok
  reason <- if (!scope_ok) {
    "not_published_incomplete_cohort_transition_or_module_scope"
  } else if (!primary_ok) {
    "not_published_primary_eight_cohort_edge_models_not_complete_and_clean"
  } else if (!family_ok) {
    "not_published_primary_four_edge_families_not_structurally_complete"
  } else {
    "eligible_for_explicit_canonical_publish"
  }
  list(eligible = eligible, reason = reason)
}

publish_canonical_outputs <- function(catalog, runs, estimates, diagnostics, fdr_audit,
                                      completion_manifest) {
  eligibility <- canonical_publish_eligibility(catalog, runs, fdr_audit)
  status <- data.table(
    requested = as.integer(config$publish_canonical),
    eligible = as.integer(eligibility$eligible),
    published = 0L,
    run_fingerprint = run_fingerprint,
    reason = if (!config$publish_canonical) "not_requested_set_GA_COX_PUBLISH_CANONICAL_true" else eligibility$reason,
    canonical_bundle = NA_character_
  )
  if (!config$publish_canonical || !eligibility$eligible) return(status)

  canonical_root <- file.path(output_root, "results")
  bundle_root <- file.path(canonical_root, "cox_canonical_bundles", run_fingerprint, run_instance_id)
  dir.create(bundle_root, recursive = TRUE, showWarnings = FALSE)
  atomic_fwrite(estimates, file.path(bundle_root, "cox_transition_estimates.csv"))
  atomic_fwrite(diagnostics, file.path(bundle_root, "cox_model_diagnostics.csv"))
  atomic_fwrite(fdr_audit, file.path(bundle_root, "cox_primary_fdr_gate.csv"))
  atomic_fwrite(completion_manifest, file.path(bundle_root, "COMPLETION_MANIFEST.csv"))

  lock_dir <- file.path(canonical_root, ".cox_canonical_publish.lock")
  if (!dir.create(lock_dir, showWarnings = FALSE)) {
    status[, reason := "not_published_canonical_lock_unavailable"]
    return(status)
  }
  on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)
  manifest_path <- file.path(canonical_root, "cox_canonical_manifest.rds")
  existing <- if (file.exists(manifest_path)) tryCatch(readRDS(manifest_path), error = function(e) NULL) else NULL
  if (!is.null(existing) && !identical(existing$run_fingerprint, run_fingerprint) &&
      !config$force_publish_canonical) {
    status[, reason := "not_published_existing_canonical_has_different_fingerprint_use_explicit_force"]
    return(status)
  }

  atomic_fwrite(estimates, file.path(canonical_root, "cox_transition_estimates.csv"))
  atomic_fwrite(diagnostics, file.path(canonical_root, "cox_model_diagnostics.csv"))
  manifest <- list(
    manifest_schema = "cox-canonical-v1",
    ready = TRUE,
    run_fingerprint = run_fingerprint,
    run_instance_id = run_instance_id,
    bundle_root = normalizePath(bundle_root, winslash = "/", mustWork = FALSE),
    estimates_sha256 = sha256_file(file.path(canonical_root, "cox_transition_estimates.csv")),
    diagnostics_sha256 = sha256_file(file.path(canonical_root, "cox_model_diagnostics.csv")),
    published_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  atomic_saveRDS(manifest, manifest_path, compress = "gzip")
  atomic_fwrite(as.data.table(manifest), file.path(canonical_root, "cox_canonical_READY.csv"))
  status[, `:=`(
    published = 1L,
    reason = "published_complete_fingerprinted_canonical_bundle",
    canonical_bundle = normalizePath(bundle_root, winslash = "/", mustWork = FALSE)
  )]
  status
}

publish_compatibility_catalog <- function(full_catalog, protocol_catalog_fingerprint) {
  compatibility_path <- file.path(catalog_root, "transition_cox_model_specifications.csv")
  version_root <- file.path(catalog_root, "catalog_versions", protocol_catalog_fingerprint)
  dir.create(version_root, recursive = TRUE, showWarnings = FALSE)
  version_path <- file.path(version_root, "transition_cox_model_specifications.csv")

  status <- data.table(
    protocol_catalog_fingerprint = protocol_catalog_fingerprint,
    versioned_catalog = normalizePath(version_path, winslash = "/", mustWork = FALSE),
    compatibility_catalog = normalizePath(compatibility_path, winslash = "/", mustWork = FALSE),
    compatibility_status = "not_checked"
  )
  lock_dir <- file.path(catalog_root, ".cox_catalog_publish.lock")
  if (!dir.create(lock_dir, showWarnings = FALSE)) {
    status[, compatibility_status := "lock_unavailable_run_specific_catalog_only"]
    return(status)
  }
  on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)
  if (!file.exists(version_path)) atomic_fwrite(full_catalog, version_path)

  if (!file.exists(compatibility_path)) {
    atomic_fwrite(full_catalog, compatibility_path)
    status[, compatibility_status := "created_compatibility_catalog"]
    return(status)
  }
  existing <- tryCatch(fread(compatibility_path), error = function(e) data.table())
  same_catalog <- nrow(existing) == nrow(full_catalog) &&
    "spec_hash_full" %in% names(existing) &&
    setequal(existing$spec_hash_full, full_catalog$spec_hash_full)
  if (same_catalog) {
    status[, compatibility_status := "existing_compatibility_catalog_matches"]
  } else if (config$force_publish_canonical) {
    atomic_fwrite(full_catalog, compatibility_path)
    status[, compatibility_status := "replaced_by_explicit_force"]
  } else {
    status[, compatibility_status := "different_existing_catalog_not_overwritten_versioned_catalog_written"]
  }
  status
}

main <- function() {
  full_catalog <- make_spec_catalog()
  specification_qc <- validate_spec_catalog(full_catalog)
  protocol_catalog_fingerprint <- sha256_data_table(full_catalog, "spec_hash_full")
  catalog <- select_specs(copy(full_catalog))
  catalog_fingerprint <- sha256_data_table(catalog, "spec_hash_full")
  run_fingerprint <<- sha256_named_fields(list(
    base_run_fingerprint = base_run_fingerprint,
    selected_catalog_fingerprint = catalog_fingerprint
  ))
  run_label <<- paste0(run_label_requested, "__", substr(run_fingerprint, 1L, 16L))
  result_root <<- file.path(catalog_root, run_label, run_instance_id)
  qc_root <<- file.path(output_root, "qc", "cox_sensitivity", run_label, run_instance_id)
  checkpoint_root <<- file.path(result_root, "checkpoints")

  dir.create(catalog_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)
  compatibility_catalog_status <- publish_compatibility_catalog(
    full_catalog, protocol_catalog_fingerprint
  )
  run_completed <- FALSE
  on.exit({
    if (!run_completed) {
      try(write_run_state(
        "incomplete_or_failed", expected_specs = catalog[selected == TRUE, .N],
        checkpoint_specs = length(list.files(checkpoint_root, pattern = "\\.rds$")),
        detail = "Process exited before the completion manifest was finalized",
        catalog_fingerprint = catalog_fingerprint
      ), silent = TRUE)
    }
  }, add = TRUE)

  atomic_fwrite(catalog, file.path(result_root, "transition_cox_model_specifications.csv"))
  atomic_fwrite(compatibility_catalog_status, file.path(result_root, "compatibility_catalog_status.csv"))
  atomic_fwrite(make_model_adjustment_dictionary(), file.path(result_root, "transition_cox_adjustment_sets.csv"))
  atomic_fwrite(specification_qc, file.path(result_root, "transition_cox_specification_qc.csv"))
  atomic_fwrite(input_manifest, file.path(result_root, "input_file_manifest.csv"))
  atomic_fwrite(make_run_metadata(catalog_fingerprint, protocol_catalog_fingerprint), file.path(result_root, "run_metadata.csv"))

  selected <- catalog[selected == TRUE]
  write_run_state(
    "running", expected_specs = nrow(selected), checkpoint_specs = 0L,
    detail = "Run-specific outputs are isolated by analytic fingerprint and execution instance",
    catalog_fingerprint = catalog_fingerprint
  )
  if (!nrow(selected)) {
    completion_manifest <- data.table(
      status = "completed_specification_only", run_fingerprint = run_fingerprint,
      catalog_fingerprint = catalog_fingerprint,
      protocol_catalog_fingerprint = protocol_catalog_fingerprint, expected_specifications = 0L,
      checkpoint_specifications = 0L, completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    )
    atomic_fwrite(completion_manifest, file.path(result_root, "COMPLETION_MANIFEST.csv"))
    write_run_state(
      "completed_specification_only", 0L, 0L, "No model was selected by the spec tier",
      catalog_fingerprint
    )
    run_completed <- TRUE
    message("No models selected. The complete specification catalog was written to: ", result_root)
    return(invisible(NULL))
  }

  intervals_by_cohort <- list()
  input_qc <- list()
  for (cohort_name in unique(selected$cohort)) {
    path <- file.path(input_root, paste0(tolower(cohort_name), "_multistate_intervals.rds"))
    if (!file.exists(path)) stop("Missing R/05 interval input: ", path)
    x <- normalize_intervals(readRDS(path), cohort_name)
    intervals_by_cohort[[cohort_name]] <- x
    input_qc[[length(input_qc) + 1L]] <- data.table(
      cohort = cohort_name,
      rows = nrow(x),
      participants = uniqueN(x$id),
      nonpositive_intervals = x[!is.finite(tstop - tstart) | tstop <= tstart, .N],
      reverse_transitions = x[to_state < from_state & to_state != 4L, .N],
      direct_S0_S3_observations = x[from_state == 0L & to_state == 3L & state_changed == 1L, .N],
      direct_S0_S3_target_specifications = catalog[selected == TRUE & cohort == cohort_name & from_state == 0L & to_state == 3L, .N]
    )
  }
  input_manifest_after_read <- build_input_manifest(unique(selected$cohort))
  input_manifest_check <- merge(
    input_manifest[cohort %chin% unique(selected$cohort)],
    input_manifest_after_read,
    by = "cohort", suffixes = c("_before", "_after"), all = TRUE, sort = TRUE
  )
  input_manifest_check[, digest_match := as.integer(
    same_with_na(sha256_before, sha256_after) &
      same_with_na(bytes_before, bytes_after) &
      exists_before == exists_after
  )]
  atomic_fwrite(input_manifest_check, file.path(qc_root, "input_file_fingerprint_check.csv"))
  if (input_manifest_check[digest_match != 1L | is.na(digest_match), .N]) {
    stop("One or more interval inputs changed between initial fingerprinting and read; rerun with a stable input snapshot")
  }
  input_qc <- rbindlist(input_qc, use.names = TRUE, fill = TRUE)
  atomic_fwrite(input_qc, file.path(qc_root, "input_interval_qc.csv"))
  if (input_qc[nonpositive_intervals > 0L | reverse_transitions > 0L | direct_S0_S3_target_specifications > 0L, .N]) {
    stop("Input interval QC failed; inspect input_interval_qc.csv")
  }

  edge_registry <- build_edge_registry(intervals_by_cohort, catalog)
  atomic_fwrite(edge_registry, file.path(result_root, "transition_cox_edge_registry.csv"))
  if (length(intervals_by_cohort) == 2L && nrow(edge_registry) != 16L) {
    stop("The edge registry must contain all 16 cohort-by-edge records")
  }

  mice_by_cohort <- list()
  mice_diagnostics <- list()
  if (config$run_mice && selected[basis == "mice_pooled", .N]) {
    for (cohort_name in unique(selected[basis == "mice_pooled", cohort])) {
      message("Creating ", config$mice_m, " baseline-covariate imputations for ", cohort_name, "...")
      mi <- run_baseline_mice(intervals_by_cohort[[cohort_name]], cohort_name)
      mice_by_cohort[[cohort_name]] <- mi
      mice_diagnostics[[length(mice_diagnostics) + 1L]] <- data.table(
        cohort = cohort_name,
        status = mi$status,
        message = mi$message,
        warnings = mi$warnings,
        imputations_requested = config$mice_m,
        imputations_created = length(mi$completed),
        protected_field_violations = mi$protection_violations,
        logged_events = nrow(mi$log),
        cache_status = mi$cache_status %||% "not_available",
        cache_file = mi$cache_file %||% "",
        baseline_sha256 = mi$baseline_sha256 %||% mi$cache_signature$baseline_sha256 %||% "",
        imputation_frame_sha256 = mi$imputation_frame_sha256 %||% mi$cache_signature$imputation_frame_sha256 %||% "",
        imputation_plan_sha256 = mi$imputation_plan_sha256 %||% mi$cache_signature$imputation_plan_sha256 %||% "",
        completed_sha256 = mi$completed_sha256 %||% ""
      )
      if (nrow(mi$log)) atomic_fwrite(mi$log, file.path(qc_root, paste0(tolower(cohort_name), "_mice_logged_events.csv")))
    }
  }
  mice_diagnostics <- bind_or_empty(
    mice_diagnostics,
    c("cohort", "status", "message", "warnings", "imputations_requested", "imputations_created", "protected_field_violations", "logged_events", "cache_status", "cache_file", "baseline_sha256", "imputation_frame_sha256", "imputation_plan_sha256", "completed_sha256")
  )
  atomic_fwrite(mice_diagnostics, file.path(qc_root, "mice_diagnostics.csv"))

  run_rows <- list()
  estimate_rows <- list()
  term_test_rows <- list()
  ph_rows <- list()
  ph_summary_rows <- list()
  qc_rows <- list()
  individual_rows <- list()
  model_objects <- list()

  for (i in seq_len(nrow(selected))) {
    spec <- selected[i]
    message(sprintf(
      "[%d/%d] %s %s %s %s %s",
      i, nrow(selected), spec$cohort, spec$transition_id,
      spec$analysis_id, spec$model_level, spec$basis
    ))
    result <- load_spec_checkpoint(spec)
    if (is.null(result)) {
      result <- run_spec(
        intervals_by_cohort[[spec$cohort]], spec,
        mice_by_cohort[[spec$cohort]]
      )
      result$run[, checkpoint_status := "computed_and_checkpointed"]
      save_spec_checkpoint(spec, result)
    } else {
      result$run[, checkpoint_status := "reused_valid_checkpoint"]
    }
    run_rows[[length(run_rows) + 1L]] <- result$run
    if (nrow(result$estimates)) estimate_rows[[length(estimate_rows) + 1L]] <- result$estimates
    if (nrow(result$term_tests)) term_test_rows[[length(term_test_rows) + 1L]] <- result$term_tests
    if (nrow(result$ph)) ph_rows[[length(ph_rows) + 1L]] <- result$ph
    if (nrow(result$ph_summary)) ph_summary_rows[[length(ph_summary_rows) + 1L]] <- result$ph_summary
    if (nrow(result$qc)) qc_rows[[length(qc_rows) + 1L]] <- result$qc
    if (nrow(result$individual_estimates)) {
      result$individual_estimates <- cbind(spec_metadata(spec)[rep(1L, nrow(result$individual_estimates))], result$individual_estimates)
      individual_rows[[length(individual_rows) + 1L]] <- result$individual_estimates
    }
    if (length(result$objects)) model_objects[[spec$spec_id]] <- result$objects
    if (i %% 10L == 0L || i == nrow(selected)) {
      write_run_state(
        "running", nrow(selected), i,
        paste0("Completed or restored ", i, " of ", nrow(selected), " specifications"),
        catalog_fingerprint
      )
    }
  }

  runs <- bind_or_empty(run_rows, c("spec_id", "status", "message"))
  estimates <- bind_or_empty(estimate_rows, c("spec_id", "term", "estimate", "robust_se", "p_value", "hazard_ratio", "conf_low", "conf_high"))
  term_tests <- bind_or_empty(term_test_rows, c("spec_id", "term_block", "terms", "statistic", "df", "p_value", "m_used", "method"))
  if (nrow(runs) != nrow(selected) || anyDuplicated(runs$spec_id) ||
      !setequal(runs$spec_id, selected$spec_id)) {
    stop("Checkpoint aggregation is incomplete or contains duplicate/mismatched specification IDs")
  }
  if (nrow(estimates) && anyDuplicated(estimates[, .(spec_id, term)])) {
    stop("Checkpoint aggregation produced duplicate estimate keys (spec_id, term)")
  }
  if (nrow(term_tests) && anyDuplicated(term_tests[, .(spec_id, term_block)])) {
    stop("Checkpoint aggregation produced duplicate term-test keys (spec_id, term_block)")
  }
  fdr_result <- if (nrow(estimates)) {
    add_bh_fdr(estimates, runs, term_tests)
  } else {
    list(estimates = estimates, audit = empty_table(c(
      "cohort", "confirmatory_family_id", "expected_edges", "expected_family_size",
      "run_edges_observed", "estimate_edges_observed", "joint_test_edges_observed",
      "clean_models", "family_complete", "joint_overall_gate_method", "joint_overall_p_value",
      "alpha", "joint_overall_gate_pass", "bh_applied", "inference_reason"
    )))
  }
  estimates <- fdr_result$estimates
  fdr_audit <- fdr_result$audit
  if (nrow(runs)) {
    runs[, `:=`(
      primary_inference_eligible_after_family_gate = 0L,
      primary_family_gate_reason = fifelse(
        claim_bearing_primary == 1L,
        "withheld_no_complete_four_edge_family_gate",
        "not_claim_bearing_primary"
      )
    )]
    if (nrow(fdr_audit)) {
      for (i in seq_len(nrow(fdr_audit))) {
        gate <- fdr_audit[i]
        runs[
          cohort == gate$cohort & claim_bearing_primary == 1L,
          `:=`(
            primary_inference_eligible_after_family_gate = as.integer(gate$joint_overall_gate_pass),
            primary_family_gate_reason = gate$inference_reason
          )
        ]
      }
    }
  }
  ph <- bind_or_empty(ph_rows, c("spec_id", "imputation", "ph_term", "chisq", "df", "p_value", "message"))
  ph_summary <- bind_or_empty(ph_summary_rows, c("spec_id", "ph_term", "imputations_with_test", "p_min", "p_median", "proportion_p_lt_0_05", "interpretation"))
  fit_qc <- bind_or_empty(qc_rows, c("spec_id", "imputation", "status", "reason", "n_events", "n_parameters", "events_per_parameter"))
  individual_estimates <- bind_or_empty(individual_rows, c("spec_id", "imputation", "term", "estimate", "within_variance"))
  failures <- runs[!status %chin% c("success", "success_with_warning")]
  primary_inference_withheld <- runs[
    claim_bearing_primary == 1L & primary_inference_eligible_after_family_gate != 1L
  ]
  cox_diagnostics <- make_canonical_cox_diagnostics(edge_registry, runs, ph)

  event_parameter_qc <- if (nrow(fit_qc)) fit_qc[, .(
    imputations_or_fits = .N,
    fits_estimable = sum(status %chin% c("success", "success_with_warning")),
    events_min = finite_min(n_events),
    events_max = finite_max(n_events),
    parameters_max = finite_max(n_parameters),
    epp_min = finite_min(events_per_parameter),
    statuses = collapse_messages(status),
    reasons = collapse_messages(reason)
  ), by = .(spec_id, cohort, transition_id, transition, analysis_id, model_level, basis, weighted)] else empty_table(
    c("spec_id", "cohort", "transition_id", "transition", "analysis_id", "model_level", "basis", "weighted", "imputations_or_fits", "fits_estimable", "events_min", "events_max", "parameters_max", "epp_min", "statuses", "reasons")
  )

  atomic_fwrite(runs, file.path(result_root, "transition_cox_model_runs.csv"))
  atomic_fwrite(estimates, file.path(result_root, "transition_cox_estimates.csv"))
  atomic_fwrite(term_tests, file.path(result_root, "transition_cox_term_tests.csv"))
  atomic_fwrite(fdr_audit, file.path(result_root, "transition_cox_primary_fdr_gate.csv"))
  atomic_fwrite(ph, file.path(result_root, "transition_cox_ph_tests.csv"))
  atomic_fwrite(ph_summary, file.path(result_root, "transition_cox_ph_summary.csv"))
  atomic_fwrite(cox_diagnostics, file.path(result_root, "transition_cox_diagnostics.csv"))
  atomic_fwrite(failures, file.path(qc_root, "transition_cox_model_failures.csv"))
  atomic_fwrite(primary_inference_withheld, file.path(qc_root, "transition_cox_primary_inference_withheld.csv"))
  atomic_fwrite(fit_qc, file.path(qc_root, "transition_cox_fit_qc.csv"))
  atomic_fwrite(event_parameter_qc, file.path(qc_root, "transition_cox_event_parameter_qc.csv"))
  if (config$save_imputation_estimates) {
    atomic_fwrite(individual_estimates, file.path(result_root, "transition_cox_imputation_estimates.csv"))
  }
  result_bundle <- list(
    metadata = make_run_metadata(catalog_fingerprint, protocol_catalog_fingerprint),
    run_fingerprint = run_fingerprint,
    catalog_fingerprint = catalog_fingerprint,
    protocol_catalog_fingerprint = protocol_catalog_fingerprint,
    protocol_catalog_fingerprint = protocol_catalog_fingerprint,
    input_manifest = input_manifest,
    configuration = config,
    specification_qc = specification_qc,
    specification_catalog = catalog,
    selected_specifications = selected,
    input_qc = input_qc,
    edge_registry = edge_registry,
    mice_diagnostics = mice_diagnostics,
    model_runs = runs,
    estimates = estimates,
    cox_diagnostics = cox_diagnostics,
    term_tests = term_tests,
    primary_fdr_gate = fdr_audit,
    primary_inference_withheld = primary_inference_withheld,
    proportional_hazards_tests = ph,
    proportional_hazards_summary = ph_summary,
    fit_qc = fit_qc,
    event_parameter_qc = event_parameter_qc,
    model_objects = model_objects
  )
  atomic_saveRDS(result_bundle, file.path(result_root, "transition_cox_results.rds"))

  key_outputs <- c(
    "transition_cox_model_runs.csv", "transition_cox_estimates.csv",
    "transition_cox_term_tests.csv", "transition_cox_primary_fdr_gate.csv",
    "transition_cox_diagnostics.csv", "transition_cox_results.rds"
  )
  output_manifest <- rbindlist(lapply(key_outputs, function(name) {
    path <- file.path(result_root, name)
    data.table(file = name, bytes = as.numeric(file.info(path)$size), sha256 = sha256_file(path))
  }))
  completion_manifest <- data.table(
    status = "completed",
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    run_fingerprint = run_fingerprint,
    catalog_fingerprint = catalog_fingerprint,
    expected_specifications = nrow(selected),
    checkpoint_specifications = length(unique(runs$spec_id)),
    run_records = nrow(runs),
    estimate_records = nrow(estimates),
    primary_families_complete = if (nrow(fdr_audit)) sum(fdr_audit$family_complete == 1L) else 0L,
    primary_families_bh_applied = if (nrow(fdr_audit)) sum(fdr_audit$bh_applied == 1L) else 0L,
    key_output_manifest_sha256 = sha256_object(output_manifest)
  )
  atomic_fwrite(output_manifest, file.path(result_root, "OUTPUT_FILE_MANIFEST.csv"))
  atomic_fwrite(completion_manifest, file.path(result_root, "COMPLETION_MANIFEST.csv"))
  write_run_state(
    "completed", nrow(selected), length(unique(runs$spec_id)),
    "All selected specifications reached a terminal checkpoint and final outputs were fingerprinted",
    catalog_fingerprint
  )
  run_completed <- TRUE

  canonical_status <- tryCatch(
    publish_canonical_outputs(catalog, runs, estimates, cox_diagnostics, fdr_audit, completion_manifest),
    error = function(e) data.table(
      requested = as.integer(config$publish_canonical), eligible = NA_integer_, published = 0L,
      run_fingerprint = run_fingerprint, reason = paste0("canonical_publish_failed: ", conditionMessage(e)),
      canonical_bundle = NA_character_
    )
  )
  atomic_fwrite(canonical_status, file.path(result_root, "canonical_publish_status.csv"))

  message("Transition-specific Cox sensitivity analysis completed.")
  message("Results: ", result_root)
  message("QC: ", qc_root)
  invisible(result_bundle)
}

if (sys.nframe() == 0L) main()
