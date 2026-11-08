#!/usr/bin/env Rscript

# Continuous-time panel-observed Markov models for the HRS-ELSA insomnia/CVD
# study.  This script deliberately keeps model fitting separate by cohort and
# mortality-coverage window.  It is restartable: every fitted model and every
# result table is written independently before the next specification starts.

suppressPackageStartupMessages({
  library(data.table)
  library(GlobalAgingR)
  library(msm)
})

options(stringsAsFactors = FALSE, width = 160)

script_version <- "2026-09-08-v3"
stage_root <- Sys.getenv("GA_STAGE_ROOT", unset = "private/stage")
multistate_root <- Sys.getenv(
  "GA_MULTISTATE_OUTPUT_ROOT",
  unset = file.path(stage_root, "multistate_output")
)
input_root <- file.path(multistate_root, "data_multistate")
result_root <- Sys.getenv(
  "GA_MARKOV_RESULT_ROOT",
  unset = file.path(multistate_root, "results", "markov")
)
model_root <- Sys.getenv(
  "GA_MARKOV_MODEL_ROOT",
  unset = file.path(multistate_root, "models", "markov")
)
log_root <- Sys.getenv(
  "GA_MARKOV_LOG_ROOT",
  unset = file.path(multistate_root, "logs", "markov")
)
by_model_root <- file.path(result_root, "by_model")
qc_root <- file.path(multistate_root, "qc", "markov")
mi_root <- file.path(multistate_root, "data_multistate", "mi_baseline")

for (path in c(result_root, model_root, log_root, by_model_root, qc_root, mi_root)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

env_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = if (default) "1" else "0")))
  value %in% c("1", "true", "yes", "y", "on")
}

env_number <- function(name, default = 0) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value)) default else value
}

run_mode <- tolower(Sys.getenv("GA_MARKOV_RUN_MODE", unset = "full"))
if (!run_mode %in% c("full", "prototype", "spec")) {
  stop("GA_MARKOV_RUN_MODE must be one of: full, prototype, spec")
}
model_regex <- Sys.getenv("GA_MARKOV_MODEL_REGEX", unset = ".")
analysis_regex <- Sys.getenv("GA_MARKOV_ANALYSIS_REGEX", unset = ".")
force_refit <- env_flag("GA_MARKOV_FORCE_REFIT", FALSE)
run_globalagingr_baseline <- env_flag("GA_RUN_GLOBALAGINGR_BASELINE", TRUE)
max_seconds <- env_number("GA_MARKOV_MAX_SECONDS_PER_MODEL", 0)
max_panel_rows <- as.integer(env_number("GA_MARKOV_MAX_PANEL_ROWS", 0))
max_models <- as.integer(env_number("GA_MARKOV_MAX_MODELS", 0))
max_iterations <- as.integer(env_number("GA_MARKOV_MAX_ITER", 1000))
elsa_death_cutoff <- env_number("GA_ELSA_DEATH_CUTOFF", 2012.999999)
mi_mode <- tolower(Sys.getenv("GA_MARKOV_MI_MODE", unset = "run"))
if (!mi_mode %in% c("prepare", "run", "off")) {
  stop("GA_MARKOV_MI_MODE must be one of: prepare, run, off")
}
mi_m <- as.integer(env_number("GA_MARKOV_MI_M", 20))
mi_maxit <- as.integer(env_number("GA_MARKOV_MI_MAXIT", 10))
mi_seed <- as.integer(env_number("GA_MARKOV_MI_SEED", 20260908))
mi_start <- as.integer(env_number("GA_MARKOV_MI_START", 1))
mi_end <- as.integer(env_number("GA_MARKOV_MI_END", mi_m))
mi_pool <- env_flag("GA_MARKOV_MI_POOL", TRUE)
mi_start <- max(1L, min(mi_start, mi_m))
mi_end <- max(mi_start, min(mi_end, mi_m))
if (mi_mode == "run" && mi_pool && (mi_start != 1L || mi_end != mi_m)) {
  stop(
    "MI pooling requires the complete imputation range 1:", mi_m,
    ". Run shards with GA_MARKOV_MI_POOL=0, then make one 1:", mi_m,
    " cache-reading pooling pass."
  )
}
nuisance_constraint_mode <- tolower(Sys.getenv(
  "GA_MARKOV_NUISANCE_CONSTRAINT",
  unset = "common"
))
if (!nuisance_constraint_mode %in% c("common", "none")) {
  stop("GA_MARKOV_NUISANCE_CONSTRAINT must be 'common' or 'none'")
}

failure_columns <- c(
  "timestamp", "analysis_id", "model_id", "engine", "stage", "status",
  "elapsed_seconds", "message", "warnings", "log_file"
)
failure_rows <- list()
current_run_table_files <- new.env(parent = emptyenv())
current_mi_prefixes <- character()

register_run_table <- function(suffix, path) {
  existing <- if (exists(suffix, envir = current_run_table_files, inherits = FALSE)) {
    get(suffix, envir = current_run_table_files, inherits = FALSE)
  } else character()
  assign(suffix, unique(c(existing, normalizePath(path, winslash = "/", mustWork = FALSE))),
    envir = current_run_table_files
  )
  invisible(path)
}

register_mi_prefix <- function(prefix) {
  current_mi_prefixes <<- unique(c(
    current_mi_prefixes,
    normalizePath(prefix, winslash = "/", mustWork = FALSE)
  ))
  invisible(prefix)
}

fingerprint_table <- function(x, columns) {
  columns <- unique(intersect(columns, names(x)))
  if (!length(columns)) stop("Cannot fingerprint a table without named columns")
  value <- copy(as.data.table(x)[, ..columns])
  digest::digest(value, algo = "xxhash64", serialize = TRUE)
}

record_failure <- function(
    analysis_id, model_id, engine, stage, status, message,
    warnings = character(), elapsed_seconds = NA_real_, log_file = NA_character_) {
  failure_rows[[length(failure_rows) + 1L]] <<- data.table(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    analysis_id = analysis_id,
    model_id = model_id,
    engine = engine,
    stage = stage,
    status = status,
    elapsed_seconds = as.numeric(elapsed_seconds),
    message = as.character(message),
    warnings = paste(unique(warnings), collapse = " | "),
    log_file = as.character(log_file)
  )
}

flush_failures <- function() {
  current <- if (length(failure_rows)) {
    rbindlist(failure_rows, use.names = TRUE, fill = TRUE)
  } else {
    as.data.table(setNames(replicate(length(failure_columns), character(), simplify = FALSE), failure_columns))
  }
  path <- file.path(log_root, "markov_failures.csv")
  if (file.exists(path)) {
    old <- tryCatch(fread(path), error = function(e) NULL)
    if (!is.null(old) && nrow(old)) current <- rbindlist(list(old, current), fill = TRUE)
  }
  if (nrow(current)) current <- unique(current)
  fwrite(current, path, na = "")
  invisible(current)
}

safe_run <- function(fun, log_file, analysis_id, model_id, engine, stage) {
  started <- Sys.time()
  warnings <- character()
  error_message <- NULL
  value <- NULL
  if (max_seconds > 0) {
    setTimeLimit(cpu = Inf, elapsed = max_seconds, transient = TRUE)
  }
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
  output <- capture.output({
    value <- tryCatch(
      withCallingHandlers(
        fun(),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    )
  }, type = "output")
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  writeLines(
    c(
      paste("analysis_id:", analysis_id),
      paste("model_id:", model_id),
      paste("engine:", engine),
      paste("stage:", stage),
      paste("started:", format(started, "%Y-%m-%d %H:%M:%S %z")),
      paste("elapsed_seconds:", sprintf("%.3f", elapsed)),
      if (length(warnings)) c("warnings:", paste0("- ", unique(warnings))) else "warnings: none",
      if (!is.null(error_message)) paste("error:", error_message) else "error: none",
      "captured_output:",
      output
    ),
    con = log_file,
    useBytes = TRUE
  )
  if (!is.null(error_message)) {
    record_failure(
      analysis_id, model_id, engine, stage, "FAILED", error_message,
      warnings, elapsed, log_file
    )
  }
  list(
    value = value,
    warnings = unique(warnings),
    error = error_message,
    elapsed_seconds = elapsed,
    log_file = log_file
  )
}

locf_down <- function(x) {
  if (!length(x)) return(x)
  last <- NA
  for (i in seq_along(x)) {
    if (!is.na(x[i])) last <- x[i] else if (!is.na(last)) x[i] <- last
  }
  x
}

state_name <- function(index) paste0("S", as.integer(index) - 1L)

allowed_edges <- function(n_states) {
  if (n_states == 5L) {
    edges <- data.table(
      from_index = c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L),
      to_index = c(2L, 3L, 5L, 4L, 5L, 4L, 5L, 5L)
    )
  } else if (n_states == 4L) {
    edges <- data.table(
      from_index = c(1L, 1L, 2L, 3L),
      to_index = c(2L, 3L, 4L, 4L)
    )
  } else {
    stop("Only the protocol-defined four- or five-state models are supported")
  }
  setorder(edges, from_index, to_index)
  edges[, `:=`(
    from_state = from_index - 1L,
    to_state = to_index - 1L,
    transition = paste0(state_name(from_index), "->", state_name(to_index)),
    death_transition = as.integer(n_states == 5L & to_index == 5L)
  )]
  edges[]
}

make_q_initial <- function(empirical_q, edges, n_states, minimum = 1e-4) {
  q <- matrix(
    0, nrow = n_states, ncol = n_states,
    dimnames = list(state_name(seq_len(n_states)), state_name(seq_len(n_states)))
  )
  for (i in seq_len(nrow(edges))) {
    from <- edges$from_index[i]
    to <- edges$to_index[i]
    value <- suppressWarnings(as.numeric(empirical_q[from, to]))
    if (!is.finite(value) || value <= 0) value <- minimum
    q[from, to] <- max(value, minimum)
  }
  diag(q) <- 0
  if (n_states == 5L) q[5L, ] <- 0
  q
}

matrix_long <- function(mat, analysis_id, matrix_type, allowed = NULL) {
  idx <- expand.grid(
    from_index = seq_len(nrow(mat)),
    to_index = seq_len(ncol(mat)),
    KEEP.OUT.ATTRS = FALSE
  )
  out <- as.data.table(idx)
  out[, value := as.numeric(mat[cbind(from_index, to_index)])]
  out[, `:=`(
    analysis_id = analysis_id,
    matrix_type = matrix_type,
    from_state = from_index - 1L,
    to_state = to_index - 1L,
    transition = paste0(state_name(from_index), "->", state_name(to_index))
  )]
  if (is.null(allowed)) {
    out[, allowed_intensity := NA_integer_]
  } else {
    keys <- paste(allowed$from_index, allowed$to_index)
    out[, allowed_intensity := as.integer(paste(from_index, to_index) %chin% keys)]
  }
  setcolorder(out, c(
    "analysis_id", "matrix_type", "from_index", "to_index", "from_state",
    "to_state", "transition", "allowed_intensity", "value"
  ))
  out[]
}

analysis_configs <- list(
  list(
    analysis_id = "hrs_5state_full",
    cohort = "HRS",
    input_file = file.path(input_root, "hrs_msm_panel.rds"),
    n_states = 5L,
    mortality_scope = "HRS follow-up through 2020; death time from last interview plus RAND months-to-death when available",
    death_exact_state = 5L,
    data_rule = "all"
  ),
  list(
    analysis_id = "elsa_4state_w4_w10",
    cohort = "ELSA",
    input_file = file.path(input_root, "elsa_msm_panel.rds"),
    n_states = 4L,
    mortality_scope = "W4-W10 living-state external validation; death excluded",
    death_exact_state = NULL,
    data_rule = "nonfatal"
  ),
  list(
    analysis_id = "elsa_5state_through_2012",
    cohort = "ELSA",
    input_file = file.path(input_root, "elsa_msm_panel.rds"),
    n_states = 5L,
    mortality_scope = paste0(
      "Mortality-linkage-limited validation; records through calendar year ",
      floor(elsa_death_cutoff)
    ),
    death_exact_state = NULL,
    data_rule = "through_2012"
  )
)

prepare_analysis_data <- function(config) {
  if (!file.exists(config$input_file)) {
    stop("Missing input: ", config$input_file, ". Run R/05_prepare_multistate_data.R first.")
  }
  x <- as.data.table(readRDS(config$input_file))
  required <- c(
    "id", "cohort", "state", "event_time", "exposure_entry_time", "event_type",
    "wave", "insomnia_score_primary"
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Input is missing columns: ", paste(missing, collapse = ", "))

  # Recreate these variables explicitly rather than trusting serialized copies.
  x[, state_msm := as.integer(state) + 1L]
  x[, time_msm := as.numeric(event_time) - as.numeric(exposure_entry_time)]
  setorder(x, id, time_msm, state_msm, na.last = TRUE)
  x[, insomnia_score_primary := locf_down(as.numeric(insomnia_score_primary)), by = id]

  if (config$data_rule == "nonfatal") {
    x <- x[event_type == "visit" & wave %in% 4:10 & state_msm %in% 1:4]
  } else if (config$data_rule == "through_2012") {
    x <- x[event_time <= elsa_death_cutoff & state_msm %in% 1:5]
  } else {
    x <- x[state_msm %in% 1:5]
  }
  x <- x[!is.na(id) & is.finite(time_msm) & time_msm >= 0]
  setorder(x, id, time_msm, state_msm)
  eligible_ids <- x[, .N, by = id][N >= 2L, id]
  x <- x[id %in% eligible_ids]

  # Stable scales and explicit reference categories aid optimisation and make
  # cross-cohort model specifications auditable.
  if ("race_bl" %in% names(x)) {
    if (config$cohort == "HRS") {
      x[, race_bl := factor(
        race_bl,
        levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other")
      )]
    } else {
      x[, race_bl := factor(race_bl, levels = c("White", "Non-White"))]
    }
  }
  if ("education_level_bl" %in% names(x)) {
    x[, education_level_bl := factor(education_level_bl, levels = 1:3)]
  }
  if ("wealth_quintile_bl" %in% names(x)) {
    x[, wealth_quintile_bl := factor(wealth_quintile_bl, levels = 1:5)]
  }
  x[, age_bl_5 := (as.numeric(age_bl) - 65) / 5]
  x[, bmi_bl_5 := (as.numeric(bmi_bl) - 25) / 5]
  x[, depression_score_bl_1 := as.numeric(depression_score_bl)]
  x[, insomnia_burden_factor := factor(
    as.integer(insomnia_burden_code),
    levels = 0:2,
    labels = c("0 symptoms", "1 symptom", ">=2 symptoms")
  )]
  x[]
}

analysis_qc <- function(x, config, edges) {
  transitions <- copy(x[, .(id, time_msm, state_msm)])
  setorder(transitions, id, time_msm)
  transitions[, `:=`(
    next_state = shift(state_msm, type = "lead"),
    next_time = shift(time_msm, type = "lead")
  ), by = id]
  pairs <- transitions[!is.na(next_state)]
  allowed_keys <- paste(edges$from_index, edges$to_index)
  pair_keys <- paste(pairs$state_msm, pairs$next_state)
  multistep_key <- if (config$n_states == 5L) "1 4" else "1 4"
  admissible <- pairs$state_msm == pairs$next_state |
    pair_keys %chin% allowed_keys | pair_keys == multistep_key
  death_not_final <- if (config$n_states == 5L) {
    x[, any(state_msm == 5L & seq_len(.N) != .N), by = id][V1 == TRUE, .N]
  } else 0L
  checks <- data.table(
    analysis_id = config$analysis_id,
    check = c(
      "state_msm explicitly equals state + 1",
      "time_msm explicitly uses actual elapsed years",
      "strictly increasing times within participant",
      "at least two panel records per participant",
      "states are inside the analysis state space",
      "death rows retain LOCF interval-start insomnia exposure",
      "death is final and absorbing in observed records",
      "observed pairs obey the protocol graph or S0-to-S3 multistep jump",
      "S0-to-S3 direct intensity is absent from the main Q"
    ),
    violations = c(
      x[state_msm != as.integer(state) + 1L, .N],
      x[abs(time_msm - (as.numeric(event_time) - as.numeric(exposure_entry_time))) > 1e-8, .N],
      x[, any(diff(time_msm) <= 0), by = id][V1 == TRUE, .N],
      x[, .N, by = id][N < 2L, .N],
      x[!state_msm %in% seq_len(config$n_states), .N],
      if (config$n_states == 5L) x[state_msm == 5L & is.na(insomnia_score_primary), .N] else 0L,
      death_not_final,
      sum(!admissible),
      edges[from_index == 1L & to_index == 4L, .N]
    )
  )
  checks[, critical := as.integer(check != "S0-to-S3 direct intensity is absent from the main Q")]
  checks[, status := fifelse(violations == 0L, "PASS", "FAIL")]
  checks[]
}

primary_m4_adjustment <- c(
  "age_bl_5", "female_bl", "race_bl", "education_level_bl",
  "wealth_quintile_bl", "married_bl", "current_smoker_bl",
  "current_drinker_bl", "inactive_bl", "bmi_bl_5", "diabetes_bl",
  "lung_disease_bl", "frequent_pain_bl"
)

model_definitions <- list(
  M1 = c("insomnia_score_primary", "age_bl_5", "female_bl"),
  M2 = c(
    "insomnia_score_primary", "age_bl_5", "female_bl", "race_bl",
    "education_level_bl", "wealth_quintile_bl", "married_bl"
  ),
  M3 = c(
    "insomnia_score_primary", "age_bl_5", "female_bl", "race_bl",
    "education_level_bl", "wealth_quintile_bl", "married_bl",
    "current_smoker_bl", "current_drinker_bl", "inactive_bl", "bmi_bl_5"
  ),
  M4 = c(
    "insomnia_score_primary", "age_bl_5", "female_bl", "race_bl",
    "education_level_bl", "wealth_quintile_bl", "married_bl",
    "current_smoker_bl", "current_drinker_bl", "inactive_bl", "bmi_bl_5",
    "diabetes_bl", "lung_disease_bl", "frequent_pain_bl"
  ),
  M4_CESD_sensitivity = c(
    "insomnia_score_primary", "age_bl_5", "female_bl", "race_bl",
    "education_level_bl", "wealth_quintile_bl", "married_bl",
    "current_smoker_bl", "current_drinker_bl", "inactive_bl", "bmi_bl_5",
    "diabetes_bl", "lung_disease_bl", "frequent_pain_bl", "depression_score_bl_1"
  ),
  M4_burden_categorical = c("insomnia_burden_factor", primary_m4_adjustment),
  M4_onset_symptom = c("insomnia_onset_primary", primary_m4_adjustment),
  M4_maintenance_symptom = c("insomnia_maintenance_primary", primary_m4_adjustment),
  M4_nonrestorative_symptom = c("insomnia_nonrestorative_primary", primary_m4_adjustment)
)

model_exposures <- c(
  M1 = "insomnia_score_primary",
  M2 = "insomnia_score_primary",
  M3 = "insomnia_score_primary",
  M4 = "insomnia_score_primary",
  M4_CESD_sensitivity = "insomnia_score_primary",
  M4_burden_categorical = "insomnia_burden_factor",
  M4_onset_symptom = "insomnia_onset_primary",
  M4_maintenance_symptom = "insomnia_maintenance_primary",
  M4_nonrestorative_symptom = "insomnia_nonrestorative_primary"
)

model_labels <- c(
  M1 = "M1: age and sex",
  M2 = "M2: M1 plus race/ethnicity, education, wealth quintile and marital status",
  M3 = "M3: M2 plus smoking, drinking, physical inactivity and BMI",
  M4 = "M4 primary complete-case prototype: M3 plus diabetes, lung disease and pain; CES-D omitted",
  M4_CESD_sensitivity = "M4 depression-adjusted sensitivity: M4 plus harmonized CES-D summary",
  M4_burden_categorical = "M4 categorical dose-response: 0, 1, or at least 2 symptoms",
  M4_onset_symptom = "M4 symptom-domain model: sleep-onset difficulty",
  M4_maintenance_symptom = "M4 symptom-domain model: sleep-maintenance difficulty",
  M4_nonrestorative_symptom = "M4 symptom-domain model: nonrestorative sleep"
)

make_specifications <- function() {
  rows <- list()
  for (config in analysis_configs) {
    rows[[length(rows) + 1L]] <- data.table(
      spec_key = paste(config$analysis_id, "GA_M0", sep = "__"),
      analysis_id = config$analysis_id,
      cohort = config$cohort,
      mortality_scope = config$mortality_scope,
      n_states = config$n_states,
      deathexact_state = NA_integer_,
      observation_scheme = "Panel-observed package smoke model; GlobalAgingR wrapper does not expose deathexact",
      engine = "GlobalAgingR::msm_transition_model",
      model_id = "GA_M0",
      model_label = "GlobalAgingR package-path baseline: time-varying insomnia score",
      exposure = "insomnia_score_primary (0-3, per one symptom)",
      covariates = "insomnia_score_primary",
      insomnia_effect = "Transition-specific under msm default parameterization",
      nuisance_constraint = "Not applicable",
      missing_data = "Complete records; exposure/state/time never imputed",
      status = "PLANNED"
    )
    for (model_id in names(model_definitions)) {
      exposure_var <- unname(model_exposures[model_id])
      rows[[length(rows) + 1L]] <- data.table(
        spec_key = paste(config$analysis_id, model_id, sep = "__"),
        analysis_id = config$analysis_id,
        cohort = config$cohort,
        mortality_scope = config$mortality_scope,
        n_states = config$n_states,
        deathexact_state = if (is.null(config$death_exact_state)) NA_integer_ else config$death_exact_state,
        observation_scheme = if (is.null(config$death_exact_state)) {
          "Panel-observed continuous-time Markov model"
        } else {
          paste0("Panel-observed model with exact entry into state ", config$death_exact_state)
        },
        engine = "msm::msm",
        model_id = model_id,
        model_label = unname(model_labels[model_id]),
        exposure = exposure_var,
        covariates = paste(model_definitions[[model_id]], collapse = " + "),
        insomnia_effect = paste0(
          "Transition-specific across all allowed intensities: constraint indices 1:",
          nrow(allowed_edges(config$n_states))
        ),
        nuisance_constraint = if (nuisance_constraint_mode == "common") {
          "Each non-insomnia model-matrix coefficient constrained equal across allowed transitions"
        } else {
          "All covariate effects transition-specific"
        },
        missing_data = "Participant-level complete-case prototype",
        status = "PLANNED"
      )
    }
    rows[[length(rows) + 1L]] <- data.table(
      spec_key = paste(config$analysis_id, "M4_MI20", sep = "__"),
      analysis_id = config$analysis_id,
      cohort = config$cohort,
      mortality_scope = config$mortality_scope,
      n_states = config$n_states,
      deathexact_state = if (is.null(config$death_exact_state)) NA_integer_ else config$death_exact_state,
      observation_scheme = if (is.null(config$death_exact_state)) {
        "Panel-observed continuous-time Markov model"
      } else {
        paste0("Panel-observed model with exact entry into state ", config$death_exact_state)
      },
      engine = "mice::mice + msm::msm + Rubin pooling",
      model_id = "M4_MI20",
      model_label = paste0(
        "Primary M4 continuous-exposure model after ", mi_m,
        " baseline-covariate imputations"
      ),
      exposure = "insomnia_score_primary (0-3, never imputed)",
      covariates = paste(c("insomnia_score_primary", primary_m4_adjustment), collapse = " + "),
      insomnia_effect = paste0(
        "Transition-specific across all allowed intensities: constraint indices 1:",
        nrow(allowed_edges(config$n_states))
      ),
      nuisance_constraint = "Each imputed nuisance coefficient constrained equal across allowed transitions",
      missing_data = paste0(
        "MICE m=", mi_m,
        "; baseline covariates only; exposure/state/time/cohort/weight protected"
      ),
      status = if (mi_mode == "run") "PLANNED" else if (mi_mode == "prepare") {
        "NOT_RUN_COMPUTE_LIMIT_MI_DATA_PREPARATION_PLANNED"
      } else {
        "DISABLED"
      }
    )
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

specifications <- make_specifications()
fwrite(specifications, file.path(result_root, "markov_model_specifications.csv"), na = "")

prepare_model_sample <- function(x, requested_covariates, exposure) {
  missing_columns <- setdiff(requested_covariates, names(x))
  if (length(missing_columns)) {
    stop(
      "The prespecified model cannot be silently reduced; missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  available <- intersect(requested_covariates, names(x))
  if (!length(available) || !exposure %in% available) {
    stop("Requested insomnia exposure is unavailable: ", exposure)
  }
  required <- unique(c("id", "time_msm", "state_msm", available))
  d <- copy(x[complete.cases(x[, ..required])])
  eligible_ids <- d[, .N, by = id][N >= 2L, id]
  d <- d[id %in% eligible_ids]
  setorder(d, id, time_msm, state_msm)
  if ("race_bl" %in% names(d)) d[, race_bl := droplevels(race_bl)]
  if ("education_level_bl" %in% names(d)) d[, education_level_bl := droplevels(education_level_bl)]
  if ("wealth_quintile_bl" %in% names(d)) d[, wealth_quintile_bl := droplevels(wealth_quintile_bl)]

  constant <- available[vapply(available, function(v) {
    z <- d[[v]]
    length(unique(z[!is.na(z)])) < 2L
  }, logical(1))]
  if (exposure %in% constant) {
    stop("Requested insomnia exposure has no variation in the complete-case sample: ", exposure)
  }
  used <- setdiff(available, constant)
  list(
    data = d,
    requested = requested_covariates,
    used = used,
    missing_columns = missing_columns,
    constant_columns = constant
  )
}

prepare_mi_baseline <- function(x, config) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required to prepare the 20-imputation primary M4 analysis")
  }
  target_variables <- primary_m4_adjustment
  missing_targets <- setdiff(target_variables, names(x))
  if (length(missing_targets)) {
    stop("MI target columns are unavailable: ", paste(missing_targets, collapse = ", "))
  }
  setorder(x, id, time_msm)
  baseline <- x[, .SD[1L], by = id, .SDcols = unique(c(
    target_variables, "insomnia_score_primary", "state_msm", "weight_bl"
  ))]
  auxiliary <- x[, .(
    followup_years = max(time_msm),
    n_panel_records = .N,
    n_state_changes = sum(diff(state_msm) != 0),
    death_observed = as.integer(any(state_msm == 5L)),
    last_state_msm = state_msm[.N]
  ), by = id]
  frame <- merge(baseline, auxiliary, by = "id", all.x = TRUE, sort = FALSE)
  setorder(frame, id)
  frame[, state_msm := factor(state_msm, levels = seq_len(config$n_states))]
  frame[, last_state_msm := factor(last_state_msm, levels = seq_len(config$n_states))]
  binary_targets <- intersect(c(
    "female_bl", "married_bl", "current_smoker_bl", "current_drinker_bl",
    "inactive_bl", "diabetes_bl", "lung_disease_bl", "frequent_pain_bl"
  ), names(frame))
  for (name in binary_targets) {
    set(frame, j = name, value = factor(frame[[name]], levels = 0:1))
  }
  if ("education_level_bl" %in% names(frame)) {
    frame[, education_level_bl := ordered(education_level_bl, levels = 1:3)]
  }
  if ("wealth_quintile_bl" %in% names(frame)) {
    frame[, wealth_quintile_bl := ordered(wealth_quintile_bl, levels = 1:5)]
  }
  if ("race_bl" %in% names(frame)) frame[, race_bl := droplevels(race_bl)]

  prefix <- file.path(mi_root, paste0(config$analysis_id, "__M4_MI", mi_m))
  register_mi_prefix(prefix)
  cache_manifest_path <- paste0(prefix, "__cache_manifest.rds")
  frame_fingerprint <- fingerprint_table(frame, names(frame))
  cache_signature <- list(
    analysis_id = config$analysis_id,
    participants = nrow(frame),
    ids = sort(as.character(frame$id)),
    missing = vapply(target_variables, function(name) sum(is.na(frame[[name]])), integer(1L)),
    m = mi_m,
    maxit = mi_maxit,
    seed = mi_seed,
    frame_fingerprint = frame_fingerprint,
    mice_version = as.character(packageVersion("mice")),
    cache_version = "2026-09-08-v2"
  )
  if (file.exists(cache_manifest_path) && file.exists(paste0(prefix, "__mids.rds"))) {
    cached <- tryCatch(readRDS(cache_manifest_path), error = function(e) NULL)
    if (!is.null(cached) && identical(cached$signature, cache_signature)) {
      imp <- tryCatch(readRDS(paste0(prefix, "__mids.rds")), error = function(e) NULL)
      if (!is.null(imp) && identical(imp$m, mi_m)) {
        audit <- tryCatch(fread(paste0(prefix, "__audit.csv")), error = function(e) data.table())
        return(list(
          mids = imp,
          baseline_long = NULL,
          audit = audit,
          predictor_matrix = cached$predictor_matrix,
          methods = cached$methods,
          target_variables = target_variables,
          protected_variables = cached$protected_variables,
          prefix = prefix,
          cache_status = "reused"
        ))
      }
    }
  }

  methods <- mice::make.method(as.data.frame(frame))
  methods[] <- ""
  for (name in target_variables) {
    if (!anyNA(frame[[name]])) next
    value <- frame[[name]]
    methods[name] <- if (is.ordered(value)) {
      "polr"
    } else if (is.factor(value) && nlevels(value) == 2L) {
      "logreg"
    } else if (is.factor(value)) {
      "polyreg"
    } else {
      "pmm"
    }
  }
  predictors <- mice::make.predictorMatrix(as.data.frame(frame))
  predictors[,] <- 0L
  predictor_candidates <- setdiff(
    names(frame),
    c("id", "weight_bl")
  )
  predictor_candidates <- predictor_candidates[vapply(predictor_candidates, function(name) {
    value <- frame[[name]]
    length(unique(value[!is.na(value)])) > 1L
  }, logical(1L))]
  for (target in names(methods)[methods != ""]) {
    predictors[target, setdiff(predictor_candidates, target)] <- 1L
  }
  protected_variables <- c(
    "id", "insomnia_score_primary", "state_msm", "last_state_msm",
    "followup_years", "n_panel_records", "n_state_changes", "death_observed",
    "weight_bl"
  )
  if (any(methods[intersect(protected_variables, names(methods))] != "")) {
    stop("Protected exposure/outcome/time/weight variables were incorrectly assigned an imputation method")
  }

  imp <- mice::mice(
    data = as.data.frame(frame),
    m = mi_m,
    maxit = mi_maxit,
    method = methods,
    predictorMatrix = predictors,
    seed = mi_seed + match(config$analysis_id, vapply(analysis_configs, `[[`, character(1), "analysis_id")),
    printFlag = FALSE,
    remove.collinear = FALSE,
    remove.constant = FALSE
  )
  completed_long <- as.data.table(mice::complete(imp, action = "long", include = FALSE))
  audit <- rbindlist(lapply(names(frame), function(name) {
    data.table(
      analysis_id = config$analysis_id,
      variable = name,
      target_baseline_covariate = as.integer(name %in% target_variables),
      protected_from_imputation = as.integer(name %in% protected_variables),
      method = unname(methods[name]),
      missing_before = sum(is.na(frame[[name]])),
      missing_fraction_before = mean(is.na(frame[[name]])),
      missing_after_across_all_imputations = if (name %in% names(completed_long)) {
        sum(is.na(completed_long[[name]]))
      } else NA_integer_,
      predictor_count = if (name %in% rownames(predictors)) sum(predictors[name, ] != 0L) else 0L,
      observed_levels_or_range = if (is.factor(frame[[name]])) {
        paste(levels(frame[[name]]), collapse = " | ")
      } else {
        paste(range(as.numeric(frame[[name]]), na.rm = TRUE), collapse = " to ")
      }
    )
  }))
  predictor_long <- as.data.table(as.table(predictors))
  setnames(predictor_long, c("target", "predictor", "used"))
  predictor_long[, analysis_id := config$analysis_id]
  logged_events <- if (is.null(imp$loggedEvents)) data.table() else as.data.table(imp$loggedEvents)
  if (nrow(logged_events)) logged_events[, analysis_id := config$analysis_id]

  saveRDS(imp, paste0(prefix, "__mids.rds"), compress = "gzip")
  saveRDS(completed_long, paste0(prefix, "__baseline_long.rds"), compress = "gzip")
  fwrite(audit, paste0(prefix, "__audit.csv"), na = "")
  fwrite(predictor_long, paste0(prefix, "__predictor_matrix.csv"), na = "")
  fwrite(logged_events, paste0(prefix, "__logged_events.csv"), na = "")
  saveRDS(
    list(
      signature = cache_signature,
      predictor_matrix = predictors,
      methods = methods,
      protected_variables = protected_variables
    ),
    cache_manifest_path,
    compress = "gzip"
  )
  list(
    mids = imp,
    baseline_long = completed_long,
    audit = audit,
    predictor_matrix = predictors,
    methods = methods,
    target_variables = target_variables,
    protected_variables = protected_variables,
    prefix = prefix,
    cache_status = "created"
  )
}

rubin_pool_log_hr <- function(x, expected_m) {
  x <- x[
    is_requested_insomnia_exposure == 1L & hazard_ratio > 0 &
      is.finite(hazard_ratio) & standard_error_log_hr > 0 &
      is.finite(standard_error_log_hr)
  ]
  if (!nrow(x)) return(data.table())
  by_columns <- intersect(c(
    "analysis_id", "cohort", "transition", "from_state", "to_state",
    "term", "exposure_variable", "is_requested_insomnia_exposure", "scale"
  ), names(x))
  x[, {
    theta <- log(hazard_ratio)
    within <- standard_error_log_hr^2
    m_used <- .N
    theta_bar <- mean(theta)
    u_bar <- mean(within)
    between <- if (m_used > 1L) var(theta) else 0
    total <- u_bar + (1 + 1 / m_used) * between
    df <- if (between <= .Machine$double.eps || m_used <= 1L) {
      Inf
    } else {
      (m_used - 1) * (1 + u_bar / ((1 + 1 / m_used) * between))^2
    }
    critical <- if (is.finite(df)) qt(0.975, df = df) else qnorm(0.975)
    statistic <- theta_bar / sqrt(total)
    p_value <- if (is.finite(df)) 2 * pt(-abs(statistic), df = df) else 2 * pnorm(-abs(statistic))
    list(
      m_expected = expected_m,
      m_fitted = m_used,
      pooling_status = if (m_used == expected_m) "RUBIN_POOLED_COMPLETE" else "PARTIAL_MI_NOT_PRIMARY",
      log_hazard_ratio = theta_bar,
      hazard_ratio = exp(theta_bar),
      conf_low = exp(theta_bar - critical * sqrt(total)),
      conf_high = exp(theta_bar + critical * sqrt(total)),
      standard_error_log_hr = sqrt(total),
      within_imputation_variance = u_bar,
      between_imputation_variance = between,
      fraction_missing_information = if (total > 0) ((1 + 1 / m_used) * between) / total else 0,
      degrees_freedom = df,
      statistic = statistic,
      p_value = p_value
    )
  }, by = by_columns]
}


sample_diagnostics <- function(original, sample, edges, analysis_id, model_id, covariates) {
  d <- copy(sample)
  setorder(d, id, time_msm)
  d[, next_state := shift(state_msm, type = "lead"), by = id]
  d[, next_time := shift(time_msm, type = "lead"), by = id]
  pairs <- d[!is.na(next_state) & next_time > time_msm]
  edge_counts <- copy(edges)
  edge_counts[, `:=`(
    risk_intervals = vapply(seq_len(.N), function(i) {
      pairs[state_msm == edge_counts$from_index[i], .N]
    }, integer(1)),
    risk_participants = vapply(seq_len(.N), function(i) {
      pairs[state_msm == edge_counts$from_index[i], uniqueN(id)]
    }, integer(1)),
    person_years = vapply(seq_len(.N), function(i) {
      pairs[state_msm == edge_counts$from_index[i], sum(next_time - time_msm)]
    }, numeric(1)),
    events = vapply(seq_len(.N), function(i) {
      pairs[state_msm == edge_counts$from_index[i] & next_state == edge_counts$to_index[i], .N]
    }, integer(1))
  )]
  edge_counts[, `:=`(
    analysis_id = analysis_id,
    model_id = model_id,
    n_sample_rows = nrow(sample),
    n_sample_participants = uniqueN(sample$id)
  )]
  summary <- data.table(
    analysis_id = analysis_id,
    model_id = model_id,
    n_input_rows = nrow(original),
    n_input_participants = uniqueN(original$id),
    n_complete_case_rows = nrow(sample),
    n_complete_case_participants = uniqueN(sample$id),
    n_observed_changes = pairs[state_msm != next_state, .N],
    n_observed_deaths = pairs[next_state == 5L, .N],
    n_observed_s0_s3_multistep = pairs[state_msm == 1L & next_state == 4L, .N],
    min_allowed_edge_events = if (nrow(edge_counts)) min(edge_counts$events) else NA_integer_,
    max_followup_years = if (nrow(sample)) max(sample$time_msm) else NA_real_,
    insomnia_mean = if (nrow(sample)) mean(sample$insomnia_score_primary) else NA_real_,
    insomnia_sd = if (nrow(sample)) sd(sample$insomnia_score_primary) else NA_real_,
    covariates_used = paste(covariates, collapse = " + ")
  )
  list(summary = summary, edge_counts = edge_counts)
}

make_nuisance_constraints <- function(covariate_formula, d, edges, exposure) {
  mm <- model.matrix(covariate_formula, data = d)
  coefficient_names <- setdiff(colnames(mm), "(Intercept)")
  exposure_names <- coefficient_names[startsWith(coefficient_names, exposure)]
  nuisance_names <- setdiff(coefficient_names, exposure_names)
  constraints <- list()
  # Explicit 1:K indices preserve one insomnia coefficient per allowed edge.
  for (name in exposure_names) constraints[[name]] <- seq_len(nrow(edges))
  for (name in nuisance_names) {
    constraints[[name]] <- if (nuisance_constraint_mode == "common") {
      rep(1L, nrow(edges))
    } else {
      seq_len(nrow(edges))
    }
  }
  if (!length(constraints)) NULL else constraints
}

fit_metrics <- function(
    fit, analysis_id, model_id, engine, sample, elapsed_seconds,
    warnings, requested_covariates, used_covariates, omitted_covariates) {
  ll <- tryCatch(logLik(fit), error = function(e) NA)
  ll_value <- suppressWarnings(as.numeric(ll))
  df <- suppressWarnings(as.numeric(attr(ll, "df")))
  n_subjects <- uniqueN(sample$id)
  n_records <- nrow(sample)
  aic <- if (is.finite(ll_value) && is.finite(df)) -2 * ll_value + 2 * df else NA_real_
  bic_subject <- if (is.finite(ll_value) && is.finite(df) && n_subjects > 0) {
    -2 * ll_value + log(n_subjects) * df
  } else NA_real_
  bic_record <- if (is.finite(ll_value) && is.finite(df) && n_records > 0) {
    -2 * ll_value + log(n_records) * df
  } else NA_real_

  optim_code <- tryCatch(as.integer(fit$opt$convergence), error = function(e) NA_integer_)
  if (!length(optim_code)) optim_code <- NA_integer_ else optim_code <- optim_code[1L]
  optim_message <- tryCatch(as.character(fit$opt$message), error = function(e) NA_character_)
  if (!length(optim_message)) optim_message <- NA_character_ else optim_message <- optim_message[1L]
  found_se <- isTRUE(fit$foundse)
  hessian <- NULL
  if (is.matrix(fit$opt$hessian)) hessian <- fit$opt$hessian
  if (is.null(hessian) && is.matrix(fit$hessian)) hessian <- fit$hessian
  hessian_available <- is.matrix(hessian)
  hessian_finite <- hessian_available && all(is.finite(hessian))
  hessian_min_eigen <- NA_real_
  hessian_positive_definite <- FALSE
  if (hessian_finite) {
    hessian_min_eigen <- tryCatch(
      min(eigen((hessian + t(hessian)) / 2, symmetric = TRUE, only.values = TRUE)$values),
      error = function(e) NA_real_
    )
    hessian_positive_definite <- is.finite(hessian_min_eigen) && hessian_min_eigen > 1e-8
  }
  covmat_finite <- is.matrix(fit$covmat) && all(is.finite(fit$covmat))
  estimated_q <- tryCatch(msm::qmatrix.msm(fit, ci = "none"), error = function(e) NULL)
  q_row_sum_max_abs <- if (is.matrix(estimated_q)) {
    max(abs(rowSums(estimated_q)))
  } else NA_real_
  absorbing_row_max_abs <- if (is.matrix(estimated_q)) {
    max(abs(estimated_q[nrow(estimated_q), ]))
  } else NA_real_
  death_row_max_abs <- if (is.matrix(estimated_q) && nrow(estimated_q) == 5L) {
    max(abs(estimated_q[5L, ]))
  } else NA_real_
  absorbing_state_ok <- is.finite(absorbing_row_max_abs) && absorbing_row_max_abs < 1e-10
  model_data_fingerprint <- fingerprint_table(
    sample,
    unique(c("id", "time_msm", "state_msm", used_covariates))
  )
  status <- if (identical(optim_code, 0L) && found_se && covmat_finite) {
    "FITTED"
  } else {
    "FITTED_WITH_INFERENCE_WARNING"
  }
  data.table(
    script_version = script_version,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    analysis_id = analysis_id,
    cohort = if ("cohort" %in% names(sample)) as.character(sample$cohort[1L]) else NA_character_,
    model_id = model_id,
    engine = engine,
    status = status,
    n_records = n_records,
    n_subjects = n_subjects,
    n_parameters = df,
    log_likelihood = ll_value,
    minus2_log_likelihood = -2 * ll_value,
    AIC = aic,
    BIC_subjects = bic_subject,
    BIC_records = bic_record,
    optimizer_convergence_code = optim_code,
    optimizer_message = optim_message,
    standard_errors_found = as.integer(found_se),
    covariance_matrix_finite = as.integer(covmat_finite),
    hessian_available = as.integer(hessian_available),
    hessian_finite = as.integer(hessian_finite),
    hessian_positive_definite = as.integer(hessian_positive_definite),
    hessian_min_eigenvalue = hessian_min_eigen,
    q_row_sum_max_abs = q_row_sum_max_abs,
    death_row_max_abs = death_row_max_abs,
    absorbing_row_max_abs = absorbing_row_max_abs,
    absorbing_state_ok = as.integer(absorbing_state_ok),
    model_data_fingerprint = model_data_fingerprint,
    elapsed_seconds = elapsed_seconds,
    warning_count = length(warnings),
    requested_covariates = paste(requested_covariates, collapse = " + "),
    used_covariates = paste(used_covariates, collapse = " + "),
    omitted_covariates = paste(omitted_covariates, collapse = " + ")
  )
}

extract_q_table <- function(fit, edges, analysis_id, model_id, engine) {
  q <- msm::qmatrix.msm(fit, covariates = "mean", ci = "delta")
  get_value <- function(component, i, j) {
    value <- q[[component]]
    if (is.null(value)) NA_real_ else as.numeric(value[i, j])
  }
  out <- copy(edges)
  out[, `:=`(
    analysis_id = analysis_id,
    model_id = model_id,
    engine = engine,
    covariate_profile = "sample means",
    intensity = vapply(seq_len(.N), function(k) get_value("estimates", from_index[k], to_index[k]), numeric(1)),
    standard_error = vapply(seq_len(.N), function(k) get_value("SE", from_index[k], to_index[k]), numeric(1)),
    conf_low = vapply(seq_len(.N), function(k) get_value("L", from_index[k], to_index[k]), numeric(1)),
    conf_high = vapply(seq_len(.N), function(k) get_value("U", from_index[k], to_index[k]), numeric(1))
  )]
  setcolorder(out, c(
    "analysis_id", "model_id", "engine", "transition", "from_state", "to_state",
    "from_index", "to_index", "death_transition", "covariate_profile",
    "intensity", "standard_error", "conf_low", "conf_high"
  ))
  out[]
}

extract_hazard_table <- function(fit, edges, analysis_id, model_id, engine, exposure) {
  hazards <- msm::hazard.msm(fit)
  if (!length(hazards)) return(data.table())
  rows <- list()
  for (term in names(hazards)) {
    mat <- as.matrix(hazards[[term]])
    n_take <- min(nrow(mat), nrow(edges))
    if (!n_take) next
    edge <- edges[seq_len(n_take)]
    get_col <- function(name, fallback) {
      if (name %in% colnames(mat)) as.numeric(mat[seq_len(n_take), name]) else rep(fallback, n_take)
    }
    hr <- get_col("HR", NA_real_)
    low <- get_col("L", NA_real_)
    high <- get_col("U", NA_real_)
    se_log <- (log(high) - log(low)) / (2 * qnorm(0.975))
    z <- log(hr) / se_log
    p <- 2 * pnorm(-abs(z))
    p[!is.finite(p)] <- NA_real_
    rows[[length(rows) + 1L]] <- data.table(
      analysis_id = analysis_id,
      model_id = model_id,
      engine = engine,
      transition = edge$transition,
      from_state = edge$from_state,
      to_state = edge$to_state,
      term = term,
      exposure_variable = exposure,
      is_requested_insomnia_exposure = as.integer(startsWith(term, exposure)),
      scale = fifelse(
        term == "insomnia_score_primary", "per one additional symptom",
        fifelse(
          startsWith(term, "insomnia_burden_factor"), "category versus 0 symptoms",
          fifelse(
            term %chin% c(
              "insomnia_onset_primary", "insomnia_maintenance_primary",
              "insomnia_nonrestorative_primary"
            ),
            "symptom present versus absent",
            fifelse(term == "age_bl_5", "per 5 years", fifelse(term == "bmi_bl_5", "per 5 kg/m2", "as coded"))
          )
        )
      ),
      hazard_ratio = hr,
      conf_low = low,
      conf_high = high,
      standard_error_log_hr = se_log,
      z_value = z,
      p_value = p
    )
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

confirmatory_parameter_block <- function(fit, sample, edges, covariates, exposure) {
  mm <- model.matrix(reformulate(covariates), data = sample)
  coefficient_names <- setdiff(colnames(mm), "(Intercept)")
  exposure_columns <- which(startsWith(coefficient_names, exposure))
  if (!length(exposure_columns)) stop("No fitted coefficient belongs to the requested exposure")
  n_edges <- nrow(edges)
  nonfatal_positions <- which(edges$death_transition == 0L)
  indices <- unlist(lapply(exposure_columns, function(column_index) {
    n_edges + (column_index - 1L) * n_edges + nonfatal_positions
  }), use.names = FALSE)
  theta_all <- as.numeric(fit$estimates)
  covariance_all <- as.matrix(fit$covmat)
  if (any(indices > length(theta_all)) || any(indices > nrow(covariance_all))) {
    stop("Unable to map the confirmatory exposure block to the msm covariance matrix")
  }
  list(
    theta = theta_all[indices],
    covariance = covariance_all[indices, indices, drop = FALSE],
    parameter_labels = unlist(lapply(exposure_columns, function(column_index) {
      paste0(coefficient_names[column_index], "__", edges$transition[nonfatal_positions])
    }), use.names = FALSE),
    exposure_terms = coefficient_names[exposure_columns],
    transitions = edges$transition[nonfatal_positions]
  )
}

wald_from_block <- function(block) {
  covariance <- (block$covariance + t(block$covariance)) / 2
  if (!all(is.finite(block$theta)) || !all(is.finite(covariance))) {
    stop("Confirmatory Wald block contains non-finite estimates or covariance")
  }
  decomposition <- eigen(covariance, symmetric = TRUE)
  tolerance <- max(abs(decomposition$values)) * max(dim(covariance)) * .Machine$double.eps^0.5
  keep <- decomposition$values > tolerance
  rank <- sum(keep)
  if (rank < 1L) stop("Confirmatory Wald covariance has rank zero")
  projected <- as.numeric(crossprod(decomposition$vectors[, keep, drop = FALSE], block$theta))
  statistic <- sum(projected^2 / decomposition$values[keep])
  list(statistic = statistic, df = rank, p_value = pchisq(statistic, df = rank, lower.tail = FALSE))
}

extract_confirmatory_wald <- function(
    fit, sample, edges, covariates, exposure, analysis_id, model_id, engine) {
  block <- confirmatory_parameter_block(fit, sample, edges, covariates, exposure)
  test <- wald_from_block(block)
  data.table(
    analysis_id = analysis_id,
    model_id = model_id,
    engine = engine,
    exposure_variable = exposure,
    tested_terms = paste(block$exposure_terms, collapse = " | "),
    tested_transitions = paste(unique(block$transitions), collapse = " | "),
    statistic = test$statistic,
    df = test$df,
    p_value = test$p_value,
    method = "Joint Wald test using the fitted msm covariance matrix",
    confirmatory_edge_count = uniqueN(block$transitions)
  )
}

profile_at_nuisance_means <- function(fit, exposure, value) {
  covariate_means <- fit$qcmodel$covmeans
  if (is.null(covariate_means) || !length(covariate_means)) {
    stop("The fitted model does not expose covariate means for adjusted prediction")
  }
  profile <- as.list(setNames(as.numeric(covariate_means), names(covariate_means)))
  exposure_columns <- names(profile)[startsWith(names(profile), exposure)]
  if (!length(exposure_columns)) {
    stop("Exposure is absent from the fitted model matrix: ", exposure)
  }
  if (length(exposure_columns) == 1L && identical(exposure_columns, exposure)) {
    profile[[exposure]] <- as.numeric(value)
  } else {
    profile[exposure_columns] <- as.list(rep(0, length(exposure_columns)))
    if (!identical(as.character(value), "0 symptoms")) {
      target <- paste0(exposure, as.character(value))
      if (!target %in% exposure_columns) {
        stop(
          "Unable to map factor exposure level '", value,
          "' to fitted model-matrix columns for ", exposure
        )
      }
      profile[[target]] <- 1
    }
  }
  profile
}

make_prediction_profiles <- function(fit, exposure) {
  profiles <- list(sample_means = "mean")
  if (exposure == "insomnia_score_primary") {
    for (value in 0:3) {
      profiles[[paste0("insomnia_", value)]] <- profile_at_nuisance_means(
        fit, exposure, value
      )
    }
  } else if (exposure == "insomnia_burden_factor") {
    values <- c("0 symptoms", "1 symptom", ">=2 symptoms")
    labels <- c("burden_0", "burden_1", "burden_ge2")
    for (i in seq_along(values)) {
      profiles[[labels[i]]] <- profile_at_nuisance_means(fit, exposure, values[i])
    }
  } else {
    profiles$symptom_absent <- profile_at_nuisance_means(fit, exposure, 0)
    profiles$symptom_present <- profile_at_nuisance_means(fit, exposure, 1)
  }
  profiles
}

extract_probability_table <- function(fit, edges, analysis_id, model_id, engine, max_followup, exposure) {
  profiles <- make_prediction_profiles(fit, exposure)
  allowed_keys <- paste(edges$from_index, edges$to_index)
  rows <- list()
  for (profile_name in names(profiles)) {
    for (time in c(2, 5, 10)) {
      p <- msm::pmatrix.msm(
        fit, t = time, covariates = profiles[[profile_name]], ci = "none"
      )
      index <- expand.grid(
        from_index = seq_len(nrow(p)),
        to_index = seq_len(ncol(p)),
        KEEP.OUT.ATTRS = FALSE
      )
      z <- as.data.table(index)
      z[, probability := as.numeric(p[cbind(from_index, to_index)])]
      z[, `:=`(
        analysis_id = analysis_id,
        model_id = model_id,
        engine = engine,
        exposure_variable = exposure,
        profile = profile_name,
        time_years = time,
        from_state = from_index - 1L,
        to_state = to_index - 1L,
        transition = paste0(state_name(from_index), "->", state_name(to_index)),
        extrapolated_beyond_observed_followup = as.integer(time > max_followup),
        nuisance_profile = "All fitted nuisance model-matrix covariates held at sample means",
        uncertainty_scope = "Point estimate only; ci='none'; not a 95% confidence interval",
        quantity_type = "Finite-horizon endpoint-state probability",
        direct_q_edge = as.integer(paste(from_index, to_index) %chin% allowed_keys),
        pathway_note = fifelse(
          from_index == 1L & to_index == 4L,
          "S0-to-S3 endpoint probability arises only through S1 or S2; direct Q intensity is fixed to zero",
          "Endpoint probability may include zero, one, or multiple latent transitions"
        )
      )]
      row_error <- z[, .(row_sum_error = abs(sum(probability) - 1)), by = from_index]
      z <- merge(z, row_error, by = "from_index", all.x = TRUE, sort = FALSE)
      rows[[length(rows) + 1L]] <- z
    }
  }
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  setcolorder(out, c(
    "analysis_id", "model_id", "engine", "exposure_variable", "profile", "time_years", "from_state",
    "to_state", "from_index", "to_index", "transition", "probability",
    "row_sum_error", "extrapolated_beyond_observed_followup", "nuisance_profile",
    "uncertainty_scope", "quantity_type", "direct_q_edge", "pathway_note"
  ))
  out[]
}

extract_sojourn_table <- function(fit, analysis_id, model_id, engine) {
  s <- as.data.table(msm::sojourn.msm(fit, covariates = "mean", ci = "delta"), keep.rownames = "state_raw")
  setnames(s, old = intersect(names(s), c("estimates", "SE", "L", "U")), new = c(
    estimates = "sojourn_years", SE = "standard_error", L = "conf_low", U = "conf_high"
  )[intersect(names(s), c("estimates", "SE", "L", "U"))])
  s[, state_index := seq_len(.N)]
  s[, state := state_index - 1L]
  s[, state_label := state_name(state_index)]
  s[, `:=`(analysis_id = analysis_id, model_id = model_id, engine = engine)]
  setcolorder(s, c(
    "analysis_id", "model_id", "engine", "state", "state_index", "state_label",
    setdiff(names(s), c("analysis_id", "model_id", "engine", "state", "state_index", "state_label"))
  ))
  s[]
}

extract_expected_state_time <- function(
    fit, analysis_id, model_id, engine, max_followup, exposure, horizon = 10) {
  profiles <- make_prediction_profiles(fit, exposure)
  rows <- list()
  for (profile_name in names(profiles)) {
    time_in_state <- as.numeric(msm::totlos.msm(
      fit, start = 1L, tot = horizon, covariates = profiles[[profile_name]], ci = "none"
    ))
    n_states <- length(time_in_state)
    for (j in seq_len(n_states)) {
      rows[[length(rows) + 1L]] <- data.table(
        analysis_id = analysis_id,
        model_id = model_id,
        engine = engine,
        exposure_variable = exposure,
        profile = profile_name,
        start_state = "S0",
        horizon_years = horizon,
        state = j - 1L,
        state_label = state_name(j),
        expected_years = time_in_state[j],
        summary_measure = "state-specific expected years",
        extrapolated_beyond_observed_followup = as.integer(horizon > max_followup),
        nuisance_profile = "All fitted nuisance model-matrix covariates held at sample means",
        uncertainty_scope = "Point estimate only; ci='none'; not a 95% confidence interval"
      )
    }
    # In the five-state model this excludes both multimorbidity (S3) and death
    # (S4). In the four-state nonfatal model, mortality is not represented.
    eligible <- seq_len(min(3L, n_states))
    s3_free_label <- if (n_states == 5L) {
      "Alive without S3"
    } else {
      "Without S3 (mortality not represented)"
    }
    s3_free_measure <- if (n_states == 5L) {
      "expected alive and S3-free years"
    } else {
      "expected S3-free years in the nonfatal model"
    }
    rows[[length(rows) + 1L]] <- data.table(
      analysis_id = analysis_id,
      model_id = model_id,
      engine = engine,
      exposure_variable = exposure,
      profile = profile_name,
      start_state = "S0",
      horizon_years = horizon,
      state = NA_integer_,
      state_label = s3_free_label,
      expected_years = sum(time_in_state[eligible]),
      summary_measure = s3_free_measure,
      extrapolated_beyond_observed_followup = as.integer(horizon > max_followup),
      nuisance_profile = "All fitted nuisance model-matrix covariates held at sample means",
      uncertainty_scope = "Point estimate only; ci='none'; not a 95% confidence interval"
    )
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

extract_prevalence_table <- function(fit, sample, analysis_id, model_id, engine, max_followup) {
  times <- unique(c(0, c(2, 5, 10)[c(2, 5, 10) <= max_followup], max_followup))
  times <- sort(times[is.finite(times) & times >= 0])
  d <- copy(as.data.table(sample)[, .(id, time_msm, state_msm)])
  setorder(d, id, time_msm, state_msm)
  subject_span <- d[, .(
    entry_time = time_msm[1L],
    entry_state = state_msm[1L],
    last_time = time_msm[.N],
    last_state = state_msm[.N]
  ), by = id]
  fitted_q <- msm::qmatrix.msm(fit, ci = "none")
  n_states <- nrow(fitted_q)
  absorbing_states <- which(rowSums(abs(fitted_q)) < 1e-12)
  timezero <- min(subject_span$entry_time)
  initial_states <- subject_span[abs(entry_time - timezero) < 1e-10, entry_state]
  initial_count <- tabulate(initial_states, nbins = n_states)
  initial_probability <- initial_count / sum(initial_count)
  rows <- list()
  for (time in times) {
    # This follows msm::observed.msm's start-of-interval LOCF convention:
    # participants are censored after their final transient-state record, while
    # any model-defined absorbing state persists after its observation.
    snap <- d[time_msm <= time, .SD[.N], by = id]
    snap <- merge(snap, subject_span, by = "id", all.x = TRUE, sort = FALSE)
    snap <- snap[last_time >= time | state_msm %in% absorbing_states]
    if (!nrow(snap)) next
    observed_count <- tabulate(snap$state_msm, nbins = n_states)
    transition_matrix <- msm::pmatrix.msm(
      fit, t = max(time - timezero, 0),
      covariates = "mean", ci = "none"
    )
    expected_probability <- as.numeric(initial_probability %*% transition_matrix)
    total_observed <- sum(observed_count)
    total_expected <- total_observed
    for (j in seq_len(n_states)) {
      rows[[length(rows) + 1L]] <- data.table(
        analysis_id = analysis_id,
        model_id = model_id,
        engine = engine,
        diagnostic_method = paste(
          "LOCF observed risk set with absorbing-state carry-forward versus",
          "fixed baseline state mix and mean-covariate Markov expectation"
        ),
        calibration_scope = paste(
          "Approximate prevalence diagnostic using current risk-set size;",
          "not an exposure-standardized causal occupancy estimate"
        ),
        uncertainty_scope = "Point estimate only; ci='none'; not a 95% confidence interval",
        time_years = time,
        state = j - 1L,
        state_index = j,
        state_label = state_name(j),
        observed_count = observed_count[j],
        expected_count = expected_probability[j] * total_expected,
        observed_probability = observed_count[j] / total_observed,
        expected_probability = expected_probability[j],
        total_observed = as.numeric(total_observed),
        total_expected = as.numeric(total_expected)
      )
    }
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

write_model_table <- function(prefix, suffix, value) {
  if (is.null(value)) return(invisible(NULL))
  path <- file.path(by_model_root, paste0(prefix, "_", suffix, ".csv"))
  fwrite(as.data.table(value), path, na = "")
  register_run_table(suffix, path)
  invisible(path)
}

extract_and_write <- function(
    fit, analysis_id, model_id, engine, edges, sample, sample_info,
    elapsed_seconds, warnings, requested_covariates, used_covariates,
    omitted_covariates, log_file, exposure) {
  prefix <- paste(analysis_id, model_id, sep = "__")
  metrics <- fit_metrics(
    fit, analysis_id, model_id, engine, sample, elapsed_seconds, warnings,
    requested_covariates, used_covariates, omitted_covariates
  )
  write_model_table(prefix, "fit_metrics", metrics)
  write_model_table(prefix, "sample_summary", sample_info$summary)
  write_model_table(prefix, "edge_counts", sample_info$edge_counts)

  extractors <- list(
    q_intensities = function() extract_q_table(fit, edges, analysis_id, model_id, engine),
    hazard_ratios = function() extract_hazard_table(
      fit, edges, analysis_id, model_id, engine, exposure
    ),
    confirmatory_wald = function() extract_confirmatory_wald(
      fit, sample, edges, used_covariates, exposure,
      analysis_id, model_id, engine
    ),
    transition_probabilities = function() extract_probability_table(
      fit, edges, analysis_id, model_id, engine, max(sample$time_msm), exposure
    ),
    sojourn_times = function() extract_sojourn_table(fit, analysis_id, model_id, engine),
    expected_state_time = function() extract_expected_state_time(
      fit, analysis_id, model_id, engine, max(sample$time_msm), exposure
    ),
    observed_expected_prevalence = function() extract_prevalence_table(
      fit, sample, analysis_id, model_id, engine, max(sample$time_msm)
    )
  )
  for (name in names(extractors)) {
    result <- tryCatch(extractors[[name]](), error = function(e) e)
    if (inherits(result, "error")) {
      record_failure(
        analysis_id, model_id, engine, paste0("extract_", name), "PARTIAL_FAILURE",
        conditionMessage(result), warnings, elapsed_seconds, log_file
      )
    } else {
      if (!"cohort" %in% names(result)) {
        result[, cohort := if ("cohort" %in% names(sample)) {
          as.character(sample$cohort[1L])
        } else NA_character_]
      }
      if (name %in% c("q_intensities", "hazard_ratios") && nrow(result)) {
        support <- unique(sample_info$edge_counts[, .(
          transition, risk_intervals, risk_participants, person_years, events
        )])
        result <- merge(result, support, by = "transition", all.x = TRUE, sort = FALSE)
        result[, `:=`(
          cohort = if ("cohort" %in% names(sample)) as.character(sample$cohort[1L]) else NA_character_,
          n_parameters = metrics$n_parameters[1L],
          events_per_total_model_parameter = events / metrics$n_parameters[1L],
          events_per_parameter = events / metrics$n_parameters[1L],
          events_per_parameter_scope = paste(
            "Events divided by total fitted-model degrees of freedom;",
            "a conservative model-level support flag, not transition-specific EPP"
          )
        )]
      }
      write_model_table(prefix, name, result)
    }
  }
  metrics$status[1L]
}

rubin_pool_positive <- function(x, estimate_column, se_column, by_columns, expected_m) {
  x <- x[
    get(estimate_column) > 0 & is.finite(get(estimate_column)) &
      get(se_column) > 0 & is.finite(get(se_column))
  ]
  if (!nrow(x)) return(data.table())
  x[, {
    estimate <- get(estimate_column)
    se <- get(se_column)
    theta <- log(estimate)
    within <- (se / estimate)^2
    m_used <- .N
    theta_bar <- mean(theta)
    u_bar <- mean(within)
    between <- if (m_used > 1L) var(theta) else 0
    total <- u_bar + (1 + 1 / m_used) * between
    df <- if (between <= .Machine$double.eps || m_used <= 1L) {
      Inf
    } else {
      (m_used - 1) * (1 + u_bar / ((1 + 1 / m_used) * between))^2
    }
    critical <- if (is.finite(df)) qt(0.975, df = df) else qnorm(0.975)
    list(
      m_expected = expected_m,
      m_fitted = m_used,
      pooling_status = if (m_used == expected_m) "RUBIN_POOLED_COMPLETE" else "PARTIAL_MI_NOT_PRIMARY",
      pooled_estimate = exp(theta_bar),
      pooled_standard_error_log = sqrt(total),
      conf_low = exp(theta_bar - critical * sqrt(total)),
      conf_high = exp(theta_bar + critical * sqrt(total)),
      within_imputation_variance = u_bar,
      between_imputation_variance = between,
      degrees_freedom = df
    )
  }, by = by_columns]
}

fit_mi_primary <- function(mi_object, x, config, edges, q_initial) {
  analysis_id <- config$analysis_id
  model_id <- "M4_MI20"
  engine <- "mice::mice + msm::msm + Rubin pooling"
  requested <- c("insomnia_score_primary", primary_m4_adjustment)
  component_metrics <- list()
  component_q <- list()
  component_hr <- list()
  component_p <- list()
  component_sojourn <- list()
  component_prevalence <- list()
  component_expected_time <- list()
  component_confirmatory_blocks <- list()
  fitted_files <- character()
  mi_indices <- seq.int(mi_start, mi_end)

  for (imputation in mi_indices) {
    imputation_id <- sprintf("imp%02d", imputation)
    component_model_id <- paste(model_id, imputation_id, sep = "__")
    completed <- as.data.table(mice::complete(mi_object$mids, action = imputation))
    imputed_variables <- intersect(primary_m4_adjustment, names(completed))
    d <- copy(x)
    d[, (imputed_variables) := NULL]
    d <- merge(
      d,
      completed[, c("id", imputed_variables), with = FALSE],
      by = "id", all.x = TRUE, sort = FALSE
    )
    sample_object <- tryCatch(
      prepare_model_sample(d, requested, exposure = "insomnia_score_primary"),
      error = function(e) e
    )
    if (inherits(sample_object, "error")) {
      record_failure(
        analysis_id, component_model_id, engine, "mi_preflight", "NOT_ESTIMABLE",
        conditionMessage(sample_object)
      )
      next
    }
    d <- sample_object$data
    covariates <- sample_object$used
    cov_formula <- reformulate(covariates)
    constraints <- make_nuisance_constraints(
      cov_formula, d, edges, exposure = "insomnia_score_primary"
    )
    model_data_fingerprint <- fingerprint_table(
      d, unique(c("id", "time_msm", "state_msm", covariates))
    )
    spec_hash <- digest::digest(list(
      script_version, analysis_id, model_id, imputation, covariates,
      q_initial, nuisance_constraint_mode, config$death_exact_state,
      max_iterations, model_data_fingerprint
    ), algo = "xxhash64")
    model_file <- file.path(model_root, paste0(
      analysis_id, "__", model_id, "__", imputation_id, ".rds"
    ))
    envelope <- NULL
    if (!force_refit && file.exists(model_file)) {
      cached <- tryCatch(readRDS(model_file), error = function(e) NULL)
      if (!is.null(cached) && identical(cached$spec_hash, spec_hash)) envelope <- cached
    }
    if (is.null(envelope)) {
      log_file <- file.path(log_root, paste0(
        analysis_id, "__", model_id, "__", imputation_id, ".log"
      ))
      run <- safe_run(
        fun = function() msm::msm(
          formula = state_msm ~ time_msm,
          subject = id,
          data = as.data.frame(d),
          qmatrix = q_initial,
          covariates = cov_formula,
          constraint = constraints,
          deathexact = config$death_exact_state,
          exacttimes = FALSE,
          center = TRUE,
          hessian = TRUE,
          control = list(fnscale = 10000, maxit = max_iterations, trace = 0)
        ),
        log_file = log_file,
        analysis_id = analysis_id,
        model_id = component_model_id,
        engine = engine,
        stage = "mi_fit"
      )
      if (is.null(run$value)) next
      envelope <- list(
        spec_hash = spec_hash,
        fit = run$value,
        imputation = imputation,
        model_data_fingerprint = model_data_fingerprint,
        used_covariates = covariates,
        constraints = constraints,
        warnings = run$warnings,
        elapsed_seconds = run$elapsed_seconds,
        log_file = log_file
      )
      saveRDS(envelope, model_file, compress = "xz")
    }
    fit <- envelope$fit
    fitted_files <- c(fitted_files, model_file)
    metric <- fit_metrics(
      fit, analysis_id, component_model_id, engine, d,
      envelope$elapsed_seconds, envelope$warnings, requested, covariates,
      c(sample_object$missing_columns, sample_object$constant_columns)
    )
    metric[, imputation := imputation]
    component_metrics[[length(component_metrics) + 1L]] <- metric
    collect <- function(fun, target) {
      value <- tryCatch(fun(), error = function(e) e)
      if (inherits(value, "error")) {
        record_failure(
          analysis_id, component_model_id, engine, target, "PARTIAL_FAILURE",
          conditionMessage(value), envelope$warnings, envelope$elapsed_seconds,
          envelope$log_file
        )
        return(NULL)
      }
      value[, `:=`(imputation = imputation, cohort = config$cohort)]
      value
    }
    component_q[[length(component_q) + 1L]] <- collect(
      function() extract_q_table(fit, edges, analysis_id, component_model_id, engine),
      "mi_extract_q"
    )
    component_hr[[length(component_hr) + 1L]] <- collect(
      function() extract_hazard_table(
        fit, edges, analysis_id, component_model_id, engine, "insomnia_score_primary"
      ),
      "mi_extract_hr"
    )
    component_p[[length(component_p) + 1L]] <- collect(
      function() extract_probability_table(
        fit, edges, analysis_id, component_model_id, engine,
        max(d$time_msm), "insomnia_score_primary"
      ),
      "mi_extract_probability"
    )
    component_sojourn[[length(component_sojourn) + 1L]] <- collect(
      function() extract_sojourn_table(fit, analysis_id, component_model_id, engine),
      "mi_extract_sojourn"
    )
    component_prevalence[[length(component_prevalence) + 1L]] <- collect(
      function() extract_prevalence_table(
        fit, d, analysis_id, component_model_id, engine, max(d$time_msm)
      ),
      "mi_extract_prevalence"
    )
    component_expected_time[[length(component_expected_time) + 1L]] <- collect(
      function() extract_expected_state_time(
        fit, analysis_id, component_model_id, engine,
        max(d$time_msm), "insomnia_score_primary"
      ),
      "mi_extract_expected_state_time"
    )
    block <- tryCatch(
      confirmatory_parameter_block(
        fit, d, edges, covariates, "insomnia_score_primary"
      ),
      error = function(e) e
    )
    if (inherits(block, "error")) {
      record_failure(
        analysis_id, component_model_id, engine, "mi_confirmatory_block",
        "PARTIAL_FAILURE", conditionMessage(block), envelope$warnings,
        envelope$elapsed_seconds, envelope$log_file
      )
    } else {
      block$imputation <- imputation
      component_confirmatory_blocks[[length(component_confirmatory_blocks) + 1L]] <- block
    }
    rm(d, completed, fit)
    gc(verbose = FALSE)
  }

  bind_nonnull <- function(values) {
    values <- Filter(Negate(is.null), values)
    if (length(values)) rbindlist(values, use.names = TRUE, fill = TRUE) else data.table()
  }
  metrics <- bind_nonnull(component_metrics)
  q_values <- bind_nonnull(component_q)
  hr_values <- bind_nonnull(component_hr)
  p_values <- bind_nonnull(component_p)
  sojourn_values <- bind_nonnull(component_sojourn)
  prevalence_values <- bind_nonnull(component_prevalence)
  expected_time_values <- bind_nonnull(component_expected_time)
  prefix <- paste(analysis_id, model_id, sep = "__")
  component_prefix <- paste0(
    prefix, "__imp", sprintf("%02d", min(mi_indices)), "_", sprintf("%02d", max(mi_indices))
  )
  write_model_table(component_prefix, "mi_components_metrics", metrics)
  write_model_table(component_prefix, "mi_components_q", q_values)
  write_model_table(component_prefix, "mi_components_hr", hr_values)
  write_model_table(component_prefix, "mi_components_pmatrix", p_values)
  write_model_table(component_prefix, "mi_components_sojourn", sojourn_values)
  write_model_table(component_prefix, "mi_components_prevalence", prevalence_values)
  write_model_table(component_prefix, "mi_components_expected_time", expected_time_values)

  fitted_m <- if (nrow(metrics)) uniqueN(metrics$imputation) else 0L
  valid_imputations <- if (nrow(metrics)) {
    metrics[
      status == "FITTED" & optimizer_convergence_code == 0L &
        standard_errors_found == 1L & covariance_matrix_finite == 1L &
        hessian_positive_definite == 1L,
      sort(unique(imputation))
    ]
  } else integer()
  valid_m <- length(valid_imputations)
  if (!mi_pool) {
    component_status <- if (valid_m == length(mi_indices)) {
      "MI_COMPONENTS_FITTED_AWAITING_POOL"
    } else {
      "PARTIAL_MI_COMPONENTS_AWAITING_POOL"
    }
    saveRDS(
      list(
        analysis_id = analysis_id,
        model_id = model_id,
        component_status = component_status,
        requested_imputations = mi_indices,
        fitted_imputations = sort(unique(metrics$imputation)),
        valid_imputations = valid_imputations,
        model_files = fitted_files,
        mi_data_prefix = mi_object$prefix
      ),
      file.path(model_root, paste0(component_prefix, "__manifest.rds")),
      compress = "xz"
    )
    return(component_status)
  }
  pooling_status <- if (valid_m == mi_m) "RUBIN_POOLED_COMPLETE" else "PARTIAL_MI_NOT_PRIMARY"
  if (valid_m != mi_m) {
    record_failure(
      analysis_id, model_id, engine, "mi_pool", pooling_status,
      paste(
        "Only", valid_m, "of", mi_m,
        "imputation-specific MSMs passed convergence, covariance, SE and Hessian checks;",
        "canonical pooled primary estimates were not written"
      )
    )
    saveRDS(
      list(
        analysis_id = analysis_id,
        model_id = model_id,
        pooling_status = pooling_status,
        expected_imputations = mi_m,
        fitted_imputations = if (nrow(metrics)) sort(unique(metrics$imputation)) else integer(),
        valid_imputations = valid_imputations,
        model_files = fitted_files,
        mi_data_prefix = mi_object$prefix
      ),
      file.path(model_root, paste0(component_prefix, "__invalid_for_pool_manifest.rds")),
      compress = "xz"
    )
    return(pooling_status)
  }

  keep_valid <- function(value) {
    if (nrow(value) && "imputation" %in% names(value)) {
      value[imputation %in% valid_imputations]
    } else value
  }
  metrics <- keep_valid(metrics)
  q_values <- keep_valid(q_values)
  hr_values <- keep_valid(hr_values)
  p_values <- keep_valid(p_values)
  sojourn_values <- keep_valid(sojourn_values)
  prevalence_values <- keep_valid(prevalence_values)
  expected_time_values <- keep_valid(expected_time_values)
  component_confirmatory_blocks <- Filter(
    function(block) block$imputation %in% valid_imputations,
    component_confirmatory_blocks
  )

  pooled_hr <- rubin_pool_log_hr(hr_values, mi_m)
  if (!nrow(pooled_hr) || !"transition" %in% names(pooled_hr)) {
    record_failure(
      analysis_id, model_id, engine, "mi_pool_hazard_ratios", "NOT_ESTIMABLE",
      "No complete set of finite insomnia hazard-ratio estimates was available"
    )
    return("NOT_ESTIMABLE")
  }
  pooled_hr[, `:=`(model_id = model_id, engine = engine, cohort = config$cohort)]
  primary_edges <- c("S0->S1", "S0->S2", "S1->S3", "S2->S3")
  primary_rows <- pooled_hr[
    term == "insomnia_score_primary" & transition %chin% primary_edges
  ]
  if (
    nrow(primary_rows) != length(primary_edges) ||
      !setequal(primary_rows$transition, primary_edges) ||
      any(primary_rows$m_fitted != mi_m)
  ) {
    record_failure(
      analysis_id, model_id, engine, "mi_pool_hazard_ratios", "NOT_ESTIMABLE",
      "All four confirmatory insomnia effects were not available from all imputations"
    )
    return("NOT_ESTIMABLE")
  }
  support_sample <- sample_diagnostics(
    x, x, edges, analysis_id, model_id, requested
  )$edge_counts[, .(transition, risk_intervals, risk_participants, person_years, events)]
  pooled_hr <- merge(pooled_hr, support_sample, by = "transition", all.x = TRUE, sort = FALSE)
  pooled_hr[, `:=`(
    events_per_total_model_parameter = events / mean(metrics$n_parameters, na.rm = TRUE),
    events_per_parameter = events / mean(metrics$n_parameters, na.rm = TRUE),
    events_per_parameter_scope = paste(
      "Events divided by total fitted-model degrees of freedom;",
      "a conservative model-level support flag, not transition-specific EPP"
    )
  )]
  write_model_table(prefix, "hazard_ratios", pooled_hr)

  confirmatory_blocks_complete <- length(component_confirmatory_blocks) == mi_m
  confirmatory_labels_match <- FALSE
  if (confirmatory_blocks_complete) {
    reference_labels <- component_confirmatory_blocks[[1L]]$parameter_labels
    confirmatory_labels_match <- all(vapply(
      component_confirmatory_blocks,
      function(block) identical(block$parameter_labels, reference_labels),
      logical(1L)
    ))
  }
  if (confirmatory_blocks_complete && confirmatory_labels_match) {
    theta_matrix <- do.call(rbind, lapply(component_confirmatory_blocks, `[[`, "theta"))
    covariance_array <- simplify2array(lapply(component_confirmatory_blocks, `[[`, "covariance"))
    q_bar <- colMeans(theta_matrix)
    u_bar <- apply(covariance_array, c(1L, 2L), mean)
    between <- stats::cov(theta_matrix)
    m_block <- nrow(theta_matrix)
    total <- u_bar + (1 + 1 / m_block) * between
    template_block <- component_confirmatory_blocks[[1L]]
    pooled_block <- list(theta = q_bar, covariance = total)
    pooled_test <- wald_from_block(pooled_block)
    pooled_wald <- data.table(
      analysis_id = analysis_id,
      cohort = config$cohort,
      model_id = model_id,
      engine = engine,
      exposure_variable = "insomnia_score_primary",
      tested_terms = paste(template_block$exposure_terms, collapse = " | "),
      tested_transitions = paste(unique(template_block$transitions), collapse = " | "),
      statistic = pooled_test$statistic,
      df = pooled_test$df,
      p_value = pooled_test$p_value,
      method = paste(
        "Rubin total-covariance joint Wald with a large-sample chi-square reference;",
        "not a finite-m D1/F test"
      ),
      confirmatory_edge_count = uniqueN(template_block$transitions),
      m_expected = mi_m,
      m_fitted = m_block,
      pooling_status = if (m_block == mi_m) "RUBIN_POOLED_COMPLETE" else "PARTIAL_MI_NOT_PRIMARY"
    )
    write_model_table(prefix, "confirmatory_wald", pooled_wald)
  } else {
    record_failure(
      analysis_id, model_id, engine, "mi_pool_confirmatory_wald", "NOT_ESTIMABLE",
      paste(
        "Compatible confirmatory parameter blocks were unavailable:",
        length(component_confirmatory_blocks), "of", mi_m, "blocks;",
        "consistent labels/order =", confirmatory_labels_match
      )
    )
  }

  if (nrow(q_values)) {
    q_by <- intersect(c(
      "analysis_id", "transition", "from_state", "to_state", "from_index",
      "to_index", "death_transition", "covariate_profile"
    ), names(q_values))
    pooled_q <- rubin_pool_positive(q_values, "intensity", "standard_error", q_by, mi_m)
    if (nrow(pooled_q)) {
      setnames(pooled_q, "pooled_estimate", "intensity")
      pooled_q[, `:=`(model_id = model_id, engine = engine)]
      pooled_q <- merge(pooled_q, support_sample, by = "transition", all.x = TRUE, sort = FALSE)
      pooled_q[, `:=`(
        cohort = config$cohort,
        n_parameters = mean(metrics$n_parameters, na.rm = TRUE),
        events_per_total_model_parameter = events / mean(metrics$n_parameters, na.rm = TRUE),
        events_per_parameter = events / mean(metrics$n_parameters, na.rm = TRUE),
        events_per_parameter_scope = paste(
          "Events divided by total fitted-model degrees of freedom;",
          "a conservative model-level support flag, not transition-specific EPP"
        )
      )]
      write_model_table(prefix, "q_intensities", pooled_q)
    }
  }
  if (nrow(sojourn_values)) {
    sojourn_by <- intersect(c(
      "analysis_id", "state", "state_index", "state_label"
    ), names(sojourn_values))
    pooled_sojourn <- rubin_pool_positive(
      sojourn_values, "sojourn_years", "standard_error", sojourn_by, mi_m
    )
    if (nrow(pooled_sojourn)) {
      setnames(pooled_sojourn, "pooled_estimate", "sojourn_years")
      pooled_sojourn[, `:=`(model_id = model_id, engine = engine, cohort = config$cohort)]
      write_model_table(prefix, "sojourn_times", pooled_sojourn)
    }
  }
  if (nrow(p_values)) {
    p_by <- intersect(c(
      "analysis_id", "cohort", "exposure_variable", "profile", "time_years", "from_state",
      "to_state", "from_index", "to_index", "transition",
      "extrapolated_beyond_observed_followup", "nuisance_profile", "quantity_type",
      "direct_q_edge", "pathway_note"
    ), names(p_values))
    pooled_p <- p_values[, .(
      m_expected = mi_m,
      m_fitted = .N,
      pooling_status = if (.N == mi_m) "MI_MEAN_COMPLETE" else "PARTIAL_MI_NOT_PRIMARY",
      probability = mean(probability),
      between_imputation_sd = sd(probability),
      between_imputation_p2_5 = quantile(probability, 0.025),
      between_imputation_p97_5 = quantile(probability, 0.975),
      row_sum_error = max(row_sum_error)
    ), by = p_by]
    pooled_p[, `:=`(
      model_id = model_id,
      engine = engine,
      cohort = config$cohort,
      uncertainty_scope = paste(
        "Between-imputation distribution of point estimates only;",
        "p2.5/p97.5 are not 95% confidence limits"
      )
    )]
    write_model_table(prefix, "transition_probabilities", pooled_p)
  }
  if (nrow(prevalence_values)) {
    prevalence_by <- intersect(c(
      "analysis_id", "cohort", "time_years", "state", "state_index", "state_label",
      "diagnostic_method", "calibration_scope"
    ), names(prevalence_values))
    pooled_prevalence <- prevalence_values[, .(
      m_expected = mi_m,
      m_fitted = .N,
      pooling_status = if (.N == mi_m) "MI_MEAN_COMPLETE" else "PARTIAL_MI_NOT_PRIMARY",
      observed_count = mean(observed_count),
      expected_count = mean(expected_count),
      observed_probability = mean(observed_probability),
      expected_probability = mean(expected_probability),
      between_imputation_sd_expected_probability = sd(expected_probability),
      total_observed = mean(total_observed),
      total_expected = mean(total_expected)
    ), by = prevalence_by]
    pooled_prevalence[, `:=`(
      model_id = model_id,
      engine = engine,
      cohort = config$cohort,
      uncertainty_scope = paste(
        "Between-imputation variability summarized by SD only;",
        "no 95% confidence interval"
      )
    )]
    write_model_table(prefix, "observed_expected_prevalence", pooled_prevalence)
  }
  if (nrow(expected_time_values)) {
    time_by <- intersect(c(
      "analysis_id", "cohort", "exposure_variable", "profile", "start_state",
      "horizon_years", "state", "state_label", "summary_measure",
      "extrapolated_beyond_observed_followup", "nuisance_profile"
    ), names(expected_time_values))
    pooled_time <- expected_time_values[, .(
      m_expected = mi_m,
      m_fitted = .N,
      pooling_status = if (.N == mi_m) "MI_MEAN_COMPLETE" else "PARTIAL_MI_NOT_PRIMARY",
      expected_years = mean(expected_years),
      between_imputation_sd = sd(expected_years),
      between_imputation_p2_5 = quantile(expected_years, 0.025),
      between_imputation_p97_5 = quantile(expected_years, 0.975)
    ), by = time_by]
    pooled_time[, `:=`(
      model_id = model_id,
      engine = engine,
      cohort = config$cohort,
      uncertainty_scope = paste(
        "Between-imputation distribution of point estimates only;",
        "p2.5/p97.5 are not 95% confidence limits"
      )
    )]
    write_model_table(prefix, "expected_state_time", pooled_time)
  }
  pooled_metrics <- data.table(
    script_version = script_version,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    analysis_id = analysis_id,
    cohort = config$cohort,
    model_id = model_id,
    engine = engine,
    status = pooling_status,
    m_expected = mi_m,
    m_fitted = fitted_m,
    m_valid_for_primary_pooling = valid_m,
    n_records = mean(metrics$n_records, na.rm = TRUE),
    n_subjects = mean(metrics$n_subjects, na.rm = TRUE),
    n_parameters = mean(metrics$n_parameters, na.rm = TRUE),
    AIC_mean_not_pooled_inference = mean(metrics$AIC, na.rm = TRUE),
    BIC_subjects_mean_not_pooled_inference = mean(metrics$BIC_subjects, na.rm = TRUE),
    optimizer_all_zero = as.integer(all(metrics$optimizer_convergence_code == 0L)),
    standard_errors_all_found = as.integer(all(metrics$standard_errors_found == 1L)),
    q_row_sum_max_abs = max(metrics$q_row_sum_max_abs, na.rm = TRUE),
    death_row_max_abs = if (any(is.finite(metrics$death_row_max_abs))) {
      max(metrics$death_row_max_abs[is.finite(metrics$death_row_max_abs)])
    } else NA_real_,
    absorbing_row_max_abs = max(metrics$absorbing_row_max_abs, na.rm = TRUE),
    absorbing_state_ok = as.integer(all(metrics$absorbing_state_ok == 1L)),
    note = "AIC/BIC are summarized across imputations and are not Rubin-pooled inferential statistics"
  )
  write_model_table(prefix, "fit_metrics", pooled_metrics)
  saveRDS(
    list(
      analysis_id = analysis_id,
      model_id = model_id,
      pooling_status = pooling_status,
      expected_imputations = mi_m,
      fitted_imputations = fitted_m,
      model_files = fitted_files,
      mi_data_prefix = mi_object$prefix
    ),
    file.path(model_root, paste0(prefix, "__manifest.rds")),
    compress = "xz"
  )
  pooling_status
}

write_skip_metric <- function(config, model_id, engine, status, message) {
  prefix <- paste(config$analysis_id, model_id, sep = "__")
  metric <- data.table(
    script_version = script_version,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    analysis_id = config$analysis_id,
    cohort = config$cohort,
    model_id = model_id,
    engine = engine,
    status = status,
    message = message
  )
  write_model_table(prefix, "fit_metrics", metric)
  invisible(status)
}

aggregate_model_tables <- function(suffix, output_name) {
  files <- if (exists(suffix, envir = current_run_table_files, inherits = FALSE)) {
    get(suffix, envir = current_run_table_files, inherits = FALSE)
  } else character()
  files <- unique(files[file.exists(files)])
  value <- if (length(files)) {
    rbindlist(lapply(files, function(path) tryCatch(fread(path), error = function(e) NULL)), fill = TRUE)
  } else data.table()
  fwrite(value, file.path(result_root, output_name), na = "")
  invisible(value)
}

q_initialization_rows <- list()
analysis_qc_rows <- list()
analysis_summary_rows <- list()
models_attempted <- 0L

for (config in analysis_configs) {
  analysis_id <- config$analysis_id
  if (!grepl(analysis_regex, analysis_id)) {
    specifications[analysis_id == config$analysis_id, status := "FILTERED_ANALYSIS"]
    next
  }
  message("\n=== Preparing ", analysis_id, " ===")
  prepared <- tryCatch(prepare_analysis_data(config), error = function(e) e)
  if (inherits(prepared, "error")) {
    record_failure(
      analysis_id, "ALL", "data", "prepare_analysis", "FAILED",
      conditionMessage(prepared)
    )
    specifications[analysis_id == config$analysis_id, status := "FAILED_DATA_PREPARATION"]
    next
  }
  x <- prepared
  edges <- allowed_edges(config$n_states)
  qc <- analysis_qc(x, config, edges)
  analysis_qc_rows[[length(analysis_qc_rows) + 1L]] <- qc
  critical_failures <- qc[critical == 1L & violations > 0L, .N]
  analysis_summary_rows[[length(analysis_summary_rows) + 1L]] <- data.table(
    analysis_id = analysis_id,
    cohort = config$cohort,
    mortality_scope = config$mortality_scope,
    n_states = config$n_states,
    n_records = nrow(x),
    n_participants = uniqueN(x$id),
    n_deaths = x[state_msm == 5L, .N],
    minimum_time = min(x$time_msm),
    maximum_time = max(x$time_msm),
    critical_qc_failures = critical_failures
  )
  if (critical_failures > 0L) {
    record_failure(
      analysis_id, "ALL", "data", "analysis_qc", "FAILED",
      "Critical analysis-data QC failed; models were not attempted"
    )
    specifications[analysis_id == config$analysis_id, status := "FAILED_DATA_QC"]
    next
  }

  q_log <- file.path(log_root, paste0(analysis_id, "__q_initialization.log"))
  q_run <- safe_run(
    fun = function() GlobalAgingR::estimate_Q_matrix(
      data = as.data.frame(x[, .(id, time_msm, state_msm)]),
      method = "empirical",
      n_states = config$n_states,
      scale_factor = 0.1,
      smooth_factor = 0.01,
      validate = TRUE,
      verbose = TRUE,
      id_col = "id",
      time_col = "time_msm",
      state_col = "state_msm"
    ),
    log_file = q_log,
    analysis_id = analysis_id,
    model_id = "Q_INIT",
    engine = "GlobalAgingR::estimate_Q_matrix",
    stage = "q_initialization"
  )
  empirical_q <- q_run$value
  if (is.null(empirical_q) || !is.matrix(empirical_q) || any(dim(empirical_q) != config$n_states)) {
    empirical_q <- matrix(0.05, config$n_states, config$n_states)
    diag(empirical_q) <- -0.05 * (config$n_states - 1L)
    record_failure(
      analysis_id, "Q_INIT", "GlobalAgingR::estimate_Q_matrix",
      "q_initialization_fallback", "FALLBACK",
      "GlobalAgingR empirical Q unavailable; deterministic positive seeds used before structural masking",
      q_run$warnings, q_run$elapsed_seconds, q_log
    )
  }
  q_initial <- make_q_initial(empirical_q, edges, config$n_states)
  q_initialization_rows[[length(q_initialization_rows) + 1L]] <- matrix_long(
    empirical_q, analysis_id, "GlobalAgingR empirical unconstrained information", edges
  )
  q_initialization_rows[[length(q_initialization_rows) + 1L]] <- matrix_long(
    q_initial, analysis_id, "Protocol-constrained initial Q", edges
  )
  saveRDS(
    list(
      analysis_id = analysis_id,
      empirical_unconstrained = empirical_q,
      protocol_constrained = q_initial,
      allowed_edges = edges,
      main_q_has_direct_s0_s3 = FALSE,
      absorbing_death = config$n_states == 5L
    ),
    file.path(model_root, paste0(analysis_id, "__q_initialization.rds")),
    compress = "xz"
  )

  run_spec <- function(model_id, engine) {
    if (!grepl(model_regex, paste(analysis_id, model_id, sep = "__"))) return(FALSE)
    if (run_mode == "spec") return(FALSE)
    if (run_mode == "prototype" && !model_id %in% c("GA_M0", "M1")) return(FALSE)
    if (max_models > 0L && models_attempted >= max_models) return(FALSE)
    TRUE
  }

  # Required package-path smoke/baseline model.  It is intentionally minimal;
  # the formal models below expose constraints and diagnostics unavailable in
  # the GlobalAgingR wrapper.
  if (run_globalagingr_baseline && run_spec("GA_M0", "GlobalAgingR::msm_transition_model")) {
    models_attempted <- models_attempted + 1L
    requested <- c("insomnia_score_primary")
    sample_object <- tryCatch(
      prepare_model_sample(x, requested, exposure = "insomnia_score_primary"),
      error = function(e) e
    )
    if (inherits(sample_object, "error")) {
      status <- write_skip_metric(config, "GA_M0", "GlobalAgingR::msm_transition_model", "FAILED_PREFLIGHT", conditionMessage(sample_object))
      record_failure(analysis_id, "GA_M0", "GlobalAgingR::msm_transition_model", "preflight", status, conditionMessage(sample_object))
      specifications[spec_key == paste(analysis_id, "GA_M0", sep = "__"), status := ..status]
    } else {
      d <- sample_object$data
      sample_info <- sample_diagnostics(x, d, edges, analysis_id, "GA_M0", sample_object$used)
      prefix <- paste(analysis_id, "GA_M0", sep = "__")
      model_file <- file.path(model_root, paste0(prefix, ".rds"))
      ga_data_fingerprint <- fingerprint_table(
        d, c("id", "time_msm", "state_msm", requested)
      )
      spec_hash <- digest::digest(list(
        script_version, analysis_id, "GA_M0", requested, q_initial,
        ga_data_fingerprint
      ), algo = "xxhash64")
      envelope <- NULL
      ga_loaded_from_cache <- FALSE
      if (!force_refit && file.exists(model_file)) {
        cached <- tryCatch(readRDS(model_file), error = function(e) NULL)
        if (!is.null(cached) && identical(cached$spec_hash, spec_hash)) {
          envelope <- cached
          ga_loaded_from_cache <- TRUE
        }
      }
      if (is.null(envelope)) {
        log_file <- file.path(log_root, paste0(prefix, ".log"))
        run <- safe_run(
          fun = function() GlobalAgingR::msm_transition_model(
            data = as.data.frame(d),
            state_var = "state_msm",
            time_var = "time_msm",
            id_var = "id",
            covariates = requested,
            qmatrix = q_initial,
            death_state = if (config$n_states == 5L) 5L else NULL
          ),
          log_file = log_file,
          analysis_id = analysis_id,
          model_id = "GA_M0",
          engine = "GlobalAgingR::msm_transition_model",
          stage = "fit"
        )
        if (is.null(run$value)) {
          message_text <- if (is.null(run$error)) "GlobalAgingR wrapper returned NULL" else run$error
          status <- write_skip_metric(config, "GA_M0", "GlobalAgingR::msm_transition_model", "FAILED", message_text)
          if (is.null(run$error)) record_failure(
            analysis_id, "GA_M0", "GlobalAgingR::msm_transition_model", "fit", "FAILED",
            message_text, run$warnings, run$elapsed_seconds, log_file
          )
          specifications[spec_key == prefix, status := ..status]
        } else {
          envelope <- list(
            spec_hash = spec_hash,
            engine = "GlobalAgingR::msm_transition_model",
            fit = run$value,
            model_data_fingerprint = ga_data_fingerprint,
            sample_ids = unique(d$id),
            warnings = run$warnings,
            elapsed_seconds = run$elapsed_seconds,
            log_file = log_file
          )
          saveRDS(envelope, model_file, compress = "xz")
        }
      }
      if (!is.null(envelope)) {
        fit <- envelope$fit$model
        status <- extract_and_write(
          fit, analysis_id, "GA_M0", "GlobalAgingR::msm_transition_model",
          edges, d, sample_info, envelope$elapsed_seconds, envelope$warnings,
          requested, sample_object$used,
          c(sample_object$missing_columns, sample_object$constant_columns),
          envelope$log_file,
          exposure = "insomnia_score_primary"
        )
        specifications[spec_key == prefix, status := if (ga_loaded_from_cache) {
          paste0(..status, "_CACHED")
        } else ..status]
      }
    }
  } else {
    specifications[spec_key == paste(analysis_id, "GA_M0", sep = "__"), status := if (!run_globalagingr_baseline) "DISABLED" else "FILTERED_OR_DEFERRED"]
  }

  for (model_id in names(model_definitions)) {
    prefix <- paste(analysis_id, model_id, sep = "__")
    engine <- "msm::msm"
    if (!run_spec(model_id, engine)) {
      specifications[spec_key == prefix, status := "FILTERED_OR_DEFERRED"]
      next
    }
    models_attempted <- models_attempted + 1L
    requested <- model_definitions[[model_id]]
    exposure <- unname(model_exposures[model_id])
    sample_object <- tryCatch(
      prepare_model_sample(x, requested, exposure = exposure),
      error = function(e) e
    )
    if (inherits(sample_object, "error")) {
      status <- write_skip_metric(config, model_id, engine, "FAILED_PREFLIGHT", conditionMessage(sample_object))
      record_failure(analysis_id, model_id, engine, "preflight", status, conditionMessage(sample_object))
      specifications[spec_key == prefix, status := ..status]
      next
    }
    d <- sample_object$data
    if (nrow(d) < 2L || uniqueN(d$id) < 2L) {
      status <- write_skip_metric(config, model_id, engine, "SKIPPED_TOO_SMALL", "Fewer than two participants or records after complete-case filtering")
      specifications[spec_key == prefix, status := ..status]
      next
    }
    if (max_panel_rows > 0L && nrow(d) > max_panel_rows) {
      message_text <- paste0(
        "Complete-case panel has ", nrow(d), " rows, exceeding GA_MARKOV_MAX_PANEL_ROWS=",
        max_panel_rows, ". No subsampling was performed."
      )
      status <- write_skip_metric(config, model_id, engine, "SKIPPED_ROW_LIMIT", message_text)
      record_failure(analysis_id, model_id, engine, "resource_guard", status, message_text)
      specifications[spec_key == prefix, status := ..status]
      next
    }
    covariates <- sample_object$used
    cov_formula <- reformulate(covariates)
    constraints <- make_nuisance_constraints(
      cov_formula, d, edges, exposure = exposure
    )
    sample_info <- sample_diagnostics(x, d, edges, analysis_id, model_id, covariates)
    write_model_table(prefix, "sample_summary", sample_info$summary)
    write_model_table(prefix, "edge_counts", sample_info$edge_counts)

    model_data_fingerprint <- fingerprint_table(
      d, unique(c("id", "time_msm", "state_msm", covariates))
    )
    spec_hash <- digest::digest(list(
      script_version, analysis_id, model_id, covariates, q_initial,
      nuisance_constraint_mode, config$death_exact_state,
      max_iterations, model_data_fingerprint
    ), algo = "xxhash64")
    model_file <- file.path(model_root, paste0(prefix, ".rds"))
    envelope <- NULL
    loaded_from_cache <- FALSE
    if (!force_refit && file.exists(model_file)) {
      cached <- tryCatch(readRDS(model_file), error = function(e) NULL)
      if (!is.null(cached) && identical(cached$spec_hash, spec_hash)) {
        envelope <- cached
        loaded_from_cache <- TRUE
      }
    }
    if (is.null(envelope)) {
      log_file <- file.path(log_root, paste0(prefix, ".log"))
      run <- safe_run(
        fun = function() msm::msm(
          formula = state_msm ~ time_msm,
          subject = id,
          data = as.data.frame(d),
          qmatrix = q_initial,
          covariates = cov_formula,
          constraint = constraints,
          deathexact = config$death_exact_state,
          exacttimes = FALSE,
          center = TRUE,
          hessian = TRUE,
          control = list(fnscale = 10000, maxit = max_iterations, trace = 0)
        ),
        log_file = log_file,
        analysis_id = analysis_id,
        model_id = model_id,
        engine = engine,
        stage = "fit"
      )
      if (is.null(run$value)) {
        message_text <- if (is.null(run$error)) "msm::msm returned NULL" else run$error
        status <- write_skip_metric(config, model_id, engine, "FAILED", message_text)
        if (is.null(run$error)) record_failure(
          analysis_id, model_id, engine, "fit", "FAILED", message_text,
          run$warnings, run$elapsed_seconds, log_file
        )
        specifications[spec_key == prefix, status := ..status]
        next
      }
      envelope <- list(
        spec_hash = spec_hash,
        engine = engine,
        fit = run$value,
        requested_covariates = requested,
        used_covariates = covariates,
        omitted_covariates = c(sample_object$missing_columns, sample_object$constant_columns),
        constraints = constraints,
        model_data_fingerprint = model_data_fingerprint,
        sample_ids = unique(d$id),
        warnings = run$warnings,
        elapsed_seconds = run$elapsed_seconds,
        log_file = log_file
      )
      saveRDS(envelope, model_file, compress = "xz")
    }
    status <- extract_and_write(
      envelope$fit, analysis_id, model_id, engine, edges, d, sample_info,
      envelope$elapsed_seconds, envelope$warnings, requested, covariates,
      c(sample_object$missing_columns, sample_object$constant_columns),
      envelope$log_file,
      exposure = exposure
    )
    specifications[spec_key == prefix, status := if (loaded_from_cache) paste0(..status, "_CACHED") else ..status]
  }

  # The protocol-defined primary M4 analysis uses 20 imputations of baseline
  # adjustment covariates. Exposure, state, time, cohort and sampling weight
  # are protected from imputation. Each imputed data set is fitted separately
  # and the transition-specific insomnia effects are Rubin pooled.
  mi_key <- paste(analysis_id, "M4_MI20", sep = "__")
  mi_selected <- grepl(model_regex, mi_key)
  mi_resource_available <- max_models <= 0L || models_attempted < max_models
  if (mi_mode == "off") {
    specifications[spec_key == mi_key, status := "DISABLED"]
  } else if (!mi_selected) {
    specifications[spec_key == mi_key, status := "FILTERED_OR_DEFERRED"]
  } else if (run_mode == "spec") {
    specifications[spec_key == mi_key, status := "SPECIFICATION_ONLY"]
  } else if (run_mode == "prototype") {
    specifications[spec_key == mi_key, status := "DEFERRED_IN_PROTOTYPE"]
  } else if (!mi_resource_available) {
    specifications[spec_key == mi_key, status := "DEFERRED_BY_MODEL_LIMIT"]
  } else {
    message("--- Preparing primary multiple imputation for ", analysis_id, " ---")
    mi_object <- tryCatch(prepare_mi_baseline(copy(x), config), error = function(e) e)
    if (inherits(mi_object, "error")) {
      specifications[spec_key == mi_key, status := "FAILED_MI_PREPARATION"]
      record_failure(
        analysis_id, "M4_MI20", "mice::mice + msm::msm + Rubin pooling",
        "mi_prepare", "FAILED", conditionMessage(mi_object)
      )
    } else if (mi_mode == "prepare") {
      specifications[spec_key == mi_key, status := "MI_DATA_PREPARED_NOT_FITTED"]
    } else {
      models_attempted <- models_attempted + 1L
      mi_status <- fit_mi_primary(mi_object, x, config, edges, q_initial)
      specifications[spec_key == mi_key, status := mi_status]
    }
  }
  rm(x)
  gc(verbose = FALSE)
}

initialization_table <- if (length(q_initialization_rows)) {
  rbindlist(q_initialization_rows, use.names = TRUE, fill = TRUE)
} else data.table()
qc_table <- if (length(analysis_qc_rows)) {
  rbindlist(analysis_qc_rows, use.names = TRUE, fill = TRUE)
} else data.table()
analysis_summary <- if (length(analysis_summary_rows)) {
  rbindlist(analysis_summary_rows, use.names = TRUE, fill = TRUE)
} else data.table()

fwrite(initialization_table, file.path(result_root, "q_initialization_matrices.csv"), na = "")
fwrite(qc_table, file.path(qc_root, "markov_analysis_data_qc.csv"), na = "")
fwrite(analysis_summary, file.path(result_root, "markov_analysis_samples.csv"), na = "")
fwrite(specifications, file.path(result_root, "markov_model_specifications.csv"), na = "")

aggregate_model_tables("fit_metrics", "markov_fit_metrics.csv")
aggregate_model_tables("fit_metrics", "markov_model_diagnostics.csv")
aggregate_model_tables("sample_summary", "markov_model_samples.csv")
aggregate_model_tables("edge_counts", "markov_model_edge_counts.csv")
aggregate_model_tables("q_intensities", "markov_q_intensities.csv")
wald_table <- aggregate_model_tables("confirmatory_wald", "markov_confirmatory_joint_wald.csv")
hazard_table <- aggregate_model_tables("hazard_ratios", "markov_hazard_ratios.csv")
if (nrow(hazard_table)) {
  confirmatory_edges <- c("S0->S1", "S0->S2", "S1->S3", "S2->S3")
  if (!"pooling_status" %in% names(hazard_table)) hazard_table[, pooling_status := NA_character_]
  # Complete-case/individual-imputation tables carry this flag directly from
  # extract_hazard_ratios().  Rubin-pooled tables are reconstructed from the
  # requested exposure rows and therefore need the same flag restored before
  # the confirmatory-family annotation below.
  if (!"is_requested_insomnia_exposure" %in% names(hazard_table)) {
    hazard_table[, is_requested_insomnia_exposure := as.integer(
      !is.na(term) & !is.na(exposure_variable) &
        (term == exposure_variable | startsWith(term, exposure_variable))
    )]
  }
  if (!"events_per_total_model_parameter" %in% names(hazard_table)) {
    hazard_table[, events_per_total_model_parameter := events_per_parameter]
  }
  hazard_table[, `:=`(
    confirmatory_nonfatal = as.integer(transition %chin% confirmatory_edges),
    primary_confirmatory_family = as.integer(
      model_id == "M4_MI20" & term == "insomnia_score_primary" &
        is_requested_insomnia_exposure == 1L & transition %chin% confirmatory_edges &
        !is.na(pooling_status) & pooling_status == "RUBIN_POOLED_COMPLETE"
    ),
    event_parameter_risk = fifelse(
      is.na(events_per_total_model_parameter), "not assessed",
      fifelse(events_per_total_model_parameter < 5, "sparse: events/total-model-df <5",
        fifelse(events_per_total_model_parameter < 10,
          "caution: events/total-model-df 5-<10",
          "adequate: events/total-model-df >=10"))
    ),
    p_fdr_bh_confirmatory = NA_real_,
    fdr_family_size = NA_integer_,
    joint_wald_p_value = NA_real_,
    joint_wald_gate_passed = NA_integer_
  )]
  if (nrow(wald_table) && all(c("analysis_id", "model_id", "p_value") %in% names(wald_table))) {
    if (!"pooling_status" %in% names(wald_table)) wald_table[, pooling_status := NA_character_]
    gate <- unique(wald_table[
      model_id == "M4_MI20" & pooling_status == "RUBIN_POOLED_COMPLETE",
      .(
        analysis_id,
        model_id,
        joint_wald_p_value_from_test = p_value,
        joint_wald_gate_passed_from_test = as.integer(is.finite(p_value) & p_value < 0.05)
      )
    ])
    if (nrow(gate)) {
      hazard_table <- merge(hazard_table, gate,
        by = c("analysis_id", "model_id"), all.x = TRUE, sort = FALSE
      )
      hazard_table[!is.na(joint_wald_p_value_from_test), `:=`(
        joint_wald_p_value = joint_wald_p_value_from_test,
        joint_wald_gate_passed = joint_wald_gate_passed_from_test
      )]
      hazard_table[, c(
        "joint_wald_p_value_from_test", "joint_wald_gate_passed_from_test"
      ) := NULL]
    }
  }
  hazard_table[
    primary_confirmatory_family == 1L & joint_wald_gate_passed == 1L,
    `:=`(
      p_fdr_bh_confirmatory = p.adjust(p_value, method = "BH"),
      fdr_family_size = sum(is.finite(p_value))
    ),
    by = .(analysis_id, cohort, model_id, term)
  ]
  hazard_table[, `:=`(
    fdr_gate_status = fifelse(
      primary_confirmatory_family != 1L,
      "Not in the primary M4_MI20 four-edge family",
      fifelse(
        is.na(joint_wald_gate_passed),
        "Primary joint Wald gate unavailable; no edge-level confirmatory claim",
        fifelse(
          joint_wald_gate_passed == 0L,
          "Primary joint Wald gate not passed; no edge-level confirmatory claim",
          "Primary joint Wald gate passed; interpret four-edge BH-FDR"
        )
      )
    ),
    edge_level_confirmatory_claim = as.integer(
      primary_confirmatory_family == 1L & joint_wald_gate_passed == 1L &
        fdr_family_size == 4L &
        is.finite(p_fdr_bh_confirmatory) & p_fdr_bh_confirmatory < 0.05
    )
  )]
  fwrite(hazard_table, file.path(result_root, "markov_hazard_ratios.csv"), na = "")
}
aggregate_model_tables("transition_probabilities", "markov_transition_probabilities_2_5_10y.csv")
aggregate_model_tables("sojourn_times", "markov_sojourn_times.csv")
aggregate_model_tables("expected_state_time", "markov_expected_state_time_10y.csv")
aggregate_model_tables("observed_expected_prevalence", "markov_observed_expected_prevalence.csv")

mi_audit_files <- unique(paste0(current_mi_prefixes, "__audit.csv"))
mi_audit_files <- mi_audit_files[file.exists(mi_audit_files)]
mi_audit <- if (length(mi_audit_files)) {
  rbindlist(lapply(mi_audit_files, function(path) tryCatch(fread(path), error = function(e) NULL)), fill = TRUE)
} else data.table()
fwrite(mi_audit, file.path(result_root, "multiple_imputation_audit.csv"), na = "")

mi_event_files <- unique(paste0(current_mi_prefixes, "__logged_events.csv"))
mi_event_files <- mi_event_files[file.exists(mi_event_files)]
mi_events <- if (length(mi_event_files)) {
  rbindlist(lapply(mi_event_files, function(path) tryCatch(fread(path), error = function(e) NULL)), fill = TRUE)
} else data.table()
fwrite(mi_events, file.path(result_root, "multiple_imputation_logged_events.csv"), na = "")
failures <- flush_failures()

writeLines(
  c(
    paste("script_version:", script_version),
    paste("run_timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste("run_mode:", run_mode),
    paste("model_regex:", model_regex),
    paste("analysis_regex:", analysis_regex),
    paste("force_refit:", force_refit),
    paste("nuisance_constraint_mode:", nuisance_constraint_mode),
    paste("mi_mode:", mi_mode),
    paste("mi_range:", paste0(mi_start, "-", mi_end)),
    paste("mi_pool:", mi_pool),
    paste("models_attempted:", models_attempted),
    paste("failures_or_partial_failures_this_log:", nrow(failures)),
    "",
    capture.output(sessionInfo())
  ),
  con = file.path(log_root, "markov_session_info.txt"),
  useBytes = TRUE
)

cat("Markov model workflow finished without fabricating skipped or failed estimates.\n")
cat("Run mode:", run_mode, "\n")
cat("Specifications:", file.path(result_root, "markov_model_specifications.csv"), "\n")
cat("Fit metrics:", file.path(result_root, "markov_fit_metrics.csv"), "\n")
cat("Failure log:", file.path(log_root, "markov_failures.csv"), "\n")
if (nrow(qc_table) && qc_table[critical == 1L & violations > 0L, .N] > 0L) {
  quit(save = "no", status = 2L)
}
