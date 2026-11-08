#!/usr/bin/env Rscript

# Collate restartable Markov shards, run-labelled Cox outputs and protocol
# sensitivity analyses into one audited result bundle.
#
# This script is intentionally estimate-conservative: it only copies or combines
# numerical estimates that already exist upstream.  When a prespecified output
# is absent, the canonical table receives an explicitly non-completed audit row
# with missing estimate fields and a reason; no numerical result is invented.

suppressPackageStartupMessages(library(data.table))

options(stringsAsFactors = FALSE, width = 180)

script_version <- "2026-09-08-v3"
stage_root <- Sys.getenv("GA_STAGE_ROOT", unset = "private/stage")
output_root <- Sys.getenv(
  "GA_MULTISTATE_OUTPUT_ROOT",
  unset = file.path(stage_root, "multistate_output")
)
result_root <- file.path(output_root, "results")
report_root <- file.path(output_root, "reports")
log_root <- file.path(output_root, "logs")
qc_root <- file.path(output_root, "qc")

for (path in c(result_root, report_root)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
artifact_sources <- list()
written_artifacts <- character()

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

normalise_name <- function(x) {
  gsub("(^_+|_+$)", "", gsub("[^a-z0-9]+", "_", tolower(clean_text(x))))
}

truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & x != 0)
  normalise_name(x) %chin% c("1", "true", "yes", "y", "selected", "prespecified")
}

first_existing_column <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit)) hit[1L] else NA_character_
}

cohort_from_analysis <- function(x) {
  z <- tolower(clean_text(x))
  fifelse(startsWith(z, "hrs"), "HRS", fifelse(startsWith(z, "elsa"), "ELSA", NA_character_))
}

mortality_scope_from_analysis <- function(analysis_id, cohort = cohort_from_analysis(analysis_id)) {
  analysis <- tolower(clean_text(analysis_id))
  cohort <- toupper(clean_text(cohort))
  fcase(
    grepl("elsa_5state|through_2012|death.*2012", analysis),
    "ELSA limited mortality validation through 2012; death ascertainment ends in 2012",
    grepl("elsa_4state|w4_w10", analysis),
    "ELSA W4-W10 nonfatal follow-up; death state excluded from this analysis",
    cohort == "ELSA",
    "ELSA nonfatal follow-up; mortality analyses are limited through 2012",
    cohort == "HRS",
    "HRS mortality follow-up through 2020; death timing follows the R/05 definition",
    default = "Scope not recoverable from upstream identifier"
  )
}

terminal_status <- function(x) {
  z <- normalise_name(x)
  fcase(
    grepl("^(fitted|success|completed|complete|done)$", z) |
      grepl("^fitted_with_inference_warning$|rubin_pooled_complete|mi_mean_complete", z),
    "completed",
    grepl("not_estim|not_run|skip|disabled|planned|prepare|partial_mi_not_primary", z),
    "not_estimable",
    grepl("fail|error|partial_failure|timeout|interrupt", z),
    "failed",
    z %chin% c("completed", "not_estimable", "failed"), z,
    default = NA_character_
  )
}

status_rank <- function(x) {
  fcase(x == "completed", 3L, x == "not_estimable", 2L, x == "failed", 1L, default = 0L)
}

atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path)), tmpdir = dirname(path), fileext = ".tmp")
  on.exit(unlink(tmp, force = TRUE), add = TRUE)
  fwrite(as.data.table(x), tmp, na = "")
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not write canonical output: ", path)
  written_artifacts <<- unique(c(written_artifacts, normalizePath(path, winslash = "/", mustWork = FALSE)))
  invisible(path)
}

relative_to_output <- function(path) {
  root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
  value <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(tolower(root), "/")
  if (startsWith(tolower(value), prefix)) substring(value, nchar(root) + 2L) else value
}

existing_files <- function(paths) {
  paths <- unique(normalizePath(paths, winslash = "/", mustWork = FALSE))
  paths[file.exists(paths) & !dir.exists(paths)]
}

source_files_from_rows <- function(x) {
  if (!nrow(x) || !"__source_file" %in% names(x)) return(character())
  relative_paths <- unique(clean_text(x[["__source_file"]]))
  relative_paths <- relative_paths[nzchar(relative_paths)]
  if (!length(relative_paths)) return(character())
  existing_files(file.path(output_root, relative_paths))
}

parse_completion_times <- function(x) {
  if (inherits(x, "POSIXt")) return(as.numeric(x))
  values <- clean_text(x)
  values <- sub("Z$", "+0000", values)
  values <- sub("([+-][0-9]{2}):([0-9]{2})$", "\\1\\2", values)
  parsed <- rep(NA_real_, length(values))
  for (format in c("%Y-%m-%dT%H:%M:%S%z", "%Y-%m-%d %H:%M:%S%z", "%Y-%m-%d %H:%M:%S")) {
    missing <- !is.finite(parsed) & nzchar(values)
    if (!any(missing)) break
    parsed[missing] <- suppressWarnings(as.numeric(strptime(values[missing], format = format, tz = "UTC")))
  }
  parsed
}

safe_read_many <- function(paths, source_priority = NULL) {
  paths <- existing_files(paths)
  rows <- list()
  used <- character()
  if (!length(paths)) {
    value <- data.table()
    attr(value, "source_files") <- character()
    return(value)
  }
  for (path in paths) {
    info <- file.info(path)
    if (!is.finite(info$size) || info$size <= 0L) next
    value <- tryCatch(
      suppressWarnings(fread(path, showProgress = FALSE)),
      error = function(e) NULL
    )
    if (is.null(value) || !nrow(value) || !ncol(value)) next
    priority <- if (is.null(source_priority)) {
      if (grepl("/by_model/", gsub("\\\\", "/", path), fixed = TRUE)) 1L else 2L
    } else if (length(source_priority) == 1L) {
      as.integer(source_priority)
    } else {
      as.integer(source_priority[match(path, paths)])
    }
    value[, `:=`(
      `__source_file` = relative_to_output(path),
      `__source_priority` = priority,
      `__source_mtime` = as.numeric(info$mtime)
    )]
    rows[[length(rows) + 1L]] <- value
    used <- c(used, path)
  }
  value <- if (length(rows)) rbindlist(rows, use.names = TRUE, fill = TRUE) else data.table()
  attr(value, "source_files") <- unique(used)
  value
}

remember_sources <- function(artifact_name, values = list(), paths = character()) {
  sources <- unique(c(
    paths,
    unlist(lapply(values, function(value) attr(value, "source_files")), use.names = FALSE)
  ))
  artifact_sources[[artifact_name]] <<- existing_files(sources)
  invisible(sources)
}

nonmissing_score <- function(x, columns) {
  columns <- intersect(columns, names(x))
  if (!length(columns)) return(rep(0L, nrow(x)))
  rowSums(as.data.frame(lapply(x[, ..columns], function(z) {
    !is.na(z) & (if (is.character(z) || is.factor(z)) nzchar(clean_text(z)) else TRUE)
  })))
}

deduplicate_rows <- function(x, keys, value_columns = character(), prefer_latest_source = FALSE) {
  if (!nrow(x)) return(x)
  keys <- intersect(keys, names(x))
  if (!length(keys)) keys <- "__row_key"
  if (identical(keys, "__row_key")) x[, `__row_key` := seq_len(.N)]
  if (!"__source_priority" %in% names(x)) x[, `__source_priority` := 9L]
  if (!"__source_mtime" %in% names(x)) x[, `__source_mtime` := 0]
  if (!"status" %in% names(x)) x[, status := NA_character_]
  x[, `__status_rank` := status_rank(status)]
  x[, `__value_rank` := nonmissing_score(.SD, value_columns)]
  if (isTRUE(prefer_latest_source)) {
    order_columns <- c(keys, "__source_mtime", "__status_rank", "__value_rank", "__source_priority")
    order_direction <- c(rep(1L, length(keys)), -1L, -1L, -1L, 1L)
  } else {
    order_columns <- c(keys, "__status_rank", "__value_rank", "__source_priority", "__source_mtime")
    order_direction <- c(rep(1L, length(keys)), -1L, -1L, 1L, -1L)
  }
  setorderv(x, order_columns, order_direction)
  value <- x[, .SD[1L], by = keys]
  if ("__row_key" %in% names(value)) value[, `__row_key` := NULL]
  value[]
}

strip_internal_columns <- function(x) {
  internal <- grep("^__", names(x), value = TRUE)
  if (length(internal)) x[, (internal) := NULL]
  x[]
}

ensure_context <- function(x) {
  if (!nrow(x)) return(x)
  if (!"analysis_id" %in% names(x)) x[, analysis_id := NA_character_]
  if (!"cohort" %in% names(x)) x[, cohort := cohort_from_analysis(analysis_id)]
  x[is.na(cohort) | !nzchar(clean_text(cohort)), cohort := cohort_from_analysis(analysis_id)]
  x[, cohort := toupper(clean_text(cohort))]
  if (!"model_id" %in% names(x)) x[, model_id := NA_character_]
  if (!"mortality_scope" %in% names(x)) {
    x[, mortality_scope := mortality_scope_from_analysis(analysis_id, cohort)]
  } else {
    x[is.na(mortality_scope) | !nzchar(clean_text(mortality_scope)),
      mortality_scope := mortality_scope_from_analysis(analysis_id, cohort)]
  }
  x[cohort == "ELSA" & grepl("5state|through_2012", tolower(analysis_id)),
    mortality_scope := "ELSA limited mortality validation through 2012; death ascertainment ends in 2012"]
  x[]
}

markov_spec_id <- function(analysis_id, model_id) {
  paste("MKV", clean_text(analysis_id), clean_text(model_id), sep = "::")
}

markov_roots <- c(file.path(result_root, "markov"))
shard_parent <- file.path(result_root, "markov_shards")
if (dir.exists(shard_parent)) {
  markov_roots <- c(markov_roots, list.dirs(shard_parent, recursive = FALSE, full.names = TRUE))
}
markov_roots <- unique(markov_roots[dir.exists(markov_roots)])
markov_by_model_roots <- file.path(markov_roots, "by_model")
markov_by_model_roots <- markov_by_model_roots[dir.exists(markov_by_model_roots)]

find_markov_table_files <- function(suffix, aggregate_names = character()) {
  by_model <- unlist(lapply(markov_by_model_roots, function(path) {
    list.files(path, pattern = paste0("_", suffix, "[.]csv$"), full.names = TRUE)
  }), use.names = FALSE)
  if (is.null(by_model)) by_model <- character()
  # Component-level MI estimates are inputs to Rubin pooling, not final estimates.
  by_model <- by_model[!grepl("_mi_component_", basename(by_model), fixed = TRUE)]
  aggregate <- unlist(lapply(markov_roots, function(path) {
    file.path(path, aggregate_names)
  }), use.names = FALSE)
  if (is.null(aggregate)) aggregate <- character()
  existing_files(c(by_model, aggregate))
}

collect_markov <- function(suffix, aggregate_names, keys, values) {
  files <- find_markov_table_files(suffix, aggregate_names)
  x <- safe_read_many(files)
  sources <- attr(x, "source_files")
  if (nrow(x)) x <- deduplicate_rows(x, keys, values)
  attr(x, "source_files") <- sources
  x
}

# -----------------------------------------------------------------------------
# Markov specifications, terminal fit states and failure evidence
# -----------------------------------------------------------------------------

markov_spec_files <- existing_files(file.path(markov_roots, "markov_model_specifications.csv"))
markov_specs <- safe_read_many(markov_spec_files)
markov_spec_sources <- attr(markov_specs, "source_files")
if (nrow(markov_specs)) {
  markov_specs <- ensure_context(markov_specs)
  markov_specs[, spec_id := markov_spec_id(analysis_id, model_id)]
  markov_specs[, prespecified := TRUE]
  if ("status" %in% names(markov_specs)) setnames(markov_specs, "status", "upstream_spec_status")
  markov_specs <- deduplicate_rows(
    markov_specs,
    c("cohort", "spec_id"),
    c("model_label", "exposure", "covariates", "missing_data", "mortality_scope")
  )
}

markov_fit <- collect_markov(
  "fit_metrics", "markov_fit_metrics.csv",
  c("analysis_id", "model_id"),
  c("log_likelihood", "AIC", "q_row_sum_max_abs", "death_row_max_abs", "n_subjects")
)
markov_fit_sources <- attr(markov_fit, "source_files")
if (nrow(markov_fit)) {
  markov_fit <- ensure_context(markov_fit)
  if (!"status" %in% names(markov_fit)) markov_fit[, status := NA_character_]
  markov_fit[, upstream_status := clean_text(status)]
  markov_fit[, status := terminal_status(upstream_status)]
  markov_fit[is.na(status), status := "failed"]
  reason_column <- first_existing_column(markov_fit, c("reason", "message", "optimizer_message", "note"))
  if (is.na(reason_column)) {
    markov_fit[, reason := ""]
  } else if (reason_column != "reason") {
    markov_fit[, reason := clean_text(get(reason_column))]
  } else {
    markov_fit[, reason := clean_text(reason)]
  }
  markov_fit[status == "completed", reason := ""]
  markov_fit[status != "completed" & !nzchar(reason),
    reason := paste0("Upstream Markov terminal status: ", upstream_status)]
  markov_fit[, spec_id := markov_spec_id(analysis_id, model_id)]
  markov_fit <- deduplicate_rows(
    markov_fit,
    c("analysis_id", "model_id"),
    c("log_likelihood", "AIC", "q_row_sum_max_abs", "death_row_max_abs", "n_subjects")
  )
}

markov_failure_files <- character()
for (candidate_root in c(
  file.path(log_root, "markov"),
  if (dir.exists(file.path(log_root, "markov_shards"))) {
    list.dirs(file.path(log_root, "markov_shards"), recursive = FALSE, full.names = TRUE)
  } else character()
)) {
  markov_failure_files <- c(markov_failure_files, file.path(candidate_root, "markov_failures.csv"))
}
markov_failures <- safe_read_many(markov_failure_files, source_priority = 3L)
markov_failure_sources <- attr(markov_failures, "source_files")
if (nrow(markov_failures)) {
  if (!"analysis_id" %in% names(markov_failures)) markov_failures[, analysis_id := NA_character_]
  if (!"model_id" %in% names(markov_failures)) markov_failures[, model_id := NA_character_]
  if (!"message" %in% names(markov_failures)) markov_failures[, message := ""]
  if (!"status" %in% names(markov_failures)) markov_failures[, status := "failed"]
  markov_failures[, failure_terminal_status := terminal_status(status)]
  markov_failures[is.na(failure_terminal_status), failure_terminal_status := "failed"]
  markov_failures[, failure_reason := clean_text(message)]
  markov_failures[!nzchar(failure_reason), failure_reason := paste0("Upstream Markov failure status: ", status)]
  markov_failures <- deduplicate_rows(
    markov_failures,
    c("analysis_id", "model_id"),
    c("failure_reason", "stage", "timestamp")
  )
}

if (!nrow(markov_specs) && nrow(markov_fit)) {
  markov_specs <- unique(markov_fit[, .(
    spec_id, analysis_id, cohort, model_id, engine, mortality_scope,
    model_label = model_id, exposure = NA_character_, covariates = NA_character_,
    missing_data = NA_character_, upstream_spec_status = "RECOVERED_FROM_FIT",
    prespecified = TRUE
  )], by = c("cohort", "spec_id"))
}

markov_registry <- copy(markov_specs)
if (nrow(markov_registry)) {
  fit_state <- if (nrow(markov_fit)) unique(markov_fit[, .(
    analysis_id, model_id,
    fit_status = status,
    fit_reason = reason,
    fit_upstream_status = upstream_status
  )], by = c("analysis_id", "model_id")) else data.table(
    analysis_id = character(), model_id = character(), fit_status = character(),
    fit_reason = character(), fit_upstream_status = character()
  )
  markov_registry <- merge(markov_registry, fit_state, by = c("analysis_id", "model_id"), all.x = TRUE, sort = FALSE)
  if (nrow(markov_failures)) {
    failure_state <- unique(markov_failures[, .(
      analysis_id, model_id, failure_terminal_status, failure_reason
    )], by = c("analysis_id", "model_id"))
    markov_registry <- merge(markov_registry, failure_state, by = c("analysis_id", "model_id"), all.x = TRUE, sort = FALSE)
  } else {
    markov_registry[, `:=`(failure_terminal_status = NA_character_, failure_reason = NA_character_)]
  }
  markov_registry[, status := fcase(
    !is.na(fit_status), fit_status,
    !is.na(failure_terminal_status), failure_terminal_status,
    default = "not_estimable"
  )]
  markov_registry[, reason := fcase(
    status == "completed", "",
    !is.na(fit_reason) & nzchar(clean_text(fit_reason)), clean_text(fit_reason),
    !is.na(failure_reason) & nzchar(clean_text(failure_reason)), clean_text(failure_reason),
    grepl("disabled", normalise_name(upstream_spec_status)),
    "Prespecified Markov model was disabled in the available run configuration; no pooled estimate is present",
    default = "No completed fit or terminal fit record was present at collation time; the prespecified model was not run or is unavailable"
  )]
  markov_registry[, module := "markov"]
  if (!"analysis_family" %in% names(markov_registry)) {
    markov_registry[, analysis_family := fcase(
      model_id %chin% c("M1", "M2", "M3", "M4", "M4_MI20"), "main",
      model_id == "GA_M0", "package_reference",
      default = "sensitivity"
    )]
  }
  markov_registry[, mortality_scope := mortality_scope_from_analysis(analysis_id, cohort)]
  markov_registry <- unique(markov_registry, by = c("cohort", "spec_id"))
}

markov_status_map <- if (nrow(markov_registry)) markov_registry[, .(
  analysis_id, model_id, spec_id, registry_status = status,
  registry_reason = reason, cohort, mortality_scope
)] else data.table(
  analysis_id = character(), model_id = character(), spec_id = character(),
  registry_status = character(), registry_reason = character(), cohort = character(),
  mortality_scope = character()
)

add_markov_terminal_context <- function(x, actual_result = TRUE) {
  if (!nrow(x)) return(x)
  x <- ensure_context(x)
  existing_context <- intersect(c("spec_id", "registry_status", "registry_reason"), names(x))
  if (length(existing_context)) x[, (existing_context) := NULL]
  x <- merge(x, markov_status_map, by = c("analysis_id", "model_id"), all.x = TRUE, sort = FALSE, suffixes = c("", ".registry"))
  if ("cohort.registry" %in% names(x)) {
    x[is.na(cohort) | !nzchar(cohort), cohort := cohort.registry]
    x[, cohort.registry := NULL]
  }
  if ("mortality_scope.registry" %in% names(x)) {
    x[is.na(mortality_scope) | !nzchar(mortality_scope), mortality_scope := mortality_scope.registry]
    x[, mortality_scope.registry := NULL]
  }
  x[is.na(spec_id) | !nzchar(spec_id), spec_id := markov_spec_id(analysis_id, model_id)]
  x[, status := fcase(
    !is.na(registry_status), registry_status,
    actual_result, "completed",
    default = "not_estimable"
  )]
  x[, reason := fcase(
    status == "completed", "",
    !is.na(registry_reason) & nzchar(clean_text(registry_reason)), clean_text(registry_reason),
    default = "No completed upstream result was available for this record"
  )]
  x[, c("registry_status", "registry_reason") := NULL]
  x[]
}

markov_q <- collect_markov(
  "q_intensities", "markov_q_intensities.csv",
  c("analysis_id", "model_id", "transition", "from_state", "to_state", "covariate_profile"),
  c("intensity", "standard_error", "conf_low", "conf_high")
)
markov_q_sources <- attr(markov_q, "source_files")
markov_q <- add_markov_terminal_context(markov_q)

markov_hr <- collect_markov(
  "hazard_ratios", "markov_hazard_ratios.csv",
  c("analysis_id", "model_id", "transition", "from_state", "to_state", "term"),
  c("hazard_ratio", "conf_low", "conf_high", "p_value")
)
markov_hr_sources <- attr(markov_hr, "source_files")
markov_hr <- add_markov_terminal_context(markov_hr)
if (nrow(markov_hr)) {
  if (!"is_requested_insomnia_exposure" %in% names(markov_hr)) markov_hr[, is_requested_insomnia_exposure := 0L]
  markov_hr[, confirmatory_nonfatal := as.integer(
    as.integer(is_requested_insomnia_exposure) == 1L &
      transition %chin% c("S0->S1", "S0->S2", "S1->S3", "S2->S3")
  )]
  if (!"events_per_parameter" %in% names(markov_hr)) markov_hr[, events_per_parameter := NA_real_]
  markov_hr[, event_parameter_risk := fcase(
    is.finite(events_per_parameter) & events_per_parameter >= 10, "adequate: EPP >=10",
    is.finite(events_per_parameter) & events_per_parameter >= 5, "caution: EPP 5-<10",
    is.finite(events_per_parameter), "sparse: EPP <5",
    default = "not assessed"
  )]
  markov_hr[, `:=`(p_fdr_bh_confirmatory = NA_real_, fdr_family_size = NA_integer_)]
  if ("p_value" %in% names(markov_hr)) {
    markov_hr[
      status == "completed" & confirmatory_nonfatal == 1L & is.finite(p_value),
      `:=`(
        p_fdr_bh_confirmatory = p.adjust(p_value, method = "BH"),
        fdr_family_size = .N
      ),
      by = .(analysis_id, model_id, term)
    ]
  }
}

markov_probability <- collect_markov(
  "transition_probabilities", "markov_transition_probabilities_2_5_10y.csv",
  c("analysis_id", "model_id", "exposure_variable", "profile", "time_years", "from_state", "to_state"),
  c("probability", "conf_low", "conf_high", "row_sum_error")
)
markov_probability_sources <- attr(markov_probability, "source_files")
markov_probability <- add_markov_terminal_context(markov_probability)
if (nrow(markov_probability)) {
  # R/09 recognises `time` as the prediction-horizon grouping key.
  markov_probability[, time := if ("time_years" %in% names(markov_probability)) {
    as.numeric(time_years)
  } else NA_real_]
}

markov_sojourn <- collect_markov(
  "sojourn_times", "markov_sojourn_times.csv",
  c("analysis_id", "model_id", "state", "state_index"),
  c("sojourn_years", "sojourn_time", "conf_low", "conf_high")
)
markov_sojourn_sources <- attr(markov_sojourn, "source_files")
markov_sojourn <- add_markov_terminal_context(markov_sojourn)
if (nrow(markov_sojourn) && !"sojourn_time" %in% names(markov_sojourn)) {
  source_col <- first_existing_column(markov_sojourn, c("sojourn_years", "estimate", "years"))
  markov_sojourn[, sojourn_time := if (!is.na(source_col)) as.numeric(get(source_col)) else NA_real_]
}

markov_expected_time <- collect_markov(
  "expected_state_time", "markov_expected_state_time_10y.csv",
  c("analysis_id", "model_id", "profile", "horizon_years", "state", "summary_measure"),
  c("expected_years", "conf_low", "conf_high", "empirical_low", "empirical_high")
)
markov_expected_time_sources <- attr(markov_expected_time, "source_files")
markov_expected_time <- add_markov_terminal_context(markov_expected_time)

markov_prevalence <- collect_markov(
  "observed_expected_prevalence", "markov_observed_expected_prevalence.csv",
  c("analysis_id", "model_id", "time_years", "state", "state_index"),
  c("observed_probability", "expected_probability", "observed_count", "expected_count")
)
markov_prevalence_sources <- attr(markov_prevalence, "source_files")
markov_prevalence <- add_markov_terminal_context(markov_prevalence)

markov_wald <- collect_markov(
  "confirmatory_wald", "markov_confirmatory_joint_wald.csv",
  c("analysis_id", "model_id", "exposure_variable", "tested_transitions"),
  c("statistic", "df", "p_value", "confirmatory_edge_count")
)
markov_wald_sources <- attr(markov_wald, "source_files")
markov_wald <- add_markov_terminal_context(markov_wald)

placeholder_for_missing_models <- function(actual, registry, artifact, schema) {
  if (!nrow(registry)) return(actual)
  present <- if (nrow(actual) && all(c("analysis_id", "model_id") %in% names(actual))) {
    unique(paste(actual$analysis_id, actual$model_id, sep = "::"))
  } else character()
  missing <- registry[!paste(analysis_id, model_id, sep = "::") %chin% present]
  if (!nrow(missing)) return(actual)
  placeholder <- missing[, .(
    analysis_id, model_id, cohort, spec_id, mortality_scope,
    status = fifelse(status == "completed", "failed", status),
    reason = fifelse(
      status == "completed",
      paste0("The model fit completed, but no ", artifact, " result table was found"),
      reason
    )
  )]
  for (name in names(schema)) placeholder[, (name) := schema[[name]]]
  rbindlist(list(actual, placeholder), use.names = TRUE, fill = TRUE)
}

markov_q <- placeholder_for_missing_models(markov_q, markov_registry, "Q-intensity", list(
  transition = NA_character_, from_state = NA_integer_, to_state = NA_integer_,
  intensity = NA_real_, standard_error = NA_real_, conf_low = NA_real_, conf_high = NA_real_
))
markov_hr <- placeholder_for_missing_models(markov_hr, markov_registry, "intensity-ratio", list(
  transition = NA_character_, from_state = NA_integer_, to_state = NA_integer_, term = NA_character_,
  hazard_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_
))
markov_probability <- placeholder_for_missing_models(markov_probability, markov_registry, "transition-probability", list(
  transition = NA_character_, from_state = NA_integer_, to_state = NA_integer_,
  profile = NA_character_, time_years = NA_real_, probability = NA_real_
))
markov_sojourn <- placeholder_for_missing_models(markov_sojourn, markov_registry, "sojourn-time", list(
  state = NA_integer_, state_index = NA_integer_, sojourn_time = NA_real_,
  conf_low = NA_real_, conf_high = NA_real_
))
markov_prevalence <- placeholder_for_missing_models(markov_prevalence, markov_registry, "state-occupation", list(
  state = NA_integer_, state_index = NA_integer_, time_years = NA_real_,
  observed_probability = NA_real_, expected_probability = NA_real_
))
markov_expected_time <- placeholder_for_missing_models(markov_expected_time, markov_registry, "expected-state-time", list(
  state = NA_integer_, profile = NA_character_, horizon_years = NA_real_, expected_years = NA_real_
))
markov_wald <- placeholder_for_missing_models(markov_wald, markov_registry, "confirmatory joint-Wald", list(
  exposure_variable = NA_character_, tested_transitions = NA_character_, statistic = NA_real_,
  df = NA_real_, p_value = NA_real_
))

markov_occupation <- copy(markov_prevalence)
if (nrow(markov_occupation)) {
  if (!"expected_probability" %in% names(markov_occupation)) markov_occupation[, expected_probability := NA_real_]
  markov_occupation[, `:=`(
    probability = as.numeric(expected_probability),
    occupation_probability = as.numeric(expected_probability),
    profile = "sample_means",
    # R/09 recognises `time` as the occupation-curve grouping key.
    time = if ("time_years" %in% names(markov_occupation)) as.numeric(time_years) else NA_real_
  )]
}

# The verifier treats the canonical intensity-estimate table as transition
# intensity ratios (HR/CI/p), so it is an alias of the final hazard-ratio table.
# Absolute Q estimates remain available separately in markov_q_intensities.csv.
markov_intensity <- copy(markov_hr)
if (nrow(markov_intensity)) markov_intensity[, estimate_type := "hazard_ratio"]

# -----------------------------------------------------------------------------
# Cox run-labelled results, diagnostics and specification coverage
# -----------------------------------------------------------------------------

cox_root <- file.path(result_root, "cox_sensitivity")
find_named_recursive <- function(root, filename) {
  if (!dir.exists(root)) return(character())
  list.files(root, pattern = paste0("^", gsub("[.]", "[.]", filename), "$"), recursive = TRUE, full.names = TRUE)
}

normalise_directory <- function(path) {
  tolower(normalizePath(path, winslash = "/", mustWork = FALSE))
}

# R/06 v2 isolates every execution below
#   results/cox_sensitivity/<fingerprinted-run-label>/<run-instance>/
# and writes COMPLETION_MANIFEST.csv only after all selected specifications have
# reached terminal checkpoints and all final CSVs have been written.  A recursive
# filename search alone can see files from a process that stopped between those
# writes, so claim-bearing tables must be restricted to validated completed run
# directories.  Legacy pre-v2 folders remain readable only when no stable-ID
# catalog is available.
cox_completion_files <- find_named_recursive(cox_root, "COMPLETION_MANIFEST.csv")
cox_completion_manifests <- safe_read_many(cox_completion_files, source_priority = 1L)
cox_completion_sources <- attr(cox_completion_manifests, "source_files")
cox_completed_run_directories <- character()
cox_valid_completion_sources <- character()
cox_completion_directory_map <- data.table()
if (nrow(cox_completion_manifests)) {
  # __source_file is relative to output_root; recover an absolute directory for
  # exact path matching without relying on run-label naming conventions.
  cox_completion_manifests[, completion_directory := normalise_directory(
    file.path(output_root, dirname(get("__source_file")))
  )]
  for (column in c(
    "expected_specifications", "checkpoint_specifications", "run_records"
  )) {
    if (!column %in% names(cox_completion_manifests)) {
      cox_completion_manifests[, (column) := NA_real_]
    } else {
      cox_completion_manifests[, (column) := suppressWarnings(as.numeric(get(column)))]
    }
  }
  if (!"run_fingerprint" %in% names(cox_completion_manifests)) {
    cox_completion_manifests[, run_fingerprint := NA_character_]
  }
  if (!"status" %in% names(cox_completion_manifests)) {
    cox_completion_manifests[, status := NA_character_]
  }
  cox_completion_manifests[, manifest_valid :=
    normalise_name(status) == "completed" &
      nzchar(clean_text(run_fingerprint)) &
      is.finite(expected_specifications) & expected_specifications > 0 &
      is.finite(checkpoint_specifications) &
      is.finite(run_records) &
      expected_specifications == checkpoint_specifications &
      expected_specifications == run_records
  ]
  cox_completed_run_directories <- unique(
    cox_completion_manifests[manifest_valid == TRUE, completion_directory]
  )
  valid_manifests <- cox_completion_manifests[manifest_valid == TRUE]
  if (nrow(valid_manifests)) {
    if (!"completed_at" %in% names(valid_manifests)) valid_manifests[, completed_at := NA_character_]
    valid_manifests[, `__completion_time_parsed` := parse_completion_times(completed_at)]
    valid_manifests[, `__completion_time` := fifelse(
      is.finite(get("__completion_time_parsed")),
      get("__completion_time_parsed"), as.numeric(get("__source_mtime"))
    )]
    valid_manifests[, `__completion_time_source` := fifelse(
      is.finite(get("__completion_time_parsed")), "manifest_completed_at", "manifest_mtime_fallback"
    )]
    setorderv(valid_manifests, c("completion_directory", "__completion_time", "__source_mtime"),
      c(1L, -1L, -1L), na.last = TRUE)
    cox_completion_directory_map <- valid_manifests[, .SD[1L], by = "completion_directory"][, .(
      `__run_directory` = completion_directory,
      `__implementation_completed_at` = clean_text(completed_at),
      `__implementation_completed_time` = as.numeric(get("__completion_time")),
      `__implementation_completion_time_source` = clean_text(get("__completion_time_source")),
      `__implementation_completion_source_file` = clean_text(get("__source_file"))
    )]
    cox_valid_completion_sources <- source_files_from_rows(valid_manifests)
  }
}
if (length(cox_completion_files) && !length(cox_completed_run_directories)) {
  stop("Cox completion manifests were found, but none passed completion validation")
}

is_stable_cox_spec_id <- function(x) {
  grepl("^CXS_[0-9a-fA-F]{20}$", clean_text(x))
}

completed_cox_artifact_files <- function(filename, stable_mode) {
  paths <- find_named_recursive(cox_root, filename)
  if (!length(paths)) return(character())
  if (!stable_mode) return(paths)
  if (!length(cox_completed_run_directories)) return(character())
  paths[normalise_directory(dirname(paths)) %chin% cox_completed_run_directories]
}

cox_catalog_paths <- find_named_recursive(cox_root, "transition_cox_model_specifications.csv")
if (length(cox_completion_files)) {
  completed_catalog_paths <- cox_catalog_paths[
    normalise_directory(dirname(cox_catalog_paths)) %chin% cox_completed_run_directories
  ]
  # Once completion manifests exist, never fall back to an incomplete catalog.
  # An empty completed set is a provenance error, not permission to fail open.
  if (!length(completed_catalog_paths)) {
    stop("Validated completed Cox runs were found, but none contained a specification catalog")
  }
  cox_catalog_paths <- completed_catalog_paths
}
cox_catalog <- safe_read_many(cox_catalog_paths)
cox_catalog_all_sources <- attr(cox_catalog, "source_files")
cox_catalog_sources <- cox_catalog_all_sources
cox_run_evidence <- safe_read_many(completed_cox_artifact_files(
  "transition_cox_model_runs.csv", TRUE
))
cox_run_evidence_sources <- attr(cox_run_evidence, "source_files")
cox_run_evidence_summary <- data.table()
if (nrow(cox_run_evidence) && "spec_id" %in% names(cox_run_evidence)) {
  cox_run_evidence <- cox_run_evidence[is_stable_cox_spec_id(spec_id)]
  if (!"status" %in% names(cox_run_evidence)) cox_run_evidence[, status := NA_character_]
  cox_run_evidence[, `__run_directory` := normalise_directory(file.path(
    output_root, dirname(get("__source_file"))
  ))]
  if (nrow(cox_completion_directory_map)) {
    cox_run_evidence <- merge(
      cox_run_evidence, cox_completion_directory_map,
      by = "__run_directory", all.x = TRUE, sort = FALSE
    )
  }
  if (!"__implementation_completed_time" %in% names(cox_run_evidence)) {
    cox_run_evidence[, `__implementation_completed_time` := as.numeric(get("__source_mtime"))]
  }
  cox_run_evidence[!is.finite(get("__implementation_completed_time")),
    `__implementation_completed_time` := as.numeric(get("__source_mtime"))]
  for (column in c(
    "__implementation_completed_at", "__implementation_completion_time_source",
    "__implementation_completion_source_file"
  )) {
    if (!column %in% names(cox_run_evidence)) cox_run_evidence[, (column) := ""]
    cox_run_evidence[is.na(get(column)), (column) := ""]
  }
  cox_run_evidence[, `:=`(
    `__run_terminal_status` = terminal_status(status),
    `__run_terminal_rank` = status_rank(terminal_status(status))
  )]
  setorderv(
    cox_run_evidence,
    c("spec_id", "__implementation_completed_time", "__run_terminal_rank", "__source_mtime"),
    c(1L, -1L, -1L, -1L), na.last = TRUE
  )
  cox_run_evidence_summary <- cox_run_evidence[, .SD[1L], by = "spec_id"][, .(
    spec_id,
    `__completed_run_present` = 1L,
    `__implementation_completed_at` = clean_text(get("__implementation_completed_at")),
    `__implementation_completed_time` = as.numeric(get("__implementation_completed_time")),
    `__implementation_completion_time_source` = clean_text(get("__implementation_completion_time_source")),
    `__implementation_completion_source_file` = clean_text(get("__implementation_completion_source_file")),
    `__implementation_run_mtime` = as.numeric(get("__source_mtime")),
    `__implementation_run_status` = clean_text(get("__run_terminal_status")),
    `__implementation_run_status_rank` = as.integer(get("__run_terminal_rank")),
    `__implementation_run_source_file` = clean_text(get("__source_file"))
  )]
}
cox_specification_supersession_audit <- data.table()
if (nrow(cox_catalog)) {
  if (!"spec_id" %in% names(cox_catalog)) stop("Cox specification catalog has no spec_id column")
  # R/06 v2 uses a stable content hash (`CXS_` plus 20 hexadecimal digits).
  # If both the legacy sequential catalog and the v2 catalog are present in a
  # restart directory, retain only the stable universe so old CXS00001-style
  # identifiers cannot be mistaken for additional prespecified analyses.
  stable_spec <- is_stable_cox_spec_id(cox_catalog$spec_id)
  if (any(stable_spec) && !length(cox_completion_files)) {
    stop("A stable-ID Cox catalog was found without any completion manifest; refusing incomplete provenance")
  }
  if (any(stable_spec)) cox_catalog <- cox_catalog[stable_spec]
  if (!"cohort" %in% names(cox_catalog)) cox_catalog[, cohort := NA_character_]
  cox_catalog[, cohort := toupper(clean_text(cohort))]
  cox_catalog[, prespecified := TRUE]
  # A corrected implementation can legitimately change formula_rhs and hence
  # the content-hashed spec_id while retaining the same prespecified analytic
  # question (for example, replacing a duplicated diabetes coding in an
  # interaction model).  Select the newest completed implementation by a
  # formula-independent semantic key, and retain an explicit supersession audit
  # so an obsolete failed specification cannot survive beside its correction.
  cox_semantic_fields <- c(
    "cohort", "transition_id", "from_state", "to_state", "module",
    "analysis_family", "analysis_id", "scenario", "exposure_var",
    "exposure_type", "interaction_var", "model_level", "basis", "weighted"
  )
  missing_semantic_fields <- setdiff(cox_semantic_fields, names(cox_catalog))
  if (length(missing_semantic_fields)) {
    stop("Cox catalog lacks semantic identity fields: ", paste(missing_semantic_fields, collapse = ", "))
  }
  cox_catalog[, `__cox_semantic_key` := do.call(
    paste,
    c(lapply(.SD, function(value) fifelse(is.na(value), "<NA>", clean_text(value))), sep = "||")
  ), .SDcols = cox_semantic_fields]
  cox_catalog_versions <- deduplicate_rows(
    copy(cox_catalog), c("__cox_semantic_key", "spec_id"),
    c("analysis_id", "analysis_label", "transition", "formula_rhs", "mortality_scope"),
    prefer_latest_source = TRUE
  )
  if (nrow(cox_run_evidence_summary)) {
    cox_catalog_versions <- merge(
      cox_catalog_versions, cox_run_evidence_summary,
      by = "spec_id", all.x = TRUE, sort = FALSE
    )
  }
  for (column in c("__completed_run_present", "__implementation_run_status_rank")) {
    if (!column %in% names(cox_catalog_versions)) cox_catalog_versions[, (column) := 0L]
    cox_catalog_versions[is.na(get(column)), (column) := 0L]
  }
  for (column in c("__implementation_completed_time", "__implementation_run_mtime")) {
    if (!column %in% names(cox_catalog_versions)) cox_catalog_versions[, (column) := -Inf]
    cox_catalog_versions[is.na(get(column)), (column) := -Inf]
  }
  for (column in c(
    "__implementation_completed_at", "__implementation_completion_time_source",
    "__implementation_completion_source_file", "__implementation_run_status",
    "__implementation_run_source_file"
  )) {
    if (!column %in% names(cox_catalog_versions)) cox_catalog_versions[, (column) := ""]
    cox_catalog_versions[is.na(get(column)), (column) := ""]
  }
  # Prefer evidence from an actually selected, validated completed run.  Within
  # a semantic question, the most recent completed implementation supersedes
  # older content hashes; terminal rank only breaks an exact-time tie.
  setorderv(
    cox_catalog_versions,
    c("__cox_semantic_key", "__completed_run_present", "__implementation_completed_time",
      "__implementation_run_status_rank", "__implementation_run_mtime", "__source_mtime"),
    c(1L, -1L, -1L, -1L, -1L, -1L), na.last = TRUE
  )
  cox_catalog <- cox_catalog_versions[, .SD[1L], by = "__cox_semantic_key"]
  semantic_selection <- cox_catalog[, .(
    `__cox_semantic_key`, selected_spec_id = spec_id,
    selected_formula_rhs = clean_text(formula_rhs),
    selected_completed_at = clean_text(get("__implementation_completed_at")),
    selected_source_file = clean_text(get("__source_file")),
    selected_source_mtime = as.numeric(get("__source_mtime"))
  )]
  superseded_keys <- cox_catalog_versions[, .(spec_id_versions = uniqueN(spec_id)),
    by = "__cox_semantic_key"][spec_id_versions > 1L, `__cox_semantic_key`]
  if (length(superseded_keys)) {
    cox_specification_supersession_audit <- merge(
      cox_catalog_versions[get("__cox_semantic_key") %chin% superseded_keys],
      semantic_selection,
      by = "__cox_semantic_key", all.x = TRUE, sort = FALSE
    )
    cox_specification_supersession_audit[, `:=`(
      implementation_catalog_source_file = clean_text(get("__source_file")),
      implementation_catalog_source_mtime = as.numeric(get("__source_mtime")),
      implementation_completed_at = clean_text(get("__implementation_completed_at")),
      implementation_completed_time = as.numeric(get("__implementation_completed_time")),
      implementation_completion_time_source = clean_text(get("__implementation_completion_time_source")),
      implementation_completion_source_file = clean_text(get("__implementation_completion_source_file")),
      implementation_run_source_file = clean_text(get("__implementation_run_source_file")),
      implementation_run_mtime = as.numeric(get("__implementation_run_mtime")),
      implementation_run_status = clean_text(get("__implementation_run_status")),
      superseded = as.integer(spec_id != selected_spec_id),
      supersession_reason = fifelse(
        spec_id == selected_spec_id,
        "Most recent validated completed implementation retained",
        "Older content-hashed implementation replaced by a more recent validated completed semantic equivalent"
      )
    )]
    # `selected` in an upstream catalog only describes that source shard.  It is
    # not the final cross-run decision, so expose it under an unambiguous name
    # and publish the actual retained/superseded decision separately.
    if ("selected" %in% names(cox_specification_supersession_audit)) {
      setnames(cox_specification_supersession_audit, "selected", "catalog_selected_at_source")
    }
    if ("selection_reason" %in% names(cox_specification_supersession_audit)) {
      setnames(
        cox_specification_supersession_audit, "selection_reason",
        "catalog_selection_reason_at_source"
      )
    }
    cox_specification_supersession_audit[, implementation_retained := as.integer(superseded == 0L)]
  }
  if ("__source_file" %in% names(cox_catalog)) {
    cox_catalog_sources <- existing_files(file.path(
      output_root, unique(clean_text(cox_catalog[["__source_file"]]))
    ))
  }
}
cox_stable_mode <- nrow(cox_catalog) && any(is_stable_cox_spec_id(cox_catalog$spec_id))

cox_runs <- safe_read_many(completed_cox_artifact_files(
  "transition_cox_model_runs.csv", cox_stable_mode
))
cox_run_sources <- attr(cox_runs, "source_files")
if (nrow(cox_runs)) {
  if (!"spec_id" %in% names(cox_runs)) stop("Cox model-run table has no spec_id column")
  stable_spec_universe <- if (nrow(cox_catalog)) unique(cox_catalog$spec_id) else character()
  if (length(stable_spec_universe)) cox_runs <- cox_runs[spec_id %chin% stable_spec_universe]
  if (!"status" %in% names(cox_runs)) cox_runs[, status := NA_character_]
  cox_runs[, upstream_status := clean_text(status)]
  cox_runs[, status := terminal_status(upstream_status)]
  cox_runs[is.na(status), status := fifelse(grepl("not", normalise_name(upstream_status)), "not_estimable", "failed")]
  if (!"message" %in% names(cox_runs)) cox_runs[, message := ""]
  cox_runs[, reason := clean_text(message)]
  cox_runs[status == "completed", reason := ""]
  cox_runs[status != "completed" & !nzchar(reason), reason := paste0("Upstream Cox terminal status: ", upstream_status)]
  if (cox_stable_mode) {
    cox_runs[, `__run_directory` := normalise_directory(file.path(
      output_root, dirname(get("__source_file"))
    ))]
    if (nrow(cox_completion_directory_map)) {
      cox_runs <- merge(
        cox_runs, cox_completion_directory_map,
        by = "__run_directory", all.x = TRUE, sort = FALSE
      )
    }
    if (!"__implementation_completed_time" %in% names(cox_runs)) {
      cox_runs[, `__implementation_completed_time` := as.numeric(get("__source_mtime"))]
    }
    cox_runs[!is.finite(get("__implementation_completed_time")),
      `__implementation_completed_time` := as.numeric(get("__source_mtime"))]
    cox_runs[, `__run_selection_status_rank` := status_rank(status)]
    setorderv(
      cox_runs,
      c("spec_id", "__implementation_completed_time", "__run_selection_status_rank", "__source_mtime"),
      c(1L, -1L, -1L, -1L), na.last = TRUE
    )
    cox_runs <- cox_runs[, .SD[1L], by = "spec_id"]
  } else {
    cox_runs <- deduplicate_rows(
      cox_runs, "spec_id",
      c("imputations_successful", "minimum_events", "minimum_events_per_parameter", "message")
    )
  }
}
cox_run_sources <- source_files_from_rows(cox_runs)
if (nrow(cox_catalog)) {
  selected_cox_specs <- if (nrow(cox_runs)) unique(clean_text(cox_runs$spec_id)) else character()
  cox_catalog[, selected := spec_id %chin% selected_cox_specs]
  cox_catalog[, selection_reason := fifelse(
    selected,
    "Selected in a validated completed Cox run",
    "Not selected in any validated completed Cox run"
  )]
}

cox_run_source_map <- data.table()
cox_selected_completion_sources <- character()
if (nrow(cox_runs) && "__source_file" %in% names(cox_runs)) {
  cox_runs[, `__run_directory` := normalise_directory(file.path(
    output_root, dirname(get("__source_file"))
  ))]
  cox_run_source_map <- unique(cox_runs[, .(
    spec_id, `__selected_run_directory` = get("__run_directory")
  )], by = "spec_id")
  selected_run_directories <- unique(cox_run_source_map[["__selected_run_directory"]])
  if (nrow(cox_completion_manifests) && "manifest_valid" %in% names(cox_completion_manifests)) {
    cox_selected_completion_sources <- source_files_from_rows(
      cox_completion_manifests[
        manifest_valid == TRUE & completion_directory %chin% selected_run_directories
      ]
    )
  }
}

filter_to_selected_cox_run <- function(x, id_column) {
  if (!cox_stable_mode || !nrow(x)) return(x)
  if (!nrow(cox_run_source_map) ||
      !id_column %in% names(x) ||
      !"__source_file" %in% names(x)) return(x[0L])
  selected <- copy(cox_run_source_map)
  setnames(selected, "spec_id", id_column)
  x[, `__artifact_run_directory` := normalise_directory(file.path(
    output_root, dirname(get("__source_file"))
  ))]
  x <- merge(x, selected, by = id_column, all = FALSE, sort = FALSE)
  x[
    get("__artifact_run_directory") == get("__selected_run_directory")
  ]
}

cox_estimates <- safe_read_many(completed_cox_artifact_files(
  "transition_cox_estimates.csv", cox_stable_mode
))
cox_estimate_sources <- attr(cox_estimates, "source_files")
if (nrow(cox_estimates)) {
  if (!"spec_id" %in% names(cox_estimates)) stop("Cox estimate table has no spec_id column")
  stable_spec_universe <- if (nrow(cox_catalog)) unique(cox_catalog$spec_id) else character()
  if (length(stable_spec_universe)) cox_estimates <- cox_estimates[spec_id %chin% stable_spec_universe]
  if (cox_stable_mode) {
    if (nrow(cox_run_source_map) && "__source_file" %in% names(cox_estimates)) {
      cox_estimates[, `__estimate_run_directory` := normalise_directory(file.path(
        output_root, dirname(get("__source_file"))
      ))]
      cox_estimates <- merge(
        cox_estimates, cox_run_source_map,
        by = "spec_id", all = FALSE, sort = FALSE
      )
      cox_estimates <- cox_estimates[
        get("__estimate_run_directory") == get("__selected_run_directory")
      ]
    } else {
      cox_estimates <- cox_estimates[0L]
    }
  }
  run_state <- if (nrow(cox_runs)) unique(cox_runs[, .(
    spec_id, run_status = status, run_reason = reason
  )], by = "spec_id") else data.table(spec_id = character(), run_status = character(), run_reason = character())
  cox_estimates <- merge(cox_estimates, run_state, by = "spec_id", all.x = TRUE, sort = FALSE)
  cox_estimates[, status := fifelse(is.na(run_status), "completed", run_status)]
  cox_estimates[, reason := fifelse(status == "completed", "", clean_text(run_reason))]
  if (!"term" %in% names(cox_estimates)) cox_estimates[, term := NA_character_]
  cox_estimates <- deduplicate_rows(
    cox_estimates, c("spec_id", "term"),
    c("hazard_ratio", "conf_low", "conf_high", "p_value", "m_used")
  )
}
cox_estimate_sources <- source_files_from_rows(cox_estimates)

estimated_cox_specs <- if (nrow(cox_estimates)) unique(cox_estimates[
  status == "completed" &
    is.finite(hazard_ratio) & is.finite(conf_low) & is.finite(conf_high) & is.finite(p_value), spec_id
]) else character()

cox_registry <- copy(cox_catalog)
if (nrow(cox_registry)) {
  run_state <- if (nrow(cox_runs)) unique(cox_runs[, .(
    spec_id, run_status = status, run_reason = reason,
    upstream_run_status = upstream_status
  )], by = "spec_id") else data.table(
    spec_id = character(), run_status = character(), run_reason = character(), upstream_run_status = character()
  )
  cox_registry <- merge(cox_registry, run_state, by = "spec_id", all.x = TRUE, sort = FALSE)
  cox_registry[, status := fcase(
    spec_id %chin% estimated_cox_specs, "completed",
    run_status == "completed", "failed",
    !is.na(run_status), run_status,
    default = "not_estimable"
  )]
  if (!"selection_reason" %in% names(cox_registry)) cox_registry[, selection_reason := ""]
  cox_registry[, reason := fcase(
    status == "completed", "",
    run_status == "completed" & !spec_id %chin% estimated_cox_specs,
    "The Cox run was marked successful, but no estimate row was present at collation time",
    !is.na(run_reason) & nzchar(clean_text(run_reason)), clean_text(run_reason),
    nzchar(clean_text(selection_reason)),
    paste0("No completed run was present at collation time (catalog selection: ", clean_text(selection_reason), ")"),
    default = "No completed run was present at collation time; this prespecified Cox sensitivity was not run or is unavailable"
  )]
  cox_registry[, module := paste0("cox_", clean_text(module))]
  cox_registry[, mortality_scope := fcase(
    cohort == "ELSA" & as.integer(to_state) == 4L,
    "ELSA limited mortality validation through 2012; death ascertainment ends in 2012",
    cohort == "ELSA",
    "ELSA W4-W10 nonfatal follow-up; mortality analyses are limited through 2012",
    default = "HRS follow-up through 2020"
  )]
  cox_registry <- unique(cox_registry, by = c("cohort", "spec_id"))
}

if (nrow(cox_registry)) {
  present <- if (nrow(cox_estimates)) unique(cox_estimates$spec_id) else character()
  absent <- cox_registry[!spec_id %chin% present]
  if (nrow(absent)) {
    placeholders <- absent[, .(
      spec_id, cohort, transition_id, transition, from_state, to_state,
      analysis_id, analysis_label, scenario, analysis_family, model_level,
      basis, weighted, exposure_var, mortality_scope, status, reason,
      term = NA_character_, estimate = NA_real_, robust_se = NA_real_,
      hazard_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
      p_value = NA_real_
    )]
    cox_estimates <- rbindlist(list(cox_estimates, placeholders), use.names = TRUE, fill = TRUE)
  }
}

cox_diag_actual <- safe_read_many(completed_cox_artifact_files(
  "transition_cox_diagnostics.csv", cox_stable_mode
), source_priority = 1L)
cox_diag_sources <- attr(cox_diag_actual, "source_files")
cox_diag_actual <- filter_to_selected_cox_run(cox_diag_actual, "reference_spec_id")
if (nrow(cox_diag_actual)) {
  if (!"status" %in% names(cox_diag_actual)) cox_diag_actual[, status := NA_character_]
  cox_diag_actual[, upstream_status := clean_text(status)]
  cox_diag_actual[, status := terminal_status(upstream_status)]
  cox_diag_actual[is.na(status), status := "not_estimable"]
  if (!"reason" %in% names(cox_diag_actual)) cox_diag_actual[, reason := ""]
  cox_diag_actual[status != "completed" & !nzchar(clean_text(reason)),
    reason := paste0("Upstream Cox diagnostic status: ", upstream_status)]
  cox_diag_actual[status == "completed", reason := ""]
}

cox_edge_registry <- safe_read_many(completed_cox_artifact_files(
  "transition_cox_edge_registry.csv", cox_stable_mode
), source_priority = 4L)
cox_edge_sources <- attr(cox_edge_registry, "source_files")
cox_edge_registry <- filter_to_selected_cox_run(cox_edge_registry, "reference_spec_id")
if (nrow(cox_edge_registry)) {
  if (!"reference_spec_id" %in% names(cox_edge_registry)) cox_edge_registry[, reference_spec_id := NA_character_]
  run_state <- if (nrow(cox_runs)) unique(cox_runs[, .(
    reference_spec_id = spec_id, reference_run_status = status,
    reference_run_reason = reason
  )], by = "reference_spec_id") else data.table(
    reference_spec_id = character(), reference_run_status = character(), reference_run_reason = character()
  )
  cox_edge_registry <- merge(cox_edge_registry, run_state, by = "reference_spec_id", all.x = TRUE, sort = FALSE)
  if (!"status" %in% names(cox_edge_registry)) cox_edge_registry[, status := NA_character_]
  cox_edge_registry[, preflight_status := clean_text(status)]
  if (!"reason" %in% names(cox_edge_registry)) cox_edge_registry[, reason := ""]
  cox_edge_registry[, status := fifelse(
    reference_run_status == "completed", "completed", "not_estimable"
  )]
  cox_edge_registry[, reason := fcase(
    status == "completed", "",
    !is.na(reference_run_reason) & nzchar(clean_text(reference_run_reason)), clean_text(reference_run_reason),
    default = "Only a Cox preflight edge audit was available; no completed reference diagnostic fit was present"
  )]
  if (!"n_events" %in% names(cox_edge_registry)) {
    cox_edge_registry[, n_events := if ("target_events" %in% names(cox_edge_registry)) as.numeric(target_events) else NA_real_]
  }
  if (!"event_parameter_risk" %in% names(cox_edge_registry)) {
    cox_edge_registry[, event_parameter_risk := fcase(
      is.finite(events_per_parameter) & events_per_parameter >= 10, "adequate",
      is.finite(events_per_parameter) & events_per_parameter >= 5, "caution",
      is.finite(events_per_parameter), "sparse_events",
      default = "not_quantifiable"
    )]
  }
  if (!"ph_global_p_value" %in% names(cox_edge_registry)) cox_edge_registry[, ph_global_p_value := NA_real_]
  cox_edge_registry[, mortality_scope := fcase(
    toupper(cohort) == "ELSA" & as.integer(to_state) == 4L,
    "ELSA limited mortality validation through 2012; death ascertainment ends in 2012",
    toupper(cohort) == "ELSA", "ELSA W4-W10 nonfatal follow-up",
    default = "HRS follow-up through 2020"
  )]
}

cox_diagnostics <- rbindlist(list(cox_diag_actual, cox_edge_registry), use.names = TRUE, fill = TRUE)
if (nrow(cox_diagnostics)) {
  cox_diagnostics <- deduplicate_rows(
    cox_diagnostics, c("cohort", "transition", "from_state", "to_state"),
    c("ph_global_p_value", "events_per_parameter", "n_parameters", "n_events")
  )
}
if (!nrow(cox_diagnostics) && nrow(cox_registry)) {
  cox_diagnostics <- unique(cox_registry[, .(
    cohort, transition_id, transition, from_state, to_state,
    reference_spec_id = spec_id, status = "not_estimable",
    reason = "No Cox diagnostic or edge-preflight output was present at collation time",
    n_events = NA_real_, n_parameters = NA_real_, events_per_parameter = NA_real_,
    event_parameter_risk = "not_quantifiable", ph_global_p_value = NA_real_,
    mortality_scope
  )], by = c("cohort", "transition"))
}
cox_diag_sources <- source_files_from_rows(cox_diagnostics)
cox_edge_sources <- character()

# -----------------------------------------------------------------------------
# Protocol-targeted sensitivities from R/07b
# -----------------------------------------------------------------------------

protocol_root <- file.path(result_root, "protocol_sensitivities")
protocol_qc_root <- file.path(qc_root, "protocol_sensitivities")

collect_protocol <- function(root, filename, keys = character(), values = character()) {
  x <- safe_read_many(find_named_recursive(root, filename))
  sources <- attr(x, "source_files")
  if (nrow(x) && length(keys)) x <- deduplicate_rows(x, keys, values)
  attr(x, "source_files") <- sources
  x
}

protocol_specs <- collect_protocol(
  protocol_root, "protocol_sensitivity_specifications.csv",
  c("cohort", "spec_id"),
  c("sensitivity_id", "scenario", "scenario_label", "transition", "mortality_scope")
)
protocol_spec_sources <- attr(protocol_specs, "source_files")

protocol_model_registry_raw <- collect_protocol(
  protocol_root, "protocol_sensitivity_model_registry.csv",
  c("cohort", "spec_id"),
  c("n_records", "n_subjects", "events", "events_per_parameter", "ph_global_p_value", "reason")
)
protocol_registry_sources <- attr(protocol_model_registry_raw, "source_files")

protocol_estimates <- collect_protocol(
  protocol_root, "protocol_sensitivity_estimates.csv",
  c("cohort", "spec_id", "term", "comparison"),
  c("hazard_ratio", "conf_low", "conf_high", "p_value", "p_fdr_bh_sensitivity")
)
protocol_estimate_sources <- attr(protocol_estimates, "source_files")

protocol_joint_tests <- collect_protocol(
  protocol_root, "protocol_sensitivity_joint_tests.csv",
  c("cohort", "spec_id", "test"),
  c("statistic", "df", "p_value", "tested_terms")
)
protocol_joint_sources <- attr(protocol_joint_tests, "source_files")

protocol_coverage <- collect_protocol(
  report_root, "protocol_sensitivity_coverage.csv",
  c("cohort", "sensitivity_id"),
  c("protocol_sensitivity", "status", "reason", "evidence_files", "completed_models")
)
protocol_coverage_sources <- attr(protocol_coverage, "source_files")

protocol_ipcw_diagnostics <- collect_protocol(
  protocol_root, "ipcw_weight_diagnostics.csv",
  "cohort", c("status", "intervals_with_ipcw", "ipcw_p01", "ipcw_median", "ipcw_p99")
)
protocol_ipcw_sources <- attr(protocol_ipcw_diagnostics, "source_files")

protocol_denial_audit <- collect_protocol(
  protocol_root, "diagnosis_denial_pattern_audit.csv",
  c("cohort", "disease"),
  c("observed_visit_rows", "positive_visit_rows", "later_negative_after_positive_rows")
)
protocol_denial_sources <- attr(protocol_denial_audit, "source_files")

protocol_expanded_hypertension_audit <- collect_protocol(
  protocol_root, "expanded_hypertension_state_audit.csv",
  "cohort",
  c("rows", "primary_state_missing", "expanded_state_missing", "changed_state_rows")
)
protocol_expanded_sources <- attr(protocol_expanded_hypertension_audit, "source_files")

protocol_data_audit <- collect_protocol(
  protocol_qc_root, "protocol_sensitivity_data_audit.csv",
  c("cohort", "audit_item"),
  c("rows", "participants", "events_on_eight_edges", "detail")
)
protocol_data_audit_sources <- attr(protocol_data_audit, "source_files")

# The model registry is specification-level.  Start from the R/07b
# specification catalog so every attempted edge remains represented even when
# no fit row is estimable, then attach the terminal fit state from its registry.
protocol_registry <- copy(protocol_specs)
if (!nrow(protocol_registry) && nrow(protocol_model_registry_raw)) {
  protocol_registry <- copy(protocol_model_registry_raw)
}
if (nrow(protocol_registry)) {
  if (!"cohort" %in% names(protocol_registry)) protocol_registry[, cohort := NA_character_]
  protocol_registry[, cohort := toupper(clean_text(cohort))]
  if (nrow(protocol_model_registry_raw)) {
    protocol_fit_state <- copy(protocol_model_registry_raw)
    if (!"status" %in% names(protocol_fit_state)) protocol_fit_state[, status := NA_character_]
    if (!"reason" %in% names(protocol_fit_state)) protocol_fit_state[, reason := ""]
    protocol_fit_state[, `:=`(
      protocol_fit_status = terminal_status(status),
      protocol_fit_reason = clean_text(reason)
    )]
    protocol_fit_state <- deduplicate_rows(
      protocol_fit_state, c("cohort", "spec_id"),
      c("n_records", "n_subjects", "events", "events_per_parameter", "protocol_fit_reason")
    )
    protocol_fit_state <- protocol_fit_state[, .(
      cohort = toupper(clean_text(cohort)), spec_id,
      protocol_fit_status, protocol_fit_reason
    )]
    protocol_registry <- merge(
      protocol_registry, protocol_fit_state,
      by = c("cohort", "spec_id"), all.x = TRUE, sort = FALSE
    )
  } else {
    protocol_registry[, `:=`(
      protocol_fit_status = NA_character_, protocol_fit_reason = NA_character_
    )]
  }
  protocol_registry[, `:=`(
    module = "protocol_sensitivity",
    analysis_family = paste0("protocol_sensitivity_", sprintf("%02d", as.integer(sensitivity_id))),
    analysis_id = paste0("PS", sprintf("%02d", as.integer(sensitivity_id)), "_", clean_text(scenario)),
    model_id = clean_text(scenario),
    model_label = clean_text(scenario_label),
    exposure = "Insomnia symptom score",
    exposure_var = "insomnia_score",
    exposure_type = fifelse(clean_text(scenario) == "rcs", "restricted_cubic_spline", "continuous_score"),
    model_level = "M4",
    basis = "complete_case",
    weighted = as.integer(clean_text(scenario) == "ipcw"),
    covariates = "M4 covariates; see R/07b_protocol_sensitivities.R",
    missing_data = "complete_case",
    prespecified = TRUE,
    selected = TRUE
  )]
  protocol_registry[, status := fcase(
    !is.na(protocol_fit_status), protocol_fit_status,
    default = "not_estimable"
  )]
  protocol_registry[, reason := fcase(
    status == "completed", "",
    !is.na(protocol_fit_reason) & nzchar(clean_text(protocol_fit_reason)), clean_text(protocol_fit_reason),
    default = "No terminal R/07b fit record was present at collation time"
  )]
  protocol_registry <- unique(protocol_registry, by = c("cohort", "spec_id"))
}

# -----------------------------------------------------------------------------
# Multiple-imputation audit: preserve source evidence and add transparent
# specification-level status rows where m=20 was prespecified but not completed.
# -----------------------------------------------------------------------------

markov_mi_files <- existing_files(c(
  file.path(markov_roots, "multiple_imputation_audit.csv"),
  if (dir.exists(file.path(output_root, "data_multistate", "mi_baseline"))) {
    list.files(file.path(output_root, "data_multistate", "mi_baseline"), pattern = "__audit[.]csv$", full.names = TRUE)
  } else character()
))
markov_mi_audit <- safe_read_many(markov_mi_files)
markov_mi_sources <- attr(markov_mi_audit, "source_files")
if (nrow(markov_mi_audit)) {
  if (!"analysis_id" %in% names(markov_mi_audit)) markov_mi_audit[, analysis_id := NA_character_]
  markov_mi_audit <- ensure_context(markov_mi_audit)
  markov_mi_audit[, model_id := "M4_MI20"]
  markov_mi_audit[, spec_id := markov_spec_id(analysis_id, model_id)]
  requested_col <- first_existing_column(markov_mi_audit, c("m_requested", "m", "mi_m", "n_imputations", "imputations_requested"))
  if (is.na(requested_col)) {
    markov_mi_audit[, m_requested := 20L]
  } else if (requested_col != "m_requested") {
    markov_mi_audit[, m_requested := as.integer(get(requested_col))]
  }
  mi_fit_state <- if (nrow(markov_fit)) markov_fit[model_id == "M4_MI20", .(
    analysis_id, mi_fit_status = status,
    m_completed = if ("m_fitted" %in% names(markov_fit)) as.numeric(m_fitted) else NA_real_,
    mi_fit_reason = reason
  )] else data.table(analysis_id = character(), mi_fit_status = character(), m_completed = numeric(), mi_fit_reason = character())
  markov_mi_audit <- merge(markov_mi_audit, mi_fit_state, by = "analysis_id", all.x = TRUE, sort = FALSE)
  markov_mi_audit[, status := fifelse(is.na(mi_fit_status), "not_estimable", mi_fit_status)]
  markov_mi_audit[, reason := fcase(
    status == "completed", "",
    !is.na(mi_fit_reason) & nzchar(clean_text(mi_fit_reason)), clean_text(mi_fit_reason),
    default = "Baseline imputation audit exists, but a completed 20-imputation pooled Markov model is not evidenced"
  )]
  markov_mi_audit[, audit_module := "markov"]
}

cox_mi_audit <- data.table()
if (nrow(cox_runs) && "basis" %in% names(cox_runs)) {
  cox_mi_audit <- cox_runs[basis == "mice_pooled", .(
    spec_id, cohort, analysis_id, transition, model_level,
    m_requested = as.integer(if ("imputations_requested" %in% names(cox_runs)) imputations_requested else 20L),
    m_completed = as.integer(if ("imputations_successful" %in% names(cox_runs)) imputations_successful else NA_integer_),
    status, reason,
    audit_module = "cox"
  )]
}

mi_audit <- rbindlist(list(markov_mi_audit, cox_mi_audit), use.names = TRUE, fill = TRUE)
mi_prespecified <- rbindlist(list(
  if (nrow(markov_registry)) markov_registry[model_id == "M4_MI20", .(
    spec_id, cohort, analysis_id, model_id, m_requested = 20L,
    m_completed = NA_integer_, status,
    reason = fifelse(status == "completed", "", reason), audit_module = "markov"
  )] else data.table(),
  if (nrow(cox_registry)) cox_registry[basis == "mice_pooled", .(
    spec_id, cohort, analysis_id, transition, model_level,
    m_requested = 20L, m_completed = NA_integer_, status, reason, audit_module = "cox"
  )] else data.table()
), use.names = TRUE, fill = TRUE)
if (nrow(mi_prespecified)) {
  already <- if (nrow(mi_audit)) unique(mi_audit$spec_id) else character()
  mi_audit <- rbindlist(list(mi_audit, mi_prespecified[!spec_id %chin% already]), use.names = TRUE, fill = TRUE)
}
if (!nrow(mi_audit)) {
  mi_audit <- data.table(
    spec_id = "MI::PRESPECIFIED::M20", cohort = "BOTH", analysis_id = "MI20",
    m_requested = 20L, m_completed = NA_integer_, status = "not_estimable",
    reason = "Multiple imputation with m=20 was prespecified, but no MI audit or pooled model result was present",
    audit_module = "cross_module"
  )
}
mi_audit[is.na(m_requested), m_requested := 20L]
mi_audit[, status := terminal_status(status)]
mi_audit[is.na(status), status := "not_estimable"]
mi_audit[status != "completed" & (!"reason" %in% names(mi_audit) | !nzchar(clean_text(reason))),
  reason := "No completed 20-imputation analysis was evidenced at collation time"]
mi_audit <- deduplicate_rows(
  mi_audit,
  intersect(c("audit_module", "cohort", "spec_id", "variable"), names(mi_audit)),
  c("m_requested", "m_completed", "missing_before", "missing_after_across_all_imputations")
)

# -----------------------------------------------------------------------------
# One truthful registry and an exact-matching coverage report
# -----------------------------------------------------------------------------

registry_columns <- unique(c(
  "spec_id", "cohort", "module", "analysis_family", "analysis_id", "analysis_label",
  "sensitivity_id", "scenario", "scenario_label", "model_id", "model_label",
  "transition_id", "transition", "from_state", "to_state",
  "exposure", "exposure_var", "exposure_type", "model_level", "basis", "weighted",
  "covariates", "formula_rhs", "missing_data", "mortality_scope", "prespecified",
  "selected", "selection_reason", "status", "reason"
))
model_registry <- rbindlist(
  list(markov_registry, cox_registry, protocol_registry),
  use.names = TRUE, fill = TRUE
)
if (!nrow(model_registry)) {
  stop("No Markov, Cox or protocol-sensitivity specification catalog was found; a truthful specification registry cannot be created")
}
for (name in setdiff(registry_columns, names(model_registry))) model_registry[, (name) := NA]
model_registry <- model_registry[, ..registry_columns]
model_registry[, cohort := toupper(clean_text(cohort))]
model_registry[, prespecified := TRUE]
model_registry[, status := terminal_status(status)]
model_registry[is.na(status), status := "not_estimable"]
model_registry[status == "completed", reason := ""]
model_registry[status != "completed" & !nzchar(clean_text(reason)),
  reason := "Prespecified analysis has no completed result at collation time"]
model_registry[cohort == "ELSA" & grepl("5state|through_2012", tolower(analysis_id)),
  mortality_scope := "ELSA limited mortality validation through 2012; death ascertainment ends in 2012"]
model_registry <- unique(model_registry, by = c("cohort", "spec_id"))
setorder(model_registry, cohort, module, analysis_id, spec_id)

coverage_matrix <- copy(model_registry)
coverage_matrix[, `:=`(
  terminal = TRUE,
  result_available = status == "completed",
  coverage_note = fcase(
    status == "completed", "Completed result is present in the consolidated outputs",
    status == "not_estimable", "No estimable/completed result is present; see reason",
    default = "Analysis failed; see reason"
  )
)]

# -----------------------------------------------------------------------------
# Canonical writes.  Internal provenance columns are removed from public tables;
# provenance remains available in results_manifest.csv.
# -----------------------------------------------------------------------------

canonical <- list(
  model_specification_registry.csv = model_registry,
  markov_q_intensities.csv = markov_q,
  markov_hazard_ratios.csv = markov_hr,
  markov_transition_intensity_estimates.csv = markov_intensity,
  markov_transition_probability_matrices.csv = markov_probability,
  markov_sojourn_times.csv = markov_sojourn,
  markov_expected_state_time_10y.csv = markov_expected_time,
  markov_observed_expected_prevalence.csv = markov_prevalence,
  markov_state_occupation_predictions.csv = markov_occupation,
  markov_fit_metrics.csv = markov_fit,
  markov_model_diagnostics.csv = markov_fit,
  markov_confirmatory_joint_wald.csv = markov_wald,
  cox_transition_estimates.csv = cox_estimates,
  cox_model_diagnostics.csv = cox_diagnostics,
  cox_run_completion_audit.csv = cox_completion_manifests,
  cox_specification_supersession_audit.csv = cox_specification_supersession_audit,
  multiple_imputation_audit.csv = mi_audit,
  protocol_sensitivity_specifications.csv = protocol_specs,
  protocol_sensitivity_model_registry.csv = protocol_model_registry_raw,
  protocol_sensitivity_estimates.csv = protocol_estimates,
  protocol_sensitivity_joint_tests.csv = protocol_joint_tests,
  protocol_sensitivity_coverage.csv = protocol_coverage,
  ipcw_weight_diagnostics.csv = protocol_ipcw_diagnostics,
  diagnosis_denial_pattern_audit.csv = protocol_denial_audit,
  expanded_hypertension_state_audit.csv = protocol_expanded_hypertension_audit,
  protocol_sensitivity_data_audit.csv = protocol_data_audit
)

source_map <- list(
  model_specification_registry.csv = c(
    markov_spec_sources, markov_fit_sources, markov_failure_sources,
    cox_catalog_sources, cox_run_sources, cox_selected_completion_sources,
    protocol_spec_sources, protocol_registry_sources
  ),
  markov_q_intensities.csv = markov_q_sources,
  markov_hazard_ratios.csv = markov_hr_sources,
  markov_transition_intensity_estimates.csv = markov_hr_sources,
  markov_transition_probability_matrices.csv = markov_probability_sources,
  markov_sojourn_times.csv = markov_sojourn_sources,
  markov_expected_state_time_10y.csv = markov_expected_time_sources,
  markov_observed_expected_prevalence.csv = markov_prevalence_sources,
  markov_state_occupation_predictions.csv = markov_prevalence_sources,
  markov_fit_metrics.csv = markov_fit_sources,
  markov_model_diagnostics.csv = markov_fit_sources,
  markov_confirmatory_joint_wald.csv = markov_wald_sources,
  cox_transition_estimates.csv = c(
    cox_estimate_sources, cox_run_sources, cox_catalog_sources,
    cox_selected_completion_sources
  ),
  cox_model_diagnostics.csv = c(
    cox_diag_sources, cox_edge_sources, cox_run_sources,
    cox_selected_completion_sources
  ),
  cox_run_completion_audit.csv = cox_completion_sources,
  cox_specification_supersession_audit.csv = c(
    cox_catalog_all_sources, cox_run_evidence_sources, cox_valid_completion_sources
  ),
  multiple_imputation_audit.csv = c(markov_mi_sources, markov_fit_sources, cox_run_sources),
  protocol_sensitivity_specifications.csv = protocol_spec_sources,
  protocol_sensitivity_model_registry.csv = protocol_registry_sources,
  protocol_sensitivity_estimates.csv = protocol_estimate_sources,
  protocol_sensitivity_joint_tests.csv = protocol_joint_sources,
  protocol_sensitivity_coverage.csv = protocol_coverage_sources,
  ipcw_weight_diagnostics.csv = protocol_ipcw_sources,
  diagnosis_denial_pattern_audit.csv = protocol_denial_sources,
  expanded_hypertension_state_audit.csv = protocol_expanded_sources,
  protocol_sensitivity_data_audit.csv = protocol_data_audit_sources
)

for (name in names(canonical)) {
  value <- strip_internal_columns(as.data.table(canonical[[name]]))
  path <- file.path(result_root, name)
  atomic_fwrite(value, path)
  remember_sources(name, paths = source_map[[name]])
}

coverage_path <- file.path(report_root, "analysis_coverage_matrix.csv")
atomic_fwrite(strip_internal_columns(coverage_matrix), coverage_path)
remember_sources("analysis_coverage_matrix.csv", values = list(model_registry), paths = source_map[["model_specification_registry.csv"]])

manifest_targets <- c(
  file.path(result_root, names(canonical)),
  coverage_path
)
manifest <- rbindlist(lapply(manifest_targets, function(path) {
  info <- file.info(path)
  name <- basename(path)
  value <- tryCatch(fread(path, showProgress = FALSE), error = function(e) data.table())
  sources <- artifact_sources[[name]]
  data.table(
    artifact_id = tools::file_path_sans_ext(name),
    artifact_role = if (identical(path, coverage_path)) "coverage_report" else "canonical_result",
    path = relative_to_output(path),
    exists = file.exists(path),
    nonempty = is.finite(info$size) && info$size > 0L,
    rows = nrow(value),
    columns = ncol(value),
    size_bytes = as.numeric(info$size),
    md5 = if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_,
    source_file_count = length(sources),
    source_files = if (length(sources)) paste(vapply(sources, relative_to_output, character(1L)), collapse = " | ") else "",
    generated_at = generated_at,
    script_version = script_version,
    estimate_policy = "Existing numerical estimates only; absent outputs use status-only rows with missing estimate fields"
  )
}), use.names = TRUE, fill = TRUE)

manifest_path <- file.path(result_root, "results_manifest.csv")
atomic_fwrite(manifest, manifest_path)

cat("Canonical multistate output collation completed.\n")
cat("Output root:", normalizePath(output_root, winslash = "/", mustWork = FALSE), "\n")
cat("Registry rows:", nrow(model_registry), "\n")
cat("Completed registry rows:", model_registry[status == "completed", .N], "\n")
cat("Not-estimable registry rows:", model_registry[status == "not_estimable", .N], "\n")
cat("Failed registry rows:", model_registry[status == "failed", .N], "\n")
cat("Validated completed Cox run directories:", length(cox_completed_run_directories), "\n")
cat("Protocol sensitivity specifications:", nrow(protocol_registry), "\n")
cat("Protocol sensitivity cohort-family coverage rows:", nrow(protocol_coverage), "\n")
cat("Manifest:", normalizePath(manifest_path, winslash = "/", mustWork = FALSE), "\n")
