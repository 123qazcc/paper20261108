#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

options(stringsAsFactors = FALSE, width = 180)

stage_root <- Sys.getenv("GA_STAGE_ROOT", unset = "private/stage")
output_root <- Sys.getenv(
  "GA_MULTISTATE_OUTPUT_ROOT",
  unset = file.path(stage_root, "multistate_output")
)
data_root <- file.path(output_root, "data_multistate")
result_root <- file.path(output_root, "results")
figure_root <- file.path(output_root, "figures")
qc_root <- file.path(output_root, "qc")
report_root <- file.path(output_root, "reports")

input_directory_exists <- c(
  data_multistate = dir.exists(data_root),
  results = dir.exists(result_root),
  figures = dir.exists(figure_root),
  qc = dir.exists(qc_root)
)

# The verifier must always leave an audit trail, including when required upstream
# artefacts are absent.  Only the two output directories are created here.
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)
dir.create(report_root, recursive = TRUE, showWarnings = FALSE)

verification_rows <- list()
manifest_rows <- list()

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

normalise_name <- function(x) {
  gsub("(^_+|_+$)", "", gsub("[^a-z0-9]+", "_", tolower(x)))
}

normalise_status <- function(x) {
  y <- normalise_name(clean_text(x))
  y[y %chin% c("complete", "success", "successful", "done")] <- "completed"
  y[y %chin% c("notestimable", "not_estimateable", "not_estimatable", "inestimable")] <- "not_estimable"
  y[y %chin% c("failure", "error")] <- "failed"
  y
}

truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & x != 0)
  normalise_name(x) %chin% c("1", "true", "yes", "y", "prespecified", "pre_specified")
}

numeric_value <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  suppressWarnings(as.numeric(gsub("^[[:space:]]*[<=>~]+[[:space:]]*", "", clean_text(x))))
}

compact_detail <- function(...) {
  z <- unlist(list(...), recursive = TRUE, use.names = FALSE)
  z <- clean_text(z)
  z <- z[nzchar(z)]
  if (!length(z)) "" else paste(z, collapse = "; ")
}

relative_path <- function(path) {
  if (!length(path) || is.na(path) || !nzchar(path)) return(NA_character_)
  root_norm <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root_norm, "/")
  if (startsWith(tolower(path_norm), tolower(prefix))) {
    substring(path_norm, nchar(prefix) + 1L)
  } else {
    path_norm
  }
}

add_check <- function(
    domain,
    cohort,
    check_id,
    check,
    passed,
    observed = "",
    expected = "",
    detail = "",
    source_file = NA_character_) {
  ok <- isTRUE(passed)
  verification_rows[[length(verification_rows) + 1L]] <<- data.table(
    domain = as.character(domain),
    cohort = as.character(cohort),
    check_id = as.character(check_id),
    check = as.character(check),
    status = if (ok) "PASS" else "FAIL",
    observed = as.character(observed),
    expected = as.character(expected),
    detail = as.character(detail),
    source_file = ifelse(is.na(source_file), NA_character_, relative_path(source_file))
  )
  invisible(ok)
}

file_metadata <- function(path) {
  exists <- length(path) == 1L && !is.na(path) && file.exists(path)
  info <- if (exists) file.info(path) else NULL
  size <- if (exists) as.numeric(info$size) else NA_real_
  md5 <- if (exists && is.finite(size) && size > 0) {
    tryCatch(unname(tools::md5sum(path)), error = function(e) NA_character_)
  } else {
    NA_character_
  }
  list(
    exists = exists,
    nonempty = exists && is.finite(size) && size > 0,
    size = size,
    modified = if (exists) format(info$mtime, "%Y-%m-%dT%H:%M:%S%z") else NA_character_,
    md5 = md5
  )
}

register_artifact <- function(
    artifact_id,
    role,
    candidates,
    required = TRUE,
    cohort = "BOTH",
    require_nonempty = TRUE,
    check_file = TRUE) {
  candidates <- unique(as.character(candidates))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  selected <- if (length(existing)) existing[1L] else if (length(candidates)) candidates[1L] else NA_character_
  meta <- file_metadata(selected)
  passed <- meta$exists && (!require_nonempty || meta$nonempty)
  manifest_rows[[length(manifest_rows) + 1L]] <<- data.table(
    artifact_id = artifact_id,
    artifact_role = role,
    cohort = cohort,
    required = required,
    expected_path = if (length(candidates)) paste(vapply(candidates, relative_path, character(1L)), collapse = " | ") else NA_character_,
    actual_path = if (meta$exists) relative_path(selected) else NA_character_,
    exists = meta$exists,
    nonempty = meta$nonempty,
    size_bytes = meta$size,
    modified = meta$modified,
    md5 = meta$md5,
    verification_status = if (passed || !required) "PASS" else "FAIL",
    note = if (meta$exists) "Resolved" else "Required artefact not found"
  )
  if (check_file) {
    add_check(
      domain = "FILES",
      cohort = cohort,
      check_id = paste0("FILE_", toupper(normalise_name(artifact_id))),
      check = paste("Required artefact exists and is non-empty:", artifact_id),
      passed = if (required) passed else TRUE,
      observed = if (meta$exists) paste0("size_bytes=", meta$size) else "missing",
      expected = if (require_nonempty) "present; size_bytes > 0" else "present",
      source_file = selected
    )
  }
  if (passed) selected else NA_character_
}

find_col <- function(dt, candidates = character(), pattern = NULL) {
  if (is.null(dt) || !length(names(dt))) return(NA_character_)
  nm <- names(dt)
  nn <- normalise_name(nm)
  cand <- normalise_name(candidates)
  hit <- match(cand, nn, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) return(nm[hit[1L]])
  if (!is.null(pattern)) {
    idx <- grep(pattern, nn, perl = TRUE)
    if (length(idx)) return(nm[idx[1L]])
  }
  NA_character_
}

read_csv_artifact <- function(path, artifact_id, require_rows = TRUE) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  out <- tryCatch(
    fread(path, na.strings = c("", "NA", "N/A", "NaN", "."), showProgress = FALSE),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    add_check(
      "FILES", "BOTH", paste0("READ_", toupper(normalise_name(artifact_id))),
      paste("CSV is readable:", artifact_id), FALSE,
      observed = conditionMessage(out), expected = "readable CSV", source_file = path
    )
    return(NULL)
  }
  passed <- ncol(out) > 0L && (!require_rows || nrow(out) > 0L)
  add_check(
    "FILES", "BOTH", paste0("READ_", toupper(normalise_name(artifact_id))),
    paste("CSV has a valid non-empty schema:", artifact_id), passed,
    observed = paste0("rows=", nrow(out), ", columns=", ncol(out)),
    expected = if (require_rows) ">=1 row and >=1 column" else ">=1 column",
    source_file = path
  )
  out
}

read_data_artifact <- function(path, artifact_id, cohort) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  out <- tryCatch({
    if (grepl("\\.rds$", path, ignore.case = TRUE)) {
      as.data.table(readRDS(path))
    } else {
      fread(path, na.strings = c("", "NA", "N/A", "NaN", "."), showProgress = FALSE)
    }
  }, error = function(e) e)
  if (inherits(out, "error")) {
    add_check(
      "DATA", cohort, paste0("READ_", toupper(normalise_name(artifact_id))),
      paste("Analysis data are readable:", artifact_id), FALSE,
      observed = conditionMessage(out), expected = "readable RDS/CSV", source_file = path
    )
    return(NULL)
  }
  add_check(
    "DATA", cohort, paste0("READ_", toupper(normalise_name(artifact_id))),
    paste("Analysis data contain observations:", artifact_id),
    nrow(out) > 0L && ncol(out) > 0L,
    observed = paste0("rows=", nrow(out), ", columns=", ncol(out)),
    expected = ">=1 row and >=1 column", source_file = path
  )
  out
}

require_columns <- function(dt, cols, domain, cohort, check_id, source_file) {
  missing <- setdiff(cols, names(dt))
  add_check(
    domain, cohort, check_id,
    paste("Required fields are present:", paste(cols, collapse = ", ")),
    length(missing) == 0L,
    observed = if (length(missing)) paste("missing", paste(missing, collapse = ", ")) else "all present",
    expected = "all required fields", source_file = source_file
  )
  length(missing) == 0L
}

parse_transition_label <- function(x) {
  z <- toupper(gsub("[[:space:]]+", "", clean_text(x)))
  z <- gsub("STATE", "S", z, fixed = TRUE)
  z <- gsub("_TO_", "->", z, fixed = TRUE)
  z <- gsub("TO", "->", z, fixed = TRUE)
  pattern <- "S?([0-4])(?:→|->|–|—|>|-)S?([0-4])"
  matches <- regexec(pattern, z, perl = TRUE)
  pieces <- regmatches(z, matches)
  from <- rep(NA_integer_, length(z))
  to <- rep(NA_integer_, length(z))
  ok <- lengths(pieces) >= 3L
  if (any(ok)) {
    from[ok] <- as.integer(vapply(pieces[ok], `[[`, character(1L), 2L))
    to[ok] <- as.integer(vapply(pieces[ok], `[[`, character(1L), 3L))
  }
  data.table(from_state = from, to_state = to)
}

normalise_state_pair <- function(from, to, from_name = "", to_name = "") {
  from_chr <- clean_text(from)
  to_chr <- clean_text(to)
  labelled <- grepl("^[sS][0-4]$", from_chr) | grepl("^[sS][0-4]$", to_chr)
  from_num <- suppressWarnings(as.integer(gsub("[^0-9-]", "", from_chr)))
  to_num <- suppressWarnings(as.integer(gsub("[^0-9-]", "", to_chr)))
  observed <- c(from_num, to_num)
  use_msm_codes <- any(grepl("msm", normalise_name(c(from_name, to_name)))) ||
    (!any(labelled, na.rm = TRUE) && any(observed == 5L, na.rm = TRUE) && !any(observed == 0L, na.rm = TRUE))
  if (any(use_msm_codes)) {
    from_num <- from_num - 1L
    to_num <- to_num - 1L
  }
  data.table(from_state = from_num, to_state = to_num)
}

extract_edges <- function(dt) {
  if (is.null(dt) || !nrow(dt)) return(NULL)
  from_col <- find_col(dt, c("from_state", "state_from", "from", "origin_state", "origin", "from_state_msm"))
  to_col <- find_col(dt, c("to_state", "state_to", "to", "destination_state", "destination", "to_state_msm"))
  transition_col <- find_col(dt, c("transition", "transition_label", "edge", "transition_name"))
  if (!is.na(transition_col)) {
    parsed <- parse_transition_label(dt[[transition_col]])
  } else if (!is.na(from_col) && !is.na(to_col)) {
    parsed <- normalise_state_pair(dt[[from_col]], dt[[to_col]], from_col, to_col)
  } else {
    return(NULL)
  }
  out <- copy(dt)
  out[, `.from_state` := parsed$from_state]
  out[, `.to_state` := parsed$to_state]
  out[, `.edge` := fifelse(
    !is.na(.from_state) & !is.na(.to_state),
    paste0(.from_state, "->", .to_state), NA_character_
  )]
  out
}

row_text <- function(dt) {
  if (is.null(dt) || !nrow(dt)) return(character())
  character_cols <- names(dt)[vapply(dt, function(x) is.character(x) || is.factor(x), logical(1L))]
  if (!length(character_cols)) return(rep("", nrow(dt)))
  apply(dt[, ..character_cols], 1L, function(z) paste(clean_text(z), collapse = " | "))
}

status_mask <- function(dt) {
  status_col <- find_col(dt, c("model_status", "analysis_status", "fit_status", "status"))
  if (is.na(status_col)) return(rep(TRUE, nrow(dt)))
  status <- normalise_status(dt[[status_col]])
  if (!any(status %chin% terminal_statuses)) return(rep(TRUE, nrow(dt)))
  status == "completed"
}

reason_valid_for_terminal_status <- function(dt, domain, check_prefix, source_file) {
  status_col <- find_col(dt, c("status", "model_status", "analysis_status", "fit_status"))
  reason_col <- find_col(dt, c("reason", "status_reason", "failure_reason", "note", "notes", "explanation"))
  if (is.na(status_col)) {
    add_check(domain, "BOTH", paste0(check_prefix, "_STATUS_FIELD"),
              "Terminal status field is present", FALSE, "missing status field",
              "completed/not_estimable/failed", source_file = source_file)
    return(rep(FALSE, nrow(dt)))
  }
  status <- normalise_status(dt[[status_col]])
  allowed <- c("completed", "not_estimable", "failed")
  valid_status <- nzchar(status) & status %chin% allowed
  add_check(
    domain, "BOTH", paste0(check_prefix, "_STATUS_VALUES"),
    "All rows use an allowed terminal analysis status", all(valid_status),
    observed = paste(capture.output(print(table(status, useNA = "ifany"))), collapse = " "),
    expected = paste(allowed, collapse = ", "), source_file = source_file
  )
  need_reason <- status %chin% c("not_estimable", "failed")
  reason_ok <- if (is.na(reason_col)) !need_reason else !need_reason | nzchar(clean_text(dt[[reason_col]]))
  add_check(
    domain, "BOTH", paste0(check_prefix, "_NONCOMPLETED_REASON"),
    "Every not_estimable/failed row has a reason", all(reason_ok),
    observed = paste0("missing_reason_rows=", sum(!reason_ok)),
    expected = "0", source_file = source_file
  )
  valid_status & reason_ok
}

validate_optional_terminal_status <- function(dt, domain, check_prefix, source_file) {
  status_col <- find_col(dt, c("model_status", "analysis_status", "fit_status", "status"))
  if (is.na(status_col)) return(rep(TRUE, nrow(dt)))
  status <- normalise_status(dt[[status_col]])
  if (!any(status %chin% terminal_statuses)) return(rep(TRUE, nrow(dt)))
  reason_valid_for_terminal_status(dt, domain, check_prefix, source_file)
}

expected_edges <- c("0->1", "0->2", "0->4", "1->3", "1->4", "2->3", "2->4", "3->4")
admissible_panel_edges <- c(expected_edges, "0->3")
terminal_statuses <- c("completed", "not_estimable", "failed")
required_main_figures <- 1:7
q_network_expected <- list(
  hrs_5state_full = expected_edges,
  elsa_4state_w4_w10 = c("0->1", "0->2", "1->3", "2->3"),
  elsa_5state_through_2012 = expected_edges
)

q_network_contract <- function(dt) {
  if (is.null(dt) || !nrow(dt)) {
    return(list(ok = FALSE, detail = "empty input", selected_rows = 0L))
  }
  record_col <- find_col(dt, c("record_type"))
  if (!is.na(record_col)) {
    keep <- normalise_name(dt[[record_col]]) == "transition_intensity"
    dt <- dt[keep]
  }
  edges <- extract_edges(dt)
  if (is.null(edges) || !nrow(edges)) {
    return(list(ok = FALSE, detail = "transition rows unavailable", selected_rows = 0L))
  }
  analysis_col <- find_col(edges, c("analysis_id"))
  model_col <- find_col(edges, c("model_id"))
  profile_col <- find_col(edges, c("covariate_profile"))
  status_col <- find_col(edges, c("status"))
  pooling_col <- find_col(edges, c("pooling_status"))
  m_expected_col <- find_col(edges, c("m_expected"))
  m_fitted_col <- find_col(edges, c("m_fitted"))
  q_col <- find_col(edges, c("intensity_per_person_year", "intensity"))
  low_col <- find_col(edges, c("conf_low", "lower", "ci_low"))
  high_col <- find_col(edges, c("conf_high", "upper", "ci_high"))
  required <- c(
    analysis_col, model_col, profile_col, status_col, pooling_col,
    m_expected_col, m_fitted_col, q_col, low_col, high_col
  )
  if (any(is.na(required))) {
    return(list(ok = FALSE, detail = "one or more q-contract columns are missing", selected_rows = 0L))
  }
  analysis <- clean_text(edges[[analysis_col]])
  selected_mask <- analysis %chin% names(q_network_expected) &
    clean_text(edges[[model_col]]) == "M4_MI20" &
    clean_text(edges[[profile_col]]) == "sample means" &
    normalise_status(edges[[status_col]]) == "completed"
  selected <- copy(edges[selected_mask])
  selected[, `.analysis_verify` := clean_text(get(analysis_col))]
  q_value <- numeric_value(selected[[q_col]])
  q_low <- numeric_value(selected[[low_col]])
  q_high <- numeric_value(selected[[high_col]])
  pooled_ok <- nrow(selected) == 20L && all(
    clean_text(selected[[pooling_col]]) == "RUBIN_POOLED_COMPLETE" &
      numeric_value(selected[[m_expected_col]]) == 20 &
      numeric_value(selected[[m_fitted_col]]) == 20
  )
  numeric_ok <- nrow(selected) == 20L && all(
    is.finite(q_value) & q_value > 0 &
      is.finite(q_low) & q_low > 0 &
      is.finite(q_high) & q_high > 0 & q_low <= q_value & q_value <= q_high
  )
  edge_ok <- nrow(selected) == 20L && all(vapply(names(q_network_expected), function(id) {
    z <- selected[.analysis_verify == id]
    expected <- q_network_expected[[id]]
    nrow(z) == length(expected) && uniqueN(z$.edge) == length(expected) && setequal(z$.edge, expected)
  }, logical(1L)))
  duplicate_rows <- if (nrow(selected)) {
    selected[, .N, by = .(.analysis_verify, .edge)][N != 1L, .N]
  } else {
    0L
  }
  forbidden_direct <- sum(selected$.edge == "0->3", na.rm = TRUE)
  ok <- pooled_ok && numeric_ok && edge_ok && duplicate_rows == 0L && forbidden_direct == 0L
  list(
    ok = ok,
    selected_rows = nrow(selected),
    detail = paste0(
      "selected_rows=", nrow(selected), ", pooled20_ok=", pooled_ok,
      ", numeric_ci_ok=", numeric_ok, ", edge_sets_ok=", edge_ok,
      ", duplicate_keys=", duplicate_rows, ", direct_S0_to_S3=", forbidden_direct
    )
  )
}

# -----------------------------------------------------------------------------
# Directory and required-file contract
# -----------------------------------------------------------------------------

for (nm in names(input_directory_exists)) {
  path <- list(
    data_multistate = data_root,
    results = result_root,
    figures = figure_root,
    qc = qc_root
  )[[nm]]
  add_check(
    "FILES", "BOTH", paste0("DIR_", toupper(nm)),
    paste("Required output directory existed before verification:", nm), isTRUE(input_directory_exists[[nm]]),
    observed = if (isTRUE(input_directory_exists[[nm]])) "present" else "missing before verifier created its audit directories",
    expected = "present before verification",
    source_file = path
  )
}

data_paths <- list()
for (cohort in c("hrs", "elsa")) {
  cohort_label <- toupper(cohort)
  data_paths[[paste0(cohort, "_panel")]] <- register_artifact(
    paste0(cohort, "_msm_panel"), "analysis_data",
    c(
      file.path(data_root, paste0(cohort, "_msm_panel.rds")),
      file.path(data_root, paste0(cohort, "_msm_panel.csv.gz")),
      file.path(data_root, paste0(cohort, "_msm_panel.csv"))
    ), cohort = cohort_label
  )
  data_paths[[paste0(cohort, "_intervals")]] <- register_artifact(
    paste0(cohort, "_multistate_intervals"), "analysis_data",
    c(
      file.path(data_root, paste0(cohort, "_multistate_intervals.rds")),
      file.path(data_root, paste0(cohort, "_multistate_intervals.csv.gz")),
      file.path(data_root, paste0(cohort, "_multistate_intervals.csv"))
    ), cohort = cohort_label
  )
}

canonical_files <- list(
  transition_catalog = file.path(result_root, "transition_catalog.csv"),
  transition_counts = file.path(result_root, "transition_counts_rates.csv"),
  observed_panel_counts = file.path(result_root, "observed_panel_transition_counts.csv"),
  model_registry = file.path(result_root, "model_specification_registry.csv"),
  markov_intensity = file.path(result_root, "markov_transition_intensity_estimates.csv"),
  markov_q = file.path(result_root, "markov_q_intensities.csv"),
  markov_probability = file.path(result_root, "markov_transition_probability_matrices.csv"),
  markov_sojourn = file.path(result_root, "markov_sojourn_times.csv"),
  markov_diagnostics = file.path(result_root, "markov_model_diagnostics.csv"),
  markov_occupation = file.path(result_root, "markov_state_occupation_predictions.csv"),
  cox_estimates = file.path(result_root, "cox_transition_estimates.csv"),
  cox_diagnostics = file.path(result_root, "cox_model_diagnostics.csv"),
  mi_audit = file.path(result_root, "multiple_imputation_audit.csv"),
  coverage_matrix = file.path(report_root, "analysis_coverage_matrix.csv"),
  figure_qa = file.path(qc_root, "figure_qa.csv")
)

# Compatible aliases are accepted, but the frozen canonical name remains first.
aliases <- list(
  model_registry = c("analysis_specification_registry.csv", "analysis_registry.csv", "model_registry.csv"),
  markov_intensity = c("markov_intensity_estimates.csv", "msm_transition_intensity_estimates.csv", "qmatrix_estimates.csv"),
  markov_q = c("markov_absolute_q_intensities.csv"),
  markov_probability = c("markov_transition_probabilities.csv", "msm_transition_probability_matrices.csv", "pmatrix_estimates.csv"),
  markov_sojourn = c("msm_sojourn_times.csv", "sojourn_times.csv"),
  markov_diagnostics = c("msm_model_diagnostics.csv", "markov_diagnostics.csv"),
  markov_occupation = c("state_occupation_predictions.csv", "msm_state_occupation_predictions.csv"),
  cox_estimates = c("transition_cox_estimates.csv", "cox_transition_results.csv"),
  cox_diagnostics = c("transition_cox_diagnostics.csv", "cox_diagnostics.csv"),
  mi_audit = c("mi_audit.csv", "multiple_imputation_status.csv"),
  coverage_matrix = c("analysis_specification_coverage.csv", "analysis_coverage.csv"),
  figure_qa = c("figure_qa_audit.csv", "figures_qa.csv")
)

resolved <- list()
tables <- list()
for (artifact_id in names(canonical_files)) {
  primary <- canonical_files[[artifact_id]]
  parent <- dirname(primary)
  alias_names <- aliases[[artifact_id]]
  if (is.null(alias_names)) alias_names <- character()
  candidates <- c(primary, file.path(parent, alias_names))
  resolved[[artifact_id]] <- register_artifact(
    artifact_id, if (artifact_id %chin% c("coverage_matrix")) "report" else if (artifact_id == "figure_qa") "quality_audit" else "result",
    candidates, required = TRUE, cohort = "BOTH"
  )
  tables[[artifact_id]] <- read_csv_artifact(resolved[[artifact_id]], artifact_id, require_rows = TRUE)
}

# -----------------------------------------------------------------------------
# Cohort panel and interval integrity
# -----------------------------------------------------------------------------

panel_objects <- list()
interval_objects <- list()
forward_fields <- c(
  "insomnia_score_primary", "insomnia_burden_code", "insomnia_high_primary",
  "insomnia_score_cumavg", "insomnia_score_prev_valid",
  "age_bl", "female_bl", "education_level_bl", "wealth_quintile_bl", "married_bl",
  "current_smoker_bl", "current_drinker_bl", "inactive_bl", "bmi_bl",
  "diabetes_bl", "lung_disease_bl", "frequent_pain_bl",
  "bmi_tv", "current_smoker_tv", "current_drinker_tv", "inactive_tv",
  "diabetes_tv", "lung_disease_tv", "frequent_pain_tv", "depression_score_tv"
)

values_equal <- function(a, b) {
  if (length(a) != 1L || length(b) != 1L || is.na(a) || is.na(b)) return(FALSE)
  if (is.numeric(a) || is.integer(a) || is.numeric(b) || is.integer(b)) {
    aa <- suppressWarnings(as.numeric(a))
    bb <- suppressWarnings(as.numeric(b))
    is.finite(aa) && is.finite(bb) && abs(aa - bb) <= 1e-8 * max(1, abs(aa), abs(bb))
  } else {
    identical(as.character(a), as.character(b))
  }
}

for (cohort in c("hrs", "elsa")) {
  cohort_label <- toupper(cohort)
  panel_path <- data_paths[[paste0(cohort, "_panel")]]
  interval_path <- data_paths[[paste0(cohort, "_intervals")]]
  panel <- read_data_artifact(panel_path, paste0(cohort, "_msm_panel"), cohort_label)
  intervals <- read_data_artifact(interval_path, paste0(cohort, "_multistate_intervals"), cohort_label)
  panel_objects[[cohort]] <- panel
  interval_objects[[cohort]] <- intervals

  if (!is.null(panel)) {
    core_ok <- require_columns(
      panel, c("id", "time_msm", "state_msm"), "DATA", cohort_label,
      paste0(cohort_label, "_PANEL_CORE_FIELDS"), panel_path
    )
    if (core_ok) {
      bad_id_time <- panel[is.na(id) | !is.finite(as.numeric(time_msm)), .N]
      duplicate_id_time <- nrow(panel) - uniqueN(panel, by = c("id", "time_msm"))
      add_check(
        "DATA", cohort_label, paste0(cohort_label, "_ID_TIME_UNIQUE"),
        "Participant-time rows are unique and finite",
        bad_id_time == 0L && duplicate_id_time == 0L,
        observed = paste0("missing_or_nonfinite=", bad_id_time, ", duplicates=", duplicate_id_time),
        expected = "0; 0", source_file = panel_path
      )
      ordered <- copy(panel)
      setorderv(ordered, c("id", "time_msm"))
      timing <- ordered[, .(
        negative_time = any(as.numeric(time_msm) < 0),
        non_increasing = .N > 1L && any(diff(as.numeric(time_msm)) <= 0)
      ), by = id]
      bad_timing <- timing[negative_time | non_increasing, .N]
      add_check(
        "DATA", cohort_label, paste0(cohort_label, "_TIME_ORDER"),
        "Within-person MSM time is nonnegative and strictly increasing",
        bad_timing == 0L,
        observed = paste0("participants_with_violation=", bad_timing), expected = "0",
        source_file = panel_path
      )
      bad_state <- panel[is.na(state_msm) | !as.integer(state_msm) %in% 1:5, .N]
      add_check(
        "STATE", cohort_label, paste0(cohort_label, "_STATE_CODES"),
        "MSM state codes are exactly within 1-5 (S0-S4)", bad_state == 0L,
        observed = paste0("violating_rows=", bad_state), expected = "0",
        source_file = panel_path
      )

      ordered[, `.row_in_id` := seq_len(.N), by = id]
      ordered[, `.n_in_id` := .N, by = id]
      death_rows <- ordered[as.integer(state_msm) == 5L]
      death_not_final <- death_rows[.row_in_id != .n_in_id, .N]
      repeated_death <- death_rows[, .N, by = id][N > 1L, .N]
      add_check(
        "DEATH", cohort_label, paste0(cohort_label, "_DEATH_FINAL"),
        "Death is a single final record within participant",
        death_not_final == 0L && repeated_death == 0L,
        observed = paste0("nonfinal_death_rows=", death_not_final, ", participants_with_repeated_death=", repeated_death),
        expected = "0; 0", source_file = panel_path
      )

      field_presence <- intersect(forward_fields, names(ordered))
      missing_forward_fields <- setdiff(forward_fields, names(ordered))
      add_check(
        "DEATH", cohort_label, paste0(cohort_label, "_FORWARD_FIELDS_PRESENT"),
        "Required interval-start exposure/covariate fields exist on the MSM panel",
        length(missing_forward_fields) == 0L,
        observed = if (length(missing_forward_fields)) paste(missing_forward_fields, collapse = ", ") else "all present",
        expected = paste(forward_fields, collapse = ", "), source_file = panel_path
      )
      primary_fields <- intersect(c("insomnia_score_primary", "insomnia_burden_code"), names(ordered))
      death_primary_missing <- if (length(primary_fields) == 2L) {
        death_rows[, sum(!complete.cases(.SD)), .SDcols = primary_fields]
      } else {
        nrow(death_rows)
      }
      add_check(
        "DEATH", cohort_label, paste0(cohort_label, "_DEATH_EXPOSURE_AVAILABLE"),
        "Primary exposure is available on every death row",
        length(primary_fields) == 2L && death_primary_missing == 0L,
        observed = paste0("death_rows=", nrow(death_rows), ", rows_missing_primary_exposure=", death_primary_missing),
        expected = "0 missing", source_file = panel_path
      )

      # A missing value is allowed only when that field was never observed before death.
      # Whenever a prior nonmissing value exists, the death-row value must equal the
      # most recent prior value (last observation carried forward).
      forward_violations <- 0L
      forward_comparisons <- 0L
      if (nrow(death_rows) && length(field_presence)) {
        death_ids <- unique(death_rows$id)
        for (id_value in death_ids) {
          subject <- ordered[id == id_value]
          death_index <- which(as.integer(subject$state_msm) == 5L)[1L]
          if (is.na(death_index) || death_index <= 1L) next
          prior <- subject[seq_len(death_index - 1L)]
          death <- subject[death_index]
          for (field in field_presence) {
            prior_values <- prior[[field]]
            prior_nonmissing <- which(!is.na(prior_values))
            if (!length(prior_nonmissing)) next
            forward_comparisons <- forward_comparisons + 1L
            prior_value <- prior_values[tail(prior_nonmissing, 1L)]
            death_value <- death[[field]][1L]
            if (!values_equal(prior_value, death_value)) forward_violations <- forward_violations + 1L
          }
        }
      }
      add_check(
        "DEATH", cohort_label, paste0(cohort_label, "_DEATH_FORWARD_FILL"),
        "Death-row exposures/covariates equal the latest prior observed value",
        forward_violations == 0L,
        observed = paste0("comparisons=", forward_comparisons, ", violations=", forward_violations),
        expected = "0 violations", source_file = panel_path
      )

      # Attained age is intentionally recomputed at each event time, including
      # death, and therefore must not be tested as a last-observation-carried-
      # forward covariate.  Check its temporal direction separately.
      age_decrease_violations <- 0L
      age_comparisons <- 0L
      if (nrow(death_rows) && "age_tv" %in% names(ordered)) {
        for (id_value in unique(death_rows$id)) {
          subject <- ordered[id == id_value]
          death_index <- which(as.integer(subject$state_msm) == 5L)[1L]
          if (is.na(death_index) || death_index <= 1L) next
          prior_age <- subject[["age_tv"]][seq_len(death_index - 1L)]
          prior_nonmissing <- which(!is.na(prior_age))
          death_age <- subject[["age_tv"]][death_index]
          if (!length(prior_nonmissing) || is.na(death_age)) next
          age_comparisons <- age_comparisons + 1L
          last_age <- as.numeric(prior_age[tail(prior_nonmissing, 1L)])
          if (!is.finite(last_age) || !is.finite(as.numeric(death_age)) || as.numeric(death_age) < last_age) {
            age_decrease_violations <- age_decrease_violations + 1L
          }
        }
      }
      add_check(
        "DEATH", cohort_label, paste0(cohort_label, "_DEATH_ATTAINED_AGE"),
        "Attained age at death is no lower than the last observed living age",
        age_decrease_violations == 0L,
        observed = paste0("comparisons=", age_comparisons, ", violations=", age_decrease_violations),
        expected = "0 violations", source_file = panel_path
      )
    }
  }

  if (!is.null(intervals)) {
    core_ok <- require_columns(
      intervals,
      c("id", "tstart", "tstop", "interval_years", "from_state", "to_state"),
      "DATA", cohort_label, paste0(cohort_label, "_INTERVAL_CORE_FIELDS"), interval_path
    )
    if (core_ok) {
      duplicate_intervals <- nrow(intervals) - uniqueN(intervals, by = c("id", "tstart", "tstop"))
      tstart <- as.numeric(intervals$tstart)
      tstop <- as.numeric(intervals$tstop)
      duration <- as.numeric(intervals$interval_years)
      bad_duration <- sum(
        is.na(tstart) | is.na(tstop) | is.na(duration) |
          !is.finite(tstart) | !is.finite(tstop) | !is.finite(duration) |
          tstop <= tstart | duration <= 0,
        na.rm = TRUE
      )
      duration_mismatch <- sum(
        is.finite(tstart) & is.finite(tstop) & is.finite(duration) &
          abs((tstop - tstart) - duration) > 1e-8 * pmax(1, abs(duration)),
        na.rm = TRUE
      )
      add_check(
        "DATA", cohort_label, paste0(cohort_label, "_INTERVAL_UNIQUE_DURATION"),
        "Intervals are participant-time unique and have positive, internally consistent duration",
        duplicate_intervals == 0L && bad_duration == 0L && duration_mismatch == 0L,
        observed = paste0(
          "duplicates=", duplicate_intervals, ", nonpositive_or_nonfinite=", bad_duration,
          ", duration_mismatch=", duration_mismatch
        ), expected = "0; 0; 0", source_file = interval_path
      )

      interval_exposure_fields <- c("insomnia_score_primary", "insomnia_burden_code")
      missing_interval_exposure_fields <- setdiff(interval_exposure_fields, names(intervals))
      missing_interval_exposure <- if (!length(missing_interval_exposure_fields)) {
        intervals[, sum(!complete.cases(.SD)), .SDcols = interval_exposure_fields]
      } else {
        nrow(intervals)
      }
      add_check(
        "DATA", cohort_label, paste0(cohort_label, "_INTERVAL_START_EXPOSURE"),
        "Every risk interval has the latest valid primary exposure at its start",
        length(missing_interval_exposure_fields) == 0L && missing_interval_exposure == 0L,
        observed = compact_detail(
          paste0("rows_missing_exposure=", missing_interval_exposure),
          if (length(missing_interval_exposure_fields)) {
            paste0("missing_fields=", paste(missing_interval_exposure_fields, collapse = ","))
          } else {
            ""
          }
        ),
        expected = "0 missing interval-start exposure rows", source_file = interval_path
      )

      intervals[, `.edge_verify` := paste0(as.integer(from_state), "->", as.integer(to_state))]
      changed <- as.integer(intervals$from_state) != as.integer(intervals$to_state)
      bad_observed_edge <- sum(changed & !intervals$.edge_verify %chin% admissible_panel_edges, na.rm = TRUE)
      add_check(
        "STATE", cohort_label, paste0(cohort_label, "_PANEL_EDGE_SET"),
        "Observed changes use only eight model edges or the S0-to-S3 panel multistep jump",
        bad_observed_edge == 0L,
        observed = paste0("forbidden_changed_intervals=", bad_observed_edge), expected = "0",
        source_file = interval_path
      )

      jump <- intervals[.edge_verify == "0->3" & changed]
      jump_flag_col <- find_col(intervals, c("interval_censored_multistep_jump", "panel_multistep_jump"))
      intensity_flag_col <- find_col(intervals, c("allowed_intensity_edge", "model_intensity_edge"))
      jump_flags_ok <- !is.na(jump_flag_col) && !is.na(intensity_flag_col) &&
        (nrow(jump) == 0L || (
          all(as.integer(jump[[jump_flag_col]]) == 1L, na.rm = FALSE) &&
            all(as.integer(jump[[intensity_flag_col]]) == 0L, na.rm = FALSE)
        ))
      add_check(
        "STATE", cohort_label, paste0(cohort_label, "_S0_S3_PANEL_ONLY"),
        "S0-to-S3 is flagged only as a panel-observed multistep jump, never a main intensity",
        jump_flags_ok,
        observed = paste0(
          "jump_rows=", nrow(jump), ", multistep_field=", ifelse(is.na(jump_flag_col), "missing", jump_flag_col),
          ", intensity_field=", ifelse(is.na(intensity_flag_col), "missing", intensity_flag_col)
        ),
        expected = "multistep=1 and allowed_intensity=0 for every S0->S3 row",
        source_file = interval_path
      )
      intervals[, `.edge_verify` := NULL]
    }
  }
}

# -----------------------------------------------------------------------------
# Frozen transition graph and empirical event counts
# -----------------------------------------------------------------------------

transition_catalog <- tables$transition_catalog
if (!is.null(transition_catalog)) {
  cat_edges <- extract_edges(transition_catalog)
  if (is.null(cat_edges)) {
    add_check("STATE", "BOTH", "TRANSITION_CATALOG_SCHEMA", "Transition catalog identifies from/to states",
              FALSE, "unrecognised edge fields", "from_state/to_state or transition", source_file = resolved$transition_catalog)
  } else {
    observed_edges <- sort(unique(cat_edges$.edge[!is.na(cat_edges$.edge)]))
    exact_graph <- setequal(observed_edges, expected_edges) &&
      sum(!is.na(cat_edges$.edge)) == length(expected_edges) && !anyDuplicated(cat_edges$.edge[!is.na(cat_edges$.edge)])
    add_check(
      "STATE", "BOTH", "TRANSITION_CATALOG_EIGHT_EDGES",
      "Transition catalog is exactly the frozen eight-edge Q graph",
      exact_graph,
      observed = paste(observed_edges, collapse = ", "),
      expected = paste(expected_edges, collapse = ", "),
      source_file = resolved$transition_catalog
    )
    add_check(
      "STATE", "BOTH", "TRANSITION_CATALOG_NO_S0_S3",
      "The main transition catalog excludes a direct S0-to-S3 intensity",
      !"0->3" %chin% observed_edges,
      observed = paste0("S0->S3 entries=", sum(cat_edges$.edge == "0->3", na.rm = TRUE)),
      expected = "0", source_file = resolved$transition_catalog
    )
  }
}

transition_counts <- tables$transition_counts
counts_standard <- NULL
if (!is.null(transition_counts)) {
  count_edges <- extract_edges(transition_counts)
  cohort_col <- find_col(transition_counts, c("cohort", "study"))
  event_col <- find_col(transition_counts, c("events", "event_count", "n_events"))
  risk_col <- find_col(transition_counts, c("risk_intervals", "at_risk_intervals", "n_at_risk", "risk_set"))
  people_col <- find_col(transition_counts, c("risk_participants", "participants_at_risk", "n_participants"))
  py_col <- find_col(transition_counts, c("person_years", "person_time", "personyears"))
  schema_ok <- !is.null(count_edges) && !any(is.na(c(cohort_col, event_col, risk_col, people_col, py_col)))
  add_check(
    "EVENTS", "BOTH", "TRANSITION_COUNTS_SCHEMA",
    "Transition-count table contains cohort, edge, events, risk set and person-time",
    schema_ok,
    observed = compact_detail(
      paste0("cohort=", cohort_col), paste0("events=", event_col), paste0("risk=", risk_col),
      paste0("participants=", people_col), paste0("person_years=", py_col)
    ),
    expected = "all fields present", source_file = resolved$transition_counts
  )
  if (schema_ok) {
    counts_standard <- data.table(
      cohort = toupper(clean_text(transition_counts[[cohort_col]])),
      edge = count_edges$.edge,
      events = numeric_value(transition_counts[[event_col]]),
      risk_intervals = numeric_value(transition_counts[[risk_col]]),
      risk_participants = numeric_value(transition_counts[[people_col]]),
      person_years = numeric_value(transition_counts[[py_col]])
    )
    required_grid <- CJ(cohort = c("ELSA", "HRS"), edge = expected_edges, sorted = TRUE)
    observed_grid <- unique(counts_standard[cohort %chin% c("HRS", "ELSA") & edge %chin% expected_edges, .(cohort, edge)])
    missing_grid <- fsetdiff(required_grid, observed_grid)
    duplicate_grid <- counts_standard[cohort %chin% c("HRS", "ELSA") & edge %chin% expected_edges, .N, by = .(cohort, edge)][N != 1L]
    add_check(
      "EVENTS", "BOTH", "TRANSITION_COUNTS_COMPLETE",
      "Every frozen edge has exactly one empirical count row in each cohort",
      nrow(missing_grid) == 0L && nrow(duplicate_grid) == 0L,
      observed = paste0("missing=", nrow(missing_grid), ", nonunique=", nrow(duplicate_grid)),
      expected = "0; 0", source_file = resolved$transition_counts
    )
    invalid_counts <- counts_standard[
      is.na(events) | events < 0 | events != floor(events) |
        is.na(risk_intervals) | risk_intervals < events | risk_intervals <= 0 |
        is.na(risk_participants) | risk_participants < 0 |
        is.na(person_years) | !is.finite(person_years) | person_years <= 0,
      .N
    ]
    add_check(
      "EVENTS", "BOTH", "TRANSITION_COUNTS_VALUES",
      "Event counts, risk sets and person-years are numerically valid",
      invalid_counts == 0L,
      observed = paste0("invalid_rows=", invalid_counts), expected = "0",
      source_file = resolved$transition_counts
    )
  }
}

observed_panel_counts <- tables$observed_panel_counts
if (!is.null(observed_panel_counts)) {
  observed_edges <- extract_edges(observed_panel_counts)
  if (is.null(observed_edges)) {
    add_check("STATE", "BOTH", "OBSERVED_PANEL_EDGE_SCHEMA", "Observed-panel count table identifies edges",
              FALSE, "unrecognised edge fields", "from_state/to_state or transition", source_file = resolved$observed_panel_counts)
  } else {
    main_flag_col <- find_col(observed_panel_counts, c("model_intensity_edge", "allowed_intensity_edge"))
    jump_flag_col <- find_col(observed_panel_counts, c("interval_censored_multistep_jump", "panel_multistep_jump"))
    jump_rows <- which(observed_edges$.edge == "0->3")
    flags_ok <- !is.na(main_flag_col) && !is.na(jump_flag_col) &&
      (!length(jump_rows) || (
        all(as.integer(observed_panel_counts[[main_flag_col]][jump_rows]) == 0L, na.rm = FALSE) &&
          all(as.integer(observed_panel_counts[[jump_flag_col]][jump_rows]) == 1L, na.rm = FALSE)
      ))
    add_check(
      "STATE", "BOTH", "OBSERVED_PANEL_S0_S3_CLASSIFICATION",
      "Observed-count output classifies S0-to-S3 as panel multistep and not a fitted edge",
      flags_ok,
      observed = paste0("S0->S3 rows=", length(jump_rows)),
      expected = "model_intensity_edge=0; interval_censored_multistep_jump=1",
      source_file = resolved$observed_panel_counts
    )
  }
}

# -----------------------------------------------------------------------------
# Specification registry, coverage matrix and ELSA mortality-scope labelling
# -----------------------------------------------------------------------------

registry <- tables$model_registry
registry_prespecified <- NULL
registry_id_col <- registry_cohort_col <- registry_status_col <- registry_reason_col <- NA_character_
if (!is.null(registry)) {
  registry_id_col <- find_col(registry, c("spec_id", "specification_id", "analysis_id", "model_id"))
  registry_cohort_col <- find_col(registry, c("cohort", "study"))
  registry_status_col <- find_col(registry, c("status", "analysis_status", "model_status", "fit_status"))
  registry_reason_col <- find_col(registry, c("reason", "status_reason", "failure_reason", "note", "notes"))
  prespecified_col <- find_col(registry, c("prespecified", "pre_specified", "is_prespecified"))
  schema_ok <- !any(is.na(c(registry_id_col, registry_cohort_col, registry_status_col, registry_reason_col)))
  add_check(
    "COVERAGE", "BOTH", "REGISTRY_SCHEMA",
    "Model registry identifies specification, cohort, terminal status and reason",
    schema_ok,
    observed = compact_detail(
      paste0("id=", registry_id_col), paste0("cohort=", registry_cohort_col),
      paste0("status=", registry_status_col), paste0("reason=", registry_reason_col),
      paste0("prespecified=", prespecified_col)
    ), expected = "spec_id, cohort, status, reason", source_file = resolved$model_registry
  )
  if (schema_ok) {
    registry_prespecified <- copy(registry)
    if (!is.na(prespecified_col)) registry_prespecified <- registry_prespecified[truthy(get(prespecified_col))]
    add_check(
      "COVERAGE", "BOTH", "REGISTRY_HAS_PRESET_ROWS",
      "Registry contains at least one prespecified analysis",
      nrow(registry_prespecified) > 0L,
      observed = paste0("prespecified_rows=", nrow(registry_prespecified)), expected = ">0",
      source_file = resolved$model_registry
    )
    reason_valid_for_terminal_status(registry_prespecified, "COVERAGE", "REGISTRY", resolved$model_registry)
    registry_prespecified[, `.spec_key` := paste(
      toupper(clean_text(get(registry_cohort_col))), clean_text(get(registry_id_col)), sep = "::"
    )]
    duplicate_specs <- registry_prespecified[, .N, by = .spec_key][N > 1L, .N]
    cohorts_present <- sort(unique(toupper(clean_text(registry_prespecified[[registry_cohort_col]]))))
    add_check(
      "COVERAGE", "BOTH", "REGISTRY_UNIQUE_SPECIFICATIONS",
      "Each prespecified cohort-specification has one registry row",
      duplicate_specs == 0L,
      observed = paste0("duplicate_keys=", duplicate_specs), expected = "0",
      source_file = resolved$model_registry
    )
    add_check(
      "COVERAGE", "BOTH", "REGISTRY_BOTH_COHORTS",
      "Registry covers both HRS and ELSA",
      all(c("HRS", "ELSA") %chin% cohorts_present),
      observed = paste(cohorts_present, collapse = ", "), expected = "HRS, ELSA",
      source_file = resolved$model_registry
    )
  }
}

coverage <- tables$coverage_matrix
if (!is.null(coverage)) {
  coverage_id_col <- find_col(coverage, c("spec_id", "specification_id", "analysis_id", "model_id"))
  coverage_cohort_col <- find_col(coverage, c("cohort", "study"))
  coverage_status_col <- find_col(coverage, c("status", "analysis_status", "model_status", "fit_status"))
  coverage_reason_col <- find_col(coverage, c("reason", "status_reason", "failure_reason", "note", "notes"))
  schema_ok <- !any(is.na(c(coverage_id_col, coverage_cohort_col, coverage_status_col, coverage_reason_col)))
  add_check(
    "COVERAGE", "BOTH", "COVERAGE_MATRIX_SCHEMA",
    "Coverage matrix identifies specification, cohort, terminal status and reason",
    schema_ok,
    observed = compact_detail(
      paste0("id=", coverage_id_col), paste0("cohort=", coverage_cohort_col),
      paste0("status=", coverage_status_col), paste0("reason=", coverage_reason_col)
    ), expected = "spec_id, cohort, status, reason", source_file = resolved$coverage_matrix
  )
  if (schema_ok) {
    reason_valid_for_terminal_status(coverage, "COVERAGE", "COVERAGE_MATRIX", resolved$coverage_matrix)
    coverage[, `.spec_key_verify` := paste(
      toupper(clean_text(get(coverage_cohort_col))), clean_text(get(coverage_id_col)), sep = "::"
    )]
    if (!is.null(registry_prespecified) && ".spec_key" %in% names(registry_prespecified)) {
      missing_coverage <- setdiff(unique(registry_prespecified$.spec_key), unique(coverage$.spec_key_verify))
      extra_coverage <- setdiff(unique(coverage$.spec_key_verify), unique(registry_prespecified$.spec_key))
      add_check(
        "COVERAGE", "BOTH", "COVERAGE_MATCHES_REGISTRY",
        "Coverage matrix contains every prespecified registry entry",
        length(missing_coverage) == 0L,
        observed = paste0("missing=", length(missing_coverage), ", extra=", length(extra_coverage)),
        expected = "missing=0", detail = if (length(missing_coverage)) paste(head(missing_coverage, 10L), collapse = ", ") else "",
        source_file = resolved$coverage_matrix
      )
    }
    coverage[, `.spec_key_verify` := NULL]
  }
}

mortality_sources <- Filter(Negate(is.null), list(registry, coverage, tables$markov_diagnostics, tables$markov_intensity))
mortality_text <- unlist(lapply(mortality_sources, function(dt) {
  cohort_col <- find_col(dt, c("cohort", "study"))
  text <- row_text(dt)
  if (!is.na(cohort_col)) {
    text[toupper(clean_text(dt[[cohort_col]])) == "ELSA"]
  } else {
    text[grepl("ELSA", text, ignore.case = TRUE)]
  }
}), use.names = FALSE)
elsa_limited_label <- grepl(
  "limited[ _-]*mortality[ _-]*validation|mortality.{0,40}(through|to|until|ends?[ _-]*in)[ _-]*2012|death.{0,40}2012|2012.{0,40}(mortality|death)|limited.{0,30}(mortality|death)|截止.{0,10}2012|有限.{0,20}死亡",
  mortality_text, ignore.case = TRUE, perl = TRUE
)
add_check(
  "DEATH", "ELSA", "ELSA_MORTALITY_COVERAGE_LABEL",
  "ELSA death analyses are explicitly labelled as limited mortality validation through 2012",
  any(elsa_limited_label),
  observed = paste0("matching_rows=", sum(elsa_limited_label)),
  expected = "at least one explicit ELSA limited-mortality/2012 label in registry, coverage or Markov outputs",
  source_file = resolved$model_registry
)

# -----------------------------------------------------------------------------
# Multiple-imputation audit (m = 20)
# -----------------------------------------------------------------------------

mi_audit <- tables$mi_audit
if (!is.null(mi_audit)) {
  m_requested_col <- find_col(
    mi_audit,
    c("m", "mi_m", "n_imputations", "m_requested", "imputations_requested", "number_of_imputations")
  )
  m_completed_col <- find_col(mi_audit, c("m_completed", "imputations_completed", "n_imputations_completed"))
  status_col <- find_col(mi_audit, c("status", "analysis_status", "model_status", "fit_status"))
  marker_ok <- !is.na(m_requested_col)
  requested_values <- if (marker_ok) numeric_value(mi_audit[[m_requested_col]]) else numeric()
  requested_ok <- marker_ok && length(requested_values) > 0L && all(is.finite(requested_values) & requested_values == 20)
  add_check(
    "MI", "BOTH", "MI_M_EQUALS_20",
    "Multiple-imputation audit explicitly records m = 20",
    requested_ok,
    observed = if (marker_ok) paste(sort(unique(requested_values)), collapse = ", ") else "m field missing",
    expected = "20", source_file = resolved$mi_audit
  )
  if (!is.na(m_completed_col)) {
    completed_mask <- if (!is.na(status_col)) normalise_status(mi_audit[[status_col]]) == "completed" else rep(TRUE, nrow(mi_audit))
    m_completed <- numeric_value(mi_audit[[m_completed_col]])
    completed_ok <- !any(completed_mask & (!is.finite(m_completed) | m_completed != 20))
    add_check(
      "MI", "BOTH", "MI_COMPLETED_M_EQUALS_20",
      "Every completed imputed analysis used all 20 imputations",
      completed_ok,
      observed = paste0("completed_rows_not_using_20=", sum(completed_mask & (!is.finite(m_completed) | m_completed != 20))),
      expected = "0", source_file = resolved$mi_audit
    )
  }
}

# -----------------------------------------------------------------------------
# Markov Q/intensity structure and model diagnostics
# -----------------------------------------------------------------------------

markov_intensity <- tables$markov_intensity
if (!is.null(markov_intensity)) {
  validate_optional_terminal_status(markov_intensity, "MARKOV", "MARKOV_INTENSITY", resolved$markov_intensity)
  intensity_edges <- extract_edges(markov_intensity)
  if (is.null(intensity_edges)) {
    add_check("MARKOV", "BOTH", "MARKOV_INTENSITY_EDGE_SCHEMA", "Markov intensity output identifies transitions",
              FALSE, "unrecognised edge fields", "from_state/to_state or transition", source_file = resolved$markov_intensity)
  } else {
    completed <- status_mask(intensity_edges)
    text <- row_text(intensity_edges)
    alternative_q <- grepl(
      "alternative.{0,20}q|alternative_q|direct.{0,20}s0.{0,10}s3|s0.{0,10}s3.{0,20}(open|direct)|nine.{0,10}edge|9.{0,5}edge",
      text, ignore.case = TRUE, perl = TRUE
    )
    main_rows <- completed & !alternative_q & !is.na(intensity_edges$.edge)
    forbidden <- intensity_edges$.edge[main_rows & !intensity_edges$.edge %chin% expected_edges &
                                         intensity_edges$.from_state != intensity_edges$.to_state]
    direct_main <- sum(main_rows & intensity_edges$.edge == "0->3", na.rm = TRUE)
    add_check(
      "MARKOV", "BOTH", "MARKOV_MAIN_EDGE_STRUCTURE",
      "Completed main Markov estimates use only the frozen eight intensity edges",
      length(forbidden) == 0L && direct_main == 0L,
      observed = paste0("forbidden_rows=", length(forbidden), ", main_S0->S3_rows=", direct_main),
      expected = "0; 0", detail = paste(unique(forbidden), collapse = ", "),
      source_file = resolved$markov_intensity
    )
    intensity_group_candidates <- c(
      "cohort", "spec_id", "specification_id", "analysis_id", "model_id", "model",
      "exposure", "imputation", "weighted"
    )
    intensity_group_cols <- names(intensity_edges)[
      normalise_name(names(intensity_edges)) %chin% normalise_name(intensity_group_candidates)
    ]
    if (!length(intensity_group_cols)) {
      intensity_edges[, `.intensity_group_verify` := "all"]
    } else {
      intensity_edges[, `.intensity_group_verify` := do.call(paste, c(.SD, sep = "::")),
                      .SDcols = intensity_group_cols]
    }
    terminal_by_intensity_group <- intensity_edges[
      completed & !is.na(.from_state) & !is.na(.to_state),
      .(terminal_state = max(c(.from_state, .to_state))),
      by = .intensity_group_verify
    ]
    intensity_with_terminal <- merge(
      intensity_edges, terminal_by_intensity_group,
      by = ".intensity_group_verify", all.x = TRUE, sort = FALSE
    )
    terminal_outgoing <- sum(
      status_mask(intensity_with_terminal) &
        is.finite(intensity_with_terminal$terminal_state) &
        intensity_with_terminal$.from_state == intensity_with_terminal$terminal_state &
        intensity_with_terminal$.to_state != intensity_with_terminal$terminal_state,
      na.rm = TRUE
    )
    add_check(
      "MARKOV", "BOTH", "MARKOV_DEATH_NO_OUTGOING_INTENSITY",
      "Each fitted model's terminal absorbing state has no outgoing estimated intensity",
      terminal_outgoing == 0L,
      observed = paste0("terminal_outgoing_rows=", terminal_outgoing), expected = "0",
      source_file = resolved$markov_intensity
    )
  }
}

markov_q <- tables$markov_q
if (!is.null(markov_q)) {
  q_contract <- q_network_contract(markov_q)
  add_check(
    "MARKOV", "BOTH", "MARKOV_Q_NETWORK_8_4_8_CONTRACT",
    paste(
      "Canonical absolute-q output contains the exact HRS 8-edge, ELSA waves 4-10",
      "4-edge, and ELSA through-2012 8-edge M4-MI20 networks"
    ),
    q_contract$ok,
    observed = q_contract$detail,
    expected = paste(
      "20 unique positive q rows with valid 95% CIs; Rubin pooled 20/20;",
      "edge counts 8/4/8; no direct S0->S3"
    ),
    source_file = resolved$markov_q
  )
}

markov_diag <- tables$markov_diagnostics
if (!is.null(markov_diag)) {
  validate_optional_terminal_status(markov_diag, "MARKOV", "MARKOV_DIAGNOSTICS", resolved$markov_diagnostics)
  completed <- status_mask(markov_diag)
  row_sum_col <- find_col(
    markov_diag,
    c("q_row_sum_max_abs", "q_max_abs_row_sum", "max_abs_q_row_sum", "max_q_row_sum_error"),
    pattern = "q.*row.*sum|row.*sum.*q"
  )
  # The ELSA four-state wave 4-10 nonfatal model uses S3 as its terminal
  # absorbing state; HRS and the ELSA through-2012 mortality model use S4.
  # Prefer the generic diagnostic and retain the death-specific field as a
  # backward-compatible row-wise fallback.
  absorbing_metric_cols <- unique(c(
    intersect(
      c("absorbing_row_max_abs", "q_absorbing_row_max_abs", "death_row_max_abs"),
      names(markov_diag)
    ),
    grep(
      "(death|absorbing).*row.*(max|abs|error)",
      names(markov_diag), value = TRUE, ignore.case = TRUE, perl = TRUE
    )
  ))
  absorbing_ok_col <- find_col(
    markov_diag,
    c("death_absorbing", "absorbing_state_ok", "death_state_absorbing", "q_absorbing_ok"),
    pattern = "(death|absorbing).*(ok|valid|pass)"
  )
  if (!any(completed)) {
    q_rows_ok <- TRUE
    absorb_ok <- TRUE
    q_detail <- "No completed Markov models; terminal registry reasons govern"
    absorb_detail <- q_detail
  } else {
    q_values <- if (!is.na(row_sum_col)) numeric_value(markov_diag[[row_sum_col]]) else numeric()
    q_rows_ok <- !is.na(row_sum_col) && !any(completed & (!is.finite(q_values) | abs(q_values) > 1e-6))
    q_detail <- if (is.na(row_sum_col)) "Q row-sum metric missing" else paste0("max_abs=", max(abs(q_values[completed]), na.rm = TRUE))
    if (length(absorbing_metric_cols)) {
      a_values <- rep(NA_real_, nrow(markov_diag))
      for (metric_col in absorbing_metric_cols) {
        candidate <- numeric_value(markov_diag[[metric_col]])
        fill <- !is.finite(a_values) & is.finite(candidate)
        a_values[fill] <- candidate[fill]
      }
      a_status <- if (!is.na(absorbing_ok_col)) {
        truthy(markov_diag[[absorbing_ok_col]])
      } else {
        rep(FALSE, nrow(markov_diag))
      }
      row_ok <- (is.finite(a_values) & abs(a_values) <= 1e-8) |
        (!is.finite(a_values) & a_status)
      absorb_ok <- !any(completed & !row_ok)
      finite_completed <- a_values[completed & is.finite(a_values)]
      absorb_detail <- paste0(
        "max_abs=",
        if (length(finite_completed)) max(abs(finite_completed)) else NA_real_,
        ", failed_completed_rows=", sum(completed & !row_ok),
        ", metrics=", paste(absorbing_metric_cols, collapse = "|")
      )
    } else if (!is.na(absorbing_ok_col)) {
      a_ok <- truthy(markov_diag[[absorbing_ok_col]])
      absorb_ok <- !any(completed & !a_ok)
      absorb_detail <- paste0("failed_completed_rows=", sum(completed & !a_ok))
    } else {
      absorb_ok <- FALSE
      absorb_detail <- "absorbing-state metric/status missing"
    }
  }
  add_check(
    "MARKOV", "BOTH", "Q_ROW_SUM_ZERO",
    "Every completed fitted Q matrix has rows summing to zero",
    q_rows_ok, observed = q_detail, expected = "maximum absolute row sum <= 1e-6",
    source_file = resolved$markov_diagnostics
  )
  add_check(
    "MARKOV", "BOTH", "Q_DEATH_ABSORBING",
    "Every completed fitted Q matrix has an all-zero terminal absorbing-state row",
    absorb_ok, observed = absorb_detail, expected = "maximum absolute terminal-row entry <= 1e-8",
    source_file = resolved$markov_diagnostics
  )
}

# -----------------------------------------------------------------------------
# Transition probabilities and state-occupation probabilities
# -----------------------------------------------------------------------------

validate_probability_matrix <- function(dt, source_file) {
  edges <- extract_edges(dt)
  prob_col <- find_col(dt, c("probability", "transition_probability", "prob", "p", "estimate"))
  lower_col <- find_col(dt, c("conf_low", "ci_low", "lower", "lcl", "lower_95", "lower_0_95"))
  upper_col <- find_col(dt, c("conf_high", "ci_high", "upper", "ucl", "upper_95", "upper_0_95"))
  schema_ok <- !is.null(edges) && !is.na(prob_col)
  add_check(
    "MARKOV", "BOTH", "PMATRIX_SCHEMA",
    "Transition-probability output identifies matrix cells and probabilities",
    schema_ok,
    observed = compact_detail(paste0("probability=", prob_col), paste0("lower=", lower_col), paste0("upper=", upper_col)),
    expected = "from_state, to_state, probability", source_file = source_file
  )
  if (!schema_ok) return(invisible(FALSE))
  edges[, `.prob_verify` := numeric_value(get(prob_col))]
  completed <- status_mask(edges)
  active <- completed & !is.na(edges$.from_state) & !is.na(edges$.to_state)
  if (!any(active)) {
    add_check(
      "MARKOV", "BOTH", "PMATRIX_VALUES",
      "Completed transition probabilities are bounded and finite",
      TRUE, observed = "completed_matrix_rows=0", expected = "valid for every completed row",
      detail = "No model was forced to complete; terminal status/reason is checked separately.", source_file = source_file
    )
    add_check(
      "MARKOV", "BOTH", "PMATRIX_ROW_SUM",
      "Completed transition-probability matrices have rows summing to one",
      TRUE, observed = "completed_matrices=0", expected = "valid for every completed matrix",
      detail = "No model was forced to complete; terminal status/reason is checked separately.", source_file = source_file
    )
    return(invisible(TRUE))
  }
  bad_value <- sum(active & (!is.finite(edges$.prob_verify) | edges$.prob_verify < -1e-8 | edges$.prob_verify > 1 + 1e-8))
  interval_bad <- 0L
  if (!is.na(lower_col) && !is.na(upper_col)) {
    lower <- numeric_value(edges[[lower_col]])
    upper <- numeric_value(edges[[upper_col]])
    interval_bad <- sum(active & (
      !is.finite(lower) | !is.finite(upper) | lower < -1e-8 | upper > 1 + 1e-8 |
        lower > edges$.prob_verify + 1e-8 | edges$.prob_verify > upper + 1e-8
    ))
  }
  add_check(
    "MARKOV", "BOTH", "PMATRIX_VALUES",
    "Completed transition probabilities and confidence limits are within [0,1]",
    bad_value == 0L && interval_bad == 0L,
    observed = paste0("invalid_probabilities=", bad_value, ", invalid_intervals=", interval_bad),
    expected = "0; 0", source_file = source_file
  )

  candidate_keys <- c(
    "cohort", "spec_id", "specification_id", "analysis_id", "model_id", "model",
    "time", "horizon", "years", "prediction_time", "exposure", "exposure_value",
    "insomnia_burden", "insomnia_burden_code", "profile", "imputation", "weighted"
  )
  key_cols <- names(edges)[normalise_name(names(edges)) %chin% normalise_name(candidate_keys)]
  if (!length(key_cols)) {
    edges[, `.matrix_group_verify` := "all"]
  } else {
    edges[, `.matrix_group_verify` := do.call(paste, c(.SD, sep = "::")), .SDcols = key_cols]
  }
  matrix_rows <- edges[active]
  duplicate_cells <- matrix_rows[, .N, by = .(.matrix_group_verify, .from_state, .to_state)][N > 1L, .N]
  row_checks <- matrix_rows[, .(
    n_to = uniqueN(.to_state),
    row_sum = sum(.prob_verify)
  ), by = .(.matrix_group_verify, .from_state)]
  group_states <- matrix_rows[, .(
    n_states = uniqueN(c(.from_state, .to_state)),
    n_from = uniqueN(.from_state),
    n_to = uniqueN(.to_state)
  ), by = .matrix_group_verify]
  row_checks <- merge(row_checks, group_states[, .(.matrix_group_verify, n_states)], by = ".matrix_group_verify", all.x = TRUE)
  bad_rows <- row_checks[!is.finite(row_sum) | abs(row_sum - 1) > 1e-6 | n_to != n_states, .N]
  nonsquare <- group_states[n_from != n_states | n_to != n_states, .N]
  add_check(
    "MARKOV", "BOTH", "PMATRIX_ROW_SUM",
    "Each completed transition-probability matrix is square with every row summing to one",
    duplicate_cells == 0L && bad_rows == 0L && nonsquare == 0L,
    observed = paste0("duplicate_cells=", duplicate_cells, ", bad_rows=", bad_rows, ", nonsquare_matrices=", nonsquare),
    expected = "0; 0; 0", source_file = source_file
  )
  terminal_by_matrix <- matrix_rows[, .(
    terminal_state = max(c(.from_state, .to_state))
  ), by = .matrix_group_verify]
  terminal_rows <- merge(
    matrix_rows, terminal_by_matrix,
    by = ".matrix_group_verify", all.x = TRUE, sort = FALSE
  )[
    .from_state == terminal_state
  ]
  terminal_absorb_bad <- terminal_rows[
    (.to_state == terminal_state & abs(.prob_verify - 1) > 1e-6) |
      (.to_state != terminal_state & abs(.prob_verify) > 1e-8),
    .N
  ]
  add_check(
    "MARKOV", "BOTH", "PMATRIX_DEATH_ABSORBING",
    "Each transition-probability matrix's terminal-state row remains absorbing",
    terminal_absorb_bad == 0L,
    observed = paste0("invalid_terminal_cells=", terminal_absorb_bad), expected = "0",
    source_file = source_file
  )
  invisible(TRUE)
}

if (!is.null(tables$markov_probability)) {
  validate_optional_terminal_status(tables$markov_probability, "MARKOV", "PMATRIX", resolved$markov_probability)
  validate_probability_matrix(tables$markov_probability, resolved$markov_probability)
}

validate_occupation <- function(dt, source_file) {
  prob_col <- find_col(dt, c("probability", "state_probability", "occupation_probability", "estimate", "prob"))
  state_col <- find_col(dt, c("state", "state_msm", "to_state", "state_code"))
  schema_ok <- !any(is.na(c(prob_col, state_col)))
  add_check(
    "MARKOV", "BOTH", "STATE_OCCUPATION_SCHEMA",
    "State-occupation output identifies state and probability", schema_ok,
    observed = compact_detail(paste0("state=", state_col), paste0("probability=", prob_col)),
    expected = "state and probability", source_file = source_file
  )
  if (!schema_ok) return(invisible(FALSE))
  completed <- status_mask(dt)
  p <- numeric_value(dt[[prob_col]])
  active <- completed
  bad <- sum(active & (!is.finite(p) | p < -1e-8 | p > 1 + 1e-8))
  add_check(
    "MARKOV", "BOTH", "STATE_OCCUPATION_VALUES",
    "Completed state-occupation probabilities lie within [0,1]",
    bad == 0L, observed = paste0("invalid_rows=", bad), expected = "0", source_file = source_file
  )
  if (!any(active)) return(invisible(TRUE))
  candidate_keys <- c(
    "cohort", "spec_id", "specification_id", "analysis_id", "model_id", "model",
    "time", "horizon", "years", "prediction_time", "from_state", "start_state",
    "exposure", "exposure_value", "insomnia_burden", "insomnia_burden_code", "profile"
  )
  key_cols <- names(dt)[normalise_name(names(dt)) %chin% normalise_name(candidate_keys)]
  work <- copy(dt[active])
  work[, `.p_verify` := p[active]]
  if (!length(key_cols)) work[, `.group_verify` := "all"] else {
    work[, `.group_verify` := do.call(paste, c(.SD, sep = "::")), .SDcols = key_cols]
  }
  sums <- work[, .(probability_sum = sum(.p_verify), n_states = uniqueN(get(state_col))), by = .group_verify]
  bad_sums <- sums[!is.finite(probability_sum) | abs(probability_sum - 1) > 1e-6 | n_states < 2L, .N]
  add_check(
    "MARKOV", "BOTH", "STATE_OCCUPATION_SUM",
    "Completed state-occupation probabilities sum to one within each prediction profile/time",
    bad_sums == 0L,
    observed = paste0("bad_profile_times=", bad_sums), expected = "0", source_file = source_file
  )
  invisible(TRUE)
}

if (!is.null(tables$markov_occupation)) {
  validate_optional_terminal_status(tables$markov_occupation, "MARKOV", "STATE_OCCUPATION", resolved$markov_occupation)
  validate_occupation(tables$markov_occupation, resolved$markov_occupation)
}

# -----------------------------------------------------------------------------
# HR/CI/p-value legality for Markov and Cox estimates
# -----------------------------------------------------------------------------

parse_ci_text <- function(x) {
  z <- clean_text(x)
  matches <- regexec("([0-9.eE+-]+)[[:space:]]*[(\\[][[:space:]]*([0-9.eE+-]+)[[:space:]]*[,;][[:space:]]*([0-9.eE+-]+)", z, perl = TRUE)
  pieces <- regmatches(z, matches)
  estimate <- lower <- upper <- rep(NA_real_, length(z))
  ok <- lengths(pieces) >= 4L
  if (any(ok)) {
    estimate[ok] <- as.numeric(vapply(pieces[ok], `[[`, character(1L), 2L))
    lower[ok] <- as.numeric(vapply(pieces[ok], `[[`, character(1L), 3L))
    upper[ok] <- as.numeric(vapply(pieces[ok], `[[`, character(1L), 4L))
  }
  list(estimate = estimate, lower = lower, upper = upper)
}

validate_effect_table <- function(dt, label, source_file) {
  status_ok <- validate_optional_terminal_status(dt, "EFFECTS", paste0(toupper(label), "_EFFECTS"), source_file)
  completed <- status_mask(dt)
  hr_col <- find_col(dt, c("hazard_ratio", "hr", "exp_estimate", "effect_ratio"))
  lower_col <- find_col(dt, c("conf_low", "ci_low", "lower", "lcl", "hr_lower", "lower_95"))
  upper_col <- find_col(dt, c("conf_high", "ci_high", "upper", "ucl", "hr_upper", "upper_95"))
  p_col <- find_col(dt, c("p_value", "pvalue", "p", "p_val", "pvalue_two_sided"))
  ci_col <- find_col(dt, c("hr_ci", "estimate_ci", "effect_ci", "hazard_ratio_95_ci"))
  type_col <- find_col(dt, c("estimate_type", "effect_type", "scale", "parameter_type"))
  estimate_col <- find_col(dt, c("estimate", "effect", "value"))

  effect_rows <- rep(TRUE, nrow(dt))
  if (!is.na(type_col)) {
    type_text <- clean_text(dt[[type_col]])
    identified <- grepl("hazard|ratio|(^|[^a-z])hr([^a-z]|$)|covariate", type_text, ignore.case = TRUE, perl = TRUE)
    if (any(identified)) effect_rows <- identified
  }
  if (is.na(hr_col) && !is.na(estimate_col) && !is.na(type_col) && any(effect_rows)) hr_col <- estimate_col

  parsed <- if (!is.na(ci_col)) parse_ci_text(dt[[ci_col]]) else NULL
  hr <- if (!is.na(hr_col)) numeric_value(dt[[hr_col]]) else if (!is.null(parsed)) parsed$estimate else rep(NA_real_, nrow(dt))
  lower <- if (!is.na(lower_col)) numeric_value(dt[[lower_col]]) else if (!is.null(parsed)) parsed$lower else rep(NA_real_, nrow(dt))
  upper <- if (!is.na(upper_col)) numeric_value(dt[[upper_col]]) else if (!is.null(parsed)) parsed$upper else rep(NA_real_, nrow(dt))
  p <- if (!is.na(p_col)) numeric_value(dt[[p_col]]) else rep(NA_real_, nrow(dt))

  active <- completed & effect_rows
  schema_ok <- !is.na(p_col) && (!is.na(hr_col) || !is.na(ci_col)) &&
    ((!is.na(lower_col) && !is.na(upper_col)) || !is.na(ci_col))
  if (!any(active)) {
    # Zero completed models is not itself a verification failure; all rows still
    # require terminal status/reason and the output file must exist.
    add_check(
      "EFFECTS", "BOTH", paste0(toupper(label), "_HR_CI_P"),
      paste(label, "completed HR/CI/p-value rows are legal"),
      all(status_ok), observed = "completed_effect_rows=0",
      expected = "valid values for every completed row",
      detail = "No model was forced to complete; registry and failure reasons remain mandatory.",
      source_file = source_file
    )
  } else {
    invalid <- active & (
      !is.finite(hr) | hr <= 0 | !is.finite(lower) | lower <= 0 |
        !is.finite(upper) | upper <= 0 | lower > hr | hr > upper |
        !is.finite(p) | p < 0 | p > 1
    )
    add_check(
      "EFFECTS", "BOTH", paste0(toupper(label), "_HR_CI_P"),
      paste(label, "completed HR/CI/p-value rows are legal"),
      schema_ok && !any(invalid),
      observed = paste0("completed_effect_rows=", sum(active), ", invalid_rows=", sum(invalid)),
      expected = "HR/CI > 0; lower <= HR <= upper; 0 <= p <= 1",
      detail = if (!schema_ok) "Required HR, lower CI, upper CI or p-value field is missing" else "",
      source_file = source_file
    )
  }

  edges <- extract_edges(dt)
  if (!is.null(edges)) {
    text <- row_text(edges)
    alt_q <- grepl(
      "alternative.{0,20}q|alternative_q|direct.{0,20}s0.{0,10}s3|s0.{0,10}s3.{0,20}(open|direct)|nine.{0,10}edge|9.{0,5}edge",
      text, ignore.case = TRUE, perl = TRUE
    )
    main_direct <- sum(active & !alt_q & edges$.edge == "0->3", na.rm = TRUE)
    forbidden <- sum(active & !alt_q & !is.na(edges$.edge) &
                       !edges$.edge %chin% expected_edges & edges$.from_state != edges$.to_state, na.rm = TRUE)
    add_check(
      "EFFECTS", "BOTH", paste0(toupper(label), "_EDGE_CONTRACT"),
      paste(label, "main-effect rows respect the frozen eight-edge graph"),
      main_direct == 0L && forbidden == 0L,
      observed = paste0("main_S0->S3_rows=", main_direct, ", other_forbidden_rows=", forbidden),
      expected = "0; 0", source_file = source_file
    )
  }
  invisible(TRUE)
}

if (!is.null(tables$markov_intensity)) {
  validate_effect_table(tables$markov_intensity, "MARKOV", resolved$markov_intensity)
}
if (!is.null(tables$cox_estimates)) {
  validate_effect_table(tables$cox_estimates, "COX", resolved$cox_estimates)
}

# Sojourn estimates are not HRs, but completed transient-state estimates must be positive.
if (!is.null(tables$markov_sojourn)) {
  validate_optional_terminal_status(tables$markov_sojourn, "MARKOV", "SOJOURN", resolved$markov_sojourn)
  completed <- status_mask(tables$markov_sojourn)
  estimate_col <- find_col(tables$markov_sojourn, c("sojourn_time", "estimate", "mean_sojourn_time", "years"))
  state_col <- find_col(tables$markov_sojourn, c("state", "state_msm", "state_code"))
  if (is.na(estimate_col)) {
    add_check("MARKOV", "BOTH", "SOJOURN_VALUES", "Completed transient-state sojourn estimates are positive",
              FALSE, "estimate field missing", ">0 for completed transient-state rows", source_file = resolved$markov_sojourn)
  } else {
    estimate <- numeric_value(tables$markov_sojourn[[estimate_col]])
    transient <- rep(TRUE, nrow(tables$markov_sojourn))
    if (!is.na(state_col)) {
      state_raw <- clean_text(tables$markov_sojourn[[state_col]])
      state_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", state_raw)))
      if (any(state_num == 5L, na.rm = TRUE) && !any(state_num == 0L, na.rm = TRUE)) state_num <- state_num - 1L
      transient <- is.na(state_num) | state_num != 4L
    }
    bad <- sum(completed & transient & (!is.finite(estimate) | estimate <= 0))
    add_check(
      "MARKOV", "BOTH", "SOJOURN_VALUES",
      "Completed transient-state sojourn estimates are positive",
      bad == 0L, observed = paste0("invalid_rows=", bad), expected = "0",
      source_file = resolved$markov_sojourn
    )
  }
}

# Cox diagnostic rows must themselves have terminal status/reason; PH p-values,
# where supplied, must be valid probabilities.
if (!is.null(tables$cox_diagnostics)) {
  validate_optional_terminal_status(tables$cox_diagnostics, "COX", "COX_DIAGNOSTICS", resolved$cox_diagnostics)
  completed <- status_mask(tables$cox_diagnostics)
  p_cols <- names(tables$cox_diagnostics)[grepl("(^p$|p_value|pvalue|ph_p|schoenfeld.*p)", normalise_name(names(tables$cox_diagnostics)), perl = TRUE)]
  invalid_p <- 0L
  assessed <- 0L
  for (p_col in p_cols) {
    p <- numeric_value(tables$cox_diagnostics[[p_col]])
    finite_rows <- completed & !is.na(p)
    assessed <- assessed + sum(finite_rows)
    invalid_p <- invalid_p + sum(finite_rows & (!is.finite(p) | p < 0 | p > 1))
  }
  add_check(
    "COX", "BOTH", "COX_DIAGNOSTIC_P_VALUES",
    "Reported Cox diagnostic p values are within [0,1]",
    length(p_cols) > 0L && invalid_p == 0L,
    observed = paste0("p_fields=", length(p_cols), ", assessed_values=", assessed, ", invalid=", invalid_p),
    expected = ">=1 p field; 0 invalid", source_file = resolved$cox_diagnostics
  )
}

# -----------------------------------------------------------------------------
# Event-per-parameter risk audit
# -----------------------------------------------------------------------------

event_audit_sources <- list(
  transition_counts = transition_counts,
  model_registry = registry,
  markov_diagnostics = markov_diag,
  cox_diagnostics = tables$cox_diagnostics,
  markov_intensity = markov_intensity,
  cox_estimates = tables$cox_estimates
)
event_audit_rows <- list()
for (source_name in names(event_audit_sources)) {
  dt <- event_audit_sources[[source_name]]
  if (is.null(dt) || !nrow(dt)) next
  edges <- extract_edges(dt)
  if (is.null(edges)) next
  cohort_col <- find_col(dt, c("cohort", "study"))
  event_col <- find_col(dt, c("events", "event_count", "n_events", "transition_events"))
  parameter_col <- find_col(dt, c("n_parameters", "parameter_count", "n_params", "parameters_per_transition"))
  epp_col <- find_col(dt, c("events_per_parameter", "event_per_parameter", "epp", "epv"))
  risk_col <- find_col(dt, c("event_parameter_risk", "epp_risk", "sparse_event_risk", "risk_flag", "event_risk"))
  if (is.na(cohort_col) || (is.na(parameter_col) && is.na(epp_col) && is.na(risk_col))) next
  event_audit_rows[[length(event_audit_rows) + 1L]] <- data.table(
    source = source_name,
    cohort = toupper(clean_text(dt[[cohort_col]])),
    edge = edges$.edge,
    events_reported = if (!is.na(event_col)) numeric_value(dt[[event_col]]) else NA_real_,
    n_parameters = if (!is.na(parameter_col)) numeric_value(dt[[parameter_col]]) else NA_real_,
    events_per_parameter = if (!is.na(epp_col)) numeric_value(dt[[epp_col]]) else NA_real_,
    risk = if (!is.na(risk_col)) clean_text(dt[[risk_col]]) else ""
  )
}

if (!length(event_audit_rows)) {
  add_check(
    "EVENTS", "BOTH", "EVENT_PARAMETER_AUDIT_PRESENT",
    "Every cohort-edge has an events-per-parameter or explicit sparse-event risk audit",
    FALSE, observed = "no parsable edge-level audit fields",
    expected = "cohort + edge + (n_parameters/events_per_parameter/risk)",
    source_file = resolved$model_registry
  )
} else {
  event_audit <- rbindlist(event_audit_rows, use.names = TRUE, fill = TRUE)
  if (!is.null(counts_standard)) {
    event_audit <- merge(
      event_audit,
      counts_standard[, .(cohort, edge, empirical_events = events)],
      by = c("cohort", "edge"), all.x = TRUE
    )
    event_audit[is.na(events_reported), events_reported := empirical_events]
  } else {
    event_audit[, empirical_events := NA_real_]
  }
  event_audit[, derived_epp := fifelse(
    is.finite(events_reported) & is.finite(n_parameters) & n_parameters > 0,
    events_reported / n_parameters, NA_real_
  )]
  event_audit[, audit_available := (
    is.finite(events_per_parameter) & events_per_parameter >= 0
  ) | (
    is.finite(derived_epp) & derived_epp >= 0
  ) | nzchar(risk)]
  required_grid <- CJ(cohort = c("ELSA", "HRS"), edge = expected_edges, sorted = TRUE)
  audit_grid <- unique(event_audit[audit_available & cohort %chin% c("HRS", "ELSA") & edge %chin% expected_edges, .(cohort, edge)])
  missing_audit <- fsetdiff(required_grid, audit_grid)
  invalid_parameters <- event_audit[!is.na(n_parameters) & (!is.finite(n_parameters) | n_parameters <= 0), .N]
  inconsistent_epp <- event_audit[
    is.finite(events_per_parameter) & is.finite(derived_epp) &
      abs(events_per_parameter - derived_epp) > 1e-6 * pmax(1, abs(derived_epp)),
    .N
  ]
  add_check(
    "EVENTS", "BOTH", "EVENT_PARAMETER_AUDIT_PRESENT",
    "Every cohort-edge has an events-per-parameter or explicit sparse-event risk audit",
    nrow(missing_audit) == 0L && invalid_parameters == 0L && inconsistent_epp == 0L,
    observed = paste0(
      "missing_edge_audits=", nrow(missing_audit), ", invalid_parameter_counts=", invalid_parameters,
      ", inconsistent_EPP=", inconsistent_epp
    ), expected = "0; 0; 0",
    detail = if (nrow(missing_audit)) paste(paste(missing_audit$cohort, missing_audit$edge, sep = ":"), collapse = ", ") else "",
    source_file = resolved$model_registry
  )
}

# -----------------------------------------------------------------------------
# Figure 1-7 delivery contract and figure QA
# -----------------------------------------------------------------------------

all_figure_files <- if (dir.exists(figure_root)) {
  list.files(figure_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
} else {
  character()
}

matches_figure_number <- function(path, number) {
  stem <- tools::file_path_sans_ext(basename(path))
  grepl(paste0("^fig(?:ure)?[ _-]*0?", number, "(?:[^0-9]|$)"), stem, ignore.case = TRUE, perl = TRUE)
}

for (figure_number in required_main_figures) {
  for (format in c("svg", "pdf", "tiff", "png")) {
    extensions <- if (format == "tiff") c("tiff", "tif") else format
    candidates <- all_figure_files[
      matches_figure_number(all_figure_files, figure_number) &
        tolower(tools::file_ext(all_figure_files)) %chin% extensions
    ]
    canonical_pattern <- file.path(figure_root, paste0("Figure", figure_number, "_*.", format))
    selected <- if (length(candidates)) candidates[1L] else canonical_pattern
    meta <- if (length(candidates)) lapply(candidates, file_metadata) else list(file_metadata(selected))
    all_nonempty <- length(candidates) > 0L && all(vapply(meta, `[[`, logical(1L), "nonempty"))
    manifest_rows[[length(manifest_rows) + 1L]] <- data.table(
      artifact_id = paste0("Figure", figure_number, "_", format),
      artifact_role = "figure",
      cohort = "BOTH",
      required = TRUE,
      expected_path = relative_path(canonical_pattern),
      actual_path = if (length(candidates)) paste(vapply(candidates, relative_path, character(1L)), collapse = " | ") else NA_character_,
      exists = length(candidates) > 0L,
      nonempty = all_nonempty,
      size_bytes = if (length(candidates)) sum(vapply(meta, `[[`, numeric(1L), "size"), na.rm = TRUE) else NA_real_,
      modified = if (length(candidates)) paste(vapply(meta, `[[`, character(1L), "modified"), collapse = " | ") else NA_character_,
      md5 = if (length(candidates)) paste(vapply(meta, `[[`, character(1L), "md5"), collapse = " | ") else NA_character_,
      verification_status = if (all_nonempty) "PASS" else "FAIL",
      note = paste0("matched_files=", length(candidates))
    )
    add_check(
      "FIGURES", "BOTH", paste0("FIGURE", figure_number, "_", toupper(format)),
      paste0("Figure ", figure_number, " has a non-empty ", toupper(format), " file"),
      all_nonempty,
      observed = paste0("matched_files=", length(candidates)), expected = ">=1; all size_bytes > 0",
      source_file = if (length(candidates)) candidates[1L] else selected
    )
  }

  source_primary <- file.path(figure_root, "source_data", paste0("Figure", figure_number, "_source_data.csv"))
  source_candidates <- c(
    source_primary,
    all_figure_files[
      matches_figure_number(all_figure_files, figure_number) &
        tolower(tools::file_ext(all_figure_files)) == "csv" &
        grepl("source|data", basename(all_figure_files), ignore.case = TRUE)
    ]
  )
  source_path <- register_artifact(
    paste0("Figure", figure_number, "_source_data"), "figure_source_data",
    unique(source_candidates), required = TRUE, cohort = "BOTH"
  )
  source_dt <- read_csv_artifact(source_path, paste0("Figure", figure_number, "_source_data"), require_rows = TRUE)
  if (!is.null(source_dt)) {
    all_empty <- nrow(source_dt) == 0L || ncol(source_dt) == 0L || all(vapply(source_dt, function(x) all(is.na(x) | !nzchar(clean_text(x))), logical(1L)))
    add_check(
      "FIGURES", "BOTH", paste0("FIGURE", figure_number, "_SOURCE_CONTENT"),
      paste0("Figure ", figure_number, " source CSV contains quantitative/nonblank content"),
      !all_empty,
      observed = paste0("rows=", nrow(source_dt), ", columns=", ncol(source_dt)),
      expected = "non-empty data", source_file = source_path
    )
    if (figure_number == 7L) {
      q_source_contract <- q_network_contract(source_dt)
      add_check(
        "FIGURES", "BOTH", "FIGURE7_Q_NETWORK_SOURCE_CONTRACT",
        "Figure 7 Source Data reproduces the exact 8/4/8 absolute-q network contract",
        q_source_contract$ok,
        observed = q_source_contract$detail,
        expected = paste(
          "20 unique positive q rows with valid 95% CIs; Rubin pooled 20/20;",
          "edge counts 8/4/8; no direct S0->S3"
        ),
        source_file = source_path
      )
    }
  }
}

figure_qa <- tables$figure_qa
if (!is.null(figure_qa)) {
  figure_col <- find_col(figure_qa, c("figure", "figure_id", "figure_number", "artifact", "file"))
  status_col <- find_col(figure_qa, c("status", "qa_status", "audit_status", "result", "check_status"))
  schema_ok <- !any(is.na(c(figure_col, status_col)))
  add_check(
    "QA", "BOTH", "FIGURE_QA_SCHEMA",
    "Figure QA audit identifies figure and audit status", schema_ok,
    observed = compact_detail(paste0("figure=", figure_col), paste0("status=", status_col)),
    expected = "figure and status", source_file = resolved$figure_qa
  )
  if (schema_ok) {
    fig_text <- clean_text(figure_qa[[figure_col]])
    main_figure_row <- !grepl("^\\s*(?:Extended|Ext)", fig_text, ignore.case = TRUE, perl = TRUE)
    fig_num <- suppressWarnings(as.integer(sub(
      "^\\s*[Ff]ig(?:ure)?[ _-]*0?([0-9]+)(?:[^0-9].*)?$", "\\1",
      fifelse(main_figure_row, fig_text, NA_character_), perl = TRUE
    )))
    numeric_direct <- suppressWarnings(as.integer(fig_text))
    direct_main <- is.na(fig_num) & main_figure_row & numeric_direct %in% required_main_figures
    fig_num[direct_main] <- numeric_direct[direct_main]
    qa_status <- normalise_name(figure_qa[[status_col]])
    accepted <- qa_status %chin% c("pass", "passed", "ok", "approved", "completed")
    figures_covered <- sort(unique(fig_num[fig_num %in% required_main_figures]))
    missing_figures <- setdiff(required_main_figures, figures_covered)
    failed_rows <- sum(!accepted)
    add_check(
      "QA", "BOTH", "FIGURE_QA_ALL_PASS",
      "Figure QA audit covers Figure 1-7 and every audit row passes",
      length(missing_figures) == 0L && failed_rows == 0L,
      observed = paste0("covered=", paste(figures_covered, collapse = ","), ", nonpass_rows=", failed_rows),
      expected = "figures 1-7; 0 non-PASS rows",
      detail = if (length(missing_figures)) paste("missing figures", paste(missing_figures, collapse = ", ")) else "",
      source_file = resolved$figure_qa
    )
  }
}

# Existing upstream QC tables are also part of the audit boundary.  Model
# registries legitimately contain failed/not-estimable statuses and are excluded;
# dedicated QC tables must not contain a FAIL result or a positive violation count.
if (dir.exists(qc_root)) {
  qc_files <- list.files(qc_root, pattern = "\\.csv(?:\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  qc_files <- qc_files[basename(qc_files) != "multistate_output_verification.csv"]
  qc_failures <- character()
  qc_examined <- 0L
  for (qc_file in qc_files) {
    qc_dt <- tryCatch(fread(qc_file, showProgress = FALSE), error = function(e) NULL)
    if (is.null(qc_dt) || !nrow(qc_dt)) next
    status_col <- find_col(qc_dt, c("status", "qa_status", "check_status", "result"))
    violation_col <- find_col(qc_dt, c("violations", "violation_count", "n_violations"))
    if (is.na(status_col) && is.na(violation_col)) next
    qc_examined <- qc_examined + 1L
    file_failed <- FALSE
    if (!is.na(status_col)) {
      s <- normalise_name(qc_dt[[status_col]])
      file_failed <- file_failed || any(s %chin% c("fail", "failed", "error", "invalid"))
    }
    if (!is.na(violation_col)) {
      v <- numeric_value(qc_dt[[violation_col]])
      file_failed <- file_failed || any(is.finite(v) & v > 0)
    }
    if (file_failed) qc_failures <- c(qc_failures, basename(qc_file))
  }
  add_check(
    "QA", "BOTH", "UPSTREAM_QC_STATUS",
    "All readable upstream QC audits report no failures or positive violations",
    qc_examined > 0L && length(qc_failures) == 0L,
    observed = paste0("audits_examined=", qc_examined, ", failing_files=", length(qc_failures)),
    expected = ">=1 audit; 0 failing files", detail = paste(qc_failures, collapse = ", "),
    source_file = qc_root
  )
}

# -----------------------------------------------------------------------------
# Write verification and complete manifest before raising any failure
# -----------------------------------------------------------------------------

verification <- rbindlist(verification_rows, use.names = TRUE, fill = TRUE)
setorder(verification, domain, cohort, check_id)
verification_path <- file.path(qc_root, "multistate_output_verification.csv")
fwrite(verification, verification_path, na = "")

# Add every unregistered actual file so the manifest is both a contract and a
# complete inventory.  Do not hash the manifest itself recursively.
registered_actual <- unique(unlist(lapply(manifest_rows, function(x) clean_text(x$actual_path)), use.names = FALSE))
registered_actual <- registered_actual[nzchar(registered_actual)]
scan_roots <- c(data_root, result_root, figure_root, qc_root)
actual_files <- unique(unlist(lapply(scan_roots[dir.exists(scan_roots)], function(x) {
  list.files(x, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
}), use.names = FALSE))
actual_rel <- vapply(actual_files, relative_path, character(1L))
unregistered <- actual_files[!actual_rel %chin% registered_actual]
for (path in unregistered) {
  meta <- file_metadata(path)
  manifest_rows[[length(manifest_rows) + 1L]] <- data.table(
    artifact_id = paste0("additional_", normalise_name(relative_path(path))),
    artifact_role = "additional_output",
    cohort = "BOTH",
    required = FALSE,
    expected_path = NA_character_,
    actual_path = relative_path(path),
    exists = meta$exists,
    nonempty = meta$nonempty,
    size_bytes = meta$size,
    modified = meta$modified,
    md5 = meta$md5,
    verification_status = if (meta$nonempty) "PASS" else "FAIL",
    note = "Discovered during recursive output inventory"
  )
}

verification_meta <- file_metadata(verification_path)
manifest_rows[[length(manifest_rows) + 1L]] <- data.table(
  artifact_id = "multistate_output_verification",
  artifact_role = "quality_audit",
  cohort = "BOTH",
  required = TRUE,
  expected_path = relative_path(verification_path),
  actual_path = relative_path(verification_path),
  exists = verification_meta$exists,
  nonempty = verification_meta$nonempty,
  size_bytes = verification_meta$size,
  modified = verification_meta$modified,
  md5 = verification_meta$md5,
  verification_status = if (verification_meta$nonempty) "PASS" else "FAIL",
  note = "Generated by R/09_verify_multistate_outputs.R"
)

manifest_path <- file.path(report_root, "multistate_output_manifest.csv")
manifest_rows[[length(manifest_rows) + 1L]] <- data.table(
  artifact_id = "multistate_output_manifest",
  artifact_role = "report",
  cohort = "BOTH",
  required = TRUE,
  expected_path = relative_path(manifest_path),
  actual_path = relative_path(manifest_path),
  exists = TRUE,
  nonempty = TRUE,
  size_bytes = NA_real_,
  modified = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  md5 = NA_character_,
  verification_status = "PASS",
  note = "Generated in the current run; checksum intentionally omitted to avoid self-reference"
)

manifest <- rbindlist(manifest_rows, use.names = TRUE, fill = TRUE)
manifest[, `.required_sort` := as.integer(!required)]
setorder(manifest, .required_sort, artifact_role, artifact_id)
manifest[, `.required_sort` := NULL]
fwrite(manifest, manifest_path, na = "")

print(verification)
cat("\nVerification summary:\n")
print(verification[, .N, by = .(status)][order(status)])
cat("Manifest written to:", manifest_path, "\n")

n_fail <- verification[status == "FAIL", .N]
if (n_fail > 0L) {
  stop(
    sprintf(
      "Multistate output verification failed (%d checks). Inspect %s and %s",
      n_fail, verification_path, manifest_path
    ),
    call. = FALSE
  )
}

cat("All multistate output verification checks passed.\n")
