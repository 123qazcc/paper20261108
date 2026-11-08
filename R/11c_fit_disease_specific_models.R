#!/usr/bin/env Rscript

# Post-hoc secondary disease-specific extension.
#
# This script fits continuous-time panel-observed Markov models to the D0-D7
# hypertension/heart-disease/stroke combination lattice prepared by 11a. HRS
# additionally contains absorbing death D8; ELSA is restricted to living W4-W10
# observations.  Multi-component changes between visits remain panel-observed
# endpoints and are integrated over compatible latent single-edge paths by the
# msm likelihood; no artificial ordering is inserted.

if (.Platform$OS.type == "windows" &&
    !grepl("UTF-?8|utf8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
  locale_result <- suppressWarnings(try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE))
  if (inherits(locale_result, "try-error") ||
      !grepl("UTF-?8|utf8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
    stop("A UTF-8 LC_CTYPE locale is required for reproducible paths and report text", call. = FALSE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(msm)
  library(GlobalAgingR)
})

options(
  stringsAsFactors = FALSE, width = 180,
  contrasts = c(unordered = "contr.treatment", ordered = "contr.treatment")
)

script_version <- "2026-09-10-v2"
cache_schema_version <- "disease-specific-msm-compact-v2"
run_timestamp_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

env_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = if (default) "1" else "0")))
  value %chin% c("1", "true", "yes", "y", "on")
}
env_integer <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (length(value) != 1L || is.na(value)) as.integer(default) else value
}

project_root <- Sys.getenv("GA_PROJECT_ROOT", unset = ".")
if (!dir.exists(project_root)) stop("GA_PROJECT_ROOT does not exist: ", project_root, call. = FALSE)
script_file <- file.path(project_root, "R", "11c_fit_disease_specific_models.R")
if (!file.exists(script_file)) stop("Model script not found: ", script_file, call. = FALSE)

nextsteps_root <- file.path(
  project_root, "outputs", "hrs_elsa_multistate_nextsteps_20260909"
)
frozen_root <- file.path(
  project_root, "outputs", "hrs_elsa_multistate_final_20260909"
)
data_root <- file.path(nextsteps_root, "data_disease_specific")
mi_root <- file.path(frozen_root, "data_multistate", "mi_baseline")
result_root <- file.path(nextsteps_root, "results", "disease_specific")
figure_root <- file.path(nextsteps_root, "figures", "disease_specific")
source_data_root <- file.path(figure_root, "source_data")
report_root <- file.path(nextsteps_root, "reports", "disease_specific")
qc_root <- file.path(nextsteps_root, "qc", "disease_specific", "models")
model_root <- file.path(nextsteps_root, "models", "disease_specific")
log_root <- file.path(nextsteps_root, "logs", "disease_specific")
contract_file <- file.path(project_root, "analysis", "secondary_extension_contract.md")
for (path in c(result_root, figure_root, source_data_root, report_root, qc_root,
               model_root, log_root)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

run_mode <- tolower(Sys.getenv("GA_DS_RUN_MODE", unset = "all"))
valid_run_modes <- c("audit", "pilot", "mi", "pool", "figures", "all")
if (!run_mode %chin% valid_run_modes) {
  stop("GA_DS_RUN_MODE must be one of: ", paste(valid_run_modes, collapse = ", "), call. = FALSE)
}
cohort_filter <- toupper(Sys.getenv("GA_DS_COHORT", unset = "ALL"))
if (!cohort_filter %chin% c("ALL", "HRS", "ELSA")) {
  stop("GA_DS_COHORT must be ALL, HRS or ELSA", call. = FALSE)
}
scoped_qc_path <- function(filename) {
  if (cohort_filter == "ALL") file.path(qc_root, filename) else {
    file.path(qc_root, paste0(tolower(cohort_filter), "__", filename))
  }
}
force_refit <- env_flag("GA_DS_FORCE_REFIT", FALSE)
max_iterations <- env_integer("GA_DS_MAX_ITER", 1200L)
mi_start <- max(1L, min(20L, env_integer("GA_DS_MI_START", 1L)))
mi_end <- max(mi_start, min(20L, env_integer("GA_DS_MI_END", 20L)))
expected_m <- 20L

sha256_file <- function(path) digest::digest(file = path, algo = "sha256")
script_sha256 <- sha256_file(script_file)
atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  suffix <- if (grepl("\\.gz$", path, ignore.case = TRUE)) ".csv.gz" else ".csv"
  temporary <- tempfile(".tmp_", tmpdir = dirname(path), fileext = suffix)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  fwrite(as.data.table(x), temporary, na = "")
  if (!file.copy(temporary, path, overwrite = TRUE)) stop("Failed to publish ", path)
  invisible(path)
}
atomic_save_rds <- function(x, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(".tmp_", tmpdir = dirname(path), fileext = ".rds")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(x, temporary, compress = compress)
  if (!file.copy(temporary, path, overwrite = TRUE)) stop("Failed to publish ", path)
  invisible(path)
}
write_text <- function(lines, path) {
  temporary <- tempfile(".tmp_", tmpdir = dirname(path), fileext = ".txt")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  writeLines(enc2utf8(lines), temporary, useBytes = TRUE)
  if (!file.copy(temporary, path, overwrite = TRUE)) stop("Failed to publish ", path)
  invisible(path)
}
delivery_status_file <- file.path(qc_root, "disease_specific_delivery_status.csv")
write_delivery_status <- function(status, message) {
  atomic_fwrite(data.table(
    status = status,
    message = message,
    script_version = script_version,
    cache_schema_version = cache_schema_version,
    script_sha256 = script_sha256,
    run_timestamp_utc = run_timestamp_utc,
    required_formal_components = 40L,
    analysis_role = "post-hoc secondary extension"
  ), delivery_status_file)
}
fingerprint_table <- function(x, columns) {
  columns <- intersect(columns, names(x))
  digest::digest(copy(as.data.table(x)[, ..columns]), algo = "xxhash64", serialize = TRUE)
}

m4_covariates <- c(
  "age_bl_5", "female_bl", "race_bl", "education_level_bl",
  "wealth_quintile_bl", "married_bl", "current_smoker_bl",
  "current_drinker_bl", "inactive_bl", "bmi_bl_5", "diabetes_bl",
  "lung_disease_bl", "frequent_pain_bl"
)
binary_covariates <- c(
  "female_bl", "married_bl", "current_smoker_bl", "current_drinker_bl",
  "inactive_bl", "diabetes_bl", "lung_disease_bl", "frequent_pain_bl"
)
exposure <- "insomnia_score_primary"
target_transitions <- c("D0->D2", "D0->D3", "D1->D4", "D1->D5")
target_labels <- c(
  "D0->D2" = "No disease to heart disease only",
  "D0->D3" = "No disease to stroke only",
  "D1->D4" = "Hypertension only to hypertension plus heart disease",
  "D1->D5" = "Hypertension only to hypertension plus stroke"
)
allowed_constraint_schemes <- c("etiologic_groups", "collapsed_other_living")
validate_constraint_scheme <- function(value, context) {
  if (length(value) != 1L || is.na(value) || !value %chin% allowed_constraint_schemes) {
    stop(context, " has an unknown or missing constraint scheme: ", paste(value, collapse = ";"), call. = FALSE)
  }
  value
}

configs <- list(
  HRS = list(
    cohort = "HRS", analysis_id = "hrs_9state_disease_specific",
    n_states = 9L, death_state = 9L,
    panel_file = file.path(data_root, "hrs_disease_specific_panel.rds"),
    mi_file = file.path(mi_root, "hrs_5state_full__M4_MI20__baseline_long.rds"),
    mi_manifest = file.path(mi_root, "hrs_5state_full__M4_MI20__cache_manifest.rds"),
    phenotype_note = "HRS broad heart disease cannot be reliably separated into MI/angina/heart failure; HRS stroke includes TIA."
  ),
  ELSA = list(
    cohort = "ELSA", analysis_id = "elsa_8state_w4_w10_disease_specific",
    n_states = 8L, death_state = NULL,
    panel_file = file.path(data_root, "elsa_disease_specific_panel.rds"),
    mi_file = file.path(mi_root, "elsa_4state_w4_w10__M4_MI20__baseline_long.rds"),
    mi_manifest = file.path(mi_root, "elsa_4state_w4_w10__M4_MI20__cache_manifest.rds"),
    phenotype_note = "ELSA W4-W10 is non-fatal; ELSA stroke does not explicitly include TIA."
  )
)
if (cohort_filter != "ALL") configs <- configs[cohort_filter]

required_files <- c(
  file.path(data_root, "disease_specific_allowed_q_edges.csv"),
  file.path(data_root, "disease_specific_state_dictionary.csv"),
  contract_file,
  unlist(lapply(configs, function(z) c(z$panel_file, z$mi_file, z$mi_manifest)), use.names = FALSE)
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Required input(s) missing: ", paste(missing_files, collapse = "; "), call. = FALSE)

state_dictionary <- fread(file.path(data_root, "disease_specific_state_dictionary.csv"))
all_edges <- fread(file.path(data_root, "disease_specific_allowed_q_edges.csv"))
required_state_columns <- c(
  "state", "state_code", "hypertension", "heart_disease", "stroke", "death", "state_label"
)
required_edge_columns <- c(
  "edge_id", "cohort", "analysis_id", "from_state", "to_state", "transition",
  "from_label", "to_label", "target_edge", "death_edge"
)
if (length(setdiff(required_state_columns, names(state_dictionary)))) {
  stop("State dictionary schema is incomplete", call. = FALSE)
}
if (length(setdiff(required_edge_columns, names(all_edges)))) {
  stop("Allowed-Q edge schema is incomplete", call. = FALSE)
}
setnames(all_edges, "analysis_id", "source_analysis_id")
if ("edge_type" %chin% names(all_edges)) {
  setnames(all_edges, "edge_type", "source_edge_type_legacy")
} else {
  all_edges[, source_edge_type_legacy := NA_character_]
}
all_edges[, edge_analysis_role := fifelse(
  target_edge == 1L,
  "planned target within post-hoc secondary extension; frozen before these new fits",
  fifelse(death_edge == 1L, "nuisance death process", "nuisance living process")
)]
all_edges[, `:=`(from_index = from_state + 1L, to_index = to_state + 1L)]
setorder(all_edges, cohort, from_index, to_index)

locf_down <- function(x) {
  if (!length(x)) return(x)
  last <- NA_real_
  for (i in seq_along(x)) {
    if (!is.na(x[i])) last <- x[i] else if (!is.na(last)) x[i] <- last
  }
  x
}

prepare_panel <- function(config) {
  x <- as.data.table(readRDS(config$panel_file))
  required <- c(
    "id", "cohort", "wave", "event_type", "event_time", "exposure_entry_time",
    "disease_specific_state", "hypertension", "heart_disease", "stroke", exposure
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(config$cohort, " panel lacks: ", paste(missing, collapse = ", "))
  x[, `:=`(
    state_msm = as.integer(disease_specific_state) + 1L,
    time_msm = as.numeric(event_time) - as.numeric(exposure_entry_time),
    insomnia_score_primary = as.numeric(insomnia_score_primary)
  )]
  setorder(x, id, time_msm, state_msm, na.last = TRUE)
  x[, insomnia_score_primary := locf_down(insomnia_score_primary), by = id]
  if (config$cohort == "ELSA") {
    x <- x[event_type == "visit" & wave %in% 4:10 & state_msm %in% 1:8]
  } else {
    x <- x[state_msm %in% 1:9]
  }
  x <- x[!is.na(id) & is.finite(time_msm) & time_msm >= 0]
  setorder(x, id, time_msm, state_msm)
  eligible <- x[, .N, by = id][N >= 2L, id]
  x <- x[id %in% eligible]
  setorder(x, id, time_msm, state_msm)
  x[]
}

make_pairs <- function(panel) {
  x <- panel[, .(id, from_state = state_msm - 1L, time_msm)]
  setorder(x, id, time_msm)
  x[, `:=`(
    to_state = shift(from_state, type = "lead"),
    next_time = shift(time_msm, type = "lead")
  ), by = id]
  x <- x[!is.na(to_state) & next_time > time_msm]
  x[, `:=`(
    interval_years = next_time - time_msm,
    transition = paste0("D", from_state, "->D", to_state),
    state_changed = as.integer(from_state != to_state)
  )]
  x[]
}

attach_constraint_groups <- function(edges, scheme = "etiologic_groups") {
  x <- copy(edges)
  target_map <- setNames(paste0("target_", seq_along(target_transitions)), target_transitions)
  x[, exposure_constraint_group := unname(target_map[transition])]
  x[is.na(exposure_constraint_group) & death_edge == 1L,
    exposure_constraint_group := "death_common"]
  if (scheme == "etiologic_groups") {
    dict <- state_dictionary[state <= 7L]
    for (i in which(is.na(x$exposure_constraint_group))) {
      from <- dict[state == x$from_state[i]]
      to <- dict[state == x$to_state[i]]
      added <- c(
        hypertension = to$hypertension - from$hypertension,
        heart = to$heart_disease - from$heart_disease,
        stroke = to$stroke - from$stroke
      )
      x$exposure_constraint_group[i] <- paste0("other_", names(added)[which(added == 1L)])
    }
  } else if (scheme == "collapsed_other_living") {
    x[is.na(exposure_constraint_group), exposure_constraint_group := "other_living_common"]
  } else stop("Unknown constraint scheme: ", scheme)
  group_order <- c(
    paste0("target_", 1:4), "other_hypertension", "other_heart",
    "other_stroke", "other_living_common", "death_common"
  )
  x[, exposure_constraint_id := match(exposure_constraint_group, group_order)]
  used <- unique(x$exposure_constraint_id)
  remap <- setNames(seq_along(sort(used)), sort(used))
  x[, exposure_constraint_id := unname(remap[as.character(exposure_constraint_id)])]
  x[]
}

edge_support <- function(panel, edges, config) {
  pairs <- make_pairs(panel)
  out <- copy(edges)
  out[, `:=`(
    risk_intervals = vapply(from_state, function(s) pairs[from_state == s, .N], integer(1)),
    risk_participants = vapply(from_state, function(s) pairs[from_state == s, uniqueN(id)], integer(1)),
    person_years = vapply(from_state, function(s) pairs[from_state == s, sum(interval_years)], numeric(1)),
    direct_events = mapply(function(a, b) pairs[from_state == a & to_state == b, .N], from_state, to_state)
  )]
  groups <- out[, .(
    constraint_group_events = sum(direct_events),
    constraint_group_edges = .N
  ), by = exposure_constraint_group]
  out <- merge(out, groups, by = "exposure_constraint_group", all.x = TRUE, sort = FALSE)
  out[, `:=`(
    sparse_baseline_edge = as.integer(direct_events < 15L),
    thin_baseline_edge_warning = as.integer(direct_events < 25L),
    sparse_exposure_group = as.integer(constraint_group_events < 50L),
    target_support_ok = as.integer(target_edge == 0L | direct_events >= 50L),
    baseline_rate_per_year = fifelse(person_years > 0, direct_events / person_years, NA_real_),
    analysis_id = config$analysis_id,
    analysis_role = fifelse(
      target_edge == 1L,
      "planned target edge frozen before new fits within a post-hoc secondary extension",
      "nuisance or structural edge within a post-hoc secondary extension"
    )
  )]
  setorder(out, from_index, to_index)
  attr(out, "pairs") <- pairs
  out[]
}

make_q_initial <- function(support, config, minimum = 1e-5) {
  n <- config$n_states
  q <- matrix(0, n, n, dimnames = list(paste0("D", 0:(n - 1L)), paste0("D", 0:(n - 1L))))
  for (i in seq_len(nrow(support))) {
    rate <- support$baseline_rate_per_year[i]
    if (!is.finite(rate) || rate <= 0) rate <- minimum
    q[support$from_index[i], support$to_index[i]] <- max(rate, minimum)
  }
  diag(q) <- -rowSums(q)
  if (!is.null(config$death_state)) q[config$death_state, ] <- 0
  q
}

matrix_long <- function(mat, config, type) {
  z <- as.data.table(as.table(mat))
  setnames(z, c("from_label", "to_label", "value"))
  z[, `:=`(
    cohort = config$cohort, analysis_id = config$analysis_id,
    matrix_type = type,
    from_state = as.integer(sub("D", "", from_label)),
    to_state = as.integer(sub("D", "", to_label))
  )]
  z[, transition := paste0("D", from_state, "->D", to_state)]
  z[]
}

input_fingerprints <- rbindlist(lapply(configs, function(config) {
  paths <- c(panel = config$panel_file, frozen_mi20 = config$mi_file,
             frozen_mi_manifest = config$mi_manifest)
  rbindlist(lapply(seq_along(paths), function(i) {
    path <- unname(paths[i])
    info <- file.info(path)
    data.table(
      cohort = config$cohort, input_role = names(paths)[i],
      path = enc2utf8(normalizePath(path, winslash = "/")),
      bytes = as.numeric(info$size), modified_utc = format(info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      sha256 = sha256_file(path)
    )
  }), fill = TRUE)
}), fill = TRUE)
shared_input_paths <- c(
  allowed_q_edges = file.path(data_root, "disease_specific_allowed_q_edges.csv"),
  state_dictionary = file.path(data_root, "disease_specific_state_dictionary.csv"),
  secondary_extension_contract = contract_file,
  model_script = script_file
)
shared_input_paths <- shared_input_paths[file.exists(shared_input_paths)]
input_fingerprints <- rbindlist(list(
  input_fingerprints,
  rbindlist(lapply(seq_along(shared_input_paths), function(i) {
    path <- unname(shared_input_paths[i])
    info <- file.info(path)
    data.table(
      cohort = "ALL", input_role = names(shared_input_paths)[i],
      path = enc2utf8(normalizePath(path, winslash = "/")),
      bytes = as.numeric(info$size), modified_utc = format(info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      sha256 = sha256_file(path)
    )
  }), fill = TRUE)
), fill = TRUE)

read_frozen_mi <- function(config) {
  x <- as.data.table(readRDS(config$mi_file))
  required <- c("id", ".imp", m4_covariates)
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(config$cohort, " frozen MI20 lacks: ", paste(missing, collapse = ", "))
  x[, .imp := as.integer(.imp)]
  x <- x[.imp %in% 1:expected_m]
  setorder(x, .imp, id)
  duplicate_keys <- x[, .N, by = .(.imp, id)][N != 1L, .N]
  if (duplicate_keys) stop(config$cohort, " frozen MI20 has duplicated/missing .imp-id keys")
  x[]
}

audit_frozen_mi <- function(mi, panel, config) {
  manifest <- readRDS(config$mi_manifest)
  predictor_names <- character()
  if (is.matrix(manifest$predictor_matrix)) {
    predictor_names <- unique(c(rownames(manifest$predictor_matrix), colnames(manifest$predictor_matrix)))
  }
  specific_auxiliaries <- c(
    "hypertension", "heart_disease", "stroke", "disease_specific_state",
    "state_heart_lower", "state_heart_upper"
  )
  imputations <- sort(unique(mi$.imp))
  ids_per_imp <- mi[, uniqueN(id), by = .imp]
  panel_ids <- sort(unique(panel$id))
  id_set_audit <- rbindlist(lapply(imputations, function(i) {
    imputation_ids <- sort(unique(mi[.imp == i, id]))
    data.table(
      imputation = i,
      missing_panel_ids = length(setdiff(panel_ids, imputation_ids)),
      extra_mi_ids = length(setdiff(imputation_ids, panel_ids)),
      exact_panel_id_set = as.integer(identical(panel_ids, imputation_ids))
    )
  }))
  data.table(
    cohort = config$cohort,
    analysis_id = config$analysis_id,
    imputation_source = "Frozen primary-delivery M4-MI20 baseline covariate draws",
    m_expected = expected_m,
    m_observed = length(imputations),
    imputations_complete_1_to_20 = as.integer(identical(imputations, 1:expected_m)),
    participants_per_imputation_min = min(ids_per_imp$V1),
    participants_per_imputation_max = max(ids_per_imp$V1),
    disease_panel_participants = uniqueN(panel$id),
    panel_participants_covered = uniqueN(intersect(panel$id, mi[.imp == 1L, id])),
    exact_panel_id_set_all_imputations = as.integer(all(id_set_audit$exact_panel_id_set == 1L)),
    missing_panel_ids_max_across_imputations = max(id_set_audit$missing_panel_ids),
    extra_mi_ids_max_across_imputations = max(id_set_audit$extra_mi_ids),
    missing_m4_values_across_completed_draws = mi[, sum(!complete.cases(.SD)), .SDcols = m4_covariates],
    predictor_manifest_has_composite_state_auxiliary = as.integer(
      any(c("state_msm", "last_state_msm", "death_observed") %in% predictor_names)
    ),
    predictor_manifest_has_disease_specific_auxiliary = as.integer(
      any(specific_auxiliaries %in% predictor_names)
    ),
    imputation_compatibility = paste(
      "Existing covariate draws used composite-state outcome auxiliaries and not the D0-D8 disease-specific components;",
      "reuse is the extension's frozen-input convention and has a congeniality limitation."
    ),
    new_imputation_performed = 0L
  )
}

canonicalize_covariates <- function(x, config) {
  d <- copy(x)
  if (!"age_bl_5" %in% names(d) && "age_bl" %in% names(d)) {
    d[, age_bl_5 := (as.numeric(age_bl) - 65) / 5]
  }
  if (!"bmi_bl_5" %in% names(d) && "bmi_bl" %in% names(d)) {
    d[, bmi_bl_5 := (as.numeric(bmi_bl) - 25) / 5]
  }
  for (name in intersect(binary_covariates, names(d))) {
    value <- d[[name]]
    integer_value <- if (is.factor(value) || is.character(value)) {
      suppressWarnings(as.integer(as.character(value)))
    } else suppressWarnings(as.integer(value))
    set(d, j = name, value = factor(integer_value, levels = 0:1))
  }
  if ("education_level_bl" %in% names(d)) {
    d[, education_level_bl := factor(as.integer(education_level_bl), levels = 1:3)]
  }
  if ("wealth_quintile_bl" %in% names(d)) {
    d[, wealth_quintile_bl := factor(as.integer(wealth_quintile_bl), levels = 1:5)]
  }
  if ("race_bl" %in% names(d)) {
    if (config$cohort == "HRS") {
      d[, race_bl := factor(as.character(race_bl), levels = c(
        "Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other"
      ))]
    } else {
      d[, race_bl := factor(as.character(race_bl), levels = c("White", "Non-White"))]
    }
  }
  d[]
}

build_complete_case_sample <- function(panel, config) {
  baseline <- panel[, .SD[1L], by = id, .SDcols = intersect(c(
    "age_bl", "female_bl", "race_bl", "education_level_bl", "wealth_quintile_bl",
    "married_bl", "current_smoker_bl", "current_drinker_bl", "inactive_bl",
    "bmi_bl", "diabetes_bl", "lung_disease_bl", "frequent_pain_bl"
  ), names(panel))]
  baseline <- canonicalize_covariates(baseline, config)
  keep <- c("id", m4_covariates)
  d <- merge(
    panel[, setdiff(names(panel), c(m4_covariates, "age_bl", "bmi_bl")), with = FALSE],
    baseline[, ..keep], by = "id", all.x = TRUE, sort = FALSE
  )
  required <- c("id", "state_msm", "time_msm", exposure, m4_covariates)
  d <- d[complete.cases(d[, ..required])]
  eligible <- d[, .N, by = id][N >= 2L, id]
  d <- d[id %in% eligible]
  setorder(d, id, time_msm, state_msm)
  d[]
}

build_imputed_sample <- function(panel, mi, imputation, config) {
  baseline <- copy(mi[.imp == imputation, c("id", m4_covariates), with = FALSE])
  baseline <- canonicalize_covariates(baseline, config)
  d <- merge(
    panel[, setdiff(names(panel), c(m4_covariates, "age_bl", "bmi_bl")), with = FALSE],
    baseline, by = "id", all = FALSE, sort = FALSE
  )
  required <- c("id", "state_msm", "time_msm", exposure, m4_covariates)
  d <- d[complete.cases(d[, ..required])]
  eligible <- d[, .N, by = id][N >= 2L, id]
  d <- d[id %in% eligible]
  setorder(d, id, time_msm, state_msm)
  d[]
}

make_constraints <- function(sample, edge_map) {
  formula <- reformulate(c(exposure, m4_covariates))
  mm <- model.matrix(formula, data = sample)
  coefficient_names <- setdiff(colnames(mm), "(Intercept)")
  if (!exposure %chin% coefficient_names) stop("Exposure coefficient is absent from model matrix")
  constraints <- vector("list", length(coefficient_names))
  names(constraints) <- coefficient_names
  constraints[[exposure]] <- edge_map$exposure_constraint_id
  for (name in setdiff(coefficient_names, exposure)) {
    constraints[[name]] <- rep(1L, nrow(edge_map))
  }
  list(formula = formula, coefficient_names = coefficient_names, constraints = constraints)
}

capture_fit <- function(fun, log_file, context) {
  warnings <- character()
  error_message <- NULL
  value <- NULL
  started <- Sys.time()
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
  write_text(c(
    paste("context:", context),
    paste("started:", format(started, "%Y-%m-%d %H:%M:%S %z")),
    paste("elapsed_seconds:", sprintf("%.3f", elapsed)),
    if (length(warnings)) c("warnings:", paste0("- ", unique(warnings))) else "warnings: none",
    if (is.null(error_message)) "error: none" else paste("error:", error_message),
    "captured_output:", output
  ), log_file)
  list(value = value, warnings = unique(warnings), error = error_message,
       elapsed_seconds = elapsed, log_file = log_file)
}

fit_metrics <- function(fit, sample, config, model_id, scheme, run) {
  if (is.null(fit)) {
    return(data.table(
      cohort = config$cohort, analysis_id = config$analysis_id, model_id = model_id,
      constraint_scheme = scheme, status = "FAILED", error = run$error,
      elapsed_seconds = run$elapsed_seconds
    ))
  }
  ll <- tryCatch(logLik(fit), error = function(e) NA)
  ll_value <- suppressWarnings(as.numeric(ll))[1L]
  ll_df <- suppressWarnings(as.numeric(attr(ll, "df")))[1L]
  hessian <- if (is.matrix(fit$opt$hessian)) fit$opt$hessian else NULL
  eig <- if (is.matrix(hessian) && all(is.finite(hessian))) {
    tryCatch(eigen((hessian + t(hessian)) / 2, symmetric = TRUE, only.values = TRUE)$values,
             error = function(e) numeric())
  } else numeric()
  min_eig <- if (length(eig)) min(eig) else NA_real_
  max_eig <- if (length(eig)) max(eig) else NA_real_
  convergence <- suppressWarnings(as.integer(fit$opt$convergence))[1L]
  if (!length(convergence)) convergence <- NA_integer_
  cov_finite <- is.matrix(fit$covmat) && all(is.finite(fit$covmat))
  accepted <- identical(convergence, 0L) && isTRUE(fit$foundse) && cov_finite &&
    is.finite(min_eig) && min_eig > 1e-8 &&
    is.finite(ll_value) && is.finite(ll_df) && ll_df > 0
  data.table(
    cohort = config$cohort, analysis_id = config$analysis_id, model_id = model_id,
    constraint_scheme = scheme,
    status = if (accepted) "ACCEPTED" else "REJECTED_DIAGNOSTIC",
    n_records = nrow(sample), n_participants = uniqueN(sample$id),
    n_parameters = ll_df,
    log_likelihood = ll_value,
    optimizer_convergence_code = convergence,
    optimizer_message = paste(fit$opt$message, collapse = " | "),
    standard_errors_found = as.integer(isTRUE(fit$foundse)),
    covariance_matrix_finite = as.integer(cov_finite),
    hessian_finite = as.integer(length(eig) > 0L),
    hessian_positive_definite = as.integer(is.finite(min_eig) && min_eig > 1e-8),
    hessian_min_eigenvalue = min_eig,
    hessian_max_eigenvalue = max_eig,
    hessian_condition_number = if (is.finite(min_eig) && min_eig > 0) max_eig / min_eig else Inf,
    warning_count = length(run$warnings), warnings = paste(run$warnings, collapse = " | "),
    error = if (is.null(run$error)) "" else run$error,
    elapsed_seconds = run$elapsed_seconds,
    independent_unit = "participant",
    observation_scheme = "continuous-time panel-observed Markov likelihood"
  )
}

extract_exposure_block <- function(fit, edges) {
  k <- nrow(edges)
  if (!identical(as.integer(fit$qmodel$npars), as.integer(k))) {
    stop("Fitted Q parameter count does not match allowed-edge count")
  }
  covariate_names <- fit$qcmodel$covlabels
  exposure_column <- match(exposure, covariate_names)
  if (is.na(exposure_column)) stop("Unable to locate exposure in msm qcmodel")
  indices <- k + (exposure_column - 1L) * k + seq_len(k)
  if (max(indices) > length(fit$estimates) || max(indices) > nrow(fit$covmat)) {
    stop("Expanded exposure parameter indices exceed fitted msm parameter dimensions")
  }
  theta <- as.numeric(fit$paramdata$params[indices])
  covariance <- as.matrix(fit$covmat[indices, indices, drop = FALSE])
  coefficient_matrices <- stats::coef(fit)
  theta_from_matrix <- vapply(seq_len(k), function(i) {
    coefficient_matrices[[exposure]][edges$from_index[i], edges$to_index[i]]
  }, numeric(1))
  if (!isTRUE(all.equal(theta, theta_from_matrix, tolerance = 1e-10,
                        check.attributes = FALSE))) {
    stop("Exposure coefficient-to-edge mapping failed")
  }
  if (any(!is.finite(theta)) || any(!is.finite(covariance))) {
    stop("Exposure block contains non-finite estimates or covariance")
  }
  covariance <- (covariance + t(covariance)) / 2
  if (any(diag(covariance) < 0)) stop("Exposure covariance has a negative diagonal element")
  covariance_eigen <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  if (min(covariance_eigen) < -1e-8) {
    stop("Exposure covariance is not positive semidefinite")
  }
  dimnames(covariance) <- list(edges$transition, edges$transition)
  names(theta) <- edges$transition
  list(theta = theta, covariance = covariance, indices = indices)
}

compact_fit <- function(fit, sample, config, model_id, scheme, edges, run, spec_hash,
                        model_data_fingerprint, fit_contract, initialization_record,
                        warm_start_source = "structural observed-rate initialization") {
  metrics <- fit_metrics(fit, sample, config, model_id, scheme, run)
  block <- if (!is.null(fit)) tryCatch(extract_exposure_block(fit, edges), error = function(e) e) else NULL
  if (inherits(block, "error")) {
    metrics[, `:=`(status = "REJECTED_EXTRACTION", error = conditionMessage(block))]
    block <- NULL
  }
  q_estimate <- if (!is.null(fit)) {
    tryCatch(msm::qmatrix.msm(fit, covariates = "mean", ci = "none"), error = function(e) NULL)
  } else NULL
  covariate_warm_starts <- if (!is.null(fit)) {
    tryCatch({
      k <- nrow(edges)
      labels <- fit$qcmodel$covlabels
      setNames(lapply(seq_along(labels), function(j) {
        as.numeric(fit$paramdata$params[k + (j - 1L) * k + seq_len(k)])
      }), labels)
    }, error = function(e) NULL)
  } else NULL
  q_valid <- is.matrix(q_estimate) && all(dim(q_estimate) == c(config$n_states, config$n_states)) &&
    all(is.finite(q_estimate)) && all(vapply(seq_len(nrow(edges)), function(i) {
      q_estimate[edges$from_index[i], edges$to_index[i]] > 0
    }, logical(1)))
  warm_valid <- is.list(covariate_warm_starts) &&
    identical(names(covariate_warm_starts), fit$qcmodel$covlabels) &&
    all(vapply(covariate_warm_starts, function(value) {
      length(value) == nrow(edges) && all(is.finite(value))
    }, logical(1)))
  if (!is.null(fit) && !q_valid) {
    metrics[, `:=`(
      status = "REJECTED_Q_EXTRACTION",
      error = paste(trimws(error), "Q-at-means extraction was missing, non-finite, dimensionally invalid, or non-positive on an allowed edge")
    )]
  }
  if (!is.null(fit) && !warm_valid) {
    metrics[, `:=`(
      status = "REJECTED_WARM_PARAMETER_EXTRACTION",
      error = paste(trimws(error), "Expanded covariate parameter extraction failed schema/finite checks")
    )]
  }
  metrics[, warm_start_source := warm_start_source]
  list(
    script_version = script_version, cache_schema_version = cache_schema_version,
    script_sha256 = script_sha256, spec_hash = spec_hash, fit_contract = fit_contract,
    config = config,
    model_id = model_id, constraint_scheme = scheme,
    data_fingerprint = model_data_fingerprint, metrics = metrics,
    exposure_block = block, q_at_covariate_means = q_estimate,
    covariate_warm_starts = covariate_warm_starts,
    warm_start_source = warm_start_source,
    initialization_record = initialization_record,
    constraints = if (!is.null(fit)) fit$qcmodel$constr else NULL,
    covariate_labels = if (!is.null(fit)) fit$qcmodel$covlabels else NULL,
    warnings = run$warnings, error = run$error, log_file = run$log_file
  )
}

make_spec_signature <- function(sample, config, edges, q_initial, scheme, model_id) {
  constraint_object <- make_constraints(sample, edges)
  model_data_fingerprint <- fingerprint_table(
    sample, c("id", "time_msm", "state_msm", exposure, m4_covariates)
  )
  fit_contract <- list(
    script_version = script_version,
    cache_schema_version = cache_schema_version,
    script_sha256 = script_sha256,
    r_version = as.character(getRversion()),
    msm_version = as.character(utils::packageVersion("msm")),
    globalagingr_version = as.character(utils::packageVersion("GlobalAgingR")),
    model = "continuous-time panel-observed msm; exact death only for HRS D8",
    formula = paste(deparse(constraint_object$formula), collapse = " "),
    center = TRUE, exacttimes = FALSE, hessian = TRUE,
    control = list(fnscale = 10000, maxit = max_iterations, trace = 0, reltol = 1e-9),
    contrasts = getOption("contrasts")
  )
  spec_hash <- digest::digest(list(
    fit_contract, config$analysis_id, model_id, scheme, edges[, .(
      from_state, to_state, exposure_constraint_id
    )], q_initial, constraint_object$coefficient_names, constraint_object$constraints,
    config$death_state, max_iterations, model_data_fingerprint,
    input_fingerprints[cohort %chin% c(config$cohort, "ALL"), sha256]
  ), algo = "sha256")
  list(
    spec_hash = spec_hash, constraint_object = constraint_object,
    model_data_fingerprint = model_data_fingerprint, fit_contract = fit_contract
  )
}

validate_compact_content <- function(object, edges, config, expected_fit_contract = NULL,
                                     expected_covariates = NULL) {
  required_fields <- c(
    "script_version", "cache_schema_version", "script_sha256", "spec_hash", "fit_contract",
    "metrics", "exposure_block", "q_at_covariate_means", "covariate_warm_starts",
    "initialization_record", "covariate_labels", "constraints"
  )
  fail <- function(reason) list(ok = FALSE, reason = reason)
  if (is.null(object) || !all(required_fields %in% names(object))) return(fail("CACHE_SCHEMA_FIELDS_MISSING"))
  metrics <- object$metrics
  required_metric_columns <- c(
    "status", "optimizer_convergence_code", "standard_errors_found",
    "covariance_matrix_finite", "hessian_positive_definite", "hessian_min_eigenvalue",
    "log_likelihood", "n_parameters"
  )
  if (!is.data.frame(metrics) || nrow(metrics) != 1L ||
      !all(required_metric_columns %in% names(metrics))) {
    return(fail("FIT_METRIC_SCHEMA_FIELDS_MISSING"))
  }
  diagnostic_pass <- isTRUE(metrics$status[1L] == "ACCEPTED") &&
    isTRUE(metrics$optimizer_convergence_code[1L] == 0L) &&
    isTRUE(metrics$standard_errors_found[1L] == 1L) &&
    isTRUE(metrics$covariance_matrix_finite[1L] == 1L) &&
    isTRUE(metrics$hessian_positive_definite[1L] == 1L) &&
    isTRUE(is.finite(metrics$hessian_min_eigenvalue[1L])) &&
    isTRUE(metrics$hessian_min_eigenvalue[1L] > 1e-8) &&
    isTRUE(is.finite(metrics$log_likelihood[1L])) &&
    isTRUE(is.finite(metrics$n_parameters[1L])) &&
    isTRUE(metrics$n_parameters[1L] > 0)
  if (!diagnostic_pass) {
    return(fail("FIT_DIAGNOSTIC_GATE_FAILED"))
  }
  if (!is.null(expected_fit_contract) && !identical(object$fit_contract, expected_fit_contract)) {
    return(fail("FIT_CONTRACT_MISMATCH"))
  }
  if (!is.null(expected_covariates) && !identical(object$covariate_labels, expected_covariates)) {
    return(fail("COVARIATE_LABEL_ORDER_MISMATCH"))
  }
  k <- nrow(edges)
  block <- object$exposure_block
  covariance <- if (is.list(block)) block$covariance else NULL
  expected_dimnames <- list(edges$transition, edges$transition)
  if (!is.list(block) || length(block$theta) != k || !identical(names(block$theta), edges$transition) ||
      any(!is.finite(block$theta)) || !is.matrix(covariance) || any(dim(covariance) != c(k, k)) ||
      !identical(dimnames(covariance), expected_dimnames) || any(!is.finite(covariance)) ||
      max(abs(covariance - t(covariance))) > 1e-8 || any(diag(covariance) < 0)) {
    return(fail("EXPOSURE_BLOCK_SCHEMA_OR_FINITE_GATE_FAILED"))
  }
  covariance_eigen <- tryCatch(
    eigen((covariance + t(covariance)) / 2, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  if (any(!is.finite(covariance_eigen)) || min(covariance_eigen) < -1e-8) {
    return(fail("EXPOSURE_COVARIANCE_NOT_POSITIVE_SEMIDEFINITE"))
  }
  q_value <- object$q_at_covariate_means
  if (!is.matrix(q_value) || any(dim(q_value) != c(config$n_states, config$n_states)) ||
      any(!is.finite(q_value)) || any(!vapply(seq_len(k), function(i) {
        q_value[edges$from_index[i], edges$to_index[i]] > 0
      }, logical(1)))) {
    return(fail("Q_AT_MEANS_SCHEMA_OR_FINITE_GATE_FAILED"))
  }
  warm <- object$covariate_warm_starts
  if (!is.list(warm) || !identical(names(warm), object$covariate_labels) ||
      any(!vapply(warm, function(value) length(value) == k && all(is.finite(value)), logical(1)))) {
    return(fail("WARM_PARAMETER_SCHEMA_OR_FINITE_GATE_FAILED"))
  }
  initialization <- object$initialization_record
  initialization_hash <- if (is.list(initialization)) digest::digest(
    initialization[c("source", "qmatrix", "covinits")], algo = "sha256"
  ) else NA_character_
  init_q <- if (is.list(initialization)) initialization$qmatrix else NULL
  init_cov <- if (is.list(initialization)) initialization$covinits else NULL
  if (!is.list(initialization) || !identical(initialization$sha256, initialization_hash) ||
      !is.matrix(init_q) || any(dim(init_q) != c(config$n_states, config$n_states)) ||
      any(!is.finite(init_q)) || !is.list(init_cov) ||
      !identical(names(init_cov), object$covariate_labels) ||
      any(!vapply(init_cov, function(value) length(value) == k && all(is.finite(value)), logical(1)))) {
    return(fail("INITIALIZATION_RECORD_GATE_FAILED"))
  }
  list(ok = TRUE, reason = "CONTENT_VALID")
}

fit_one <- function(sample, config, edges, q_initial, scheme, model_id, cache_file,
                    warm_start = NULL) {
  signature <- make_spec_signature(sample, config, edges, q_initial, scheme, model_id)
  constraint_object <- signature$constraint_object
  model_data_fingerprint <- signature$model_data_fingerprint
  spec_hash <- signature$spec_hash
  if (!force_refit && file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    content_check <- validate_compact_content(
      cached, edges, config, signature$fit_contract,
      signature$constraint_object$coefficient_names
    )
    reusable <- isTRUE(content_check$ok) &&
      identical(cached$script_version, script_version) &&
      identical(cached$cache_schema_version, cache_schema_version) &&
      identical(cached$script_sha256, script_sha256) &&
      identical(cached$spec_hash, spec_hash) && cached$metrics$status[1L] == "ACCEPTED"
    if (reusable) {
      cached$metrics[, cache_status := "REUSED_MATCHING_SPEC"]
      return(cached)
    }
  }
  if (file.exists(cache_file)) {
    archive_dir <- file.path(dirname(cache_file), "superseded_or_rejected")
    dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
    archive_file <- file.path(
      archive_dir,
      paste0(
        basename(cache_file), ".", format(Sys.time(), "%Y%m%dT%H%M%S"), ".",
        substr(sha256_file(cache_file), 1L, 10L), ".rds"
      )
    )
    if (!file.exists(archive_file) && !file.copy(cache_file, archive_file, overwrite = FALSE)) {
      stop("Unable to preserve non-reusable cache before refit: ", cache_file)
    }
  }
  log_file <- file.path(log_root, paste0(config$analysis_id, "__", model_id, "__", scheme, ".log"))
  q_start <- q_initial
  covariate_starts <- setNames(
    lapply(constraint_object$coefficient_names, function(name) numeric(nrow(edges))),
    constraint_object$coefficient_names
  )
  warm_start_source <- "structural observed-rate initialization; zero covariate effects"
  if (!is.null(warm_start)) {
    if (is.matrix(warm_start$qmatrix) && all(dim(warm_start$qmatrix) == dim(q_initial))) {
      q_start <- warm_start$qmatrix
    }
    if (is.list(warm_start$covinits) && length(warm_start$covinits)) {
      for (name in intersect(names(warm_start$covinits), constraint_object$coefficient_names)) {
        value <- as.numeric(warm_start$covinits[[name]])
        if (length(value) == nrow(edges) && all(is.finite(value))) covariate_starts[[name]] <- value
      }
    }
    warm_start_source <- warm_start$source
  }
  initialization_record <- list(
    source = warm_start_source,
    qmatrix = q_start,
    covinits = covariate_starts
  )
  initialization_record$sha256 <- digest::digest(
    initialization_record[c("source", "qmatrix", "covinits")], algo = "sha256"
  )
  run <- capture_fit(
    function() msm::msm(
      formula = state_msm ~ time_msm,
      subject = id,
      data = as.data.frame(sample),
      qmatrix = q_start,
      covariates = constraint_object$formula,
      covinits = covariate_starts,
      constraint = constraint_object$constraints,
      deathexact = config$death_state,
      exacttimes = FALSE,
      center = TRUE,
      hessian = TRUE,
      control = list(fnscale = 10000, maxit = max_iterations, trace = 0, reltol = 1e-9)
    ),
    log_file, paste(config$analysis_id, model_id, scheme, warm_start_source, sep = " | ")
  )
  result <- compact_fit(
    run$value, sample, config, model_id, scheme, edges, run, spec_hash,
    model_data_fingerprint, signature$fit_contract, initialization_record,
    warm_start_source
  )
  result$metrics[, cache_status := "CREATED"]
  atomic_save_rds(result, cache_file, compress = "xz")
  result
}

warm_start_from_compact <- function(object, edges, config, source) {
  if (is.null(object) || is.null(object$q_at_covariate_means) ||
      is.null(object$exposure_block) || !is.data.frame(object$metrics) || nrow(object$metrics) != 1L ||
      !isTRUE(object$metrics$status[1L] == "ACCEPTED")) {
    return(NULL)
  }
  raw_q <- as.matrix(object$q_at_covariate_means)
  if (any(dim(raw_q) != c(config$n_states, config$n_states))) return(NULL)
  q_start <- matrix(
    0, config$n_states, config$n_states,
    dimnames = list(paste0("D", 0:(config$n_states - 1L)), paste0("D", 0:(config$n_states - 1L)))
  )
  for (i in seq_len(nrow(edges))) {
    value <- suppressWarnings(as.numeric(raw_q[edges$from_index[i], edges$to_index[i]]))
    if (!is.finite(value) || value <= 0) return(NULL)
    q_start[edges$from_index[i], edges$to_index[i]] <- value
  }
  diag(q_start) <- -rowSums(q_start)
  if (!is.null(config$death_state)) q_start[config$death_state, ] <- 0
  covinits <- object$covariate_warm_starts
  if (is.null(covinits) || !length(covinits)) {
    covinits <- setNames(list(as.numeric(object$exposure_block$theta)), exposure)
  }
  covinits <- lapply(covinits, function(value) {
    value <- as.numeric(value)
    if (length(value) != nrow(edges) || any(!is.finite(value))) numeric(nrow(edges)) else value
  })
  list(qmatrix = q_start, covinits = covinits, source = source)
}

pilot_warm_start <- function(cohort, scheme, edges, config) {
  cache_file <- file.path(
    model_root, cohort,
    paste0(config$analysis_id, "__CC_PILOT__", scheme, ".rds")
  )
  if (!file.exists(cache_file)) return(NULL)
  object <- tryCatch(readRDS(cache_file), error = function(e) NULL)
  source_type <- if (!is.null(object) && identical(object$cache_schema_version, cache_schema_version)) {
    "current-schema accepted complete-case pilot"
  } else "exact-hash legacy pilot; Q plus insomnia exposure starts only; all other covariate starts explicitly zero"
  warm_start_from_compact(
    object, edges, config,
    paste0(source_type, ": ", normalizePath(cache_file, winslash = "/"))
  )
}

panels <- lapply(configs, prepare_panel)
mi_objects <- Map(read_frozen_mi, configs)
mi_compatibility <- rbindlist(Map(audit_frozen_mi, mi_objects, panels, configs), fill = TRUE)
formal_sample_mapping_audit <- rbindlist(lapply(names(configs), function(cohort_name) {
  config <- configs[[cohort_name]]
  panel <- panels[[cohort_name]]
  rbindlist(lapply(1:expected_m, function(imputation) {
    sample <- build_imputed_sample(panel, mi_objects[[cohort_name]], imputation, config)
    ordered_covariates <- m4_covariates[vapply(sample[, ..m4_covariates], is.ordered, logical(1))]
    row <- data.table(
      cohort = cohort_name, analysis_id = config$analysis_id, imputation = imputation,
      expected_panel_records = nrow(panel), model_sample_records = nrow(sample),
      expected_panel_participants = uniqueN(panel$id), model_sample_participants = uniqueN(sample$id),
      exact_id_set = as.integer(setequal(panel$id, sample$id)),
      missing_model_values = sample[, sum(!complete.cases(.SD)),
                                    .SDcols = c("state_msm", "time_msm", exposure, m4_covariates)],
      ordered_formal_covariates = if (length(ordered_covariates)) {
        paste(ordered_covariates, collapse = ";")
      } else "none",
      mapping_status = if (
        nrow(sample) == nrow(panel) && uniqueN(sample$id) == uniqueN(panel$id) &&
          setequal(panel$id, sample$id) && sample[, sum(!complete.cases(.SD)),
                                                 .SDcols = c("state_msm", "time_msm", exposure, m4_covariates)] == 0L &&
          !length(ordered_covariates)
      ) "PASS" else "FAIL"
    )
    rm(sample)
    row
  }), fill = TRUE)
}), fill = TRUE)

primary_edges <- lapply(names(configs), function(cohort_name) {
  attach_constraint_groups(all_edges[cohort == cohort_name], "etiologic_groups")
})
names(primary_edges) <- names(configs)
fallback_edges <- lapply(names(configs), function(cohort_name) {
  attach_constraint_groups(all_edges[cohort == cohort_name], "collapsed_other_living")
})
names(fallback_edges) <- names(configs)

sparse_audit <- rbindlist(lapply(names(configs), function(cohort_name) {
  edge_support(panels[[cohort_name]], primary_edges[[cohort_name]], configs[[cohort_name]])
}), fill = TRUE)

multijump_audit <- rbindlist(lapply(names(configs), function(cohort_name) {
  pairs <- make_pairs(panels[[cohort_name]])
  changed <- pairs[state_changed == 1L & to_state <= 7L]
  changed[, components_added := mapply(function(a, b) {
    from <- state_dictionary[state == a]
    to <- state_dictionary[state == b]
    sum(c(to$hypertension - from$hypertension,
          to$heart_disease - from$heart_disease,
          to$stroke - from$stroke))
  }, from_state, to_state)]
  out <- changed[components_added > 1L, .(
    intervals = .N, participants = uniqueN(id), person_years = sum(interval_years)
  ), by = .(transition, from_state, to_state, components_added)]
  out[, cohort := cohort_name]
  setcolorder(out, "cohort")
  out[]
}), fill = TRUE)

q_initials <- lapply(names(configs), function(cohort_name) {
  make_q_initial(sparse_audit[cohort == cohort_name], configs[[cohort_name]])
})
names(q_initials) <- names(configs)

q_initial_long <- rbindlist(lapply(names(configs), function(cohort_name) {
  matrix_long(q_initials[[cohort_name]], configs[[cohort_name]], "structural observed-rate initialization")
}), fill = TRUE)

state_dictionary_qc <- data.table(
  cohort = "ALL",
  check_id = c(
    "STATE_DICTIONARY_UNIQUE_COMPLETE_D0_D8",
    "LIVING_STATE_COMPONENT_BINARY",
    "LIVING_STATE_COMPONENT_CODE_UNIQUE",
    "D8_IS_DEATH_ONLY"
  ),
  observed = c(
    paste(sort(unique(state_dictionary$state)), collapse = ";"),
    state_dictionary[state <= 7L, sum(!hypertension %in% 0:1 | !heart_disease %in% 0:1 | !stroke %in% 0:1)],
    state_dictionary[state <= 7L, uniqueN(paste0(hypertension, heart_disease, stroke))],
    state_dictionary[state == 8L & death == 1L, .N]
  ),
  expected = c("0;1;2;3;4;5;6;7;8 with one row each", "0 invalid components", "8", "1"),
  status = c(
    ifelse(nrow(state_dictionary) == 9L && !anyDuplicated(state_dictionary$state) &&
             identical(sort(state_dictionary$state), 0:8), "PASS", "FAIL"),
    ifelse(state_dictionary[state <= 7L, sum(!hypertension %in% 0:1 | !heart_disease %in% 0:1 | !stroke %in% 0:1)] == 0L, "PASS", "FAIL"),
    ifelse(state_dictionary[state <= 7L, uniqueN(paste0(hypertension, heart_disease, stroke))] == 8L, "PASS", "FAIL"),
    ifelse(state_dictionary[state == 8L & death == 1L, .N] == 1L, "PASS", "FAIL")
  ),
  critical = 1L
)

model_structure_qc <- rbindlist(c(list(state_dictionary_qc), lapply(names(configs), function(cohort_name) {
  config <- configs[[cohort_name]]
  edges <- all_edges[cohort == cohort_name]
  panel <- panels[[cohort_name]]
  q <- q_initials[[cohort_name]]
  living_dict <- state_dictionary[state <= 7L][order(state)]
  living_edges <- edges[to_state <= 7L]
  increment_ok <- if (nrow(living_edges)) mapply(function(from, to) {
    from_bits <- as.integer(unlist(living_dict[state == from, .(hypertension, heart_disease, stroke)], use.names = FALSE))
    to_bits <- as.integer(unlist(living_dict[state == to, .(hypertension, heart_disease, stroke)], use.names = FALSE))
    delta <- to_bits - from_bits
    length(delta) == 3L && all(delta %in% 0:1) && sum(delta) == 1L
  }, living_edges$from_state, living_edges$to_state) else logical()
  pairs <- make_pairs(panel)
  observed_living <- pairs[to_state <= 7L]
  reachable <- if (nrow(observed_living)) mapply(function(from, to) {
    from_bits <- as.integer(unlist(living_dict[state == from, .(hypertension, heart_disease, stroke)], use.names = FALSE))
    to_bits <- as.integer(unlist(living_dict[state == to, .(hypertension, heart_disease, stroke)], use.names = FALSE))
    length(from_bits) == 3L && length(to_bits) == 3L && all(to_bits - from_bits >= 0L)
  }, observed_living$from_state, observed_living$to_state) else logical()
  allowed_mask <- matrix(FALSE, config$n_states, config$n_states)
  for (i in seq_len(nrow(edges))) allowed_mask[edges$from_index[i], edges$to_index[i]] <- TRUE
  off_diagonal <- q
  diag(off_diagonal) <- 0
  forbidden_nonzero <- sum(abs(off_diagonal[!allowed_mask]) > 1e-14)
  allowed_nonpositive <- sum(off_diagonal[allowed_mask] <= 0 | !is.finite(off_diagonal[allowed_mask]))
  diagonal_error <- max(abs(diag(q) + rowSums(off_diagonal)))
  expected_source_id <- if (cohort_name == "HRS") "hrs_9state_full" else "elsa_8state_w4_w10_nonfatal"
  expected_edge_count <- if (cohort_name == "HRS") 20L else 12L
  expected_death_edges <- if (cohort_name == "HRS") 8L else 0L
  expected_from_labels <- state_dictionary$state_label[match(edges$from_state, state_dictionary$state)]
  expected_to_labels <- state_dictionary$state_label[match(edges$to_state, state_dictionary$state)]
  state_label_mismatches <- sum(
    is.na(expected_from_labels) | is.na(expected_to_labels) |
      edges$from_label != expected_from_labels | edges$to_label != expected_to_labels
  )
  elsa_death_observations <- if (cohort_name == "ELSA") pairs[to_state == 8L, .N] else 0L
  data.table(
    cohort = cohort_name,
    check_id = c(
      "EDGE_KEYS_UNIQUE", "EDGE_COUNT_COMPLETE_HYPERCUBE", "TARGET_EDGE_SET_EXACT",
      "SINGLE_COMPONENT_LIVING_INCREMENT_ONLY", "DEATH_EDGE_SCOPE",
      "SOURCE_ANALYSIS_ID_DOCUMENTED", "TRANSITION_CODE_MATCHES_STATE_INDICES",
      "EDGE_STATE_LABELS_MATCH_DICTIONARY",
      "POSTHOC_TARGET_ROLE_OVERRIDES_LEGACY_SOURCE_DESCRIPTOR",
      "PANEL_COHORT_MATCH", "OBSERVED_TRANSITIONS_MONOTONE_REACHABLE",
      "Q_ALLOWED_EDGES_POSITIVE", "Q_FORBIDDEN_EDGES_ZERO", "Q_DIAGONAL_ROW_SUM",
      "Q_HRS_DEATH_ROW_ZERO"
    ),
    observed = c(
      edges[, sum(duplicated(edge_id) | duplicated(paste(from_state, to_state)))],
      nrow(edges),
      paste(sort(edges[target_edge == 1L, transition]), collapse = ";"),
      sum(!increment_ok),
      edges[death_edge == 1L, .N],
      paste(sort(unique(edges$source_analysis_id)), collapse = ";"),
      edges[transition != paste0("D", from_state, "->D", to_state), .N],
      state_label_mismatches,
      edges[target_edge == 1L & edge_analysis_role !=
              "planned target within post-hoc secondary extension; frozen before these new fits", .N],
      paste(sort(unique(panel$cohort)), collapse = ";"),
      sum(!reachable),
      allowed_nonpositive,
      forbidden_nonzero,
      format(diagonal_error, scientific = TRUE),
      if (cohort_name == "HRS") sum(abs(q[config$death_state, ]) > 1e-14) else 0L
    ),
    expected = c(
      "0", as.character(expected_edge_count), paste(sort(target_transitions), collapse = ";"),
      "0", as.character(expected_death_edges), expected_source_id, "0", "0", "0", cohort_name,
      "0", "0", "0", "<=1e-12", "0"
    ),
    status = c(
      ifelse(edges[, sum(duplicated(edge_id) | duplicated(paste(from_state, to_state)))] == 0L, "PASS", "FAIL"),
      ifelse(nrow(edges) == expected_edge_count, "PASS", "FAIL"),
      ifelse(identical(sort(edges[target_edge == 1L, transition]), sort(target_transitions)) &&
               edges[target_edge == 1L, .N] == 4L, "PASS", "FAIL"),
      ifelse(length(increment_ok) == 12L && all(increment_ok), "PASS", "FAIL"),
      ifelse(edges[death_edge == 1L, .N] == expected_death_edges &&
               (cohort_name != "HRS" || setequal(edges[death_edge == 1L, from_state], 0:7)) &&
               edges[death_edge == 1L & to_state != 8L, .N] == 0L, "PASS", "FAIL"),
      ifelse(identical(sort(unique(edges$source_analysis_id)), expected_source_id), "PASS", "FAIL"),
      ifelse(edges[transition != paste0("D", from_state, "->D", to_state), .N] == 0L, "PASS", "FAIL"),
      ifelse(state_label_mismatches == 0L, "PASS", "FAIL"),
      ifelse(edges[target_edge == 1L & edge_analysis_role !=
                     "planned target within post-hoc secondary extension; frozen before these new fits", .N] == 0L,
             "PASS", "FAIL"),
      ifelse(identical(sort(unique(panel$cohort)), cohort_name), "PASS", "FAIL"),
      ifelse(all(reachable) && elsa_death_observations == 0L, "PASS", "FAIL"),
      ifelse(allowed_nonpositive == 0L, "PASS", "FAIL"),
      ifelse(forbidden_nonzero == 0L, "PASS", "FAIL"),
      ifelse(is.finite(diagonal_error) && diagonal_error <= 1e-12, "PASS", "FAIL"),
      ifelse(cohort_name != "HRS" || sum(abs(q[config$death_state, ]) > 1e-14) == 0L, "PASS", "FAIL")
    ),
    critical = 1L
  )
})), fill = TRUE)

# GlobalAgingR is retained as the package-path empirical-Q audit.  The fitted
# model uses the structurally masked, direct-rate initialization above because
# the package empirical matrix is unconstrained and observed multi-jumps must
# not become direct Q edges.
ga_q_audit <- rbindlist(lapply(names(configs), function(cohort_name) {
  config <- configs[[cohort_name]]
  panel <- panels[[cohort_name]]
  value <- tryCatch(
    GlobalAgingR::estimate_Q_matrix(
      data = as.data.frame(panel[, .(id, time_msm, state_msm)]),
      method = "empirical", n_states = config$n_states,
      scale_factor = 0.1, smooth_factor = 0.01, validate = TRUE,
      verbose = FALSE, id_col = "id", time_col = "time_msm", state_col = "state_msm"
    ),
    error = function(e) e
  )
  if (inherits(value, "error") || !is.matrix(value) || any(dim(value) != config$n_states)) {
    return(data.table(
      cohort = cohort_name, analysis_id = config$analysis_id,
      status = "UNAVAILABLE_NOT_USED_FOR_FIT",
      message = if (inherits(value, "error")) conditionMessage(value) else "invalid matrix dimensions"
    ))
  }
  z <- data.table(
    from_index = rep(seq_len(nrow(value)), times = ncol(value)),
    to_index = rep(seq_len(ncol(value)), each = nrow(value)),
    value = as.numeric(value)
  )
  z[, `:=`(
    cohort = cohort_name, analysis_id = config$analysis_id,
    status = "AVAILABLE_AUDIT_ONLY",
    message = "Unconstrained empirical package matrix; not used to create forbidden direct edges",
    from_index = as.integer(from_index), to_index = as.integer(to_index)
  )]
  z[]
}), fill = TRUE)

input_qc <- rbindlist(lapply(names(configs), function(cohort_name) {
  config <- configs[[cohort_name]]
  panel <- panels[[cohort_name]]
  edges <- sparse_audit[cohort == cohort_name]
  mi_row <- mi_compatibility[cohort == cohort_name]
  mapping_rows <- formal_sample_mapping_audit[cohort == cohort_name]
  pairs <- make_pairs(panel)
  death_not_final <- if (cohort_name == "HRS") {
    panel[, any(state_msm == 9L & seq_len(.N) != .N), by = id][V1 == TRUE, .N]
  } else 0L
  data.table(
    cohort = cohort_name,
    analysis_id = config$analysis_id,
    check_id = c(
      "INPUT_STATE_RANGE", "STRICT_TIME_ORDER", "EXPOSURE_COMPLETE",
      "Q_EDGE_SUPPORT", "TARGET_EDGE_SUPPORT", "EXPOSURE_GROUP_SUPPORT",
      "MULTIJUMPS_PRESERVED", "DEATH_SCOPE", "MI20_COMPLETE",
      "MI20_M4_COMPLETE", "MI20_ID_SET_EXACT_EACH_IMPUTATION",
      "FORMAL_SAMPLE_MAPPING_PRESERVED", "MI_COMPONENT_AUXILIARY_LIMITATION_RECORDED",
      "THIN_EDGE_WARNING", "FORMAL_TREATMENT_FACTOR_PARAMETERIZATION"
    ),
    observed = c(
      panel[!state_msm %in% seq_len(config$n_states), .N],
      panel[, any(diff(time_msm) <= 0), by = id][V1 == TRUE, .N],
      panel[is.na(insomnia_score_primary), .N],
      edges[sparse_baseline_edge == 1L, .N],
      edges[target_support_ok != 1L, .N],
      edges[sparse_exposure_group == 1L, .N],
      multijump_audit[cohort == cohort_name, sum(intervals)],
      death_not_final,
      mi_row$m_observed,
      mi_row$missing_m4_values_across_completed_draws,
      paste0(mi_row$exact_panel_id_set_all_imputations, "; missing_max=", mi_row$missing_panel_ids_max_across_imputations,
             "; extra_max=", mi_row$extra_mi_ids_max_across_imputations),
      mapping_rows[mapping_status != "PASS", .N],
      mi_row$predictor_manifest_has_disease_specific_auxiliary,
      edges[thin_baseline_edge_warning == 1L, .N],
      paste(unique(mapping_rows$ordered_formal_covariates), collapse = ";")
    ),
    expected = c(
      "0", "0", "0", "0 edges with <15 direct events", "0 unsupported targets",
      "0 groups with <50 events", ">0 retained panel multijumps", "0",
      "20", "0", "1; missing_max=0; extra_max=0", "0 failed imputations",
      "0; explicit limitation required", "informational count of edges with <25 direct events",
      "none; all categorical adjustments use treatment contrasts"
    ),
    status = c(
      ifelse(panel[!state_msm %in% seq_len(config$n_states), .N] == 0L, "PASS", "FAIL"),
      ifelse(panel[, any(diff(time_msm) <= 0), by = id][V1 == TRUE, .N] == 0L, "PASS", "FAIL"),
      ifelse(panel[is.na(insomnia_score_primary), .N] == 0L, "PASS", "FAIL"),
      ifelse(edges[sparse_baseline_edge == 1L, .N] == 0L, "PASS", "FAIL"),
      ifelse(edges[target_support_ok != 1L, .N] == 0L, "PASS", "FAIL"),
      ifelse(edges[sparse_exposure_group == 1L, .N] == 0L, "PASS", "FAIL"),
      ifelse(multijump_audit[cohort == cohort_name, sum(intervals)] > 0L, "PASS", "FAIL"),
      ifelse(death_not_final == 0L, "PASS", "FAIL"),
      ifelse(mi_row$m_observed == expected_m, "PASS", "FAIL"),
      ifelse(mi_row$missing_m4_values_across_completed_draws == 0L, "PASS", "FAIL"),
      ifelse(mi_row$exact_panel_id_set_all_imputations == 1L &&
               mi_row$missing_panel_ids_max_across_imputations == 0L &&
               mi_row$extra_mi_ids_max_across_imputations == 0L, "PASS", "FAIL"),
      ifelse(nrow(mapping_rows) == expected_m && all(mapping_rows$mapping_status == "PASS"), "PASS", "FAIL"),
      ifelse(mi_row$predictor_manifest_has_disease_specific_auxiliary == 0L, "PASS_WITH_LIMITATION", "PASS"),
      ifelse(edges[thin_baseline_edge_warning == 1L, .N] > 0L, "PASS_WITH_WARNING", "PASS"),
      ifelse(all(mapping_rows$ordered_formal_covariates == "none"), "PASS", "FAIL")
    ),
    critical = c(rep(1L, 12L), 0L, 0L, 1L)
  )
}), fill = TRUE)

atomic_fwrite(input_fingerprints, scoped_qc_path("disease_specific_model_input_fingerprints.csv"))
atomic_fwrite(mi_compatibility, scoped_qc_path("disease_specific_mi20_compatibility.csv"))
atomic_fwrite(
  formal_sample_mapping_audit,
  scoped_qc_path("disease_specific_mi20_formal_sample_mapping_audit.csv")
)
atomic_fwrite(sparse_audit, scoped_qc_path("disease_specific_sparse_edge_audit.csv"))
atomic_fwrite(multijump_audit, scoped_qc_path("disease_specific_multijump_model_audit.csv"))
atomic_fwrite(q_initial_long, scoped_qc_path("disease_specific_q_initialization.csv"))
atomic_fwrite(ga_q_audit, scoped_qc_path("disease_specific_globalagingr_q_audit.csv"))
atomic_fwrite(input_qc, scoped_qc_path("disease_specific_model_input_qc.csv"))
atomic_fwrite(model_structure_qc, scoped_qc_path("disease_specific_model_structure_qc.csv"))

critical_input_failures <- rbindlist(list(
  input_qc[critical == 1L & status == "FAIL"],
  model_structure_qc[critical == 1L & status == "FAIL"]
), fill = TRUE)
if (nrow(critical_input_failures)) {
  print(critical_input_failures)
  stop("Disease-specific model input fail-fast gate failed", call. = FALSE)
}
if (run_mode == "audit") {
  cat("Disease-specific audit completed; no model requested.\n")
  quit(save = "no", status = 0L)
}

run_pilot_for_cohort <- function(cohort) {
  config <- configs[[cohort]]
  panel <- panels[[cohort]]
  sample <- build_complete_case_sample(panel, config)
  schemes <- c("etiologic_groups", "collapsed_other_living")
  attempts <- list()
  accepted_scheme <- NA_character_
  for (scheme in schemes) {
    edges <- if (scheme == "etiologic_groups") primary_edges[[cohort]] else fallback_edges[[cohort]]
    cache_file <- file.path(model_root, cohort, paste0(config$analysis_id, "__CC_PILOT__", scheme, ".rds"))
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    fit <- fit_one(
      sample, config, edges, q_initials[[cohort]], scheme,
      paste0("CC_PILOT__", scheme), cache_file
    )
    attempts[[length(attempts) + 1L]] <- fit$metrics
    if (fit$metrics$status[1L] == "ACCEPTED") {
      accepted_scheme <- scheme
      break
    }
  }
  metrics <- rbindlist(attempts, fill = TRUE)
  decision <- data.table(
    cohort = cohort, analysis_id = config$analysis_id,
    pilot_sample_records = nrow(sample), pilot_sample_participants = uniqueN(sample$id),
    selected_constraint_scheme = accepted_scheme,
    pilot_status = if (is.na(accepted_scheme)) "NO_ACCEPTED_PILOT" else "ACCEPTED",
    primary_scheme_attempted_first = 1L,
    fallback_rule = "Only used if the etiologic-group CC pilot fails convergence/finite covariance/positive-definite Hessian",
    post_hoc_secondary_extension = 1L,
    pilot_cache_validation = "ACCEPTED_MATCHING_CURRENT_SPEC",
    formal_component_cache_schema = cache_schema_version,
    formal_component_parameterization = "unordered factor treatment contrasts; pilot is never pooled"
  )
  atomic_fwrite(metrics, file.path(qc_root, paste0(tolower(cohort), "_disease_specific_cc_pilot_metrics.csv")))
  atomic_fwrite(decision, file.path(qc_root, paste0(tolower(cohort), "_disease_specific_pilot_decision.csv")))
  list(metrics = metrics, decision = decision, sample = sample)
}

legacy_equivalent_pilot_signature <- function(panel, config, edges, q_initial, scheme, model_id) {
  # The first CC screen was fitted before the formal-component cache schema was
  # frozen and used polynomial rather than treatment contrasts for education
  # and wealth.  Both are full-rank bases for the same categorical adjustment
  # space.  Reuse is permitted only for pilot scheme selection/warm starts and
  # only after reconstructing and matching its exact legacy specification hash.
  sample <- build_complete_case_sample(panel, config)
  sample[, education_level_bl := ordered(as.integer(education_level_bl), levels = 1:3)]
  sample[, wealth_quintile_bl := ordered(as.integer(wealth_quintile_bl), levels = 1:5)]
  previous_contrasts <- options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
  on.exit(options(previous_contrasts), add = TRUE)
  constraint_object <- make_constraints(sample, edges)
  model_data_fingerprint <- fingerprint_table(
    sample, c("id", "time_msm", "state_msm", exposure, m4_covariates)
  )
  digest::digest(list(
    "2026-09-09-v1", config$analysis_id, model_id, scheme, edges[, .(
      from_state, to_state, exposure_constraint_id
    )], q_initial, constraint_object$coefficient_names, constraint_object$constraints,
    config$death_state, max_iterations, model_data_fingerprint,
    input_fingerprints[cohort == config$cohort, sha256]
  ), algo = "sha256")
}

read_or_run_pilot <- function(cohort) {
  path <- file.path(qc_root, paste0(tolower(cohort), "_disease_specific_pilot_decision.csv"))
  if (file.exists(path) && !force_refit) {
    decision <- tryCatch(fread(path), error = function(e) NULL)
    required_decision_columns <- c(
      "cohort", "analysis_id", "selected_constraint_scheme", "pilot_status"
    )
    decision_schema_valid <- is.data.frame(decision) && nrow(decision) == 1L &&
      all(required_decision_columns %chin% names(decision)) &&
      isTRUE(decision$pilot_status[1L] == "ACCEPTED") &&
      isTRUE(decision$selected_constraint_scheme[1L] %chin% allowed_constraint_schemes)
    if (decision_schema_valid) {
      scheme <- validate_constraint_scheme(
        decision$selected_constraint_scheme[1L], paste0(cohort, " pilot decision")
      )
      config <- configs[[cohort]]
      edges <- if (scheme == "etiologic_groups") primary_edges[[cohort]] else fallback_edges[[cohort]]
      sample <- build_complete_case_sample(panels[[cohort]], config)
      model_id <- paste0("CC_PILOT__", scheme)
      cache_file <- file.path(model_root, cohort, paste0(config$analysis_id, "__", model_id, ".rds"))
      object <- if (file.exists(cache_file)) tryCatch(readRDS(cache_file), error = function(e) NULL) else NULL
      expected <- make_spec_signature(sample, config, edges, q_initials[[cohort]], scheme, model_id)
      current_content <- validate_compact_content(
        object, edges, config, expected$fit_contract,
        expected$constraint_object$coefficient_names
      )
      current_valid <- isTRUE(current_content$ok) && identical(object$spec_hash, expected$spec_hash) &&
        identical(object$script_version, script_version) &&
        identical(object$cache_schema_version, cache_schema_version) &&
        identical(object$script_sha256, script_sha256)
      legacy_hash <- legacy_equivalent_pilot_signature(
        panels[[cohort]], config, edges, q_initials[[cohort]], scheme, model_id
      )
      legacy_valid <- !is.null(object) && is.data.frame(object$metrics) && nrow(object$metrics) == 1L &&
        "status" %chin% names(object$metrics) && isTRUE(object$metrics$status[1L] == "ACCEPTED") &&
        identical(object$script_version, "2026-09-09-v1") && identical(object$spec_hash, legacy_hash)
      if (current_valid || legacy_valid) {
        decision[, `:=`(
          pilot_cache_validation = if (current_valid) {
            "ACCEPTED_MATCHING_CURRENT_SPEC"
          } else "ACCEPTED_MATCHING_EXACT_LEGACY_SPEC_EQUIVALENT_CATEGORICAL_BASIS",
          formal_component_cache_schema = cache_schema_version,
          formal_component_parameterization = "unordered factor treatment contrasts; legacy pilot is never pooled"
        )]
        atomic_fwrite(decision, path)
        metrics_path <- file.path(qc_root, paste0(tolower(cohort), "_disease_specific_cc_pilot_metrics.csv"))
        metrics <- if (file.exists(metrics_path)) {
          tryCatch(fread(metrics_path), error = function(e) copy(object$metrics))
        } else copy(object$metrics)
        if (!file.exists(metrics_path)) atomic_fwrite(metrics, metrics_path)
        return(list(decision = decision, metrics = metrics, sample = NULL))
      }
    }
  }
  run_pilot_for_cohort(cohort)
}

pilot_results <- lapply(names(configs), read_or_run_pilot)
names(pilot_results) <- names(configs)
if (run_mode == "pilot") {
  print(rbindlist(lapply(pilot_results, `[[`, "decision"), fill = TRUE))
  quit(save = "no", status = ifelse(any(vapply(pilot_results, function(z) z$decision$pilot_status != "ACCEPTED", logical(1))), 2L, 0L))
}

fit_mi_components <- function(cohort) {
  config <- configs[[cohort]]
  if (pilot_results[[cohort]]$decision$pilot_status[1L] != "ACCEPTED") {
    return(data.table(cohort = cohort, imputation = mi_start:mi_end, status = "NOT_RUN_NO_ACCEPTED_PILOT"))
  }
  scheme <- validate_constraint_scheme(
    pilot_results[[cohort]]$decision$selected_constraint_scheme[1L], paste0(cohort, " MI component fit")
  )
  edges <- if (scheme == "etiologic_groups") primary_edges[[cohort]] else fallback_edges[[cohort]]
  warm_start <- pilot_warm_start(cohort, scheme, edges, config)
  inventory_path <- file.path(qc_root, sprintf(
    "%s_disease_specific_mi_inventory_imp%02d_%02d__%s.csv",
    tolower(cohort), mi_start, mi_end, substr(script_sha256, 1L, 12L)
  ))
  inventory <- list()
  for (imputation in mi_start:mi_end) {
    message("Fitting ", cohort, " disease-specific MI ", imputation, "/", expected_m)
    sample <- build_imputed_sample(panels[[cohort]], mi_objects[[cohort]], imputation, config)
    model_id <- sprintf("M4_MI20__imp%02d", imputation)
    cache_file <- file.path(model_root, cohort, paste0(config$analysis_id, "__", model_id, "__", scheme, ".rds"))
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    fit <- fit_one(
      sample, config, edges, q_initials[[cohort]], scheme, model_id, cache_file,
      warm_start = warm_start
    )
    row <- copy(fit$metrics)
    row[, `:=`(
      imputation = imputation, cache_file = normalizePath(cache_file, winslash = "/"),
      spec_hash = fit$spec_hash, script_sha256 = fit$script_sha256,
      cache_schema_version = fit$cache_schema_version, recorded_utc = run_timestamp_utc
    )]
    inventory[[length(inventory) + 1L]] <- row
    current_inventory_snapshot <- rbindlist(inventory, fill = TRUE)
    setorder(current_inventory_snapshot, cohort, imputation)
    atomic_fwrite(current_inventory_snapshot, inventory_path)
    if (fit$metrics$status[1L] == "ACCEPTED") {
      warm_start <- warm_start_from_compact(
        fit, edges, config,
        paste0("previous accepted MI component cache: ", normalizePath(cache_file, winslash = "/"))
      )
    }
    rm(sample, fit)
    gc(verbose = FALSE)
  }
  current_inventory <- rbindlist(inventory, fill = TRUE)
  shard_status <- current_inventory[, .(
    mi_start = min(imputation), mi_end = max(imputation), components_requested = .N,
    components_accepted = sum(status == "ACCEPTED"),
    status = if (all(status == "ACCEPTED")) "SHARD_COMPLETE_ACCEPTED" else "SHARD_HAS_REJECTED_COMPONENTS",
    script_sha256 = unique(script_sha256), cache_schema_version = unique(cache_schema_version)
  ), by = .(cohort, analysis_id)]
  atomic_fwrite(
    shard_status,
    file.path(qc_root, sprintf(
      "%s_disease_specific_mi_shard_status_imp%02d_%02d__%s.csv",
      tolower(cohort), mi_start, mi_end, substr(script_sha256, 1L, 12L)
    ))
  )
  current_inventory
}

component_inventories <- list()
if (run_mode %chin% c("mi", "all")) {
  component_inventories <- lapply(names(configs), fit_mi_components)
  names(component_inventories) <- names(configs)
}
if (run_mode == "mi") {
  component_summary <- rbindlist(component_inventories, fill = TRUE)
  print(component_summary[, .N, by = .(cohort, status)])
  quit(save = "no", status = ifelse(any(component_summary$status != "ACCEPTED"), 3L, 0L))
}

rubin_scalar <- function(estimates, within_variances, m_expected = expected_m) {
  keep <- is.finite(estimates) & is.finite(within_variances) & within_variances >= 0
  estimates <- estimates[keep]
  within_variances <- within_variances[keep]
  m <- length(estimates)
  if (!m) return(list(m = 0L, estimate = NA_real_, variance = NA_real_, df = NA_real_))
  qbar <- mean(estimates)
  ubar <- mean(within_variances)
  b <- if (m > 1L) var(estimates) else 0
  total <- ubar + (1 + 1 / m) * b
  r <- if (ubar > 0) (1 + 1 / m) * b / ubar else if (b > 0) Inf else 0
  df <- if (m <= 1L) NA_real_ else if (r <= .Machine$double.eps) Inf else (m - 1) * (1 + 1 / r)^2
  se <- sqrt(total)
  statistic <- if (se > 0) qbar / se else NA_real_
  p <- if (!is.finite(statistic)) NA_real_ else if (is.finite(df)) {
    2 * pt(-abs(statistic), df = df)
  } else 2 * pnorm(-abs(statistic))
  critical <- if (is.finite(df)) qt(0.975, df) else qnorm(0.975)
  list(
    m_expected = m_expected, m = m, estimate = qbar,
    within_variance = ubar, between_variance = b, total_variance = total,
    standard_error = se, relative_increase_variance = r,
    rubin_missing_variance_fraction = if (total > 0) ((1 + 1 / m) * b) / total else 0,
    degrees_freedom = df, statistic = statistic, p_value = p,
    conf_low = qbar - critical * se, conf_high = qbar + critical * se
  )
}

strict_wald <- function(theta, covariance, label) {
  covariance <- (covariance + t(covariance)) / 2
  if (any(!is.finite(theta)) || any(!is.finite(covariance))) {
    return(list(
      status = "NOT_ESTIMABLE_NONFINITE", statistic = NA_real_, df = length(theta),
      p_value = NA_real_, rank = NA_integer_, rank_tolerance = NA_real_
    ))
  }
  eigenvalues <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- max(abs(eigenvalues)) * length(theta) * .Machine$double.eps^0.5
  rank <- sum(eigenvalues > tolerance)
  if (rank != length(theta)) {
    return(list(
      status = paste0("NOT_ESTIMABLE_RANK_", rank, "_OF_", length(theta)),
      statistic = NA_real_, df = length(theta), p_value = NA_real_,
      min_eigenvalue = min(eigenvalues), rank = rank, rank_tolerance = tolerance,
      label = label
    ))
  }
  solved <- tryCatch(solve(covariance, theta), error = function(e) NULL)
  if (is.null(solved)) {
    return(list(
      status = "NOT_ESTIMABLE_SOLVE_FAILED", statistic = NA_real_,
      df = length(theta), p_value = NA_real_, rank = rank,
      rank_tolerance = tolerance, min_eigenvalue = min(eigenvalues)
    ))
  }
  statistic <- as.numeric(crossprod(theta, solved))
  list(
    status = "ESTIMABLE_FULL_RANK", statistic = statistic, df = length(theta),
    p_value = pchisq(statistic, df = length(theta), lower.tail = FALSE),
    min_eigenvalue = min(eigenvalues), rank = rank,
    rank_tolerance = tolerance, label = label
  )
}

validated_component_cache <- function(cohort, imputation, scheme) {
  config <- configs[[cohort]]
  edges <- if (scheme == "etiologic_groups") primary_edges[[cohort]] else fallback_edges[[cohort]]
  model_id <- sprintf("M4_MI20__imp%02d", imputation)
  cache_file <- file.path(model_root, cohort, paste0(config$analysis_id, "__", model_id, "__", scheme, ".rds"))
  if (!file.exists(cache_file)) return(list(status = "MISSING_CACHE", cache_file = cache_file))
  object <- tryCatch(readRDS(cache_file), error = function(e) e)
  if (inherits(object, "error")) {
    return(list(status = "UNREADABLE_CACHE", cache_file = cache_file, message = conditionMessage(object)))
  }
  sample <- build_imputed_sample(panels[[cohort]], mi_objects[[cohort]], imputation, config)
  expected <- make_spec_signature(sample, config, edges, q_initials[[cohort]], scheme, model_id)
  content_check <- validate_compact_content(
    object, edges, config, expected$fit_contract,
    expected$constraint_object$coefficient_names
  )
  status <- if (!identical(object$cache_schema_version, cache_schema_version)) {
    "CACHE_SCHEMA_VERSION_MISMATCH"
  } else if (!identical(object$script_sha256, script_sha256)) {
    "SCRIPT_SHA256_MISMATCH"
  } else if (!identical(object$spec_hash, expected$spec_hash)) {
    "SPEC_HASH_MISMATCH"
  } else if (!identical(object$script_version, script_version)) {
    "SCRIPT_VERSION_MISMATCH"
  } else if (!isTRUE(content_check$ok)) {
    paste0("CONTENT_", content_check$reason)
  } else "ACCEPTED_MATCHING_SPEC"
  list(
    status = status, cache_file = cache_file, object = object,
    expected_spec_hash = expected$spec_hash, imputation = imputation
  )
}

pool_cohort <- function(cohort) {
  config <- configs[[cohort]]
  if (pilot_results[[cohort]]$decision$pilot_status[1L] != "ACCEPTED") {
    return(list(status = data.table(cohort = cohort, status = "NO_ACCEPTED_PILOT")))
  }
  scheme <- validate_constraint_scheme(
    pilot_results[[cohort]]$decision$selected_constraint_scheme[1L], paste0(cohort, " pooling")
  )
  edges <- if (scheme == "etiologic_groups") primary_edges[[cohort]] else fallback_edges[[cohort]]
  cache_checks <- lapply(1:expected_m, function(i) validated_component_cache(cohort, i, scheme))
  inventory <- rbindlist(lapply(cache_checks, function(z) data.table(
    cohort = cohort, imputation = z$imputation %||% NA_integer_, status = z$status,
    cache_file = normalizePath(z$cache_file, winslash = "/", mustWork = FALSE),
    expected_spec_hash = z$expected_spec_hash %||% NA_character_,
    observed_spec_hash = if (!is.null(z$object)) z$object$spec_hash else NA_character_
  )), fill = TRUE)
  atomic_fwrite(inventory, file.path(qc_root, paste0(tolower(cohort), "_disease_specific_pool_cache_validation.csv")))
  accepted <- Filter(function(z) identical(z$status, "ACCEPTED_MATCHING_SPEC"), cache_checks)
  if (length(accepted) != expected_m) {
    return(list(
      status = data.table(
        cohort = cohort, analysis_id = config$analysis_id,
        status = "PARTIAL_MI_NOT_REPORTABLE", m_expected = expected_m,
        m_accepted = length(accepted), message = "All 20 matching, accepted fits are required; no pooled inference written"
      ), inventory = inventory
    ))
  }
  objects <- lapply(accepted, `[[`, "object")
  component_metrics <- rbindlist(lapply(objects, `[[`, "metrics"), fill = TRUE)
  edge_labels <- edges$transition
  if (!all(vapply(objects, function(z) identical(names(z$exposure_block$theta), edge_labels), logical(1)))) {
    stop(cohort, " exposure edge ordering differs across imputations")
  }
  theta_matrix <- do.call(rbind, lapply(objects, function(z) z$exposure_block$theta))
  covariance_list <- lapply(objects, function(z) z$exposure_block$covariance)
  covariance_schema_ok <- vapply(covariance_list, function(z) {
    is.matrix(z) && all(dim(z) == c(length(edge_labels), length(edge_labels))) &&
      identical(dimnames(z), list(edge_labels, edge_labels)) && all(is.finite(z)) &&
      all(diag(z) >= 0)
  }, logical(1))
  if (any(!is.finite(theta_matrix)) || any(!covariance_schema_ok)) {
    stop(cohort, " accepted caches contain non-finite exposure estimates/covariances")
  }
  qbar <- colMeans(theta_matrix)
  ubar <- Reduce(`+`, covariance_list) / expected_m
  between <- cov(theta_matrix)
  total <- ubar + (1 + 1 / expected_m) * between
  dimnames(ubar) <- dimnames(between) <- dimnames(total) <- list(edge_labels, edge_labels)

  per_imputation <- rbindlist(lapply(seq_len(expected_m), function(i) {
    data.table(
      cohort = cohort, analysis_id = config$analysis_id, imputation = i,
      transition = edge_labels, from_state = edges$from_state, to_state = edges$to_state,
      exposure_constraint_group = edges$exposure_constraint_group,
      target_edge = edges$target_edge, death_edge = edges$death_edge,
      log_hazard_ratio = theta_matrix[i, ], hazard_ratio = exp(theta_matrix[i, ]),
      standard_error_log_hr = sqrt(diag(covariance_list[[i]]))
    )
  }), fill = TRUE)
  all_pooled <- rbindlist(lapply(seq_along(edge_labels), function(j) {
    scalar <- rubin_scalar(theta_matrix[, j], vapply(covariance_list, function(z) z[j, j], numeric(1)))
    data.table(
      cohort = cohort, analysis_id = config$analysis_id,
      transition = edge_labels[j], from_state = edges$from_state[j], to_state = edges$to_state[j],
      exposure_constraint_group = edges$exposure_constraint_group[j],
      target_edge = edges$target_edge[j], death_edge = edges$death_edge[j],
      m_expected = scalar$m_expected, m_fitted = scalar$m,
      pooling_status = if (scalar$m == expected_m) "RUBIN_POOLED_COMPLETE" else "PARTIAL_MI_NOT_REPORTABLE",
      log_hazard_ratio = scalar$estimate, standard_error_log_hr = scalar$standard_error,
      hazard_ratio = exp(scalar$estimate), conf_low = exp(scalar$conf_low), conf_high = exp(scalar$conf_high),
      degrees_freedom = scalar$degrees_freedom, statistic = scalar$statistic,
      p_value = scalar$p_value, within_imputation_variance = scalar$within_variance,
      between_imputation_variance = scalar$between_variance,
      total_variance = scalar$total_variance,
      rubin_missing_variance_fraction = scalar$rubin_missing_variance_fraction,
      scale = "per one additional insomnia symptom (0-3)",
      analysis_role = "post-hoc secondary extension; target edges frozen before these new fits"
    )
  }), fill = TRUE)
  if (any(all_pooled$m_fitted != expected_m) ||
      any(all_pooled$pooling_status != "RUBIN_POOLED_COMPLETE")) {
    stop(cohort, " scalar Rubin pooling dropped one or more accepted imputations")
  }

  target_positions <- match(target_transitions, edge_labels)
  if (anyNA(target_positions)) stop(cohort, " target edge missing from fitted edge map")
  target_theta <- qbar[target_positions]
  target_ubar <- ubar[target_positions, target_positions, drop = FALSE]
  target_between <- between[target_positions, target_positions, drop = FALSE]
  target_total <- total[target_positions, target_positions, drop = FALSE]
  per_imputation_target_covariance <- rbindlist(lapply(seq_len(expected_m), function(i) {
    block <- covariance_list[[i]][target_positions, target_positions, drop = FALSE]
    z <- as.data.table(as.table(block))
    setnames(z, c("row_transition", "column_transition", "covariance"))
    z[, `:=`(
      cohort = cohort, analysis_id = config$analysis_id, imputation = i,
      covariance_source = "within-imputation msm covariance"
    )]
    z[]
  }), fill = TRUE)
  target <- copy(all_pooled[match(target_transitions, transition)])
  target[, `:=`(
    target_label = unname(target_labels[transition]),
    target_order = seq_len(.N),
    target_outcome = rep(c("Heart disease", "Stroke"), 2L),
    origin_state_label = rep(c("No disease", "Hypertension only"), each = 2L),
    phenotype_note = config$phenotype_note
  )]

  contrast_matrix <- rbind(
    no_disease_heart_minus_stroke = c(1, -1, 0, 0),
    hypertension_heart_minus_stroke = c(0, 0, 1, -1)
  )
  contrast_rows <- rbindlist(lapply(seq_len(nrow(contrast_matrix)), function(i) {
    cvec <- contrast_matrix[i, ]
    per_imp_est <- as.numeric(theta_matrix[, target_positions, drop = FALSE] %*% cvec)
    per_imp_var <- vapply(covariance_list, function(z) {
      block <- z[target_positions, target_positions, drop = FALSE]
      as.numeric(crossprod(cvec, block %*% cvec))
    }, numeric(1))
    scalar <- rubin_scalar(per_imp_est, per_imp_var)
    data.table(
      cohort = cohort, analysis_id = config$analysis_id,
      contrast_id = rownames(contrast_matrix)[i],
      origin_state_label = c("No disease", "Hypertension only")[i],
      contrast = "heart-disease coefficient minus stroke coefficient",
      m_expected = scalar$m_expected, m_fitted = scalar$m,
      pooling_status = if (scalar$m == expected_m) "RUBIN_POOLED_COMPLETE" else "PARTIAL_MI_NOT_REPORTABLE",
      log_hr_ratio = scalar$estimate, standard_error = scalar$standard_error,
      hr_ratio_heart_vs_stroke = exp(scalar$estimate),
      conf_low = exp(scalar$conf_low), conf_high = exp(scalar$conf_high),
      degrees_freedom = scalar$degrees_freedom, statistic = scalar$statistic,
      p_value = scalar$p_value, within_imputation_variance = scalar$within_variance,
      between_imputation_variance = scalar$between_variance,
      total_variance = scalar$total_variance,
      interpretation = "HR ratio >1 indicates a larger insomnia coefficient for broad heart disease than stroke; a non-significant contrast is not equivalence",
      analysis_role = "planned contrast within a post-hoc secondary extension"
    )
  }), fill = TRUE)
  contrast_rows[, p_bh_two_contrasts_within_cohort := p.adjust(p_value, method = "BH")]
  if (any(contrast_rows$m_fitted != expected_m) ||
      any(contrast_rows$pooling_status != "RUBIN_POOLED_COMPLETE")) {
    stop(cohort, " contrast Rubin pooling dropped one or more accepted imputations")
  }

  target_joint <- strict_wald(target_theta, target_total, "four target coefficients jointly equal zero")
  contrast_theta <- as.numeric(contrast_matrix %*% target_theta)
  contrast_total <- contrast_matrix %*% target_total %*% t(contrast_matrix)
  contrast_joint <- strict_wald(contrast_theta, contrast_total, "two heart-versus-stroke contrasts jointly equal zero")
  joint <- rbindlist(list(
    data.table(
      cohort = cohort, analysis_id = config$analysis_id,
      test_id = "target_edges_4df", hypothesis = "four target insomnia coefficients jointly equal zero",
      statistic = target_joint$statistic, df = target_joint$df, p_value = target_joint$p_value,
      covariance_rank_status = target_joint$status,
      covariance_rank = target_joint$rank, covariance_rank_tolerance = target_joint$rank_tolerance,
      method = "Rubin total-covariance asymptotic Wald chi-square; not finite-m D1/F",
      m_expected = expected_m, m_fitted = expected_m
    ),
    data.table(
      cohort = cohort, analysis_id = config$analysis_id,
      test_id = "heart_vs_stroke_2df", hypothesis = "both within-origin heart-versus-stroke coefficient contrasts equal zero",
      statistic = contrast_joint$statistic, df = contrast_joint$df, p_value = contrast_joint$p_value,
      covariance_rank_status = contrast_joint$status,
      covariance_rank = contrast_joint$rank, covariance_rank_tolerance = contrast_joint$rank_tolerance,
      method = "Rubin total-covariance asymptotic Wald chi-square; not finite-m D1/F",
      m_expected = expected_m, m_fitted = expected_m
    )
  ), fill = TRUE)

  covariance_long <- rbindlist(lapply(c("within" = "within", "between" = "between", "total" = "total"), function(type) {
    mat <- switch(type, within = target_ubar, between = target_between, total = target_total)
    z <- as.data.table(as.table(mat))
    setnames(z, c("row_transition", "column_transition", "covariance"))
    z[, `:=`(cohort = cohort, analysis_id = config$analysis_id, covariance_component = type)]
    z[]
  }), fill = TRUE)

  q_values <- lapply(objects, `[[`, "q_at_covariate_means")
  q_long <- rbindlist(lapply(seq_along(q_values), function(i) {
    mat <- q_values[[i]]
    if (!is.matrix(mat)) return(NULL)
    z <- as.data.table(as.table(mat))
    setnames(z, c("from_label", "to_label", "intensity"))
    z[, `:=`(cohort = cohort, analysis_id = config$analysis_id, imputation = i)]
    z[]
  }), fill = TRUE)
  q_summary <- if (nrow(q_long)) q_long[, .(
    m_fitted = .N, intensity_mean = mean(intensity), intensity_between_sd = sd(intensity),
    intensity_across_imputation_p2_5_descriptive = quantile(intensity, 0.025),
    intensity_across_imputation_p97_5_descriptive = quantile(intensity, 0.975),
    quantile_note = "Descriptive distribution of 20 point estimates; not a confidence interval"
  ), by = .(cohort, analysis_id, from_label, to_label)] else data.table()

  covariance_gate_pass <- identical(target_joint$status, "ESTIMABLE_FULL_RANK") &&
    identical(contrast_joint$status, "ESTIMABLE_FULL_RANK")
  status <- data.table(
    cohort = cohort, analysis_id = config$analysis_id,
    status = if (covariance_gate_pass) "RUBIN_POOLED_COMPLETE" else "POOL_COVARIANCE_RANK_FAILURE",
    m_expected = expected_m, m_accepted = expected_m,
    selected_constraint_scheme = scheme,
    target_covariance_rank = target_joint$rank,
    target_covariance_rank_tolerance = target_joint$rank_tolerance,
    contrast_covariance_rank = contrast_joint$rank,
    contrast_covariance_rank_tolerance = contrast_joint$rank_tolerance,
    message = if (covariance_gate_pass) {
      "All 20 caches matched the current specification; all fit gates and full-rank target/contrast covariance gates passed"
    } else {
      "All 20 fit caches passed, but a target or contrast covariance block was rank deficient; no pseudoinverse result accepted"
    }
  )
  list(
    status = status, inventory = inventory, per_imputation = per_imputation,
    per_imputation_target_covariance = per_imputation_target_covariance,
    all_pooled = all_pooled, target = target, contrasts = contrast_rows,
    joint = joint, target_covariance = covariance_long, q_summary = q_summary,
    component_metrics = component_metrics,
    rubin_matrices = list(
      transition_order = target_transitions,
      qbar = target_theta,
      within = target_ubar,
      between = target_between,
      total = target_total,
      contrast_matrix = contrast_matrix,
      contrast_total = contrast_total
    )
  )
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

pooled <- list()
if (run_mode %chin% c("pool", "figures", "all")) {
  # A canonical all-cohort pooling pass is required for final inference.
  if (cohort_filter != "ALL") {
    stop("Pooling/figures require GA_DS_COHORT=ALL so multiplicity families are complete", call. = FALSE)
  }
  write_delivery_status(
    "IN_PROGRESS_NOT_FOR_INFERENCE",
    "Current-spec caches are being validated and pooled outputs have not passed the final inference/figure gates"
  )
  pooled <- lapply(names(configs), pool_cohort)
  names(pooled) <- names(configs)
  pool_status <- rbindlist(lapply(pooled, `[[`, "status"), fill = TRUE)
  atomic_fwrite(pool_status, file.path(qc_root, "disease_specific_mi20_pool_status.csv"))
  complete_pool <- nrow(pool_status) == length(configs) && all(pool_status$status == "RUBIN_POOLED_COMPLETE")
  if (complete_pool) {
    all_hr <- rbindlist(lapply(pooled, `[[`, "all_pooled"), fill = TRUE)
    target_hr <- rbindlist(lapply(pooled, `[[`, "target"), fill = TRUE)
    contrasts <- rbindlist(lapply(pooled, `[[`, "contrasts"), fill = TRUE)
    joint_tests <- rbindlist(lapply(pooled, `[[`, "joint"), fill = TRUE)
    target_covariance <- rbindlist(lapply(pooled, `[[`, "target_covariance"), fill = TRUE)
    per_imputation_target_covariance <- rbindlist(
      lapply(pooled, `[[`, "per_imputation_target_covariance"), fill = TRUE
    )
    per_imputation <- rbindlist(lapply(pooled, `[[`, "per_imputation"), fill = TRUE)
    q_summary <- rbindlist(lapply(pooled, `[[`, "q_summary"), fill = TRUE)
    component_metrics <- rbindlist(lapply(pooled, `[[`, "component_metrics"), fill = TRUE)
  }
}

markdown_table <- function(x) {
  d <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!ncol(d)) return(character())
  d[] <- lapply(d, function(value) {
    value <- ifelse(is.na(value), "NA", as.character(value))
    gsub("\\|", "\\\\|", value)
  })
  header <- paste0("| ", paste(names(d), collapse = " | "), " |")
  divider <- paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |")
  rows <- apply(d, 1L, function(value) paste0("| ", paste(value, collapse = " | "), " |"))
  c(header, divider, rows)
}

format_p <- function(x) {
  fifelse(is.na(x), "NA", fifelse(x < 0.001, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = 3)))
}

atomic_export_plot <- function(plot, stem, width, height) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("svglite", quietly = TRUE) ||
      !requireNamespace("ragg", quietly = TRUE)) {
    stop("ggplot2, svglite and ragg are required for the R-only figure exports")
  }
  exporters <- list(
    pdf = function(path) ggplot2::ggsave(
      path, plot = plot, device = grDevices::cairo_pdf,
      width = width, height = height, units = "in", bg = "white"
    ),
    svg = function(path) ggplot2::ggsave(
      path, plot = plot, device = svglite::svglite,
      width = width, height = height, units = "in", bg = "white"
    ),
    tiff = function(path) ggplot2::ggsave(
      path, plot = plot, device = ragg::agg_tiff,
      width = width, height = height, units = "in", res = 600,
      compression = "lzw", bg = "white"
    ),
    png = function(path) ggplot2::ggsave(
      path, plot = plot, device = ragg::agg_png,
      width = width, height = height, units = "in", res = 600, bg = "white"
    )
  )
  paths <- character()
  for (extension in names(exporters)) {
    destination <- file.path(figure_root, paste0(stem, ".", extension))
    temporary <- tempfile(".tmp_figure_", tmpdir = figure_root, fileext = paste0(".", extension))
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    exporters[[extension]](temporary)
    if (!file.exists(temporary) || file.info(temporary)$size <= 0L) {
      stop("Figure export failed: ", destination)
    }
    if (!file.copy(temporary, destination, overwrite = TRUE)) stop("Failed to publish ", destination)
    unlink(temporary, force = TRUE)
    paths <- c(paths, destination)
  }
  paths
}

if (run_mode %chin% c("pool", "figures", "all")) {
  if (!exists("complete_pool") || !isTRUE(complete_pool)) {
    write_delivery_status(
      "REJECTED_INCOMPLETE_MI20",
      "Fewer than 20/20 accepted current-spec component caches were available for at least one cohort; any older result artifacts are superseded and not valid for this run"
    )
    stop(
      "Disease-specific MI20 pooling was not complete. See disease_specific_mi20_pool_status.csv; no pooled secondary-extension figures/report were generated.",
      call. = FALSE
    )
  }

  model_specification <- rbindlist(lapply(names(configs), function(cohort_name) {
    config <- configs[[cohort_name]]
    scheme <- validate_constraint_scheme(
      pilot_results[[cohort_name]]$decision$selected_constraint_scheme[1L],
      paste0(cohort_name, " final model specification")
    )
    edges <- if (scheme == "etiologic_groups") primary_edges[[cohort_name]] else fallback_edges[[cohort_name]]
    support <- sparse_audit[cohort == cohort_name]
    out <- merge(
      edges[, .(
        cohort, edge_id, from_state, to_state, transition, from_label, to_label,
        target_edge, death_edge, exposure_constraint_id, exposure_constraint_group,
        source_edge_type_legacy, edge_analysis_role
      )],
      support[, .(
        cohort, transition, direct_events, risk_intervals, risk_participants,
        person_years, thin_baseline_edge_warning, baseline_rate_per_year
      )],
      by = c("cohort", "transition"), all.x = TRUE, sort = FALSE
    )
    out[, `:=`(
      analysis_id = config$analysis_id,
      selected_constraint_scheme = scheme,
      states = config$n_states,
      death_scope = if (cohort_name == "HRS") "absorbing D8 death" else "no death in W4-W10 secondary-extension model",
      insomnia_effect_role = fifelse(
        target_edge == 1L,
        "transition-specific planned target within post-hoc secondary extension",
        fifelse(death_edge == 1L, "shared across all death edges", paste0("shared within ", exposure_constraint_group))
      ),
      nuisance_covariate_constraint = "each M4 nuisance coefficient shared across all allowed Q edges",
      phenotype_note = config$phenotype_note,
      pilot_cache_validation = pilot_results[[cohort_name]]$decision$pilot_cache_validation[1L],
      pilot_role = "constraint-screen and initialization provenance only; pilot estimates are never pooled",
      formal_initialization_note = "initial values may inherit an accepted pilot or prior MI component; acceptance is evaluated independently under the frozen formal likelihood"
    )]
    setorder(out, from_state, to_state)
    out[]
  }), fill = TRUE)
  pool_inventory <- rbindlist(lapply(pooled, `[[`, "inventory"), fill = TRUE)
  atomic_fwrite(pool_inventory, file.path(qc_root, "disease_specific_pool_cache_validation_all.csv"))

  probability_ok <- function(x) {
    length(x) > 0L && all(is.finite(x) & x >= 0 & x <= 1)
  }
  positive_df_ok <- function(x) {
    length(x) > 0L && all(!is.na(x) & x > 0)
  }
  ordered_interval_ok <- function(point, lower, upper) {
    length(point) > 0L && all(
      is.finite(point) & is.finite(lower) & is.finite(upper) &
        lower > 0 & lower <= point & point <= upper
    )
  }
  pooled_scalar_numeric_ok <- function(x, point_column, se_column) {
    point <- x[[point_column]]
    standard_error <- x[[se_column]]
    all(
      is.finite(x$log_hazard_ratio %||% x$log_hr_ratio),
      is.finite(point), point > 0,
      is.finite(standard_error), standard_error > 0,
      ordered_interval_ok(point, x$conf_low, x$conf_high),
      positive_df_ok(x$degrees_freedom),
      is.finite(x$statistic), probability_ok(x$p_value),
      is.finite(x$within_imputation_variance), x$within_imputation_variance >= 0,
      is.finite(x$between_imputation_variance), x$between_imputation_variance >= 0,
      is.finite(x$total_variance), x$total_variance > 0
    )
  }
  audit_covariance_long <- function(x, group_columns, expected_blocks, expected_rows) {
    z <- copy(x)
    z[, `:=`(
      row_transition = as.character(row_transition),
      column_transition = as.character(column_transition)
    )]
    block_counts <- z[, .(
      rows = .N,
      unique_pairs = uniqueN(paste(row_transition, column_transition, sep = "|")),
      row_set_ok = setequal(unique(row_transition), target_transitions),
      column_set_ok = setequal(unique(column_transition), target_transitions),
      finite = all(is.finite(covariance))
    ), by = group_columns]
    reverse <- copy(z)
    temporary_transition <- reverse$row_transition
    reverse[, row_transition := column_transition]
    reverse[, column_transition := temporary_transition]
    setnames(reverse, "covariance", "reverse_covariance")
    symmetry <- merge(
      z, reverse,
      by = c(group_columns, "row_transition", "column_transition"),
      all = TRUE, sort = FALSE
    )
    symmetry_error <- if (nrow(symmetry) &&
                          all(is.finite(symmetry$covariance)) &&
                          all(is.finite(symmetry$reverse_covariance))) {
      max(abs(symmetry$covariance - symmetry$reverse_covariance))
    } else Inf
    psd_blocks <- z[, {
      mat <- matrix(
        NA_real_, length(target_transitions), length(target_transitions),
        dimnames = list(target_transitions, target_transitions)
      )
      row_index <- match(row_transition, target_transitions)
      column_index <- match(column_transition, target_transitions)
      complete_block <- .N == 16L && uniqueN(paste(row_transition, column_transition, sep = "|")) == 16L &&
        !anyNA(row_index) && !anyNA(column_index) && all(is.finite(covariance))
      if (complete_block) mat[cbind(row_index, column_index)] <- covariance
      eigenvalues <- if (complete_block && all(is.finite(mat))) {
        tryCatch(
          eigen((mat + t(mat)) / 2, symmetric = TRUE, only.values = TRUE)$values,
          error = function(e) NA_real_
        )
      } else NA_real_
      scale <- if (all(is.finite(eigenvalues))) max(1, max(abs(eigenvalues))) else NA_real_
      .(
        min_eigenvalue = if (all(is.finite(eigenvalues))) min(eigenvalues) else NA_real_,
        psd_tolerance = if (is.finite(scale)) scale * 1e-8 else NA_real_,
        psd_ok = all(is.finite(eigenvalues)) && min(eigenvalues) >= -scale * 1e-8
      )
    }, by = group_columns]
    psd_failures <- psd_blocks[psd_ok != TRUE, .N]
    minimum_eigenvalue <- if (nrow(psd_blocks) && any(is.finite(psd_blocks$min_eigenvalue))) {
      min(psd_blocks$min_eigenvalue[is.finite(psd_blocks$min_eigenvalue)])
    } else NA_real_
    diagonal_failures <- z[
      row_transition == column_transition & (!is.finite(covariance) | covariance < -1e-10), .N
    ]
    ok <- nrow(z) == expected_rows && nrow(block_counts) == expected_blocks &&
      all(block_counts$rows == 16L) && all(block_counts$unique_pairs == 16L) &&
      all(block_counts$row_set_ok) && all(block_counts$column_set_ok) &&
      all(block_counts$finite) && symmetry_error <= 1e-8 && diagonal_failures == 0L &&
      psd_failures == 0L
    list(
      ok = ok,
      observed = paste0(
        "rows=", nrow(z), "; blocks=", nrow(block_counts),
        "; symmetry_error=", format(symmetry_error, scientific = TRUE),
        "; negative_diagonal=", diagonal_failures,
        "; psd_failures=", psd_failures,
        "; min_eigenvalue=", format(minimum_eigenvalue, scientific = TRUE)
      )
    )
  }

  expected_cohorts <- c("HRS", "ELSA")
  expected_model_ids <- sprintf("M4_MI20__imp%02d", 1:expected_m)
  target_count_by_cohort <- target_hr[, .N, by = cohort]
  target_identity_ok <- setequal(target_count_by_cohort$cohort, expected_cohorts) &&
    all(target_count_by_cohort$N == 4L) &&
    all(target_hr[, .(
      ok = .N == 4L && uniqueN(transition) == 4L && setequal(transition, target_transitions)
    ), by = cohort]$ok)
  target_pool_ok <- nrow(target_hr) == 8L &&
    all(target_hr$m_expected == expected_m) && all(target_hr$m_fitted == expected_m) &&
    all(target_hr$pooling_status == "RUBIN_POOLED_COMPLETE") &&
    pooled_scalar_numeric_ok(target_hr, "hazard_ratio", "standard_error_log_hr") &&
    probability_ok(target_hr$rubin_missing_variance_fraction)

  contrast_ids <- c("no_disease_heart_minus_stroke", "hypertension_heart_minus_stroke")
  contrast_count_by_cohort <- contrasts[, .N, by = cohort]
  contrast_identity_ok <- nrow(contrasts) == 4L &&
    setequal(contrast_count_by_cohort$cohort, expected_cohorts) &&
    all(contrast_count_by_cohort$N == 2L) &&
    all(contrasts[, .(
      ok = .N == 2L && uniqueN(contrast_id) == 2L && setequal(contrast_id, contrast_ids)
    ), by = cohort]$ok)
  contrast_bh_recalculation <- contrasts[, .(
    ok = isTRUE(all.equal(
      p_bh_two_contrasts_within_cohort,
      p.adjust(p_value, method = "BH"),
      tolerance = 1e-14, check.attributes = FALSE
    ))
  ), by = cohort]
  contrast_pool_ok <- contrast_identity_ok &&
    all(contrasts$m_expected == expected_m) && all(contrasts$m_fitted == expected_m) &&
    all(contrasts$pooling_status == "RUBIN_POOLED_COMPLETE") &&
    pooled_scalar_numeric_ok(contrasts, "hr_ratio_heart_vs_stroke", "standard_error") &&
    probability_ok(contrasts$p_bh_two_contrasts_within_cohort) &&
    all(contrasts$p_bh_two_contrasts_within_cohort + 1e-15 >= contrasts$p_value) &&
    nrow(contrast_bh_recalculation) == 2L && all(contrast_bh_recalculation$ok)

  joint_ids <- c("target_edges_4df", "heart_vs_stroke_2df")
  expected_joint_df <- c(target_edges_4df = 4L, heart_vs_stroke_2df = 2L)
  joint_count_by_cohort <- joint_tests[, .N, by = cohort]
  joint_ok <- nrow(joint_tests) == 4L &&
    setequal(joint_count_by_cohort$cohort, expected_cohorts) &&
    all(joint_count_by_cohort$N == 2L) &&
    all(joint_tests[, .(
      ok = .N == 2L && uniqueN(test_id) == 2L && setequal(test_id, joint_ids)
    ), by = cohort]$ok) &&
    all(joint_tests$m_expected == expected_m) && all(joint_tests$m_fitted == expected_m) &&
    all(joint_tests$covariance_rank_status == "ESTIMABLE_FULL_RANK") &&
    all(joint_tests$df == unname(expected_joint_df[joint_tests$test_id])) &&
    all(joint_tests$covariance_rank == joint_tests$df) &&
    all(is.finite(joint_tests$covariance_rank_tolerance) & joint_tests$covariance_rank_tolerance >= 0) &&
    all(is.finite(joint_tests$statistic) & joint_tests$statistic >= -1e-10) &&
    probability_ok(joint_tests$p_value)

  all_edge_counts <- all_hr[, .N, by = cohort]
  all_edge_ok <- nrow(all_hr) == 32L && setequal(all_edge_counts$cohort, expected_cohorts) &&
    all_edge_counts[cohort == "HRS", N] == 20L && all_edge_counts[cohort == "ELSA", N] == 12L &&
    all(all_hr$m_expected == expected_m) && all(all_hr$m_fitted == expected_m) &&
    all(all_hr$pooling_status == "RUBIN_POOLED_COMPLETE") &&
    pooled_scalar_numeric_ok(all_hr, "hazard_ratio", "standard_error_log_hr") &&
    probability_ok(all_hr$rubin_missing_variance_fraction)

  per_imputation_counts <- per_imputation[, .N, by = .(cohort, imputation)]
  per_imputation_ok <- nrow(per_imputation) == 640L &&
    setequal(unique(per_imputation$cohort), expected_cohorts) &&
    all(per_imputation_counts[cohort == "HRS", N] == 20L) &&
    all(per_imputation_counts[cohort == "ELSA", N] == 12L) &&
    all(per_imputation_counts[, setequal(imputation, 1:expected_m), by = cohort]$V1) &&
    all(is.finite(per_imputation$log_hazard_ratio)) &&
    all(is.finite(per_imputation$hazard_ratio) & per_imputation$hazard_ratio > 0) &&
    all(is.finite(per_imputation$standard_error_log_hr) & per_imputation$standard_error_log_hr >= 0)

  per_imputation_covariance_audit <- audit_covariance_long(
    per_imputation_target_covariance,
    c("cohort", "analysis_id", "imputation"), expected_blocks = 40L, expected_rows = 640L
  )
  pooled_covariance_audit <- audit_covariance_long(
    target_covariance,
    c("cohort", "analysis_id", "covariance_component"), expected_blocks = 6L, expected_rows = 96L
  )
  pooled_covariance_components_ok <- setequal(
    unique(target_covariance$covariance_component), c("within", "between", "total")
  )

  rubin_scalar_matrix_consistency <- rbindlist(lapply(expected_cohorts, function(cohort_name) {
    matrices <- pooled[[cohort_name]]$rubin_matrices
    cohort_target <- target_hr[cohort == cohort_name][match(target_transitions, transition)]
    cohort_contrast <- contrasts[cohort == cohort_name][
      match(rownames(matrices$contrast_matrix), contrast_id)
    ]
    expected_contrast_estimate <- as.numeric(matrices$contrast_matrix %*% matrices$qbar)
    expected_contrast_within <- diag(
      matrices$contrast_matrix %*% matrices$within %*% t(matrices$contrast_matrix)
    )
    expected_contrast_between <- diag(
      matrices$contrast_matrix %*% matrices$between %*% t(matrices$contrast_matrix)
    )
    expected_contrast_total <- diag(matrices$contrast_total)
    errors <- c(
      target_estimate = max(abs(cohort_target$log_hazard_ratio - as.numeric(matrices$qbar))),
      target_within = max(abs(cohort_target$within_imputation_variance - diag(matrices$within))),
      target_between = max(abs(cohort_target$between_imputation_variance - diag(matrices$between))),
      target_total = max(abs(cohort_target$total_variance - diag(matrices$total))),
      contrast_estimate = max(abs(cohort_contrast$log_hr_ratio - expected_contrast_estimate)),
      contrast_within = max(abs(cohort_contrast$within_imputation_variance - expected_contrast_within)),
      contrast_between = max(abs(cohort_contrast$between_imputation_variance - expected_contrast_between)),
      contrast_total = max(abs(cohort_contrast$total_variance - expected_contrast_total))
    )
    data.table(
      cohort = cohort_name,
      maximum_absolute_error = if (all(is.finite(errors))) max(errors) else Inf,
      all_components_finite = as.integer(all(is.finite(errors)))
    )
  }), fill = TRUE)
  rubin_scalar_matrix_ok <- nrow(rubin_scalar_matrix_consistency) == 2L &&
    all(rubin_scalar_matrix_consistency$all_components_finite == 1L) &&
    all(rubin_scalar_matrix_consistency$maximum_absolute_error <= 1e-12)

  component_key_audit <- component_metrics[, .(
    ok = .N == expected_m && uniqueN(model_id) == expected_m &&
      setequal(model_id, expected_model_ids)
  ), by = cohort]
  component_diagnostics_ok <- nrow(component_metrics) == 40L &&
    setequal(component_key_audit$cohort, expected_cohorts) && all(component_key_audit$ok) &&
    all(component_metrics$status == "ACCEPTED") &&
    all(component_metrics$optimizer_convergence_code == 0L) &&
    all(component_metrics$standard_errors_found == 1L) &&
    all(component_metrics$covariance_matrix_finite == 1L) &&
    all(component_metrics$hessian_finite == 1L) &&
    all(component_metrics$hessian_positive_definite == 1L) &&
    all(is.finite(component_metrics$hessian_min_eigenvalue) & component_metrics$hessian_min_eigenvalue > 1e-8) &&
    all(is.finite(component_metrics$hessian_max_eigenvalue)) &&
    all(is.finite(component_metrics$hessian_condition_number) & component_metrics$hessian_condition_number >= 1) &&
    all(is.finite(component_metrics$log_likelihood)) &&
    all(is.finite(component_metrics$n_parameters) & component_metrics$n_parameters > 0) &&
    all(component_metrics$n_records > 0L & component_metrics$n_participants > 0L)

  cache_inventory_audit <- pool_inventory[, .(
    ok = .N == expected_m && uniqueN(imputation) == expected_m &&
      setequal(imputation, 1:expected_m) && all(status == "ACCEPTED_MATCHING_SPEC")
  ), by = cohort]
  cache_inventory_ok <- nrow(pool_inventory) == 40L &&
    setequal(cache_inventory_audit$cohort, expected_cohorts) && all(cache_inventory_audit$ok)

  mapping_key_audit <- formal_sample_mapping_audit[, .(
    ok = .N == expected_m && uniqueN(imputation) == expected_m &&
      setequal(imputation, 1:expected_m) && all(mapping_status == "PASS")
  ), by = cohort]
  mapping_ok <- nrow(formal_sample_mapping_audit) == 40L &&
    setequal(mapping_key_audit$cohort, expected_cohorts) && all(mapping_key_audit$ok)
  input_gate_ok <- setequal(unique(input_qc$cohort), expected_cohorts) &&
    input_qc[critical == 1L & status != "PASS", .N] == 0L
  structure_gate_ok <- all(c("ALL", expected_cohorts) %chin% unique(model_structure_qc$cohort)) &&
    model_structure_qc[critical == 1L & status != "PASS", .N] == 0L
  fingerprint_gate_ok <- nrow(input_fingerprints) == 10L &&
    all(is.finite(input_fingerprints$bytes) & input_fingerprints$bytes > 0) &&
    all(nchar(input_fingerprints$sha256) == 64L) &&
    setequal(
      input_fingerprints[cohort == "ALL", input_role],
      c("allowed_q_edges", "state_dictionary", "secondary_extension_contract", "model_script")
    )
  model_specification_counts <- model_specification[, .N, by = cohort]
  model_specification_ok <- nrow(model_specification) == 32L &&
    model_specification_counts[cohort == "HRS", N] == 20L &&
    model_specification_counts[cohort == "ELSA", N] == 12L &&
    all(model_specification$pilot_cache_validation %chin% c(
      "ACCEPTED_MATCHING_CURRENT_SPEC",
      "ACCEPTED_MATCHING_EXACT_LEGACY_SPEC_EQUIVALENT_CATEGORICAL_BASIS"
    ))
  q_summary_ok <- nrow(q_summary) == 145L && all(q_summary$m_fitted == expected_m) &&
    all(is.finite(q_summary$intensity_mean)) &&
    all(is.finite(q_summary$intensity_between_sd) & q_summary$intensity_between_sd >= 0) &&
    all(is.finite(q_summary$intensity_across_imputation_p2_5_descriptive)) &&
    all(is.finite(q_summary$intensity_across_imputation_p97_5_descriptive)) &&
    all(q_summary$intensity_across_imputation_p2_5_descriptive <=
          q_summary$intensity_across_imputation_p97_5_descriptive)

  final_qc <- rbindlist(list(
    pool_status[, .(
      cohort, check_id = "MI20_ACCEPTED_MATCHING_SPEC",
      observed = paste0(m_accepted, "/", m_expected), expected = "20/20",
      status = fifelse(status == "RUBIN_POOLED_COMPLETE" & m_accepted == 20L, "PASS", "FAIL"), critical = 1L
    )],
    pool_status[, .(
      cohort, check_id = "TARGET_RUBIN_COVARIANCE_FULL_RANK",
      observed = as.character(target_covariance_rank), expected = "4",
      status = fifelse(target_covariance_rank == 4L, "PASS", "FAIL"), critical = 1L
    )],
    pool_status[, .(
      cohort, check_id = "CONTRAST_RUBIN_COVARIANCE_FULL_RANK",
      observed = as.character(contrast_covariance_rank), expected = "2",
      status = fifelse(contrast_covariance_rank == 2L, "PASS", "FAIL"), critical = 1L
    )],
    data.table(
      cohort = "ALL", check_id = "ALL_40_COMPONENT_OPTIMIZERS_AND_COVARIANCES_ACCEPTED",
      observed = paste0(sum(component_metrics$status == "ACCEPTED"), "/", nrow(component_metrics)),
      expected = "40/40 with finite likelihood, finite PD Hessian and complete model keys",
      status = ifelse(component_diagnostics_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "ALL_ALLOWED_EDGE_POOLED_RESULTS_COMPLETE_FINITE",
      observed = paste0("rows=", nrow(all_hr)), expected = "32 complete finite Rubin rows",
      status = ifelse(all_edge_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "FOUR_TARGETS_PER_COHORT_EXACT",
      observed = paste(target_hr[, .N, by = cohort][order(cohort), paste0(cohort, "=", N)], collapse = "; "),
      expected = "ELSA=4; HRS=4",
      status = ifelse(target_identity_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "TARGET_RUBIN_ROWS_COMPLETE_FINITE_ORDERED",
      observed = paste0("rows=", nrow(target_hr)), expected = "8 rows; m=20; finite estimates/variances/P/CI",
      status = ifelse(target_pool_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "TWO_CONTRASTS_PER_COHORT_COMPLETE_FINITE_BH",
      observed = paste0("rows=", nrow(contrasts)), expected = "4 rows; m=20; finite estimates/P; exact within-cohort two-test BH",
      status = ifelse(contrast_pool_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "JOINT_TESTS_FULL_RANK_COMPLETE_FINITE",
      observed = paste0("rows=", nrow(joint_tests)), expected = "4 rows; correct 4df/2df ranks; finite statistics/P",
      status = ifelse(joint_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "PER_IMPUTATION_EXPOSURE_OUTPUT_COMPLETE_FINITE",
      observed = paste0("rows=", nrow(per_imputation)), expected = "640 rows; exact cohort-imputation-edge grid",
      status = ifelse(per_imputation_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "PER_IMPUTATION_TARGET_COVARIANCE_COMPLETE_SYMMETRIC",
      observed = per_imputation_covariance_audit$observed,
      expected = "640 rows; 40 complete symmetric 4x4 blocks; nonnegative diagonal",
      status = ifelse(per_imputation_covariance_audit$ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "RUBIN_TARGET_COVARIANCE_COMPLETE_SYMMETRIC",
      observed = pooled_covariance_audit$observed,
      expected = "96 rows; 6 complete symmetric 4x4 within/between/total blocks",
      status = ifelse(pooled_covariance_audit$ok && pooled_covariance_components_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "RUBIN_SCALAR_FULL_MATRIX_INTERNAL_CONSISTENCY",
      observed = paste(
        rubin_scalar_matrix_consistency[order(cohort),
          paste0(cohort, "=", format(maximum_absolute_error, scientific = TRUE))],
        collapse = "; "
      ),
      expected = "maximum absolute error <=1e-12 in each cohort",
      status = ifelse(rubin_scalar_matrix_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "CACHE_INVENTORY_40_MATCH_CURRENT_SPEC",
      observed = paste0(sum(pool_inventory$status == "ACCEPTED_MATCHING_SPEC"), "/", nrow(pool_inventory)),
      expected = "40/40 exact current-spec caches", status = ifelse(cache_inventory_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "FORMAL_MI20_MAPPING_40_EXACT",
      observed = paste0(sum(formal_sample_mapping_audit$mapping_status == "PASS"), "/", nrow(formal_sample_mapping_audit)),
      expected = "40/40 exact mappings", status = ifelse(mapping_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "CANONICAL_INPUT_AND_STRUCTURE_GATES_PASS",
      observed = paste0(
        "input_fail=", input_qc[critical == 1L & status != "PASS", .N],
        "; structure_fail=", model_structure_qc[critical == 1L & status != "PASS", .N]
      ),
      expected = "0;0", status = ifelse(input_gate_ok && structure_gate_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "FROZEN_INPUT_FINGERPRINT_SET_COMPLETE",
      observed = paste0("rows=", nrow(input_fingerprints)), expected = "10 rows including Q/state/contract/script",
      status = ifelse(fingerprint_gate_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "MODEL_SPECIFICATION_EDGE_SET_COMPLETE",
      observed = paste0("rows=", nrow(model_specification)), expected = "32 allowed edges with accepted pilot provenance",
      status = ifelse(model_specification_ok, "PASS", "FAIL"), critical = 1L
    ),
    data.table(
      cohort = "ALL", check_id = "Q_INTENSITY_SUMMARY_COMPLETE_FINITE",
      observed = paste0("rows=", nrow(q_summary)), expected = "145 HRS/ELSA matrix cells, each m=20",
      status = ifelse(q_summary_ok, "PASS", "FAIL"), critical = 1L
    ),
    mi_compatibility[, .(
      cohort, check_id = "FROZEN_MI_AUXILIARY_CONGENIALITY_RECORDED",
      observed = as.character(predictor_manifest_has_disease_specific_auxiliary), expected = "0 with explicit limitation",
      status = "PASS_WITH_LIMITATION", critical = 0L
    )],
    sparse_audit[thin_baseline_edge_warning == 1L, .(
      cohort, check_id = paste0("THIN_EDGE_", transition), observed = as.character(direct_events),
      expected = "informational: <25 direct observed transitions", status = "PASS_WITH_WARNING", critical = 0L
    )]
  ), fill = TRUE)
  atomic_fwrite(final_qc, file.path(qc_root, "disease_specific_final_inference_qc.csv"))
  if (final_qc[critical == 1L & status == "FAIL", .N]) {
    write_delivery_status(
      "REJECTED_FINAL_INFERENCE_QC",
      "All claim-bearing result, figure and report artifacts were withheld because at least one critical final QC gate failed"
    )
    stop("Final disease-specific inference QC failed; figures/report not generated", call. = FALSE)
  }

  # Claim-bearing artifacts are published only after every current-run MI20,
  # covariance, mapping, structure and numerical gate above has passed.  If a
  # later graphics/report step fails, delivery status remains explicitly
  # IN_PROGRESS_NOT_FOR_INFERENCE rather than endorsing a partial publication.
  atomic_fwrite(model_specification, file.path(result_root, "disease_specific_model_specification.csv"))
  atomic_fwrite(all_hr, file.path(result_root, "disease_specific_all_edge_insomnia_hr_mi20.csv"))
  atomic_fwrite(target_hr, file.path(result_root, "disease_specific_target_edge_hr_mi20.csv"))
  atomic_fwrite(contrasts, file.path(result_root, "disease_specific_heart_vs_stroke_contrasts_mi20.csv"))
  atomic_fwrite(joint_tests, file.path(result_root, "disease_specific_joint_tests_mi20.csv"))
  atomic_fwrite(target_covariance, file.path(result_root, "disease_specific_target_rubin_covariance.csv"))
  atomic_fwrite(
    per_imputation_target_covariance,
    file.path(result_root, "disease_specific_per_imputation_target_covariance.csv.gz")
  )
  atomic_fwrite(
    per_imputation,
    file.path(result_root, "disease_specific_per_imputation_exposure_estimates.csv.gz")
  )
  atomic_fwrite(q_summary, file.path(result_root, "disease_specific_q_intensity_mi_summary.csv"))
  atomic_fwrite(component_metrics, file.path(qc_root, "disease_specific_mi20_component_fit_metrics.csv"))
  atomic_save_rds(
    lapply(pooled, `[[`, "rubin_matrices"),
    file.path(result_root, "disease_specific_target_rubin_matrices.rds")
  )

  target_figure_data <- copy(target_hr)
  target_figure_data[, row_label := paste0(
    cohort, "  |  ", origin_state_label, " -> ", target_outcome
  )]
  target_figure_data[, row_order := match(cohort, c("HRS", "ELSA")) * 10L + target_order]
  target_figure_data[, row_label := factor(row_label, levels = rev(row_label[order(row_order)]))]
  target_figure_data[, result_text := sprintf(
    "%.2f (%.2f-%.2f); P=%s", hazard_ratio, conf_low, conf_high, format_p(p_value)
  )]
  atomic_fwrite(
    target_figure_data,
    file.path(source_data_root, "figure_ds1_target_transition_forest_source.csv")
  )

  contrast_figure_data <- copy(contrasts)
  contrast_figure_data[, row_label := paste0(cohort, "  |  ", origin_state_label)]
  contrast_figure_data[, row_order := match(cohort, c("HRS", "ELSA")) * 10L + match(origin_state_label, c("No disease", "Hypertension only"))]
  contrast_figure_data[, row_label := factor(row_label, levels = rev(row_label[order(row_order)]))]
  contrast_figure_data[, result_text := sprintf(
    "%.2f (%.2f-%.2f); P=%s", hr_ratio_heart_vs_stroke, conf_low, conf_high, format_p(p_value)
  )]
  atomic_fwrite(
    contrast_figure_data,
    file.path(source_data_root, "figure_ds2_heart_vs_stroke_contrast_source.csv")
  )
  figure_participant_note <- paste(
    formal_sample_mapping_audit[imputation == 1L][order(cohort),
      paste0(cohort, " n=", model_sample_participants)],
    collapse = "; "
  )

  base_theme <- ggplot2::theme_classic(base_family = "Arial", base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "bold", margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "#374151", margin = ggplot2::margin(b = 7)),
      plot.caption = ggplot2::element_text(size = 6.5, colour = "#4B5563", hjust = 0, margin = ggplot2::margin(t = 7)),
      axis.title = ggplot2::element_text(size = 8.5),
      axis.text = ggplot2::element_text(size = 7.5, colour = "#111827"),
      legend.title = ggplot2::element_blank(), legend.position = "top",
      legend.text = ggplot2::element_text(size = 7.5),
      panel.grid.major.x = ggplot2::element_line(colour = "#E5E7EB", linewidth = 0.3),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )

  target_plot <- ggplot2::ggplot(
    target_figure_data,
    ggplot2::aes(x = hazard_ratio, y = row_label, colour = target_outcome)
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, colour = "#6B7280") +
    ggplot2::geom_segment(
      ggplot2::aes(x = conf_low, xend = conf_high, yend = row_label),
      linewidth = 0.8, lineend = "round"
    ) +
    ggplot2::geom_point(size = 2.4, shape = 15) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c("Heart disease" = "#2166AC", "Stroke" = "#B2182B")) +
    ggplot2::labs(
      title = "Disease-specific target transitions",
      subtitle = "Adjusted M4-MI20 panel Markov estimates per one additional insomnia symptom",
      x = "Hazard ratio (95% CI, logarithmic scale)", y = NULL,
      caption = paste(
        "Points are MI20 Rubin-pooled hazard ratios; horizontal lines are 95% finite-imputation t confidence intervals; participant is the independent unit.",
        paste0("Analytic participants: ", figure_participant_note, "."),
        "Post-hoc secondary extension; four target edges were frozen before these new fits.",
        "HRS stroke includes TIA; ELSA stroke does not explicitly include TIA and its W4-W10 model is non-fatal.",
        sep = "\n"
      )
    ) + base_theme

  contrast_plot <- ggplot2::ggplot(
    contrast_figure_data,
    ggplot2::aes(x = hr_ratio_heart_vs_stroke, y = row_label, colour = origin_state_label)
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, colour = "#6B7280") +
    ggplot2::geom_segment(
      ggplot2::aes(x = conf_low, xend = conf_high, yend = row_label),
      linewidth = 0.8, lineend = "round"
    ) +
    ggplot2::geom_point(size = 2.4, shape = 16) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c("No disease" = "#0072B2", "Hypertension only" = "#D55E00")) +
    ggplot2::labs(
      title = "Within-origin comparison of heart-disease and stroke associations",
      subtitle = "HR ratio = exp(beta[heart disease] - beta[stroke]) = HR[heart disease] / HR[stroke]",
      x = "Heart-disease versus stroke HR ratio (95% CI, logarithmic scale)", y = NULL,
      caption = paste(
        "Points are MI20 Rubin-pooled HR ratios; horizontal lines are 95% finite-imputation t confidence intervals; participant is the independent unit.",
        paste0("Analytic participants: ", figure_participant_note, "."),
        "Ratios above 1 indicate a larger insomnia coefficient for broad heart disease than stroke.",
        "A non-significant contrast does not establish equivalence; phenotype definitions differ between cohorts.",
        sep = "\n"
      )
    ) + base_theme

  target_paths <- atomic_export_plot(target_plot, "figure_ds1_target_transition_forest", 7.2, 5.2)
  contrast_paths <- atomic_export_plot(contrast_plot, "figure_ds2_heart_vs_stroke_contrasts", 7.2, 4.2)
  figure_manifest <- rbindlist(list(
    data.table(
      figure_id = "DS1", title = "Disease-specific target transition forest",
      artifact_path = normalizePath(target_paths, winslash = "/"),
      source_data = normalizePath(file.path(source_data_root, "figure_ds1_target_transition_forest_source.csv"), winslash = "/"),
      backend = "R ggplot2", analysis_role = "post-hoc secondary extension"
    ),
    data.table(
      figure_id = "DS2", title = "Heart-versus-stroke coefficient contrasts",
      artifact_path = normalizePath(contrast_paths, winslash = "/"),
      source_data = normalizePath(file.path(source_data_root, "figure_ds2_heart_vs_stroke_contrast_source.csv"), winslash = "/"),
      backend = "R ggplot2", analysis_role = "post-hoc secondary extension"
    )
  ), fill = TRUE)
  figure_manifest[, automated_qa_status := "PENDING_INDEPENDENT_POST_EXPORT_AUDIT; statistics are frozen before graphics QA"]
  atomic_fwrite(figure_manifest, file.path(figure_root, "disease_specific_figure_manifest.csv"))

  report_target <- copy(target_hr)[, .(
    Cohort = cohort,
    Origin = origin_state_label,
    Outcome = target_outcome,
    `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", hazard_ratio, conf_low, conf_high),
    `P` = format_p(p_value),
    `Rubin missing-variance fraction` = sprintf("%.3f", rubin_missing_variance_fraction)
  )]
  report_contrasts <- copy(contrasts)[, .(
    Cohort = cohort,
    Origin = origin_state_label,
    `Heart/stroke HR ratio (95% CI)` = sprintf("%.2f (%.2f-%.2f)", hr_ratio_heart_vs_stroke, conf_low, conf_high),
    `P` = format_p(p_value),
    `BH P (two planned contrasts within cohort)` = format_p(p_bh_two_contrasts_within_cohort)
  )]
  report_joint <- copy(joint_tests)[, .(
    Cohort = cohort, Test = test_id,
    `Wald chi-square` = sprintf("%.3f", statistic), df,
    `P` = format_p(p_value),
    Rank = covariance_rank_status
  )]
  report_diagnostics <- component_metrics[, .(
    `Accepted MI fits` = paste0(sum(status == "ACCEPTED"), "/", .N),
    `Elapsed minutes, median (range)` = sprintf(
      "%.1f (%.1f-%.1f)", median(elapsed_seconds) / 60, min(elapsed_seconds) / 60, max(elapsed_seconds) / 60
    ),
    `Min Hessian eigenvalue` = format(min(hessian_min_eigenvalue), scientific = TRUE, digits = 3),
    `Max Hessian condition number` = format(max(hessian_condition_number), scientific = TRUE, digits = 3),
    `Fits with warnings` = sum(warning_count > 0L)
  ), by = .(Cohort = cohort)]
  report_schemes <- unique(model_specification[, .(
    Cohort = cohort,
    `Selected exposure-constraint scheme` = selected_constraint_scheme,
    `Pilot cache validation` = pilot_cache_validation
  )])
  setorder(report_schemes, Cohort)
  report_initialization <- component_metrics[, .(
    `MI components` = .N
  ), by = .(Cohort = cohort, `Recorded initialization source` = warm_start_source)]
  setorder(report_initialization, Cohort, `Recorded initialization source`)
  thin_report <- sparse_audit[thin_baseline_edge_warning == 1L, .(
    Cohort = cohort, Transition = transition, `Direct observed transitions` = direct_events,
    `Exposure constraint group` = exposure_constraint_group
  )]
  report_lines <- c(
    "# Disease-specific HRS/ELSA panel Markov extension",
    "",
    paste0("Generated: ", run_timestamp_utc),
    "",
    "## Analysis status",
    "",
    paste(
      "This is a **post-hoc secondary extension** performed after the primary composite-state analysis.",
      "The four disease-specific target edges and two within-origin heart-versus-stroke contrasts were frozen",
      "in the secondary-extension contract before the models reported here were fitted; they were not preregistered confirmatory aims."
    ),
    "",
    "All reported pooled estimates require 20/20 matching-specification component fits, optimizer convergence code 0,",
    "finite standard errors and covariance matrices, positive-definite Hessians, and full-rank Rubin target/contrast covariance blocks.",
    "No pseudoinverse or reduced degrees of freedom was used to rescue a rank-deficient joint test.",
    "",
    "## Model",
    "",
    paste(
      "The HRS model contains D0-D7 living component combinations plus absorbing D8 death.",
      "The ELSA post-hoc secondary extension contains D0-D7 living states at W4-W10 and no death state.",
      "Allowed Q edges add one disease component at a time; observed same-wave multi-component jumps were retained as interval endpoints",
      "and integrated through compatible latent paths by the continuous-time panel likelihood rather than assigned an artificial order."
    ),
    "",
    paste(
      "Insomnia symptom score (0-3) has four transition-specific target coefficients:",
      "D0->D2 heart disease, D0->D3 stroke, D1->D4 hypertension plus heart disease, and D1->D5 hypertension plus stroke.",
      "For each cohort, the non-target living-effect constraint scheme was selected by the complete-case pilot:",
      "the etiologic grouping when accepted, otherwise the predeclared collapsed fallback;",
      "HRS death effects are common across death edges. Each M4 nuisance covariate coefficient is shared across all allowed Q edges."
    ),
    "",
    markdown_table(report_schemes),
    "",
    paste(
      "The pilot is used only to screen the constraint scheme and provide optional starting values; no pilot estimate enters Rubin pooling.",
      "A validated legacy pilot may use a polynomial rather than treatment-contrast basis for ordered-looking category codes,",
      "but only because both are full-rank parameterizations of the same categorical nuisance space and its exact legacy hash must match.",
      "All 40 formal MI components use unordered treatment contrasts and independently pass the current likelihood, covariance and cache gates."
    ),
    "",
    markdown_table(report_initialization),
    "",
    paste(
      "Initialization provenance is recorded for reproducibility but is not part of the estimand:",
      "a formal component may start from an accepted pilot or the preceding accepted MI component, while acceptance is evaluated independently under the frozen formal specification."
    ),
    "",
    "M4 adjusts for age, sex, race, education, wealth, marital status, smoking, drinking, inactivity, BMI, diabetes, lung disease and frequent pain.",
    "Scalar intervals and contrasts use finite-imputation Rubin t inference; joint tests use the full Rubin total covariance with an explicitly labelled asymptotic Wald chi-square.",
    "Only the two planned heart-versus-stroke contrasts within each cohort form a BH multiplicity family; target-edge and joint-test P values are reported unadjusted.",
    "",
    "## Target-edge estimates",
    "",
    markdown_table(report_target),
    "",
    "## Heart-versus-stroke contrasts",
    "",
    markdown_table(report_contrasts),
    "",
    "A ratio above 1 denotes a larger insomnia coefficient for broad heart disease than for stroke within the same origin state.",
    "A non-significant contrast is not evidence of equivalence.",
    "",
    "## Joint tests",
    "",
    markdown_table(report_joint),
    "",
    "## Component-fit diagnostics",
    "",
    markdown_table(report_diagnostics),
    "",
    "## Sparse-edge audit",
    "",
    if (nrow(thin_report)) markdown_table(thin_report) else "No allowed edge had fewer than 25 directly observed transitions.",
    "",
    "Edges with 15-24 direct observed transitions passed the hard structural-support gate but are explicitly flagged as thin;",
    "exposure effects for non-target edges are constrained to supported groups and are not interpreted as disease-specific target estimates.",
    "",
    "## Phenotype and missing-data limitations",
    "",
    paste(
      "HRS broad heart disease cannot be reliably separated into myocardial infarction, angina and heart failure in this harmonized panel;",
      "HRS stroke includes TIA. ELSA stroke does not explicitly include TIA, and the W4-W10 secondary extension is non-fatal.",
      "Accordingly, cross-cohort agreement is evaluated at the broad heart-disease/stroke level, not as identical clinical phenotypes."
    ),
    "",
    paste(
      "The frozen primary-delivery M4-MI20 covariate draws included composite-state outcome auxiliaries but not hypertension, heart disease and stroke",
      "as separate imputation auxiliaries. Their reuse preserves the primary record mapping and analysis convention, but it is not fully congenial",
      "to the disease-specific outcome decomposition. Results must therefore be described as a frozen-MI secondary sensitivity/extension, and this limitation cannot be hidden."
    ),
    "",
    paste(
      "The authoritative Q-edge preparation file contains the legacy source descriptor 'confirmatory_disease_specific' for target rows.",
      "That descriptor is retained only as input provenance and has no inferential meaning here; this report and the model specification override it with",
      "'planned target within post-hoc secondary extension; frozen before these new fits'."
    ),
    "",
    "## Reproducibility outputs",
    "",
    "Component caches are accepted only when their current data/model specification hash matches. Per-imputation estimates and full target covariance matrices,",
    "Rubin within/between/total matrices, Q summaries, model specification, input fingerprints, and all QC gates are exported with this report.",
    "The figure manifest is marked PENDING_INDEPENDENT_POST_EXPORT_AUDIT at generation; independent collision/text/export QA is required before delivery status can be upgraded."
  )
  write_text(report_lines, file.path(report_root, "disease_specific_model_report.md"))
  write_text(capture.output(sessionInfo()), file.path(qc_root, "disease_specific_session_info.txt"))

  artifact_roots <- c(result_root, figure_root, report_root, qc_root, model_root, log_root)
  artifact_files <- unique(unlist(lapply(artifact_roots, function(path) {
    list.files(path, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  }), use.names = FALSE))
  artifact_files <- artifact_files[!grepl(
    "disease_specific_(artifact_manifest|delivery_status)\\.csv$", artifact_files
  )]
  project_abs <- normalizePath(project_root, winslash = "/")
  artifact_manifest <- rbindlist(lapply(artifact_files, function(path) {
    normalized <- normalizePath(path, winslash = "/")
    info <- file.info(path)
    data.table(
      relative_path = substring(normalized, nchar(project_abs) + 2L),
      bytes = as.numeric(info$size), modified_utc = format(info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      sha256 = sha256_file(path)
    )
  }), fill = TRUE)
  setorder(artifact_manifest, relative_path)
  atomic_fwrite(artifact_manifest, file.path(qc_root, "disease_specific_artifact_manifest.csv"))
  write_delivery_status(
    "GENERATED_PENDING_INDEPENDENT_FIGURE_QA",
    "All 40 formal components and final inference gates passed; claim-bearing artifacts were generated, but figure delivery remains pending independent post-export QA"
  )
  cat("Disease-specific MI20 modeling, pooling, figures and report completed.\n")
}
