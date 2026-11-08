#!/usr/bin/env Rscript

# Post-hoc secondary extension: formal comparison of three insomnia symptom
# domains in the HRS/ELSA continuous-time multistate models.  Contrasts and
# their hierarchy were frozen before these joint models were fitted.
#
# IMPORTANT: inputs are read exclusively from the immutable final-delivery
# copy under outputs/hrs_elsa_multistate_final_20260909/data_multistate.

suppressPackageStartupMessages({
  library(data.table)
  library(msm)
  library(mice)
  library(digest)
})

options(stringsAsFactors = FALSE, width = 180)

script_version <- "2026-09-09-joint-domains-v1"
mi_seed <- 20260908L
default_m <- 20L
default_maxit <- 10L
max_iterations <- 1000L

# Keep operational paths relative when launched from the project root.  This
# avoids Windows/R locale conversion failures on the Chinese project folder,
# while still pointing to the exact final-delivery copy.
project_root <- "."
delivery_root <- file.path("outputs", "hrs_elsa_multistate_final_20260909")
input_root <- file.path(delivery_root, "data_multistate")
default_output_root <- file.path("outputs", "hrs_elsa_multistate_nextsteps_20260909")
output_root <- Sys.getenv("JOINT_DOMAIN_OUTPUT_ROOT", unset = default_output_root)

result_root <- file.path(output_root, "results", "symptom_domains")
figure_root <- file.path(output_root, "figures", "symptom_domains")
report_root <- file.path(output_root, "reports", "symptom_domains")
qc_root <- file.path(output_root, "qc", "symptom_domains")
model_root <- file.path(result_root, "model_cache")
mi_root <- file.path(result_root, "mi_cache")
log_root <- file.path(result_root, "logs")
source_root <- file.path(figure_root, "source_data")
for (path in c(result_root, figure_root, report_root, qc_root, model_root, mi_root, log_root, source_root)) {
  if (!dir.exists(path)) {
    try(dir.create(path, recursive = TRUE, showWarnings = FALSE), silent = TRUE)
    if (!dir.exists(path)) stop("Unable to create output directory: ", path)
  }
}

if (!dir.exists(input_root)) stop("Final-delivery input directory is missing: ", input_root)
if (grepl("ga_stage", input_root, fixed = TRUE)) {
  stop("Refusing to use a staging-data path")
}

env_flag <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = if (default) "1" else "0")))
  value %in% c("1", "true", "yes", "y", "on")
}
env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value)) as.integer(default) else value
}

run_stage <- tolower(Sys.getenv("JOINT_DOMAIN_STAGE", unset = "full"))
if (!run_stage %in% c("full", "mi", "fit", "pool")) stop("JOINT_DOMAIN_STAGE must be full, mi, fit, or pool")
analysis_regex <- Sys.getenv("JOINT_DOMAIN_ANALYSIS_REGEX", unset = ".")
mi_m <- env_int("JOINT_DOMAIN_M", default_m)
mi_maxit <- env_int("JOINT_DOMAIN_MI_MAXIT", default_maxit)
imp_start <- max(1L, env_int("JOINT_DOMAIN_IMP_START", 1L))
imp_end <- min(mi_m, max(imp_start, env_int("JOINT_DOMAIN_IMP_END", mi_m)))
force_mi <- env_flag("JOINT_DOMAIN_FORCE_MI", FALSE)
force_fit <- env_flag("JOINT_DOMAIN_FORCE_REFIT", FALSE)
enable_hrs4 <- env_flag("JOINT_DOMAIN_ENABLE_HRS4", TRUE)
use_warm_start <- env_flag("JOINT_DOMAIN_WARM_START", TRUE)

domain_vars <- c(
  "insomnia_onset_sync",
  "insomnia_maintenance_sync",
  "insomnia_nonrestorative_sync"
)
domain_labels <- c(
  insomnia_onset_sync = "Sleep-onset difficulty",
  insomnia_maintenance_sync = "Sleep-maintenance difficulty",
  insomnia_nonrestorative_sync = "Nonrestorative sleep"
)
domain_labels_zh <- c(
  insomnia_onset_sync = "入睡困难",
  insomnia_maintenance_sync = "睡眠维持困难",
  insomnia_nonrestorative_sync = "非恢复性睡眠"
)
primary_m4_adjustment <- c(
  "age_bl_5", "female_bl", "race_bl", "education_level_bl",
  "wealth_quintile_bl", "married_bl", "current_smoker_bl",
  "current_drinker_bl", "inactive_bl", "bmi_bl_5", "diabetes_bl",
  "lung_disease_bl", "frequent_pain_bl"
)

state_labels <- c(
  S0 = "No target disease",
  S1 = "Hypertension only",
  S2 = "One CVD without hypertension",
  S3 = "Multimorbidity",
  S4 = "Death"
)

locf_down <- function(x) {
  if (!length(x)) return(x)
  last <- NA
  for (i in seq_along(x)) {
    if (!is.na(x[i])) last <- x[i] else if (!is.na(last)) x[i] <- last
  }
  x
}

file_sha256 <- function(path) digest::digest(file = path, algo = "sha256")
table_fingerprint <- function(x, columns) {
  columns <- unique(intersect(columns, names(x)))
  digest::digest(as.data.table(x)[, ..columns], algo = "xxhash64", serialize = TRUE)
}

allowed_edges <- function(n_states) {
  if (n_states == 5L) {
    out <- data.table(
      from_index = c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L),
      to_index = c(2L, 3L, 5L, 4L, 5L, 4L, 5L, 5L)
    )
  } else if (n_states == 4L) {
    out <- data.table(
      from_index = c(1L, 1L, 2L, 3L),
      to_index = c(2L, 3L, 4L, 4L)
    )
  } else stop("Only four- and five-state protocol graphs are supported")
  out[, `:=`(from_state = from_index - 1L, to_state = to_index - 1L)]
  out[, `:=`(
    transition = paste0("S", from_state, "->S", to_state),
    death_transition = as.integer(n_states == 5L & to_index == 5L)
  )]
  out[]
}

analysis_configs <- list(
  hrs_5state_joint_domains = list(
    analysis_id = "hrs_5state_joint_domains", cohort = "HRS", input = "hrs_msm_panel.rds",
    n_states = 5L, rule = "all", death_exact_state = 5L, mi_id = "hrs_joint_domains"
  ),
  elsa_4state_w4_w10_joint_domains = list(
    analysis_id = "elsa_4state_w4_w10_joint_domains", cohort = "ELSA", input = "elsa_msm_panel.rds",
    n_states = 4L, rule = "elsa_nonfatal", death_exact_state = NULL, mi_id = "elsa_joint_domains"
  ),
  hrs_4state_death_censored_joint_domains = list(
    analysis_id = "hrs_4state_death_censored_joint_domains", cohort = "HRS", input = "hrs_msm_panel.rds",
    n_states = 4L, rule = "hrs_nonfatal", death_exact_state = NULL, mi_id = "hrs_joint_domains"
  )
)
if (!enable_hrs4) analysis_configs$hrs_4state_death_censored_joint_domains <- NULL

derive_synchronised_domains <- function(x, cohort) {
  x <- copy(as.data.table(x))
  setorder(x, id, event_time, wave, state, na.last = TRUE)
  if (cohort == "HRS") {
    x[, insomnia_onset_observed_joint := fifelse(fallslp %in% 1:3, as.integer(fallslp == 1), NA_integer_)]
    x[, insomnia_maintenance_observed_joint := fifelse(wakent %in% 1:3, as.integer(wakent == 1), NA_integer_)]
    x[, insomnia_nonrestorative_observed_joint := fifelse(rested %in% 1:3, as.integer(rested == 3), NA_integer_)]
  } else {
    needed <- c("insomnia_onset_observed", "insomnia_maintenance_observed", "insomnia_nonrestorative_observed")
    if (length(setdiff(needed, names(x)))) stop("ELSA observed symptom columns are missing")
    x[, insomnia_onset_observed_joint := as.integer(insomnia_onset_observed)]
    x[, insomnia_maintenance_observed_joint := as.integer(insomnia_maintenance_observed)]
    x[, insomnia_nonrestorative_observed_joint := as.integer(insomnia_nonrestorative_observed)]
  }
  obs <- c(
    "insomnia_onset_observed_joint", "insomnia_maintenance_observed_joint",
    "insomnia_nonrestorative_observed_joint"
  )
  x[, insomnia_triplet_complete_derived := as.integer(complete.cases(.SD)), .SDcols = obs]
  x[, insomnia_onset_sync := fifelse(insomnia_triplet_complete_derived == 1L, insomnia_onset_observed_joint, NA_integer_)]
  x[, insomnia_maintenance_sync := fifelse(insomnia_triplet_complete_derived == 1L, insomnia_maintenance_observed_joint, NA_integer_)]
  x[, insomnia_nonrestorative_sync := fifelse(insomnia_triplet_complete_derived == 1L, insomnia_nonrestorative_observed_joint, NA_integer_)]
  x[, (domain_vars) := lapply(.SD, locf_down), by = id, .SDcols = domain_vars]
  x[, insomnia_score_sync := fifelse(
    complete.cases(.SD),
    insomnia_onset_sync + insomnia_maintenance_sync + insomnia_nonrestorative_sync,
    NA_integer_
  ), .SDcols = domain_vars]
  async_vars <- c("insomnia_onset_primary", "insomnia_maintenance_primary", "insomnia_nonrestorative_primary")
  x[, insomnia_async_sum := fifelse(
    complete.cases(.SD),
    insomnia_onset_primary + insomnia_maintenance_primary + insomnia_nonrestorative_primary,
    NA_integer_
  ), .SDcols = async_vars]
  x[, sync_score_disagrees_primary := as.integer(
    !is.na(insomnia_score_sync) & !is.na(insomnia_score_primary) & insomnia_score_sync != insomnia_score_primary
  )]
  x[, async_score_disagrees_primary := as.integer(
    !is.na(insomnia_async_sum) & !is.na(insomnia_score_primary) & insomnia_async_sum != insomnia_score_primary
  )]
  x[, any_sync_async_domain_difference := as.integer(
    (!is.na(insomnia_onset_sync) & !is.na(insomnia_onset_primary) & insomnia_onset_sync != insomnia_onset_primary) |
      (!is.na(insomnia_maintenance_sync) & !is.na(insomnia_maintenance_primary) & insomnia_maintenance_sync != insomnia_maintenance_primary) |
      (!is.na(insomnia_nonrestorative_sync) & !is.na(insomnia_nonrestorative_primary) & insomnia_nonrestorative_sync != insomnia_nonrestorative_primary)
  )]
  x[]
}

load_base_panels <- function() {
  out <- list()
  for (cohort in c("HRS", "ELSA")) {
    filename <- paste0(tolower(cohort), "_msm_panel.rds")
    path <- file.path(input_root, filename)
    if (!file.exists(path)) stop("Required final-delivery input is missing: ", path)
    raw <- readRDS(path)
    out[[cohort]] <- derive_synchronised_domains(raw, cohort)
  }
  out
}

write_sync_audit <- function(panels) {
  summaries <- list()
  waves <- list()
  sources <- list()
  correlations <- list()
  correlation_summaries <- list()
  for (cohort in names(panels)) {
    x <- panels[[cohort]]
    cohort_label <- cohort
    summaries[[cohort]] <- data.table(
      cohort = cohort,
      n_rows = nrow(x), n_participants = uniqueN(x$id),
      complete_triplet_rows = sum(x$insomnia_triplet_complete_derived == 1L),
      rows_without_prior_complete_triplet = sum(!complete.cases(x[, ..domain_vars])),
      participants_without_sync_at_any_row = uniqueN(x[!complete.cases(x[, ..domain_vars]), id]),
      async_score_mismatch_rows = sum(x$async_score_disagrees_primary == 1L),
      async_score_mismatch_participants = uniqueN(x[async_score_disagrees_primary == 1L, id]),
      sync_score_mismatch_rows = sum(x$sync_score_disagrees_primary == 1L),
      sync_score_mismatch_participants = uniqueN(x[sync_score_disagrees_primary == 1L, id]),
      sync_async_domain_difference_rows = sum(x$any_sync_async_domain_difference == 1L),
      sync_async_domain_difference_participants = uniqueN(x[any_sync_async_domain_difference == 1L, id]),
      assessment_flag_disagreement_rows = if ("insomnia_assessment_complete" %in% names(x)) {
        sum(!is.na(x$insomnia_assessment_complete) & x$insomnia_assessment_complete != x$insomnia_triplet_complete_derived)
      } else NA_integer_
    )
    waves[[cohort]] <- x[, .(
      n_rows = .N,
      n_participants = uniqueN(id),
      complete_triplet_rows = sum(insomnia_triplet_complete_derived == 1L),
      rows_without_prior_complete_triplet = sum(!complete.cases(.SD)),
      async_score_mismatch_rows = sum(async_score_disagrees_primary == 1L),
      sync_score_mismatch_rows = sum(sync_score_disagrees_primary == 1L),
      sync_async_domain_difference_rows = sum(any_sync_async_domain_difference == 1L)
    ), by = .(cohort, wave), .SDcols = domain_vars]
    path <- file.path(input_root, paste0(tolower(cohort), "_msm_panel.rds"))
    sources[[cohort]] <- data.table(
      cohort = cohort, source_role = "final-delivery multistate panel",
      path = gsub("\\", "/", path, fixed = TRUE),
      sha256 = file_sha256(path), bytes = file.info(path)$size
    )
    correlation_sample <- x[complete.cases(x[, ..domain_vars])]
    correlation_matrix <- stats::cor(correlation_sample[, ..domain_vars])
    correlation_long <- as.data.table(as.table(correlation_matrix))
    setnames(correlation_long, c("domain_a", "domain_b", "correlation"))
    correlation_long[, `:=`(
      cohort = cohort_label,
      sample = "all synchronized panel rows"
    )]
    setcolorder(correlation_long, c("cohort", "sample", "domain_a", "domain_b", "correlation"))
    correlations[[cohort]] <- correlation_long
    off_diagonal <- correlation_matrix[row(correlation_matrix) != col(correlation_matrix)]
    eigenvalues <- eigen(correlation_matrix, symmetric = TRUE, only.values = TRUE)$values
    correlation_summaries[[cohort]] <- data.table(
      cohort = cohort_label,
      n_complete_synchronized_rows = nrow(correlation_sample),
      max_absolute_pairwise_correlation = max(abs(off_diagonal)),
      correlation_matrix_condition_number = max(eigenvalues) / min(eigenvalues)
    )
  }
  summary <- rbindlist(summaries, use.names = TRUE, fill = TRUE)
  by_wave <- rbindlist(waves, use.names = TRUE, fill = TRUE)
  source <- rbindlist(sources, use.names = TRUE, fill = TRUE)
  correlation <- rbindlist(correlations, use.names = TRUE, fill = TRUE)
  correlation_summary <- rbindlist(correlation_summaries, use.names = TRUE, fill = TRUE)
  fwrite(summary, file.path(qc_root, "symptom_domain_sync_audit_summary.csv"), na = "")
  fwrite(by_wave, file.path(qc_root, "symptom_domain_sync_audit_by_wave.csv"), na = "")
  fwrite(source, file.path(qc_root, "input_source_manifest.csv"), na = "")
  fwrite(correlation, file.path(qc_root, "symptom_domain_correlations.csv"), na = "")
  fwrite(correlation_summary, file.path(qc_root, "symptom_domain_correlation_summary.csv"), na = "")
  list(
    summary = summary, by_wave = by_wave, source = source,
    correlation = correlation, correlation_summary = correlation_summary
  )
}

prepare_analysis_data <- function(panel, config) {
  x <- copy(panel)
  x[, state_msm := as.integer(state) + 1L]
  x[, time_msm := as.numeric(event_time) - as.numeric(exposure_entry_time)]
  if (config$rule == "elsa_nonfatal") {
    x <- x[event_type == "visit" & wave %in% 4:10 & state_msm %in% 1:4]
  } else if (config$rule == "hrs_nonfatal") {
    x <- x[event_type == "visit" & state_msm %in% 1:4]
  } else {
    x <- x[state_msm %in% 1:5]
  }
  x <- x[!is.na(id) & is.finite(time_msm) & time_msm >= 0]
  x <- x[complete.cases(x[, ..domain_vars])]
  setorder(x, id, time_msm, state_msm)
  eligible <- x[, .N, by = id][N >= 2L, id]
  x <- x[id %in% eligible]
  if (config$cohort == "HRS") {
    x[, race_bl := factor(race_bl, levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other"))]
  } else {
    x[, race_bl := factor(race_bl, levels = c("White", "Non-White"))]
  }
  x[, education_level_bl := factor(education_level_bl, levels = 1:3)]
  x[, wealth_quintile_bl := factor(wealth_quintile_bl, levels = 1:5)]
  x[, age_bl_5 := (as.numeric(age_bl) - 65) / 5]
  x[, bmi_bl_5 := (as.numeric(bmi_bl) - 25) / 5]
  x[]
}

make_analysis_sets <- function(panels) {
  out <- list()
  for (config in analysis_configs) {
    if (!grepl(analysis_regex, config$analysis_id)) next
    out[[config$analysis_id]] <- prepare_analysis_data(panels[[config$cohort]], config)
  }
  out
}

edge_support <- function(x, config) {
  edges <- allowed_edges(config$n_states)
  d <- copy(x[, .(id, time_msm, state_msm, insomnia_onset_sync, insomnia_maintenance_sync, insomnia_nonrestorative_sync)])
  setorder(d, id, time_msm, state_msm)
  d[, next_state := shift(state_msm, type = "lead"), by = id]
  d[, next_time := shift(time_msm, type = "lead"), by = id]
  pairs <- d[!is.na(next_state) & next_time > time_msm]
  out <- copy(edges)
  out[, `:=`(
    risk_intervals = vapply(from_index, function(s) pairs[state_msm == s, .N], integer(1)),
    risk_participants = vapply(from_index, function(s) pairs[state_msm == s, uniqueN(id)], integer(1)),
    person_years = vapply(from_index, function(s) pairs[state_msm == s, sum(next_time - time_msm)], numeric(1)),
    events = mapply(function(a, b) pairs[state_msm == a & next_state == b, .N], from_index, to_index),
    analysis_id = config$analysis_id,
    n_records = nrow(x), n_participants = uniqueN(x$id)
  )]
  out[]
}

symptom_pattern_support <- function(x, config) {
  d <- copy(x[, c("id", "time_msm", "state_msm", domain_vars), with = FALSE])
  setorder(d, id, time_msm, state_msm)
  d[, next_state := shift(state_msm, type = "lead"), by = id]
  d[, next_time := shift(time_msm, type = "lead"), by = id]
  d <- d[!is.na(next_state) & next_time > time_msm]
  d[, pattern := paste0(insomnia_onset_sync, insomnia_maintenance_sync, insomnia_nonrestorative_sync)]
  out <- d[, .(
    risk_intervals = .N, risk_participants = uniqueN(id), person_years = sum(next_time - time_msm),
    observed_state_changes = sum(next_state != state_msm)
  ), by = .(origin_state = state_msm - 1L, pattern)]
  out[, analysis_id := config$analysis_id]
  setcolorder(out, c("analysis_id", "origin_state", "pattern"))
  out[]
}

make_mi_frame <- function(x, config) {
  setorder(x, id, time_msm, state_msm)
  baseline <- x[, .SD[1L], by = id, .SDcols = unique(c(primary_m4_adjustment, domain_vars, "state_msm", "weight_bl"))]
  auxiliary <- x[, .(
    followup_years = max(time_msm) - min(time_msm),
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
  for (name in binary_targets) set(frame, j = name, value = factor(frame[[name]], levels = 0:1))
  frame[, education_level_bl := ordered(education_level_bl, levels = 1:3)]
  frame[, wealth_quintile_bl := ordered(wealth_quintile_bl, levels = 1:5)]
  frame[, race_bl := droplevels(race_bl)]
  for (name in domain_vars) set(frame, j = name, value = as.integer(frame[[name]]))
  if (anyNA(frame[, ..domain_vars])) stop("Protected symptom domains are missing in the MI frame")
  if ("insomnia_score_primary" %in% names(frame)) stop("Deterministic total score must not enter the joint-domain MI frame")
  frame[]
}

mi_cache_paths <- function(mi_id) {
  prefix <- file.path(mi_root, paste0(mi_id, "__M4_joint_domains_MI", mi_m))
  list(
    prefix = prefix, mids = paste0(prefix, "__mids.rds"), long = paste0(prefix, "__baseline_long.rds"),
    manifest = paste0(prefix, "__cache_manifest.rds"), audit = paste0(prefix, "__audit.csv"),
    predictors = paste0(prefix, "__predictor_matrix.csv"), logged = paste0(prefix, "__logged_events.csv"),
    convergence = paste0(prefix, "__convergence.csv")
  )
}

prepare_mi <- function(x, config) {
  paths <- mi_cache_paths(config$mi_id)
  frame <- make_mi_frame(copy(x), config)
  seed <- mi_seed + ifelse(config$cohort == "HRS", 1L, 2L)
  signature <- list(
    script_version = script_version, mi_id = config$mi_id, m = mi_m, maxit = mi_maxit, seed = seed,
    ids = sort(as.character(frame$id)), frame_fingerprint = table_fingerprint(frame, names(frame)),
    input_sha256 = file_sha256(file.path(input_root, config$input))
  )
  if (!force_mi && file.exists(paths$manifest) && file.exists(paths$mids)) {
    old <- tryCatch(readRDS(paths$manifest), error = function(e) NULL)
    imp <- tryCatch(readRDS(paths$mids), error = function(e) NULL)
    if (!is.null(old) && identical(old$signature, signature) && !is.null(imp) && imp$m == mi_m) {
      message("Reusing MI cache: ", config$mi_id)
      return(list(mids = imp, frame = frame, signature = signature, paths = paths, cache_status = "reused"))
    }
  }
  methods <- mice::make.method(as.data.frame(frame))
  methods[] <- ""
  for (name in primary_m4_adjustment) {
    if (!anyNA(frame[[name]])) next
    value <- frame[[name]]
    methods[name] <- if (is.ordered(value)) "polr" else if (is.factor(value) && nlevels(value) == 2L) {
      "logreg"
    } else if (is.factor(value)) "polyreg" else "pmm"
  }
  predictors <- mice::make.predictorMatrix(as.data.frame(frame))
  predictors[,] <- 0L
  candidates <- setdiff(names(frame), c("id", "weight_bl"))
  candidates <- candidates[vapply(candidates, function(name) length(unique(frame[[name]][!is.na(frame[[name]])])) > 1L, logical(1))]
  for (target in names(methods)[methods != ""]) predictors[target, setdiff(candidates, target)] <- 1L
  protected <- c("id", domain_vars, "state_msm", "last_state_msm", "followup_years", "n_panel_records", "n_state_changes", "death_observed", "weight_bl")
  if (any(methods[intersect(protected, names(methods))] != "")) stop("A protected variable was assigned an imputation method")
  started <- Sys.time()
  imp <- mice::mice(
    data = as.data.frame(frame), m = mi_m, maxit = mi_maxit,
    method = methods, predictorMatrix = predictors, seed = seed,
    printFlag = FALSE, remove.collinear = FALSE, remove.constant = FALSE
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  completed_long <- as.data.table(mice::complete(imp, action = "long", include = FALSE))
  audit <- rbindlist(lapply(names(frame), function(name) data.table(
    mi_id = config$mi_id, variable = name,
    target_baseline_covariate = as.integer(name %in% primary_m4_adjustment),
    protected_from_imputation = as.integer(name %in% protected),
    method = unname(methods[name]), missing_before = sum(is.na(frame[[name]])),
    missing_fraction_before = mean(is.na(frame[[name]])),
    missing_after_across_all_imputations = sum(is.na(completed_long[[name]])),
    predictor_count = if (name %in% rownames(predictors)) sum(predictors[name, ] != 0L) else 0L
  )))
  predictor_long <- as.data.table(as.table(predictors))
  setnames(predictor_long, c("target", "predictor", "used"))
  logged <- if (is.null(imp$loggedEvents)) data.table() else as.data.table(imp$loggedEvents)
  convergence <- tryCatch(as.data.table(mice::convergence(imp)), error = function(e) data.table(error = conditionMessage(e)))
  saveRDS(imp, paths$mids, compress = "gzip")
  saveRDS(completed_long, paths$long, compress = "gzip")
  fwrite(audit, paths$audit, na = "")
  fwrite(predictor_long, paths$predictors, na = "")
  fwrite(logged, paths$logged, na = "")
  fwrite(convergence, paths$convergence, na = "")
  saveRDS(list(signature = signature, methods = methods, predictor_matrix = predictors, protected = protected, elapsed_seconds = elapsed), paths$manifest, compress = "gzip")
  message(sprintf("Prepared MI cache %s in %.1f seconds", config$mi_id, elapsed))
  list(mids = imp, frame = frame, signature = signature, paths = paths, cache_status = "created")
}

make_q_initial <- function(x, edges, n_states, minimum = 1e-4) {
  d <- copy(x[, .(id, time_msm, state_msm)])
  setorder(d, id, time_msm, state_msm)
  d[, next_state := shift(state_msm, type = "lead"), by = id]
  d[, next_time := shift(time_msm, type = "lead"), by = id]
  d <- d[!is.na(next_state) & next_time > time_msm]
  q <- matrix(0, n_states, n_states, dimnames = list(paste0("S", 0:(n_states - 1L)), paste0("S", 0:(n_states - 1L))))
  for (i in seq_len(nrow(edges))) {
    a <- edges$from_index[i]
    b <- edges$to_index[i]
    py <- d[state_msm == a, sum(next_time - time_msm)]
    ev <- d[state_msm == a & next_state == b, .N]
    value <- if (is.finite(py) && py > 0) ev / py else minimum
    q[a, b] <- max(minimum, value)
  }
  diag(q) <- 0
  if (n_states == 5L) q[5L, ] <- 0
  q
}

prepare_fit_sample <- function(x, completed) {
  imputed <- intersect(primary_m4_adjustment, names(completed))
  d <- copy(x)
  d[, (imputed) := NULL]
  d <- merge(d, completed[, c("id", imputed), with = FALSE], by = "id", all.x = TRUE, sort = FALSE)
  required <- c("id", "time_msm", "state_msm", domain_vars, primary_m4_adjustment)
  d <- d[complete.cases(d[, ..required])]
  eligible <- d[, .N, by = id][N >= 2L, id]
  d <- d[id %in% eligible]
  setorder(d, id, time_msm, state_msm)
  d[, race_bl := droplevels(race_bl)]
  d[, education_level_bl := droplevels(education_level_bl)]
  d[, wealth_quintile_bl := droplevels(wealth_quintile_bl)]
  d[]
}

make_constraints <- function(formula, d, edges) {
  mm <- model.matrix(formula, data = d)
  coefficient_names <- setdiff(colnames(mm), "(Intercept)")
  missing_domains <- setdiff(domain_vars, coefficient_names)
  if (length(missing_domains)) stop("Joint-domain model matrix is missing: ", paste(missing_domains, collapse = ", "))
  constraints <- list()
  for (name in coefficient_names) {
    constraints[[name]] <- if (name %in% domain_vars) seq_len(nrow(edges)) else rep(1L, nrow(edges))
  }
  list(constraints = constraints, coefficient_names = coefficient_names)
}

extract_domain_block <- function(fit, d, edges, formula) {
  coefficient_names <- setdiff(colnames(model.matrix(formula, data = d)), "(Intercept)")
  columns <- match(domain_vars, coefficient_names)
  if (anyNA(columns)) stop("Unable to locate all symptom-domain coefficients")
  k <- nrow(edges)
  indices <- unlist(lapply(columns, function(j) k + (j - 1L) * k + seq_len(k)), use.names = FALSE)
  estimates <- as.numeric(fit$estimates)
  covariance <- as.matrix(fit$covmat)
  if (max(indices) > length(estimates) || max(indices) > nrow(covariance)) stop("Domain block indices exceed msm parameter dimensions")
  labels <- unlist(lapply(domain_vars, function(name) paste0(name, "__", edges$transition)), use.names = FALSE)
  theta <- estimates[indices]
  names(theta) <- labels
  dimnames(covariance) <- list(names(estimates), names(estimates))
  block_cov <- covariance[indices, indices, drop = FALSE]
  dimnames(block_cov) <- list(labels, labels)
  list(theta = theta, covariance = block_cov, labels = labels, indices = indices, coefficient_names = coefficient_names)
}

model_diagnostics <- function(fit, block, d, config, imputation, elapsed, warnings, cache_status) {
  hessian <- if (is.matrix(fit$opt$hessian)) fit$opt$hessian else fit$hessian
  hessian_values <- if (is.matrix(hessian) && all(is.finite(hessian))) eigen((hessian + t(hessian)) / 2, symmetric = TRUE, only.values = TRUE)$values else NA_real_
  cov_values <- if (all(is.finite(block$covariance))) eigen((block$covariance + t(block$covariance)) / 2, symmetric = TRUE, only.values = TRUE)$values else NA_real_
  convergence_code <- suppressWarnings(as.integer(fit$opt$convergence)[1L])
  coefficient_finite <- all(is.finite(as.numeric(fit$estimates)))
  full_covariance_finite <- is.matrix(fit$covmat) && all(is.finite(fit$covmat))
  hessian_pd <- !all(is.na(hessian_values)) && min(hessian_values, na.rm = TRUE) > max(abs(hessian_values), na.rm = TRUE) * sqrt(.Machine$double.eps)
  domain_covariance_pd <- !all(is.na(cov_values)) && min(cov_values, na.rm = TRUE) > max(abs(cov_values), na.rm = TRUE) * sqrt(.Machine$double.eps)
  accepted <- identical(convergence_code, 0L) && isTRUE(fit$foundse) && coefficient_finite &&
    full_covariance_finite && hessian_pd && domain_covariance_pd
  rejection_reasons <- c(
    if (!identical(convergence_code, 0L)) paste0("convergence_code=", convergence_code),
    if (!isTRUE(fit$foundse)) "foundse_false",
    if (!coefficient_finite) "nonfinite_coefficients",
    if (!full_covariance_finite) "nonfinite_full_covariance",
    if (!hessian_pd) "free_hessian_not_positive_definite",
    if (!domain_covariance_pd) "domain_covariance_not_positive_definite"
  )
  data.table(
    analysis_id = config$analysis_id, cohort = config$cohort, imputation = imputation,
    n_records = nrow(d), n_participants = uniqueN(d$id),
    n_free_parameters = suppressWarnings(as.numeric(attr(logLik(fit), "df"))),
    log_likelihood = suppressWarnings(as.numeric(logLik(fit))),
    convergence_code = convergence_code,
    found_se = as.integer(isTRUE(fit$foundse)), coefficient_finite = as.integer(coefficient_finite),
    full_covariance_finite = as.integer(full_covariance_finite),
    domain_covariance_finite = as.integer(all(is.finite(block$covariance))),
    domain_covariance_min_eigen = if (all(is.na(cov_values))) NA_real_ else min(cov_values, na.rm = TRUE),
    domain_covariance_condition_number = if (all(is.na(cov_values)) || min(cov_values, na.rm = TRUE) <= 0) NA_real_ else max(cov_values, na.rm = TRUE) / min(cov_values, na.rm = TRUE),
    hessian_finite = as.integer(!all(is.na(hessian_values))), hessian_positive_definite = as.integer(hessian_pd),
    domain_covariance_positive_definite = as.integer(domain_covariance_pd),
    hessian_min_eigen = if (all(is.na(hessian_values))) NA_real_ else min(hessian_values, na.rm = TRUE),
    elapsed_seconds = elapsed, warning_count = length(warnings), warnings = paste(unique(warnings), collapse = " | "),
    cache_status = cache_status, accepted = as.integer(accepted),
    rejection_reasons = paste(rejection_reasons, collapse = " | ")
  )
}

make_fit_spec <- function(x, config, mi_object, imputation, q_initial) {
  edges <- allowed_edges(config$n_states)
  completed <- as.data.table(mice::complete(mi_object$mids, action = imputation))
  d <- prepare_fit_sample(x, completed)
  formula <- reformulate(c(domain_vars, primary_m4_adjustment))
  constraint_object <- make_constraints(formula, d, edges)
  data_hash <- table_fingerprint(d, c("id", "time_msm", "state_msm", domain_vars, primary_m4_adjustment))
  spec_hash <- digest::digest(list(
    script_version, config$analysis_id, imputation, mi_object$signature, q_initial,
    constraint_object$constraints, max_iterations, data_hash
  ), algo = "xxhash64")
  list(
    edges = edges, d = d, formula = formula, constraint_object = constraint_object,
    data_hash = data_hash, spec_hash = spec_hash
  )
}

make_warm_start <- function(config, imputation, edges, coefficient_names) {
  none <- function(reason) list(covinits = NULL, mode = paste0("zero_default:", reason), source_spec_hash = NA_character_, max_abs = 0)
  if (!use_warm_start) return(none("disabled"))
  if (imputation <= 1L) return(none("reference_imputation"))
  tag <- sprintf("%s__MI%02d__imp01", config$analysis_id, mi_m)
  source_model <- file.path(model_root, paste0(tag, ".rds"))
  source_component <- file.path(model_root, paste0(tag, "__component.rds"))
  if (!file.exists(source_model) || !file.exists(source_component)) return(none("reference_cache_missing"))
  envelope <- tryCatch(readRDS(source_model), error = function(e) NULL)
  component <- tryCatch(readRDS(source_component), error = function(e) NULL)
  if (is.null(envelope) || is.null(component) || !isTRUE(envelope$accepted) || !isTRUE(component$accepted)) return(none("reference_not_accepted"))
  if (!identical(envelope$script_version, script_version) || !identical(envelope$analysis_id, config$analysis_id) ||
      !identical(envelope$coefficient_names, coefficient_names) || !identical(envelope$spec_hash, component$spec_hash)) {
    return(none("reference_incompatible"))
  }
  k <- nrow(edges)
  estimates <- as.numeric(envelope$fit$estimates)
  expected <- k * (1L + length(coefficient_names))
  if (length(estimates) != expected || any(!is.finite(estimates))) return(none("reference_parameter_shape"))
  qcov <- estimates[k + seq_len(k * length(coefficient_names))]
  covinits <- split(qcov, rep(coefficient_names, each = k))
  covinits <- covinits[coefficient_names]
  list(
    covinits = covinits, mode = "accepted_imp01_mle_covinits",
    source_spec_hash = envelope$spec_hash, max_abs = max(abs(qcov))
  )
}

fit_one_imputation <- function(x, config, mi_object, imputation, q_initial) {
  spec <- make_fit_spec(x, config, mi_object, imputation, q_initial)
  edges <- spec$edges
  d <- spec$d
  formula <- spec$formula
  constraint_object <- spec$constraint_object
  data_hash <- spec$data_hash
  spec_hash <- spec$spec_hash
  warm_start <- make_warm_start(config, imputation, edges, constraint_object$coefficient_names)
  tag <- sprintf("%s__MI%02d__imp%02d", config$analysis_id, mi_m, imputation)
  model_file <- file.path(model_root, paste0(tag, ".rds"))
  component_file <- file.path(model_root, paste0(tag, "__component.rds"))
  if (!force_fit && file.exists(model_file) && file.exists(component_file)) {
    envelope <- tryCatch(readRDS(model_file), error = function(e) NULL)
    component <- tryCatch(readRDS(component_file), error = function(e) NULL)
    if (!is.null(envelope) && identical(envelope$spec_hash, spec_hash) && !is.null(component) && identical(component$spec_hash, spec_hash)) {
      message(if (isTRUE(component$accepted)) "Reusing accepted model cache: " else "Reusing rejected model cache: ", tag)
      component$diagnostics[, cache_status := if (isTRUE(component$accepted)) "reused_accepted" else "reused_rejected"]
      return(component)
    }
  }
  warnings <- character()
  started <- Sys.time()
  fit <- withCallingHandlers(
    msm::msm(
      formula = state_msm ~ time_msm, subject = id, data = as.data.frame(d),
      qmatrix = q_initial, covariates = formula, constraint = constraint_object$constraints,
      covinits = warm_start$covinits,
      deathexact = config$death_exact_state, exacttimes = FALSE, center = TRUE, hessian = TRUE,
      control = list(fnscale = 10000, maxit = max_iterations, trace = 0)
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  block <- extract_domain_block(fit, d, edges, formula)
  diagnostics <- model_diagnostics(fit, block, d, config, imputation, elapsed, warnings, "created")
  diagnostics[, `:=`(
    optimizer_initialization = warm_start$mode,
    warm_start_source_spec_hash = warm_start$source_spec_hash,
    warm_start_max_abs = warm_start$max_abs
  )]
  accepted <- diagnostics$accepted[1L] == 1L
  envelope <- list(
    spec_hash = spec_hash, script_version = script_version, analysis_id = config$analysis_id,
    imputation = imputation, data_fingerprint = data_hash, constraints = constraint_object$constraints,
    coefficient_names = constraint_object$coefficient_names, fit = fit, diagnostics = diagnostics,
    accepted = accepted, optimizer_initialization = warm_start
  )
  component <- list(
    spec_hash = spec_hash, analysis_id = config$analysis_id, cohort = config$cohort,
    imputation = imputation, edges = edges, theta = block$theta, covariance = block$covariance,
    labels = block$labels, diagnostics = diagnostics, data_fingerprint = data_hash,
    script_version = script_version, mi_m = mi_m, accepted = accepted
  )
  saveRDS(envelope, model_file, compress = "gzip")
  saveRDS(component, component_file, compress = "gzip")
  message(sprintf("Fitted %s in %.1f seconds; acceptance=%s", tag, elapsed, accepted))
  component
}

fit_analysis <- function(x, config, mi_object) {
  edges <- allowed_edges(config$n_states)
  q_initial <- make_q_initial(x, edges, config$n_states)
  q_long <- as.data.table(as.table(q_initial))
  setnames(q_long, c("from_state", "to_state", "initial_intensity"))
  q_long[, analysis_id := config$analysis_id]
  fwrite(q_long, file.path(qc_root, paste0(config$analysis_id, "__q_initial.csv")), na = "")
  components <- list()
  for (imputation in seq.int(imp_start, imp_end)) {
    components[[as.character(imputation)]] <- tryCatch(
      fit_one_imputation(x, config, mi_object, imputation, q_initial),
      error = function(e) {
        message("FAILED ", config$analysis_id, " imputation ", imputation, ": ", conditionMessage(e))
        data.table(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), analysis_id = config$analysis_id, imputation = imputation, error = conditionMessage(e)) |>
          fwrite(file.path(log_root, sprintf("%s__imp%02d__failure.csv", config$analysis_id, imputation)), na = "")
        NULL
      }
    )
  }
  completed_components <- Filter(function(z) is.list(z) && !is.null(z$diagnostics), components)
  if (length(completed_components)) {
    shard_tag <- sprintf("imp%02d_%02d", imp_start, imp_end)
    partial_diagnostics <- rbindlist(lapply(completed_components, function(z) z$diagnostics), use.names = TRUE, fill = TRUE)
    fwrite(partial_diagnostics, file.path(qc_root, paste0(config$analysis_id, "__", shard_tag, "__fit_diagnostics_partial.csv")), na = "")
    inventory <- rbindlist(lapply(completed_components, function(z) {
      tag <- sprintf("%s__MI%02d__imp%02d", config$analysis_id, mi_m, z$imputation)
      model_file <- file.path(model_root, paste0(tag, ".rds"))
      component_file <- file.path(model_root, paste0(tag, "__component.rds"))
      data.table(
        analysis_id = config$analysis_id, imputation = z$imputation,
        spec_hash = z$spec_hash, accepted = as.integer(isTRUE(z$accepted)),
        model_file = model_file, component_file = component_file,
        model_sha256 = if (file.exists(model_file)) file_sha256(model_file) else NA_character_,
        component_sha256 = if (file.exists(component_file)) file_sha256(component_file) else NA_character_
      )
    }))
    fwrite(inventory, file.path(qc_root, paste0(config$analysis_id, "__", shard_tag, "__model_inventory_partial.csv")), na = "")
  }
  invisible(components)
}

component_paths <- function(analysis_id) {
  file.path(model_root, sprintf("%s__MI%02d__imp%02d__component.rds", analysis_id, mi_m, seq_len(mi_m)))
}

read_components <- function(config, x, mi_object) {
  paths <- component_paths(config$analysis_id)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stop("Missing ", length(missing), " of ", mi_m, " fitted component caches for ", config$analysis_id)
  out <- lapply(paths, readRDS)
  edges <- allowed_edges(config$n_states)
  q_initial <- make_q_initial(x, edges, config$n_states)
  for (i in seq_along(out)) {
    expected <- make_fit_spec(x, config, mi_object, i, q_initial)$spec_hash
    z <- out[[i]]
    if (!identical(z$script_version, script_version) || !identical(z$mi_m, mi_m) ||
        !identical(z$analysis_id, config$analysis_id) || !identical(z$imputation, i) ||
        !identical(z$spec_hash, expected)) {
      stop("Stale or incompatible component cache for ", config$analysis_id, " imputation ", i)
    }
    if (!isTRUE(z$accepted) || z$diagnostics$accepted[1L] != 1L) {
      stop("Component failed the frozen post-hoc-extension acceptance gate for ", config$analysis_id, " imputation ", i,
           ": ", z$diagnostics$rejection_reasons[1L])
    }
  }
  labels <- lapply(out, function(z) z$labels)
  if (!all(vapply(labels, identical, logical(1), labels[[1L]]))) stop("Coefficient labels/order differ across imputations for ", config$analysis_id)
  out
}

require_pd_inverse <- function(x, label = "covariance", tolerance = NULL) {
  x <- (x + t(x)) / 2
  if (!all(is.finite(x))) stop(label, " contains non-finite values")
  e <- eigen(x, symmetric = TRUE)
  if (is.null(tolerance)) tolerance <- max(abs(e$values)) * max(dim(x)) * sqrt(.Machine$double.eps)
  if (any(e$values <= tolerance)) {
    stop(label, " is not full-rank positive definite; min eigenvalue=", signif(min(e$values), 6),
         ", tolerance=", signif(tolerance, 6))
  }
  list(inverse = e$vectors %*% diag(1 / e$values, nrow = length(e$values)) %*% t(e$vectors), rank = length(e$values), eigenvalues = e$values)
}

rubin_pool_block <- function(components) {
  theta <- do.call(cbind, lapply(components, function(z) z$theta))
  covariance <- simplify2array(lapply(components, function(z) z$covariance))
  qbar <- rowMeans(theta)
  ubar <- apply(covariance, c(1L, 2L), mean)
  between <- stats::cov(t(theta))
  total <- ubar + (1 + 1 / ncol(theta)) * between
  list(theta_m = theta, covariance_m = covariance, qbar = qbar, ubar = ubar, between = between, total = total, m = ncol(theta))
}

d1_test <- function(pool, contrast) {
  contrast <- as.matrix(contrast)
  qhat <- contrast %*% pool$theta_m
  uhat <- array(NA_real_, dim = c(nrow(contrast), nrow(contrast), pool$m))
  for (i in seq_len(pool$m)) uhat[, , i] <- contrast %*% pool$covariance_m[, , i] %*% t(contrast)
  qbar <- rowMeans(qhat)
  ubar <- apply(uhat, c(1L, 2L), mean)
  b <- if (nrow(contrast) == 1L) matrix(stats::var(as.numeric(qhat)), 1L, 1L) else stats::cov(t(qhat))
  k <- nrow(contrast)
  inv_u <- require_pd_inverse(ubar, "D1 within-imputation covariance")
  total_c <- contrast %*% pool$total %*% t(contrast)
  inv_total <- require_pd_inverse(total_c, "Rubin total contrast covariance")
  if (k == 1L) {
    u_scalar <- ubar[1L, 1L]
    b_scalar <- b[1L, 1L]
    r <- (1 + 1 / pool$m) * b_scalar / u_scalar
    r_for_df <- max(r, .Machine$double.eps)
    v <- (pool$m - 1) * (1 + 1 / r_for_df)^2
    estimate <- qbar[1L]
    se <- sqrt(total_c[1L, 1L])
    t_stat <- estimate / se
    return(list(
      estimate = qbar, covariance = total_c, F = t_stat^2, numerator_df = 1L,
      denominator_df = v, p_value = 2 * stats::pt(-abs(t_stat), df = v),
      relative_increase_variance = r, chi_square = t_stat^2, chi_square_df = 1L,
      chi_square_p = stats::pchisq(t_stat^2, 1L, lower.tail = FALSE),
      method = "Scalar Rubin t with old degrees of freedom"
    ))
  }
  r <- (1 + 1 / pool$m) * sum(diag(b %*% inv_u$inverse)) / k
  r <- max(r, .Machine$double.eps)
  tdf <- k * (pool$m - 1)
  v <- if (tdf > 4) {
    4 + (tdf - 4) * (1 + (1 - 2 / tdf) / r)^2
  } else {
    tdf * (1 + 1 / k) * (1 + 1 / r)^2 / 2
  }
  ttilde <- (1 + r) * ubar
  inv_ttilde <- require_pd_inverse(ttilde, "D1 scaled within-imputation covariance")
  f_stat <- as.numeric(t(qbar) %*% inv_ttilde$inverse %*% qbar / k)
  chi <- as.numeric(t(qbar) %*% inv_total$inverse %*% qbar)
  list(
    estimate = qbar, covariance = total_c, F = f_stat, numerator_df = k, denominator_df = v,
    p_value = stats::pf(f_stat, k, v, lower.tail = FALSE), relative_increase_variance = r,
    chi_square = chi, chi_square_df = inv_total$rank,
    chi_square_p = stats::pchisq(chi, inv_total$rank, lower.tail = FALSE),
    method = "Multivariate D1 F; Rubin total-covariance chi-square sensitivity"
  )
}

contrast_row <- function(labels, positive, negative = NULL) {
  out <- numeric(length(labels)); names(out) <- labels
  out[positive] <- 1
  if (!is.null(negative)) out[negative] <- -1
  out
}

build_tests <- function(config, pool) {
  edges <- allowed_edges(config$n_states)
  labels <- rownames(pool$total)
  if (is.null(labels)) labels <- names(pool$qbar)
  nonfatal <- edges[death_transition == 0L, transition]
  death <- edges[death_transition == 1L, transition]
  families <- list(common_nonfatal = nonfatal)
  if (length(death)) families$hrs_death_exploratory <- death
  joint_rows <- list(); pair_rows <- list(); estimate_rows <- list()
  for (domain in domain_vars) {
    for (transition in edges$transition) {
      label <- paste0(domain, "__", transition)
      cmat <- matrix(contrast_row(labels, label), nrow = 1)
      test <- d1_test(pool, cmat)
      se <- sqrt(test$covariance[1, 1])
      df <- test$denominator_df
      estimate_rows[[length(estimate_rows) + 1L]] <- data.table(
        analysis_id = config$analysis_id, cohort = config$cohort,
        transition = transition, death_transition = as.integer(transition %in% death),
        domain = domain, domain_label = unname(domain_labels[domain]),
        beta = pool$qbar[label], standard_error = se,
        transition_intensity_ratio = exp(pool$qbar[label]),
        conf_low = exp(pool$qbar[label] - qt(.975, df) * se),
        conf_high = exp(pool$qbar[label] + qt(.975, df) * se),
        denominator_df = df, p_value = test$p_value, m = pool$m,
        inference_method = test$method,
        covariance_gate = "full-rank positive definite; no pseudoinverse"
      )
    }
  }
  pair_def <- list(
    onset_vs_maintenance = c(domain_vars[1], domain_vars[2]),
    onset_vs_nonrestorative = c(domain_vars[1], domain_vars[3]),
    maintenance_vs_nonrestorative = c(domain_vars[2], domain_vars[3])
  )
  for (family in names(families)) {
    transitions <- families[[family]]
    global_contrasts <- list()
    for (transition in transitions) {
      global_contrasts[[length(global_contrasts) + 1L]] <- contrast_row(labels, paste0(domain_vars[1], "__", transition), paste0(domain_vars[2], "__", transition))
      global_contrasts[[length(global_contrasts) + 1L]] <- contrast_row(labels, paste0(domain_vars[1], "__", transition), paste0(domain_vars[3], "__", transition))
    }
    global_c <- do.call(rbind, global_contrasts)
    global_test <- d1_test(pool, global_c)
    global_pass <- is.finite(global_test$p_value) && global_test$p_value < .05
    joint_rows[[length(joint_rows) + 1L]] <- data.table(
      analysis_id = config$analysis_id, cohort = config$cohort, family = family,
      test_scope = "global_domain_equality", transition = paste(transitions, collapse = " | "), domain = NA_character_,
      F = global_test$F, numerator_df = global_test$numerator_df, denominator_df = global_test$denominator_df,
      p_value = global_test$p_value, chi_square = global_test$chi_square,
      chi_square_df = global_test$chi_square_df, chi_square_p = global_test$chi_square_p,
      relative_increase_variance = global_test$relative_increase_variance,
      global_gate_pass = as.integer(global_pass),
      inference_method = global_test$method,
      covariance_gate = "full-rank positive definite; no pseudoinverse",
      status = if (family == "common_nonfatal") "post-hoc secondary; contrasts frozen before new fit" else "exploratory death family"
    )
    for (transition in transitions) {
      edge_c <- rbind(
        contrast_row(labels, paste0(domain_vars[1], "__", transition), paste0(domain_vars[2], "__", transition)),
        contrast_row(labels, paste0(domain_vars[1], "__", transition), paste0(domain_vars[3], "__", transition))
      )
      edge_test <- d1_test(pool, edge_c)
      joint_rows[[length(joint_rows) + 1L]] <- data.table(
        analysis_id = config$analysis_id, cohort = config$cohort, family = family,
        test_scope = "per_transition_domain_equality", transition = transition, domain = NA_character_,
        F = edge_test$F, numerator_df = edge_test$numerator_df, denominator_df = edge_test$denominator_df,
        p_value = edge_test$p_value, chi_square = edge_test$chi_square,
        chi_square_df = edge_test$chi_square_df, chi_square_p = edge_test$chi_square_p,
        relative_increase_variance = edge_test$relative_increase_variance,
        global_gate_pass = as.integer(global_pass),
        inference_method = edge_test$method,
        covariance_gate = "full-rank positive definite; no pseudoinverse",
        status = if (global_pass) "gate open" else "descriptive only; global gate closed"
      )
      for (pair_name in names(pair_def)) {
        pair <- pair_def[[pair_name]]
        cmat <- matrix(contrast_row(labels, paste0(pair[1], "__", transition), paste0(pair[2], "__", transition)), nrow = 1)
        pair_test <- d1_test(pool, cmat)
        estimate <- pair_test$estimate[1]
        se <- sqrt(pair_test$covariance[1, 1])
        pair_rows[[length(pair_rows) + 1L]] <- data.table(
          analysis_id = config$analysis_id, cohort = config$cohort, family = family,
          transition = transition, contrast = pair_name,
          domain_a = pair[1], domain_b = pair[2],
          beta_difference = estimate, standard_error = se,
          ratio_of_transition_intensity_ratios = exp(estimate),
          conf_low = exp(estimate - qt(.975, pair_test$denominator_df) * se),
          conf_high = exp(estimate + qt(.975, pair_test$denominator_df) * se),
          F = pair_test$F, numerator_df = 1L, denominator_df = pair_test$denominator_df,
          p_value = pair_test$p_value, chi_square = pair_test$chi_square,
          chi_square_p = pair_test$chi_square_p, global_gate_pass = as.integer(global_pass),
          inference_method = pair_test$method,
          covariance_gate = "full-rank positive definite; no pseudoinverse"
        )
      }
    }
  }
  for (domain in domain_vars) {
    cmat <- do.call(rbind, lapply(nonfatal, function(transition) contrast_row(labels, paste0(domain, "__", transition))))
    test <- d1_test(pool, cmat)
    joint_rows[[length(joint_rows) + 1L]] <- data.table(
      analysis_id = config$analysis_id, cohort = config$cohort, family = "common_nonfatal",
      test_scope = "per_domain_joint_association", transition = paste(nonfatal, collapse = " | "), domain = domain,
      F = test$F, numerator_df = test$numerator_df, denominator_df = test$denominator_df,
      p_value = test$p_value, chi_square = test$chi_square, chi_square_df = test$chi_square_df,
      chi_square_p = test$chi_square_p, relative_increase_variance = test$relative_increase_variance,
      global_gate_pass = NA_integer_,
      inference_method = test$method,
      covariance_gate = "full-rank positive definite; no pseudoinverse",
      status = "association test; not a between-domain difference test"
    )
  }
  estimates <- rbindlist(estimate_rows, use.names = TRUE, fill = TRUE)
  joint <- rbindlist(joint_rows, use.names = TRUE, fill = TRUE)
  pairwise <- rbindlist(pair_rows, use.names = TRUE, fill = TRUE)
  pairwise[, p_bh := p.adjust(p_value, method = "BH"), by = .(analysis_id, family)]
  pairwise[, inferential_status := fifelse(
    global_gate_pass == 1L & p_bh < .05, "gate open and BH < 0.05",
    fifelse(global_gate_pass == 1L, "gate open; BH >= 0.05", "descriptive only; global gate closed")
  )]
  list(estimates = estimates, joint = joint, pairwise = pairwise)
}

pool_analysis <- function(config, x, mi_object) {
  components <- read_components(config, x, mi_object)
  pool <- rubin_pool_block(components)
  tests <- build_tests(config, pool)
  diagnostics <- rbindlist(lapply(components, function(z) z$diagnostics), use.names = TRUE, fill = TRUE)
  coefficients <- rbindlist(lapply(components, function(z) data.table(
    analysis_id = config$analysis_id, cohort = config$cohort, imputation = z$imputation,
    parameter = names(z$theta), beta = as.numeric(z$theta)
  )))
  covariance <- rbindlist(lapply(seq_len(pool$m), function(i) {
    value <- as.data.table(as.table(pool$covariance_m[, , i]))
    setnames(value, c("parameter_row", "parameter_col", "within_covariance"))
    value[, `:=`(analysis_id = config$analysis_id, imputation = i)]
    value
  }))
  pooled_cov <- as.data.table(as.table(pool$total)); setnames(pooled_cov, c("parameter_row", "parameter_col", "total_covariance"))
  ubar <- as.data.table(as.table(pool$ubar)); setnames(ubar, c("parameter_row", "parameter_col", "within_covariance_mean"))
  between <- as.data.table(as.table(pool$between)); setnames(between, c("parameter_row", "parameter_col", "between_covariance"))
  pooled_cov <- Reduce(function(a, b) merge(a, b, by = c("parameter_row", "parameter_col"), all = TRUE), list(pooled_cov, ubar, between))
  pooled_cov[, analysis_id := config$analysis_id]
  inventory_expected <- rbindlist(lapply(components, function(z) data.table(
    analysis_id = z$analysis_id, imputation = z$imputation, spec_hash = z$spec_hash,
    accepted = as.integer(isTRUE(z$accepted))
  )))
  list(
    tests = tests, diagnostics = diagnostics, coefficients = coefficients,
    covariance = covariance, pooled_covariance = pooled_cov,
    inventory_expected = inventory_expected
  )
}

render_raster_staged <- function(plot, target, format = c("png", "tiff"), width, height, res = 600) {
  format <- match.arg(format)
  staged <- tempfile(fileext = paste0(".", format))
  on.exit(unlink(staged), add = TRUE)
  device_open <- FALSE
  tryCatch({
    if (format == "png") {
      ragg::agg_png(staged, width = width, height = height, units = "in", res = res)
    } else {
      ragg::agg_tiff(staged, width = width, height = height, units = "in", res = res, compression = "lzw")
    }
    device_open <- TRUE
    print(plot)
  }, finally = {
    if (device_open) grDevices::dev.off()
  })
  if (!file.exists(staged) || !is.finite(file.info(staged)$size) || file.info(staged)$size <= 0) {
    stop("R raster staging failed for ", target)
  }
  if (!file.copy(staged, target, overwrite = TRUE)) stop("Unable to copy staged R raster to ", target)
  invisible(target)
}

export_figure_set <- function(plot, base, width, height, res = 600) {
  svglite::svglite(paste0(base, ".svg"), width = width, height = height)
  print(plot)
  grDevices::dev.off()
  grDevices::pdf(
    paste0(base, ".pdf"), width = width, height = height,
    family = "Helvetica", useDingbats = FALSE
  )
  print(plot)
  grDevices::dev.off()
  render_raster_staged(plot, paste0(base, ".tiff"), "tiff", width, height, res)
  render_raster_staged(plot, paste0(base, ".png"), "png", width, height, res)
  invisible(base)
}

# patchwork stores panel tracks as "null" units.  With ggplot2 4.0 these units
# convert to zero outside their parent layout, so resolve the remaining device
# space explicitly before passing the geometry to the strict Nature audit.
write_patchwork_panel_layout_resolved <- function(
  plot, manifest_path, width_in, height_in, panel_ids = NULL,
  row_groups = NULL, column_groups = NULL, exemptions = list()
) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required for alignment QA")
  grob <- patchwork::patchworkGrob(plot)
  panel_rows <- grob$layout[grepl("^panel(?:-[0-9]+)?$", grob$layout$name, perl = TRUE), , drop = FALSE]
  panel_rows <- panel_rows[order(panel_rows$t, panel_rows$l), , drop = FALSE]
  if (nrow(panel_rows) < 2L) stop("At least two flat patchwork panel cells are required for alignment QA")
  if (is.null(panel_ids)) panel_ids <- letters[seq_len(nrow(panel_rows))]
  if (length(panel_ids) != nrow(panel_rows) || anyDuplicated(panel_ids)) stop("Invalid alignment panel IDs")

  probe_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(probe_path, width = width_in, height = height_in, useDingbats = FALSE)
  device_open <- TRUE
  on.exit({
    if (device_open) grDevices::dev.off()
    unlink(probe_path)
  }, add = TRUE)
  grid::grid.newpage(); grid::grid.draw(grob); grid::grid.force()
  resolve_tracks <- function(units, total_pt, axis = c("width", "height")) {
    axis <- match.arg(axis)
    types <- grid::unitType(units)
    values <- if (axis == "width") {
      grid::convertWidth(units, "pt", valueOnly = TRUE)
    } else {
      grid::convertHeight(units, "pt", valueOnly = TRUE)
    }
    null <- types == "null"
    if (any(null)) {
      weights <- as.numeric(units[null])
      remaining <- total_pt - sum(values[!null])
      if (!is.finite(remaining) || remaining <= 0 || any(!is.finite(weights)) || sum(weights) <= 0) {
        stop("Unable to resolve patchwork null tracks for alignment QA")
      }
      values[null] <- remaining * weights / sum(weights)
    }
    if (any(!is.finite(values)) || any(values < 0)) stop("Invalid resolved patchwork track geometry")
    values
  }
  widths_pt <- resolve_tracks(grob$widths, width_in * 72, "width")
  heights_pt <- resolve_tracks(grob$heights, height_in * 72, "height")
  x_edges <- c(0, cumsum(widths_pt))
  top_edges <- c(0, cumsum(heights_pt))
  total_height_pt <- sum(heights_pt)
  panels <- lapply(seq_len(nrow(panel_rows)), function(index) {
    row <- panel_rows[index, ]
    list(
      id = as.character(panel_ids[index]),
      bbox_pt = unname(c(
        x_edges[row$l], total_height_pt - top_edges[row$b + 1L],
        x_edges[row$r + 1L], total_height_pt - top_edges[row$t]
      )),
      grid_id = "patchwork-grid-1",
      row_start = as.integer(row$t - 1L), row_stop = as.integer(row$b),
      col_start = as.integer(row$l - 1L), col_stop = as.integer(row$r)
    )
  })
  manifest <- list(
    schema_version = 1L, backend = "r-patchwork-resolved-null-tracks",
    figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
    panels = panels, exemptions = exemptions
  )
  if (!is.null(row_groups)) manifest$row_groups <- lapply(row_groups, as.character)
  if (!is.null(column_groups)) manifest$column_groups <- lapply(column_groups, as.character)
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  grDevices::dev.off(); device_open <- FALSE; unlink(probe_path)
  invisible(manifest)
}

write_model_inventory <- function(expected) {
  partial_files <- list.files(qc_root, pattern = "model_inventory_partial[.]csv$", full.names = TRUE)
  if (!length(partial_files)) stop("No model inventory shards were found")
  partial <- rbindlist(lapply(partial_files, fread), use.names = TRUE, fill = TRUE)
  partial <- partial[accepted == 1L & nzchar(spec_hash)]
  matched <- merge(
    expected, partial,
    by = c("analysis_id", "imputation", "spec_hash", "accepted"), all.x = TRUE,
    allow.cartesian = TRUE
  )
  if (anyNA(matched$model_file) || anyNA(matched$component_file)) {
    stop("At least one current accepted component lacks a matching inventory shard")
  }
  duplicates <- matched[, .N, by = .(analysis_id, imputation)][N != 1L]
  if (nrow(duplicates)) stop("Current model inventory is duplicated or incomplete")
  if (nrow(matched) != nrow(expected) || uniqueN(matched, by = c("analysis_id", "imputation")) != nrow(expected)) {
    stop("Current model inventory does not contain exactly the expected accepted components")
  }
  missing_paths <- c(matched$model_file[!file.exists(matched$model_file)], matched$component_file[!file.exists(matched$component_file)])
  if (length(missing_paths)) stop("Inventoried model artifacts are missing: ", paste(missing_paths, collapse = "; "))
  matched[, component_sha256_current := vapply(component_file, file_sha256, character(1))]
  matched[, component_hash_matches_fit_inventory := as.integer(component_sha256_current == component_sha256)]
  if (any(matched$component_hash_matches_fit_inventory != 1L)) stop("A component hash changed after model fitting")
  setorder(matched, analysis_id, imputation)
  fwrite(matched, file.path(qc_root, "symptom_domain_model_inventory.csv"), na = "")
  matched[]
}

write_reproducibility_bundle <- function(model_inventory) {
  package_names <- c("R", "data.table", "msm", "mice", "digest", "ggplot2", "patchwork", "svglite", "ragg", "coda")
  package_versions <- data.table(
    package = package_names,
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      vapply(package_names[-1L], function(name) {
        if (requireNamespace(name, quietly = TRUE)) as.character(utils::packageVersion(name)) else NA_character_
      }, character(1))
    )
  )
  fwrite(package_versions, file.path(qc_root, "software_versions.csv"), na = "")
  writeLines(capture.output(utils::sessionInfo()), file.path(qc_root, "session_info.txt"), useBytes = TRUE)

  specification <- rbindlist(lapply(analysis_configs, function(config) {
    edges <- allowed_edges(config$n_states)
    data.table(
      script_version = script_version,
      analysis_status = if (config$analysis_id == "hrs_5state_joint_domains" || config$analysis_id == "elsa_4state_w4_w10_joint_domains") {
        "post-hoc secondary extension; contrasts frozen before new joint-model fitting"
      } else "structural-comparability sensitivity",
      analysis_id = config$analysis_id, cohort = config$cohort,
      input_file = file.path(input_root, config$input),
      input_sha256 = file_sha256(file.path(input_root, config$input)),
      n_states = config$n_states, allowed_edges = paste(edges$transition, collapse = " | "),
      exact_death_state = if (is.null(config$death_exact_state)) NA_integer_ else config$death_exact_state,
      mi_id = config$mi_id, imputations = mi_m, fcs_iterations = mi_maxit,
      mi_seed = mi_seed + ifelse(config$cohort == "HRS", 1L, 2L),
      domain_variables = paste(domain_vars, collapse = " | "),
      domain_constraint = "transition-specific",
      nuisance_variables = paste(primary_m4_adjustment, collapse = " | "),
      nuisance_constraint = "common across allowed transitions",
      optimizer_max_iterations = max_iterations,
      accepted_components = model_inventory[analysis_id == config$analysis_id, .N]
    )
  }), use.names = TRUE, fill = TRUE)
  fwrite(specification, file.path(qc_root, "model_specification_manifest.csv"), na = "")

  manifest_target <- file.path(qc_root, "symptom_domain_output_file_manifest.csv")
  artifact_roots <- c(result_root, figure_root, report_root, qc_root)
  files <- unique(unlist(lapply(artifact_roots, function(path) {
    list.files(path, recursive = TRUE, full.names = TRUE)
  }), use.names = FALSE))
  files <- files[file.info(files)$isdir %in% FALSE]
  normalised <- gsub("\\\\", "/", files)
  excluded <- grepl("/(model_cache|mi_cache)/", normalised) |
    normalizePath(files, winslash = "/", mustWork = FALSE) == normalizePath(manifest_target, winslash = "/", mustWork = FALSE)
  files <- files[!excluded]
  file_manifest <- data.table(
    path = substring(gsub("\\\\", "/", files), nchar(gsub("\\\\", "/", output_root)) + 2L),
    bytes = as.numeric(file.info(files)$size),
    sha256 = vapply(files, file_sha256, character(1))
  )
  setorder(file_manifest, path)
  fwrite(file_manifest, manifest_target, na = "")
  invisible(list(packages = package_versions, specification = specification, files = file_manifest))
}

classical_chain_psrf <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (n < 2L || ncol(x) < 2L || any(!is.finite(x))) return(NA_real_)
  within <- mean(apply(x, 2L, stats::var))
  between <- n * stats::var(colMeans(x))
  if (!is.finite(within) || within <= 0 || !is.finite(between)) return(NA_real_)
  sqrt((((n - 1) / n) * within + between / n) / within)
}

coda_chain_psrf <- function(x) {
  if (!requireNamespace("coda", quietly = TRUE)) return(c(point = NA_real_, upper = NA_real_))
  result <- tryCatch(
    coda::gelman.diag(
      coda::mcmc.list(lapply(seq_len(ncol(x)), function(j) coda::mcmc(x[, j]))),
      autoburnin = FALSE, multivariate = FALSE
    )$psrf[1L, ],
    error = function(e) c(NA_real_, NA_real_)
  )
  setNames(as.numeric(result), c("point", "upper"))
}

write_mi_diagnostics <- function(mi_objects) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("svglite", quietly = TRUE) ||
      !requireNamespace("ragg", quietly = TRUE)) {
    stop("ggplot2, svglite and ragg are required for R-only MI diagnostic figures")
  }
  cohort_by_mi <- c(hrs_joint_domains = "HRS", elsa_joint_domains = "ELSA")
  trace_rows <- list(); summary_rows <- list(); event_rows <- list(); event_count_rows <- list()
  for (mi_id in names(mi_objects)) {
    object <- mi_objects[[mi_id]]
    imp <- object$mids
    manifest <- readRDS(object$paths$manifest)
    methods <- manifest$methods
    variables <- names(methods)[methods != ""]
    cohort <- unname(cohort_by_mi[mi_id])
    for (variable in variables) {
      mean_matrix <- imp$chainMean[variable, , , drop = TRUE]
      sd_matrix <- sqrt(imp$chainVar[variable, , , drop = TRUE])
      if (!is.matrix(mean_matrix)) mean_matrix <- matrix(mean_matrix, nrow = imp$iteration, ncol = imp$m)
      if (!is.matrix(sd_matrix)) sd_matrix <- matrix(sd_matrix, nrow = imp$iteration, ncol = imp$m)
      iterations <- seq_len(nrow(mean_matrix))
      psrf_by_iteration <- vapply(iterations, function(last) {
        if (last < 2L) NA_real_ else classical_chain_psrf(mean_matrix[seq_len(last), , drop = FALSE])
      }, numeric(1))
      adjacent_ac <- c(NA_real_, vapply(2:nrow(mean_matrix), function(iteration) {
        value <- suppressWarnings(stats::cor(
          mean_matrix[iteration - 1L, ], mean_matrix[iteration, ], use = "pairwise.complete.obs"
        ))
        if (is.finite(value)) value else 0
      }, numeric(1)))
      cumulative_ac <- c(NA_real_, cumsum(adjacent_ac[-1L]) / seq_along(adjacent_ac[-1L]))
      long <- as.data.table(as.table(mean_matrix))
      setnames(long, c("iteration", "imputation", "chain_mean"))
      long[, `:=`(
        iteration = as.integer(iteration),
        imputation = as.integer(sub("[^0-9]+", "", as.character(imputation))),
        cohort = cohort, mi_id = mi_id, variable = variable,
        method = unname(methods[variable])
      )]
      long[, chain_sd := as.numeric(sd_matrix[cbind(iteration, imputation)])]
      long[, `:=`(
        classical_psrf = psrf_by_iteration[iteration],
        cumulative_adjacent_chain_correlation = cumulative_ac[iteration]
      )]
      trace_rows[[length(trace_rows) + 1L]] <- long
      coda_value <- coda_chain_psrf(mean_matrix)
      summary_rows[[length(summary_rows) + 1L]] <- data.table(
        cohort = cohort, mi_id = mi_id, variable = variable,
        method = unname(methods[variable]), n_imputations = ncol(mean_matrix),
        iterations = nrow(mean_matrix), finite_chain_mean_values = sum(is.finite(mean_matrix)),
        finite_chain_sd_values = sum(is.finite(sd_matrix)),
        classical_psrf_final = tail(psrf_by_iteration, 1L),
        coda_psrf_point_full_trace = coda_value["point"],
        coda_psrf_upper_full_trace = coda_value["upper"],
        cumulative_adjacent_chain_correlation_final = tail(cumulative_ac, 1L)
      )
    }
    logged <- imp$loggedEvents
    event_count_rows[[length(event_count_rows) + 1L]] <- data.table(
      cohort = cohort, mi_id = mi_id,
      total_logged_notices = if (is.null(logged)) 0L else nrow(logged)
    )
    if (!is.null(logged) && nrow(logged)) {
      logged <- as.data.table(logged)
      event_rows[[length(event_rows) + 1L]] <- logged[, .(
        notices = .N
      ), by = .(iteration = it, target = dep, method = meth, removed_predictors = out)][
        , `:=`(cohort = cohort, mi_id = mi_id)
      ]
    }
  }
  trace <- rbindlist(trace_rows, fill = TRUE)
  summary <- rbindlist(summary_rows, fill = TRUE)
  events <- rbindlist(event_rows, fill = TRUE)
  event_counts <- rbindlist(event_count_rows, fill = TRUE)
  fwrite(trace, file.path(qc_root, "mi_chain_mean_trace.csv"), na = "")
  fwrite(summary, file.path(qc_root, "mi_chain_diagnostics_summary.csv"), na = "")
  fwrite(events, file.path(qc_root, "mi_logged_event_summary.csv"), na = "")
  fwrite(event_counts, file.path(qc_root, "mi_logged_event_counts.csv"), na = "")
  fwrite(trace, file.path(source_root, "figure_mi_chain_mean_trace_source.csv"), na = "")

  plot_data <- copy(trace)
  plot_data[, variable_label := gsub("_bl$", "", gsub("_", " ", variable))]
  plot_data[, cohort := factor(cohort, levels = c("HRS", "ELSA"))]
  plot_data[, panel_label := paste(cohort, variable_label, sep = " | ")]
  plot_data[, panel_label := factor(panel_label, levels = unique(panel_label))]
  figure <- ggplot2::ggplot(plot_data, ggplot2::aes(x = iteration, y = chain_mean, group = imputation)) +
    ggplot2::geom_line(linewidth = 0.22, alpha = 0.3, colour = "#2F6B9A") +
    ggplot2::stat_summary(ggplot2::aes(group = 1L), fun = mean, geom = "line", linewidth = 0.65, colour = "#B64A3A") +
    ggplot2::facet_wrap(~panel_label, scales = "free_y", ncol = 2) +
    ggplot2::scale_x_continuous(breaks = c(1, 5, 10)) +
    ggplot2::labs(
      x = "FCS iteration", y = "Stored imputation-chain mean",
      title = "Multiple-imputation chain-mean diagnostics",
      subtitle = "Thin lines: each imputation chain; red line: across-chain mean"
    ) +
    ggplot2::theme_classic(base_size = 7.2, base_family = "Helvetica") +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold", size = 6.5),
      axis.line = ggplot2::element_line(linewidth = 0.35), axis.ticks = ggplot2::element_line(linewidth = 0.35),
      plot.title = ggplot2::element_text(face = "bold", size = 8.2),
      plot.subtitle = ggplot2::element_text(size = 6.5), panel.spacing = grid::unit(2.5, "mm")
    )
  base <- file.path(figure_root, "fig_mi_chain_mean_diagnostics")
  width_in <- 183 / 25.4; height_in <- 250 / 25.4
  export_figure_set(figure, base, width_in, height_in, 600)
  writeLines(c(
    "Supplementary QC figure | Multiple-imputation chain-mean diagnostics.",
    "Each thin line is the stored mean imputed value in one FCS chain at each iteration; the red line is the across-chain mean. Ordered and binary factors use their internal numeric codes, so levels should not be interpreted as prevalence. These are imputation-process diagnostics, not MCMC diagnostics for the multistate-model parameters. Classical PSRF and coda upper confidence limits are supplied in the accompanying source table because the installed mice build requires the unavailable optional rstan package for mice::convergence().",
    "",
    "补充质控图 | 多重插补链均值诊断。",
    "细线表示各FCS插补链每轮保存的插补均值，红线表示20条链的均值。有序变量和二元因子使用内部数值编码，不能把纵轴直接解释为患病率。本图诊断插补过程，不是多状态模型参数的MCMC诊断。由于本机mice版本的convergence()依赖未安装的可选rstan包，经典PSRF及coda上限另见配套源表。"
  ), paste0(base, "_caption.md"), useBytes = TRUE)
  list(trace = trace, summary = summary, events = events, event_counts = event_counts)
}

format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

analysis_display <- c(
  hrs_5state_joint_domains = "HRS five-state model",
  elsa_4state_w4_w10_joint_domains = "ELSA W4-W10 four-state model",
  hrs_4state_death_censored_joint_domains = "HRS four-state death-censored sensitivity"
)
transition_display <- c(
  "S0->S1" = "No disease -> hypertension only",
  "S0->S2" = "No disease -> one CVD without hypertension",
  "S0->S4" = "No disease -> death",
  "S1->S3" = "Hypertension only -> multimorbidity",
  "S1->S4" = "Hypertension only -> death",
  "S2->S3" = "One CVD without hypertension -> multimorbidity",
  "S2->S4" = "One CVD without hypertension -> death",
  "S3->S4" = "Multimorbidity -> death"
)
transition_display_zh <- c(
  "S0->S1" = "无目标疾病→单纯高血压",
  "S0->S2" = "无目标疾病→无高血压的单种心血管疾病",
  "S0->S4" = "无目标疾病→死亡",
  "S1->S3" = "单纯高血压→心血管多病共存",
  "S1->S4" = "单纯高血压→死亡",
  "S2->S3" = "无高血压的单种心血管疾病→心血管多病共存",
  "S2->S4" = "无高血压的单种心血管疾病→死亡",
  "S3->S4" = "心血管多病共存→死亡"
)
pairwise_display_zh <- c(
  onset_vs_maintenance = "入睡困难/睡眠维持困难",
  onset_vs_nonrestorative = "入睡困难/非恢复性睡眠",
  maintenance_vs_nonrestorative = "睡眠维持困难/非恢复性睡眠"
)

write_reports <- function(estimates, joint, pairwise, diagnostics, sync_summary, correlation_summary, edge_table, mi_qc) {
  main_ids <- c("hrs_5state_joint_domains", "elsa_4state_w4_w10_joint_domains")
  en <- c(
    "# Joint insomnia symptom-domain multistate analysis",
    "",
    "## Status and scope",
    "",
    "This is a **post-hoc secondary analysis** conducted after inspection of the original multistate results. The three domains, the four shared non-fatal transitions, and the hierarchical contrast family were prospectively fixed for this post-hoc extension before inspecting estimates from the new joint models. They were not preregistered confirmatory analyses.",
    "",
    "A statistically significant estimate for one domain and a non-significant estimate for another is not evidence that the domains differ. Between-domain claims below are governed by direct covariance-aware contrasts from the same joint model.",
    "",
    "## Methods",
    "",
    paste0("Twenty imputations (10 FCS iterations; seeds 20260909 for HRS and 20260910 for ELSA) were generated for missing baseline M4 adjustment covariates. The synchronised binary indicators `insomnia_onset_sync`, `insomnia_maintenance_sync`, and `insomnia_nonrestorative_sync` were protected from imputation and entered as predictors; their deterministic total score was excluded. Each imputed data set was analysed using one continuous-time panel-observed `msm` model containing all three domains. Domain effects were transition-specific, whereas non-domain M4 nuisance effects were constrained to be common across allowed transitions, matching the original primary-model convention. HRS used the five-state model with exact death times; ELSA W4-W10 used the four-state living-state model. A death-censored HRS four-state model was treated only as a structural-comparability sensitivity."),
    "",
    "The full domain coefficient vector and within-imputation covariance matrix were pooled with Rubin's rules. Scalar estimates and pairwise contrasts use Rubin t inference. Multivariate tests use the finite-imputation D1 F approximation; a Rubin total-covariance large-sample Wald chi-square is retained as a sensitivity statistic. Every inferential covariance had to be finite, full-rank and positive definite; no pseudoinverse-based test was permitted. The primary secondary-extension gate is the 8-df equality test across all three domains and four shared non-fatal transitions. If that gate is open, four 2-df per-transition tests and 12 frozen pairwise comparisons are interpreted; the pairwise family uses Benjamini-Hochberg adjustment within cohort. Per-domain 4-df association tests do not establish differences between domains.",
    "",
    "Model and pooling conventions follow the [msm documentation](https://chjackson.github.io/msm/reference/msm.html), [mice pooling documentation](https://amices.org/mice/reference/pool.html), and [mice D1 documentation](https://amices.org/mice/reference/D1.html).",
    "",
    "## Synchronised-exposure audit",
    ""
  )
  for (i in seq_len(nrow(sync_summary))) {
    z <- sync_summary[i]
    en <- c(en, sprintf(
      "- %s: %s rows among %s participants; the existing asynchronous three-item sum disagreed with the primary total score on %s rows (%s participants), whereas the synchronised triplet disagreed on %s rows (%s participants). %s rows had no earlier complete triplet.",
      z$cohort, format(z$n_rows, big.mark = ","), format(z$n_participants, big.mark = ","),
      format(z$async_score_mismatch_rows, big.mark = ","), format(z$async_score_mismatch_participants, big.mark = ","),
      format(z$sync_score_mismatch_rows, big.mark = ","), format(z$sync_score_mismatch_participants, big.mark = ","),
      format(z$rows_without_prior_complete_triplet, big.mark = ",")
    ))
  }
  en <- c(en, "", "## Collinearity and analysis-support audit", "")
  for (i in seq_len(nrow(correlation_summary))) {
    z <- correlation_summary[i]
    en <- c(en, sprintf(
      "- %s: maximum absolute row-level correlation among the synchronised binary domains %.3f; three-domain correlation-matrix condition number %.2f.",
      z$cohort, z$max_absolute_pairwise_correlation, z$correlation_matrix_condition_number
    ))
  }
  for (id in main_ids) {
    z <- edge_table[analysis_id == id & death_transition == 0L]
    if (!nrow(z)) next
    event_text <- paste(sprintf("%s=%s", z$transition, format(z$events, big.mark = ",")), collapse = "; ")
    en <- c(en, sprintf(
      "- %s: %s records from %s participants; observed transition counts %s.",
      unname(analysis_display[id]), format(z$n_records[1], big.mark = ","),
      format(z$n_participants[1], big.mark = ","), event_text
    ))
  }
  en <- c(en, "", "## Multiple-imputation diagnostics", "")
  for (cohort_label in c("HRS", "ELSA")) {
    z <- mi_qc$summary[cohort == cohort_label]
    count <- mi_qc$event_counts[cohort == cohort_label, total_logged_notices][1]
    en <- c(en, sprintf(
      "- %s: %d imputed variables, 20 chains and 10 iterations; maximum coda PSRF point estimate %.3f (upper limit %.3f); %s logged predictor-removal notices.",
      cohort_label, nrow(z), max(z$coda_psrf_point_full_trace, na.rm = TRUE),
      max(z$coda_psrf_upper_full_trace, na.rm = TRUE), format(count, big.mark = ",")
    ))
  }
  if (nrow(mi_qc$events)) {
    event_detail <- mi_qc$events[, .(notices = sum(notices)), by = .(cohort, removed_predictors)][order(cohort, -notices)]
    en <- c(en, paste0(
      "- HRS notices were local removals of deterministic/redundant outcome auxiliaries: ",
      paste(sprintf("%s (%s)", event_detail$removed_predictors, format(event_detail$notices, big.mark = ",")), collapse = "; "),
      ". All targeted covariates were complete after imputation."
    ))
  }
  en <- c(en,
    "- The installed `mice::convergence()` attempted to call the unavailable optional `rstan` package; the original error CSVs are retained. The supplied trace figure and classical/coda PSRF table are FCS chain-mean diagnostics, not MCMC diagnostics for multistate-model parameters."
  )
  en <- c(en, "", "## Formal domain-difference results", "")
  for (id in main_ids[main_ids %in% joint$analysis_id]) {
    global <- joint[analysis_id == id & family == "common_nonfatal" & test_scope == "global_domain_equality"][1]
    sig <- pairwise[analysis_id == id & family == "common_nonfatal" & global_gate_pass == 1L & p_bh < .05]
    model_name <- unname(analysis_display[id])
    if (!nrow(global)) next
    if (global$global_gate_pass == 1L) {
      sentence <- sprintf(
        "- **%s:** the global domain-equality gate was open (D1 F(%d, %.1f) = %.2f, P = %s).",
        model_name, global$numerator_df, global$denominator_df, global$F, format_p(global$p_value)
      )
      if (nrow(sig)) {
        detail <- paste(vapply(seq_len(nrow(sig)), function(i) sprintf(
          "%s, %s: ratio of transition-intensity ratios %.2f (95%% CI %.2f-%.2f; BH q=%s)",
          unname(transition_display[sig$transition[i]]), gsub("_", " ", sig$contrast[i]),
          sig$ratio_of_transition_intensity_ratios[i], sig$conf_low[i], sig$conf_high[i], format_p(sig$p_bh[i])
        ), character(1)), collapse = "; ")
        sentence <- paste(sentence, "Gate-supported pairwise findings:", detail)
      } else sentence <- paste(sentence, "No pairwise contrast remained below BH q=0.05.")
    } else {
      sentence <- sprintf(
        "- **%s:** there was no evidence that domain effects differed jointly across the four shared non-fatal transitions (D1 F(%d, %.1f) = %.2f, P = %s). Per-edge and pairwise estimates are descriptive because the global gate was closed.",
        model_name, global$numerator_df, global$denominator_df, global$F, format_p(global$p_value)
      )
    }
    en <- c(en, sentence)
  }
  hrs_global <- joint[analysis_id == "hrs_5state_joint_domains" & family == "common_nonfatal" & test_scope == "global_domain_equality"]
  elsa_global <- joint[analysis_id == "elsa_4state_w4_w10_joint_domains" & family == "common_nonfatal" & test_scope == "global_domain_equality"]
  if (nrow(hrs_global) && nrow(elsa_global) && hrs_global$global_gate_pass[1] == 1L && elsa_global$global_gate_pass[1] == 0L) {
    en <- c(en,
      "",
      "Overall, HRS supported transition-dependent conditional domain heterogeneity, with lower sleep-maintenance coefficients than sleep-onset or nonrestorative coefficients on the gate-supported S0->S1 and S2->S3 contrasts. ELSA did not provide corresponding domain-heterogeneity evidence. This pattern does not identify a single cross-cohort symptom driver, and separate cohort P values are not a formal cohort-by-domain interaction test."
    )
  }
  en <- c(en, "", "## Domain-association tests", "")
  assoc <- joint[analysis_id %in% main_ids & test_scope == "per_domain_joint_association"]
  if (nrow(assoc)) {
    for (i in seq_len(nrow(assoc))) {
      en <- c(en, sprintf(
        "- %s, %s: D1 F(%d, %.1f) = %.2f, P = %s. This tests association across transitions, not superiority over another domain.",
        unname(analysis_display[assoc$analysis_id[i]]), unname(domain_labels[assoc$domain[i]]),
        assoc$numerator_df[i], assoc$denominator_df[i], assoc$F[i], format_p(assoc$p_value[i])
      ))
    }
  }
  elsa_maintenance <- assoc[
    analysis_id == "elsa_4state_w4_w10_joint_domains" & domain == "insomnia_maintenance_sync"
  ]
  if (nrow(elsa_maintenance) && nrow(elsa_global) && elsa_global$global_gate_pass[1] == 0L) {
    en <- c(en, sprintf(
      "- ELSA illustrates the distinction: sleep-maintenance difficulty was associated jointly across transitions (4-df P = %s), but the direct domain-equality gate was closed (8-df P = %s). This cannot be reported as superiority over the other domains.",
      format_p(elsa_maintenance$p_value[1]), format_p(elsa_global$p_value[1])
    ))
  }
  en <- c(en, "", "## Exploratory and sensitivity analyses", "")
  death <- joint[analysis_id == "hrs_5state_joint_domains" & family == "hrs_death_exploratory" & test_scope == "global_domain_equality"]
  if (nrow(death)) en <- c(en, sprintf(
    "- HRS death-edge domain heterogeneity was exploratory: D1 F(%d, %.1f) = %.2f, P = %s.",
    death$numerator_df[1], death$denominator_df[1], death$F[1], format_p(death$p_value[1])
  ))
  death_sig <- pairwise[
    analysis_id == "hrs_5state_joint_domains" & family == "hrs_death_exploratory" &
      global_gate_pass == 1L & p_bh < .05
  ]
  if (nrow(death_sig)) {
    death_detail <- paste(vapply(seq_len(nrow(death_sig)), function(i) sprintf(
      "%s, %s %.2f (95%% CI %.2f-%.2f; BH q=%s)",
      unname(transition_display[death_sig$transition[i]]), gsub("_", " ", death_sig$contrast[i]),
      death_sig$ratio_of_transition_intensity_ratios[i], death_sig$conf_low[i], death_sig$conf_high[i],
      format_p(death_sig$p_bh[i])
    ), character(1)), collapse = "; ")
    en <- c(en, paste0("  Exploratory gate/BH-supported death-edge contrasts: ", death_detail, "."))
  }
  sens <- joint[analysis_id == "hrs_4state_death_censored_joint_domains" & family == "common_nonfatal" & test_scope == "global_domain_equality"]
  if (nrow(sens)) {
    sens_sig <- pairwise[
      analysis_id == "hrs_4state_death_censored_joint_domains" & family == "common_nonfatal" &
        global_gate_pass == 1L & p_bh < .05
    ]
    en <- c(en, sprintf(
      "- The HRS four-state death-censored sensitivity gave a global D1 F(%d, %.1f) = %.2f, P = %s; %s. It is a structural-comparability sensitivity and does not replace the HRS five-state analysis. Cross-model attenuation was assessed descriptively, not by a formal paired model-difference test.",
      sens$numerator_df[1], sens$denominator_df[1], sens$F[1], format_p(sens$p_value[1]),
      if (nrow(sens_sig)) sprintf("%d pairwise contrast(s) remained below BH q=0.05", nrow(sens_sig)) else "no pairwise contrast remained below BH q=0.05"
    ))
  }
  en <- c(en, "", "## Fit acceptance and interpretation limits", "")
  diag_summary <- diagnostics[, .(
    accepted = sum(accepted == 1L), total = .N,
    maximum_elapsed_minutes = max(elapsed_seconds, na.rm = TRUE) / 60,
    maximum_domain_covariance_condition_number = max(domain_covariance_condition_number, na.rm = TRUE)
  ), by = analysis_id]
  for (i in seq_len(nrow(diag_summary))) en <- c(en, sprintf(
    "- %s: %d/%d accepted fits; longest fit %.1f min; maximum domain-covariance condition number %.1f.",
    unname(analysis_display[diag_summary$analysis_id[i]]), diag_summary$accepted[i], diag_summary$total[i],
    diag_summary$maximum_elapsed_minutes[i], diag_summary$maximum_domain_covariance_condition_number[i]
  ))
  en <- c(en,
    "",
    "All fitted components passed convergence, standard-error, finite-coefficient, full free-Hessian, and domain-covariance positive-definiteness gates. Four-state fits produced the expected warning for repeated observations after the living multimorbidity state became absorbing; all fits also noted that `msm` output helpers do not display polynomial contrasts for ordered nuisance factors. The domain block was extracted directly from fitted coefficients and covariance matrices, so these warnings did not alter domain inference.",
    "",
    "Conditional domain coefficients compare symptom presence versus absence while holding the other two domains fixed and assume additivity on the log-intensity scale. They are not marginal total effects. In particular, an inverse conditional coefficient should not be interpreted as protective; suppression or phenotype differences remain possible. HRS and ELSA symptom questions use harmonised binary indicators but are not measurement-invariant. ELSA W4-W10 lacks comparable death ascertainment, so cross-cohort consistency is judged only on the four shared non-fatal transitions. A non-significant contrast is not evidence of equivalence, and this analysis does not replace the cumulative symptom-burden analysis.",
    "",
    "All participant-level repeated observations were analysed within the multistate likelihood; the participant is the independent unit. Baseline-covariate imputation uses outcome and follow-up summaries but remains an approximation for a time-varying exposure analysis."
  )

  zh <- c(
    "# 失眠症状域联合多状态分析",
    "",
    "## 分析定位",
    "",
    "本分析是在查看原始多状态结果后开展的**事后次级扩展**。三个症状域、四条两队列共有的非致死转移及层级检验族，均在查看本次新联合模型估计值之前冻结；它们不是原研究预注册的确认性分析。",
    "",
    "某一症状域显著而另一症状域不显著，并不等于两域效应存在差异。域间结论只依据同一联合模型、使用完整协方差的直接对比。",
    "",
    "## 方法",
    "",
    "对主 M4 基线协变量缺失进行20重插补（10轮FCS；HRS种子20260909，ELSA种子20260910）。同步构造的 `insomnia_onset_sync`（入睡困难）、`insomnia_maintenance_sync`（睡眠维持困难）和 `insomnia_nonrestorative_sync`（非恢复性睡眠）均作为受保护预测变量而不被插补，且不把三域的确定性总分同时纳入插补模型。每个插补数据集拟合同一个连续时间、面板观测的 `msm` 模型：三域均采用转移特异系数，其他 M4 协变量按原主模型约束为各允许边共用。HRS主分析为含死亡的五状态模型；ELSA W4-W10为四状态非致死模型；HRS死亡删失四状态模型仅作结构可比性敏感性分析。",
    "",
    "完整域系数向量和插补内协方差按 Rubin 规则池化。标量估计和两两对比使用 Rubin t 推断；多参数检验使用有限插补 D1 F 近似，并保留 Rubin 总协方差的大样本 Wald 卡方作为敏感性统计量。所有推断协方差必须有限、满秩且正定；不允许用伪逆继续检验。事后次级扩展的首要门控为三域在四条共有非致死边上的8自由度整体等同性检验；只有整体门开启时才解释每边2自由度检验和冻结的12个两两比较，队列内两两比较采用BH校正。每域跨四边的4自由度检验只回答“该域是否与任一转移相关”，不能证明域间不同。",
    "",
    "方法依据：[msm文档](https://chjackson.github.io/msm/reference/msm.html)、[mice池化文档](https://amices.org/mice/reference/pool.html)和[mice D1文档](https://amices.org/mice/reference/D1.html)。",
    "",
    "## 同步暴露审计",
    ""
  )
  for (i in seq_len(nrow(sync_summary))) {
    z <- sync_summary[i]
    zh <- c(zh, sprintf(
      "- %s：共%s行、%s人；原异步三域之和与主总分不一致%s行（%s人），同步三联体与主总分不一致%s行（%s人）；%s行此前没有完整三联体。",
      z$cohort, format(z$n_rows, big.mark = ","), format(z$n_participants, big.mark = ","),
      format(z$async_score_mismatch_rows, big.mark = ","), format(z$async_score_mismatch_participants, big.mark = ","),
      format(z$sync_score_mismatch_rows, big.mark = ","), format(z$sync_score_mismatch_participants, big.mark = ","),
      format(z$rows_without_prior_complete_triplet, big.mark = ",")
    ))
  }
  zh <- c(zh, "", "## 共线性与分析信息量审计", "")
  for (i in seq_len(nrow(correlation_summary))) {
    z <- correlation_summary[i]
    zh <- c(zh, sprintf(
      "- %s：同步二元症状域的行级两两相关绝对值最大为%.3f，三域相关矩阵条件数为%.2f。",
      z$cohort, z$max_absolute_pairwise_correlation, z$correlation_matrix_condition_number
    ))
  }
  for (id in main_ids) {
    z <- edge_table[analysis_id == id & death_transition == 0L]
    if (!nrow(z)) next
    event_text <- paste(sprintf("%s=%s", z$transition, format(z$events, big.mark = ",")), collapse = "；")
    model_name <- if (id == "hrs_5state_joint_domains") "HRS五状态模型" else "ELSA W4-W10四状态模型"
    zh <- c(zh, sprintf(
      "- %s：%s条记录、%s名参与者；观察到的转移数为%s。",
      model_name, format(z$n_records[1], big.mark = ","),
      format(z$n_participants[1], big.mark = ","), event_text
    ))
  }
  zh <- c(zh, "", "## 多重插补诊断", "")
  for (cohort_label in c("HRS", "ELSA")) {
    z <- mi_qc$summary[cohort == cohort_label]
    count <- mi_qc$event_counts[cohort == cohort_label, total_logged_notices][1]
    zh <- c(zh, sprintf(
      "- %s：插补变量%d个、20条链、10轮；coda PSRF点估计最大值%.3f（上限%.3f）；预测变量移除提示%s条。",
      cohort_label, nrow(z), max(z$coda_psrf_point_full_trace, na.rm = TRUE),
      max(z$coda_psrf_upper_full_trace, na.rm = TRUE), format(count, big.mark = ",")
    ))
  }
  if (nrow(mi_qc$events)) {
    event_detail <- mi_qc$events[, .(notices = sum(notices)), by = .(cohort, removed_predictors)][order(cohort, -notices)]
    zh <- c(zh, paste0(
      "- HRS提示均为在局部插补方程中移除确定性或冗余结局辅助变量：",
      paste(sprintf("%s（%s条）", event_detail$removed_predictors, format(event_detail$notices, big.mark = ",")), collapse = "；"),
      "。插补后所有目标协变量均无缺失。"
    ))
  }
  zh <- c(zh,
    "- 本机安装的 `mice::convergence()` 会调用未安装的可选包 `rstan`；原错误CSV已保留。另行提供的轨迹图、经典PSRF和coda PSRF是FCS链均值诊断，不是多状态模型参数的MCMC诊断。"
  )
  zh <- c(zh, "", "## 域间正式比较", "")
  for (id in main_ids[main_ids %in% joint$analysis_id]) {
    global <- joint[analysis_id == id & family == "common_nonfatal" & test_scope == "global_domain_equality"][1]
    sig <- pairwise[analysis_id == id & family == "common_nonfatal" & global_gate_pass == 1L & p_bh < .05]
    model_name <- if (id == "hrs_5state_joint_domains") "HRS五状态模型" else "ELSA W4-W10四状态模型"
    if (!nrow(global)) next
    if (global$global_gate_pass == 1L) {
      sentence <- sprintf("- **%s：**整体域等同性门开启（D1 F(%d, %.1f)=%.2f，P=%s）。", model_name, global$numerator_df, global$denominator_df, global$F, format_p(global$p_value))
      if (nrow(sig)) {
        detail <- paste(vapply(seq_len(nrow(sig)), function(i) sprintf(
          "%s，%s：转移强度比之比%.2f（95%%CI %.2f-%.2f；BH q=%s）",
          unname(transition_display_zh[sig$transition[i]]), unname(pairwise_display_zh[sig$contrast[i]]),
          sig$ratio_of_transition_intensity_ratios[i], sig$conf_low[i], sig$conf_high[i], format_p(sig$p_bh[i])
        ), character(1)), collapse = "；")
        sentence <- paste(sentence, "门控支持的两两结果：", detail)
      } else sentence <- paste(sentence, "没有两两对比达到BH q<0.05。")
    } else {
      sentence <- sprintf("- **%s：**未发现三域在四条共有非致死转移上的关联存在整体差异（D1 F(%d, %.1f)=%.2f，P=%s）。因整体门关闭，每边和两两估计仅作描述。", model_name, global$numerator_df, global$denominator_df, global$F, format_p(global$p_value))
    }
    zh <- c(zh, sentence)
  }
  if (nrow(hrs_global) && nrow(elsa_global) && hrs_global$global_gate_pass[1] == 1L && elsa_global$global_gate_pass[1] == 0L) {
    zh <- c(zh,
      "",
      "总体而言，HRS支持随转移而变化的条件性域异质性：在通过门控的S0→S1和S2→S3对比中，睡眠维持困难系数低于入睡困难或非恢复性睡眠系数；ELSA未提供相应的域异质性证据。因此，本分析没有识别出一个跨队列一致的单一症状驱动因素；分别比较两个队列的P值也不构成正式的队列×症状域交互检验。"
    )
  }
  zh <- c(zh, "", "## 各域跨转移关联检验", "")
  if (nrow(assoc)) {
    for (i in seq_len(nrow(assoc))) {
      model_name <- if (assoc$analysis_id[i] == "hrs_5state_joint_domains") "HRS五状态模型" else "ELSA W4-W10四状态模型"
      zh <- c(zh, sprintf(
        "- %s，%s：D1 F(%d, %.1f)=%.2f，P=%s。该检验回答跨转移关联，不证明该域优于其他域。",
        model_name, unname(domain_labels_zh[assoc$domain[i]]), assoc$numerator_df[i],
        assoc$denominator_df[i], assoc$F[i], format_p(assoc$p_value[i])
      ))
    }
  }
  if (nrow(elsa_maintenance) && nrow(elsa_global) && elsa_global$global_gate_pass[1] == 0L) {
    zh <- c(zh, sprintf(
      "- ELSA清楚说明了两类检验的区别：睡眠维持困难的四边联合关联检验P=%s，但域等同性整体门P=%s并未开启，不能据此声称其优于另外两个域。",
      format_p(elsa_maintenance$p_value[1]), format_p(elsa_global$p_value[1])
    ))
  }
  zh <- c(zh, "", "## 探索与敏感性分析", "")
  if (nrow(death)) zh <- c(zh, sprintf("- HRS死亡边域异质性仅作探索：D1 F(%d, %.1f)=%.2f，P=%s。", death$numerator_df[1], death$denominator_df[1], death$F[1], format_p(death$p_value[1])))
  if (nrow(death_sig)) {
    death_detail <- paste(vapply(seq_len(nrow(death_sig)), function(i) sprintf(
      "%s，%s %.2f（95%%CI %.2f-%.2f；BH q=%s）",
      unname(transition_display_zh[death_sig$transition[i]]), unname(pairwise_display_zh[death_sig$contrast[i]]),
      death_sig$ratio_of_transition_intensity_ratios[i], death_sig$conf_low[i], death_sig$conf_high[i],
      format_p(death_sig$p_bh[i])
    ), character(1)), collapse = "；")
    zh <- c(zh, paste0("  探索性死亡边中经门控和BH支持的对比：", death_detail, "。"))
  }
  if (nrow(sens)) {
    sens_sig <- pairwise[
      analysis_id == "hrs_4state_death_censored_joint_domains" & family == "common_nonfatal" &
        global_gate_pass == 1L & p_bh < .05
    ]
    zh <- c(zh, sprintf(
      "- HRS死亡删失四状态敏感性分析：整体D1 F(%d, %.1f)=%.2f，P=%s；%s。该结果不替代HRS五状态主分析，模型间效应减弱仅作描述，未进行正式配对模型差异检验。",
      sens$numerator_df[1], sens$denominator_df[1], sens$F[1], format_p(sens$p_value[1]),
      if (nrow(sens_sig)) sprintf("%d个两两对比达到BH q<0.05", nrow(sens_sig)) else "没有两两对比达到BH q<0.05"
    ))
  }
  zh <- c(zh,
    "",
    "## 拟合验收与解释边界",
    ""
  )
  for (i in seq_len(nrow(diag_summary))) zh <- c(zh, sprintf(
    "- %s：%d/%d个拟合通过；单次最长%.1f分钟；域协方差条件数最大%.1f。",
    if (diag_summary$analysis_id[i] == "hrs_5state_joint_domains") "HRS五状态模型" else if (diag_summary$analysis_id[i] == "elsa_4state_w4_w10_joint_domains") "ELSA W4-W10四状态模型" else "HRS死亡删失四状态敏感性模型",
    diag_summary$accepted[i], diag_summary$total[i], diag_summary$maximum_elapsed_minutes[i],
    diag_summary$maximum_domain_covariance_condition_number[i]
  ))
  zh <- c(zh,
    "",
    "所有组件均通过收敛、标准误、有限系数、完整自由参数Hessian正定和域协方差正定门。四状态模型中“吸收态后重复观察”的警告来自S3作为活体吸收态；所有模型还提示 `msm` 输出辅助函数不展示有序混杂因子的多项式对比。域系数和协方差均直接从拟合对象提取，因此这些警告不改变域推断。",
    "",
    "联合模型中的域系数是在另外两个域保持不变时，症状存在相对于不存在的条件关联，并假定在对数转移强度尺度上可加；它不是边际总效应。尤其不能把反向的条件系数解释为保护作用，抑制效应或表型差异仍可能存在。HRS与ELSA的二元症状变量虽已协调编码，但题干和阈值并非严格测量等值。ELSA W4-W10没有可比的死亡随访，因此跨队列一致性只评价四条共有非致死边。非显著对比不能解释为等效，本分析也不替代累积失眠症状负担分析。",
    "",
    "参与者是独立分析单位；同一参与者的重复观察由多状态似然统一处理。基线协变量插补纳入结局和随访摘要，但对时间变化暴露模型仍属于近似处理。"
  )
  writeLines(en, file.path(report_root, "joint_symptom_domain_report_en.md"), useBytes = TRUE)
  writeLines(zh, file.path(report_root, "joint_symptom_domain_report_zh.md"), useBytes = TRUE)
}

render_figures <- function(estimates, joint, pairwise) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE) ||
      !requireNamespace("svglite", quietly = TRUE) || !requireNamespace("ragg", quietly = TRUE)) {
    stop("ggplot2, patchwork, svglite and ragg are required for R-only figure export")
  }
  library(ggplot2)
  library(patchwork)
  palette <- c(
    insomnia_onset_sync = "#2F6B9A",
    insomnia_maintenance_sync = "#D08726",
    insomnia_nonrestorative_sync = "#2A8C82"
  )
  theme_joint <- function() {
    theme_classic(base_size = 7.2, base_family = "Helvetica") +
      theme(
        axis.line = element_line(linewidth = .35, colour = "#222222"),
        axis.ticks = element_line(linewidth = .35, colour = "#222222"),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7),
        legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
        legend.justification = "left", legend.title = element_text(face = "bold", size = 6.5),
        legend.text = element_text(size = 6.5), legend.spacing.y = grid::unit(1, "mm"),
        plot.title = element_text(face = "bold", size = 8.2), plot.subtitle = element_text(size = 6.5, colour = "#444444"),
        plot.tag = element_text(face = "bold", size = 9), panel.spacing = grid::unit(4, "mm")
      )
  }
  main_ids <- c("hrs_5state_joint_domains", "elsa_4state_w4_w10_joint_domains")
  est <- copy(estimates[analysis_id %in% main_ids & death_transition == 0L])
  est[, model := factor(unname(analysis_display[analysis_id]), levels = unname(analysis_display[main_ids]))]
  est[, transition_label := factor(unname(transition_display[transition]), levels = rev(unname(transition_display[c("S0->S1", "S0->S2", "S1->S3", "S2->S3")]))) ]
  est[, domain_label := factor(domain_label, levels = unname(domain_labels))]
  fwrite(est, file.path(source_root, "figure_joint_domain_estimates_source.csv"), na = "")
  p_a <- ggplot(est, aes(x = transition_intensity_ratio, y = transition_label, colour = domain, shape = domain)) +
    geom_vline(xintercept = 1, linewidth = .35, linetype = "22", colour = "#777777") +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0, linewidth = .45, position = position_dodge(width = .55), orientation = "y") +
    geom_point(size = 1.8, position = position_dodge(width = .55)) +
    facet_wrap(~model, ncol = 1) +
    scale_x_log10() +
    scale_colour_manual(values = palette, breaks = names(domain_labels), labels = unname(domain_labels), name = "Symptom domain") +
    scale_shape_manual(values = c(16, 17, 15), breaks = names(domain_labels), labels = unname(domain_labels), name = "Symptom domain") +
    labs(
      x = "Conditional transition-intensity ratio (95% CI)", y = NULL,
      title = "Joint symptom-domain associations",
      subtitle = "All three domains entered simultaneously; four shared non-fatal transitions"
    ) + theme_joint()

  pw <- copy(pairwise[analysis_id %in% main_ids & family == "common_nonfatal"])
  pw[, model := factor(unname(analysis_display[analysis_id]), levels = unname(analysis_display[main_ids]))]
  pw[, transition_label := factor(unname(transition_display[transition]), levels = rev(unname(transition_display[c("S0->S1", "S0->S2", "S1->S3", "S2->S3")]))) ]
  pair_labels <- c(
    onset_vs_maintenance = "Onset / maintenance",
    onset_vs_nonrestorative = "Onset / nonrestorative",
    maintenance_vs_nonrestorative = "Maintenance / nonrestorative"
  )
  pw[, contrast_label := factor(unname(pair_labels[contrast]), levels = unname(pair_labels))]
  pw[, supported := factor(
    fifelse(global_gate_pass == 1L & p_bh < .05, "Gate open; BH q < 0.05", "Not gate-supported"),
    levels = c("Gate open; BH q < 0.05", "Not gate-supported")
  )]
  fwrite(pw, file.path(source_root, "figure_pairwise_domain_contrasts_source.csv"), na = "")
  p_b <- ggplot(pw, aes(x = ratio_of_transition_intensity_ratios, y = transition_label, colour = contrast_label, shape = supported)) +
    geom_vline(xintercept = 1, linewidth = .35, linetype = "22", colour = "#777777") +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0, linewidth = .45, position = position_dodge(width = .58), orientation = "y") +
    geom_point(size = 1.7, position = position_dodge(width = .58)) +
    facet_wrap(~model, ncol = 1) + scale_x_log10() +
    scale_colour_manual(values = c("#6A51A3", "#D95F59", "#4D7770"), name = "Direct contrast") +
    scale_shape_manual(
      values = c(16, 1), limits = c("Gate open; BH q < 0.05", "Not gate-supported"),
      drop = FALSE, name = "Gate/BH status"
    ) +
    labs(
      x = "Ratio of transition-intensity ratios (95% CI)", y = NULL,
      title = "Direct pairwise domain contrasts",
      subtitle = "Inference uses covariance from the same joint model; BH applied to 12 frozen contrasts per cohort"
    ) + theme_joint()

  sens <- copy(estimates[analysis_id %in% c("hrs_5state_joint_domains", "hrs_4state_death_censored_joint_domains") & death_transition == 0L])
  sens[, model := factor(
    fifelse(analysis_id == "hrs_5state_joint_domains", "Five-state (death competing state)", "Four-state (death censored)"),
    levels = c("Five-state (death competing state)", "Four-state (death censored)")
  )]
  sens[, transition_label := factor(unname(transition_display[transition]), levels = rev(unname(transition_display[c("S0->S1", "S0->S2", "S1->S3", "S2->S3")]))) ]
  sens[, domain_label := factor(domain_label, levels = unname(domain_labels))]
  fwrite(sens, file.path(source_root, "figure_hrs4_sensitivity_source.csv"), na = "")
  p_c <- ggplot(sens, aes(x = transition_intensity_ratio, y = transition_label, colour = domain, shape = model)) +
    geom_vline(xintercept = 1, linewidth = .35, linetype = "22", colour = "#777777") +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0, linewidth = .4, position = position_dodge(width = .62), orientation = "y") +
    geom_point(size = 1.65, position = position_dodge(width = .62)) +
    scale_x_log10() +
    scale_colour_manual(values = palette, breaks = names(domain_labels), labels = unname(domain_labels), name = "Symptom domain") +
    scale_shape_manual(values = c(16, 1), name = "HRS specification") +
    labs(
      x = "Conditional transition-intensity ratio (95% CI)", y = NULL,
      title = "HRS structural-comparability sensitivity",
      subtitle = "Death-censored four-state results do not replace the five-state analysis"
    ) + guides(colour = "none") + theme_joint()

  p_a_hrs <- p_a + est[analysis_id == "hrs_5state_joint_domains"] +
    facet_null() + labs(title = "HRS joint symptom-domain associations", subtitle = "Five-state model; four non-fatal transitions shown")
  p_a_elsa <- p_a + est[analysis_id == "elsa_4state_w4_w10_joint_domains"] +
    facet_null() + labs(title = "ELSA joint symptom-domain associations", subtitle = "W4-W10 four-state model")
  p_b_hrs <- p_b + pw[analysis_id == "hrs_5state_joint_domains"] +
    facet_null() + labs(title = "HRS direct pairwise domain contrasts", subtitle = "Filled symbols: global gate open and BH q < 0.05")
  p_b_elsa <- p_b + pw[analysis_id == "elsa_4state_w4_w10_joint_domains"] +
    facet_null() + guides(shape = "none") +
    labs(title = "ELSA direct pairwise domain contrasts", subtitle = "Global domain-equality gate closed; all contrasts descriptive")
  figure <- wrap_plots(
    p_a_hrs, p_a_elsa, p_b_hrs, p_b_elsa, p_c,
    ncol = 1, guides = "collect", heights = c(1, 1, 1, 1, .9)
  ) + plot_annotation(tag_levels = list(c("A1", "A2", "B1", "B2", "C"))) &
    theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left")
  base <- file.path(figure_root, "fig_joint_symptom_domain_comparison")
  width_in <- 183 / 25.4
  height_in <- 250 / 25.4
  alignment_script <- "private/render_tools/nature-figure/scripts/panel_alignment.R"
  alignment_audit <- "private/render_tools/nature-figure/scripts/audit_panel_alignment.py"
  alignment_python <- "python"
  if (file.exists(alignment_script)) {
    source(alignment_script, local = environment())
    write_patchwork_panel_layout <- write_patchwork_panel_layout_resolved
    require_patchwork_panel_alignment(
      figure, manifest_path = paste0(base, ".alignment-layout.json"),
      report_path = paste0(base, ".alignment.json"), overlay_svg = paste0(base, ".alignment.svg"),
      audit_script = alignment_audit, python = alignment_python,
      width_in = width_in, height_in = height_in, tolerance_pt = 1.5,
      gutter_tolerance_pt = 1.5, strict = TRUE
    )
  }
  export_figure_set(figure, base, width_in, height_in, 600)
  caption_en <- c(
    "Fig. | Joint comparison of insomnia symptom domains across cardiovascular multistate transitions.",
    "(A1-A2) Conditional transition-intensity ratios (95% confidence intervals) from HRS and ELSA models simultaneously containing all three binary symptom domains. (B1-B2) Direct ratios of transition-intensity ratios for the frozen pairwise domain contrasts. Filled symbols identify contrasts supported after the cohort-level global gate and Benjamini-Hochberg adjustment; open symbols are not gate-supported and should not be interpreted from separate significance patterns. (C) HRS five-state estimates compared descriptively with the death-censored four-state structural-comparability sensitivity. Estimates combine 20 imputations using Rubin's rules; scalar intervals use Rubin t degrees of freedom. This figure reports a post-hoc secondary extension whose contrasts were fixed before the new joint models were fitted."
  )
  caption_zh <- c(
    "图 | 三个失眠症状域与心血管多状态转移的联合比较。",
    "A1-A2，HRS和ELSA中同时纳入三个二元症状域的联合模型所估计的条件转移强度比及95%置信区间。B1-B2，冻结的域间两两对比，效应量为转移强度比之比；实心符号表示整体门控开启且经Benjamini-Hochberg校正后得到支持，空心符号表示未获门控支持，不能依据分别显著或不显著作域间判断。C，HRS五状态模型与死亡删失四状态结构可比性敏感性分析的描述性比较。结果按Rubin规则合并20重插补，标量区间采用Rubin t自由度。本图属于事后次级扩展，所有对比均在本次新联合模型拟合前冻结。"
  )
  writeLines(c(caption_en, "", caption_zh), file.path(figure_root, "fig_joint_symptom_domain_comparison_caption.md"), useBytes = TRUE)

  death_est <- copy(estimates[
    analysis_id == "hrs_5state_joint_domains" & death_transition == 1L
  ])
  death_pw <- copy(pairwise[
    analysis_id == "hrs_5state_joint_domains" & family == "hrs_death_exploratory"
  ])
  if (nrow(death_est) && nrow(death_pw)) {
    death_order <- c("S0->S4", "S1->S4", "S2->S4", "S3->S4")
    death_est[, transition_label := factor(
      unname(transition_display[transition]), levels = rev(unname(transition_display[death_order]))
    )]
    death_pw[, transition_label := factor(
      unname(transition_display[transition]), levels = rev(unname(transition_display[death_order]))
    )]
    death_pw[, contrast_label := factor(unname(pair_labels[contrast]), levels = unname(pair_labels))]
    death_pw[, supported := factor(
      fifelse(global_gate_pass == 1L & p_bh < .05, "Gate open; BH q < 0.05", "Not gate-supported"),
      levels = c("Gate open; BH q < 0.05", "Not gate-supported")
    )]
    fwrite(death_est, file.path(source_root, "figure_hrs_death_domain_estimates_source.csv"), na = "")
    fwrite(death_pw, file.path(source_root, "figure_hrs_death_pairwise_contrasts_source.csv"), na = "")

    p_d <- ggplot(death_est, aes(
      x = transition_intensity_ratio, y = transition_label, colour = domain, shape = domain
    )) +
      geom_vline(xintercept = 1, linewidth = .35, linetype = "22", colour = "#777777") +
      geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0, linewidth = .45, position = position_dodge(width = .55), orientation = "y") +
      geom_point(size = 1.8, position = position_dodge(width = .55)) +
      scale_x_log10() +
      scale_colour_manual(values = palette, breaks = names(domain_labels), labels = unname(domain_labels), name = "Symptom domain") +
      scale_shape_manual(values = c(16, 17, 15), breaks = names(domain_labels), labels = unname(domain_labels), name = "Symptom domain") +
      labs(
        x = "Conditional transition-intensity ratio (95% CI)", y = NULL,
        title = "HRS death-transition associations (exploratory)",
        subtitle = "All three domains entered simultaneously; no comparable ELSA death analysis"
      ) + theme_joint()

    p_e <- ggplot(death_pw, aes(
      x = ratio_of_transition_intensity_ratios, y = transition_label,
      colour = contrast_label, shape = supported
    )) +
      geom_vline(xintercept = 1, linewidth = .35, linetype = "22", colour = "#777777") +
      geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0, linewidth = .45, position = position_dodge(width = .58), orientation = "y") +
      geom_point(size = 1.7, position = position_dodge(width = .58)) +
      scale_x_log10() +
      scale_colour_manual(values = c("#6A51A3", "#D95F59", "#4D7770"), name = "Direct contrast") +
      scale_shape_manual(
        values = c(16, 1), limits = c("Gate open; BH q < 0.05", "Not gate-supported"),
        drop = FALSE, name = "Gate/BH status"
      ) +
      labs(
        x = "Ratio of transition-intensity ratios (95% CI)", y = NULL,
        title = "HRS direct death-edge domain contrasts (exploratory)",
        subtitle = "BH adjustment applied within the frozen 12-contrast death-edge family"
      ) + theme_joint()

    death_figure <- p_d / p_e + plot_annotation(tag_levels = "A") +
      plot_layout(guides = "collect", heights = c(1, 1)) &
      theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left")
    death_base <- file.path(figure_root, "fig_hrs_death_symptom_domain_exploratory")
    death_width <- 183 / 25.4
    death_height <- 180 / 25.4
    if (file.exists(alignment_script)) {
      require_patchwork_panel_alignment(
        death_figure, manifest_path = paste0(death_base, ".alignment-layout.json"),
        report_path = paste0(death_base, ".alignment.json"), overlay_svg = paste0(death_base, ".alignment.svg"),
        audit_script = alignment_audit, python = alignment_python,
        width_in = death_width, height_in = death_height, tolerance_pt = 1.5,
        gutter_tolerance_pt = 1.5, strict = TRUE
      )
    }
    export_figure_set(death_figure, death_base, death_width, death_height, 600)
    writeLines(c(
      "Supplementary Fig. | Exploratory HRS symptom-domain results for transitions to death.",
      "(A) Mutually adjusted conditional transition-intensity ratios and 95% confidence intervals for the four HRS death transitions. (B) Direct ratios of transition-intensity ratios for the three frozen pairwise domain contrasts on each death edge. Filled symbols passed the exploratory death-family global gate and BH correction; open symbols were not gate-supported. The analysis combines 20 imputations using Rubin's rules. Death-edge findings were not part of the cross-cohort primary secondary-extension family and cannot be validated against ELSA W4-W10 because comparable mortality follow-up was unavailable.",
      "",
      "补充图 | HRS转入死亡各边的探索性症状域结果。",
      "A，HRS四条死亡边中三个症状域相互调整后的条件转移强度比及95%置信区间。B，每条死亡边上三个冻结的域间两两对比，效应量为转移强度比之比；实心符号表示探索性死亡边整体门控开启且通过BH校正，空心符号表示未获门控支持。结果按Rubin规则合并20重插补。死亡边不属于跨队列事后次级扩展的主要检验族，且ELSA W4-W10缺乏可比死亡随访，不能作跨队列验证。"
    ), paste0(death_base, "_caption.md"), useBytes = TRUE)
  }
}

panels <- load_base_panels()
sync_audit <- write_sync_audit(panels)
analysis_sets <- make_analysis_sets(panels)

support <- list(); patterns <- list()
for (config in analysis_configs) {
  if (!config$analysis_id %in% names(analysis_sets)) next
  support[[config$analysis_id]] <- edge_support(analysis_sets[[config$analysis_id]], config)
  patterns[[config$analysis_id]] <- symptom_pattern_support(analysis_sets[[config$analysis_id]], config)
}
if (length(support)) fwrite(rbindlist(support, fill = TRUE), file.path(qc_root, "joint_domain_edge_support.csv"), na = "")
if (length(patterns)) fwrite(rbindlist(patterns, fill = TRUE), file.path(qc_root, "joint_domain_symptom_pattern_support.csv"), na = "")

mi_objects <- list()
mi_configs <- list(
  hrs_joint_domains = analysis_configs$hrs_5state_joint_domains,
  elsa_joint_domains = analysis_configs$elsa_4state_w4_w10_joint_domains
)
for (mi_id in names(mi_configs)) {
  config <- mi_configs[[mi_id]]
  if (is.null(config)) next
  needed <- any(vapply(analysis_configs, function(z) identical(z$mi_id, mi_id) && z$analysis_id %in% names(analysis_sets), logical(1)))
  if (!needed) next
  base_id <- if (mi_id == "hrs_joint_domains") "hrs_5state_joint_domains" else "elsa_4state_w4_w10_joint_domains"
  base_x <- prepare_analysis_data(panels[[config$cohort]], config)
  mi_objects[[mi_id]] <- prepare_mi(base_x, config)
}

if (run_stage %in% c("full", "fit")) {
  for (config in analysis_configs) {
    if (!config$analysis_id %in% names(analysis_sets)) next
    fit_analysis(analysis_sets[[config$analysis_id]], config, mi_objects[[config$mi_id]])
  }
}

if (run_stage %in% c("full", "pool")) {
  pooled <- list()
  pooling_errors <- character()
  for (config in analysis_configs) {
    if (!config$analysis_id %in% names(analysis_sets)) next
    pooled[[config$analysis_id]] <- tryCatch(pool_analysis(config, analysis_sets[[config$analysis_id]], mi_objects[[config$mi_id]]), error = function(e) {
      pooling_errors[config$analysis_id] <<- conditionMessage(e)
      message("POOLING REJECTED for ", config$analysis_id, ": ", conditionMessage(e)); NULL
    })
  }
  pooled <- Filter(Negate(is.null), pooled)
  expected_analysis_ids <- names(analysis_configs)[names(analysis_configs) %in% names(analysis_sets)]
  if (!setequal(names(pooled), expected_analysis_ids)) {
    missing_ids <- setdiff(expected_analysis_ids, names(pooled))
    stop(
      "Complete pooling requires every analysis and all ", mi_m, " accepted components. Missing: ",
      paste(missing_ids, collapse = ", "), ". Details: ",
      paste(sprintf("%s=%s", names(pooling_errors), pooling_errors), collapse = " | ")
    )
  }
  if (length(pooled)) {
    estimates <- rbindlist(lapply(pooled, function(x) x$tests$estimates), fill = TRUE)
    joint <- rbindlist(lapply(pooled, function(x) x$tests$joint), fill = TRUE)
    pairwise <- rbindlist(lapply(pooled, function(x) x$tests$pairwise), fill = TRUE)
    diagnostics <- rbindlist(lapply(pooled, function(z) z$diagnostics), fill = TRUE)
    coefficients <- rbindlist(lapply(pooled, function(z) z$coefficients), fill = TRUE)
    covariance <- rbindlist(lapply(pooled, function(z) z$covariance), fill = TRUE)
    pooled_covariance <- rbindlist(lapply(pooled, function(z) z$pooled_covariance), fill = TRUE)
    fwrite(estimates, file.path(result_root, "symptom_domain_joint_estimates.csv"), na = "")
    fwrite(joint, file.path(result_root, "symptom_domain_joint_tests.csv"), na = "")
    fwrite(pairwise, file.path(result_root, "symptom_domain_pairwise_contrasts.csv"), na = "")
    fwrite(diagnostics, file.path(qc_root, "symptom_domain_model_diagnostics.csv"), na = "")
    fwrite(coefficients, file.path(result_root, "symptom_domain_mi_component_coefficients.csv"), na = "")
    fwrite(covariance, file.path(result_root, "symptom_domain_mi_component_covariance.csv.gz"), na = "")
    fwrite(pooled_covariance, file.path(result_root, "symptom_domain_pooled_covariance.csv"), na = "")
    edge_table <- if (length(support)) rbindlist(support, fill = TRUE) else data.table()
    model_inventory <- write_model_inventory(rbindlist(lapply(pooled, function(z) z$inventory_expected), fill = TRUE))
    mi_qc <- write_mi_diagnostics(mi_objects)
    write_reports(
      estimates, joint, pairwise, diagnostics, sync_audit$summary,
      sync_audit$correlation_summary, edge_table, mi_qc
    )
    render_figures(estimates, joint, pairwise)
    write_reproducibility_bundle(model_inventory)
  }
}

cat(sprintf("11b stage '%s' complete at %s\n", run_stage, format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")))
