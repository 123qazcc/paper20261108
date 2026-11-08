#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

stage_root <- Sys.getenv("GA_STAGE_ROOT", unset = "private/stage")
output_root <- file.path(stage_root, "multistate_output")
result_root <- file.path(output_root, "results", "markov_shards", "hrs_mi20_pool_final")
model_root <- file.path(output_root, "models", "markov_shards", "hrs_mi20_final")
log_root <- file.path(output_root, "logs", "markov_shards", "hrs_mi20_pool_final")
by_model <- file.path(result_root, "by_model")

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  fread(path)
}

fit <- read_required(file.path(result_root, "markov_fit_metrics.csv"))
hr <- read_required(file.path(result_root, "markov_hazard_ratios.csv"))
component_fit <- read_required(file.path(
  by_model, "hrs_5state_full__M4_MI20__imp01_20_mi_components_metrics.csv"
))
component_hr <- read_required(file.path(
  by_model, "hrs_5state_full__M4_MI20__imp01_20_mi_components_hr.csv"
))
source_fit <- read_required(file.path(
  result_root, "source_audit", "hrs_mi20_per_imputation_fit_metrics.csv"
))
source_sample <- read_required(file.path(
  result_root, "source_audit", "hrs_mi20_per_imputation_sample.csv"
))
source_edges <- read_required(file.path(
  result_root, "source_audit", "hrs_mi20_per_imputation_edge_counts.csv"
))
copy_manifest <- read_required(file.path(model_root, "hrs_mi20_source_copy_manifest.csv"))
destination_sha256 <- vapply(
  copy_manifest$destination_file,
  digest::digest,
  character(1L),
  algo = "sha256",
  file = TRUE,
  serialize = FALSE
)
pool_manifest <- readRDS(file.path(
  model_root, "hrs_5state_full__M4_MI20__manifest.rds"
))
failure_log <- read_required(file.path(log_root, "markov_failures.csv"))

requested_component_hr <- component_hr[is_requested_insomnia_exposure == 1L]
requested_pooled_hr <- hr[is_requested_insomnia_exposure == 1L]
allowed_transitions <- c(
  "S0->S1", "S0->S2", "S0->S4", "S1->S3",
  "S1->S4", "S2->S3", "S2->S4", "S3->S4"
)

checks <- list()
add_check <- function(check, passed, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = check,
    status = if (isTRUE(passed)) "PASS" else "FAIL",
    observed = as.character(observed),
    expected = as.character(expected)
  )
}

add_check("source_models_accepted", nrow(source_fit) == 20L && all(source_fit$validation_status == "ACCEPTED"),
  paste(nrow(source_fit), sum(source_fit$validation_status == "ACCEPTED"), sep = "/"), "20/20")
add_check("source_samples_reconstructed", nrow(source_sample) == 20L && uniqueN(source_sample$imputation) == 20L,
  paste(nrow(source_sample), uniqueN(source_sample$imputation), sep = "/"), "20 rows/20 imputations")
add_check("source_edges_reconstructed", nrow(source_edges) == 160L && uniqueN(source_edges$transition) == 8L,
  paste(nrow(source_edges), uniqueN(source_edges$transition), sep = "/"), "160 rows/8 edges")
add_check("copied_model_hashes_preserved", nrow(copy_manifest) == 20L &&
    all(file.exists(copy_manifest$destination_file)) &&
    identical(tolower(copy_manifest$sha256), tolower(unname(destination_sha256))),
  paste(nrow(copy_manifest), sum(tolower(copy_manifest$sha256) == tolower(destination_sha256)), sep = "/"),
  "20 files/20 matching SHA256 hashes")
add_check("pool_manifest_complete", identical(pool_manifest$pooling_status, "RUBIN_POOLED_COMPLETE") &&
    identical(as.integer(pool_manifest$expected_imputations), 20L) &&
    identical(as.integer(pool_manifest$fitted_imputations), 20L) &&
    length(pool_manifest$model_files) == 20L,
  paste(pool_manifest$pooling_status, pool_manifest$expected_imputations,
    pool_manifest$fitted_imputations, length(pool_manifest$model_files), sep = "/"),
  "RUBIN_POOLED_COMPLETE/20/20/20")
add_check("pooled_fit_counts", nrow(fit) == 1L && fit$m_expected == 20L && fit$m_fitted == 20L &&
    fit$m_valid_for_primary_pooling == 20L,
  if (nrow(fit)) paste(fit$m_expected, fit$m_fitted, fit$m_valid_for_primary_pooling, sep = "/") else "no row",
  "20/20/20")
add_check("component_fit_counts", nrow(component_fit) == 20L && uniqueN(component_fit$imputation) == 20L,
  paste(nrow(component_fit), uniqueN(component_fit$imputation), sep = "/"), "20/20")
add_check("component_fit_inference_valid", all(
  component_fit$status == "FITTED" & component_fit$optimizer_convergence_code == 0L &
    component_fit$standard_errors_found == 1L & component_fit$covariance_matrix_finite == 1L &
    component_fit$hessian_positive_definite == 1L
), sum(component_fit$status == "FITTED"), "20 valid component fits")
add_check("failure_log_empty", nrow(failure_log) == 0L, nrow(failure_log), "0 failures")
add_check("component_requested_hr_complete", nrow(requested_component_hr) == 160L &&
    uniqueN(requested_component_hr$imputation) == 20L &&
    setequal(unique(requested_component_hr$transition), allowed_transitions),
  paste(nrow(requested_component_hr), uniqueN(requested_component_hr$imputation),
    uniqueN(requested_component_hr$transition), sep = "/"), "160/20/8")
add_check("component_requested_hr_finite", all(
  is.finite(requested_component_hr$hazard_ratio) & requested_component_hr$hazard_ratio > 0 &
    is.finite(requested_component_hr$conf_low) & is.finite(requested_component_hr$conf_high) &
    is.finite(requested_component_hr$standard_error_log_hr) &
    requested_component_hr$standard_error_log_hr > 0
), sum(is.finite(requested_component_hr$hazard_ratio)), "160 finite estimable HRs")
add_check("pooled_requested_hr_complete", nrow(requested_pooled_hr) == 8L &&
    setequal(unique(requested_pooled_hr$transition), allowed_transitions) &&
    all(requested_pooled_hr$m_expected == 20L & requested_pooled_hr$m_fitted == 20L &
      requested_pooled_hr$pooling_status == "RUBIN_POOLED_COMPLETE"),
  paste(nrow(requested_pooled_hr), uniqueN(requested_pooled_hr$transition), sep = "/"), "8/8 pooled edges")
add_check("pooled_requested_hr_finite", all(
  is.finite(requested_pooled_hr$hazard_ratio) & requested_pooled_hr$hazard_ratio > 0 &
    is.finite(requested_pooled_hr$conf_low) & is.finite(requested_pooled_hr$conf_high) &
    is.finite(requested_pooled_hr$standard_error_log_hr) & requested_pooled_hr$standard_error_log_hr > 0
), sum(is.finite(requested_pooled_hr$hazard_ratio)), "8 finite estimable HRs")

qc <- rbindlist(checks)
fwrite(qc, file.path(result_root, "hrs_mi20_pool_validation.csv"), na = "")
fwrite(requested_pooled_hr, file.path(result_root, "hrs_mi20_pooled_exposure_hr_validation.csv"), na = "")

if (any(qc$status != "PASS")) {
  print(qc)
  stop("HRS MI20 pooled-output validation failed")
}

cat("Validation checks passed:", nrow(qc), "\n")
cat("Pooled requested-exposure edges:", nrow(requested_pooled_hr), "\n")
cat("Failure rows:", nrow(failure_log), "\n")
cat("Result root:", normalizePath(result_root, winslash = "/", mustWork = TRUE), "\n")
