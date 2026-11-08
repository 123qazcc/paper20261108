#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(msm)
})

stage_root <- Sys.getenv("GA_STAGE_ROOT", unset = "private/stage")
output_root <- file.path(stage_root, "multistate_output")
model_shards <- file.path(output_root, "models", "markov_shards")
result_shards <- file.path(output_root, "results", "markov_shards")
audit_root <- Sys.getenv(
  "GA_HRS_MI20_AUDIT_ROOT",
  unset = file.path(result_shards, "hrs_mi20_pool_final", "source_audit")
)
dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)

candidate_dirs <- list.dirs(model_shards, recursive = FALSE, full.names = TRUE)
candidate_dirs <- candidate_dirs[grepl("^v3_hrs_mi(_|$)", basename(candidate_dirs))]
candidate_files <- unlist(lapply(
  candidate_dirs,
  list.files,
  pattern = "^hrs_5state_full__M4_MI20__imp[0-9]{2}\\.rds$",
  full.names = TRUE
), use.names = FALSE)

inventory <- data.table(source_model_file = normalizePath(
  candidate_files, winslash = "/", mustWork = TRUE
))
inventory[, imputation := as.integer(sub(
  ".*__imp([0-9]{2})\\.rds$", "\\1", source_model_file
))]
inventory[, `:=`(
  source_shard = basename(dirname(source_model_file)),
  bytes = file.info(source_model_file)$size,
  modified_utc = format(
    file.info(source_model_file)$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
  )
)]
setorder(inventory, imputation, source_model_file)

duplicate_counts <- inventory[, .N, by = imputation]
if (!identical(duplicate_counts$imputation, 1:20) || any(duplicate_counts$N != 1L)) {
  fwrite(inventory, file.path(audit_root, "hrs_mi20_candidate_inventory.csv"), na = "")
  stop("Expected exactly one candidate for each HRS MI imputation 1:20")
}

allowed_edges <- data.table(
  from_state = c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L),
  to_state = c(2L, 3L, 5L, 4L, 5L, 4L, 5L, 5L),
  transition = c("T01", "T02", "T04", "T13", "T14", "T23", "T24", "T34")
)

fit_rows <- vector("list", 20L)
sample_rows <- vector("list", 20L)
edge_rows <- vector("list", 20L)
component_rows <- vector("list", 20L)

as_state_integer <- function(x) {
  if (is.factor(x)) suppressWarnings(as.integer(as.character(x))) else as.integer(x)
}

for (i in seq_len(nrow(inventory))) {
  item <- inventory[i]
  imp <- item$imputation
  envelope <- tryCatch(readRDS(item$source_model_file), error = function(e) e)
  if (inherits(envelope, "error")) {
    fit_rows[[i]] <- data.table(
      imputation = imp, object_readable = 0L, validation_status = "REJECTED",
      validation_message = conditionMessage(envelope)
    )
    next
  }
  fit <- envelope$fit
  if (!inherits(fit, "msm")) {
    fit_rows[[i]] <- data.table(
      imputation = imp, object_readable = 1L, validation_status = "REJECTED",
      validation_message = "Envelope does not contain an msm fit"
    )
    next
  }

  hessian <- if (is.matrix(fit$opt$hessian)) fit$opt$hessian else fit$hessian
  hessian_finite <- is.matrix(hessian) && all(is.finite(hessian))
  hessian_min_eigenvalue <- if (hessian_finite) {
    min(eigen((hessian + t(hessian)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  } else NA_real_
  q <- tryCatch(msm::qmatrix.msm(fit, covariates = "mean", ci = "none"), error = function(e) NULL)
  q_finite <- is.matrix(q) && all(is.finite(q))
  q_row_sum_max_abs <- if (q_finite) max(abs(rowSums(q))) else NA_real_
  death_row_max_abs <- if (q_finite && nrow(q) == 5L) max(abs(q[5L, ])) else NA_real_

  mf <- fit$data$mf
  sample_ok <- is.data.frame(mf) && all(c("(subject)", "(state)", "(time)") %in% names(mf))
  if (sample_ok) {
    panel <- data.table(
      id = as.character(mf[["(subject)"]]),
      state = as_state_integer(mf[["(state)"]]),
      time = as.numeric(mf[["(time)"]]),
      original_order = seq_len(nrow(mf))
    )
    setorder(panel, id, time, original_order)
    pairs <- panel[, if (.N >= 2L) .(
      from_state = state[-.N],
      to_observed = state[-1L],
      interval_years = time[-1L] - time[-.N]
    ) else NULL, by = id]
    sample_rows[[i]] <- panel[, .(
      n_records = .N,
      n_subjects = uniqueN(id),
      n_death_records = sum(state == 5L, na.rm = TRUE),
      n_observed_intervals = nrow(pairs),
      min_time_years = min(time, na.rm = TRUE),
      max_time_years = max(time, na.rm = TRUE),
      nonfinite_time_count = sum(!is.finite(time)),
      invalid_state_count = sum(!state %in% 1:5 | is.na(state)),
      negative_interval_count = sum(pairs$interval_years < 0, na.rm = TRUE)
    )][, imputation := imp]
    edge_rows[[i]] <- allowed_edges[, {
      source_pairs <- pairs[from_state == .BY$from_state]
      .(
        risk_intervals = nrow(source_pairs),
        risk_participants = uniqueN(source_pairs$id),
        person_years = sum(source_pairs$interval_years, na.rm = TRUE),
        events = sum(source_pairs$to_observed == .BY$to_state, na.rm = TRUE)
      )
    }, by = .(from_state, to_state, transition)][, imputation := imp]
  }

  source_result_dir <- file.path(result_shards, item$source_shard)
  component_prefix <- sprintf(
    "hrs_5state_full__M4_MI20__imp%02d_%02d_mi_components_", imp, imp
  )
  expected_component_suffixes <- c(
    "metrics", "q", "hr", "pmatrix", "sojourn", "prevalence", "expected_time"
  )
  component_files <- file.path(
    source_result_dir, "by_model",
    paste0(component_prefix, expected_component_suffixes, ".csv")
  )
  original_metrics <- component_files[expected_component_suffixes == "metrics"]
  original_metric <- if (file.exists(original_metrics)) {
    tryCatch(fread(original_metrics), error = function(e) data.table())
  } else data.table()
  metric_matches_object <- if (nrow(original_metric) == 1L) {
    identical(as.integer(original_metric$optimizer_convergence_code), as.integer(fit$opt$convergence)) &&
      identical(as.integer(original_metric$n_records), as.integer(nrow(mf))) &&
      identical(as.character(original_metric$model_data_fingerprint), as.character(envelope$model_data_fingerprint))
  } else NA
  component_rows[[i]] <- data.table(
    imputation = imp,
    source_shard = item$source_shard,
    original_component_files_present = sum(file.exists(component_files)),
    original_component_files_expected = length(component_files),
    original_metrics_present = as.integer(nrow(original_metric) == 1L),
    original_metrics_match_object = as.integer(metric_matches_object),
    original_manifest_present = as.integer(any(file.exists(file.path(
      dirname(item$source_model_file),
      sprintf("hrs_5state_full__M4_MI20__imp%02d_%02d__manifest.rds", imp, imp)
    ))))
  )

  convergence <- suppressWarnings(as.integer(fit$opt$convergence)[1L])
  found_se <- isTRUE(fit$foundse)
  covmat_finite <- is.matrix(fit$covmat) && all(is.finite(fit$covmat))
  valid <- identical(as.integer(envelope$imputation), imp) &&
    identical(convergence, 0L) && found_se && covmat_finite &&
    hessian_finite && is.finite(hessian_min_eigenvalue) && hessian_min_eigenvalue > 1e-8 &&
    q_finite && is.finite(q_row_sum_max_abs) && q_row_sum_max_abs < 1e-10 &&
    is.finite(death_row_max_abs) && death_row_max_abs < 1e-10 && sample_ok &&
    nrow(sample_rows[[i]]) == 1L && sample_rows[[i]]$nonfinite_time_count == 0L &&
    sample_rows[[i]]$invalid_state_count == 0L && sample_rows[[i]]$negative_interval_count == 0L &&
    nrow(edge_rows[[i]]) == 8L && all(edge_rows[[i]]$risk_intervals > 0L)
  fit_rows[[i]] <- data.table(
    imputation = imp,
    source_shard = item$source_shard,
    source_model_file = item$source_model_file,
    bytes = item$bytes,
    object_readable = 1L,
    envelope_imputation = as.integer(envelope$imputation),
    model_class = paste(class(fit), collapse = "|"),
    optimizer_convergence_code = convergence,
    standard_errors_found = as.integer(found_se),
    covariance_matrix_finite = as.integer(covmat_finite),
    hessian_finite = as.integer(hessian_finite),
    hessian_positive_definite = as.integer(
      is.finite(hessian_min_eigenvalue) && hessian_min_eigenvalue > 1e-8
    ),
    hessian_min_eigenvalue = hessian_min_eigenvalue,
    q_finite = as.integer(q_finite),
    q_row_sum_max_abs = q_row_sum_max_abs,
    death_row_max_abs = death_row_max_abs,
    model_data_fingerprint = as.character(envelope$model_data_fingerprint),
    spec_hash = as.character(envelope$spec_hash),
    warning_count = length(envelope$warnings),
    elapsed_seconds = as.numeric(envelope$elapsed_seconds),
    validation_status = if (valid) "ACCEPTED" else "REJECTED",
    validation_message = if (valid) {
      "Readable converged msm fit with finite inference, positive-definite Hessian, valid five-state sample, and all eight supported edges"
    } else {
      "One or more strict fit/sample/edge checks failed"
    }
  )
  rm(envelope, fit, mf, q)
  gc(verbose = FALSE)
}

fit_audit <- rbindlist(fit_rows, use.names = TRUE, fill = TRUE)
sample_audit <- rbindlist(sample_rows, use.names = TRUE, fill = TRUE)
edge_audit <- rbindlist(edge_rows, use.names = TRUE, fill = TRUE)
component_audit <- rbindlist(component_rows, use.names = TRUE, fill = TRUE)
setcolorder(sample_audit, c("imputation", setdiff(names(sample_audit), "imputation")))
setcolorder(edge_audit, c("imputation", setdiff(names(edge_audit), "imputation")))
setorder(sample_audit, imputation)
setorder(edge_audit, imputation, transition)

inventory <- merge(
  inventory,
  fit_audit[, .(imputation, validation_status)],
  by = "imputation", all.x = TRUE, sort = TRUE
)
fwrite(inventory, file.path(audit_root, "hrs_mi20_candidate_inventory.csv"), na = "")
fwrite(fit_audit, file.path(audit_root, "hrs_mi20_per_imputation_fit_metrics.csv"), na = "")
fwrite(sample_audit, file.path(audit_root, "hrs_mi20_per_imputation_sample.csv"), na = "")
fwrite(edge_audit, file.path(audit_root, "hrs_mi20_per_imputation_edge_counts.csv"), na = "")
fwrite(component_audit, file.path(audit_root, "hrs_mi20_source_component_outputs.csv"), na = "")

if (nrow(fit_audit) != 20L || any(fit_audit$validation_status != "ACCEPTED")) {
  stop("HRS MI20 source audit rejected at least one candidate model")
}
if (nrow(sample_audit) != 20L || uniqueN(sample_audit$imputation) != 20L) {
  stop("HRS MI20 source audit did not reconstruct all sample summaries")
}
if (nrow(edge_audit) != 160L || uniqueN(edge_audit$imputation) != 20L) {
  stop("HRS MI20 source audit did not reconstruct all 8 edges for all imputations")
}

cat("Accepted HRS MI models:", nrow(fit_audit), "\n")
cat("Per-imputation sample rows:", nrow(sample_audit), "\n")
cat("Per-imputation edge rows:", nrow(edge_audit), "\n")
cat("Audit root:", normalizePath(audit_root, winslash = "/", mustWork = TRUE), "\n")
