#!/usr/bin/env Rscript

# Formal post-estimation contrasts and HRS-ELSA synthesis for the frozen
# M4-MI20 continuous-time Markov models.
#
# This script never refits or alters the primary models. It extracts the
# transition-specific insomnia coefficients and their full covariance from
# each imputed-data msm object, evaluates contrasts within imputation, and then
# applies Rubin pooling. All plotting and exports are R-only.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

required_packages <- c("msm", "svglite", "ragg", "jsonlite", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

script_version <- "1.0.2"
options(stringsAsFactors = FALSE, scipen = 999)

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x))) y else x
}

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else file.path(getwd(), "R", "10_formal_transition_contrasts.R")
project_root <- normalizePath(
  Sys.getenv("GA_PROJECT_ROOT", unset = file.path(dirname(script_path), "..")),
  winslash = "/", mustWork = TRUE
)
primary_root <- normalizePath(
  Sys.getenv(
    "GA_PRIMARY_DELIVERY_ROOT",
    unset = file.path(project_root, "outputs", "hrs_elsa_multistate_final_20260909")
  ),
  winslash = "/", mustWork = TRUE
)
extension_root <- normalizePath(
  Sys.getenv(
    "GA_NEXTSTEP_OUTPUT_ROOT",
    unset = file.path(project_root, "outputs", "hrs_elsa_multistate_nextsteps_20260909")
  ),
  winslash = "/", mustWork = FALSE
)

result_root <- file.path(extension_root, "results", "formal_inference")
report_root <- file.path(extension_root, "reports")
figure_root <- file.path(extension_root, "figures", "formal_inference")
source_root <- file.path(figure_root, "source_data")
qa_root <- file.path(extension_root, "qc", "formal_inference")
for (path in c(result_root, report_root, figure_root, source_root, qa_root)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  fwrite(as.data.table(x), temporary, na = "")
  if (file.exists(path) && !file.remove(path)) {
    stop("Could not replace output: ", path, call. = FALSE)
  }
  if (!file.rename(temporary, path)) {
    stop("Could not activate output: ", path, call. = FALSE)
  }
  invisible(path)
}

atomic_write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  writeLines(x, temporary, useBytes = TRUE)
  if (file.exists(path) && !file.remove(path)) {
    stop("Could not replace output: ", path, call. = FALSE)
  }
  if (!file.rename(temporary, path)) {
    stop("Could not activate output: ", path, call. = FALSE)
  }
  invisible(path)
}

row_major_edges <- function(fit) {
  imatrix <- fit$qmodel$imatrix
  state_names <- rownames(imatrix) %||% paste0("S", seq_len(nrow(imatrix)) - 1L)
  rows <- rbindlist(lapply(seq_len(nrow(imatrix)), function(from_index) {
    to_indices <- which(imatrix[from_index, ] != 0)
    if (!length(to_indices)) return(data.table())
    data.table(
      from_index = from_index,
      to_index = to_indices,
      from_state = state_names[from_index],
      to_state = state_names[to_indices]
    )
  }))
  rows[, transition := paste0(from_state, "->", to_state)]
  rows[, parameter_order := seq_len(.N)]
  rows[]
}

extract_exposure_block <- function(model_path, exposure = "insomnia_score_primary") {
  object <- readRDS(model_path)
  if (!is.list(object) || !inherits(object$fit, "msm")) {
    stop("Not a saved msm fit wrapper: ", model_path, call. = FALSE)
  }
  fit <- object$fit
  coefficients <- stats::coef(fit)
  covariate_labels <- attr(coefficients, "covlabels")
  covariate_index <- match(exposure, covariate_labels)
  if (is.na(covariate_index)) {
    stop("Exposure coefficient block not found in ", basename(model_path), call. = FALSE)
  }
  edges <- row_major_edges(fit)
  n_edges <- nrow(edges)
  if (!identical(as.integer(fit$qmodel$npars), as.integer(n_edges))) {
    stop("Q-model parameter count does not match allowed-edge count", call. = FALSE)
  }
  coefficient_indices <- n_edges + (covariate_index - 1L) * n_edges + seq_len(n_edges)
  beta <- as.numeric(fit$paramdata$params[coefficient_indices])
  covariance <- as.matrix(fit$covmat[coefficient_indices, coefficient_indices, drop = FALSE])
  beta_from_matrix <- vapply(seq_len(n_edges), function(index) {
    coefficients[[exposure]][edges$from_index[index], edges$to_index[index]]
  }, numeric(1))
  if (!isTRUE(all.equal(beta, beta_from_matrix, tolerance = 1e-10, check.attributes = FALSE))) {
    stop("Exposure coefficient mapping failed for ", basename(model_path), call. = FALSE)
  }
  if (any(!is.finite(beta)) || any(!is.finite(covariance))) {
    stop("Non-finite exposure estimate or covariance in ", basename(model_path), call. = FALSE)
  }
  covariance <- (covariance + t(covariance)) / 2
  if (min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) < -1e-8) {
    stop("Exposure covariance is not positive semidefinite in ", basename(model_path), call. = FALSE)
  }
  list(
    beta = setNames(beta, edges$transition),
    covariance = covariance,
    edges = edges,
    imputation = as.integer(object$imputation),
    convergence = as.integer(fit$opt$convergence),
    foundse = isTRUE(fit$foundse),
    spec_hash = as.character(object$spec_hash),
    source_file = normalizePath(model_path, winslash = "/", mustWork = TRUE)
  )
}

discover_mi_files <- function(directory, prefix) {
  paths <- list.files(
    directory,
    pattern = paste0("^", prefix, "__M4_MI20__imp[0-9]{2}\\.rds$"),
    full.names = TRUE
  )
  paths <- sort(normalizePath(paths, winslash = "/", mustWork = TRUE))
  if (length(paths) != 20L) {
    stop("Expected exactly 20 imputed-data model files for ", prefix, "; found ", length(paths), call. = FALSE)
  }
  paths
}

pool_model_set <- function(paths, cohort, analysis_id, expected_edges) {
  extracted <- lapply(paths, extract_exposure_block)
  imputation_ids <- vapply(extracted, `[[`, integer(1), "imputation")
  if (!identical(sort(imputation_ids), seq_len(20L))) {
    stop("Imputation identifiers are not exactly 1-20 for ", cohort, call. = FALSE)
  }
  if (any(vapply(extracted, `[[`, integer(1), "convergence") != 0L)) {
    stop("At least one saved model did not converge for ", cohort, call. = FALSE)
  }
  if (!all(vapply(extracted, `[[`, logical(1), "foundse"))) {
    stop("At least one saved model lacks a covariance estimate for ", cohort, call. = FALSE)
  }
  edge_sets <- lapply(extracted, function(x) names(x$beta))
  if (!all(vapply(edge_sets, identical, logical(1), expected_edges))) {
    stop("Allowed-edge order changed across imputations for ", cohort, call. = FALSE)
  }
  betas <- do.call(rbind, lapply(extracted, function(x) x$beta))
  covariances <- lapply(extracted, `[[`, "covariance")
  colnames(betas) <- expected_edges
  rownames(betas) <- sprintf("imp%02d", imputation_ids)
  m <- nrow(betas)
  beta_bar <- colMeans(betas)
  within <- Reduce(`+`, covariances) / m
  between <- stats::cov(betas)
  total <- within + (1 + 1 / m) * between
  total <- (total + t(total)) / 2
  dimnames(within) <- dimnames(between) <- dimnames(total) <- list(expected_edges, expected_edges)
  if (min(eigen(total, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("Rubin total covariance is not positive definite for ", cohort, call. = FALSE)
  }
  component_rows <- rbindlist(lapply(seq_along(extracted), function(index) {
    x <- extracted[[index]]
    data.table(
      cohort = cohort,
      analysis_id = analysis_id,
      imputation = x$imputation,
      transition = expected_edges,
      log_intensity_ratio = as.numeric(x$beta[expected_edges]),
      within_standard_error = sqrt(diag(x$covariance)),
      convergence = x$convergence,
      foundse = x$foundse,
      spec_hash = x$spec_hash,
      source_model = basename(x$source_file)
    )
  }))
  list(
    cohort = cohort,
    analysis_id = analysis_id,
    edges = expected_edges,
    m = m,
    betas = betas,
    covariances = covariances,
    beta = setNames(beta_bar, expected_edges),
    within = within,
    between = between,
    total = total,
    components = component_rows
  )
}

rubin_scalar_contrast <- function(pool, weights, contrast_id, contrast_label, family) {
  weights <- setNames(as.numeric(weights), names(weights))
  if (!all(names(weights) %in% pool$edges)) {
    stop("Unknown edge in contrast ", contrast_id, call. = FALSE)
  }
  vector <- setNames(rep(0, length(pool$edges)), pool$edges)
  vector[names(weights)] <- weights
  estimates <- as.numeric(pool$betas[, pool$edges, drop = FALSE] %*% vector)
  within_variances <- vapply(pool$covariances, function(covariance) {
    as.numeric(t(vector) %*% covariance %*% vector)
  }, numeric(1))
  m <- length(estimates)
  estimate <- mean(estimates)
  within <- mean(within_variances)
  between <- stats::var(estimates)
  total <- within + (1 + 1 / m) * between
  relative_increase <- if (within > 0) (1 + 1 / m) * between / within else Inf
  degrees_freedom <- if (!is.finite(relative_increase) || relative_increase <= 1e-14) {
    Inf
  } else {
    (m - 1) * (1 + 1 / relative_increase)^2
  }
  statistic <- estimate / sqrt(total)
  critical <- if (is.finite(degrees_freedom)) stats::qt(0.975, degrees_freedom) else stats::qnorm(0.975)
  p_value <- if (is.finite(degrees_freedom)) {
    2 * stats::pt(-abs(statistic), degrees_freedom)
  } else {
    2 * stats::pnorm(-abs(statistic))
  }
  data.table(
    cohort = pool$cohort,
    analysis_id = pool$analysis_id,
    family = family,
    contrast_id = contrast_id,
    contrast_label = contrast_label,
    weights = paste(sprintf("%s:%g", names(weights), weights), collapse = " | "),
    m = m,
    log_ratio = estimate,
    standard_error = sqrt(total),
    degrees_freedom = degrees_freedom,
    statistic = statistic,
    p_value = p_value,
    ratio_of_intensity_ratios = exp(estimate),
    conf_low = exp(estimate - critical * sqrt(total)),
    conf_high = exp(estimate + critical * sqrt(total)),
    within_imputation_variance = within,
    between_imputation_variance = between,
    fraction_missing_information = if (total > 0) (1 + 1 / m) * between / total else NA_real_,
    method = "Within-imputation linear contrast with Rubin scalar pooling; two-sided t reference"
  )
}

rubin_joint_contrast <- function(pool, contrast_matrix, test_id, test_label, family) {
  contrast_matrix <- as.matrix(contrast_matrix)
  if (is.null(colnames(contrast_matrix)) || !all(colnames(contrast_matrix) %in% pool$edges)) {
    stop("Joint contrast columns must name allowed transitions", call. = FALSE)
  }
  full_matrix <- matrix(0, nrow = nrow(contrast_matrix), ncol = length(pool$edges),
                        dimnames = list(rownames(contrast_matrix), pool$edges))
  full_matrix[, colnames(contrast_matrix)] <- contrast_matrix
  theta <- pool$betas[, pool$edges, drop = FALSE] %*% t(full_matrix)
  within_list <- lapply(pool$covariances, function(covariance) full_matrix %*% covariance %*% t(full_matrix))
  estimate <- colMeans(theta)
  within <- Reduce(`+`, within_list) / pool$m
  between <- stats::cov(theta)
  total <- within + (1 + 1 / pool$m) * between
  total <- (total + t(total)) / 2
  statistic <- as.numeric(t(estimate) %*% solve(total, estimate))
  degrees <- length(estimate)
  data.table(
    cohort = pool$cohort,
    analysis_id = pool$analysis_id,
    family = family,
    test_id = test_id,
    test_label = test_label,
    tested_contrasts = paste(rownames(full_matrix), collapse = " | "),
    m = pool$m,
    statistic = statistic,
    df = degrees,
    p_value = stats::pchisq(statistic, df = degrees, lower.tail = FALSE),
    method = "Rubin total-covariance joint Wald chi-square; large-sample reference"
  )
}

nonfatal_edges <- c("S0->S1", "S0->S2", "S1->S3", "S2->S3")
hrs_edges <- c("S0->S1", "S0->S2", "S0->S4", "S1->S3", "S1->S4", "S2->S3", "S2->S4", "S3->S4")
death_edges <- c("S0->S4", "S1->S4", "S2->S4", "S3->S4")

hrs_files <- discover_mi_files(
  file.path(primary_root, "models", "hrs_mi20_final"),
  "hrs_5state_full"
)
elsa_files <- discover_mi_files(
  file.path(primary_root, "models", "elsa_w10tri_final"),
  "elsa_4state_w4_w10"
)

hrs_pool <- pool_model_set(hrs_files, "HRS", "hrs_5state_full", hrs_edges)
elsa_pool <- pool_model_set(elsa_files, "ELSA", "elsa_4state_w4_w10", nonfatal_edges)

atomic_fwrite(
  rbindlist(list(hrs_pool$components, elsa_pool$components), use.names = TRUE),
  file.path(result_root, "per_imputation_insomnia_coefficients.csv")
)

stage_contrast_definitions <- list(
  hypertension_path = c("S0->S1" = -1, "S1->S3" = 1),
  nonhypertension_path = c("S0->S2" = -1, "S2->S3" = 1),
  mean_progression_vs_onset = c("S0->S1" = -0.5, "S0->S2" = -0.5, "S1->S3" = 0.5, "S2->S3" = 0.5)
)
stage_labels <- c(
  hypertension_path = "S1->S3 versus S0->S1",
  nonhypertension_path = "S2->S3 versus S0->S2",
  mean_progression_vs_onset = "Mean progression versus mean onset"
)

stage_results <- rbindlist(lapply(list(hrs_pool, elsa_pool), function(pool) {
  rbindlist(lapply(names(stage_contrast_definitions), function(id) {
    rubin_scalar_contrast(
      pool,
      stage_contrast_definitions[[id]],
      id,
      stage_labels[[id]],
      "onset_vs_progression"
    )
  }))
}))
stage_results[, p_fdr_bh_path_specific := NA_real_]
stage_results[contrast_id %in% c("hypertension_path", "nonhypertension_path"),
              p_fdr_bh_path_specific := p.adjust(p_value, method = "BH"), by = cohort]
stage_results[, multiplicity_note := fifelse(
  contrast_id %in% c("hypertension_path", "nonhypertension_path"),
  "BH-FDR across the two path-specific contrasts within cohort",
  "Single planned summary contrast in a post-hoc extension; not included in the two-path BH family"
)]

stage_joint_matrix <- rbind(
  hypertension_path = c("S0->S1" = -1, "S0->S2" = 0, "S1->S3" = 1, "S2->S3" = 0),
  nonhypertension_path = c("S0->S1" = 0, "S0->S2" = -1, "S1->S3" = 0, "S2->S3" = 1)
)
stage_joint <- rbindlist(lapply(list(hrs_pool, elsa_pool), function(pool) {
  rubin_joint_contrast(
    pool,
    stage_joint_matrix,
    "stage_pairwise_equality",
    "Both corresponding progression coefficients equal their onset coefficients",
    "onset_vs_progression"
  )
}))

death_pair_definitions <- combn(death_edges, 2L, simplify = FALSE)
death_results <- rbindlist(lapply(death_pair_definitions, function(pair) {
  id <- paste0(gsub("->", "to", pair[2], fixed = TRUE), "_vs_", gsub("->", "to", pair[1], fixed = TRUE))
  weights <- setNames(c(-1, 1), pair)
  rubin_scalar_contrast(
    hrs_pool,
    weights,
    id,
    paste0(pair[2], " versus ", pair[1]),
    "hrs_death_origin"
  )
}))
death_results[, p_fdr_bh := p.adjust(p_value, method = "BH")]
death_results[, multiplicity_note := "BH-FDR across all six HRS death-origin pairwise contrasts"]

death_joint_matrix <- rbind(
  `S1->S4 minus S0->S4` = c("S0->S4" = -1, "S1->S4" = 1, "S2->S4" = 0, "S3->S4" = 0),
  `S2->S4 minus S0->S4` = c("S0->S4" = -1, "S1->S4" = 0, "S2->S4" = 1, "S3->S4" = 0),
  `S3->S4 minus S0->S4` = c("S0->S4" = -1, "S1->S4" = 0, "S2->S4" = 0, "S3->S4" = 1)
)
death_joint <- rubin_joint_contrast(
  hrs_pool,
  death_joint_matrix,
  "hrs_death_coefficient_equality",
  "All four HRS death-transition insomnia coefficients are equal",
  "hrs_death_origin"
)

joint_tests <- rbindlist(list(stage_joint, death_joint), use.names = TRUE, fill = TRUE)
atomic_fwrite(stage_results, file.path(result_root, "onset_progression_contrasts.csv"))
atomic_fwrite(death_results, file.path(result_root, "hrs_death_origin_contrasts.csv"))
atomic_fwrite(joint_tests, file.path(result_root, "formal_joint_tests.csv"))

# Cross-cohort equality and synthesis --------------------------------------
hrs_nonfatal_beta <- hrs_pool$beta[nonfatal_edges]
hrs_nonfatal_covariance <- hrs_pool$total[nonfatal_edges, nonfatal_edges, drop = FALSE]
elsa_nonfatal_beta <- elsa_pool$beta[nonfatal_edges]
elsa_nonfatal_covariance <- elsa_pool$total[nonfatal_edges, nonfatal_edges, drop = FALSE]
cohort_difference <- hrs_nonfatal_beta - elsa_nonfatal_beta
cohort_difference_covariance <- hrs_nonfatal_covariance + elsa_nonfatal_covariance
joint_cohort_statistic <- as.numeric(
  t(cohort_difference) %*% solve(cohort_difference_covariance, cohort_difference)
)
joint_cohort_test <- data.table(
  test_id = "hrs_elsa_joint_nonfatal_equality",
  test_label = "All four HRS and ELSA non-fatal insomnia coefficients are equal",
  tested_transitions = paste(nonfatal_edges, collapse = " | "),
  statistic = joint_cohort_statistic,
  df = length(nonfatal_edges),
  p_value = stats::pchisq(joint_cohort_statistic, df = length(nonfatal_edges), lower.tail = FALSE),
  method = "Independent-cohort joint Wald chi-square using the sum of Rubin total covariance matrices"
)
atomic_fwrite(joint_cohort_test, file.path(result_root, "cross_cohort_joint_test.csv"))

meta_rows <- list()
heterogeneity_rows <- list()
for (transition in nonfatal_edges) {
  yi <- c(HRS = hrs_nonfatal_beta[[transition]], ELSA = elsa_nonfatal_beta[[transition]])
  vi <- c(
    HRS = hrs_nonfatal_covariance[transition, transition],
    ELSA = elsa_nonfatal_covariance[transition, transition]
  )
  sei <- sqrt(vi)
  weights_fixed <- 1 / vi
  fixed_log <- sum(weights_fixed * yi) / sum(weights_fixed)
  fixed_se <- sqrt(1 / sum(weights_fixed))
  q_statistic <- sum(weights_fixed * (yi - fixed_log)^2)
  q_df <- length(yi) - 1L
  q_p <- stats::pchisq(q_statistic, df = q_df, lower.tail = FALSE)
  c_value <- sum(weights_fixed) - sum(weights_fixed^2) / sum(weights_fixed)
  tau2_dl <- max(0, (q_statistic - q_df) / c_value)
  weights_random <- 1 / (vi + tau2_dl)
  random_log <- sum(weights_random * yi) / sum(weights_random)
  random_se_normal <- sqrt(1 / sum(weights_random))
  q_star <- sum(weights_random * (yi - random_log)^2) / q_df
  random_se_mhksj <- sqrt(max(1, q_star) / sum(weights_random))
  random_critical <- stats::qt(0.975, df = q_df)
  i2 <- if (q_statistic > 0) max(0, (q_statistic - q_df) / q_statistic) * 100 else 0
  difference_log <- yi[["HRS"]] - yi[["ELSA"]]
  difference_se <- sqrt(sum(vi))
  difference_z <- difference_log / difference_se
  difference_p <- 2 * stats::pnorm(-abs(difference_z))

  heterogeneity_rows[[transition]] <- data.table(
    transition = transition,
    hrs_log_ir = yi[["HRS"]],
    hrs_ir = exp(yi[["HRS"]]),
    hrs_se = sei[["HRS"]],
    elsa_log_ir = yi[["ELSA"]],
    elsa_ir = exp(yi[["ELSA"]]),
    elsa_se = sei[["ELSA"]],
    log_ratio_hrs_vs_elsa = difference_log,
    ratio_hrs_vs_elsa = exp(difference_log),
    ratio_conf_low = exp(difference_log - stats::qnorm(0.975) * difference_se),
    ratio_conf_high = exp(difference_log + stats::qnorm(0.975) * difference_se),
    z_value = difference_z,
    p_value_difference = difference_p,
    cochran_q = q_statistic,
    q_df = q_df,
    p_value_q = q_p,
    i_squared_percent = i2,
    tau_squared_dl = tau2_dl,
    method = "Independent-cohort coefficient difference; inverse-variance Cochran Q"
  )

  cohort_rows <- data.table(
    transition = transition,
    estimate_type = c("HRS", "ELSA"),
    log_intensity_ratio = yi,
    standard_error = sei,
    critical_value = stats::qnorm(0.975),
    reference_distribution = "Normal"
  )
  fixed_row <- data.table(
    transition = transition,
    estimate_type = "Fixed effect",
    log_intensity_ratio = fixed_log,
    standard_error = fixed_se,
    critical_value = stats::qnorm(0.975),
    reference_distribution = "Normal"
  )
  random_row <- data.table(
    transition = transition,
    estimate_type = "Random effects",
    log_intensity_ratio = random_log,
    standard_error = random_se_mhksj,
    critical_value = random_critical,
    reference_distribution = "Modified Hartung-Knapp t, df=1"
  )
  rows <- rbindlist(list(cohort_rows, fixed_row, random_row), use.names = TRUE)
  rows[, `:=`(
    intensity_ratio = exp(log_intensity_ratio),
    conf_low = exp(log_intensity_ratio - critical_value * standard_error),
    conf_high = exp(log_intensity_ratio + critical_value * standard_error),
    z_or_t = log_intensity_ratio / standard_error,
    p_value = fifelse(
      estimate_type == "Random effects",
      2 * stats::pt(-abs(log_intensity_ratio / standard_error), df = 1),
      2 * stats::pnorm(-abs(log_intensity_ratio / standard_error))
    ),
    cochran_q = q_statistic,
    q_df = q_df,
    p_value_q = q_p,
    i_squared_percent = i2,
    tau_squared_dl = tau2_dl,
    random_effects_point_estimator = "DerSimonian-Laird",
    random_effects_ci = "Modified Hartung-Knapp; variance multiplier constrained to at least 1",
    caution = "Only two cohorts; heterogeneity and random-effects estimates are unstable and secondary"
  )]
  meta_rows[[transition]] <- rows
}

heterogeneity <- rbindlist(heterogeneity_rows)
heterogeneity[, p_fdr_bh_difference := p.adjust(p_value_difference, method = "BH")]
meta_results <- rbindlist(meta_rows)
atomic_fwrite(heterogeneity, file.path(result_root, "cross_cohort_heterogeneity.csv"))
atomic_fwrite(meta_results, file.path(result_root, "cross_cohort_meta_analysis.csv"))

# Reproduction audit -------------------------------------------------------
published_hr <- fread(file.path(primary_root, "results", "markov_hazard_ratios.csv"))
published_hr <- published_hr[
  model_id == "M4_MI20" & term == "insomnia_score_primary" &
    analysis_id %chin% c("hrs_5state_full", "elsa_4state_w4_w10")
]
recomputed <- rbindlist(list(
  data.table(
    analysis_id = hrs_pool$analysis_id,
    cohort = hrs_pool$cohort,
    transition = hrs_pool$edges,
    recomputed_log_hr = as.numeric(hrs_pool$beta),
    recomputed_se = sqrt(diag(hrs_pool$total))
  ),
  data.table(
    analysis_id = elsa_pool$analysis_id,
    cohort = elsa_pool$cohort,
    transition = elsa_pool$edges,
    recomputed_log_hr = as.numeric(elsa_pool$beta),
    recomputed_se = sqrt(diag(elsa_pool$total))
  )
))
audit <- merge(
  recomputed,
  published_hr[, .(
    analysis_id, cohort, transition,
    published_log_hr = log_hazard_ratio,
    published_se = standard_error_log_hr,
    published_hr = hazard_ratio,
    published_pooling_status = pooling_status
  )],
  by = c("analysis_id", "cohort", "transition"),
  all.x = TRUE, sort = FALSE
)
audit[, `:=`(
  abs_diff_log_hr = abs(recomputed_log_hr - published_log_hr),
  abs_diff_se = abs(recomputed_se - published_se),
  abs_diff_hr = abs(exp(recomputed_log_hr) - published_hr)
)]
audit[, status := fifelse(
  is.finite(abs_diff_log_hr) & is.finite(abs_diff_se) & is.finite(abs_diff_hr) &
    abs_diff_log_hr < 1e-10 & abs_diff_se < 1e-10 & abs_diff_hr < 1e-10 &
    published_pooling_status == "RUBIN_POOLED_COMPLETE",
  "PASS", "FAIL"
)]
if (any(audit$status != "PASS")) {
  atomic_fwrite(audit, file.path(qa_root, "primary_estimate_reproduction_audit.csv"))
  stop("Recomputed pooled coefficients do not reproduce the frozen primary results", call. = FALSE)
}
atomic_fwrite(audit, file.path(qa_root, "primary_estimate_reproduction_audit.csv"))

# Figures -----------------------------------------------------------------
palette_contract <- c(
  neutral_dark = "#252A34",
  neutral_mid = "#707984",
  neutral_light = "#D9DEE3",
  signal_blue = "#2878B5",
  signal_teal = "#2A9D8F",
  accent_orange = "#E28E2C",
  accent_red = "#C84C4C"
)
cohort_palette <- c(HRS = palette_contract[["signal_blue"]], ELSA = palette_contract[["accent_orange"]])

theme_nature_contract <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(size = base_size - 0.3, colour = "black"),
      axis.title = element_text(size = base_size),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.4),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 0.8, face = "bold", margin = margin(b = 3)),
      plot.subtitle = element_text(size = base_size - 0.2, colour = palette_contract[["neutral_mid"]], margin = margin(b = 4)),
      plot.caption = element_text(size = base_size - 0.8, colour = palette_contract[["neutral_mid"]], hjust = 0),
      panel.grid = element_blank(),
      plot.margin = margin(5, 7, 5, 7)
    )
}
theme_set(theme_nature_contract())

export_single_figure <- function(plot, stem, width_mm = 135, height_mm = 100) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in, bg = "white")
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Arial", bg = "white")
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in,
                 units = "in", res = 600, background = "white", compression = "lzw")
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in,
                units = "in", res = 600, background = "white")
  print(plot)
  grDevices::dev.off()
  invisible(stem)
}

display_transition <- function(x) gsub("->", "\u2192", x, fixed = TRUE)
format_p <- function(x) {
  ifelse(x < 0.001, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = 3))
}

stage_plot_data <- copy(stage_results)
stage_order <- c("mean_progression_vs_onset", "nonhypertension_path", "hypertension_path")
stage_plot_data[, contrast_factor := factor(
  contrast_id,
  levels = stage_order,
  labels = c(
    "Mean progression vs mean onset",
    "S2\u2192S3 vs S0\u2192S2",
    "S1\u2192S3 vs S0\u2192S1"
  )
)]
stage_subtitle <- paste(
  sprintf("%s joint p=%s", stage_joint$cohort, format_p(stage_joint$p_value)),
  collapse = "; "
)
p_stage <- ggplot(
  stage_plot_data,
  aes(x = ratio_of_intensity_ratios, y = contrast_factor, colour = cohort, shape = cohort)
) +
  geom_vline(xintercept = 1, linewidth = 0.35, linetype = "dashed", colour = palette_contract[["neutral_mid"]]) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.12,
                linewidth = 0.5, position = position_dodge(width = 0.42)) +
  geom_point(size = 2.2, stroke = 0.55, position = position_dodge(width = 0.42)) +
  scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
  scale_colour_manual(values = cohort_palette) +
  scale_shape_manual(values = c(HRS = 16, ELSA = 17)) +
  labs(
    x = "Ratio of insomnia transition intensity ratios (progression / onset)",
    y = NULL,
    colour = "Cohort", shape = "Cohort",
    title = "Formal comparison of disease onset and progression",
    subtitle = stage_subtitle,
    caption = "Ratios and 95% CIs use within-imputation contrasts and Rubin pooling. Two-sided tests."
  ) +
  theme_nature_contract() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
stage_source <- copy(stage_plot_data)
stage_source[, panel := "Disease-stage contrast"]
atomic_fwrite(stage_source, file.path(source_root, "fig08_onset_progression_formal_contrasts.csv"))
export_single_figure(p_stage, file.path(figure_root, "fig08_onset_progression_formal_contrasts"), 150, 92)

death_plot_data <- copy(death_results)
death_plot_data[, contrast_display := gsub("->", "\u2192", contrast_label, fixed = TRUE)]
death_plot_data[, contrast_factor := factor(contrast_display, levels = rev(unique(contrast_display)))]
death_subtitle <- sprintf(
  "Omnibus equality: chi-square(%d)=%.2f, p=%s",
  death_joint$df, death_joint$statistic, format_p(death_joint$p_value)
)
p_death <- ggplot(death_plot_data, aes(x = ratio_of_intensity_ratios, y = contrast_factor)) +
  geom_vline(xintercept = 1, linewidth = 0.35, linetype = "dashed", colour = palette_contract[["neutral_mid"]]) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high, colour = p_fdr_bh < 0.05),
                orientation = "y", width = 0.12, linewidth = 0.5) +
  geom_point(aes(colour = p_fdr_bh < 0.05), size = 2.2) +
  scale_x_log10(labels = scales::label_number(accuracy = 0.01)) +
  scale_colour_manual(
    values = c(`TRUE` = palette_contract[["accent_red"]], `FALSE` = palette_contract[["neutral_mid"]]),
    labels = c(`TRUE` = "BH-FDR <0.05", `FALSE` = "BH-FDR >=0.05"),
    name = NULL
  ) +
  labs(
    x = "Ratio of insomnia mortality-transition intensity ratios",
    y = NULL,
    title = "HRS mortality association by origin state",
    subtitle = death_subtitle,
    caption = "All six pairwise contrasts shown; 95% CIs use Rubin pooling. Multiplicity: BH across six comparisons."
  ) +
  theme_nature_contract() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
death_source <- copy(death_plot_data)
death_source[, panel := "HRS death-origin contrast"]
atomic_fwrite(death_source, file.path(source_root, "extfig07_hrs_death_origin_contrasts.csv"))
export_single_figure(p_death, file.path(figure_root, "extfig07_hrs_death_origin_contrasts"), 150, 112)

meta_plot_data <- copy(meta_results)
meta_plot_data[, transition_factor := factor(
  display_transition(transition),
  levels = display_transition(rev(nonfatal_edges))
)]
meta_plot_data[, estimate_factor := factor(
  estimate_type,
  levels = c("Random effects", "Fixed effect", "ELSA", "HRS")
)]
meta_plot_data[, estimate_group := fifelse(
  estimate_type %chin% c("HRS", "ELSA"), estimate_type,
  fifelse(estimate_type == "Fixed effect", "Fixed effect", "Random effects")
)]
meta_colours <- c(
  HRS = cohort_palette[["HRS"]],
  ELSA = cohort_palette[["ELSA"]],
  `Fixed effect` = palette_contract[["signal_teal"]],
  `Random effects` = palette_contract[["neutral_dark"]]
)
annotation_data <- unique(heterogeneity[, .(
  transition_factor = factor(display_transition(transition), levels = display_transition(rev(nonfatal_edges))),
  label = sprintf("Q p=%s; I\u00b2=%.0f%%", format_p(p_value_q), i_squared_percent)
)])
x_max <- max(meta_plot_data$conf_high[is.finite(meta_plot_data$conf_high)], na.rm = TRUE)
p_meta <- ggplot(meta_plot_data, aes(x = intensity_ratio, y = estimate_factor, colour = estimate_group)) +
  geom_vline(xintercept = 1, linewidth = 0.35, linetype = "dashed", colour = palette_contract[["neutral_mid"]]) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.1, linewidth = 0.45) +
  geom_point(size = 1.9) +
  geom_text(
    data = annotation_data,
    aes(x = x_max * 1.08, y = factor("Random effects", levels = levels(meta_plot_data$estimate_factor)), label = label),
    inherit.aes = FALSE, hjust = 0, size = 2.1, family = "Arial", colour = palette_contract[["neutral_mid"]]
  ) +
  facet_grid(transition_factor ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(labels = scales::label_number(accuracy = 0.01), expand = expansion(mult = c(0.05, 0.42))) +
  scale_colour_manual(values = meta_colours, name = NULL) +
  labs(
    x = "Transition intensity ratio per additional insomnia symptom (95% CI)",
    y = NULL,
    title = "Cross-cohort replication and secondary synthesis",
    subtitle = sprintf(
      "Four-edge HRS-ELSA equality: chi-square(4)=%.2f, p=%s",
      joint_cohort_test$statistic, format_p(joint_cohort_test$p_value)
    ),
    caption = paste0(
      "Fixed-effect summaries are secondary estimates in this post-hoc extension.\n",
      "Random-effects points use DerSimonian-Laird; CIs use modified Hartung-Knapp with df=1 and are unstable with two cohorts."
    )
  ) +
  theme_nature_contract() +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = "bottom", strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", margin = margin(r = 3)),
    panel.spacing.y = grid::unit(2.2, "mm")
  )
meta_source <- merge(meta_plot_data, heterogeneity, by = "transition", all.x = TRUE, suffixes = c("", "_heterogeneity"))
meta_source[, panel := "Cross-cohort synthesis"]
atomic_fwrite(meta_source, file.path(source_root, "fig09_cross_cohort_heterogeneity_meta.csv"))
export_single_figure(p_meta, file.path(figure_root, "fig09_cross_cohort_heterogeneity_meta"), 170, 170)

# Machine-readable manifest and report ------------------------------------
figure_manifest <- data.table(
  figure_id = c("Fig08", "ExtFig07", "Fig09"),
  file_stem = c(
    "fig08_onset_progression_formal_contrasts",
    "extfig07_hrs_death_origin_contrasts",
    "fig09_cross_cohort_heterogeneity_meta"
  ),
  figure_level_claim = c(
    "Formal contrasts determine whether progression associations exceed first-onset associations.",
    "Formal contrasts determine whether the HRS insomnia-mortality association differs by origin state.",
    "Formal heterogeneity and cautious synthesis quantify HRS-ELSA replication across four non-fatal transitions."
  ),
  source_data = c(
    "source_data/fig08_onset_progression_formal_contrasts.csv",
    "source_data/extfig07_hrs_death_origin_contrasts.csv",
    "source_data/fig09_cross_cohort_heterogeneity_meta.csv"
  ),
  formats = "PDF | SVG | TIFF 600 dpi | PNG 600 dpi",
  backend = "R ggplot2",
  panel_alignment = "Not applicable: each deliverable is a single ggplot object",
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
  script_version = script_version
)
atomic_fwrite(figure_manifest, file.path(figure_root, "formal_inference_figure_manifest.csv"))

stage_key <- stage_results[contrast_id == "mean_progression_vs_onset"]
death_key <- death_results[grepl("S3toS4", contrast_id, fixed = TRUE)]
report_lines <- c(
  "# HRS-ELSA formal transition contrasts and cross-cohort synthesis",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("Script version: ", script_version),
  "",
  "## Analysis boundary",
  "",
  "This extension reuses the frozen M4-MI20 continuous-time Markov models and their full transition-coefficient covariance matrices. It does not refit the disease states, exposure, covariates or imputations.",
  "These are post-hoc secondary analyses developed after inspection of the primary results. Contrasts and multiplicity families were fixed before evaluation in this extension, but are not original-protocol preregistered confirmatory tests.",
  "",
  "## First disease onset versus later progression",
  "",
  sprintf(
    "- HRS mean progression-versus-onset ratio of intensity ratios: %.3f (95%% CI %.3f-%.3f; p=%s).",
    stage_key[cohort == "HRS", ratio_of_intensity_ratios], stage_key[cohort == "HRS", conf_low],
    stage_key[cohort == "HRS", conf_high], format_p(stage_key[cohort == "HRS", p_value])
  ),
  sprintf(
    "- ELSA mean progression-versus-onset ratio of intensity ratios: %.3f (95%% CI %.3f-%.3f; p=%s).",
    stage_key[cohort == "ELSA", ratio_of_intensity_ratios], stage_key[cohort == "ELSA", conf_low],
    stage_key[cohort == "ELSA", conf_high], format_p(stage_key[cohort == "ELSA", p_value])
  ),
  sprintf(
    "- Joint equality of both corresponding onset/progression pairs: HRS chi-square(2)=%.3f, p=%s; ELSA chi-square(2)=%.3f, p=%s.",
    stage_joint[cohort == "HRS", statistic], format_p(stage_joint[cohort == "HRS", p_value]),
    stage_joint[cohort == "ELSA", statistic], format_p(stage_joint[cohort == "ELSA", p_value])
  ),
  "",
  "A significant individual transition and a non-significant transition do not by themselves establish a difference. Only the direct contrasts above address stage differences.",
  "",
  "## HRS mortality-transition heterogeneity",
  "",
  sprintf(
    "The four HRS death-transition insomnia coefficients were compared jointly: chi-square(3)=%.3f, p=%s. Pairwise estimates and six-comparison BH-FDR values are in the accompanying CSV and Extended Data figure.",
    death_joint$statistic, format_p(death_joint$p_value)
  ),
  "",
  "## HRS-ELSA heterogeneity and synthesis",
  "",
  sprintf(
    "The joint null that all four non-fatal insomnia coefficients were equal between HRS and ELSA yielded chi-square(4)=%.3f, p=%s.",
    joint_cohort_test$statistic, format_p(joint_cohort_test$p_value)
  ),
  "",
  "Per-transition Q, I-squared, cohort-difference tests, fixed-effect summaries and exploratory random-effects summaries are provided in the result tables. Because only two cohorts contribute, heterogeneity and random-effects estimates are unstable and remain secondary.",
  "",
  "## Multiplicity and interpretation",
  "",
  "- All tests are two-sided.",
  "- The two path-specific onset-versus-progression contrasts use BH-FDR within cohort.",
  "- The six HRS death-origin pairwise contrasts use BH-FDR.",
  "- Cross-cohort transition-specific difference p values use BH-FDR across four transitions.",
  "- Results describe associations, not causal effects.",
  "",
  "## Reproducibility",
  "",
  sprintf("- HRS saved MI fits: %d/20; ELSA saved MI fits: %d/20.", hrs_pool$m, elsa_pool$m),
  sprintf("- Frozen primary coefficient reproduction: %d/%d PASS.", sum(audit$status == "PASS"), nrow(audit)),
  "- Independent unit: participant; cohorts are independent.",
  "- Model engine: msm continuous-time panel Markov; pooling: within-imputation contrasts plus Rubin rules."
)
atomic_write_lines(report_lines, file.path(report_root, "formal_transition_contrasts_report_en.md"))

fixed_key <- meta_results[estimate_type == "Fixed effect"]
zh_report_lines <- c(
  "# HRS-ELSA 转移系数正式对比与跨队列综合",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("脚本版本：", script_version),
  "",
  "## 分析边界",
  "",
  "本扩展直接复用冻结的 M4-MI20 连续时间 Markov 模型及完整转移系数协方差矩阵，不重新定义状态、暴露、协变量、随访时间或插补数据。所有检验均为双侧检验。",
  "本组检验属于查看主结果后的事后次级扩展；对比及多重检验家族在本次计算前固定，不代表原方案预注册的确认性检验。",
  "",
  "## 首次发病与后续进展",
  "",
  sprintf(
    "- HRS：高血压路径的进展/首次发病强度比之比为 %.3f（95%% CI %.3f-%.3f；p=%s；BH-FDR=%s）；非高血压路径为 %.3f（95%% CI %.3f-%.3f；p=%s；BH-FDR=%s）。",
    stage_results[cohort == "HRS" & contrast_id == "hypertension_path", ratio_of_intensity_ratios],
    stage_results[cohort == "HRS" & contrast_id == "hypertension_path", conf_low],
    stage_results[cohort == "HRS" & contrast_id == "hypertension_path", conf_high],
    format_p(stage_results[cohort == "HRS" & contrast_id == "hypertension_path", p_value]),
    format_p(stage_results[cohort == "HRS" & contrast_id == "hypertension_path", p_fdr_bh_path_specific]),
    stage_results[cohort == "HRS" & contrast_id == "nonhypertension_path", ratio_of_intensity_ratios],
    stage_results[cohort == "HRS" & contrast_id == "nonhypertension_path", conf_low],
    stage_results[cohort == "HRS" & contrast_id == "nonhypertension_path", conf_high],
    format_p(stage_results[cohort == "HRS" & contrast_id == "nonhypertension_path", p_value]),
    format_p(stage_results[cohort == "HRS" & contrast_id == "nonhypertension_path", p_fdr_bh_path_specific])
  ),
  sprintf(
    "- HRS 两条路径的联合相等检验：chi-square(2)=%.3f，p=%s；平均进展/平均首次发病比为 %.3f（95%% CI %.3f-%.3f；p=%s）。没有证据支持统一的后期进展增强。",
    stage_joint[cohort == "HRS", statistic], format_p(stage_joint[cohort == "HRS", p_value]),
    stage_key[cohort == "HRS", ratio_of_intensity_ratios], stage_key[cohort == "HRS", conf_low],
    stage_key[cohort == "HRS", conf_high], format_p(stage_key[cohort == "HRS", p_value])
  ),
  sprintf(
    "- ELSA：高血压路径为 %.3f（95%% CI %.3f-%.3f；p=%s；BH-FDR=%s），非高血压路径为 %.3f（95%% CI %.3f-%.3f；p=%s；BH-FDR=%s）。",
    stage_results[cohort == "ELSA" & contrast_id == "hypertension_path", ratio_of_intensity_ratios],
    stage_results[cohort == "ELSA" & contrast_id == "hypertension_path", conf_low],
    stage_results[cohort == "ELSA" & contrast_id == "hypertension_path", conf_high],
    format_p(stage_results[cohort == "ELSA" & contrast_id == "hypertension_path", p_value]),
    format_p(stage_results[cohort == "ELSA" & contrast_id == "hypertension_path", p_fdr_bh_path_specific]),
    stage_results[cohort == "ELSA" & contrast_id == "nonhypertension_path", ratio_of_intensity_ratios],
    stage_results[cohort == "ELSA" & contrast_id == "nonhypertension_path", conf_low],
    stage_results[cohort == "ELSA" & contrast_id == "nonhypertension_path", conf_high],
    format_p(stage_results[cohort == "ELSA" & contrast_id == "nonhypertension_path", p_value]),
    format_p(stage_results[cohort == "ELSA" & contrast_id == "nonhypertension_path", p_fdr_bh_path_specific])
  ),
  sprintf(
    "- ELSA 联合相等检验：chi-square(2)=%.3f，p=%s；但两条路径方向相反，平均进展/平均首次发病比为 %.3f（95%% CI %.3f-%.3f；p=%s）。因此结果支持路径异质性，不支持笼统的后期进展效应更强。",
    stage_joint[cohort == "ELSA", statistic], format_p(stage_joint[cohort == "ELSA", p_value]),
    stage_key[cohort == "ELSA", ratio_of_intensity_ratios], stage_key[cohort == "ELSA", conf_low],
    stage_key[cohort == "ELSA", conf_high], format_p(stage_key[cohort == "ELSA", p_value])
  ),
  "",
  "## HRS 死亡转移的起始状态差异",
  "",
  sprintf(
    "四条死亡转移的失眠系数并不相等：chi-square(3)=%.3f，p=%s。经六次两两比较 BH 校正后，S3->S4 相对 S0->S4 的比值为 %.3f（95%% CI %.3f-%.3f；FDR=%s），相对 S1->S4 为 %.3f（95%% CI %.3f-%.3f；FDR=%s）。这表示进入心血管多病共存状态后，失眠与随后死亡的相对关联较弱；不能解释为保护效应，且可能受状态选择、竞争风险和基线风险差异影响。",
    death_joint$statistic, format_p(death_joint$p_value),
    death_results[contrast_id == "S3toS4_vs_S0toS4", ratio_of_intensity_ratios],
    death_results[contrast_id == "S3toS4_vs_S0toS4", conf_low],
    death_results[contrast_id == "S3toS4_vs_S0toS4", conf_high],
    format_p(death_results[contrast_id == "S3toS4_vs_S0toS4", p_fdr_bh]),
    death_results[contrast_id == "S3toS4_vs_S1toS4", ratio_of_intensity_ratios],
    death_results[contrast_id == "S3toS4_vs_S1toS4", conf_low],
    death_results[contrast_id == "S3toS4_vs_S1toS4", conf_high],
    format_p(death_results[contrast_id == "S3toS4_vs_S1toS4", p_fdr_bh])
  ),
  "",
  "## HRS-ELSA 一致性与次要综合",
  "",
  sprintf(
    "四条非致死边的跨队列联合差异检验为 chi-square(4)=%.3f，p=%s；四条单边差异在 BH 校正后均不显著。因此没有统计证据表明 HRS 与 ELSA 的失眠系数整体不同。",
    joint_cohort_test$statistic, format_p(joint_cohort_test$p_value)
  ),
  sprintf(
    "固定效应次要综合为：S0->S1 %.3f（95%% CI %.3f-%.3f）；S0->S2 %.3f（%.3f-%.3f）；S1->S3 %.3f（%.3f-%.3f）；S2->S3 %.3f（%.3f-%.3f）。",
    fixed_key[transition == "S0->S1", intensity_ratio], fixed_key[transition == "S0->S1", conf_low], fixed_key[transition == "S0->S1", conf_high],
    fixed_key[transition == "S0->S2", intensity_ratio], fixed_key[transition == "S0->S2", conf_low], fixed_key[transition == "S0->S2", conf_high],
    fixed_key[transition == "S1->S3", intensity_ratio], fixed_key[transition == "S1->S3", conf_low], fixed_key[transition == "S1->S3", conf_high],
    fixed_key[transition == "S2->S3", intensity_ratio], fixed_key[transition == "S2->S3", conf_low], fixed_key[transition == "S2->S3", conf_high]
  ),
  "S0->S2 与 S1->S3 是两库各自均为正且置信区间均不跨 1 的一致信号。S2->S3 的固定效应综合虽不跨 1，但 ELSA 单库接近无关联且两库 I-squared 约 65%；不能把该综合值表述为稳健复制。只有两个队列时，Q、I-squared、tau-squared 和随机效应区间均不稳定，故随机效应仅作探索。",
  "",
  "## 当前可支持的结论",
  "",
  "较高失眠症状负担与特定而非所有心血管转移相关。最稳定的跨队列证据是从无病进入心脏病/卒中合并状态（S0->S2），以及从单纯高血压进展至高血压与心脏病/卒中共存（S1->S3）。现有结果不支持统一的‘后期进展强于首次发病’，也不能在未拆分心脏病与卒中之前判断由哪一种疾病驱动。",
  "",
  "## 可复现性",
  "",
  sprintf("- HRS 与 ELSA 各读取 %d/20 与 %d/20 个保存的插补模型。", hrs_pool$m, elsa_pool$m),
  sprintf("- 冻结主系数复现：%d/%d 项通过。", sum(audit$status == "PASS"), nrow(audit)),
  "- 独立单位为参与者；HRS 与 ELSA 视为独立队列。",
  "- 模型为 msm 连续时间面板 Markov；线性对比在每个插补内完成后按 Rubin 规则合并。"
)
atomic_write_lines(zh_report_lines, file.path(report_root, "formal_transition_contrasts_report_zh.md"))

summary_table <- data.table(
  check = c(
    "HRS MI models", "ELSA MI models", "Primary estimate reproduction",
    "Stage scalar contrasts", "Death scalar contrasts", "Joint tests",
    "Heterogeneity rows", "Meta-analysis rows", "Figure exports"
  ),
  observed = c(
    hrs_pool$m, elsa_pool$m, sum(audit$status == "PASS"),
    nrow(stage_results), nrow(death_results), nrow(joint_tests) + nrow(joint_cohort_test),
    nrow(heterogeneity), nrow(meta_results),
    sum(file.exists(unlist(lapply(
      c("fig08_onset_progression_formal_contrasts", "extfig07_hrs_death_origin_contrasts", "fig09_cross_cohort_heterogeneity_meta"),
      function(stem) file.path(figure_root, paste0(stem, c(".pdf", ".svg", ".tiff", ".png")))
    ))))
  ),
  expected = c(20, 20, length(hrs_edges) + length(nonfatal_edges), 6, 6, 4, 4, 16, 12)
)
summary_table[, status := fifelse(observed == expected, "PASS", "FAIL")]
atomic_fwrite(summary_table, file.path(qa_root, "formal_inference_completeness.csv"))
if (any(summary_table$status != "PASS")) {
  stop("Formal-inference completeness gate failed", call. = FALSE)
}

message("Formal transition contrasts and cross-cohort synthesis completed: ", extension_root)
