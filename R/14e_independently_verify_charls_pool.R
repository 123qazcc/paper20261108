#!/usr/bin/env Rscript
# Independent computational audit, not a new fit and not an assessment of all
# epidemiological assumptions. Never writes into models/ or results/.
# Run from the project root; --self-test uses synthetic data only.
options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.treatment"))
suppressPackageStartupMessages({library(data.table); library(msm); library(digest)})
setDTthreads(1L)

same_number <- function(a, b, tolerance = 1e-8) {
  if (!identical(dim(a), dim(b)) || length(a) != length(b)) return(FALSE)
  a <- as.numeric(a); b <- as.numeric(b)
  if (!identical(is.na(a), is.na(b))) return(FALSE)
  keep <- !is.na(a)
  a <- a[keep]; b <- b[keep]
  if (!identical(is.infinite(a), is.infinite(b))) return(FALSE)
  infinite <- is.infinite(a)
  if (any(a[infinite] != b[infinite])) return(FALSE)
  a <- a[!infinite]; b <- b[!infinite]
  all(abs(a - b) <= tolerance * pmax(1, abs(a), abs(b)))
}
require_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}
equal_groups <- function(a, b) {
  length(a) == length(b) && identical(outer(a, a, `==`), outer(b, b, `==`))
}

# Recover covariance from the saved *free-parameter* Hessian. msm minimizes
# minus twice log likelihood, hence covariance = 2 * inverse(H). The binary
# expansion matrix maps free parameters to every constrained Q/covariate cell.
# The current contract has positive equality constraints and no fixed effects;
# reject other schemas rather than silently assume they are equivalent.
extract_from_full_fit <- function(fit, edges, exposure = "restless_bl") {
  require_true(inherits(fit, "msm"), "Full msm fit is unavailable")
  p <- fit$paramdata
  require_true(length(p$fixedpars) == 0L && length(p$auxpars) == 0L &&
                 all(p$constr > 0), "Unsupported fixed, auxiliary or negative constraint schema")
  require_true(identical(as.integer(fit$opt$convergence), 0L) && isTRUE(fit$foundse),
               "Optimizer convergence / standard-error gate failed")
  require_true(is.finite(as.numeric(logLik(fit))), "Non-finite log likelihood")
  h <- as.matrix(fit$opt$hessian)
  require_true(all(is.finite(h)) && same_number(h, t(h), 1e-10),
               "Free Hessian missing, non-finite or asymmetric")
  h <- (h + t(h)) / 2
  eigenvalues <- eigen(h, symmetric = TRUE, only.values = TRUE)$values
  require_true(min(eigenvalues) > 1e-8, "Free Hessian does not pass positive-definiteness gate")
  require_true(nrow(h) == length(p$optpars) && nrow(h) == length(fit$opt$par),
               "Free-parameter Hessian dimensions disagree")
  group_index <- match(p$constr, p$constr[p$optpars])
  require_true(!anyNA(group_index) && !anyDuplicated(p$constr[p$optpars]),
               "Cannot map all constraints to unique optimized parameters")
  expand <- matrix(0, length(p$constr), length(p$optpars))
  expand[cbind(seq_along(group_index), group_index)] <- 1
  expanded_theta <- drop(expand %*% as.numeric(fit$opt$par))
  expanded_covariance <- expand %*% (2 * chol2inv(chol(h))) %*% t(expand)
  require_true(same_number(expanded_theta, as.numeric(p$params), 1e-9),
               "Free-to-expanded parameter reconstruction disagrees")
  require_true(same_number(expanded_covariance, fit$covmat, 2e-7),
               "Hessian-based covariance reconstruction disagrees with full fit")

  # msm enumerates Q cells by row (not R's default column-major order).
  q <- fit$qmodel$imatrix
  cells <- which(q != 0, arr.ind = TRUE)
  cells <- cells[order(cells[, 1L], cells[, 2L]), , drop = FALSE]
  require_true(nrow(cells) == nrow(edges) &&
                 identical(as.integer(cells[, 1L]), as.integer(edges$from_index)) &&
                 identical(as.integer(cells[, 2L]), as.integer(edges$to_index)),
               "Full-fit Q edge ordering differs from the frozen bundle")
  labels <- fit$qcmodel$covlabels
  exposure_column <- match(exposure, labels)
  require_true(!is.na(exposure_column), "Exposure absent from the full-fit covariates")
  # Identify covariate-class parameters first, rather than assuming a fixed
  # absolute offset of 20, and cross-check the coefficient matrix independently.
  qcov_indices <- which(p$plabs == "qcov")
  require_true(length(qcov_indices) == nrow(edges) * length(labels),
               "Unexpected Q-covariate parameter layout")
  idx <- matrix(qcov_indices, nrow = nrow(edges))[, exposure_column]
  theta <- as.numeric(expanded_theta[idx])
  V <- expanded_covariance[idx, idx, drop = FALSE]
  coefficient_matrix <- coef(fit)[[exposure]]
  matrix_theta <- as.numeric(coefficient_matrix[cells])
  require_true(same_number(theta, matrix_theta, 1e-10),
               "Reconstructed exposure differs from coef(fit)")
  require_true(equal_groups(p$constr[idx], edges$exposure_group),
               "Exposure equality constraints do not match the frozen specification")
  for (j in setdiff(seq_along(labels), exposure_column)) {
    nuisance_idx <- matrix(qcov_indices, nrow = nrow(edges))[, j]
    require_true(length(unique(p$constr[nuisance_idx])) == 1L,
                 paste("Nuisance effect is not shared across edges:", labels[j]))
  }
  v_eigen <- eigen((V + t(V)) / 2, symmetric = TRUE, only.values = TRUE)$values
  require_true(all(is.finite(V)) && all(diag(V) >= 0) &&
                 min(v_eigen) >= -1e-8 * max(1, abs(v_eigen)),
               "Expanded exposure covariance failed the PSD gate")
  # The full 20-edge covariance is expected to be singular due to equality
  # constraints; positive definiteness is required only for tested contrasts.
  list(theta = theta, covariance = (V + t(V)) / 2,
       hessian_min = min(eigenvalues), hessian_condition = max(eigenvalues) / min(eigenvalues),
       likelihood = as.numeric(logLik(fit)), expanded_indices = idx)
}

independent_pool <- function(theta, covariance) {
  m <- ncol(theta); k <- nrow(theta)
  require_true(m > 1L && length(covariance) == m && all(is.finite(theta)),
               "At least two complete finite imputations are required")
  require_true(all(vapply(covariance, function(v)
    identical(dim(v), c(k, k)) && all(is.finite(v)), TRUE)), "Invalid within-imputation covariance")
  mean_theta <- rowSums(theta) / m
  centered <- sweep(theta, 1L, mean_theta, `-`)
  between <- tcrossprod(centered) / (m - 1)
  within <- matrix(0, k, k)
  for (i in seq_len(m)) within <- within + covariance[[i]] / m
  total <- within + ((m + 1) / m) * between
  list(m = m, Q = theta, qbar = mean_theta, U = within, B = between, T = total)
}

scalar_statistics <- function(pool, L) {
  estimates <- as.numeric(L %*% pool$qbar)
  within <- diag(L %*% pool$U %*% t(L))
  between <- diag(L %*% pool$B %*% t(L))
  total <- diag(L %*% pool$T %*% t(L))
  require_true(all(within > 0) && all(total > 0) && all(between >= -1e-12),
               "Non-positive variance for a reported scalar effect")
  # Contract uses Rubin's large-complete-sample df, not Barnard-Rubin df.
  # Match the explicitly documented 1e-15 numerical-zero convention.
  df <- rep(Inf, length(estimates))
  nonzero <- between > 1e-15
  relative_increase <- ((pool$m + 1) / pool$m) * between[nonzero] / within[nonzero]
  df[nonzero] <- (pool$m - 1) * ((relative_increase + 1) / relative_increase)^2
  standard_error <- sqrt(total)
  half_width <- qt(0.975, df = df) * standard_error
  data.table(log_effect = estimates, se = standard_error, df = df,
             HR = exp(estimates), lower = exp(estimates - half_width),
             upper = exp(estimates + half_width),
             p = 2 * pt(abs(estimates / standard_error), df = df, lower.tail = FALSE))
}

joint_statistic <- function(pool, L) {
  V <- L %*% pool$T %*% t(L)
  V <- (V + t(V)) / 2
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  require_true(all(is.finite(ev)) && min(ev) > 1e-12 * max(ev) && qr(V)$rank == nrow(L),
               "Joint contrast covariance is singular or poorly identified")
  # Triangular whitening is independent of the production solve(V, a) route.
  standardized <- forwardsolve(t(chol(V)), as.numeric(L %*% pool$qbar))
  statistic <- sum(standardized^2)
  data.table(chi_square = statistic, df = nrow(L),
             p = pchisq(statistic, df = nrow(L), lower.tail = FALSE))
}

self_test <- function() {
  m <- 20L; offsets <- seq(-0.19, 0.19, length.out = m)
  location <- c(0.2, 0.1, -0.05, 0.15); loading <- c(1, 2, -1, 0.5)
  theta <- matrix(location, 4L, m) + outer(loading, offsets)
  U <- diag(c(0.04, 0.09, 0.06, 0.08)); U[1, 2] <- U[2, 1] <- 0.015
  U[3, 4] <- U[4, 3] <- -0.008
  pool <- independent_pool(theta, rep(list(U), m))
  stopifnot(same_number(pool$qbar, location), same_number(pool$U, U),
            same_number(pool$B, var(offsets) * tcrossprod(loading)))
  C <- rbind(c(1, -1, 0, 0), c(0, 0, 1, -1))
  s <- scalar_statistics(pool, C)
  stopifnot(same_number(s$log_effect, c(0.1, -0.2)),
            same_number(s$se[1]^2, 0.04 + 0.09 - 2 * 0.015 + 1.05 * var(offsets)),
            all(s$df > 19), all(s$lower < s$HR), all(s$upper > s$HR))
  joint <- joint_statistic(pool, C)
  a <- as.numeric(C %*% pool$qbar); v <- C %*% pool$T %*% t(C)
  stopifnot(same_number(joint$chi_square, as.numeric(crossprod(a, solve(v, a)))))
  constant <- independent_pool(matrix(location, 4L, m), rep(list(U), m))
  stopifnot(all(is.infinite(scalar_statistics(constant, diag(4))$df)))
  stopifnot(inherits(try(joint_statistic(pool, rbind(C[1, ], C[1, ])), silent = TRUE), "try-error"))

  # A small constrained synthetic msm fit tests Hessian scaling, Q ordering,
  # nuisance constraints, and covariate expansion; never enters study results.
  set.seed(149013)
  n <- 180L
  z <- data.frame(id = sprintf("toy%03d", rep(seq_len(n), each = 3L)),
                  subject_id = rep(seq_len(n), each = 3L), time = rep(c(0, 1, 3), n),
                  state = 1L, restless_bl = rep(rep(0:3, length.out = n), each = 3L),
                  auxiliary = rep(rnorm(n), each = 3L))
  for (i in seq_len(n)) {
    rows <- (3L * i - 2L):(3L * i)
    z$state[rows[2L]] <- sample(c(1L, 2L, 3L), 1L, prob = c(0.62, 0.30, 0.08))
    pool_states <- seq.int(z$state[rows[2L]], 3L)
    z$state[rows[3L]] <- pool_states[sample.int(length(pool_states), 1L)]
  }
  z$state[c(2, 5, 8)] <- 99L
  q <- matrix(c(0, .15, .03, 0, 0, .15, 0, 0, 0), 3L, 3L, byrow = TRUE)
  fit <- suppressWarnings(msm(state ~ time, subject = subject_id, data = z, qmatrix = q,
    covariates = ~restless_bl + auxiliary,
    constraint = list(restless_bl = c(1L, 2L, 1L), auxiliary = c(1L, 1L, 1L)),
    obstype = rep(1L, nrow(z)), censor = 99L, censor.states = list(c(1L, 2L)),
    hessian = TRUE, control = list(maxit = 5000, reltol = 1e-9, fnscale = 1)))
  e <- data.frame(from_index = c(1L, 1L, 2L), to_index = c(2L, 3L, 3L),
                   exposure_group = c(1L, 2L, 1L))
  independent <- extract_from_full_fit(fit, e)
  stopifnot(same_number(independent$theta[1], independent$theta[3]),
            same_number(independent$covariance[1, ], independent$covariance[3, ]))
  message("SELF_TEST_PASS: analytic Rubin/contrasts/zero-between-variance/singular-joint rejection; ",
          "synthetic constrained msm Hessian and full covariance reconstruction. No study estimates used.")
}

run_review <- function() {
  out <- "outputs/charls_external_extension_20260913"
  qc <- file.path(out, "model_qc")
  require_true(dir.exists(qc), "CHARLS QC directory not found; run from the project root")
  checks <- list(); component_diagnostics <- list()
  add <- function(check, pass, detail = "", severity = "P0") {
    checks[[length(checks) + 1L]] <<- data.table(check = check, pass = isTRUE(pass),
      severity = severity, detail = as.character(detail))
    isTRUE(pass)
  }
  guard <- function(check, pass, detail = "") {
    if (!add(check, pass, detail)) stop(paste(check, detail, sep = ": "), call. = FALSE)
  }
  publish <- function(status, error = "") {
    results <- rbindlist(checks, fill = TRUE)
    results[, review_status := status]
    fwrite(results, file.path(qc, "independent_pool_review.csv"))
    if (length(component_diagnostics))
      fwrite(rbindlist(component_diagnostics, fill = TRUE), file.path(qc, "independent_component_diagnostics.csv"))
    failures <- if (nrow(results)) results[pass == FALSE, paste0("- ", check, ": ", detail)] else "- No checks available."
    writeLines(c("# CHARLS MI20 独立计算核验", "", paste0("状态：", status),
      paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), "",
      paste0("已通过：", sum(results$pass), "/", nrow(results), " 项。"), "",
      "## 核验范围", "",
      "从每份完整 msm 模型的自由参数和 Hessian 重建所有约束后的暴露系数与完整协方差，并与 coef(fit)、保存的 block 和已发布结果逐项对照。",
      "以参与者为独立单位核对同一组 ID、完整观测序列、插补协变量及模型记录数；20 份插补不是 20 个独立人群。",
      "独立重算 Rubin 的 U、B、T、有限插补 t 区间与双侧 p、BH(4)/BH(2)、心脏病减卒中对比和渐近卡方联合 Wald。", "",
      "## 边界", "",
      "计算核验通过不等于统计或流行病学假设已成立；本脚本不验证 Markov/时间齐性、无未测混杂、失访可忽略性、死亡完整性或 MICE 与多状态似然的同源性。",
      "Hessian 是已保存拟合的数值导数；本核验独立重建和核对协方差，但不重新拟合或独立数值微分原始似然。",
      "联合检验是使用 Rubin 总协方差的渐近卡方 Wald，不是有限插补 D1 检验；标量自由度采用 large-complete-sample Rubin 公式。",
      "心脏病/卒中差异须依据直接对比，不能用一项显著、另一项不显著来声称病种差异；结果仅支持关联而非因果。", "",
      "## 未通过或尚未就绪", "", if (length(failures)) failures else "无计算核验失败项。",
      if (nzchar(error)) c("", paste0("运行说明：", error)) else character()),
      file.path(qc, "independent_pool_review.md"), useBytes = TRUE)
    status
  }
  bundle_path <- file.path(out, "data/charls_model_bundle.rds")
  model_paths <- file.path(out, "models", sprintf("MI_%02d.rds", 1:20))
  result_paths <- file.path(out, "results", c("all_edge_MI20.csv", "four_target_MI20.csv",
    "heart_stroke_contrasts_MI20.csv", "joint_tests_MI20.csv", "rubin_matrices_MI20.rds"))
  prerequisite_paths <- c(bundle_path, model_paths, result_paths,
    file.path(qc, "model_bundle_review.csv"), "R/14b_fit_charls_disease_models.R")
  missing <- prerequisite_paths[!file.exists(prerequisite_paths)]
  if (length(missing)) {
    add("all_20_models_and_pooled_outputs_ready", FALSE, paste(missing, collapse = " | "))
    publish("NOT_READY", "Missing required files; no partial-imputation result was accepted.")
    message("NOT_READY: ", length(missing), " required files absent. No passing review generated.")
    return(2L)
  }
  status <- tryCatch({
    sha <- function(path) digest(file = path, algo = "sha256")
    bundle <- readRDS(bundle_path)
    guard("bundle_status_and_checks", identical(bundle$status, "APPROVED_FOR_PILOT_AND_MI") &&
      all(bundle$preparation_checks$pass) && all(fread(file.path(qc, "model_bundle_review.csv"))$pass))
    guard("bundle_has_20_imputations", length(bundle$imputations) == 20L)
    guard("all_signed_input_files_exist", all(file.exists(names(bundle$input_hashes))))
    guard("all_signed_input_hashes_current", identical(unname(vapply(names(bundle$input_hashes), sha, "")),
                                                       unname(bundle$input_hashes)))
    signature <- digest(list(bundle = sha(bundle_path), script = sha("R/14b_fit_charls_disease_models.R"),
                              msm = as.character(packageVersion("msm"))), algo = "sha256")
    edges <- as.data.table(bundle$edges)
    expected_edges <- rbindlist(lapply(0:7, function(from) {
      bits <- c(1L, 2L, 4L)
      to <- c(vapply(bits[bitwAnd(from, bits) == 0L], function(bit) bitwOr(from, bit), 1L), 8L)
      data.table(from_index = from + 1L, to_index = to + 1L)
    }))
    setorder(expected_edges, from_index, to_index)
    guard("9_state_20_edge_contract", nrow(edges) == 20L &&
      identical(as.integer(edges$from_index), expected_edges$from_index) &&
      identical(as.integer(edges$to_index), expected_edges$to_index) &&
      identical(as.character(edges$transition), paste(edges$from_index - 1L, edges$to_index - 1L, sep = "-")))
    targets <- match(c("0-2", "0-4", "1-3", "1-5"), edges$transition)
    guard("four_prespecified_targets_present", !anyNA(targets) && !anyDuplicated(targets))
    expected_groups <- ifelse(edges$to_index == 9L, 8L,
      match(edges$to_index - edges$from_index, c(1L, 2L, 4L)) + 4L)
    expected_groups[targets] <- 1:4
    guard("prespecified_exposure_group_structure", !anyNA(expected_groups) &&
      identical(as.integer(edges$exposure_group), as.integer(expected_groups)))
    panel <- as.data.table(copy(bundle$panel)); setorder(panel, id, time)
    ids <- unique(as.character(panel$id))
    expected_subject <- match(as.character(panel$id), ids)
    guard("participant_ids_unique_baseline", !anyDuplicated(bundle$baseline$id) &&
      setequal(ids, as.character(bundle$baseline$id)))
    guard("all_panel_observations_type1", all(panel$obstype == 1L))
    extracted <- vector("list", 20L)
    for (i in 1:20) {
      label <- sprintf("MI_%02d", i)
      saved <- readRDS(model_paths[i])
      guard(paste0(label, "_identity_signature"), identical(saved$status, "ACCEPTED") &&
        identical(as.integer(saved$imputation), i) && identical(saved$model, label) &&
        identical(saved$signature, signature))
      imp <- as.data.table(bundle$imputations[[i]])
      guard(paste0(label, "_imputation_IDs"), !anyDuplicated(imp$id) && setequal(ids, as.character(imp$id)))
      original <- saved$fit$call$data
      guard(paste0(label, "_saved_original_data_available"), is.data.frame(original) &&
        all(c("id", "subject_id", "time", "state", "obstype") %in% names(original)),
        "The contract's do.call(msm, data=z) must preserve its actual dataframe; symbolic calls are not evaluated.")
      guard(paste0(label, "_sample_IDs_and_observations"),
        identical(as.character(original$id), as.character(panel$id)) &&
        identical(as.integer(original$subject_id), as.integer(expected_subject)) &&
        same_number(original$time, panel$time, 1e-12) &&
        identical(as.integer(original$state), as.integer(panel$state)) &&
        identical(as.integer(original$obstype), as.integer(panel$obstype)) &&
        saved$n_people == length(ids) && saved$n_records == nrow(panel),
        paste0("participant_n=", length(ids), "; observation_n=", nrow(panel),
               "; ID_SHA256=", digest(ids, algo = "sha256")))
      mf <- saved$fit$data$mf
      guard(paste0(label, "_model_frame_unchanged_rows"), nrow(mf) == nrow(panel) &&
        identical(as.integer(mf[["(subject)"]]), as.integer(expected_subject)) &&
        same_number(mf[["(time)"]], panel$time, 1e-12) &&
        identical(as.integer(mf[["(state)"]]), as.integer(panel$state)) &&
        identical(as.integer(mf[["(obstype)"]]), as.integer(panel$obstype)))
      variables <- c("restless_bl", bundle$core_covariates)
      guard(paste0(label, "_covariate_formula_matches_contract"),
        setequal(all.vars(saved$fit$covariates), variables))
      expected_rows <- imp[match(as.character(panel$id), as.character(imp$id))]
      covariates_match <- function(a, b) {
        if (is.factor(a) || is.factor(b)) identical(as.character(a), as.character(b))
        else same_number(a, b, 1e-12)
      }
      guard(paste0(label, "_covariates_match_imputation"), all(vapply(variables, function(v)
        v %in% names(mf) && v %in% names(original) && !anyNA(expected_rows[[v]]) &&
        covariates_match(original[[v]], expected_rows[[v]]) &&
        covariates_match(mf[[v]], expected_rows[[v]]), TRUE)))
      extracted[[i]] <- extract_from_full_fit(saved$fit, edges)
      guard(paste0(label, "_independent_Hessian_expansion"), TRUE,
        "Validated convergence, Hessian, Q ordering, constraints, coef(fit), and full covariance.")
      guard(paste0(label, "_cached_block_matches_full_fit"),
        same_number(extracted[[i]]$theta, saved$block$theta, 1e-9) &&
        same_number(extracted[[i]]$covariance, saved$block$covariance, 2e-7))
      guard(paste0(label, "_saved_Hessian_diagnostic"),
        same_number(extracted[[i]]$hessian_min, saved$hessian_min_eigenvalue, 1e-7))
      component_diagnostics[[i]] <- data.table(model = label, imputation = i, n_people = length(ids),
        n_records = nrow(panel), hessian_min = extracted[[i]]$hessian_min,
        hessian_condition = extracted[[i]]$hessian_condition, log_likelihood = extracted[[i]]$likelihood,
        warning_count = length(saved$warnings), warnings = paste(saved$warnings, collapse = " | "))
      rm(saved, original, mf, imp, expected_rows); invisible(gc(FALSE))
    }
    theta <- do.call(cbind, lapply(extracted, `[[`, "theta"))
    pooled <- independent_pool(theta, lapply(extracted, `[[`, "covariance"))
    published <- readRDS(result_paths[5L])
    guard("pooled_signature_and_imputation_count", identical(published$signature, signature) && published$m == 20L)
    for (field in c("Q", "qbar", "U", "B", "T"))
      guard(paste0("Rubin_", field), same_number(pooled[[field]], published[[field]], 2e-7))
    C <- matrix(0, 2L, nrow(edges))
    C[1L, targets[1:2]] <- c(1, -1); C[2L, targets[3:4]] <- c(1, -1)
    guard("pooled_contrast_matrix_and_target_order", same_number(C, published$C, 0) &&
      identical(as.integer(targets), as.integer(published$targets)))
    expected_all <- scalar_statistics(pooled, diag(nrow(edges)))
    expected_target <- expected_all[targets]
    expected_target[, p_BH_4 := p.adjust(p, method = "BH", n = 4L)]
    expected_contrast <- scalar_statistics(pooled, C)
    expected_contrast[, p_BH_2 := p.adjust(p, method = "BH", n = 2L)]
    expected_joint <- rbind(joint_statistic(pooled, diag(nrow(edges))[targets, , drop = FALSE]),
                            joint_statistic(pooled, C))
    all_table <- fread(result_paths[1L]); target_table <- fread(result_paths[2L])
    contrast_table <- fread(result_paths[3L]); joint_table <- fread(result_paths[4L])
    guard("all_edge_CSV_order", identical(all_table$transition, edges$transition))
    guard("target_CSV_order", identical(target_table$transition, edges$transition[targets]))
    guard("contrast_CSV_order_and_definition", identical(contrast_table$origin, c("No disease", "HTN only")) &&
      all(contrast_table$contrast == "Heart minus stroke"))
    guard("joint_CSV_definitions", identical(joint_table$test, c("Four target coefficients = 0",
      "Two heart-stroke coefficient differences = 0")) &&
      all(joint_table$reference == "Asymptotic chi-square Wald on Rubin total covariance"))
    compare_table <- function(label, actual, expected) {
      guard(paste0(label, "_rows_and_columns"), nrow(actual) == nrow(expected) && all(names(expected) %in% names(actual)))
      for (column in names(expected))
        guard(paste0(label, "_", column), same_number(actual[[column]], expected[[column]], 2e-7))
    }
    compare_table("all_edges", all_table, expected_all)
    compare_table("four_targets", target_table, expected_target)
    compare_table("heart_stroke_contrasts", contrast_table, expected_contrast)
    compare_table("joint_tests", joint_table, expected_joint)
    guard("global_gate_and_BH_decisions", same_number(target_table$global_four_edge_p,
      rep(expected_joint$p[1L], 4L), 2e-7) && identical(as.logical(target_table$passes_global_and_BH),
        expected_joint$p[1L] < .05 & expected_target$p_BH_4 < .05))
    publish("COMPUTATIONAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_REVIEW")
    0L
  }, error = function(e) {
    add("review_execution_completed", FALSE, conditionMessage(e))
    publish("FAILED_REVIEW_REQUIRED", conditionMessage(e))
    message("FAILED_REVIEW_REQUIRED: ", conditionMessage(e))
    1L
  })
  status
}

arguments <- commandArgs(trailingOnly = TRUE)
if (identical(arguments, "--self-test")) {
  self_test()
} else {
  if (length(arguments)) stop("Use no arguments for real-data review, or --self-test", call. = FALSE)
  exit_code <- run_review()
  quit(save = "no", status = exit_code)
}
