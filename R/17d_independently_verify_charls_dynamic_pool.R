#!/usr/bin/env Rscript
# Independent reconstruction of all dynamic-pattern coefficients, full covariance,
# Rubin pooling, derived contrasts, multiplicity corrections and joint tests.
options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.treatment"))
suppressPackageStartupMessages({library(data.table); library(msm); library(digest); library(jsonlite)})
setDTthreads(1L)

root <- "outputs/charls_dynamic_sleep_20260920"
bundle_path <- file.path(root, "data/charls_dynamic_model_bundle.rds")
control <- fromJSON("analysis/charls_stage3_dynamic_sleep_control.json", simplifyVector = TRUE)
bun <- readRDS(bundle_path); edges <- as.data.table(bun$edges); k <- nrow(edges)
model_root <- control$model_directory
labels <- c("patternonset", "patternremission", "patternpersistent_high")
targets <- match(as.character(control$target_transitions), edges$transition)
qc <- file.path(root, "model_qc"); dir.create(qc, recursive = TRUE, showWarnings = FALSE)
sha <- function(p) digest(file = p, algo = "sha256")
same <- function(a, b, tol = 2e-7) {
  if (length(a) != length(b) || !identical(dim(a), dim(b)) || !identical(is.na(a), is.na(b))) return(FALSE)
  a <- as.numeric(a); b <- as.numeric(b); keep <- !is.na(a)
  all(abs(a[keep] - b[keep]) <= tol * pmax(1, abs(a[keep]), abs(b[keep])))
}
checks <- list()
add <- function(name, pass, detail = "") {
  checks[[length(checks) + 1L]] <<- data.table(check = name, pass = isTRUE(pass), detail = as.character(detail))
  if (!isTRUE(pass)) stop(paste(name, detail, sep = ": "), call. = FALSE)
}
extract <- function(fit) {
  stopifnot(inherits(fit, "msm"), fit$opt$convergence == 0L, isTRUE(fit$foundse), all(is.finite(fit$covmat)))
  cols <- match(labels, fit$qcmodel$covlabels); stopifnot(!anyNA(cols), as.integer(fit$qmodel$npars) == k)
  idx <- unlist(lapply(cols, function(j) k + (j - 1L) * k + seq_len(k)), use.names = FALSE)
  theta <- as.numeric(fit$paramdata$params[idx]); V <- fit$covmat[idx, idx, drop = FALSE]; V <- (V + t(V)) / 2
  stopifnot(all(is.finite(theta)), all(is.finite(V)), min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) >=
    -1e-8 * max(1, max(abs(V))))
  list(theta = theta, covariance = V)
}
pool <- function(Q, covs) {
  m <- ncol(Q); qbar <- rowMeans(Q); U <- Reduce(`+`, covs) / m; B <- cov(t(Q)); T <- U + (1 + 1 / m) * B
  list(Q = Q, qbar = qbar, U = U, B = B, T = T, m = m)
}
scalar <- function(pl, L) {
  est <- drop(L %*% pl$qbar); u <- diag(L %*% pl$U %*% t(L)); b <- diag(L %*% pl$B %*% t(L))
  se <- sqrt(diag(L %*% pl$T %*% t(L))); r <- (1 + 1 / pl$m) * b / u
  df <- ifelse(b > 1e-15 & u > 0, (pl$m - 1) * (1 + 1 / r)^2, Inf)
  data.table(log_effect = est, se = se, df = df, HR = exp(est),
    lower = exp(est - qt(.975, df) * se), upper = exp(est + qt(.975, df) * se),
    p = 2 * pt(-abs(est / se), df), within_U = u, between_B = b, total_variance = se^2)
}
joint <- function(pl, L) {
  V <- L %*% pl$T %*% t(L); V <- (V + t(V)) / 2; a <- drop(L %*% pl$qbar)
  data.table(chi_square = drop(t(a) %*% solve(V, a)), df = nrow(L),
    p = pchisq(drop(t(a) %*% solve(V, a)), nrow(L), lower.tail = FALSE))
}
defs <- list(onset_vs_persistent_low = c(1, 0, 0), persistent_high_vs_remission = c(0, -1, 1),
  remission_vs_persistent_low = c(0, 1, 0), persistent_high_vs_persistent_low = c(0, 0, 1))
roles <- c(onset_vs_persistent_low = "primary", persistent_high_vs_remission = "primary",
  remission_vs_persistent_low = "contextual", persistent_high_vs_persistent_low = "contextual")
build <- function() {
  meta <- CJ(comparison = names(defs), edge_order = 1:4, sorted = FALSE)
  meta[, `:=`(transition = control$target_transitions[edge_order], target_index = targets[edge_order], role = roles[comparison])]
  L <- matrix(0, nrow(meta), 3L * k)
  for (r in seq_len(nrow(meta))) for (j in 1:3) L[r, (j - 1L) * k + meta$target_index[r]] <- defs[[meta$comparison[r]]][j]
  origins <- data.table(origin = c("No disease", "HTN only"), heart_edge = targets[c(1, 3)], stroke_edge = targets[c(2, 4)])
  dm <- CJ(comparison = names(defs), origin_order = 1:2, sorted = FALSE)
  dm[, `:=`(origin = origins$origin[origin_order], role = roles[comparison])]
  DL <- matrix(0, nrow(dm), 3L * k)
  for (r in seq_len(nrow(dm))) for (j in 1:3) {
    w <- defs[[dm$comparison[r]]][j]; o <- dm$origin_order[r]
    DL[r, (j - 1L) * k + origins$heart_edge[o]] <- w
    DL[r, (j - 1L) * k + origins$stroke_edge[o]] <- -w
  }
  list(meta = meta, L = L, direct_meta = dm, direct_L = DL)
}

run <- function() {
  prep_review <- fromJSON(file.path(qc, "dynamic_bundle_review.json"), simplifyVector = TRUE)
  pooled_status <- fromJSON(file.path(qc, "pooled_status.json"), simplifyVector = TRUE)
  add("upstream_reviews_ready", identical(prep_review$status, "PREPARATION_NUMERICAL_CHECKS_PASS") &&
    identical(pooled_status$status, "POOLED_AWAITING_INDEPENDENT_REVIEW_AND_FIGURES") && pooled_status$accepted_fits == 40L)
  add("bundle_and_control_identity", identical(bun$status, "APPROVED_FOR_DYNAMIC_MODELS") &&
    nrow(bun$panel) == control$expected_panel_records && nrow(bun$baseline) == control$expected_people)
  C <- build(); summaries <- list()
  for (spec in control$models) {
    sig <- digest(list(bundle = sha(bundle_path), script = sha("R/17c_fit_charls_dynamic_models.R"),
      msm = as.character(packageVersion("msm")), spec = spec), algo = "sha256")
    theta <- matrix(NA_real_, 3L * k, 20L); covs <- vector("list", 20L)
    for (i in 1:20) {
      p <- file.path(model_root, spec, sprintf("MI_%02d.rds", i)); add(paste0(spec, "_MI", i, "_exists"), file.exists(p))
      obj <- readRDS(p)
      add(paste0(spec, "_MI", i, "_identity"), identical(obj$status, "ACCEPTED") &&
        identical(obj$signature, sig) && as.integer(obj$imputation) == i &&
        obj$n_people == control$expected_people && obj$n_records == control$expected_panel_records)
      z <- merge(bun$panel, as.data.table(bun$imputations[[i]])[, c("id", "pattern", bun$model_covariates[[spec]]), with = FALSE],
        by = "id", sort = FALSE); setorder(z, id, time); z[, subject_id := .GRP, by = id]
      add(paste0(spec, "_MI", i, "_input_hash"), identical(obj$model_input_hash, digest(as.data.frame(z), algo = "sha256")) &&
        identical(obj$draw_hash, digest(as.data.frame(bun$imputations[[i]]), algo = "sha256")))
      ex <- extract(obj$fit)
      add(paste0(spec, "_MI", i, "_saved_block"), same(ex$theta, obj$block$theta, 1e-9) &&
        same(ex$covariance, obj$block$covariance, 2e-7))
      theta[, i] <- ex$theta; covs[[i]] <- ex$covariance
    }
    pl <- pool(theta, covs)
    published <- readRDS(file.path(root, "pooled", spec, "results/rubin_matrices_MI20.rds"))
    for (v in c("Q", "qbar", "U", "B", "T")) add(paste0(spec, "_Rubin_", v), same(pl[[v]], published[[v]], 2e-7))

    target <- cbind(C$meta, scalar(pl, C$L)); target[, p_BH_family := p.adjust(p, "BH"), by = role]
    direct <- cbind(C$direct_meta, scalar(pl, C$direct_L)); direct[, p_BH_family := p.adjust(p, "BH"), by = role]
    pub_target <- fread(file.path(root, "pooled", spec, "results/four_target_dynamic_contrasts_MI20.csv"))
    pub_direct <- fread(file.path(root, "pooled", spec, "results/heart_stroke_dynamic_direct_MI20.csv"))
    for (v in c("log_effect", "se", "df", "HR", "lower", "upper", "p", "within_U", "between_B", "total_variance", "p_BH_family")) {
      add(paste0(spec, "_target_", v), same(target[[v]], pub_target[[v]], 2e-7))
      add(paste0(spec, "_direct_", v), same(direct[[v]], pub_direct[[v]], 2e-7))
    }

    omnibus <- rbindlist(lapply(seq_along(targets), function(e) {
      L <- matrix(0, 3, 3L * k); for (j in 1:3) L[j, (j - 1L) * k + targets[e]] <- 1
      cbind(data.table(transition = control$target_transitions[e]), joint(pl, L))
    }))
    omnibus[, p_BH_4 := p.adjust(p, "BH")]
    pub_omni <- fread(file.path(root, "pooled", spec, "results/four_target_pattern_omnibus_MI20.csv"))
    for (v in c("chi_square", "df", "p", "p_BH_4")) add(paste0(spec, "_omnibus_", v), same(omnibus[[v]], pub_omni[[v]], 2e-7))

    primary_L <- C$L[C$meta$role == "primary", , drop = FALSE]
    contextual_L <- C$L[C$meta$role == "contextual", , drop = FALSE]
    allL <- matrix(0, 12, 3L * k); r <- 0L
    for (j in 1:3) for (e in targets) {r <- r + 1L; allL[r, (j - 1L) * k + e] <- 1}
    recomputed_joint <- rbind(joint(pl, primary_L), joint(pl, contextual_L), joint(pl, allL),
      joint(pl, C$direct_L[C$direct_meta$role == "primary", , drop = FALSE]),
      joint(pl, C$direct_L[C$direct_meta$role == "contextual", , drop = FALSE]))
    pub_joint <- fread(file.path(root, "pooled", spec, "results/joint_tests_MI20.csv"))
    for (v in c("chi_square", "df", "p")) add(paste0(spec, "_joint_", v), same(recomputed_joint[[v]], pub_joint[[v]], 2e-7))
    summaries[[spec]] <- data.table(model = spec, people = control$expected_people, records = control$expected_panel_records,
      primary_global_p = recomputed_joint$p[1], all_pattern_global_p = recomputed_joint$p[3], primary_direct_global_p = recomputed_joint$p[4])
  }
  tab <- rbindlist(checks)
  fwrite(tab, file.path(qc, "dynamic_independent_pool_review.csv"), bom = TRUE)
  fwrite(rbindlist(summaries), file.path(qc, "dynamic_independent_pool_summary.csv"), bom = TRUE)
  write_json(list(status = "COMPUTATIONAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_REVIEW",
    checks_total = nrow(tab), checks_passed = sum(tab$pass), accepted_fits = 40L,
    models = as.character(control$models), bundle_sha256 = sha(bundle_path),
    reviewer_script_sha256 = sha("R/17d_independently_verify_charls_dynamic_pool.R"),
    substantive_assumptions_proven = FALSE), file.path(qc, "dynamic_independent_pool_review.json"),
    pretty = TRUE, auto_unbox = TRUE)
  writeLines(c("# CHARLS 动态睡眠模式独立数值核验", "",
    "状态：COMPUTATIONAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_REVIEW", paste0("检查：", nrow(tab), "/", nrow(tab), " 项通过。"), "",
    "从40个完整拟合对象重新提取60维模式系数及完整协方差，独立重算Rubin U/B/T、四路径派生比较、心脏病－卒中直接比较、BH及联合Wald，并与发布表逐项核对。", "",
    "该核验不证明MAR、动态样本选择可忽略、无反向因果、Markov/时间齐性、无未测混杂或因果效应。"),
    file.path(qc, "dynamic_independent_pool_review.md"), useBytes = TRUE)
  message("COMPUTATIONAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_REVIEW: ", nrow(tab), "/", nrow(tab), " checks.")
}

tryCatch(run(), error = function(e) {
  if (length(checks)) fwrite(rbindlist(checks, fill = TRUE), file.path(qc, "dynamic_independent_pool_review.csv"), bom = TRUE)
  message("FAILED_REVIEW_REQUIRED: ", conditionMessage(e)); quit(save = "no", status = 1L)
})
