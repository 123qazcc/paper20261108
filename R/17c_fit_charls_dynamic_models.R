#!/usr/bin/env Rscript
# CHARLS stage 3 continuous-time panel models for the four-level W1-W2 dynamic
# restless-sleep pattern. Full fits are stored separately from aggregate results.
options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.treatment"))
suppressPackageStartupMessages({library(data.table); library(msm); library(digest); library(jsonlite)})
setDTthreads(1L)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L, args[[1]] %in% c("validate", "component", "pool"))
root <- "outputs/charls_dynamic_sleep_20260920"
stage2 <- "outputs/charls_cesd9_duration_20260914"
script <- "R/17c_fit_charls_dynamic_models.R"
bundle_path <- file.path(root, "data/charls_dynamic_model_bundle.rds")
review_path <- file.path(root, "model_qc/dynamic_bundle_review.json")
control_path <- "analysis/charls_stage3_dynamic_sleep_control.json"
control <- read_json(control_path, simplifyVector = TRUE)
stopifnot(control$research_stage == 3L, isTRUE(control$stage_enabled),
  isTRUE(control$stop_after_stage), identical(control$automatic_next_stage, FALSE),
  identical(control$allow_later_stages, FALSE),
  identical(as.character(control$models), c("CORE", "BOTH")), control$required_model_fits == 40L)
model_root <- control$model_directory
for (d in c("model_qc", "results", "pooled", "logs")) dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
dir.create(model_root, recursive = TRUE, showWarnings = FALSE)

sha <- function(p) digest(file = p, algo = "sha256")
save_atomic <- function(x, p) {
  if (file.exists(p)) stop("Existing object retained; refusing overwrite: ", p)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile("charls_dynamic_fit_", tmpdir = dirname(p), fileext = ".rds")
  saveRDS(x, tmp, compress = "gzip")
  if (!file.rename(tmp, p)) stop("Atomic publication failed; retained ", tmp)
}
bun <- readRDS(bundle_path)
review <- read_json(review_path, simplifyVector = TRUE)
stopifnot(identical(bun$status, "APPROVED_FOR_DYNAMIC_MODELS"), all(bun$preparation_checks$pass),
  identical(review$status, "PREPARATION_NUMERICAL_CHECKS_PASS"), review$checks_total == review$checks_passed,
  identical(review$bundle_sha256, sha(bundle_path)), identical(review$review_script_sha256, sha("R/17b_verify_charls_dynamic_bundle.R")),
  identical(review$review_csv_sha256, sha(file.path(root, "model_qc/dynamic_bundle_review.csv"))))
stopifnot(all(file.exists(names(bun$input_hashes))),
  identical(unname(vapply(names(bun$input_hashes), sha, "")), unname(as.character(bun$input_hashes))))

edges <- as.data.table(bun$edges)
k <- nrow(edges)
stopifnot(k == 20L, identical(as.character(bun$pattern_levels), as.character(control$pattern_levels)))
target_names <- as.character(control$target_transitions)
target <- match(target_names, edges$transition)
stopifnot(!anyNA(target), length(unique(target)) == 4L)
pattern_labels <- as.character(bun$pattern_coefficient_labels)
stopifnot(identical(pattern_labels, c("patternonset", "patternremission", "patternpersistent_high")))
bundle_sha <- sha(bundle_path); script_sha <- sha(script)
signature_for <- function(spec) digest(list(bundle = bundle_sha, script = script_sha,
  msm = as.character(packageVersion("msm")), spec = spec), algo = "sha256")

make_sample <- function(spec, imp) {
  covs <- bun$model_covariates[[spec]]
  b <- as.data.table(copy(bun$imputations[[imp]]))
  stopifnot(nrow(b) == control$expected_people, !anyDuplicated(b$id), identical(b$id, bun$baseline$id),
    !anyNA(b[, c("id", "pattern", covs), with = FALSE]))
  z <- merge(bun$panel, b[, c("id", "pattern", covs), with = FALSE], by = "id", sort = FALSE)
  setorder(z, id, time)
  z[, subject_id := .GRP, by = id]
  stopifnot(nrow(z) == control$expected_panel_records, uniqueN(z$id) == control$expected_people,
    !anyNA(z[, c("id", "time", "state", "pattern", covs), with = FALSE]))
  as.data.frame(z)
}
design <- function(spec, z) {
  form <- reformulate(c("pattern", bun$model_covariates[[spec]]))
  mm <- model.matrix(form, z)
  labels <- setdiff(colnames(mm), "(Intercept)")
  stopifnot(all(is.finite(mm)), qr(mm)$rank == ncol(mm), all(pattern_labels %in% labels),
    identical(labels[match(pattern_labels, labels)], pattern_labels))
  constraints <- setNames(lapply(labels, function(v)
    if (v %in% pattern_labels) edges$exposure_group else rep(1L, k)), labels)
  list(formula = form, labels = labels, constraints = constraints)
}
warm_source <- function(spec, imp) {
  p <- file.path("private/model_store/CHARLS_CESD9_duration_20260914/models", spec, sprintf("MI_%02d.rds", imp))
  candidates <- names(bun$warm_start_hashes)
  j <- match(tolower(normalizePath(p, winslash = "/", mustWork = TRUE)),
    tolower(normalizePath(candidates, winslash = "/", mustWork = TRUE)))
  stopifnot(!is.na(j), identical(sha(p), unname(as.character(bun$warm_start_hashes[[j]]))))
  obj <- readRDS(p)
  stopifnot(identical(obj$status, "ACCEPTED"), identical(obj$spec, spec), as.integer(obj$imputation) == imp,
    obj$n_people == 12036L, obj$n_records == 53786L, inherits(obj$fit, "msm"))
  obj
}
warm_values <- function(warm, labels) {
  fit <- warm$fit
  q <- msm::qmatrix.msm(fit, covariates = "mean", ci = "none")
  diag(q) <- 0
  starts <- setNames(lapply(labels, function(v) {
    col <- match(v, fit$qcmodel$covlabels)
    if (is.na(col)) rep(0, k) else as.numeric(fit$paramdata$params[k + (col - 1L) * k + seq_len(k)])
  }), labels)
  stopifnot(all(vapply(starts, function(x) length(x) == k && all(is.finite(x)), TRUE)))
  list(q = q, starts = starts)
}
inspect_fit <- function(fit) {
  stopifnot(inherits(fit, "msm"), identical(as.integer(fit$opt$convergence), 0L), isTRUE(fit$foundse),
    is.finite(as.numeric(logLik(fit))), is.matrix(fit$covmat), all(is.finite(fit$covmat)),
    as.integer(fit$qmodel$npars) == k)
  h <- (fit$opt$hessian + t(fit$opt$hessian)) / 2
  hev <- eigen(h, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(all(is.finite(hev)), min(hev) > 1e-8)
  cols <- match(pattern_labels, fit$qcmodel$covlabels)
  stopifnot(!anyNA(cols), length(unique(cols)) == 3L)
  idx <- unlist(lapply(cols, function(col) k + (col - 1L) * k + seq_len(k)), use.names = FALSE)
  theta <- as.numeric(fit$paramdata$params[idx])
  V <- fit$covmat[idx, idx, drop = FALSE]
  stopifnot(length(theta) == 3L * k, all(is.finite(theta)), all(is.finite(V)),
    max(abs(V - t(V))) < 1e-8, all(diag(V) >= 0))
  V <- (V + t(V)) / 2
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(min(ev) >= -1e-8 * max(1, max(abs(ev))))
  for (j in seq_along(pattern_labels)) {
    cm <- coef(fit)[[pattern_labels[j]]]
    expected <- vapply(seq_len(k), function(e) cm[edges$from_index[e], edges$to_index[e]], 0.0)
    stopifnot(isTRUE(all.equal(theta[((j - 1L) * k + 1L):(j * k)], expected, tolerance = 1e-10)))
  }
  list(block = list(theta = theta, covariance = V, labels = pattern_labels,
    parameter_names = unlist(lapply(pattern_labels, function(v) paste(v, edges$transition, sep = "::")))),
    hessian_min_eigenvalue = min(hev))
}
verify_saved <- function(obj, spec, imp, z, des) {
  stopifnot(identical(obj$status, "ACCEPTED"), identical(obj$signature, signature_for(spec)),
    identical(obj$spec, spec), as.integer(obj$imputation) == imp,
    obj$n_people == control$expected_people, obj$n_records == control$expected_panel_records,
    identical(obj$constraints, des$constraints), identical(obj$model_input_hash, digest(z, algo = "sha256")),
    identical(obj$draw_hash, digest(as.data.frame(bun$imputations[[imp]]), algo = "sha256")),
    identical(obj$bundle_sha256, bundle_sha), identical(obj$script_sha256, script_sha))
  audited <- inspect_fit(obj$fit)
  stopifnot(isTRUE(all.equal(audited$block, obj$block, tolerance = 1e-10)),
    identical(as.character(obj$fit$qcmodel$covlabels), des$labels))
  invisible(TRUE)
}

fit_one <- function(spec, imp) {
  label <- sprintf("%s_MI_%02d", spec, imp)
  path <- file.path(model_root, spec, sprintf("MI_%02d.rds", imp))
  z <- make_sample(spec, imp); des <- design(spec, z)
  if (file.exists(path)) {
    obj <- readRDS(path); verify_saved(obj, spec, imp, z, des)
    message("Verified accepted cache: ", label)
    return(invisible(obj))
  }
  warm <- warm_source(spec, imp); init <- warm_values(warm, des$labels)
  begin <- Sys.time(); warnings_seen <- character(); error <- NULL; fit <- NULL; audited <- NULL
  status_path <- file.path(root, "model_qc", paste0(label, "_status.csv"))
  fwrite(data.table(model = label, status = "RUNNING", pid = Sys.getpid(), started = as.character(begin),
    signature = signature_for(spec), model_file = path), status_path)
  fitargs <- list(formula = state ~ time, subject = quote(subject_id), data = z, qmatrix = init$q,
    covariates = des$formula, constraint = des$constraints, covinits = init$starts,
    obstype = z$obstype, center = TRUE, hessian = TRUE,
    control = list(fnscale = 10000, maxit = 2000, reltol = 1e-9, trace = 0))
  if (length(bun$censor_codes)) {
    fitargs$censor <- bun$censor_codes
    fitargs$censor.states <- bun$censor_states
  }
  fit <- tryCatch(withCallingHandlers(do.call(msm::msm, fitargs), warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) { error <<- conditionMessage(e); NULL })
  if (!is.null(fit)) audited <- tryCatch(inspect_fit(fit), error = function(e) { error <<- conditionMessage(e); NULL })
  accepted <- !is.null(audited)
  obj <- list(status = if (accepted) "ACCEPTED" else "REJECTED_REVIEW_REQUIRED",
    signature = signature_for(spec), spec = spec, model = sprintf("MI_%02d", imp), imputation = imp,
    fit = fit, block = if (accepted) audited$block else NULL, constraints = des$constraints,
    q_at_means = if (accepted) msm::qmatrix.msm(fit, covariates = "mean", ci = "none") else NULL,
    model_input_hash = digest(z, algo = "sha256"),
    draw_hash = digest(as.data.frame(bun$imputations[[imp]]), algo = "sha256"),
    warm_source_sha256 = sha(file.path("private/model_store/CHARLS_CESD9_duration_20260914/models", spec, sprintf("MI_%02d.rds", imp))),
    bundle_sha256 = bundle_sha, script_sha256 = script_sha, review_sha256 = sha(review_path),
    warnings = unique(warnings_seen), error = error, n_people = uniqueN(z$id), n_records = nrow(z),
    hessian_min_eigenvalue = if (accepted) audited$hessian_min_eigenvalue else NA_real_,
    elapsed_seconds = as.numeric(difftime(Sys.time(), begin, units = "secs")) )
  if (accepted) verify_saved(obj, spec, imp, z, des)
  save_atomic(obj, path)
  fwrite(data.table(model = label, status = obj$status, pid = Sys.getpid(), n_people = obj$n_people,
    n_records = obj$n_records, elapsed_seconds = obj$elapsed_seconds,
    convergence = if (is.null(fit)) NA_integer_ else fit$opt$convergence,
    found_se = if (is.null(fit)) FALSE else fit$foundse,
    hessian_min_eigenvalue = obj$hessian_min_eigenvalue,
    warnings = paste(unique(warnings_seen), collapse = "; "), error = if (is.null(error)) "" else error,
    signature = obj$signature, model_file = path, model_sha256 = sha(path)), status_path, bom = TRUE)
  message(label, ": ", obj$status, "; seconds=", round(obj$elapsed_seconds, 1))
  if (!accepted) stop("Rejected fit retained for review; no automatic retry or omission: ", label)
  invisible(obj)
}

pool_spec <- function(spec) {
  paths <- file.path(model_root, spec, sprintf("MI_%02d.rds", 1:20))
  stopifnot(all(file.exists(paths)))
  fits <- lapply(paths, readRDS)
  stopifnot(all(vapply(seq_along(fits), function(i) identical(fits[[i]]$status, "ACCEPTED") &&
    identical(fits[[i]]$signature, signature_for(spec)) && as.integer(fits[[i]]$imputation) == i, TRUE)))
  Q <- sapply(fits, function(x) x$block$theta)
  covs <- lapply(fits, function(x) x$block$covariance)
  m <- ncol(Q); qbar <- rowMeans(Q); U <- Reduce(`+`, covs) / m; B <- cov(t(Q)); T <- U + (1 + 1 / m) * B
  scalar <- function(L) {
    est <- drop(L %*% qbar); u <- diag(L %*% U %*% t(L)); b <- diag(L %*% B %*% t(L))
    se <- sqrt(diag(L %*% T %*% t(L))); r <- (1 + 1 / m) * b / u
    df <- ifelse(b > 1e-15 & u > 0, (m - 1) * (1 + 1 / r)^2, Inf)
    data.table(log_effect = est, se = se, df = df, HR = exp(est),
      lower = exp(est - qt(.975, df) * se), upper = exp(est + qt(.975, df) * se),
      p = 2 * pt(-abs(est / se), df), within_U = u, between_B = b, total_variance = se^2)
  }
  joint <- function(L, label, family) {
    V <- L %*% T %*% t(L); V <- (V + t(V)) / 2
    ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    stopifnot(all(is.finite(ev)), min(ev) > 1e-12 * max(ev), qr(V)$rank == nrow(V))
    a <- drop(L %*% qbar); stat <- drop(t(a) %*% solve(V, a))
    data.table(test = label, family = family, chi_square = stat, df = nrow(L),
      p = pchisq(stat, nrow(L), lower.tail = FALSE),
      reference = "Asymptotic chi-square Wald on Rubin total covariance")
  }
  npar <- 3L * k
  raw_meta <- cbind(data.table(pattern_coefficient = rep(pattern_labels, each = k)),
    edges[rep(seq_len(k), times = 3L)])
  raw <- cbind(raw_meta, scalar(diag(npar)))

  defs <- list(onset_vs_persistent_low = c(1, 0, 0),
    persistent_high_vs_remission = c(0, -1, 1),
    remission_vs_persistent_low = c(0, 1, 0),
    persistent_high_vs_persistent_low = c(0, 0, 1))
  roles <- c(onset_vs_persistent_low = "primary", persistent_high_vs_remission = "primary",
    remission_vs_persistent_low = "contextual", persistent_high_vs_persistent_low = "contextual")
  build_target_L <- function(which_defs) {
    meta <- CJ(comparison = which_defs, edge_order = seq_along(target_names), sorted = FALSE)
    meta[, `:=`(transition = target_names[edge_order], target_index = target[edge_order], role = roles[comparison])]
    L <- matrix(0, nrow(meta), npar)
    for (r in seq_len(nrow(meta))) for (j in 1:3)
      L[r, (j - 1L) * k + meta$target_index[r]] <- defs[[meta$comparison[r]]][j]
    list(meta = meta, L = L)
  }
  all_defs <- names(defs)
  tc <- build_target_L(all_defs)
  target_tab <- cbind(tc$meta, scalar(tc$L))
  target_tab[, p_BH_family := p.adjust(p, "BH"), by = role]

  omnibus <- rbindlist(lapply(seq_along(target), function(e) {
    L <- matrix(0, 3, npar)
    for (j in 1:3) L[j, (j - 1L) * k + target[e]] <- 1
    x <- joint(L, paste0("Four-level pattern omnibus: ", target_names[e]), "edge_pattern_omnibus")
    x[, `:=`(transition = target_names[e], edge_order = e)]
    x
  }))
  omnibus[, p_BH_4 := p.adjust(p, "BH")]

  primary_L <- tc$L[tc$meta$role == "primary", , drop = FALSE]
  contextual_L <- tc$L[tc$meta$role == "contextual", , drop = FALSE]
  all_target_L <- matrix(0, 12L, npar); rr <- 0L
  for (j in 1:3) for (e in target) { rr <- rr + 1L; all_target_L[rr, (j - 1L) * k + e] <- 1 }
  joint_tab <- rbind(joint(primary_L, "Eight primary dynamic target contrasts = 0", "primary_target"),
    joint(contextual_L, "Eight contextual target contrasts = 0", "contextual_target"),
    joint(all_target_L, "All 12 four-level-pattern target coefficients = 0", "four_level_target"))

  origins <- data.table(origin = c("No disease", "HTN only"), heart_edge = target[c(1, 3)], stroke_edge = target[c(2, 4)])
  direct_meta <- CJ(comparison = all_defs, origin_order = 1:2, sorted = FALSE)
  direct_meta[, `:=`(origin = origins$origin[origin_order], role = roles[comparison])]
  direct_L <- matrix(0, nrow(direct_meta), npar)
  for (r in seq_len(nrow(direct_meta))) {
    o <- direct_meta$origin_order[r]
    for (j in 1:3) {
      w <- defs[[direct_meta$comparison[r]]][j]
      direct_L[r, (j - 1L) * k + origins$heart_edge[o]] <- w
      direct_L[r, (j - 1L) * k + origins$stroke_edge[o]] <- -w
    }
  }
  direct_tab <- cbind(direct_meta, scalar(direct_L))
  direct_tab[, p_BH_family := p.adjust(p, "BH"), by = role]
  direct_joint <- rbind(joint(direct_L[direct_meta$role == "primary", , drop = FALSE],
    "Four primary heart-minus-stroke dynamic contrasts = 0", "primary_direct"),
    joint(direct_L[direct_meta$role == "contextual", , drop = FALSE],
      "Four contextual heart-minus-stroke contrasts = 0", "contextual_direct"))

  outdir <- file.path(root, "pooled", spec, "results")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(Q = Q, qbar = qbar, U = U, B = B, T = T, pattern_labels = pattern_labels,
    edges = edges, target_contrast_matrix = tc$L, target_contrast_meta = tc$meta,
    direct_contrast_matrix = direct_L, direct_contrast_meta = direct_meta),
    file.path(outdir, "rubin_matrices_MI20.rds"))
  fwrite(raw, file.path(outdir, "all_pattern_coefficients_MI20.csv"), bom = TRUE)
  fwrite(target_tab, file.path(outdir, "four_target_dynamic_contrasts_MI20.csv"), bom = TRUE)
  fwrite(omnibus, file.path(outdir, "four_target_pattern_omnibus_MI20.csv"), bom = TRUE)
  fwrite(direct_tab, file.path(outdir, "heart_stroke_dynamic_direct_MI20.csv"), bom = TRUE)
  fwrite(rbind(joint_tab, direct_joint, fill = TRUE), file.path(outdir, "joint_tests_MI20.csv"), bom = TRUE)
  list(target = target_tab, omnibus = omnibus, direct = direct_tab,
    joint = rbind(joint_tab, direct_joint, fill = TRUE), paths = paths)
}

pool_all <- function() {
  manifest <- list()
  for (spec in control$models) for (i in 1:20) {
    p <- file.path(model_root, spec, sprintf("MI_%02d.rds", i)); stopifnot(file.exists(p))
    z <- make_sample(spec, i); des <- design(spec, z); obj <- readRDS(p); verify_saved(obj, spec, i, z, des)
    manifest[[length(manifest) + 1L]] <- data.table(spec = spec, imputation = i, path = p, sha256 = sha(p),
      n_people = obj$n_people, n_records = obj$n_records, elapsed_seconds = obj$elapsed_seconds)
    rm(z, obj)
  }
  pooled <- lapply(control$models, pool_spec); names(pooled) <- control$models
  for (name in c("four_target_dynamic_contrasts_MI20", "four_target_pattern_omnibus_MI20",
    "heart_stroke_dynamic_direct_MI20", "joint_tests_MI20", "all_pattern_coefficients_MI20")) {
    combined <- rbindlist(lapply(control$models, function(spec) {
      x <- fread(file.path(root, "pooled", spec, "results", paste0(name, ".csv")))
      x[, model := spec]; x
    }), fill = TRUE)
    fwrite(combined, file.path(root, "results", paste0(name, ".csv")), bom = TRUE)
  }
  core <- pooled$CORE$target; both <- pooled$BOTH$target
  stopifnot(identical(core$comparison, both$comparison), identical(core$transition, both$transition))
  fwrite(data.table(comparison = core$comparison, role = core$role, transition = core$transition,
    CORE_logHR = core$log_effect, BOTH_logHR = both$log_effect,
    change_logHR = both$log_effect - core$log_effect, CORE_HR = core$HR, BOTH_HR = both$HR,
    interpretation = "Descriptive correlated-model change; no mediation proportion or change P/CI"),
    file.path(root, "results/matched_CORE_BOTH_change_descriptive.csv"), bom = TRUE)
  fwrite(rbindlist(manifest), file.path(root, "model_qc/accepted_40_model_manifest.csv"), bom = TRUE)
  write_json(list(status = "POOLED_AWAITING_INDEPENDENT_REVIEW_AND_FIGURES", accepted_fits = 40L,
    models = as.character(control$models), m = 20L, people = control$expected_people,
    panel_records = control$expected_panel_records, automatic_next_stage = FALSE),
    file.path(root, "model_qc/pooled_status.json"), pretty = TRUE, auto_unbox = TRUE)
  message("POOLED_AWAITING_INDEPENDENT_REVIEW_AND_FIGURES: 40/40 accepted.")
}

if (args[[1]] == "validate") {
  stopifnot(length(args) == 1L)
  for (spec in control$models) for (i in 1:20) {
    z <- make_sample(spec, i); d <- design(spec, z); w <- warm_source(spec, i); init <- warm_values(w, d$labels)
    stopifnot(length(d$constraints) == length(d$labels), all(vapply(d$constraints, length, 1L) == k),
      identical(dim(init$q), c(9L, 9L)))
  }
  write_json(list(status = "FIT_PREFLIGHT_PASS", required_fits = 40L,
    model_directory = model_root, bundle_sha256 = bundle_sha, script_sha256 = script_sha,
    no_models_fitted = TRUE), file.path(root, "model_qc/fit_preflight.json"), pretty = TRUE, auto_unbox = TRUE)
  message("FIT_PREFLIGHT_PASS: 40 dynamic designs; no fit started.")
} else if (args[[1]] == "component") {
  stopifnot(length(args) == 3L, args[[2]] %in% control$models)
  imp <- as.integer(args[[3]]); stopifnot(!is.na(imp), imp %in% 1:20)
  fit_one(args[[2]], imp)
} else {
  stopifnot(length(args) == 1L)
  pool_all()
}
