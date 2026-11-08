#!/usr/bin/env Rscript
# Independent numerical/data audit of the CHARLS dynamic-pattern preparation.
# This script never fits or interprets a disease-transition association model.
options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.treatment"))
suppressPackageStartupMessages({library(data.table); library(mice); library(digest); library(jsonlite)})
setDTthreads(1L)

out <- "outputs/charls_dynamic_sleep_20260920"
bundle_path <- file.path(out, "data/charls_dynamic_model_bundle.rds")
mids_path <- file.path(out, "data/charls_dynamic_mice20.rds")
pre_path <- file.path(out, "data/charls_dynamic_preimputation.rds")
control_path <- "analysis/charls_stage3_dynamic_sleep_control.json"
observed_path <- "outputs/charls_external_extension_20260913/data/charls_observed_exposures_tristate_panel.rds"
conflict_path <- "outputs/charls_external_extension_20260913/secondary_readiness/internal/living_status_conflict_locators.rds"
qc_dir <- file.path(out, "model_qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
sha <- function(p) digest(file = p, algo = "sha256")
same <- function(a, b, tol = 1e-10) {
  if (length(a) != length(b) || !identical(dim(a), dim(b))) return(FALSE)
  if (is.factor(a) || is.factor(b) || is.character(a) || is.character(b))
    return(identical(as.character(a), as.character(b)))
  if (!identical(is.na(a), is.na(b))) return(FALSE)
  keep <- !is.na(a); a <- as.numeric(a[keep]); b <- as.numeric(b[keep])
  all(abs(a - b) <= tol * pmax(1, abs(a), abs(b)))
}
checks <- list()
add <- function(name, pass, detail = "") {
  checks[[length(checks) + 1L]] <<- data.table(check = name, pass = isTRUE(pass), detail = as.character(detail))
  if (!isTRUE(pass)) stop(paste(name, detail, sep = ": "), call. = FALSE)
}
publish <- function(status, reason = "") {
  tab <- if (length(checks)) rbindlist(checks, fill = TRUE) else data.table(check = "execution", pass = FALSE, detail = reason)
  fwrite(tab, file.path(qc_dir, "dynamic_bundle_review.csv"), bom = TRUE)
  write_json(list(status = if (all(tab$pass) && status == "PREPARATION_NUMERICAL_CHECKS_PASS") status else "NOT_READY_NO_MODELS_AUTHORIZED",
    checks_total = nrow(tab), checks_passed = sum(tab$pass), reason = reason,
    bundle_sha256 = if (file.exists(bundle_path)) sha(bundle_path) else NA_character_,
    mids_sha256 = if (file.exists(mids_path)) sha(mids_path) else NA_character_,
    preimputation_sha256 = if (file.exists(pre_path)) sha(pre_path) else NA_character_,
    review_script_sha256 = sha("R/17b_verify_charls_dynamic_bundle.R"),
    review_csv_sha256 = sha(file.path(qc_dir, "dynamic_bundle_review.csv")),
    substantive_assumptions_proven = FALSE,
    updated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
    file.path(qc_dir, "dynamic_bundle_review.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  writeLines(c("# CHARLS 动态睡眠模式准备独立核验", "", paste0("状态：", status),
    paste0("检查：", sum(tab$pass), "/", nrow(tab), " 项通过。"), if (nzchar(reason)) paste0("原因：", reason) else "", "",
    "本脚本独立重建 W1–W2 模式和 W2 landmark，核对人/记录/死亡/部分状态分母、两条争议终点剔除、MICE 原观察值保留、20份完成数据、CESD-9/时长样条及模型矩阵。未拟合或查看疾病关联效应。", "",
    "数值通过不证明 MAR、无选择偏倚、无反向因果、Markov/时间齐性、无未测混杂或因果识别。"),
    file.path(qc_dir, "dynamic_bundle_review.md"), useBytes = TRUE)
}

run <- function() {
  required <- c(bundle_path, mids_path, pre_path, control_path, observed_path, conflict_path)
  add("all_review_inputs_exist", all(file.exists(required)), paste(required[!file.exists(required)], collapse = " | "))
  initial <- setNames(vapply(required, sha, ""), required)
  control <- fromJSON(control_path, simplifyVector = TRUE)
  bun <- readRDS(bundle_path); imp <- readRDS(mids_path); pre <- readRDS(pre_path)
  add("stage_scope_check", control$research_stage == 3L && isTRUE(control$stage_enabled) &&
    isTRUE(control$stop_after_stage) && identical(control$automatic_next_stage, FALSE) &&
    identical(control$allow_later_stages, FALSE))
  add("approved_bundle_identity", identical(bun$status, "APPROVED_FOR_DYNAMIC_MODELS") && bun$research_stage == 3L &&
    bun$m == 20L && bun$seed == control$seed && all(bun$preparation_checks$pass))
  add("expected_model_batch", identical(as.character(control$models), c("CORE", "BOTH")) &&
    control$required_model_fits == 40L && identical(control$primary_model, "CORE"))

  input_hashes <- unlist(bun$input_hashes, use.names = TRUE)
  add("all_signed_inputs_current", length(input_hashes) > 0L && all(file.exists(names(input_hashes))) &&
    identical(unname(vapply(names(input_hashes), sha, "")), unname(as.character(input_hashes))))
  warm_hashes <- unlist(bun$warm_start_hashes, use.names = TRUE)
  add("all_40_warm_sources_current", length(warm_hashes) == 40L && all(file.exists(names(warm_hashes))) &&
    identical(unname(vapply(names(warm_hashes), sha, "")), unname(as.character(warm_hashes))))

  p <- as.data.table(bun$panel); lm <- as.data.table(bun$landmark); b <- as.data.table(bun$baseline)
  add("frozen_sample_dimensions", nrow(lm) == control$expected_people && nrow(p) == control$expected_panel_records &&
    uniqueN(p$id) == control$expected_people && nrow(b) == control$expected_people)
  add("ID_order_and_uniqueness", !anyDuplicated(lm$id) && identical(lm$id, b$id) && setequal(lm$id, p$id))
  add("entry_state_and_time", nrow(p[time == 0]) == nrow(lm) && all(p[time == 0, state] %in% 1:8) &&
    all(p[, min(time), by = id]$V1 == 0) && all(p[, any(time > 0), by = id]$V1))
  add("death_partial_denominators", sum(p$state == 9L) == control$expected_interval_deaths &&
    sum(p$state >= 100L) == control$expected_partial_records)
  end <- p[, .(death_n = sum(state == 9L), last_state = last(state)), by = id]
  add("death_once_terminal", all(end$death_n <= 1L & (end$death_n == 0L | end$last_state == 9L)))
  add("pattern_levels_and_counts", identical(levels(lm$pattern), as.character(control$pattern_levels)) &&
    identical(as.integer(lm[, .N, by = pattern][match(control$pattern_levels, as.character(pattern)), N]),
      as.integer(control$pattern_counts)))

  conflicts <- as.data.table(readRDS(conflict_path))
  add("disputed_W4_rows_absent", nrow(conflicts) == 2L &&
    p[conflicts[, .(id, wave)], on = .(id, wave), nomatch = 0L, .N] == 0L)
  obs <- as.data.table(readRDS(observed_path))
  sl <- obs[id %in% lm$id & wave %in% 1:2, .(id = as.character(id), wave, score = restless_score_observed)]
  sl <- dcast(sl, id ~ wave, value.var = "score"); setnames(sl, c("1", "2"), c("w1", "w2")); setorder(sl, id)
  expected_pattern <- c("persistent_low", "onset", "remission", "persistent_high")[
    1L + 2L * as.integer(sl$w1 >= 2) + as.integer(sl$w2 >= 2)]
  add("independent_pattern_reconstruction", identical(sl$id, lm$id) &&
    same(sl$w1, lm$w1_score, 0) && same(sl$w2, lm$w2_score, 0) &&
    identical(expected_pattern, as.character(lm$pattern)))

  add("mids_and_preimputation_identity", inherits(imp, "mids") && imp$m == 20L && imp$iteration == 20L &&
    identical(as.data.frame(imp$data), as.data.frame(pre$data)) &&
    identical(as.character(pre$ids), as.character(lm$id)))
  add("method_and_predictor_identity", identical(as.character(imp$method), as.character(pre$method)) &&
    identical(dimnames(imp$predictorMatrix), dimnames(pre$predictor_matrix)) &&
    same(imp$predictorMatrix, pre$predictor_matrix, 0))
  add("dynamic_fields_fixed_not_imputed", all(imp$method[c("restless_bl", "pattern", "w2_state_factor", "w2_score", "w2_landmark_time")] == "") &&
    !anyNA(imp$data[, c("restless_bl", "pattern", "w2_state_factor", "w2_score", "w2_landmark_time")]))
  add("no_mice_logged_events", is.null(imp$loggedEvents) || NROW(imp$loggedEvents) == 0L)

  core <- as.character(bun$core_covariates); items <- as.character(bun$raw_item_names)
  expected_specs <- list(CORE = core, BOTH = c(core, "cesd9_5", paste0("night_ns", 1:3)))
  add("model_covariate_sets_exact", identical(names(bun$model_covariates), names(expected_specs)) &&
    all(vapply(names(expected_specs), function(s) identical(as.character(bun$model_covariates[[s]]), expected_specs[[s]]), TRUE)))
  add("factor_labels_exact", identical(as.character(bun$pattern_coefficient_labels),
    c("patternonset", "patternremission", "patternpersistent_high")))
  add("spline_frozen", same(bun$spline_spec$knots, c(6, 7), 0) && same(bun$spline_spec$boundary, c(1, 13), 0))

  support_rows <- list(); rank_rows <- list()
  for (i in 1:20) {
    z <- as.data.frame(complete(imp, i)); draw <- as.data.table(bun$imputations[[i]])
    add(sprintf("MI_%02d_complete", i), nrow(z) == control$expected_people && !anyNA(z))
    add(sprintf("MI_%02d_observed_preserved", i), all(vapply(names(imp$data), function(v) {
      keep <- !is.na(imp$data[[v]]); same(z[[v]][keep], imp$data[[v]][keep], 0)
    }, TRUE)))
    add(sprintf("MI_%02d_ID_pattern", i), identical(draw$id, lm$id) &&
      identical(as.character(draw$pattern), as.character(lm$pattern)))
    completed_items <- as.matrix(z[, items, drop = FALSE])
    cesd <- rowSums(completed_items)
    ns <- splines::ns(z$night_hours, knots = bun$spline_spec$knots,
      Boundary.knots = bun$spline_spec$boundary, intercept = FALSE)
    add(sprintf("MI_%02d_derived_CESD9", i), all(completed_items %in% 0:3) &&
      same(draw$cesd9_5, cesd / 5, 1e-12))
    add(sprintf("MI_%02d_derived_duration_spline", i), all(vapply(1:3, function(j)
      same(draw[[paste0("night_ns", j)]], ns[, j], 1e-12), TRUE)))
    for (spec in names(expected_specs)) {
      mm <- model.matrix(reformulate(c("pattern", expected_specs[[spec]])), data = draw)
      full <- nrow(mm) == control$expected_people && all(is.finite(mm)) && qr(mm)$rank == ncol(mm)
      rank_rows[[length(rank_rows) + 1L]] <- data.table(imputation = i, model = spec,
        people = nrow(mm), columns = ncol(mm), rank = qr(mm)$rank, full_rank = full)
      add(sprintf("MI_%02d_%s_rank", i, spec), full)
    }
    for (v in names(imp$method)[nzchar(imp$method)]) {
      donor <- imp$data[[v]][!is.na(imp$data[[v]])]; fill <- z[[v]][is.na(imp$data[[v]])]
      ok <- if (imp$method[[v]] == "pmm") all(fill %in% donor) else if (imp$method[[v]] == "polyreg")
        all(as.character(fill) %in% as.character(1:4)) else all(as.character(fill) %in% c("0", "1"))
      support_rows[[length(support_rows) + 1L]] <- data.table(imputation = i, variable = v,
        method = imp$method[[v]], missing = sum(is.na(imp$data[[v]])), support_pass = ok)
      add(sprintf("MI_%02d_support_%s", i, v), ok)
    }
  }
  fwrite(rbindlist(rank_rows), file.path(qc_dir, "dynamic_model_matrix_ranks.csv"), bom = TRUE)
  fwrite(rbindlist(support_rows), file.path(qc_dir, "dynamic_imputation_support.csv"), bom = TRUE)
  add("all_review_inputs_unchanged", identical(unname(vapply(required, sha, "")), unname(initial)))
  publish("PREPARATION_NUMERICAL_CHECKS_PASS")
  message("PREPARATION_NUMERICAL_CHECKS_PASS: ", length(checks), "/", length(checks), " checks.")
}

tryCatch(run(), error = function(e) {
  publish("NOT_READY_NO_MODELS_AUTHORIZED", conditionMessage(e))
  message("NOT_READY_NO_MODELS_AUTHORIZED: ", conditionMessage(e))
  quit(save = "no", status = 1L)
})
