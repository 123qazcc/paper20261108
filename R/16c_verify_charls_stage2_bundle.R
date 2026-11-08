#!/usr/bin/env Rscript
# Independent preparation audit only. No model fitting, effect extraction or edits
# to source data, completed datasets, models or scientific results.
options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.treatment"))
suppressPackageStartupMessages({library(data.table); library(mice); library(digest); library(jsonlite)})
setDTthreads(1L)

out <- "outputs/charls_cesd9_duration_20260914"
old_out <- "outputs/charls_external_extension_20260913"
qc <- file.path(out, "model_qc")
bundle_path <- file.path(out, "data/charls_stage2_model_bundle.rds")
mids_path <- file.path(out, "data/charls_stage2_mice20.rds")
old_bundle_path <- file.path(old_out, "data/charls_model_bundle.rds")
old_mids_path <- file.path(old_out, "data/charls_baseline_mice20.rds")
observed_path <- file.path(old_out, "data/charls_observed_exposures_tristate_panel.rds")
contract_path <- "analysis/charls_stage2_cesd9_duration_contract.md"
control_path <- "analysis/charls_stage2_cesd9_duration_control.json"
sha <- function(p) digest(file = p, algo = "sha256")
canonical <- function(p) tolower(normalizePath(p, winslash = "/", mustWork = TRUE))
same_number <- function(a, b, tolerance = 1e-10) {
  if (length(a) != length(b) || !identical(dim(a), dim(b))) return(FALSE)
  a <- as.numeric(a); b <- as.numeric(b)
  if (!identical(is.na(a), is.na(b))) return(FALSE)
  keep <- !is.na(a); a <- a[keep]; b <- b[keep]
  if (!identical(is.infinite(a), is.infinite(b))) return(FALSE)
  infinite <- is.infinite(a)
  if (any(a[infinite] != b[infinite])) return(FALSE)
  a <- a[!infinite]; b <- b[!infinite]
  all(abs(a - b) <= tolerance * pmax(1, abs(a), abs(b)))
}
same_values <- function(a, b) {
  if (is.factor(a) || is.factor(b) || is.character(a) || is.character(b))
    return(identical(as.character(a), as.character(b)))
  same_number(a, b, 1e-12)
}

checks <- list()
add <- function(name, pass, detail = "") {
  checks[[length(checks) + 1L]] <<- data.table(check = name, pass = isTRUE(pass), detail = as.character(detail))
}
guard <- function(name, value, detail = "") {
  add(name, value, detail)
  if (!isTRUE(value)) stop(paste0(name, if (nzchar(detail)) paste0(": ", detail)), call. = FALSE)
}
publish <- function(status, reason = "") {
  dir.create(qc, recursive = TRUE, showWarnings = FALSE)
  review <- rbindlist(checks, fill = TRUE)
  if (!nrow(review)) review <- data.table(check = "execution", pass = FALSE, detail = reason)
  fwrite(review, file.path(qc, "model_bundle_review.csv"), bom = TRUE)
  write_json(list(status=if(all(review$pass)&&grepl("^PREPARATION_NUMERICAL_CHECKS_PASS",status))
      "PREPARATION_NUMERICAL_CHECKS_PASS" else status,
    substantive_imputation_assumptions_proven=FALSE,
    checks_total=nrow(review),checks_passed=sum(review$pass),
    bundle_sha256=if(file.exists(bundle_path))sha(bundle_path)else NA_character_,
    mids_sha256=if(file.exists(mids_path))sha(mids_path)else NA_character_,
    review_script_sha256=sha("R/16c_verify_charls_stage2_bundle.R"),
    review_csv_sha256=sha(file.path(qc,"model_bundle_review.csv")),
    updated_utc=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")),
    file.path(qc,"model_bundle_review.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
  lines <- c("# CHARLS 第二研究环节：共同插补准备独立核验", "",
    paste0("状态：", status), paste0("时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
    paste0("实际检查：", sum(review$pass), "/", nrow(review), " 项通过。"),
    if (nzchar(reason)) paste0("停止原因：", reason) else "",
    "本核验仅读取既有数据、MICE对象、签名和旧模型身份字段；未拟合任何模型，也未提取或解释疾病关联效应。", "",
    "核对共同12036人、53786条记录、原观察值/缺失标记、20份完成数据、九题总分与固定时长样条、每规格满秩及签名。参与者为独立单位，面板记录不是独立样本量。", "",
    "MICE链的有限性和数值摘要只属于数值卫生检查，不证明链已充分混合、MAR成立或与多状态似然相容。loggedEvents非空、支持范围/观察值保留/矩阵秩失败均需停止复核，不能自动改插补或模型规格。",
    "一般MICE后派生自然样条不是严格同源插补；共同数据中的模型变化也不是因果中介效应。", "",
    "nature-statistics 用于缺失处理、人/记录分母、共同比较样本及推断边界核验。",
    paste0("R ", getRversion(), "; mice ", packageVersion("mice"), "; data.table ", packageVersion("data.table"), "."))
  writeLines(lines, file.path(qc, "model_bundle_review.md"), useBytes = TRUE)
  message(status, ": ", sum(review$pass), "/", nrow(review), " checks passed.")
}

run_review <- function() {
  needed <- c(bundle_path, mids_path, old_bundle_path, old_mids_path, observed_path, contract_path, control_path)
  missing <- needed[!file.exists(needed)]
  if (length(missing)) {
    add("all_preparation_inputs_exist", FALSE, paste(missing, collapse = " | "))
    publish("NOT_READY_NO_MODELS_AUTHORIZED", "Preparation inputs are missing.")
    return(2L)
  }
  tryCatch({
    initial_hashes <- setNames(vapply(needed, sha, ""), needed)
    bun <- readRDS(bundle_path); old <- readRDS(old_bundle_path)
    imp <- readRDS(mids_path); old_imp <- readRDS(old_mids_path)
    control <- fromJSON(control_path, simplifyVector = TRUE)
    guard("current_stage_authorization_and_stop", control$research_stage == 2L &&
      isTRUE(control$stage_enabled) && isTRUE(control$stop_after_stage) &&
      identical(control$automatic_next_stage, FALSE) && identical(control$allow_later_stages, FALSE) &&
      isTRUE(control$requires_results_and_figure_review))
    models <- c("CORE", "BOTH")
    guard("control_two_models_twenty_draws_each", identical(as.character(control$models), models) &&
      identical(control$primary_extended_model, "BOTH") && control$required_imputations_per_model == 20L &&
      control$required_model_fits == 40L && control$expected_people == 12036L && control$expected_panel_records == 53786L)
    guard("approved_preparation_status", identical(bun$status, "APPROVED_FOR_STAGE2_MODELS") &&
      !is.null(bun$preparation_checks) && nrow(bun$preparation_checks) > 0L && all(bun$preparation_checks$pass))
    guard("primary_bundle_approved", identical(old$status, "APPROVED_FOR_PILOT_AND_MI") && all(old$preparation_checks$pass))

    hashes <- unlist(bun$input_hashes, use.names = TRUE)
    guard("signed_inputs_exist_and_named", length(hashes) > 0L && !is.null(names(hashes)) &&
      !anyDuplicated(names(hashes)) && all(nzchar(names(hashes))) && all(file.exists(names(hashes))))
    guard("signed_inputs_hashes_current", identical(unname(vapply(names(hashes), sha, "")), unname(as.character(hashes))))
    required_signed <- c(old_bundle_path, old_mids_path, observed_path, contract_path, control_path)
    guard("critical_preparation_sources_signed", all(vapply(required_signed, canonical, "") %in%
      vapply(names(hashes), canonical, "")))
    guard("primary_signed_inputs_current", all(file.exists(names(old$input_hashes))) &&
      identical(unname(vapply(names(old$input_hashes), sha, "")), unname(as.character(old$input_hashes))))
    original_signature <- digest(list(bundle = sha(old_bundle_path), script = sha("R/14b_fit_charls_disease_models.R"),
      msm = as.character(packageVersion("msm"))), algo = "sha256")
    guard("primary_signature_matches_control", identical(original_signature, control$primary_stage1_signature) &&
      identical(bun$primary_stage1_signature, control$primary_stage1_signature))
    guard("bundle_stage_and_model_counts_match_control", bun$research_stage == 2L && bun$m == 20L &&
      bun$seed == 20260914L && identical(as.character(bun$models), models))

    warm_paths <- file.path(old_out, "models", sprintf("MI_%02d.rds", 1:20))
    warm_hashes <- if (!is.null(bun$warm_start_hashes)) unlist(bun$warm_start_hashes, use.names = TRUE) else hashes
    guard("all_twenty_warm_sources_exist", all(file.exists(warm_paths)))
    guard("warm_hashes_are_named_paths", length(warm_hashes) > 0L && !is.null(names(warm_hashes)) &&
      all(file.exists(names(warm_hashes))) && !anyDuplicated(vapply(names(warm_hashes), canonical, "")))
    warm_index <- match(vapply(warm_paths, canonical, ""), vapply(names(warm_hashes), canonical, ""))
    guard("all_twenty_warm_source_hashes_frozen", !anyNA(warm_index) && identical(unname(vapply(warm_paths, sha, "")),
      unname(as.character(warm_hashes[warm_index]))))
    for (i in 1:20) {
      saved <- readRDS(warm_paths[i])
      guard(sprintf("warm_MI_%02d_identity_only", i), identical(saved$status, "ACCEPTED") &&
        identical(saved$signature, original_signature) && as.integer(saved$imputation) == i &&
        identical(saved$model, sprintf("MI_%02d", i)) && saved$n_people == 12036L && saved$n_records == 53786L,
        "Only identity/status/sample fields checked; no effect extraction.")
      rm(saved)
    }

    same_object <- function(a,b) {
      if(is.data.frame(a)&&is.data.frame(b))return(identical(as.data.frame(a),as.data.frame(b)))
      identical(a,b)
    }
    guard("same_primary_panel_values_types_order", same_object(bun$panel, old$panel) && nrow(bun$panel) == 53786L)
    for (field in c("edges", "q_initial", "core_covariates", "censor_codes", "censor_states"))
      guard(paste0("same_primary_", field), same_object(bun[[field]], old[[field]]))
    b <- as.data.table(bun$baseline); ob <- as.data.table(old$baseline)
    ids <- as.character(b$id); core <- as.character(old$core_covariates)
    guard("same_12036_baseline_ids_in_order", nrow(b) == 12036L && !anyDuplicated(ids) &&
      identical(ids, as.character(ob$id)) && setequal(ids, as.character(bun$panel$id)))
    guard("baseline_original_exposure_core_unchanged", all(c("id", "restless_bl", core) %in% names(b)) &&
      all(vapply(c("restless_bl", core), function(v) same_values(b[[v]], ob[[v]]), TRUE)))
    guard("fixed_exposure_complete_0_to_3", all(b$restless_bl %in% 0:3))
    guard("twenty_new_completed_datasets", length(bun$imputations) == 20L && inherits(imp, "mids") && imp$m == 20L && imp$iteration == 20L)
    guard("new_mice_seed", same_number(imp$seed, 20260914, 0))

    raw_fields <- c("botherl", "mindtsl", "depresl", "effortl", "fhopel", "fearll", "whappyl", "flonel", "goingl")
    item_fields <- paste0("cesd_", raw_fields)
    guard("nine_item_names_no_sleep_item", identical(as.character(bun$raw_item_names), item_fields) &&
      !("cesd_sleeprl" %in% item_fields))
    md <- as.data.frame(bun$original_mi_data)
    old_md <- as.data.frame(old_imp$data)
    expected_fields <- c(names(old_md), item_fields, "night_hours")
    guard("original_MI_data_has_only_prespecified_variables", nrow(md) == 12036L && identical(names(md), expected_fields) &&
      !any(c("id", "subject_id", "cesd9_5", "night_ns1", "night_ns2", "night_ns3") %in% names(md)))
    guard("core_original_missingness_and_outcome_auxiliaries_preserved", all(vapply(names(old_md), function(v)
      same_values(md[[v]], old_md[[v]]), TRUE)), "Reuses original incomplete data, not old completed draws.")
    guard("mids_data_matches_frozen_original_MI_data", identical(names(imp$data), names(md)) &&
      all(vapply(names(md), function(v) same_values(imp$data[[v]], md[[v]]), TRUE)))
    guard("baseline_ID_to_original_MI_data_order", all(vapply(c("restless_bl", core), function(v)
      same_values(md[[v]], b[[v]]), TRUE)))

    obs <- as.data.table(readRDS(observed_path)); w1 <- obs[wave == 1L]
    guard("unique_W1_observed_source_ids", !anyDuplicated(w1$id))
    w1 <- w1[match(ids, as.character(w1$id))]
    guard("W1_sources_match_every_baseline_id", identical(as.character(w1$id), ids))
    expected_items <- as.data.frame(lapply(raw_fields, function(v) {
      a <- as.numeric(w1[[v]]); a[!a %in% 1:4] <- NA_real_
      if (v %in% c("fhopel", "whappyl")) 4 - a else a - 1
    }))
    names(expected_items) <- item_fields
    guard("all_nine_items_scored_once_from_observed_source", all(vapply(item_fields, function(v)
      same_number(md[[v]], expected_items[[v]], 0), TRUE)))
    guard("night_hours_same_as_valid_observed_source", same_number(md$night_hours, w1$sleep_hours_observed, 0) &&
      all(md$night_hours[!is.na(md$night_hours)] >= 1 & md$night_hours[!is.na(md$night_hours)] <= 18))
    observed_sum <- rowSums(as.matrix(expected_items))
    guard("observed_CESD9_coverage_and_sum", sum(!is.na(observed_sum)) == 11439L &&
      same_number(observed_sum, w1$cesd9_no_sleep, 0) && sum(!is.na(md$night_hours)) == 11922L &&
      sum(!is.na(observed_sum) & !is.na(md$night_hours)) == 11362L)
    guard("baseline_observed_CESD9_derived", "cesd9_5" %in% names(b) && same_number(b$cesd9_5, observed_sum / 5, 1e-12))

    sp <- bun$spline_spec
    observed_knots <- as.numeric(quantile(md$night_hours, probs = c(1/3, 2/3), na.rm = TRUE, type = 7L))
    observed_boundary <- range(md$night_hours, na.rm = TRUE)
    guard("frozen_natural_spline_specification", same_number(sp$knots, c(6, 7), 0) &&
      same_number(sp$boundary, c(1, 13), 0) && identical(as.integer(sp$df), 3L) &&
      identical(as.integer(sp$quantile_type), 7L) && same_number(sp$knots, observed_knots, 0) &&
      same_number(sp$boundary, observed_boundary, 0) && all(diff(c(sp$boundary[1], sp$knots, sp$boundary[2])) > 0))
    spline <- function(x) as.matrix(splines::ns(x, knots = sp$knots, Boundary.knots = sp$boundary, intercept = FALSE))
    baseline_ns <- spline(md$night_hours)
    guard("baseline_observed_night_spline_derived", all(paste0("night_ns", 1:3) %in% names(b)) &&
      all(vapply(1:3, function(j) same_number(b[[paste0("night_ns", j)]], baseline_ns[, j], 1e-12), TRUE)))

    expected_covs <- list(CORE = core, CESD9 = c(core, "cesd9_5"), DURATION = c(core, paste0("night_ns", 1:3)),
      BOTH = c(core, "cesd9_5", paste0("night_ns", 1:3)))
    guard("two_model_covariate_sets", identical(names(bun$model_covariates), models) &&
      all(vapply(models, function(v) identical(as.character(bun$model_covariates[[v]]), expected_covs[[v]]), TRUE)))
    expected_methods <- setNames(rep("", ncol(md)), names(md))
    binary <- c("female", "married", "rural", "smoker", "drinker", "diabetes", "lung", "pain")
    for (v in c(core, item_fields, "night_hours")) if (anyNA(md[[v]]))
      expected_methods[v] <- if (v == "education4") "polyreg" else if (v %in% binary) "logreg" else "pmm"
    guard("mice_methods_exactly_prespecified", identical(names(imp$method), names(expected_methods)) &&
      identical(as.character(imp$method), as.character(expected_methods)))
    targets <- names(expected_methods)[nzchar(expected_methods)]
    guard("imputation_candidate_variables_exactly_prespecified",
      setequal(as.character(bun$imputation_target_variables), c(core, item_fields, "night_hours")) &&
      all(targets %in% bun$imputation_target_variables) &&
      !any(c("restless_bl", setdiff(names(old_md), core)) %in% targets),
      "Candidate variables may be fully observed; actual imputed variables are exactly the nonempty methods.")
    predictor <- matrix(1L, ncol(md), ncol(md), dimnames = list(names(md), names(md)))
    diag(predictor) <- 0L; predictor[expected_methods == "", ] <- 0L
    guard("actual_predictor_matrix_exactly_prespecified", identical(dimnames(imp$predictorMatrix), dimnames(predictor)) &&
      same_number(imp$predictorMatrix, predictor, 0))
    guard("actual_missingness_where_preserved", identical(dim(imp$where), dim(md)) &&
      identical(unname(imp$where), unname(is.na(as.matrix(md)))))
    guard("mice_logged_events_empty", is.null(imp$loggedEvents) || NROW(imp$loggedEvents) == 0L,
      "Any automatic collinearity/removal or fitting event requires explicit review.")

    model_rank_rows <- list(); support_rows <- list()
    for (i in 1:20) {
      label <- sprintf("MI_%02d", i)
      complete_data <- as.data.frame(mice::complete(imp, action = i))
      draw <- as.data.table(bun$imputations[[i]])
      guard(paste0(label, "_complete_data_fields"), identical(names(complete_data), names(md)) && nrow(complete_data) == 12036L &&
        !anyNA(complete_data))
      guard(paste0(label, "_every_observed_value_preserved"), all(vapply(names(md), function(v) {
        keep <- !is.na(md[[v]])
        same_values(complete_data[[v]][keep], md[[v]][keep])
      }, TRUE)), "Includes exposure, original core observations, nine items, duration and all fixed outcome auxiliaries.")
      for (v in targets) {
        values <- complete_data[[v]]; donors <- md[[v]][!is.na(md[[v]])]
        if (is.factor(values)) values <- as.character(values)
        if (is.factor(donors)) donors <- as.character(donors)
        support_ok <- if (expected_methods[v] == "pmm") all(values[is.na(md[[v]])] %in% donors) else
          if (v == "education4") all(as.character(values) %in% as.character(1:4)) else all(as.character(values) %in% c("0", "1"))
        support_rows[[length(support_rows) + 1L]] <- data.table(imputation = i, variable = v,
          method = expected_methods[v], missing_n = sum(is.na(md[[v]])), support_pass = support_ok)
        guard(paste0(label, "_support_", v), support_ok)
      }
      guard(paste0(label, "_same_IDs_and_complete_model_variables"), nrow(draw) == 12036L && !anyDuplicated(draw$id) &&
        identical(as.character(draw$id), ids) && all(c("restless_bl", expected_covs$BOTH) %in% names(draw)) &&
        !anyNA(draw[, c("restless_bl", expected_covs$BOTH), with = FALSE]))
      guard(paste0(label, "_core_and_exposure_same_shared_draw"), all(vapply(c("restless_bl", core), function(v)
        same_values(draw[[v]], complete_data[[v]]), TRUE)))
      guard(paste0(label, "_education_levels_fixed"), is.factor(draw$education4) && identical(levels(draw$education4), as.character(1:4)))
      completed_items <- as.matrix(complete_data[, item_fields, drop = FALSE])
      cesd <- rowSums(completed_items)
      ns <- spline(complete_data$night_hours)
      guard(paste0(label, "_CESD9_reconstructed_from_all_nine_items"), all(completed_items %in% 0:3) &&
        all(cesd %in% 0:27) && same_number(draw$cesd9_5, cesd / 5, 1e-12))
      guard(paste0(label, "_night_spline_reconstructed_fixed_knots"), identical(dim(ns), c(12036L, 3L)) &&
        all(is.finite(ns)) && all(vapply(1:3, function(j)
          same_number(draw[[paste0("night_ns", j)]], ns[, j], 1e-12), TRUE)))
      for (model in models) {
        mm <- model.matrix(reformulate(c("restless_bl", expected_covs[[model]])), data = as.data.frame(draw))
        rank <- qr(mm)$rank
        full <- nrow(mm) == 12036L && all(is.finite(mm)) && rank == ncol(mm)
        model_rank_rows[[length(model_rank_rows) + 1L]] <- data.table(imputation = i, model = model,
          people = nrow(mm), columns = ncol(mm), matrix_rank = rank, full_rank = full)
        guard(paste0(label, "_", model, "_full_rank"), full)
      }
    }
    dir.create(qc, recursive = TRUE, showWarnings = FALSE)
    fwrite(rbindlist(support_rows), file.path(qc, "independent_imputation_support.csv"), bom = TRUE)
    fwrite(rbindlist(model_rank_rows), file.path(qc, "independent_model_matrix_ranks.csv"), bom = TRUE)

    guard("mice_chain_arrays_available", is.array(imp$chainMean) && is.array(imp$chainVar) &&
      length(dim(imp$chainMean)) == 3L && identical(dim(imp$chainMean), dim(imp$chainVar)) &&
      all(targets %in% dimnames(imp$chainMean)[[1]]) && all(dim(imp$chainMean)[2:3] == 20L))
    chain_rows <- lapply(targets, function(v) {
      means <- imp$chainMean[v, , ]; variances <- imp$chainVar[v, , ]
      missing_n <- sum(is.na(md[[v]]))
      singleton <- missing_n==1L
      good <- all(is.finite(means)) && if(singleton)all(is.na(variances)) else
        all(is.finite(variances)) && all(variances >= 0)
      add(paste0("chain_numeric_hygiene_", v), good, "Finite traces do not establish convergence or adequate mixing.")
      data.table(variable = v, missing_n=missing_n,
        variance_status=if(singleton)"VAR_UNDEFINED_SINGLE_MISSING"else"FINITE_NONNEGATIVE_REQUIRED",
        iterations = nrow(means), chains = ncol(means), finite_nonnegative = good,
        chain_mean_min = min(means), chain_mean_max = max(means),
        final_iteration_mean_across_chains = mean(means[nrow(means), ]),
        final_iteration_sd_across_chains = sd(means[nrow(means), ]),
        chain_variance_min = if(singleton)NA_real_ else min(variances),
        chain_variance_max = if(singleton)NA_real_ else max(variances))
    })
    chain_summary <- rbindlist(chain_rows)
    fwrite(chain_summary, file.path(qc, "independent_chain_numeric_summary.csv"), bom = TRUE)
    guard("all_chain_numeric_hygiene_checks", all(chain_summary$finite_nonnegative),
      "This is not a convergence verdict; substantive imputation review remains required.")
    guard("all_reviewed_preparation_inputs_unchanged", identical(unname(vapply(needed, sha, "")), unname(initial_hashes)))
    guard("all_frozen_warm_models_unchanged", identical(unname(vapply(warm_paths, sha, "")), unname(as.character(warm_hashes[warm_index]))))
    guard("all_signed_sources_unchanged_at_review_end", identical(unname(vapply(names(hashes), sha, "")), unname(as.character(hashes))))
    publish("PREPARATION_NUMERICAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_IMPUTATION_REVIEW")
    0L
  }, error = function(e) {
    add("independent_preparation_review_completed", FALSE, conditionMessage(e))
    publish("FAILED_REVIEW_REQUIRED", conditionMessage(e))
    1L
  })
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments)) stop("This preparation-only reviewer takes no arguments.", call. = FALSE)
quit(save = "no", status = run_review())
