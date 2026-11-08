#!/usr/bin/env Rscript
# CHARLS stage 3 preparation: W1-W2 dynamic restless-sleep pattern, W2 landmark,
# estimand-matched MI20. No disease-transition association model is fitted here.
options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.treatment"), warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(mice)
  library(digest)
  library(jsonlite)
})
setDTthreads(1L)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, args[[1]] == "prepare")
out <- "outputs/charls_dynamic_sleep_20260920"
old <- "outputs/charls_external_extension_20260913"
stage2 <- "outputs/charls_cesd9_duration_20260914"
for (folder in c("data", "tables", "model_qc", "results", "reports", "logs"))
  dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)

bundle_path <- file.path(out, "data/charls_dynamic_model_bundle.rds")
mids_path <- file.path(out, "data/charls_dynamic_mice20.rds")
pre_path <- file.path(out, "data/charls_dynamic_preimputation.rds")
if (any(file.exists(c(bundle_path, mids_path, pre_path))))
  stop("Existing stage-3 preparation retained; refusing to overwrite.")

control_path <- "analysis/charls_stage3_dynamic_sleep_control.json"
contract_path <- "analysis/charls_stage3_dynamic_sleep_contract.md"
stage2_bundle_path <- file.path(stage2, "data/charls_stage2_model_bundle.rds")
stage2_mids_path <- file.path(stage2, "data/charls_stage2_mice20.rds")
observed_path <- file.path(old, "data/charls_observed_exposures_tristate_panel.rds")
conflict_path <- file.path(old, "secondary_readiness/internal/living_status_conflict_locators.rds")
readiness_path <- file.path(old, "secondary_readiness/W2_landmark_flow.csv")
sources <- c(control_path, contract_path, stage2_bundle_path, stage2_mids_path,
  observed_path, conflict_path, readiness_path,
  "R/13a_prepare_charls_external_extension.R", "R/14a_prepare_charls_model_bundle.R",
  "R/16a_prepare_charls_cesd9_duration.R", "R/17a_prepare_charls_dynamic_sleep.R")

sha <- function(p) digest(file = p, algo = "sha256")
utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
write_json_safe <- function(x, p) write_json(x, p, pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null")
save_new <- function(x, p) {
  if (file.exists(p)) stop("Refusing to overwrite ", p)
  tmp <- tempfile("charls_dynamic_", tmpdir = dirname(p), fileext = ".rds")
  saveRDS(x, tmp, compress = "gzip")
  if (!file.rename(tmp, p)) stop("Could not atomically publish ", p, "; retained ", tmp)
}
checks <- list()
add <- function(name, pass, detail = "") {
  checks[[length(checks) + 1L]] <<- data.table(check = name, pass = isTRUE(pass), detail = as.character(detail))
  if (!isTRUE(pass)) stop(paste(name, detail, sep = ": "), call. = FALSE)
}
status <- function(value, detail = "") {
  write_json_safe(list(status = value, stage = 3L, pid = Sys.getpid(), updated_utc = utc(),
    detail = detail, model_fits_started = FALSE, automatic_next_stage = FALSE),
    file.path(out, "model_qc/preparation_status.json"))
}
same_observed <- function(a, b) {
  keep <- !is.na(a)
  identical(as.character(a[keep]), as.character(b[keep]))
}

run <- function() {
  status("AUDITING_W2_LANDMARK_INPUTS")
  add("all_sources_exist", all(file.exists(sources)), paste(sources[!file.exists(sources)], collapse = " | "))
  hashes_before <- setNames(vapply(sources, sha, ""), sources)
  control <- read_json(control_path, simplifyVector = TRUE)
  add("stage3_scope_check", control$research_stage == 3L && isTRUE(control$stage_enabled) &&
    isTRUE(control$stop_after_stage) && identical(control$automatic_next_stage, FALSE) &&
    identical(control$allow_later_stages, FALSE))
  add("two_models_and_forty_fits", identical(as.character(control$models), c("CORE", "BOTH")) &&
    identical(control$primary_model, "CORE") && control$required_imputations_per_model == 20L &&
    control$required_model_fits == 40L)

  bun <- readRDS(stage2_bundle_path)
  old_imp <- readRDS(stage2_mids_path)
  add("approved_stage2_source", identical(bun$status, "APPROVED_FOR_STAGE2_MODELS") &&
    all(bun$preparation_checks$pass) && length(bun$imputations) == 20L &&
    inherits(old_imp, "mids") && old_imp$m == 20L)
  stage2_signed <- unlist(bun$input_hashes, use.names = TRUE)
  add("stage2_signed_sources_current_including_teaching_code",
    length(stage2_signed) > 0L && all(file.exists(names(stage2_signed))) &&
      identical(unname(vapply(names(stage2_signed), sha, "")), unname(as.character(stage2_signed))))
  add("frozen_W1_source_population", nrow(bun$baseline) == 12036L && nrow(bun$panel) == 53786L &&
    !anyDuplicated(bun$baseline$id) && setequal(bun$baseline$id, bun$panel$id))

  obs <- as.data.table(readRDS(observed_path))
  panel0 <- as.data.table(copy(bun$panel))
  base_ids <- as.character(bun$baseline$id)
  sleep <- obs[id %in% base_ids & wave %in% 1:2,
    .(id = as.character(id), wave = as.integer(wave), score = as.numeric(restless_score_observed))]
  add("unique_observed_sleep_ID_wave", !anyDuplicated(sleep[, .(id, wave)]))
  sleep <- dcast(sleep, id ~ wave, value.var = "score")
  setnames(sleep, c("1", "2"), c("w1_score", "w2_score"))
  sleep[, `:=`(w1_high = as.integer(w1_score >= 2), w2_high = as.integer(w2_score >= 2))]
  sleep[, pattern := factor(c("persistent_low", "onset", "remission", "persistent_high")
    [1L + 2L * w1_high + w2_high],
    levels = c("persistent_low", "onset", "remission", "persistent_high"))]
  pattern_candidates <- sleep[!is.na(w1_score) & !is.na(w2_score)]
  add("pattern_uses_two_observed_waves_only", nrow(pattern_candidates) == 9958L &&
    !anyNA(pattern_candidates[, .(w1_score, w2_score, pattern)]))
  add("pattern_threshold_and_mapping", all(pattern_candidates$w1_high == as.integer(pattern_candidates$w1_score >= 2)) &&
    all(pattern_candidates$w2_high == as.integer(pattern_candidates$w2_score >= 2)))

  conflicts <- as.data.table(readRDS(conflict_path))
  add("two_known_W4_status_conflicts", nrow(conflicts) == 2L && all(conflicts$wave == 4L) &&
    all(conflicts$conflicting_start_intervals == 0L) && all(conflicts$conflicting_end_intervals == 1L))
  panel <- panel0[!conflicts[, .(id, wave)], on = .(id, wave)]
  add("only_two_disputed_rows_removed", nrow(panel0) - nrow(panel) == 2L)

  w2_living <- panel[wave == 2L & record_type == "living_visit", .(id, w2_time = time, w2_state = state)]
  add("unique_W2_living_row", !anyDuplicated(w2_living$id))
  landmark <- merge(pattern_candidates, w2_living[w2_state %in% 1:8], by = "id")
  post <- merge(panel, landmark[, .(id, pattern, w1_score, w2_score, w2_time, w2_state)], by = "id")
  post <- post[time >= w2_time]
  post[, source_time_from_W1 := time]
  post[, time := time - w2_time]
  eligible <- post[, .(records = .N, positive_followup = any(time > 0)), by = id][records >= 2L & positive_followup, id]
  post <- post[id %in% eligible]
  landmark <- landmark[id %in% eligible]
  setorder(landmark, id)
  setorder(post, id, time)
  add("frozen_dynamic_people", nrow(landmark) == control$expected_people && uniqueN(post$id) == control$expected_people)
  add("frozen_dynamic_records", nrow(post) == control$expected_panel_records)
  add("one_known_W2_entry_each", nrow(post[time == 0]) == nrow(landmark) &&
    all(post[time == 0, state] %in% 1:8) && all(post[, min(time), by = id]$V1 == 0))
  add("strict_positive_post_entry_times", all(post[time > 0, time] > 0))
  add("death_and_partial_counts", sum(post$state == 9L) == control$expected_interval_deaths &&
    sum(post$state >= 100L) == control$expected_partial_records)
  end <- post[, .(death_n = sum(state == 9L), last_state = last(state)), by = id]
  add("death_once_and_terminal", all(end$death_n <= 1L & (end$death_n == 0L | end$last_state == 9L)))
  counts <- landmark[, .N, by = pattern][match(control$pattern_levels, as.character(pattern))]
  add("frozen_pattern_counts", identical(as.integer(counts$N), as.integer(control$pattern_counts)))
  add("no_pattern_or_initial_state_imputation", !anyNA(landmark[, .(pattern, w2_state)]))

  model_panel <- post[, .(id, wave, time, state, obstype, record_type)]
  used_censor <- bun$censor_codes %in% unique(model_panel$state)
  censor_codes <- bun$censor_codes[used_censor]
  censor_states <- bun$censor_states[used_censor]
  add("partial_state_sets_retained", length(censor_codes) > 0L &&
    all(model_panel[state >= 100L, state] %in% censor_codes) && all(unlist(censor_states) %in% 1:8))

  flow <- data.table(step = c("W1 main analysis population", "Reliable W2 living interview",
    "Observed W1 and W2 restless sleep", "Plus fully known W2 disease state",
    "Plus reliable information after W2; final landmark sample"),
    n = c(length(base_ids), uniqueN(panel[wave == 2L & record_type == "living_visit", id]),
      nrow(pattern_candidates), nrow(merge(pattern_candidates, w2_living[w2_state %in% 1:8], by = "id")), nrow(landmark)))
  fwrite(flow, file.path(out, "tables/dynamic_landmark_flow.csv"), bom = TRUE)
  pattern_table <- landmark[, .(people = .N, percent = 100 * .N / nrow(landmark)), by = pattern]
  pattern_table <- pattern_table[match(control$pattern_levels, as.character(pattern))]
  fwrite(pattern_table, file.path(out, "tables/dynamic_pattern_distribution.csv"), bom = TRUE)
  fwrite(landmark[, .N, by = .(pattern, w2_state)][order(w2_state, pattern)],
    file.path(out, "tables/W2_origin_state_by_pattern.csv"), bom = TRUE)

  known <- post[state %in% 1:8]
  setorder(known, id, time)
  known[, next_state := shift(state, type = "lead"), by = id]
  known[, next_time := shift(time, type = "lead"), by = id]
  pairs <- known[!is.na(next_state) & next_time > time,
    .(id, pattern, from = state - 1L, to = next_state - 1L)]
  targets <- data.table(edge_order = 1:4, from = c(0L, 0L, 1L, 1L), to = c(2L, 4L, 3L, 5L),
    transition = c("0-2", "0-4", "1-3", "1-5"))
  support <- pairs[, .(observed_known_pairs = .N, people = uniqueN(id)), by = .(pattern, from, to)]
  support <- merge(merge(CJ(pattern = control$pattern_levels, edge_order = 1:4), targets, by = "edge_order"),
    support, by = c("pattern", "from", "to"), all.x = TRUE)
  support[is.na(observed_known_pairs), `:=`(observed_known_pairs = 0L, people = 0L)]
  setorder(support, edge_order, pattern)
  add("all_pattern_target_cells_have_observed_support", all(support$observed_known_pairs > 0L))
  fwrite(support, file.path(out, "tables/four_target_observed_known_pairs_by_pattern.csv"), bom = TRUE)

  # Estimand-matched MICE: preserve the original incomplete W1 variables, subset
  # to the frozen landmark cohort, and add the dynamic exposure and W2 state as
  # fixed predictors. Disease/exposure/time fields themselves are not imputed.
  mi <- as.data.table(copy(bun$original_mi_data))
  mi[, id := base_ids]
  mi <- mi[match(landmark$id, id)]
  add("MI_rows_match_landmark_ID_order", identical(mi$id, landmark$id))
  mi[, id := NULL]
  mi[, pattern := factor(landmark$pattern, levels = control$pattern_levels)]
  mi[, w2_state_factor := factor(landmark$w2_state, levels = 1:8)]
  mi[, w2_score := as.numeric(landmark$w2_score)]
  mi[, w2_landmark_time := as.numeric(landmark$w2_time)]
  add("dynamic_fixed_predictors_complete", !anyNA(mi[, .(pattern, w2_state_factor, w2_score, w2_landmark_time)]))

  core <- as.character(bun$core_covariates)
  item_names <- as.character(bun$raw_item_names)
  extra <- c(item_names, "night_hours")
  targets_mi <- c(core, extra)
  methods <- setNames(rep("", ncol(mi)), names(mi))
  for (v in core) if (anyNA(mi[[v]])) methods[v] <- old_imp$method[[v]]
  for (v in extra) if (anyNA(mi[[v]])) methods[v] <- "pmm"
  predictor <- matrix(1L, ncol(mi), ncol(mi), dimnames = list(names(mi), names(mi)))
  diag(predictor) <- 0L
  predictor[methods == "", ] <- 0L
  fixed <- names(methods)[methods == ""]
  add("fixed_fields_complete", !anyNA(mi[, ..fixed]))
  add("dynamic_exposure_and_state_never_imputed", all(methods[c("restless_bl", "pattern", "w2_state_factor", "w2_score", "w2_landmark_time")] == ""))
  add("MI_targets_exact", setequal(names(methods)[nzchar(methods)], targets_mi[vapply(targets_mi, function(v) anyNA(mi[[v]]), TRUE)]))

  miss <- rbindlist(lapply(names(mi), function(v) data.table(variable = v, n = nrow(mi),
    missing = sum(is.na(mi[[v]])), method = methods[[v]],
    role = if (v %in% core) "core_covariate" else if (v %in% extra) "CESD_item_or_night_hours" else "fixed_exposure_or_outcome_auxiliary")))
  fwrite(miss, file.path(out, "tables/dynamic_imputation_missingness.csv"), bom = TRUE)
  save_new(list(data = as.data.frame(mi), ids = landmark$id, method = methods,
    predictor_matrix = predictor, raw_item_names = item_names,
    spline_spec = bun$spline_spec, pattern_levels = control$pattern_levels), pre_path)

  status("DYNAMIC_MICE20_RUNNING", "No transition model has started.")
  warnings_seen <- character()
  imp <- withCallingHandlers(mice(as.data.frame(mi), m = 20L, maxit = 20L, method = methods,
    predictorMatrix = predictor, seed = control$seed, printFlag = FALSE),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  save_new(imp, mids_path)
  events <- as.data.table(imp$loggedEvents)
  if (!ncol(events)) events <- data.table(it = integer(), im = integer(), dep = character(), meth = character(), out = character())
  fwrite(events, file.path(out, "model_qc/mice_logged_events.csv"), bom = TRUE)
  write_json_safe(list(warnings = as.list(unique(warnings_seen)), logged_events = nrow(events),
    m = imp$m, iterations = imp$iteration, seed = control$seed, methods = as.list(imp$method)),
    file.path(out, "model_qc/mice_runtime_review.json"))
  add("MICE_no_runtime_warnings", length(warnings_seen) == 0L, paste(unique(warnings_seen), collapse = "; "))
  add("MICE_no_logged_events", nrow(events) == 0L)
  add("MICE_20_by_20", imp$m == 20L && imp$iteration == 20L)

  knots <- as.numeric(bun$spline_spec$knots)
  boundary <- as.numeric(bun$spline_spec$boundary)
  add("frozen_duration_spline", identical(knots, c(6, 7)) && identical(boundary, c(1, 13)))
  derive <- function(z) {
    total <- rowSums(as.matrix(z[, item_names, drop = FALSE]))
    ns <- splines::ns(z$night_hours, knots = knots, Boundary.knots = boundary, intercept = FALSE)
    ans <- as.data.table(z[, c("pattern", core), drop = FALSE])
    ans[, pattern := factor(pattern, levels = control$pattern_levels)]
    for (v in setdiff(core, "education4")) if (is.factor(ans[[v]]))
      ans[, (v) := as.numeric(as.character(get(v)))]
    ans[, education4 := factor(education4, levels = 1:4)]
    ans[, cesd9_5 := as.numeric(total) / 5]
    for (j in 1:3) ans[, (paste0("night_ns", j)) := as.numeric(ns[, j])]
    ans[, id := landmark$id]
    setcolorder(ans, c("id", "pattern", core, "cesd9_5", paste0("night_ns", 1:3)))
    ans
  }
  specs <- list(CORE = core, BOTH = c(core, "cesd9_5", paste0("night_ns", 1:3)))
  draws <- vector("list", 20L)
  for (i in 1:20) {
    z <- complete(imp, i)
    add(sprintf("MI_%02d_complete", i), !anyNA(z))
    add(sprintf("MI_%02d_observed_values_preserved", i),
      all(vapply(names(mi), function(v) same_observed(mi[[v]], z[[v]]), TRUE)))
    add(sprintf("MI_%02d_fixed_fields_preserved", i),
      all(vapply(fixed, function(v) identical(as.character(mi[[v]]), as.character(z[[v]])), TRUE)))
    add(sprintf("MI_%02d_CESD_item_support", i), all(as.matrix(z[, item_names]) %in% 0:3))
    draws[[i]] <- derive(z)
    add(sprintf("MI_%02d_CESD9_range", i), all(draws[[i]]$cesd9_5 >= 0 & draws[[i]]$cesd9_5 <= 27 / 5))
    add(sprintf("MI_%02d_pattern_preserved", i), identical(as.character(draws[[i]]$pattern), as.character(landmark$pattern)))
    for (spec in names(specs)) {
      mm <- model.matrix(reformulate(c("pattern", specs[[spec]])), data = draws[[i]])
      add(sprintf("MI_%02d_%s_full_rank", i, spec), nrow(mm) == nrow(landmark) &&
        qr(mm)$rank == ncol(mm) && all(is.finite(mm)))
    }
  }
  observed_baseline <- derive(as.data.frame(mi))

  warm_paths <- unlist(lapply(c("CORE", "BOTH"), function(spec)
    file.path("private/model_store/CHARLS_CESD9_duration_20260914/models", spec, sprintf("MI_%02d.rds", 1:20))))
  add("all_40_stage2_warm_models_exist", all(file.exists(warm_paths)))
  warm_hashes <- setNames(vapply(warm_paths, sha, ""), warm_paths)
  input_hashes <- hashes_before
  add("all_read_only_sources_unchanged", identical(unname(vapply(sources, sha, "")), unname(hashes_before)))
  prep_checks <- rbindlist(checks)
  add("all_preparation_checks_pass", all(prep_checks$pass))
  prep_checks <- rbindlist(checks)

  bundle <- list(status = "APPROVED_FOR_DYNAMIC_MODELS", research_stage = 3L,
    panel = model_panel, baseline = observed_baseline, imputations = draws,
    landmark = landmark[, .(id, pattern, w1_score, w2_score, w2_time, w2_state)],
    edges = copy(bun$edges), q_initial = bun$q_initial, core_covariates = core,
    model_covariates = specs, pattern_levels = control$pattern_levels,
    pattern_coefficient_labels = c("patternonset", "patternremission", "patternpersistent_high"),
    censor_codes = censor_codes, censor_states = censor_states,
    raw_item_names = item_names, original_mi_data = as.data.frame(mi),
    imputation_target_variables = targets_mi, spline_spec = bun$spline_spec,
    input_hashes = input_hashes, warm_start_hashes = warm_hashes,
    preparation_checks = prep_checks, m = 20L, seed = control$seed,
    estimand = "W2-landmark conditional transition intensity by observed W1-W2 restless-sleep pattern",
    primary_dynamic_contrasts = control$primary_dynamic_contrasts,
    contextual_contrasts = control$contextual_contrasts,
    status_conflict_rows_excluded = nrow(conflicts))
  save_new(bundle, bundle_path)
  fwrite(prep_checks, file.path(out, "model_qc/preparation_checks.csv"), bom = TRUE)
  write_json_safe(list(status = "PREPARED_AWAITING_INDEPENDENT_REVIEW", people = nrow(landmark),
    panel_records = nrow(model_panel), interval_deaths = sum(model_panel$state == 9L),
    partial_records = sum(model_panel$state >= 100L), m = 20L, required_fits = 40L,
    bundle_sha256 = sha(bundle_path), mids_sha256 = sha(mids_path), model_fits_started = FALSE,
    automatic_next_stage = FALSE), file.path(out, "model_qc/preparation_status.json"))
  message("PREPARED_AWAITING_INDEPENDENT_REVIEW: n=", nrow(landmark), "; records=", nrow(model_panel),
    "; patterns=", paste(pattern_table$people, collapse = "/"), "; no transition model fitted.")
}

tryCatch(run(), error = function(e) {
  if (length(checks)) fwrite(rbindlist(checks, fill = TRUE), file.path(out, "model_qc/preparation_checks.csv"), bom = TRUE)
  status("FAILED_PREPARATION_REVIEW_REQUIRED", conditionMessage(e))
  message("FAILED_PREPARATION_REVIEW_REQUIRED: ", conditionMessage(e))
  quit(save = "no", status = 1L)
})
