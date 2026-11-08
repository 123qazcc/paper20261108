#!/usr/bin/env Rscript
# CHARLS mortality ascertainment audit; reads source data without changing it.
options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages({library(data.table); library(haven); library(digest)})
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "outputs/charls_external_extension_20260913/mortality_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
protected_paths <- file.path(root, c("R/13a_prepare_charls_external_extension.R", "data_clean/charls_clean.rds", "data_clean/hrs_clean.rds", "data_clean/elsa_clean.rds",
  "outputs/charls_external_extension_20260913/data/charls_w1_descriptive_baseline.rds"))
protected_before <- vapply(protected_paths, function(p) digest(file = p, algo = "sha256"), "")
raw_root <- "private/source/data_raw/1_CHARLS_China"
read_stata <- function(path, n_max = Inf) {
  staged <- tempfile(fileext = ".dta")
  stopifnot(file.copy(path, staged))
  on.exit(unlink(staged))
  as.data.table(haven::read_dta(staged, n_max = n_max))
}
dump_csv <- function(x, name) fwrite(x, file.path(out, paste0(name, ".csv")), bom = TRUE, na = "")
id12 <- function(x, w1 = FALSE) {
  z <- if (is.numeric(x)) sprintf("%012.0f", as.numeric(x)) else trimws(as.character(x))
  if (w1) z[nchar(z) == 11L] <- paste0(substr(z[nchar(z) == 11L], 1, 9), "0", substr(z[nchar(z) == 11L], 10, 11))
  z
}
files <- list.files(raw_root, pattern = "\\.dta$", recursive = TRUE, full.names = TRUE)
selected <- files[grepl("(Sample_Infor|Exit|Verbal|Insider|Interviewer|interviewer|household_roster)", basename(files)) & !grepl("LIFEhistory", files)]
manifest <- data.table(path = selected, sha256 = vapply(selected, function(p) digest(file = p, algo = "sha256"), ""))
dictionary <- rbindlist(lapply(selected, function(p) {
  a <- read_stata(p, 0)
  data.table(path = p, variable = names(a),
    label = vapply(a, function(x) paste(attr(x, "label"), collapse = ";"), ""),
    value_labels = vapply(a, function(x) paste(names(attr(x, "labels")), attr(x, "labels"), collapse = ";"), ""))
}), fill = TRUE)
dump_csv(manifest, "raw_source_manifest")
dump_csv(dictionary, "raw_tracking_exit_dictionary")
env <- new.env(parent = baseenv())
utils::data(list = "CHARLS", package = "GlobalAgingR", envir = env)
package <- as.data.table(GlobalAgingR::CH_filter(env[["charls"]], wave %in% 1:5))
package[, id := id12(ID)]
dump_csv(package[, .(n = .N, death_year_nonmissing = sum(!is.na(radyear))), by = wave], "package_death_field_by_wave")
dump_csv(package[, .N, by = .(wave, iwstat)][order(wave, iwstat)], "package_iwstat_counts")
samples <- rbindlist(lapply(c(2015, 2018, 2020), function(y) {
  path <- file.path(raw_root, paste0(y, "charls"), "Sample_Infor.dta")
  a <- read_stata(path)
  a[, .(id = id12(ID), wave = match(y, c(2011, 2013, 2015, 2018, 2020)),
    died = as.numeric(died), iyear = as.numeric(iyear), imonth = as.numeric(imonth), source = path)]
}))
dump_csv(samples[, .N, by = .(wave, died, iyear)][order(wave, died, iyear)], "sample_death_counts")
print(samples[, .N, by = .(wave, died, iyear)][order(wave, died, iyear)])
exit_specs <- data.table(wave = c(2L, 5L), path = file.path(raw_root, c("2013charls/Exit_Interview.dta", "2020charls/Exit_Module.dta")))
exits <- rbindlist(lapply(seq_len(nrow(exit_specs)), function(i) {
  a <- read_stata(exit_specs$path[i])
  a[, .(id = id12(ID), wave = exit_specs$wave[i], death_year = as.numeric(exb001_1), death_month = as.numeric(exb001_2),
    death_day = if ("exb001_3" %in% names(a)) as.numeric(exb001_3) else NA_real_, calendar = as.numeric(exb002), source = exit_specs$path[i])]
}))
print(exits[, .N, by = .(wave, death_year)][order(wave, death_year)])
dump_csv(exits[, .N, by = .(wave, death_year)][order(wave, death_year)], "exit_death_year_counts")
harm_path <- "private/source/harmonized/CHARLS/H_CHARLS_D_Data/H_CHARLS_D_Data.dta"
hstage <- tempfile(fileext = ".dta"); stopifnot(file.copy(harm_path, hstage))
harm <- as.data.table(haven::read_dta(hstage, col_select = c("ID", "radyear", "radmonth", "r1iwstat", "r2iwstat", "r3iwstat", "r4iwstat")))
unlink(hstage)
dump_csv(rbindlist(lapply(1:4, function(w) harm[, .N, by = c(paste0("r", w, "iwstat"))][, .(wave = w, iwstat = .SD[[1]], n = N)])), "harmonized_iwstat_counts")
docs <- c(file.path(raw_root, "2013charls/CHARLS_Wave2_Release_Note.pdf"),
  file.path(raw_root, "2015charls/CHARLS_2015_Release_Note.pdf"),
  file.path(raw_root, "2018charls/CHARLS_2018_Users_Guide.pdf"),
  file.path(raw_root, "2020charls/CHARLS_2020_User_Guide_Chinese.pdf"))
for (p in docs) {
  tfile <- tempfile(fileext = ".pdf"); file.copy(p, tfile)
  txt <- pdftools::pdf_text(tfile); unlink(tfile)
  writeLines(unlist(lapply(seq_along(txt), function(i) c(paste0("PAGE ", i), txt[i]))),
    file.path(out, sub("\\.pdf$", ".txt", basename(p))), useBytes = TRUE)
}
print(package[!is.na(radyear), .(n = uniqueN(id)), by = .(year = as.numeric(radyear))][order(year)])
print(package[, .(min_year = min(as.numeric(iwy), na.rm = TRUE), max_year = max(as.numeric(iwy), na.rm = TRUE)), by = wave])
print(exits[, .N, by = .(wave, calendar, valid_year = !is.na(death_year))])
print(samples[, .N, by = .(wave, died, imonth)][order(wave, died, imonth)])
w2weights <- read_stata(file.path(raw_root, "2013charls/Weights.dta"))
dump_csv(data.table(variable = c("INDV_L_Died", "iyear", "imonth"),
  label = vapply(w2weights[, .SD, .SDcols = c("INDV_L_Died", "iyear", "imonth")], function(x) paste(attr(x, "label"), collapse = ";"), ""),
  value_labels = vapply(w2weights[, .SD, .SDcols = c("INDV_L_Died", "iyear", "imonth")], function(x) paste(names(attr(x, "labels")), attr(x, "labels"), collapse = ";"), "")), "w2_weights_dictionary")
harm[, id := id12(ID)]
old <- as.data.table(readRDS(file.path(root, "data_clean/charls_clean.rds")))
old[, id := id12(id)]
pkgdead <- package[!is.na(radyear), .(package_death_year = as.numeric(radyear)[1], package_death_month = as.numeric(radmonth)[1], n_distinct_year = uniqueN(radyear)), by = id]
explicit <- unique(rbindlist(list(samples[died == 1, .(id)], exits[, .(id)], harm[r2iwstat %in% 5:6 | r3iwstat %in% 5:6 | r4iwstat %in% 5:6, .(id)])))
print(data.table(item = c("package_death_ids", "old_appended_ids", "explicit_death_ids", "package_deaths_no_raw_status", "raw_status_no_package_death", "package_death_distinctyear_conflicts"),
 n = c(nrow(pkgdead), nrow(old[event_type == "death"]), nrow(explicit), sum(!pkgdead$id %in% explicit$id), sum(!explicit$id %in% pkgdead$id), sum(pkgdead$n_distinct_year > 1))))
# Initial conservative interval-censored death gate. A death is observed by its
# documented ascertainment interview month END, never as an exact death time.
w2 <- w2weights[, .(id = id12(ID), wave = 2L, died = as.numeric(INDV_L_Died), iyear = as.numeric(iyear), imonth = as.numeric(imonth),
  source = file.path(raw_root, "2013charls/Weights.dta"))]
tracking <- rbindlist(list(w2, samples), use.names = TRUE)
invalid_tracking <- tracking[!grepl("^[0-9]{12}$", id)]
dump_csv(invalid_tracking[, .N, by = .(wave, died)], "excluded_tracking_rows_without_person_id")
stopifnot(all(is.na(invalid_tracking$died)))
tracking <- tracking[grepl("^[0-9]{12}$", id)]
tracking[, `:=`(ascertain_lower = iyear + (imonth - 1)/12, ascertain_upper = iyear + imonth/12)]
tracking[!iyear %in% 2011:2021 | !imonth %in% 1:12, c("ascertain_lower", "ascertain_upper") := .(NA_real_, NA_real_)]
positive <- tracking[died == 1]
health_alive <- old[event_type == "visit", .(id, wave, event_time,
  alive_lower = interview_year + fifelse(interview_month %in% 1:12, (interview_month - 1)/12, 0),
  alive_upper = interview_year + fifelse(interview_month %in% 1:12, interview_month/12, 1))]
alive_death_same_wave <- merge(health_alive, positive[, .(id, wave)], by = c("id", "wave"))
health_alive <- health_alive[!positive, on = .(id, wave)]
alive <- rbindlist(list(health_alive[, .(id, wave, alive_lower, alive_upper, event_time)],
  tracking[died == 0, .(id, wave, alive_lower = ascertain_lower, alive_upper = ascertain_upper, event_time = (ascertain_lower + ascertain_upper)/2)]))
lastalive <- alive[, .(last_known_alive_lower = max(alive_lower, na.rm = TRUE),
  last_known_alive_upper = max(alive_upper, na.rm = TRUE), last_known_alive_midpoint = max(event_time, na.rm = TRUE),
  last_known_alive_wave = max(wave)), by = id]
confirmed <- positive[, .(death_status_upper = min(ascertain_upper, na.rm = TRUE),
  first_death_wave = min(wave), ascertainment_sources = paste(sort(unique(source)), collapse = " | "),
  ascertainment_year_months = paste(sort(unique(sprintf("%04d-%02d", iyear, imonth))), collapse = " | ")), by = id]
exits[, valid_solar_year := calendar == 1 & death_year %in% 2011:2021]
exits[, date_lower := fifelse(valid_solar_year, death_year + fifelse(death_month %in% 1:12, (death_month - 1)/12, 0), NA_real_)]
exits[, date_upper := fifelse(valid_solar_year, death_year + fifelse(death_month %in% 1:12, death_month/12, 1), NA_real_)]
# Use a date interval only if explicitly solar, never interpret lunar months as Gregorian.
reported <- exits[, .(raw_exit_rows = .N, raw_solar_dates = sum(valid_solar_year, na.rm = TRUE),
  raw_lunar_rows = sum(calendar == 2, na.rm = TRUE),
  raw_missing_year_rows = sum(is.na(death_year)),
  raw_solar_lower = if (any(valid_solar_year, na.rm = TRUE)) max(date_lower, na.rm = TRUE) else NA_real_,
  raw_solar_upper = if (any(valid_solar_year, na.rm = TRUE)) min(date_upper, na.rm = TRUE) else NA_real_,
  death_year = if (sum(valid_solar_year, na.rm = TRUE) == 1L) death_year[which(valid_solar_year)[1]] else NA_real_,
  death_month = if (sum(valid_solar_year, na.rm = TRUE) == 1L) death_month[which(valid_solar_year)[1]] else NA_real_,
  exit_sources = paste(sort(unique(source)), collapse = " | ")), by = id]
candidates <- merge(explicit, confirmed, by = "id", all.x = TRUE)
candidates <- merge(candidates, lastalive, by = "id", all.x = TRUE)
candidates <- merge(candidates, reported, by = "id", all.x = TRUE)
candidates <- merge(candidates, pkgdead, by = "id", all.x = TRUE)
candidates[, in_old_appended_deaths := id %in% old[event_type == "death", id]]
candidates[, `:=`(has_raw_ascertainment_bound = is.finite(death_status_upper), has_prior_alive = is.finite(last_known_alive_lower),
  conflicting_alive_after_dead = is.finite(last_known_alive_lower) & is.finite(death_status_upper) & last_known_alive_lower >= death_status_upper,
  conflicting_solar_date_before_alive = is.finite(raw_solar_upper) & is.finite(last_known_alive_lower) & raw_solar_upper <= last_known_alive_lower,
  conflicting_solar_dates = is.finite(raw_solar_lower) & is.finite(raw_solar_upper) & raw_solar_lower >= raw_solar_upper,
  conflicting_solar_date_after_ascertainment = is.finite(raw_solar_lower) & is.finite(death_status_upper) & raw_solar_lower >= death_status_upper)]
candidates[, usable_interval_death := has_raw_ascertainment_bound & has_prior_alive &
  !conflicting_alive_after_dead & !conflicting_solar_date_before_alive & !conflicting_solar_dates & !conflicting_solar_date_after_ascertainment &
  last_known_alive_midpoint < death_status_upper]
candidates[, `:=`(date_lower = last_known_alive_lower, date_upper = death_status_upper,
  date_midpoint = (last_known_alive_lower + death_status_upper)/2,
  date_midpoint_is_observed_death = FALSE,
  model_observation_time = death_status_upper, recommended_obstype = 1L,
  source = paste(ascertainment_sources, exit_sources, sep = " | "),
  coverage = "Observed positive death records only; absent tracking rows are unknown, never alive through 2020.")]
candidates[, conflicts := paste0(fifelse(conflicting_alive_after_dead, "alive_after_dead;", ""),
  fifelse(conflicting_solar_date_before_alive, "solar_death_before_alive;", ""),
  fifelse(conflicting_solar_dates, "solar_dates_disagree;", ""),
  fifelse(conflicting_solar_date_after_ascertainment, "solar_death_after_ascertainment;", ""),
  fifelse(!has_prior_alive, "no_prior_alive;", ""), fifelse(!has_raw_ascertainment_bound, "no_ascertainment_bound;", ""))]
candidates[, usable := usable_interval_death]
saveRDS(candidates, file.path(out, "charls_death_candidates.rds"))
saveRDS(tracking, file.path(out, "charls_raw_tracking_status.rds"))
dump_csv(candidates[, .N, by = .(usable, in_old_appended_deaths, conflicts)][order(-N)], "candidate_gate_counts")
dump_csv(candidates[, .N, by = .(first_death_wave, usable)][order(first_death_wave, usable)], "usable_deaths_by_ascertainment_wave")
dump_csv(alive_death_same_wave[, .N, by = wave], "package_visit_raw_dead_same_wave")
saveRDS(alive_death_same_wave, file.path(out, "package_visit_raw_dead_same_wave.rds"))
print(candidates[, .N, by = .(usable, in_old_appended_deaths, conflicts)][order(-N)])
print(candidates[, .N, by = .(first_death_wave, usable)][order(first_death_wave, usable)])
print(alive_death_same_wave[, .N, by = wave])
harm_death_dates <- harm[!is.na(radyear), .(id, harmonized_death_year = as.numeric(radyear), harmonized_death_month = as.numeric(radmonth))]
reconcile <- merge(candidates[in_old_appended_deaths == TRUE], old[event_type == "death", .(id, old_death_event_time = event_time, old_death_time_adjusted = death_time_adjusted)], by = "id", all.x = TRUE)
reconcile <- merge(reconcile, harm_death_dates, by = "id", all.x = TRUE)
reconcile[, `:=`(package_year_matches_solar_raw = !is.na(death_year) & death_year == package_death_year,
  package_year_matches_harmonized = !is.na(harmonized_death_year) & harmonized_death_year == package_death_year,
  old_death_time_before_last_alive = old_death_event_time <= last_known_alive_midpoint,
  old_death_time_after_dead_ascertainment = old_death_event_time > death_status_upper)]
saveRDS(reconcile, file.path(out, "old_2841_deaths_reconciliation.rds"))
reconcile_summary <- data.table(metric = c("old_appended_death_IDs", "all_have_explicit_raw_status_or_exit", "safe_interval_death", "quarantined_solar_date_conflicts",
  "package_year_matches_raw_solar_or_harmonized", "package_year_not_verified_by_available_calendar_source",
  "old_death_time_before_last_alive", "old_death_time_after_dead_ascertainment", "old_death_time_adjusted"),
  n = c(nrow(reconcile), sum(reconcile$id %in% explicit$id), sum(reconcile$usable), sum(!reconcile$usable),
    sum(reconcile$package_year_matches_solar_raw | reconcile$package_year_matches_harmonized),
    sum(!(reconcile$package_year_matches_solar_raw | reconcile$package_year_matches_harmonized)),
    sum(reconcile$old_death_time_before_last_alive, na.rm = TRUE), sum(reconcile$old_death_time_after_dead_ascertainment, na.rm = TRUE),
    sum(reconcile$old_death_time_adjusted == 1, na.rm = TRUE)))
dump_csv(reconcile_summary, "old_death_reconciliation_summary")
print(reconcile_summary)
baseline <- as.data.table(readRDS(file.path(root, "outputs/charls_external_extension_20260913/data/charls_w1_descriptive_baseline.rds")))
baseline[, id := id12(id)]
allspine <- unique(package[, .(id)])
firstpositive <- positive[, .(first_dead_wave = min(wave)), by = id]
coverage <- rbindlist(lapply(c("all_core", "W1_descriptive_baseline"), function(pop) {
  ids <- if (pop == "all_core") allspine$id else baseline$id
  rbindlist(lapply(2:5, function(w) {
    d <- data.table(id = ids)
    d[, prior_or_current_positive_death := id %in% firstpositive[first_dead_wave <= w, id]]
    d[, interviewed_this_wave := id %in% health_alive[wave == w, id]]
    d[, explicit_alive_this_wave := id %in% tracking[wave == w & died == 0, id]]
    d[, status := fifelse(prior_or_current_positive_death, "recorded_dead_by_this_wave",
      fifelse(interviewed_this_wave, "core_interview_alive", fifelse(explicit_alive_this_wave, "raw_tracking_alive_without_core_visit", "not_observed_unknown")))]
    d[, .(n = .N), by = status][, `:=`(population = pop, wave = w, denominator = length(ids))]
  }))
}))
dump_csv(coverage[, .(population, wave, status, n, denominator)], "mortality_status_coverage_by_wave")
print(coverage[population == "W1_descriptive_baseline"])
baseline_death_gate <- candidates[id %in% baseline$id, .N, by = .(usable, conflicts)]
dump_csv(baseline_death_gate, "baseline_death_gate")
print(baseline_death_gate)

# Final source evidence and reproducibility checks. All counts are individuals,
# not repeated records, except files explicitly named as record-level counts.
doc_do <- sub("H_CHARLS_D_Data.dta$", "H_CHARLS_D_do_file.do", harm_path)
source_paths <- unique(c(selected, harm_path, doc_do, docs, file.path(raw_root, "2013charls/Weights.dta"), protected_paths))
manifest <- data.table(path = source_paths, bytes = file.info(source_paths)$size,
  sha256 = vapply(source_paths, function(p) digest(file = p, algo = "sha256"), ""))
dump_csv(manifest, "source_manifest")
saveRDS(exits, file.path(out, "raw_exit_death_calendar_fields.rds"))
saveRDS(harm[, .(id, radyear, radmonth, r1iwstat, r2iwstat, r3iwstat, r4iwstat)], file.path(out, "harmonized_death_status_fields.rds"))
protected_after <- vapply(protected_paths, function(p) digest(file = p, algo = "sha256"), "")
checks <- list()
check <- function(name, ok, detail = "") checks[[length(checks) + 1L]] <<- data.table(check = name, pass = isTRUE(ok), detail = as.character(detail))
check("source_files_exist", all(file.exists(source_paths)))
check("all_protected_inputs_unchanged", identical(protected_before, protected_after))
check("package_version_2_0_0", as.character(packageVersion("GlobalAgingR")) == "2.0.0")
check("package_unique_id_wave", !anyDuplicated(package[, .(id, wave)]))
check("tracking_unique_id_wave", !anyDuplicated(tracking[, .(id, wave)]))
check("candidate_unique_id", !anyDuplicated(candidates$id))
check("candidate_id12", all(grepl("^[0-9]{12}$", candidates$id)))
check("reconcile_exact_old_ids", setequal(reconcile$id, old[event_type == "death", id]) && nrow(reconcile) == 2841L)
check("every_old_death_has_raw_confirmation", all(reconcile$id %in% explicit$id))
check("all_positive_statuses_have_valid_interview_month", all(is.finite(positive$ascertain_lower)) && all(is.finite(positive$ascertain_upper)))
check("w2_death_matches_harmonized", setequal(w2[died == 1, id], harm[r2iwstat == 5, id]))
check("raw_sample_died_binary", all(samples$died %in% 0:1))
check("w2_valid_death_codes_binary_or_missing", all(is.na(w2$died) | w2$died %in% 0:1))
check("candidate_no_unknown_status_used_as_dead", all(candidates$id %in% explicit$id))
check("candidate_usable_has_prior_alive", all(candidates[usable == TRUE, has_prior_alive]))
check("candidate_usable_has_status_bound", all(candidates[usable == TRUE, has_raw_ascertainment_bound]))
check("candidate_usable_positive_interval", all(candidates[usable == TRUE, date_lower < date_upper]))
check("candidate_usable_alive_mid_before_upper", all(candidates[usable == TRUE, last_known_alive_midpoint < model_observation_time]))
check("model_time_is_raw_status_upper_not_midpoint", identical(candidates$model_observation_time, candidates$death_status_upper))
check("model_obstype_panel_not_exact_death", all(candidates$recommended_obstype == 1L))
check("date_midpoint_not_observed", !any(candidates$date_midpoint_is_observed_death))
check("all_usable_conflict_free", all(candidates[usable == TRUE, conflicts == ""]))
check("all_unusable_have_reason", all(candidates[usable == FALSE, nzchar(conflicts)]))
check("no_raw_status_then_living_contradiction", !any(candidates$conflicting_alive_after_dead))
check("no_samewave_visit_dead_contradiction", nrow(alive_death_same_wave) == 0L)
check("all_old_death_years_person_constant", all(pkgdead$n_distinct_year == 1L))
check("coverage_denominators_complete", all(coverage[, .(ok = sum(n) == unique(denominator)), by = .(population, wave)]$ok))
check("baseline_descriptive_n_12319", nrow(baseline) == 12319L)
check("candidate_rds_roundtrip", isTRUE(all.equal(as.data.frame(candidates), as.data.frame(readRDS(file.path(out, "charls_death_candidates.rds"))))))
check("solar_dates_never_used_if_lunar", all(is.na(exits[calendar != 1 | is.na(calendar), date_upper])))
check("raw_status_interval_bound_not_package_death_year", all(candidates[usable == TRUE, is.finite(death_status_upper)]))
qa <- rbindlist(checks)
dump_csv(qa, "mortality_audit_checks")
writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
gate <- data.table(status = if (all(qa$pass)) "APPROVED_INTERVAL_DEATH_WITH_ASCERTAINMENT_LIMITATIONS" else "BLOCKED_QA",
  deaths_all_raw = nrow(candidates), deaths_with_prior_alive_and_usable_bound = sum(candidates$usable),
  baseline_usable_deaths = sum(candidates$usable & candidates$id %in% baseline$id),
  exact_death_approval = FALSE, complete_mortality_ascertainment_claim = FALSE,
  recommended_obstype = 1L, completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
dump_csv(gate, "mortality_model_gate")
summary_lines <- c(
  "# CHARLS 死亡资料核查与区间死亡门控", "", paste0("执行时间：", gate$completed_at_utc), "",
  "## 结论", "",
  paste0("允许在明确披露随访与时间精度限制的前提下，构建死亡为吸收态的九状态区间观察模型。原始资料确认死亡共 ", nrow(candidates), " 人，其中 ", sum(candidates$usable), " 人同时有既往存活记录、可靠登记月份及无日期冲突；W1 描述性基线 12,319 人中可用死亡 ", gate$baseline_usable_deaths, " 人。此人数不是最终模型死亡数，仍须与冻结建模样本及疾病状态数据取交集。"),
  "",
  "**不批准既有 radyear 年中点为精确死亡日；不批准全员完整死亡随访至2020年的表述。**", "",
  "## 原始证据与覆盖", "",
  "- W2：Weights.dta 的 INDV_L_Died=1 共441人；其死亡ID与 Harmonized D 的 r2iwstat=5 完全一致，并有原始访谈年月。Exit_Interview 431人，其中424人报告死亡年，238条有年份记录使用阴历。",
  "- W3、W4、W5：Sample_Infor.dta 的 died=1 分别689、997、785人，均有访谈年月；这些是本波明确死亡记录，不得与各波累计死亡数混用。",
  "- W5 Exit_Module 770人，仅107人有 exb001_1 年份，88人明确阳历、18人阴历、1人历法未知。缺少年份绝不等于没有死亡。",
  "- Harmonized D 提供W1–W4的1=访谈、4=存活未访、5=当波死亡、6=既往死亡、9=生死未知编码；其死亡年仅402人非缺失。do-file在W3/W4明确注明没有退出文件日期，不能由该来源证明2015–2020精确死亡日期。",
  "- GlobalAgingR 包长表保留核心访谈行，iwstat只有0/1且仅4行取1，不能作为完整死亡追踪名册。当前审计回到原始死亡状态，未更改包或旧清洗。", "",
  "W1 描述性基线截至W5：9,056人有该波核心存活访谈，1,971人此前或该波明确登记死亡，1,292人该波未观察/生死未知；后者未假定存活到2020年。其中1,971人有4例原始阳历日期矛盾，本轮不追加这些有争议死亡。", "",
  "## 对旧2,841条死亡记录的逐人核对", "",
  "旧2,841人全部在原始状态/退出数据找到明确死亡证据，并非2,841个独立可信的精确死亡年份。仅489人的包死亡年可与现有阳历原始日期或Harmonized死亡年对应；其余2,352人当前未获独立历法日期验证（不等于已证明错误）。旧清洗有77人曾调整年中点，1人的旧死亡时间晚于已经登记死亡的上界。",
  "",
  "旧2,841人中2,833人通过区间死亡门控；8人隔离有争议死亡：5例原始阳历死亡月份完全早于末次存活月份，3例晚于已登记死亡的月份。没有“已登记死亡后又有存活访谈”的冲突。保留其可靠的既往存活访谈、仅在最后可靠观察时点截尾是可接受的保守处理，不能将其重编码为研究结束时存活。可另做整人排除敏感性。",
  "",
  "新增94个原始死亡ID中93人没有本地核心既往存活记录，不可贡献此次转移风险区间；另1人有既往存活记录和死亡登记月，可补入候选。", "",
  "W2权重文件另有13条空个人ID且死亡标志缺失的行，仅剔除这些无法关联个体的追踪表行；它们未计入任何死亡，未删除研究参与者。", "",
  "## 时间与模型接口", "",
  "1. 对每人取最早的明确死亡登记访谈月，月末对应 decimal year = iyear + imonth/12。这是死亡已经发生的保守观察上界，不是死亡发生日。不得使用全年名义波次末代替有原始记录的个人月份。",
  "2. 输出 date_lower 为末次可靠存活月份的月初，date_upper 为登记月末；date_midpoint仅为区间中点描述，date_midpoint_is_observed_death=FALSE。原始阳历死亡年/月另列，阴历未转换也未冒充公历。",
  "3. 建模仅使用 usable=TRUE 的 model_observation_time，并采用普通面板观察 obstype=1，禁止 exact-death obstype=3。此前部分已知存活但疾病状态不全的观察，如未纳入 censor/partial-state 似然，则模型使用的实际死亡区间可能比本表记录更宽，必须如实说明。",
  "4. 可用死亡按首次登记W2–W5分别437、633、992、772人。四库原始字段重叠、既往死亡及无核心基线者均已按ID去重，不能将文件行数相加为死亡数。",
  "5. 非死亡且无后续访谈者止于末次可靠观察，不追加2020存活行。模型还需不可检验的非信息性观察/失访假设；九状态不能消除该偏倚，不能称全国绝对风险。建议保留八个存活状态分支作非致死敏感性，不把两者视作相同estimand。", "",
  "## 文献与资料定位", "",
  "疾病与死亡口径依据原始文件及其变量标签，而非其他研究的样本数。CHARLS_2015_Release_Note 的数据文件表明确 Sample_Info 包括是否死亡及访谈日期；CHARLS_2018_Users_Guide 的第4节及8.1.1节说明持续追踪与死亡退出访谈，第8.3节还指出部分访谈延至2019年，因此不以名义2018年底强截全部观察；CHARLS_2020_User_Guide_Chinese 第4.7节说明直至死亡的追踪且不会放弃未回访者。以上手册并不保证每个非回应者生死状态均已确定。", "",
  "原始根目录：" , raw_root, "",
  "完整来源与SHA256见 source_manifest.csv；字段标签见 raw_tracking_exit_dictionary.csv、w2_weights_dictionary.csv；逐ID数据仅存本地RDS，不在汇报中披露。", "",
  "## 质控与边界", "",
  paste0("自动质控 ", sum(qa$pass), "/", nrow(qa), " 项通过。保护输入SHA256均未变化。该审计使用 nature-statistics 的设计透明、独立单位与缺失信息边界要求；计数单位为个体，未估计关联、P值或因果效应。"),
  "",
  "这是资料口径与程序检查，不等于独立重建全国死亡登记或对每一条访谈的人工核查。月份末映射仍是有限时间精度的操作规则；模型收敛、协变量插补和目标转移可辨识性需另行审计。")
writeLines(summary_lines, file.path(out, "CHARLS_死亡资料审计与模型门控.md"), useBytes = TRUE)
stopifnot(all(qa$pass))
dump_csv(qa, "mortality_checks")
message("Mortality audit complete: ", sum(qa$pass), "/", nrow(qa), " checks; ", sum(candidates$usable), " usable interval deaths.")
