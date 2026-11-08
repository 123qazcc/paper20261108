#!/usr/bin/env Rscript
# Phase 1 only: no effect estimation; all earlier HRS/ELSA files are read-only.
options(stringsAsFactors = FALSE, width = 120)
suppressPackageStartupMessages({library(data.table); library(digest)})
stopifnot(requireNamespace("GlobalAgingR", quietly = TRUE),
          requireNamespace("haven", quietly = TRUE))
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "outputs", "charls_external_extension_20260913")
for (folder in c("data", "tables", "qc", "reports", "code"))
  dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) fwrite(x, file.path(out, "tables", paste0(name, ".csv")), bom = TRUE, na = "")
sha <- function(p) digest(file = p, algo = "sha256")
protected_paths <- c("data_clean/charls_clean.rds", "data_clean/hrs_clean.rds", "data_clean/elsa_clean.rds",
  "analysis/secondary_extension_contract.md", "R/11c_fit_disease_specific_models.R")
protected <- data.table(path = protected_paths, before = vapply(protected_paths, sha, ""))
env <- new.env(parent = baseenv())
utils::data(list = "CHARLS", package = "GlobalAgingR", envir = env)
original <- env[["charls"]]
pkg <- as.data.table(GlobalAgingR::CH_filter(original, wave %in% 1:5))
num <- function(x) as.numeric(x)
valid_code <- function(x, codes) {x <- num(x); x[!x %in% codes] <- NA_real_; x}
id_string <- function(x) sprintf("%012.0f", num(x))
checks <- list()
check <- function(name, ok, detail = "") {
  checks[[length(checks) + 1L]] <<- data.table(check = name, pass = isTRUE(ok), detail = as.character(detail))
}
check("GlobalAgingR_version_2_0_0", as.character(packageVersion("GlobalAgingR")) == "2.0.0")
pkg[, id := id_string(ID)]
check("package_unique_ID_wave", !anyDuplicated(pkg[, .(id, wave)]))
check("ID_12_digit_strings", all(grepl("^[0-9]{12}$", pkg$id)))

# Raw question comparison uses checksum-verified ASCII-path copies to avoid
# haven's Windows non-ASCII path issue. Original files are hashed below too.
raw_specs <- data.table(wave = 1:5, item = c(rep("dc015", 4), "dc022"),
  staged_file = file.path("private/stage/charls_field_audit",
    c("w1_health.dta", "w2_health.dta", "w3_health.dta", "w4_cognition.dta", "w5_health.dta")),
  expected_sha256 = tolower(c("87BAE0F7E75899BA3C844985ACEAED42DD9E60483EDBEBC9C402D3CB866A810E",
    "8D2F3D9B9E289E614B30D011939433B8F67044585301B7F763FA4865BD7C2EFE",
    "9E4DC760DD1363BA1F41D551E16FA3AE998DE96D4D3796C88E95340C7C52AE5C",
    "F8F31143DDEBB134DBFAE0ECF3E88F0F9482709AD0E2561FFAB1781F9AD84792",
    "80E8C473F4F3D587B9B9AFC676B3BE353EE96C4B9D14C95E699D1F94C075584E")))
raw_checks <- rbindlist(lapply(seq_len(nrow(raw_specs)), function(i) {
  spec <- raw_specs[i]
  check(paste0("staged_raw_hash_W", i), sha(spec$staged_file) == spec$expected_sha256)
  a <- as.data.table(haven::read_dta(spec$staged_file, col_select = c("ID", spec$item)))
  a[, id := as.character(ID)]
  if (i == 1L) a[, id := paste0(substr(id, 1, 9), "0", substr(id, 10, 11))]
  a[, raw_score := valid_code(get(spec$item), 1:4)]
  j <- merge(pkg[wave == i, .(id, package_score = num(sleeprl))], a[, .(id, raw_score)], by = "id", all = TRUE)
  discord <- sum(xor(is.na(j$package_score), is.na(j$raw_score))) +
    sum(j$package_score != j$raw_score, na.rm = TRUE)
  check(paste0("raw_sleep_per_ID_W", i), discord == 0 && !anyDuplicated(a$id), discord)
  data.table(wave = i, raw_item = spec$item, raw_rows = nrow(a), package_rows = nrow(pkg[wave == i]),
             valid_sleep = sum(!is.na(j$package_score)), discordances = discord)
}))
save_table(raw_checks, "raw_sleep_per_ID_validation")

# Labels remain in an explicit dictionary, separate from numerical analysis fields.
cesd_vars <- c("botherl", "mindtsl", "depresl", "effortl", "fhopel", "fearll",
               "sleeprl", "whappyl", "flonel", "goingl")
needed <- unique(c("sleeprl", cesd_vars, "hibpe", "hearte", "stroke", "ragender", "raeducl",
                    "marry", "hrural", "hatotfa", "income_total", "iwstat", "radyear", "radmonth"))
dictionary <- rbindlist(lapply(needed, function(v) data.table(variable = v,
  label = paste(attr(original[[v]], "label"), collapse = ";"),
  value_labels = paste(names(attr(original[[v]], "labels")), attr(original[[v]], "labels"), collapse = ";"))))
save_table(dictionary, "package_variable_dictionary")
base <- as.data.table(readRDS(file.path(root, "data_clean/charls_clean.rds")))
visits <- copy(base[event_type == "visit"])
visits[, id := id_string(id)]
keep <- c("id", "wave", "event_time", "interview_year", "interview_month", "age", "sex_male",
  "education_code", "marital_status", "current_smoker", "current_drinker", "bmi",
  "sleep_hours_observed", "systolic_bp", "diastolic_bp", "pulse_pressure", "antihypertensive_medication",
  "diabetes_raw", "lung_disease_raw", "self_rated_health", "vigorous_activity_days",
  "moderate_activity_days", "light_activity_days", "hypertension_raw", "heart_disease_raw", "stroke_raw")
visits <- visits[, ..keep]
supplement <- pkg[, c("id", "wave", needed), with = FALSE]
supplement[, (needed) := lapply(.SD, num), .SDcols = needed]
d <- merge(visits, supplement, by = c("id", "wave"), all = TRUE, sort = FALSE)
setorder(d, id, event_time, wave)
check("package_clean_visit_keys_identical", nrow(d) == nrow(visits) && nrow(d) == nrow(pkg) && !anyNA(d$event_time))
check("positive_visit_times", all(is.finite(d$event_time)))
check("unique_visit_keys", !anyDuplicated(d[, .(id, wave)]))
check("sex_matches_package", all(d$sex_male == d$ragender, na.rm = TRUE) && identical(is.na(d$sex_male), is.na(d$ragender)))
check("HTN_raw_matches_package", identical(d$hypertension_raw, as.integer(valid_code(d$hibpe, 0:1))))
check("heart_raw_matches_package", identical(d$heart_disease_raw, as.integer(valid_code(d$hearte, 0:1))))
stroke_expected <- fifelse(d$stroke == 1, 1L, fifelse(d$stroke %in% c(0, 2), 0L, NA_integer_))
check("stroke_wave4_no2_mapping", identical(d$stroke_raw, stroke_expected))

d[, restless_frequency_code := valid_code(sleeprl, 1:4)]
d[, restless_score_observed := restless_frequency_code - 1]
d[, restless_high_observed := fifelse(is.na(restless_frequency_code), NA_integer_, as.integer(restless_frequency_code >= 3))]
d[, sleep_category := factor(restless_frequency_code, levels = 1:4,
  labels = c("<1 day/week", "1-2 days/week", "3-4 days/week", "5-7 days/week"))]
m <- as.matrix(d[, lapply(.SD, valid_code, codes = 1:4), .SDcols = cesd_vars])
m <- m - 1
m[, c("fhopel", "whappyl")] <- 3 - m[, c("fhopel", "whappyl")]
d[, cesd10_complete_items := rowSums(!is.na(m))]
d[, cesd10_recomputed := rowSums(m)]
d[, cesd9_complete_items := rowSums(!is.na(m[, setdiff(cesd_vars, "sleeprl"), drop = FALSE]))]
d[, cesd9_no_sleep := rowSums(m[, setdiff(cesd_vars, "sleeprl"), drop = FALSE])]
pkg_cesd <- num(pkg$cesd10[match(paste(d$id, d$wave), paste(pkg$id, pkg$wave))])
check("CESD10_complete_item_sum_matches_package", all(!is.na(pkg_cesd[!is.na(d$cesd10_recomputed)])) &&
  all(d$cesd10_recomputed == pkg_cesd, na.rm = TRUE))
check("CESD9_complete_only", identical(is.na(d$cesd9_no_sleep), d$cesd9_complete_items != 9))
check("CESD9_valid_range", all(d$cesd9_no_sleep >= 0 & d$cesd9_no_sleep <= 27, na.rm = TRUE))

# Never carry an observed negative into a currently unobserved disease assessment.
ever_positive_only <- function(x) {
  seen_yes <- cummax(replace(as.integer(x == 1L), is.na(x), 0L)) == 1L
  ifelse(seen_yes, 1L, ifelse(!is.na(x) & x == 0, 0L, NA_integer_))
}
for (v in c("hypertension", "heart_disease", "stroke")) {
  raw <- paste0(v, "_raw")
  dest <- paste0(v, "_tri")
  d[, (dest) := ever_positive_only(get(raw)), by = id]
}
d[, death_coded_visit := iwstat == 1]
d[, known_components := rowSums(!is.na(cbind(hypertension_tri, heart_disease_tri, stroke_tri)))]
d[, disease_state_bits := fifelse(known_components == 3L & !death_coded_visit,
  as.integer(hypertension_tri + 2 * heart_disease_tri + 4 * stroke_tri), NA_integer_)]
state_names <- c("No HTN/heart/stroke", "HTN only", "Heart only", "HTN + heart", "Stroke only",
                 "HTN + stroke", "Heart + stroke", "HTN + heart + stroke")
d[, disease_state := factor(disease_state_bits, levels = 0:7, labels = state_names)]
state_labels_zh <- c("无高血压/心脏病/卒中", "单纯高血压", "单纯心脏病", "高血压＋心脏病",
                     "单纯卒中", "高血压＋卒中", "心脏病＋卒中", "高血压＋心脏病＋卒中")
save_table(data.table(disease_state_bits = 0:7, label = state_names, label_zh = state_labels_zh), "disease_state_dictionary")
old <- copy(base[event_type == "visit"])
old[, id := id_string(id)]
old <- old[match(paste(d$id, d$wave), paste(id, wave))]
comparison <- rbindlist(lapply(c("hypertension", "heart_disease", "stroke"), function(v) {
  raw <- d[[paste0(v, "_raw")]]; tri <- d[[paste0(v, "_tri")]]
  data.table(wave = d$wave, disease = v, current_missing = is.na(raw),
    legacy_no_now_unknown = !is.na(old[[v]]) & old[[v]] == 0 & is.na(tri),
    raw_no_after_positive = !is.na(raw) & raw == 0 & !is.na(tri) & tri == 1)
}))[, lapply(.SD, sum), by = .(wave, disease)]
save_table(comparison, "disease_missingness_and_reversals")

availability <- d[, .(visit_records = .N, people = uniqueN(id),
  restless_observed = sum(!is.na(restless_score_observed)), sleep_hours_observed = sum(!is.na(sleep_hours_observed)),
  cesd9_observed = sum(!is.na(cesd9_no_sleep)), sbp_observed = sum(!is.na(systolic_bp)),
  dbp_observed = sum(!is.na(diastolic_bp)), hypertension_raw_observed = sum(!is.na(hypertension_raw)),
  heart_raw_observed = sum(!is.na(heart_disease_raw)), stroke_raw_observed = sum(!is.na(stroke_raw)),
  full_state_observed = sum(!is.na(disease_state_bits)), partial_state = sum(known_components %in% 1:2),
  death_coded_visit = sum(death_coded_visit)), by = wave]
save_table(availability, "wave_availability")
save_table(d[, .N, by = .(wave, restless_frequency_code)], "sleep_frequency_by_wave")
save_table(d[, .N, by = .(wave, disease_state_bits, disease_state)], "state_distribution_by_wave")
save_table(d[, .(complete10 = sum(cesd10_complete_items == 10), complete9 = sum(cesd9_complete_items == 9)), by = wave], "cesd_item_completeness")

# Description at W1 does not condition on future attendance or survival.
b0 <- copy(d[wave == 1])
b0[, eligible_age := !is.na(age) & age >= 50]
b0[, eligible_sleep := eligible_age & !is.na(restless_score_observed)]
b0[, eligible_baseline := eligible_sleep & !is.na(disease_state_bits)]
b <- copy(b0[eligible_baseline == TRUE])
later <- d[wave > 1 & !is.na(disease_state_bits), .(later_known_visits = .N), by = id]
b[later, on = "id", later_known_visits := i.later_known_visits]
b[is.na(later_known_visits), later_known_visits := 0L]
death <- copy(base[event_type == "death"])
death[, id := id_string(id)]
b[, legacy_death_record := id %in% death$id]
flow <- data.table(step = c("W1 respondents", "Age >=50", "Plus observed restless sleep", "Plus known baseline disease combination",
  "Subset: later known living state (not final model N)", "Subset: no later known living state", "Subset: legacy death record, coverage unverified"),
  n = c(nrow(b0), sum(b0$eligible_age), sum(b0$eligible_sleep), nrow(b),
    sum(b$later_known_visits > 0), sum(b$later_known_visits == 0), sum(b$legacy_death_record)))
save_table(flow, "baseline_eligibility_flow")
saveRDS(d, file.path(out, "data/charls_observed_exposures_tristate_panel.rds"))
saveRDS(b, file.path(out, "data/charls_w1_descriptive_baseline.rds"))

# Two-wave exposure patterns only; downstream analysis must start after W2.
dyn <- merge(b[, .(id, w1_time = event_time, w1_score = restless_score_observed, w1_high = restless_high_observed)],
  d[wave == 2, .(id, w2_time = event_time, w2_score = restless_score_observed,
                w2_high = restless_high_observed, w2_state = disease_state_bits)], by = "id", all.x = TRUE)
dyn[, elapsed := w2_time - w1_time]
dyn[, pattern := fifelse(is.na(w2_high), "Unobserved at W2", c("Persistent low", "Onset", "Remission", "Persistent high")[1 + 2*w1_high + w2_high])]
dyn[, cumulative_score_years := fifelse(elapsed > 0, (w1_score + w2_score)/2 * elapsed, NA_real_)]
dyn[, mean_two_wave_score := cumulative_score_years / elapsed]
dyn[, valid_W2_landmark := !is.na(w2_score) & !is.na(w2_state) & !is.na(elapsed) & elapsed > 0]
saveRDS(dyn, file.path(out, "data/charls_w1_w2_exposure_patterns.rds"))
save_table(dyn[, .(n = .N, valid_W2_landmark = sum(valid_W2_landmark)), by = pattern], "w1_w2_dynamic_pattern_counts")

# Consecutive known states may bridge missing visits; these are NOT q_ij estimates.
p <- d[id %in% b$id & !is.na(disease_state_bits)]
setorder(p, id, event_time)
p[, `:=`(next_time = shift(event_time, type = "lead"), next_state = shift(disease_state_bits, type = "lead"),
         next_wave = shift(wave, type = "lead")), by = id]
intervals <- p[!is.na(next_time), .(id, from_wave = wave, to_wave = next_wave, from = disease_state_bits,
  to = next_state, elapsed_years = next_time - event_time, start_restless_score = restless_score_observed)]
intervals[, compatible_increment := bitwAnd(from, to) == from]
check("strictly_increasing_known_visit_times", all(intervals$elapsed_years > 0))
check("no_loss_of_ever_positive_disease", all(intervals$compatible_increment))
counts <- intervals[, .(observed_intervals = .N, people = uniqueN(id),
  skipped_wave_intervals = sum(to_wave - from_wave > 1)), by = .(from, to)]
counts[, `:=`(from_label = state_names[from + 1], to_label = state_names[to + 1])]
setorder(counts, from, to)
targets <- data.table(from = c(0L, 0L, 1L, 1L), to = c(2L, 4L, 3L, 5L))
target_counts <- merge(targets, counts, by = c("from", "to"), all.x = TRUE, sort = FALSE)
save_table(counts, "all_observed_interval_state_pairs")
save_table(target_counts, "four_target_observed_interval_counts")
saveRDS(intervals, file.path(out, "data/charls_observed_known_state_intervals.rds"))

# Baseline descriptive table: no hypothesis tests; available-case denominators explicit.
continuous <- c(age = "年龄（岁）", bmi = "BMI", sleep_hours_observed = "夜间睡眠时长（小时）",
  cesd9_no_sleep = "CESD-9（剔除睡眠题）", systolic_bp = "收缩压（mmHg）", diastolic_bp = "舒张压（mmHg）")
categorical <- c(sex_male = "男性", current_smoker = "当前吸烟", current_drinker = "当前饮酒",
                 diabetes_raw = "糖尿病", lung_disease_raw = "慢性肺病")
groups <- c("Overall", levels(b$sleep_category))
tab <- rbindlist(lapply(groups, function(g) {
  z <- if (g == "Overall") b else b[as.character(sleep_category) == g]
  rows <- list(data.table(variable = "n", label = "人数", group = g, n_group = nrow(z), n_available = nrow(z),
    missing_n = 0L, summary = as.character(nrow(z))))
  for (v in names(continuous)) {
    a <- z[[v]]; good <- a[!is.na(a)]
    val <- if (length(good)) sprintf("%.2f (%.2f)", mean(good), sd(good)) else "NA"
    rows[[length(rows)+1L]] <- data.table(variable = v, label = continuous[[v]], group = g,
      n_group = nrow(z), n_available = length(good), missing_n = sum(is.na(a)), summary = val)
  }
  for (v in names(categorical)) {
    a <- z[[v]]; n <- sum(!is.na(a)); yes <- sum(a == 1, na.rm = TRUE)
    rows[[length(rows)+1L]] <- data.table(variable = v, label = categorical[[v]], group = g,
      n_group = nrow(z), n_available = n, missing_n = sum(is.na(a)),
      summary = if (n) sprintf("%d/%d (%.1f%%)", yes, n, 100*yes/n) else "NA")
  }
  multicat <- list(raeducl = c(`1` = "教育：低于初中（包标签）", `2` = "教育：高中/职业培训", `3` = "教育：高等教育"),
    marry = c(`1` = "已婚同住", `2` = "已婚不同住", `3` = "同居", `4` = "分居", `5` = "离异", `7` = "丧偶", `8` = "从未结婚"),
    hrural = c(`0` = "城市居住", `1` = "农村居住"))
  for (v in names(multicat)) for (code in names(multicat[[v]])) {
    a <- z[[v]]; den <- sum(!is.na(a)); n_yes <- sum(a == as.numeric(code), na.rm = TRUE)
    rows[[length(rows)+1L]] <- data.table(variable = paste(v, code, sep = "_"), label = multicat[[v]][[code]],
      group = g, n_group = nrow(z), n_available = den, missing_n = sum(is.na(a)),
      summary = if (den) sprintf("%d/%d (%.1f%%)", n_yes, den, 100*n_yes/den) else "NA")
  }
  for (state in 0:7) rows[[length(rows)+1L]] <- data.table(variable = paste0("state_", state),
    label = state_labels_zh[state+1], group = g, n_group = nrow(z), n_available = nrow(z), missing_n = 0L,
    summary = sprintf("%d (%.1f%%)", sum(z$disease_state_bits == state), 100*mean(z$disease_state_bits == state)))
  rbindlist(rows)
}))
save_table(tab, "baseline_table_long_with_denominators")
wide <- dcast(tab, variable + label ~ group, value.var = "summary")
wide <- wide[match(unique(tab$variable), variable)]
setcolorder(wide, c("variable", "label", groups))
save_table(wide, "baseline_table_by_restless_sleep")
group_labels_zh <- c("总体", "睡眠不安＜1天/周", "1–2天/周", "3–4天/周", "5–7天/周")
md <- c("# CHARLS W1 基线表（第一阶段描述性样本）", "",
  "年龄≥50岁、W1睡眠不安题有效、三种疾病组件可确定；不以未来回访或存活筛选。不是最终模型样本量。",
  "", "以下为未加调查权重的描述统计，不代表中国总体患病率。连续变量为均值（标准差）；分类变量以非缺失者为分母；缺失人数另见配套长表。不进行基线显著性检验。", "",
  paste0("| 变量 | ", paste(group_labels_zh, collapse = " | "), " |"),
  paste0("|---|", paste(rep("---", length(groups)), collapse = "|"), "|"))
for (i in seq_len(nrow(wide))) md <- c(md, paste0("| ", wide$label[i], " | ", paste(unlist(wide[i, ..groups]), collapse = " | "), " |"))
writeLines(md, file.path(out, "reports/CHARLS_第一阶段基线表.md"), useBytes = TRUE)

source_root <- Sys.getenv("CHARLS_SOURCE_ROOT", "private/source")
raw_root <- file.path(source_root, "data_raw/1_CHARLS_China")
teaching_root <- file.path(source_root, "reference_scripts")
source_files <- c(file.path(root, "data_clean/charls_clean.rds"), file.path(root, "R/13a_prepare_charls_external_extension.R"),
  file.path(root, "analysis/charls_external_extension_contract.md"),
  file.path(teaching_root, "msm_example.R"),
  file.path(raw_root, c("2011charls/health_status_and_functioning.dta", "2013charls/Health_Status_and_Functioning.dta",
    "2015charls/Health_Status_and_Functioning.dta", "2018charls/Cognition.dta", "2020charls/Health_Status_and_Functioning.dta")))
pkg_home <- system.file(package = "GlobalAgingR")
pkg_files <- c(list.files(file.path(pkg_home, "data"), full.names = TRUE),
               list.files(file.path(pkg_home, "R"), full.names = TRUE), file.path(pkg_home, c("DESCRIPTION", "NAMESPACE")))
source_files <- unique(c(source_files, raw_specs$staged_file, pkg_files))
manifest <- data.table(path = source_files, exists = file.exists(source_files))
manifest[, sha256 := vapply(path, function(p) if(file.exists(p)) sha(p) else NA_character_, "")]
check("all_source_files_present", all(manifest$exists))
for (i in 1:5) {
  original_sha <- manifest$sha256[match(source_files[4+i], manifest$path)]
  check(paste0("original_staged_raw_hash_match_W", i), identical(original_sha, raw_specs$expected_sha256[i]))
}
fwrite(manifest, file.path(out, "qc/source_manifest.csv"), bom = TRUE)
protected[, after := vapply(path, sha, "")]
check("protected_files_unchanged", all(protected$before == protected$after))
fwrite(protected, file.path(out, "qc/protected_file_hashes.csv"))
fwrite(rbindlist(checks), file.path(out, "qc/phase1_checks.csv"))
writeLines(capture.output(sessionInfo()), file.path(out, "qc/session_info.txt"))
saveRDS(list(status = "PHASE1_DESCRIPTIVE_ONLY_MODELS_NOT_RUN", package_version = as.character(packageVersion("GlobalAgingR")),
  n_baseline = nrow(b), n_followup_candidates = sum(b$later_known_visits > 0), started_from = source_files,
  target_counts = target_counts, complete_checks = all(rbindlist(checks)$pass)), file.path(out, "qc/phase1_status.rds"))
file.copy(file.path(root, "R/13a_prepare_charls_external_extension.R"), file.path(out, "code/13a_prepare_charls_external_extension.R"), overwrite = TRUE)
print(flow)
print(target_counts)
print(rbindlist(checks))
if (!all(rbindlist(checks)$pass)) stop("Phase 1 check failure: inspect qc/phase1_checks.csv")
message("CHARLS phase 1 complete. No models fitted; mortality and raw disease routing review remain open.")
