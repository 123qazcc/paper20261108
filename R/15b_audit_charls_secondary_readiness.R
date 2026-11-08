#!/usr/bin/env Rscript
# Read-only coverage audit: no model, imputation, filled exposures or person-level exports.
options(stringsAsFactors = FALSE, width = 140)
suppressPackageStartupMessages({library(data.table); library(digest)})
setDTthreads(1)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "outputs/charls_external_extension_20260913")
dest <- file.path(out, "secondary_readiness")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
files <- c(bundle = file.path(out, "data/charls_model_bundle.rds"),
  observed = file.path(out, "data/charls_observed_exposures_tristate_panel.rds"),
  raw_baseline = file.path(out, "data/charls_raw_W1_baseline_stage2.rds"),
  tracking = file.path(out, "mortality_audit/charls_raw_tracking_status.rds"),
  death_audit = file.path(out, "mortality_audit/charls_death_candidates.rds"),
  phase1_contract = file.path(root, "analysis/charls_external_extension_contract.md"),
  primary_contract = file.path(root, "analysis/charls_model_stage2_contract.md"),
  phase1_code = file.path(root, "R/13a_prepare_charls_external_extension.R"),
  prep_code = file.path(root, "R/14a_prepare_charls_model_bundle.R"),
  cleaning_code = file.path(root, "R/01_clean_three_cohorts.R"),
  audit_code = file.path(root, "R/15b_audit_charls_secondary_readiness.R"))
stopifnot(all(file.exists(files)))
sha <- function(p) digest(file = p, algo = "sha256")
before <- vapply(files, sha, "")
checks <- list()
addcheck <- function(name, pass, detail = "", severity = "HARD_VALIDITY") {
  checks[[length(checks) + 1L]] <<- data.table(check = name, pass = isTRUE(pass), detail = as.character(detail), severity = severity)
}
put <- function(x, name) {
  stopifnot(!any(names(x) %in% c("id", "ID", "participant_id", "subject_id")))
  fwrite(x, file.path(dest, paste0(name, ".csv")), bom = TRUE, na = "")
}
bun <- readRDS(files[["bundle"]])
p <- copy(as.data.table(bun$panel))
b <- copy(as.data.table(bun$baseline))
d <- as.data.table(readRDS(files[["observed"]]))
rawb <- as.data.table(readRDS(files[["raw_baseline"]]))
track <- as.data.table(readRDS(files[["tracking"]]))
death_audit <- as.data.table(readRDS(files[["death_audit"]]))
addcheck("approved_primary_bundle", identical(bun$status, "APPROVED_FOR_PILOT_AND_MI"))
addcheck("main_people_12036_unique", nrow(b) == 12036L && !anyDuplicated(b$id))
addcheck("main_panel_53786_records", nrow(p) == 53786L)
addcheck("unique_observed_person_wave", !anyDuplicated(d[, .(id, wave)]))
addcheck("living_no_registered_death_same_wave", nrow(merge(p[record_type == "living_visit"], track[died == 1], by = c("id", "wave"))) == 0L)
addcheck("primary_input_hashes_current", all(vapply(names(bun$input_hashes), sha, "") == bun$input_hashes))

# Only use package data here to recover unfilled source dates and audit the
# inherited night-sleep validity filter. Disease states come from the raw branch.
env <- new.env(parent = baseenv())
utils::data(list = "CHARLS", package = "GlobalAgingR", envir = env)
pkg <- as.data.table(GlobalAgingR::CH_filter(env[["charls"]], wave %in% 1:5))
source <- pkg[, .(id = sprintf("%012.0f", as.numeric(ID)), wave = as.integer(wave),
  source_year = as.numeric(iwy), source_month = as.numeric(iwm), night_raw = as.numeric(sleep_night))]
rm(pkg, env)
source[, `:=`(source_year_valid = is.finite(source_year) & source_year >= 1990 & source_year <= 2030,
  source_month_valid = is.finite(source_month) & source_month >= 1 & source_month <= 12)]
source[, time_quality := fcase(source_year_valid & source_month_valid, "recorded_year_month_midpoint",
  source_year_valid, "recorded_year_month_missing_midyear_proxy", default = "year_missing_nominal_wave_proxy")]
addcheck("source_date_keys_unique", !anyDuplicated(source[, .(id, wave)]))
d <- merge(d, source, by = c("id", "wave"), all.x = TRUE, sort = FALSE)
d <- d[id %in% b$id]
addcheck("all_observed_visits_have_source_record", !anyNA(d$time_quality))
valid_night <- ifelse(is.finite(d$night_raw) & d$night_raw >= 1 & d$night_raw <= 18, d$night_raw, NA_real_)
addcheck("night_filter_reproduced", identical(as.numeric(d$sleep_hours_observed), as.numeric(valid_night)))
addcheck("CESD9_requires_all_nine_items", identical(is.na(d$cesd9_no_sleep), d$cesd9_complete_items != 9L))
addcheck("CESD9_range_0_27", all(d$cesd9_no_sleep >= 0 & d$cesd9_no_sleep <= 27, na.rm = TRUE))

state_names <- c("无高血压/心脏病/卒中", "单纯高血压", "单纯心脏病", "高血压+心脏病",
  "单纯卒中", "高血压+卒中", "心脏病+卒中", "三病共存")
state_set <- function(s) {
  if (s %in% 1:8) return(as.integer(s))
  if (s == 9L) return(9L)
  z <- match(s, bun$censor_codes)
  if (is.na(z)) stop("Unknown panel state code")
  as.integer(bun$censor_states[[z]])
}
dict <- data.table(state = sort(unique(p$state)))
dict[, compatible_n := lengths(lapply(state, state_set))]
dict[, compatible_states := vapply(state, function(s) paste(state_set(s), collapse = ";"), "")]
dict[, state_class := fifelse(state == 9L, "death", fifelse(compatible_n == 1, "fully_known",
  fifelse(compatible_n == 8, "all_components_unknown_alive", "partially_known_alive")))]
dict[, label := vapply(state, function(s) if (s %in% 1:8) state_names[s] else if (s == 9L) "死亡" else paste0("部分观测代码", s), "")]
put(dict, "state_set_dictionary")

w1 <- merge(b[, .(id, restless_bl)], d[wave == 1], by = "id", all.x = TRUE, sort = FALSE)
w1 <- merge(w1, p[wave == 1 & record_type == "living_visit", .(id, primary_state = state)], by = "id", all.x = TRUE)
addcheck("baseline_joins_one_row_per_person", nrow(w1) == nrow(b) && !anyNA(w1$primary_state) && !anyDuplicated(w1$id))
addcheck("W1_state_matches_raw_branch_baseline", identical(as.integer(w1$primary_state),
  as.integer(rawb$bits[match(w1$id, rawb$id)]) + 1L))
addcheck("primary_exposure_matches_observed", identical(w1$restless_bl, as.numeric(w1$restless_score_observed)))
w1[, `:=`(cesd9_available = !is.na(cesd9_no_sleep), night_available = !is.na(sleep_hours_observed),
  binary = restless_high_observed, frequency = restless_frequency_code)]
coverage <- function(z, group, level) data.table(group = group, level = as.character(level), n = nrow(z),
  cesd9_observed = sum(z$cesd9_available), cesd9_missing = sum(!z$cesd9_available),
  night_observed = sum(z$night_available), night_missing = sum(!z$night_available),
  both_observed = sum(z$cesd9_available & z$night_available))
basecov <- coverage(w1, "overall", "all")
for (v in c("primary_state", "frequency", "binary")) {
  basecov <- rbind(basecov, rbindlist(lapply(sort(unique(w1[[v]])), function(x) coverage(w1[get(v) == x], v, x))))
}
put(basecov, "W1_secondary_coverage")
put(w1[, .N, by = .(cesd9_complete_items)][order(cesd9_complete_items)], "W1_CESD9_item_completion")
w1[, night_missing_reason := fcase(night_available, "valid_1_18_hours", is.na(night_raw), "source_missing_or_special_missing",
  !is.finite(night_raw), "nonfinite_source", default = "outside_prespecified_1_18_hours")]
put(w1[, .N, by = night_missing_reason], "W1_night_missing_reasons")

# Genuine living visits must exist in the frozen raw-state panel. No historical
# harmonized disease_state_bits is used for landmark classification.
live <- merge(p[record_type == "living_visit"], d[, .(id, wave, event_time, iwstat, restless_score_observed,
  restless_high_observed, sleep_hours_observed, cesd9_no_sleep, time_quality)], by = c("id", "wave"), all.x = TRUE, sort = FALSE)
live <- merge(live, dict[, .(state, state_class, compatible_n)], by = "state", all.x = TRUE)
addcheck("all_living_panel_visits_matched", !anyNA(live$event_time))
# The raw branch deliberately follows the raw status audit, not retrospectively
# modified package iwstat. Document disagreements instead of declaring them clear.
status_conflicts <- merge(live[iwstat != 0], track[, .(id, wave, raw_died = died)], by = c("id", "wave"), all.x = TRUE)
status_conflicts <- merge(status_conflicts, death_audit[, .(id, death_usable = usable, conflicts)], by = "id", all.x = TRUE)
addcheck("all_living_interview_flags_zero", all(live$iwstat == 0),
  paste0(nrow(status_conflicts), " known W4 raw/package disagreements; formal time-updated preparation blocked for review"),
  severity = "REVIEW_REQUIRED")
addcheck("package_status_disagreements_identified_in_death_audit", all(status_conflicts$raw_died == 0 &
  !status_conflicts$death_usable & grepl("solar_death_before_alive", status_conflicts$conflicts)))
put(status_conflicts[, .(records = .N, observed_sleep = sum(!is.na(restless_score_observed))),
  by = .(wave, iwstat, raw_died, death_usable, conflicts)], "living_package_status_conflicts")
entry <- w1[, .(id, entry_time = event_time)]
live <- merge(live, entry, by = "id", all.x = TRUE)
addcheck("observed_source_time_matches_panel", max(abs(live$event_time - live$entry_time - live$time)) < 1e-8)
lm <- merge(w1[, .(id, w1_high = restless_high_observed, w1_score = restless_score_observed, w1_state = primary_state)],
  live[wave == 2, .(id, w2_time = time, w2_score = restless_score_observed, w2_high = restless_high_observed,
    w2_state = state, w2_class = state_class, w2_compatible_n = compatible_n, time_quality)], by = "id", all.x = TRUE)
lm[, w2_living := !is.na(w2_time)]
lm[, w2_sleep := w2_living & !is.na(w2_score)]
lm[, descriptive_candidate := w2_sleep & is.finite(w2_time) & w2_time > 0]
later <- merge(p, lm[descriptive_candidate == TRUE, .(id, w2_time)], by = "id")[time > w2_time]
later[, known_or_partial_living := record_type == "living_visit" &
  !state %in% dict[state_class == "all_components_unknown_alive", state]]
follow <- later[, .(later_n = .N, later_living = sum(record_type == "living_visit"),
  later_disease_information = sum(known_or_partial_living), later_death = sum(state == 9L),
  followup_years = max(time) - first(w2_time)), by = id]
lm <- merge(lm, follow, by = "id", all.x = TRUE)
for (v in c("later_n", "later_living", "later_disease_information", "later_death")) set(lm, which(is.na(lm[[v]])), v, 0L)
lm[, followup_candidate := descriptive_candidate & later_n > 0]
lm[, complete_state_followup := followup_candidate & w2_class == "fully_known"]
lm[, partial_state_followup := followup_candidate & w2_class == "partially_known_alive"]
lm[, all_unknown_followup := followup_candidate & w2_class == "all_components_unknown_alive"]
lm[, pattern := fifelse(descriptive_candidate, c("持续低频", "新出现高频", "由高频缓解", "持续高频")[1 + 2*w1_high + w2_high], NA_character_)]
lm[, mutually_exclusive_reason := fcase(!w2_living, "no_reliable_W2_living_interview",
  !w2_sleep, "W2_sleep_unobserved", !descriptive_candidate, "W2_time_not_after_W1",
  later_n == 0, "no_reliable_information_after_W2", complete_state_followup, "ready_known_W2_state_with_followup",
  partial_state_followup, "candidate_partial_W2_state_needs_initial_state_policy", default = "candidate_unknown_W2_state_needs_initial_state_policy")]
flow <- data.table(step = c("main_W1_analysis_population", "plus_reliable_W2_living_interview", "plus_observed_W2_sleep",
  "plus_positive_W1_W2_elapsed_descriptive_candidates", "plus_any_reliable_observation_after_W2",
  "subset_fully_known_W2_state_with_followup", "subset_partial_W2_state_with_followup", "subset_unknown_W2_state_with_followup"),
  n = c(nrow(lm), sum(lm$w2_living), sum(lm$w2_sleep), sum(lm$descriptive_candidate), sum(lm$followup_candidate),
    sum(lm$complete_state_followup), sum(lm$partial_state_followup), sum(lm$all_unknown_followup)))
put(flow, "W2_landmark_flow")
put(lm[, .N, by = mutually_exclusive_reason], "W2_landmark_exclusion_reasons")
put(lm[, .N, by = .(w1_state, mutually_exclusive_reason)], "W2_landmark_reasons_by_W1_state")
pattern_counts <- lm[descriptive_candidate == TRUE, .(descriptive_candidates = .N, post_W2_followup_candidates = sum(followup_candidate),
  fully_known_W2_with_followup = sum(complete_state_followup), partial_W2_with_followup = sum(partial_state_followup),
  all_unknown_W2_with_followup = sum(all_unknown_followup), later_disease_information_people = sum(later_disease_information > 0),
  later_death_people = sum(later_death > 0), death_only_followup_people = sum(later_death > 0 & later_living == 0)), by = pattern]
put(pattern_counts, "W1_W2_pattern_readiness")
put(lm[descriptive_candidate == TRUE, .(n_candidates = .N, with_followup = sum(followup_candidate)), by = .(w2_state, w2_class, pattern)], "W2_patterns_by_raw_state")
put(lm[descriptive_candidate == TRUE, .(n_candidates = .N, with_followup = sum(followup_candidate)), by = time_quality], "W2_landmark_time_source")
addcheck("landmark_partition", sum(lm$complete_state_followup) + sum(lm$partial_state_followup) + sum(lm$all_unknown_followup) == sum(lm$followup_candidate))
addcheck("landmark_does_not_use_historical_state", all(lm$w2_state[lm$w2_living] == live[wave == 2]$state[match(lm$id[lm$w2_living], live[wave == 2]$id)]))
addcheck("all_followup_strictly_after_landmark", nrow(later) > 0 && all(later$time > later$w2_time))
without_conflicting_later <- later[!status_conflicts[, .(id, wave)], on = .(id, wave)]
addcheck("W2_followup_readiness_not_changed_by_known_W4_conflicts",
  setequal(lm[followup_candidate == TRUE, id], without_conflicting_later$id))

# Readiness for an observed-start time-updated sensitivity analysis. Aggregate
# adjacent intervals only; no future backfill and no carry-forward-filled file.
put(live[, .(living_visits = .N, people = uniqueN(id), restless_observed = sum(!is.na(restless_score_observed)),
  restless_missing = sum(is.na(restless_score_observed)), night_observed = sum(!is.na(sleep_hours_observed)),
  cesd9_observed = sum(!is.na(cesd9_no_sleep))), by = .(wave, state_class)], "time_updated_wave_state_coverage")
put(live[, .N, by = .(wave, time_quality)], "time_updated_time_source")
ip <- merge(p, live[, .(id, wave, event_time, restless_score_observed, sleep_hours_observed, state_class)],
  by = c("id", "wave"), all.x = TRUE, sort = FALSE)
setorder(ip, id, time)
ip[, `:=`(end_time = shift(time, type = "lead"), end_wave = shift(wave, type = "lead"),
  end_state = shift(state, type = "lead"), end_record_type = shift(record_type, type = "lead")), by = id]
ip <- ip[!is.na(end_time)]
ip[, `:=`(elapsed = end_time - time, start_sleep_observed = !is.na(restless_score_observed),
  start_night_observed = !is.na(sleep_hours_observed), skipped_wave = end_wave - wave > 1)]
addcheck("time_updated_all_intervals_positive", all(ip$elapsed > 0))
addcheck("no_intervals_start_after_death", all(ip$record_type == "living_visit"))
addcheck("exposure_time_from_same_start_visit", all(is.finite(ip$event_time)))
status_conflicts[, conflicting_start_intervals := ip[.SD, on = .(id, wave), .N, by = .EACHI]$N]
status_conflicts[, conflicting_end_intervals := ip[.SD, on = .(id, end_wave = wave), .N, by = .EACHI]$N]
put(status_conflicts[, .(records = .N, observed_sleep = sum(!is.na(restless_score_observed)),
  conflicting_start_intervals = sum(conflicting_start_intervals), conflicting_end_intervals = sum(conflicting_end_intervals)),
  by = .(wave, iwstat, raw_died, death_usable, conflicts)], "living_package_status_conflicts")
# Explicitly requested by the coordinating audit: private locator, never included
# in public CSV/Markdown. Keep only the fields needed to inspect source conflicts.
dir.create(file.path(dest, "internal"), showWarnings = FALSE)
saveRDS(status_conflicts[, .(id, wave, event_time, iwstat, raw_died, death_usable, conflicts,
  restless_score_observed, conflicting_start_intervals, conflicting_end_intervals)],
  file.path(dest, "internal/living_status_conflict_locators.rds"))
tv <- ip[, .(intervals = .N, people = uniqueN(id), observed_start_intervals = sum(start_sleep_observed),
  missing_start_intervals = sum(!start_sleep_observed), observed_night_start_intervals = sum(start_night_observed),
  skipped_wave_intervals = sum(skipped_wave), interval_years = sum(elapsed),
  years_with_observed_start_sleep = sum(elapsed[start_sleep_observed])), by = .(wave, state_class, end_record_type)]
put(tv, "time_updated_observed_start_interval_coverage")
person_tv <- ip[, .(has_observed_start = any(start_sleep_observed), any_missing_start = any(!start_sleep_observed),
  all_starts_observed = all(start_sleep_observed), observed_start_count = sum(start_sleep_observed),
  observed_updates_after_W1 = sum(start_sleep_observed & wave > 1)), by = id]
put(person_tv[, .N, by = .(has_observed_start, any_missing_start, all_starts_observed, observed_updates_after_W1)], "time_updated_person_coverage")

after <- vapply(files, sha, "")
addcheck("all_fixed_inputs_unchanged", identical(before, after))
qc <- rbindlist(checks)
put(qc, "audit_checks")
put(data.table(role = names(files), path = unname(files), sha256_before = unname(before), sha256_after = unname(after)), "source_manifest")
hard_ok <- all(qc[severity == "HARD_VALIDITY", pass])
gates <- data.table(module = c("fixed_input_and_coverage_integrity", "W1_secondary_exposure_preparation", "W2_landmark_data_preparation", "time_updated_exposure_preparation"),
  status = c(if (hard_ok) "PASSED_HARD_VALIDITY" else "FAILED_HARD_VALIDITY",
    "READY_FOR_DATA_PREPARATION", "READY_FOR_DATA_PREPARATION", "REVIEW_REQUIRED"),
  boundary = c("Known disagreements remain separately flagged, not converted to agreement",
    "CESD9 missingness requires a revised imputation plan; no effects fitted",
    "Known-state N=9397; partial/unknown initial-state policy and appropriate MI plan unresolved",
    "Two disputed outcome endpoints; one observed sleep value but neither row starts an interval; review before formal construction"))
put(gates, "module_gates")
stopifnot(hard_ok)
overall <- basecov[group == "overall"]
lines <- c("# CHARLS 次要暴露及敏感性分析可用性审计", "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "。仅为数据覆盖审计，不拟合模型、不插补、不推断关联。"), "",
  "## 范围与分母", "",
  paste0("固定主队列为 ", nrow(b), " 名参与者、", nrow(p), " 条既定模型观测（其中两条争议 W4 存活行另列复核）。独立单位为人；访谈/区间是重复观测，不是独立样本。采用最终原始疾病问答链及前瞻性既往阳性口径，不使用第一阶段被后续回答回改的历史疾病状态或历史动态人数。"), "",
  "## W1 次要指标", "",
  sprintf("CESD-9：九题完整 %d 人，缺失 %d 人；得分剔除睡眠题，两个正向题反向计分。不是包内不完整总分减睡眠分。", overall$cesd9_observed, overall$cesd9_missing),
  sprintf("夜间睡眠时长：1–18 小时有效观测 %d 人，缺失/超范围 %d 人；与 CESD-9 同时可用 %d 人。缺失细分见 W1_night_missing_reasons.csv。此指标是时长，不是失眠负担。", overall$night_observed, overall$night_missing, overall$both_observed),
  paste0("W1 四档频率人数（<1、1–2、3–4、5–7 天/周）：", paste(w1[, .N, by = frequency][order(frequency)]$N, collapse = "、"), "。"),
  paste0("W1 二分类人数（低频<3天、高频≥3天/周）：", paste(w1[, .N, by = binary][order(binary)]$N, collapse = "、"), "。分类和连续每级效应属于不同分析规格，不按显著性互换。"), "",
  "## W1–W2 动态 landmark", "",
  paste0("主队列 → W2 真实存活访谈 ", sum(lm$w2_living), " 人 → W2 睡眠有观测且晚于 W1 的描述候选 ", sum(lm$descriptive_candidate),
    " 人 → W2 后至少一条可靠存活或死亡观测 ", sum(lm$followup_candidate), " 人。"),
  paste0("具有 W2 后随访者中：W2 状态完整 ", sum(lm$complete_state_followup), " 人，部分组件已知 ", sum(lm$partial_state_followup),
    " 人，全部组件未知但确认存活 ", sum(lm$all_unknown_followup), " 人。后两类只列兼容状态集合候选，尚未批准初始状态分布或左条件化方法，不能直接当成已就绪模型人数。"),
  "所有动态效应的风险时间必须从 W2 开始；W1–W2 只能定义过去的暴露模式，不算作模式确定后的结局随访。无后续访谈不等于无事件；只有可靠死亡记录者也可提供后续信息。死亡不是精确死亡日。", "",
  "| 模式 | W2 描述候选 | 有 W2 后信息 | W2 状态完整且有随访 | W2 部分状态且有随访 | W2 全未知且有随访 |",
  "|---|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(pattern_counts))) {
  z <- pattern_counts[i]
  lines <- c(lines, sprintf("| %s | %d | %d | %d | %d | %d |", z$pattern, z$descriptive_candidates, z$post_W2_followup_candidates,
    z$fully_known_W2_with_followup, z$partial_W2_with_followup, z$all_unknown_W2_with_followup))
}
lines <- c(lines, "", "## 时间更新暴露：只计来源覆盖，不作填充值", "",
  sprintf("共 %d 个既定相邻模型观测区间（含下述争议终点）；起点当波睡眠实测 %d 个，起点睡眠缺失 %d 个。%d 人的所有区间起点都有实测；%d 人至少一个起点缺失。",
    nrow(ip), sum(ip$start_sleep_observed), sum(!ip$start_sleep_observed), sum(person_tv$all_starts_observed), sum(person_tv$any_missing_start)),
  "使用同一个人的当前访谈睡眠值及其来源时间进行覆盖汇总，不用下一波睡眠回填前一时段；也没有生成可直接拟合的前向填充数据。把起点值视为整个区间的暴露仍是未来模型需要明确的阶梯函数/陈旧测量假设，本审计不批准该假设。跨缺失波次区间单列；删除睡眠缺失区间可能改变风险集并产生选择偏倚，不能直接当成完整主样本的替代。",
  paste0("来源时间从包内未填充 iwy/iwm 重新核对。本次 ", sum(live$time_quality == "recorded_year_month_midpoint"), "/", nrow(live),
    " 条存活访谈有记录年月，采用月中点，未触发缺月/缺年的历史代理规则；见 time_updated_time_source.csv。记录年月不意味着已知精确访谈日。"), "",
  sprintf("已识别的死亡日期/存活状态矛盾：主分支中 %d 条 W4 行的原始追踪 died=0，但包 iwstat=1；W5 退出记录给出的死亡年月早于 W4 访谈。这些人已经被既有死亡审计标为不接受的死亡终点，主分支仍保留原始 W4 存活行。本审计不擅自改其状态；其中 %d 条有睡眠观测，涉及 %d 个后续区间起点、%d 个区间终点。W2 模式定义不受这些 W4 标记影响；观测真实性和保守截尾敏感性仍需复核，不能把相互冲突的来源写成完全一致。",
    nrow(status_conflicts), sum(!is.na(status_conflicts$restless_score_observed)), sum(status_conflicts$conflicting_start_intervals), sum(status_conflicts$conflicting_end_intervals)),
  "独立覆盖复核表明，即使不计这两条 W4 争议后续行，W2 后仍有信息的参与者集合不变；因此 W2 landmark 数据准备门控可单独放行，但这些争议记录不纳入新的模型。", "",
  "这两个 W4 行都没有后续模型区间，因此其中那个有记录的睡眠值不充当当前时间更新区间的起点暴露；争议主要涉及两个存活结局终点的真实性，属于共享的结局来源限制，而不只是时间更新睡眠测量的问题。", "",
  "## 下一步边界与未解决方法问题", "",
  "1. CESD-9 扩展调整需单独冻结缺失处理。现有 MI20 未插补 CESD-9，不得把其缺失设为 0 或声称已纳入原 MI 模型。完整 CESD-9 子集比较也需说明选择与调整差异。",
  "2. W2 动态分析可先以已知 W2 状态且有后续信息者明确数据分母；保留部分状态者需要明确初始状态分布/条件化方案并单独审核，不能只改时间零点就复用 W1 条件模型。这里只审计沿用主 W1 队列的人群，尚未评价更广的 W2 新入组人群。原 MI20 的辅助变量针对 W1 主暴露及全程结局构建，不能默认充分支持 W2 新样本、动态分类或新增 CESD-9 调整；需要重新审查或另备匹配分析目标的插补数据。",
  "3. 时长的具体分组/非线性参数、二分类与四分类的目标比较族、多重比较校正、时间更新缺失/滞后策略尚未在本次审计中冻结。必须在查看相应效应前决定。",
  "4. 全部数字未加调查权重，不是全国患病率或可推广人群的估计。广义心脏病仍不能命名为心梗/冠心病/心衰亚型；本审计无 HR、q 或 P 值。",
  "5. 公开 CSV/说明不包含参与者 ID；根据独立复核需求，另在 internal 目录保存仅用于源记录定位的受控 RDS，不用于公开报告或共享。未改主模型、签名输入、既有结果或 HRS/ELSA。", "",
  "## 分模块门控", "",
  "W1 次要指标和 W2 landmark 均为 READY_FOR_DATA_PREPARATION：允许下一步准备相应数据，绝不是立即可拟合或分析已完成。时间更新暴露为 REVIEW_REQUIRED：上述两条 W4 状态冲突（其中一条睡眠有效）需额外判断；不能因为已知来源就把不一致检查改为通过。固定主模型仍执行原死亡审计分支，本脚本不作追溯改写。", "",
  sprintf("复核：硬性有效性 %d/%d 项通过；另有 %d 项来源不一致复核检查未通过并保留 REVIEW_REQUIRED。nature-statistics 用于明确人/记录分母、缺失口径和因果/显著性边界；软件 R %s、data.table %s、GlobalAgingR %s。",
    sum(qc[severity == "HARD_VALIDITY", pass]), nrow(qc[severity == "HARD_VALIDITY"]), sum(!qc[severity == "REVIEW_REQUIRED", pass]),
    as.character(getRversion()), as.character(packageVersion("data.table")), as.character(packageVersion("GlobalAgingR"))))
writeLines(lines, file.path(dest, "CHARLS_次要分析可用性审计.md"), useBytes = TRUE)
print(basecov[group == "overall"])
print(flow)
print(pattern_counts)
message("Secondary readiness audit complete: hard checks ", sum(qc[severity == "HARD_VALIDITY", pass]), "/", nrow(qc[severity == "HARD_VALIDITY"]),
  "; time-updated status REVIEW_REQUIRED; no models or imputation performed.")
