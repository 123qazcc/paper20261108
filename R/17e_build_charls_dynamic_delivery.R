#!/usr/bin/env Rscript
# R-only Nature-style delivery for CHARLS stage 3 dynamic sleep patterns.
options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 1L, args[[1]] == "build")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(svglite)
  library(ragg); library(jsonlite); library(digest)
})
setDTthreads(1L)

root <- "outputs/charls_dynamic_sleep_20260920"
dest <- file.path(root, "stage3_delivery"); figdir <- file.path(dest, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
control <- fromJSON("analysis/charls_stage3_dynamic_sleep_control.json", simplifyVector = TRUE)
review <- fromJSON(file.path(root, "model_qc/dynamic_independent_pool_review.json"), simplifyVector = TRUE)
pooled <- fromJSON(file.path(root, "model_qc/pooled_status.json"), simplifyVector = TRUE)
stopifnot(identical(review$status, "COMPUTATIONAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_REVIEW"),
  review$checks_total == review$checks_passed, identical(pooled$status, "POOLED_AWAITING_INDEPENDENT_REVIEW_AND_FIGURES"),
  pooled$accepted_fits == 40L)

target <- fread(file.path(root, "results/four_target_dynamic_contrasts_MI20.csv"))
omnibus <- fread(file.path(root, "results/four_target_pattern_omnibus_MI20.csv"))
direct <- fread(file.path(root, "results/heart_stroke_dynamic_direct_MI20.csv"))
joint <- fread(file.path(root, "results/joint_tests_MI20.csv"))
change <- fread(file.path(root, "results/matched_CORE_BOTH_change_descriptive.csv"))
flow <- fread(file.path(root, "tables/dynamic_landmark_flow.csv"))
patterns <- fread(file.path(root, "tables/dynamic_pattern_distribution.csv"))
support <- fread(file.path(root, "tables/four_target_observed_known_pairs_by_pattern.csv"))
stopifnot(nrow(target) == 32L, nrow(omnibus) == 8L, nrow(direct) == 16L, nrow(joint) == 10L,
  all(is.finite(target$HR)), all(target$lower > 0), all(target$upper >= target$HR),
  all(is.finite(direct$HR)), all(direct$lower > 0), sum(patterns$people) == control$expected_people)

write_csv <- function(x, p) { fwrite(x, p, bom = TRUE); stopifnot(file.exists(p), file.info(p)$size > 0) }
write_csv(flow, file.path(dest, "table1_dynamic_landmark_flow.csv"))
write_csv(patterns, file.path(dest, "table2_dynamic_pattern_distribution.csv"))
write_csv(target[role == "primary"], file.path(dest, "table3_primary_dynamic_four_paths_CORE_BOTH.csv"))
write_csv(target[role == "contextual"], file.path(dest, "table4_contextual_dynamic_four_paths_CORE_BOTH.csv"))
write_csv(omnibus, file.path(dest, "table5_four_level_pattern_omnibus_CORE_BOTH.csv"))
write_csv(direct, file.path(dest, "table6_heart_stroke_direct_dynamic_CORE_BOTH.csv"))
write_csv(joint, file.path(dest, "table7_joint_tests_CORE_BOTH.csv"))
write_csv(change, file.path(dest, "table8_CORE_BOTH_change_descriptive.csv"))
write_csv(support, file.path(dest, "table9_observed_known_pair_support.csv"))

pattern_labels <- c(persistent_low = "Persistent low", onset = "Onset",
  remission = "Remission", persistent_high = "Persistent high")
comparison_labels <- c(onset_vs_persistent_low = "Onset vs persistent low",
  persistent_high_vs_remission = "Persistent high vs remission",
  remission_vs_persistent_low = "Remission vs persistent low",
  persistent_high_vs_persistent_low = "Persistent high vs persistent low")
path_labels <- c(`0-2` = "No disease → heart disease only",
  `0-4` = "No disease → stroke only",
  `1-3` = "Hypertension only → hypertension + heart disease",
  `1-5` = "Hypertension only → hypertension + stroke")
origin_labels <- c(`No disease` = "No target disease at origin", `HTN only` = "Hypertension only at origin")
model_labels <- c(CORE = "CORE", BOTH = "BOTH: CESD-9 + duration")
palette <- c(neutral_dark = "#424A52", neutral_mid = "#77828C", neutral_light = "#D9DEE3",
  signal_blue = "#2F6F9F", signal_teal = "#4C9085", accent_red = "#B64C3C", accent_orange = "#CE8629")
pattern_cols <- c(persistent_low = "#789C98", onset = palette[["signal_blue"]],
  remission = palette[["accent_orange"]], persistent_high = palette[["accent_red"]])
comparison_cols <- c(onset_vs_persistent_low = palette[["signal_blue"]],
  persistent_high_vs_remission = palette[["accent_red"]],
  remission_vs_persistent_low = palette[["accent_orange"]],
  persistent_high_vs_persistent_low = palette[["neutral_mid"]])
model_cols <- c(CORE = palette[["neutral_dark"]], BOTH = palette[["signal_blue"]])

theme_pub <- function(base_size = 6.7) theme_classic(base_family = "Arial", base_size = base_size) +
  theme(axis.line = element_line(linewidth = .32, colour = "#222222"),
    axis.ticks = element_line(linewidth = .32, colour = "#222222"),
    axis.text = element_text(size = base_size - .3, colour = "#222222"),
    axis.title = element_text(size = base_size), strip.background = element_rect(fill = "#F1F3F5", colour = NA),
    strip.text = element_text(size = base_size - .1, face = "bold", colour = "#27313B"),
    legend.title = element_blank(), legend.text = element_text(size = base_size - .4),
    plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = base_size - .1, colour = "#4C5661", hjust = 0),
    plot.caption = element_text(size = 5.4, colour = "#4C5661", hjust = 0, lineheight = 1.12),
    panel.grid.major.x = element_line(colour = "#E6E9EC", linewidth = .25),
    panel.grid.minor = element_blank(), plot.margin = margin(5, 6, 5, 5))

log_scale <- function(d, n = 5L) {
  values <- c(d$lower, d$upper, 1); lr <- range(log(values[is.finite(values) & values > 0]))
  pad <- max(diff(lr) * .08, .045); lims <- exp(lr + c(-pad, pad))
  br <- sort(unique(c(1, exp(pretty(log(lims), n = n))))); br <- br[br >= lims[1] & br <= lims[2]]
  list(limits = lims, breaks = br)
}
hr_axis <- function(s) scale_x_log10(limits = s$limits, breaks = s$breaks,
  labels = function(x) formatC(x, digits = 3, format = "g"))

# Figure contract: the main composite asks whether onset/persistence show selective
# associations across all four prespecified paths and whether direct heart-stroke
# contrasts support disease specificity. Panel a establishes the selected landmark
# cohort; panel b is the hero evidence; panel c is the direct-comparison safeguard.
patterns[, pattern := factor(pattern, levels = names(pattern_labels), labels = unname(pattern_labels))]
p_a <- ggplot(patterns, aes(x = people, y = pattern, fill = pattern)) +
  geom_col(width = .62, colour = NA) +
  geom_text(aes(label = sprintf("%s  (%.1f%%)", format(people, big.mark = ","), percent)),
    hjust = -0.06, size = 2.25, family = "Arial", colour = "#27313B") +
  scale_fill_manual(values = setNames(pattern_cols, unname(pattern_labels)), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, .28)), labels = scales::comma) +
  labs(title = "W2 landmark cohort", subtitle = "Observed W1–W2 pattern; n = 9,397", x = "Participants", y = NULL) +
  theme_pub() + theme(panel.grid.major.x = element_blank())

primary_core <- target[model == "CORE" & role == "primary"]
primary_core[, `:=`(path = factor(path_labels[transition], levels = rev(unname(path_labels))),
  comparison_label = factor(comparison_labels[comparison], levels = unname(comparison_labels[c("onset_vs_persistent_low", "persistent_high_vs_remission")])))]
s_b <- log_scale(primary_core)
p_b <- ggplot(primary_core, aes(x = HR, y = path, colour = comparison_label, shape = comparison_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .32, colour = "#6B737B") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = .18, linewidth = .48,
    position = position_dodge(width = .48)) +
  geom_point(size = 2.05, position = position_dodge(width = .48)) +
  scale_colour_manual(values = unname(comparison_cols[c("onset_vs_persistent_low", "persistent_high_vs_remission")])) +
  scale_shape_manual(values = c(16, 15)) + hr_axis(s_b) +
  labs(title = "Primary dynamic comparisons", subtitle = "CORE-adjusted transition HR (95% CI)",
    x = "Hazard ratio", y = NULL) + theme_pub() +
  theme(legend.position = "top", axis.text.y = element_text(size = 5.7))

direct_core <- direct[model == "CORE" & role == "primary"]
direct_core[, `:=`(origin_label = factor(origin_labels[origin], levels = rev(unname(origin_labels))),
  comparison_label = factor(comparison_labels[comparison], levels = unname(comparison_labels[c("onset_vs_persistent_low", "persistent_high_vs_remission")])))]
s_c <- log_scale(direct_core)
p_c <- ggplot(direct_core, aes(x = HR, y = origin_label, colour = comparison_label, shape = comparison_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .32, colour = "#6B737B") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = .16, linewidth = .48,
    position = position_dodge(width = .44)) +
  geom_point(size = 2.05, position = position_dodge(width = .44)) +
  scale_colour_manual(values = unname(comparison_cols[c("onset_vs_persistent_low", "persistent_high_vs_remission")])) +
  scale_shape_manual(values = c(16, 15)) + hr_axis(s_c) +
  labs(title = "Heart versus stroke", subtitle = "Heart-transition HR / stroke-transition HR",
    x = "Ratio of hazard ratios", y = NULL) + theme_pub() +
  theme(legend.position = "none", axis.text.y = element_text(size = 5.9))

design <- "
ABB
ACC
"
fig1 <- p_a + p_b + p_c + plot_layout(design = design, widths = c(.9, 1, 1), heights = c(1.12, .88), guides = "collect") +
  plot_annotation(tag_levels = "a", title = "Dynamic restless-sleep patterns and disease transitions",
    caption = paste(strwrap("W1–W2 patterns require two observed sleep assessments; risk time begins at W2. Points and bars are Rubin-pooled HRs and pointwise 95% finite-MI t intervals. Primary comparisons preserve W1 sleep-frequency background. Direct ratios use full within-model covariance and are not incidence-rate ratios.", width = 132), collapse = "\n")) &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "Arial"),
    plot.title = element_text(size = 8.5, face = "bold", family = "Arial"),
    plot.caption = element_text(size = 5.4, colour = "#4C5661", family = "Arial", hjust = 0))

all_plot <- copy(target)
all_plot[, `:=`(path = factor(path_labels[transition], levels = rev(unname(path_labels))),
  comparison_label = factor(comparison_labels[comparison], levels = unname(comparison_labels[names(comparison_labels)])),
  model_label = factor(model_labels[model], levels = unname(model_labels)),
  role_label = factor(role, levels = c("primary", "contextual")))]
s2 <- log_scale(all_plot)
p2 <- ggplot(all_plot, aes(x = HR, y = path, colour = role_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .28, colour = "#6B737B") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = .14, linewidth = .42) +
  geom_point(size = 1.7, aes(shape = role_label)) +
  facet_grid(model_label ~ comparison_label) +
  scale_colour_manual(values = c(primary = palette[["signal_blue"]], contextual = palette[["neutral_mid"]]),
    labels = c(primary = "Primary within-W1 comparison", contextual = "Contextual comparison")) +
  scale_shape_manual(values = c(primary = 16, contextual = 15),
    labels = c(primary = "Primary within-W1 comparison", contextual = "Contextual comparison")) +
  hr_axis(s2) + labs(title = "Complete four-pattern results across all prespecified disease paths",
    subtitle = "CORE is primary; BOTH additionally adjusts for W1 CESD-9 and sleep-duration spline",
    x = "Transition hazard ratio", y = NULL,
    caption = paste(strwrap("Each facet displays all four prespecified paths; no transition is omitted by statistical significance. Intervals are pointwise 95% finite-MI t intervals. BH-adjusted values, 3-df per-edge omnibus tests and joint tests are in the companion tables.", width = 145), collapse = "\n")) +
  theme_pub(6.25) + theme(legend.position = "top", axis.text.y = element_text(size = 5.15),
    strip.text.x = element_text(size = 5.5), strip.text.y = element_text(size = 5.7))

direct_plot <- copy(direct)
direct_plot[, `:=`(origin_label = factor(origin_labels[origin], levels = rev(unname(origin_labels))),
  comparison_label = factor(comparison_labels[comparison], levels = unname(comparison_labels[names(comparison_labels)])),
  model_label = factor(model_labels[model], levels = unname(model_labels)),
  role_label = factor(role, levels = c("primary", "contextual")))]
s3 <- log_scale(direct_plot)
p3 <- ggplot(direct_plot, aes(x = HR, y = origin_label, colour = role_label, shape = role_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .28, colour = "#6B737B") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = .14, linewidth = .42) +
  geom_point(size = 1.7) + facet_grid(model_label ~ comparison_label) +
  scale_colour_manual(values = c(primary = palette[["accent_red"]], contextual = palette[["neutral_mid"]]),
    labels = c(primary = "Primary within-W1 comparison", contextual = "Contextual comparison")) +
  scale_shape_manual(values = c(primary = 16, contextual = 15),
    labels = c(primary = "Primary within-W1 comparison", contextual = "Contextual comparison")) +
  hr_axis(s3) + labs(title = "Direct comparison of heart-disease and stroke transition associations",
    subtitle = "Ratio >1 indicates a larger sleep-pattern HR for the heart-disease transition at the same origin",
    x = "Heart-transition HR / stroke-transition HR", y = NULL,
    caption = paste(strwrap("Ratios and intervals retain full covariance between the two transition coefficients within each model. These are ratios of association HRs, not disease incidence-rate ratios; failure to reject 1 is not evidence of equivalence.", width = 145), collapse = "\n")) +
  theme_pub(6.25) + theme(legend.position = "top", axis.text.y = element_text(size = 5.4),
    strip.text.x = element_text(size = 5.5), strip.text.y = element_text(size = 5.7))

skill_root <- Sys.getenv("NATURE_FIGURE_SKILL_ROOT", path.expand("private/render_tools/nature-figure"))
source(file.path(skill_root, "scripts/panel_alignment.R"))
python_bin <- Sys.getenv("NATURE_FIGURE_PYTHON", Sys.which("python"))
stopifnot(nzchar(python_bin), file.exists(python_bin))
audit_script <- file.path(skill_root, "scripts/audit_panel_alignment.py")
sha <- function(p) digest(file = p, algo = "sha256")
export_plot <- function(plot, stem, width_mm, height_mm, multipanel = FALSE) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  if (multipanel) {
    require_patchwork_panel_alignment(plot,
      manifest_path = paste0(stem, ".alignment-layout.json"),
      report_path = paste0(stem, ".alignment.json"),
      overlay_svg = paste0(stem, ".alignment.svg"),
      width_in = w, height_in = h, tolerance_pt = 1.5, gutter_tolerance_pt = 1.5,
      strict = TRUE, audit_script = audit_script, python = python_bin)
  } else {
    write_json(list(verdict = "NOT APPLICABLE", reason = "Single ggplot object; facet geometry is managed within one rendered plot area.",
      backend = "R ggplot2", figure = basename(stem)), paste0(stem, ".alignment.json"), pretty = TRUE, auto_unbox = TRUE)
  }
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = "Arial"); print(plot); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in", res = 600); print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(plot); dev.off()
  files <- paste0(stem, c(".pdf", ".svg", ".png", ".tiff")); stopifnot(all(file.info(files)$size > 0))
  write_json(list(backend = "R", width_mm = width_mm, height_mm = height_mm, dpi_raster = 600,
    files = basename(files), sha256 = as.list(setNames(vapply(files, sha, ""), basename(files)))),
    paste0(stem, ".export.json"), pretty = TRUE, auto_unbox = TRUE)
}

s1 <- file.path(figdir, "figure1_dynamic_primary_composite")
s2stem <- file.path(figdir, "figure2_all_four_pattern_contrasts")
s3stem <- file.path(figdir, "figure3_heart_stroke_direct_comparisons")
export_plot(fig1, s1, 183, 124, multipanel = TRUE)
export_plot(p2, s2stem, 183, 137, multipanel = FALSE)
export_plot(p3, s3stem, 183, 105, multipanel = FALSE)
write_csv(patterns, paste0(s1, ".panel_a.source_data.csv"))
write_csv(primary_core, paste0(s1, ".panel_b.source_data.csv"))
write_csv(direct_core, paste0(s1, ".panel_c.source_data.csv"))
write_csv(all_plot, paste0(s2stem, ".source_data.csv"))
write_csv(direct_plot, paste0(s3stem, ".source_data.csv"))

fmtp <- function(x) ifelse(x < .0001, "<0.0001", formatC(x, digits = 4, format = "fg", flag = "#"))
fmthr <- function(hr, lo, hi) sprintf("%.2f (%.2f–%.2f)", hr, lo, hi)
core_joint <- joint[model == "CORE"]
both_joint <- joint[model == "BOTH"]
supported_rows <- function(model) target[model == ..model & role == "primary" & p_BH_family < .05]
core_supported <- supported_rows("CORE"); both_supported <- supported_rows("BOTH")
describe_supported <- function(x) {
  if (!nrow(x)) return("无主要边级比较在本家族 BH(8)<0.05。")
  paste(vapply(seq_len(nrow(x)), function(i) paste0(comparison_labels[[x$comparison[i]]], "；",
    path_labels[[x$transition[i]]], "：HR ", fmthr(x$HR[i], x$lower[i], x$upper[i]),
    "，BH(8) P=", fmtp(x$p_BH_family[i])), ""), collapse = "；")
}
direct_core_sig <- direct[model == "CORE" & role == "primary" & p_BH_family < .05]
direct_both_sig <- direct[model == "BOTH" & role == "primary" & p_BH_family < .05]
describe_direct <- function(x) if (!nrow(x)) "主要直接比较均未通过 BH(4)。" else paste(vapply(seq_len(nrow(x)), function(i)
  paste0(comparison_labels[[x$comparison[i]]], "；", origin_labels[[x$origin[i]]], "：HR比值 ",
    fmthr(x$HR[i], x$lower[i], x$upper[i]), "，BH(4) P=", fmtp(x$p_BH_family[i])), ""), collapse = "；")

report <- c("# CHARLS 第三研究环节：动态睡眠模式结果与复盘", "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "。本报告在 40/40 模型接受且独立数值核验通过后生成。"), "",
  "## 本环节回答的问题", "",
  "W1–W2 睡眠不安由低频转为高频、或持续高频，是否与 W2 后四条预设疾病转移相关；同一起点下心脏病与卒中关联是否存在直接差异。动态暴露不是临床失眠诊断。", "",
  "## 样本与时间顺序", "",
  paste0("W1 主队列 12,036 人，经 W2 真实存活访谈、两波睡眠实测、W2 三疾病组件完整及 W2 后可靠信息筛选后，固定 landmark 样本为 ",
    format(control$expected_people, big.mark = ","), " 人、", format(control$expected_panel_records, big.mark = ","), " 条 W2 后模型记录。风险时间从 W2 开始；W1–W2 仅定义过去模式。"),
  paste0("四类人数：", paste0(unname(pattern_labels), " ", patterns$people[match(unname(pattern_labels), as.character(patterns$pattern))], collapse = "；"), "。"),
  "动态模式和 W2 初始疾病状态均未插补；20份MICE仅处理调整变量缺失。", "",
  "## 主要动态比较", "",
  paste0("CORE 八项主要比较联合 Wald P=", fmtp(core_joint[family == "primary_target", p]), "。", describe_supported(core_supported)),
  paste0("BOTH 八项主要比较联合 Wald P=", fmtp(both_joint[family == "primary_target", p]), "。", describe_supported(both_supported)), "",
  "主要比较为‘新出现高频 vs 持续低频’及‘持续高频 vs 缓解’，分别保持 W1 低频和 W1 高频背景。缓解 vs 持续低频、持续高频 vs 持续低频完整报告为情境性比较，不能解释为纯改变效应。", "",
  "## 心脏病与卒中直接比较", "",
  paste0("CORE 四项主要直接比较联合 Wald P=", fmtp(core_joint[family == "primary_direct", p]), "。", describe_direct(direct_core_sig)),
  paste0("BOTH 四项主要直接比较联合 Wald P=", fmtp(both_joint[family == "primary_direct", p]), "。", describe_direct(direct_both_sig)),
  "直接量是同一模型中两个睡眠关联 HR 的比值，不是心脏病与卒中的发生率比；未显著不等于等效。", "",
  "## 与前两环节的衔接", "",
  "CORE 是与 HRS/ELSA 动态分析解释一致的主规格；BOTH 加入 W1 CESD-9 与睡眠时长样条，仅作预定稳健性检查。两模型估计相关，CORE→BOTH 的变化不作中介比例或独立变化检验。", "",
  "## 结论边界", "",
  "本结果估计的是条件于参加 W2、两波睡眠实测、W2 状态完整且有后续信息者的转移强度关联。选择进入 landmark、W2 睡眠与疾病同次测量、访谈间发病时间未知和未测混杂均可能影响结果；不得作因果、临床诊断或全国代表性解释。广义心脏病不能细化命名为心梗、冠心病或心衰。", "",
  "## 交付与下一步", "",
  "已交付样本流程、完整四模式四路径结果、主要/情境性比较、每边 omnibus、直接心脏病－卒中比较、联合检验、CORE/BOTH 描述性变化、三幅图及源数据。图件仍需最终尺寸人工检查及PDF文字/碰撞QA记录。", "",
  "本阶段不估计血压背景或纵向机制效应。")
writeLines(report, file.path(dest, "CHARLS第三环节_结果与复盘.md"), useBytes = TRUE)

writeLines(c("# Figure contract", "", "Core conclusion: Dynamic-pattern associations are displayed across every prespecified disease path without significance-based omission, with direct heart-versus-stroke tests and adjustment robustness.",
  "Results-level question: Do onset or persistence of high-frequency restless sleep predict selective W2-forward disease transitions, and are heart and stroke associations directly different?",
  "Figure archetype: quantitative grid with one hero forest panel.", "Backend: R only.", "Final size: 183 mm wide; editable PDF/SVG plus 600 dpi PNG/TIFF.",
  "Evidence hierarchy: panel b primary path estimates; panel c direct disease-specific contrasts; panel a selected landmark composition.",
  "Statistics: Rubin-pooled HR and pointwise 95% finite-MI t intervals; BH and joint Wald in companion tables.",
  "Reviewer risks: W2 landmark selection, reverse causality, interval observation, multiplicity, direct-comparison interpretation, non-equivalence of null results."),
  file.path(dest, "figure_contract.md"), useBytes = TRUE)
write_json(list(status = "DRAFT_REQUIRES_HUMAN_FIGURE_AND_SCIENTIFIC_REVIEW", backend = "R",
  accepted_fits = 40L, independent_checks = review$checks_total,
  figures = c(s1, s2stem, s3stem), human_figure_review = "PENDING", human_scientific_review = "PENDING",
  automatic_next_stage = FALSE), file.path(dest, "delivery_status.json"), pretty = TRUE, auto_unbox = TRUE)
message("STAGE3_DELIVERY_DRAFT_READY: three figures, nine tables, Chinese recap; hard stop before blood-pressure stage.")
