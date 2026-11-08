#!/usr/bin/env Rscript

# Plot-only redraw for the frozen 11b symptom-domain results.  This script reads
# only saved result CSVs and does not load, alter, or refit any msm object.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(jsonlite)
})

options(stringsAsFactors = FALSE, width = 180)

output_root <- file.path("outputs", "hrs_elsa_multistate_nextsteps_20260909")
result_root <- file.path(output_root, "results", "symptom_domains")
figure_root <- file.path(output_root, "figures", "symptom_domains")
source_root <- file.path(figure_root, "source_data")
qc_root <- file.path(output_root, "qc", "symptom_domains")
for (path in c(figure_root, source_root, qc_root)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

estimate_path <- file.path(result_root, "symptom_domain_joint_estimates.csv")
pairwise_path <- file.path(result_root, "symptom_domain_pairwise_contrasts.csv")
if (!file.exists(estimate_path) || !file.exists(pairwise_path)) stop("Frozen 11b result CSVs are missing")
estimates <- fread(estimate_path)
pairwise <- fread(pairwise_path)

analysis_display <- c(
  hrs_5state_joint_domains = "HRS five-state model",
  elsa_4state_w4_w10_joint_domains = "ELSA W4-W10 four-state model",
  hrs_4state_death_censored_joint_domains = "HRS four-state death-censored sensitivity"
)
domain_labels <- c(
  insomnia_onset_sync = "Sleep-onset difficulty",
  insomnia_maintenance_sync = "Sleep-maintenance difficulty",
  insomnia_nonrestorative_sync = "Nonrestorative sleep"
)
transition_display <- c(
  "S0->S1" = "No disease -> hypertension only",
  "S0->S2" = "No disease -> one CVD without hypertension",
  "S0->S4" = "No disease -> death",
  "S1->S3" = "Hypertension only -> multimorbidity",
  "S1->S4" = "Hypertension only -> death",
  "S2->S3" = "One CVD without hypertension -> multimorbidity",
  "S2->S4" = "One CVD without hypertension -> death",
  "S3->S4" = "Multimorbidity -> death"
)
pair_labels <- c(
  onset_vs_maintenance = "Onset / maintenance",
  onset_vs_nonrestorative = "Onset / nonrestorative",
  maintenance_vs_nonrestorative = "Maintenance / nonrestorative"
)
palette <- c(
  insomnia_onset_sync = "#2F6B9A",
  insomnia_maintenance_sync = "#D08726",
  insomnia_nonrestorative_sync = "#2A8C82"
)

association_limits <- c(0.70, 1.65)
association_breaks <- c(0.7, 0.8, 1.0, 1.2, 1.5)
association_labels <- c("0.7", "0.8", "1.0", "1.2", "1.5")
contrast_limits <- c(0.55, 2.05)
contrast_breaks <- c(0.6, 0.8, 1.0, 1.2, 1.5, 2.0)
contrast_labels <- c("0.6", "0.8", "1.0", "1.2", "1.5", "2.0")
death_association_limits <- c(0.68, 1.78)
death_association_breaks <- c(0.7, 0.8, 1.0, 1.2, 1.5)
death_association_labels <- c("0.7", "0.8", "1.0", "1.2", "1.5")
death_contrast_limits <- c(0.58, 2.05)
death_contrast_breaks <- contrast_breaks
death_contrast_labels <- contrast_labels

main_ids <- c("hrs_5state_joint_domains", "elsa_4state_w4_w10_joint_domains")
main_est <- estimates[analysis_id %in% main_ids & death_transition == 0L]
main_pw <- pairwise[analysis_id %in% main_ids & family == "common_nonfatal"]
sens_est <- estimates[
  analysis_id %in% c("hrs_5state_joint_domains", "hrs_4state_death_censored_joint_domains") &
    death_transition == 0L
]
death_est <- estimates[analysis_id == "hrs_5state_joint_domains" & death_transition == 1L]
death_pw <- pairwise[analysis_id == "hrs_5state_joint_domains" & family == "hrs_death_exploratory"]

assert_interval_inside <- function(x, limits, label) {
  if (!nrow(x) || any(!is.finite(x$conf_low)) || any(!is.finite(x$conf_high)) ||
      min(x$conf_low) < limits[1L] || max(x$conf_high) > limits[2L]) {
    stop(label, " confidence intervals do not fit the frozen display limits")
  }
}
assert_interval_inside(rbindlist(list(main_est, sens_est), fill = TRUE), association_limits, "Shared association panels")
assert_interval_inside(main_pw, contrast_limits, "Shared pairwise panels")
assert_interval_inside(death_est, death_association_limits, "Exploratory death association panel")
assert_interval_inside(death_pw, death_contrast_limits, "Exploratory death contrast panel")

scale_contract <- rbindlist(list(
  data.table(
    scale_group = "shared_association", panels = "A1 | A2 | C",
    estimand = "conditional transition-intensity ratio",
    actual_min_ci = min(c(main_est$conf_low, sens_est$conf_low)),
    actual_max_ci = max(c(main_est$conf_high, sens_est$conf_high)),
    display_min = association_limits[1L], display_max = association_limits[2L],
    breaks = paste(association_breaks, collapse = " | "), reference = 1
  ),
  data.table(
    scale_group = "shared_pairwise", panels = "B1 | B2",
    estimand = "ratio of transition-intensity ratios",
    actual_min_ci = min(main_pw$conf_low), actual_max_ci = max(main_pw$conf_high),
    display_min = contrast_limits[1L], display_max = contrast_limits[2L],
    breaks = paste(contrast_breaks, collapse = " | "), reference = 1
  ),
  data.table(
    scale_group = "death_association", panels = "exploratory A",
    estimand = "conditional transition-intensity ratio",
    actual_min_ci = min(death_est$conf_low), actual_max_ci = max(death_est$conf_high),
    display_min = death_association_limits[1L], display_max = death_association_limits[2L],
    breaks = paste(death_association_breaks, collapse = " | "), reference = 1
  ),
  data.table(
    scale_group = "death_pairwise", panels = "exploratory B",
    estimand = "ratio of transition-intensity ratios",
    actual_min_ci = min(death_pw$conf_low), actual_max_ci = max(death_pw$conf_high),
    display_min = death_contrast_limits[1L], display_max = death_contrast_limits[2L],
    breaks = paste(death_contrast_breaks, collapse = " | "), reference = 1
  )
))
fwrite(scale_contract, file.path(qc_root, "symptom_domain_figure_scale_contract.csv"), na = "")

theme_joint <- function() {
  theme_classic(base_size = 7.2, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = .35, colour = "#222222"),
      axis.ticks = element_line(linewidth = .35, colour = "#222222"),
      strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7),
      legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
      legend.justification = "left", legend.title = element_text(face = "bold", size = 6.5),
      legend.text = element_text(size = 6.5), legend.spacing.y = grid::unit(1, "mm"),
      plot.title = element_text(face = "bold", size = 8.2),
      plot.subtitle = element_text(size = 6.5, colour = "#444444"),
      plot.tag = element_text(face = "bold", size = 9), panel.spacing = grid::unit(4, "mm")
    )
}

render_raster_staged <- function(plot, target, format = c("png", "tiff"), width, height, res = 600) {
  format <- match.arg(format)
  staged <- tempfile(fileext = paste0(".", format))
  on.exit(unlink(staged), add = TRUE)
  if (format == "png") {
    ragg::agg_png(staged, width = width, height = height, units = "in", res = res)
  } else {
    ragg::agg_tiff(staged, width = width, height = height, units = "in", res = res, compression = "lzw")
  }
  print(plot); grDevices::dev.off()
  if (!file.copy(staged, target, overwrite = TRUE)) stop("Unable to copy staged raster: ", target)
}

export_figure_set <- function(plot, base, width, height, res = 600) {
  svglite::svglite(paste0(base, ".svg"), width = width, height = height)
  print(plot); grDevices::dev.off()
  staged_pdf <- tempfile(fileext = ".pdf")
  on.exit(unlink(staged_pdf), add = TRUE)
  grDevices::pdf(
    staged_pdf, width = width, height = height,
    family = "Helvetica", useDingbats = FALSE
  )
  print(plot); grDevices::dev.off()
  if (!file.copy(staged_pdf, paste0(base, ".pdf"), overwrite = TRUE)) stop("Unable to copy staged PDF")
  render_raster_staged(plot, paste0(base, ".tiff"), "tiff", width, height, res)
  render_raster_staged(plot, paste0(base, ".png"), "png", width, height, res)
}

# ggplot2 4.0 leaves patchwork null tracks unresolved when converted outside
# their parent layout; explicitly allocate the remaining device space so the
# strict backend-neutral alignment audit receives physical panel rectangles.
write_patchwork_panel_layout_resolved <- function(
  plot, manifest_path, width_in, height_in, panel_ids = NULL,
  row_groups = NULL, column_groups = NULL, exemptions = list()
) {
  grob <- patchwork::patchworkGrob(plot)
  panel_rows <- grob$layout[grepl("^panel(?:-[0-9]+)?$", grob$layout$name, perl = TRUE), , drop = FALSE]
  panel_rows <- panel_rows[order(panel_rows$t, panel_rows$l), , drop = FALSE]
  if (nrow(panel_rows) < 2L) stop("At least two flat patchwork panels are required")
  if (is.null(panel_ids)) panel_ids <- letters[seq_len(nrow(panel_rows))]
  probe <- tempfile(fileext = ".pdf")
  grDevices::pdf(probe, width = width_in, height = height_in, useDingbats = FALSE)
  open <- TRUE
  on.exit({ if (open) grDevices::dev.off(); unlink(probe) }, add = TRUE)
  grid::grid.newpage(); grid::grid.draw(grob); grid::grid.force()
  resolve_tracks <- function(units, total_pt, axis = c("width", "height")) {
    axis <- match.arg(axis)
    type <- grid::unitType(units)
    value <- if (axis == "width") grid::convertWidth(units, "pt", valueOnly = TRUE) else grid::convertHeight(units, "pt", valueOnly = TRUE)
    null <- type == "null"
    if (any(null)) {
      weight <- as.numeric(units[null])
      remaining <- total_pt - sum(value[!null])
      if (!is.finite(remaining) || remaining <= 0 || sum(weight) <= 0) stop("Unable to resolve null tracks")
      value[null] <- remaining * weight / sum(weight)
    }
    value
  }
  widths <- resolve_tracks(grob$widths, width_in * 72, "width")
  heights <- resolve_tracks(grob$heights, height_in * 72, "height")
  xedge <- c(0, cumsum(widths)); topedge <- c(0, cumsum(heights)); total_height <- sum(heights)
  panels <- lapply(seq_len(nrow(panel_rows)), function(i) {
    z <- panel_rows[i, ]
    list(
      id = as.character(panel_ids[i]),
      bbox_pt = unname(c(xedge[z$l], total_height - topedge[z$b + 1L], xedge[z$r + 1L], total_height - topedge[z$t])),
      grid_id = "patchwork-grid-1", row_start = as.integer(z$t - 1L), row_stop = as.integer(z$b),
      col_start = as.integer(z$l - 1L), col_stop = as.integer(z$r)
    )
  })
  manifest <- list(
    schema_version = 1L, backend = "r-patchwork-resolved-null-tracks",
    figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
    panels = panels, exemptions = exemptions
  )
  if (!is.null(row_groups)) manifest$row_groups <- lapply(row_groups, as.character)
  if (!is.null(column_groups)) manifest$column_groups <- lapply(column_groups, as.character)
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  grDevices::dev.off(); open <- FALSE; unlink(probe)
}

alignment_script <- "private/render_tools/nature-figure/scripts/panel_alignment.R"
alignment_audit <- "private/render_tools/nature-figure/scripts/audit_panel_alignment.py"
alignment_python <- "python"
if (!file.exists(alignment_script) || !file.exists(alignment_audit) || !file.exists(alignment_python)) {
  stop("Required strict alignment-audit runtime is missing")
}
source(alignment_script, local = environment())
write_patchwork_panel_layout <- write_patchwork_panel_layout_resolved

transition_order <- c("S0->S1", "S0->S2", "S1->S3", "S2->S3")
main_est[, transition_label := factor(unname(transition_display[transition]), levels = rev(unname(transition_display[transition_order])))]
main_pw[, transition_label := factor(unname(transition_display[transition]), levels = rev(unname(transition_display[transition_order])))]
main_pw[, contrast_label := factor(unname(pair_labels[contrast]), levels = unname(pair_labels))]
main_pw[, supported := factor(
  fifelse(global_gate_pass == 1L & p_bh < .05, "Gate open; BH q < 0.05", "Not gate-supported"),
  levels = c("Gate open; BH q < 0.05", "Not gate-supported")
)]
sens_est[, transition_label := factor(unname(transition_display[transition]), levels = rev(unname(transition_display[transition_order])))]
sens_est[, model := factor(
  fifelse(analysis_id == "hrs_5state_joint_domains", "Five-state (death competing state)", "Four-state (death censored)"),
  levels = c("Five-state (death competing state)", "Four-state (death censored)")
)]

association_scale <- scale_x_log10(
  limits = association_limits, breaks = association_breaks, labels = association_labels,
  expand = expansion(mult = c(0.01, 0.01))
)
contrast_scale <- scale_x_log10(
  limits = contrast_limits, breaks = contrast_breaks, labels = contrast_labels,
  expand = expansion(mult = c(0.01, 0.01))
)

p_association <- ggplot(main_est, aes(x = transition_intensity_ratio, y = transition_label, colour = domain, shape = domain)) +
  geom_vline(xintercept = 1, linewidth = .35, linetype = "22", colour = "#777777") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0, linewidth = .45, position = position_dodge(width = .55), orientation = "y") +
  geom_point(size = 1.8, position = position_dodge(width = .55)) + association_scale +
  scale_colour_manual(values = palette, breaks = names(domain_labels), labels = unname(domain_labels), name = "Symptom domain") +
  scale_shape_manual(values = c(16, 17, 15), breaks = names(domain_labels), labels = unname(domain_labels), name = "Symptom domain") +
  labs(x = "Conditional transition-intensity ratio (95% CI)", y = NULL) + theme_joint()

p_pairwise <- ggplot(main_pw, aes(x = ratio_of_transition_intensity_ratios, y = transition_label, colour = contrast_label, shape = supported)) +
  geom_vline(xintercept = 1, linewidth = .35, linetype = "22", colour = "#777777") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0, linewidth = .45, position = position_dodge(width = .58), orientation = "y") +
  geom_point(size = 1.7, position = position_dodge(width = .58)) + contrast_scale +
  scale_colour_manual(values = c("#6A51A3", "#D95F59", "#4D7770"), name = "Direct contrast") +
  scale_shape_manual(
    values = c(16, 1), limits = c("Gate open; BH q < 0.05", "Not gate-supported"),
    drop = FALSE, name = "Gate/BH status"
  ) + labs(x = "Ratio of transition-intensity ratios (95% CI)", y = NULL) + theme_joint()

p_c <- ggplot(sens_est, aes(x = transition_intensity_ratio, y = transition_label, colour = domain, shape = model)) +
  geom_vline(xintercept = 1, linewidth = .35, linetype = "22", colour = "#777777") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0, linewidth = .4, position = position_dodge(width = .62), orientation = "y") +
  geom_point(size = 1.65, position = position_dodge(width = .62)) + association_scale +
  scale_colour_manual(values = palette, breaks = names(domain_labels), labels = unname(domain_labels), name = "Symptom domain") +
  scale_shape_manual(values = c(16, 1), name = "HRS specification") +
  guides(colour = "none") +
  labs(
    x = "Conditional transition-intensity ratio (95% CI)", y = NULL,
    title = "HRS structural-comparability sensitivity",
    subtitle = "Death-censored four-state results do not replace the five-state analysis"
  ) + theme_joint()

p_a1 <- p_association + main_est[analysis_id == "hrs_5state_joint_domains"] +
  labs(title = "HRS joint symptom-domain associations", subtitle = "Five-state model; four non-fatal transitions shown")
p_a2 <- p_association + main_est[analysis_id == "elsa_4state_w4_w10_joint_domains"] +
  labs(title = "ELSA joint symptom-domain associations", subtitle = "W4-W10 four-state model")
p_b1 <- p_pairwise + main_pw[analysis_id == "hrs_5state_joint_domains"] +
  labs(title = "HRS direct pairwise domain contrasts", subtitle = "Filled symbols: global gate open and BH q < 0.05")
p_b2 <- p_pairwise + main_pw[analysis_id == "elsa_4state_w4_w10_joint_domains"] +
  guides(shape = "none") +
  labs(title = "ELSA direct pairwise domain contrasts", subtitle = "Global domain-equality gate closed; all contrasts descriptive")

main_figure <- wrap_plots(
  p_a1, p_a2, p_b1, p_b2, p_c, ncol = 1, guides = "collect",
  heights = c(1, 1, 1, 1, .9)
) + plot_annotation(tag_levels = list(c("A1", "A2", "B1", "B2", "C"))) &
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left")
main_base <- file.path(figure_root, "fig_joint_symptom_domain_comparison")
main_width <- 183 / 25.4; main_height <- 250 / 25.4
require_patchwork_panel_alignment(
  main_figure, manifest_path = paste0(main_base, ".alignment-layout.json"),
  report_path = paste0(main_base, ".alignment.json"), overlay_svg = paste0(main_base, ".alignment.svg"),
  audit_script = alignment_audit, python = alignment_python,
  width_in = main_width, height_in = main_height, tolerance_pt = 1.5,
  gutter_tolerance_pt = 1.5, strict = TRUE
)
export_figure_set(main_figure, main_base, main_width, main_height, 600)
writeLines(c(
  "Fig. | Joint comparison of insomnia symptom domains across cardiovascular multistate transitions.",
  "(A1-A2) Conditional transition-intensity ratios (95% confidence intervals) from HRS and ELSA models simultaneously containing all three binary symptom domains. (B1-B2) Direct ratios of transition-intensity ratios for the frozen pairwise domain contrasts. Filled symbols identify contrasts supported after the cohort-level global gate and Benjamini-Hochberg adjustment; open symbols are not gate-supported. (C) HRS five-state estimates compared descriptively with the death-censored four-state sensitivity. A1, A2 and C share the same logarithmic effect scale and limits; B1 and B2 share a separate common logarithmic contrast scale, so the reference line is directly aligned within each estimand family. Estimates combine 20 imputations using Rubin's rules. This is a post-hoc secondary extension whose contrasts were fixed before the new joint models were fitted.",
  "",
  "图 | 三个失眠症状域与心血管多状态转移的联合比较。",
  "A1-A2，HRS和ELSA中同时纳入三个二元症状域的联合模型所估计的条件转移强度比及95%置信区间。B1-B2，冻结的域间两两对比；实心符号表示整体门控开启且经BH校正后得到支持，空心符号表示未获门控支持。C，HRS五状态与死亡删失四状态敏感性结果的描述性比较。A1、A2和C使用完全相同的对数效应坐标范围，B1和B2使用另一组完全相同的对数对比坐标，因而同类估计的无效线可直接对齐比较。结果按Rubin规则合并20重插补。本图属于事后次级扩展，所有对比在新联合模型拟合前冻结。"
), paste0(main_base, "_caption.md"), useBytes = TRUE)

death_order <- c("S0->S4", "S1->S4", "S2->S4", "S3->S4")
death_est[, transition_label := factor(unname(transition_display[transition]), levels = rev(unname(transition_display[death_order])))]
death_pw[, transition_label := factor(unname(transition_display[transition]), levels = rev(unname(transition_display[death_order])))]
death_pw[, contrast_label := factor(unname(pair_labels[contrast]), levels = unname(pair_labels))]
death_pw[, supported := factor(
  fifelse(global_gate_pass == 1L & p_bh < .05, "Gate open; BH q < 0.05", "Not gate-supported"),
  levels = c("Gate open; BH q < 0.05", "Not gate-supported")
)]

p_death_a <- ggplot(death_est, aes(x = transition_intensity_ratio, y = transition_label, colour = domain, shape = domain)) +
  geom_vline(xintercept = 1, linewidth = .35, linetype = "22", colour = "#777777") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0, linewidth = .45, position = position_dodge(width = .55), orientation = "y") +
  geom_point(size = 1.8, position = position_dodge(width = .55)) +
  scale_x_log10(
    limits = death_association_limits, breaks = death_association_breaks,
    labels = death_association_labels, expand = expansion(mult = c(.01, .01))
  ) +
  scale_colour_manual(values = palette, breaks = names(domain_labels), labels = unname(domain_labels), name = "Symptom domain") +
  scale_shape_manual(values = c(16, 17, 15), breaks = names(domain_labels), labels = unname(domain_labels), name = "Symptom domain") +
  labs(
    x = "Conditional transition-intensity ratio (95% CI)", y = NULL,
    title = "HRS death-transition associations (exploratory)",
    subtitle = "All three domains entered simultaneously; no comparable ELSA death analysis"
  ) + theme_joint()

p_death_b <- ggplot(death_pw, aes(x = ratio_of_transition_intensity_ratios, y = transition_label, colour = contrast_label, shape = supported)) +
  geom_vline(xintercept = 1, linewidth = .35, linetype = "22", colour = "#777777") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0, linewidth = .45, position = position_dodge(width = .58), orientation = "y") +
  geom_point(size = 1.7, position = position_dodge(width = .58)) +
  scale_x_log10(
    limits = death_contrast_limits, breaks = death_contrast_breaks,
    labels = death_contrast_labels, expand = expansion(mult = c(.01, .01))
  ) +
  scale_colour_manual(values = c("#6A51A3", "#D95F59", "#4D7770"), name = "Direct contrast") +
  scale_shape_manual(
    values = c(16, 1), limits = c("Gate open; BH q < 0.05", "Not gate-supported"),
    drop = FALSE, name = "Gate/BH status"
  ) +
  labs(
    x = "Ratio of transition-intensity ratios (95% CI)", y = NULL,
    title = "HRS direct death-edge domain contrasts (exploratory)",
    subtitle = "BH adjustment applied within the frozen 12-contrast death-edge family"
  ) + theme_joint()

death_figure <- p_death_a / p_death_b + plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left")
death_base <- file.path(figure_root, "fig_hrs_death_symptom_domain_exploratory")
death_width <- 183 / 25.4; death_height <- 180 / 25.4
require_patchwork_panel_alignment(
  death_figure, manifest_path = paste0(death_base, ".alignment-layout.json"),
  report_path = paste0(death_base, ".alignment.json"), overlay_svg = paste0(death_base, ".alignment.svg"),
  audit_script = alignment_audit, python = alignment_python,
  width_in = death_width, height_in = death_height, tolerance_pt = 1.5,
  gutter_tolerance_pt = 1.5, strict = TRUE
)
export_figure_set(death_figure, death_base, death_width, death_height, 600)
writeLines(c(
  "Supplementary Fig. | Exploratory HRS symptom-domain results for transitions to death.",
  "(A) Mutually adjusted conditional transition-intensity ratios and 95% confidence intervals for the four HRS death transitions. (B) Direct ratios of transition-intensity ratios for the three frozen pairwise contrasts on each death edge. Filled symbols passed the exploratory death-family global gate and BH correction; open symbols were not gate-supported. Both logarithmic axes show explicit ticks above one; panels use different limits because they display different estimands. Death-edge findings cannot be validated against ELSA W4-W10 because comparable mortality follow-up was unavailable.",
  "",
  "补充图 | HRS转入死亡各边的探索性症状域结果。",
  "A，HRS四条死亡边中三个症状域相互调整后的条件转移强度比及95%置信区间。B，每条死亡边上三个冻结的域间两两对比；实心符号表示通过探索性整体门控和BH校正，空心符号表示未获门控支持。两个对数坐标均明确标出大于1的刻度；因估计对象不同，A和B使用不同坐标范围。ELSA W4-W10缺乏可比死亡随访，死亡边结果不能作跨队列验证。"
), paste0(death_base, "_caption.md"), useBytes = TRUE)

fwrite(main_est, file.path(source_root, "figure_joint_domain_estimates_source.csv"), na = "")
fwrite(main_pw, file.path(source_root, "figure_pairwise_domain_contrasts_source.csv"), na = "")
fwrite(sens_est, file.path(source_root, "figure_hrs4_sensitivity_source.csv"), na = "")
fwrite(death_est, file.path(source_root, "figure_hrs_death_domain_estimates_source.csv"), na = "")
fwrite(death_pw, file.path(source_root, "figure_hrs_death_pairwise_contrasts_source.csv"), na = "")

cat(sprintf("Plot-only 11b redraw complete at %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")))
