#!/usr/bin/env Rscript
# Presentation-only Figure 6; no fitting or pooling.
# Contract: one quantitative panel asks whether four observed sleep-pattern
# levels differ collectively within each pathway and model specification.
# It displays all eight precomputed asymptotic Wald chi-square statistics,
# df = 3, raw P and BH P. These are not effect sizes or directional evidence.
# R-only export: Arial, 183 mm, editable PDF/SVG and 600 dpi PNG/TIFF.
source("R/17j_dynamic_extra_figure_helpers.R", encoding = "UTF-8")

omnibus_input <- file.path(root, "stage3_delivery/table5_four_level_pattern_omnibus_CORE_BOTH.csv")
omnibus_hash <- sha(omnibus_input)
omnibus <- fread(omnibus_input)
stopifnot(nrow(omnibus) == 8L,
  uniqueN(omnibus, by = c("model", "transition")) == 8L,
  setequal(unique(omnibus$model), c("CORE", "BOTH")),
  setequal(unique(omnibus$transition), paths), all(omnibus$df == 3L),
  all(is.finite(omnibus$chi_square)), all(omnibus$chi_square >= 0),
  all(is.finite(omnibus$p)), all(omnibus$p > 0 & omnibus$p <= 1),
  all(is.finite(omnibus$p_BH_4)), all(omnibus$p_BH_4 > 0 & omnibus$p_BH_4 <= 1),
  all(omnibus$reference == "Asymptotic chi-square Wald on Rubin total covariance"))
# Read-only consistency checks; no source values are replaced or recomputed for plotting.
stopifnot(max(abs(pchisq(omnibus$chi_square, df = 3, lower.tail = FALSE) - omnibus$p)) < 1e-12)
for (spec in c("CORE", "BOTH")) {
  z <- omnibus[model == spec]
  stopifnot(nrow(z) == 4L, max(abs(p.adjust(z$p, method = "BH") - z$p_BH_4)) < 1e-12)
}
omnibus[, model := factor(model, levels = c("CORE", "BOTH"))]
omnibus[, plot_y := 10 - 2 * match(transition, paths) - as.integer(model == "BOTH")]
fmt_exact_p <- function(value) vapply(value, function(v) {
  if (v < .001) formatC(v, digits = 5, format = "e")
  else format(signif(v, 6), scientific = FALSE, trim = TRUE)
}, "")
omnibus[, statistic_label := sprintf("%.3f", chi_square)]
omnibus[, raw_p_label := fmt_exact_p(p)]
omnibus[, bh_p_label := fmt_exact_p(p_BH_4)]
setorder(omnibus, edge_order, model)

fig6 <- ggplot(omnibus, aes(x = chi_square, y = plot_y, colour = model, shape = model)) +
  geom_point(size = 2.3) +
  geom_text(aes(x = 31.0, label = statistic_label), colour = "#27313B",
    hjust = 0.5, size = 2.45, family = font, show.legend = FALSE) +
  geom_text(aes(x = 43.0, label = raw_p_label), colour = "#27313B",
    hjust = 0.5, size = 2.45, family = font, show.legend = FALSE) +
  geom_text(aes(x = 56.0, label = bh_p_label), colour = "#27313B",
    hjust = 0.5, size = 2.45, family = font, show.legend = FALSE) +
  annotate("text", x = 31.0, y = 9.05, label = "Wald statistic",
    hjust = .5, size = 2.45, fontface = "bold", family = font) +
  annotate("text", x = 43.0, y = 9.05, label = "Raw P",
    hjust = .5, size = 2.45, fontface = "bold", family = font) +
  annotate("text", x = 56.0, y = 9.05, label = "BH P",
    hjust = .5, size = 2.45, fontface = "bold", family = font) +
  scale_colour_manual(values = cols) + scale_shape_manual(values = shapes) +
  scale_x_continuous(limits = c(0, 63), breaks = seq(0, 25, 5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(.35, 9.65), breaks = c(7.5, 5.5, 3.5, 1.5),
    labels = path_names, expand = c(0, 0)) +
  labs(title = "Overall evidence across four dynamic sleep patterns",
    subtitle = "Joint asymptotic Wald tests; df = 3; all four target pathways and both specifications",
    x = "Wald chi-square (df = 3)", y = NULL,
    caption = paste(
      "CHARLS W2 landmark cohort: 9,397 participants; 20 imputations per specification. Rubin total covariance.",
      "BH correction: four pathway tests within each specification. P values shown to six significant digits.",
      "No test-statistic confidence intervals apply; effect-size CIs appear in Figures 1, 2 and 5. HTN: hypertension.",
      "Omnibus evidence does not establish which patterns differ, effect direction, or a CORE-BOTH difference.", sep = "\n")) +
  theme_pub() + theme(legend.position = "top", legend.justification = "left",
    axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    axis.title.x = element_text(hjust = .03, family = "Arial", size = 7),
    plot.title = element_text(family = "Arial", size = 9, face = "bold"),
    plot.subtitle = element_text(family = "Arial", size = 7),
    plot.caption = element_text(family = "Arial", size = 6.3, hjust = 0,
      lineheight = 1.2, margin = margin(t = 9)),
    plot.margin = margin(8, 9, 8, 8))

export_single(fig6, "figure6_pattern_omnibus", 130)
fwrite(omnibus, file.path(dest, "figure6_source_data.csv"), bom = TRUE)
assert_scientific_unchanged()
stopifnot(identical(sha(omnibus_input), omnibus_hash))
write_json(list(status = "RENDERED_REQUIRES_FINAL_QA", input_sha256 = omnibus_hash,
  input_rows = 8L, displayed_rows = 8L, excluded_rows = 0L,
  source = "R/17m_dynamic_omnibus_EN.R", backend = "R", font = "Arial",
  primary_statistic = "Precomputed asymptotic chi-square Wald statistic on Rubin total covariance",
  degrees_of_freedom = 3L, p_adjustment = "BH over four pathways within each specification",
  p_display = "Six significant digits; full-precision source data retained",
  uncertainty = "Not applicable to displayed test statistics; effect estimates with CIs appear elsewhere",
  altered_scientific_values = FALSE, raw_p_and_BH_crosschecks = "PASS",
  scientific_frozen_hashes_unchanged = TRUE, alignment = "NOT_APPLICABLE_SINGLE_PANEL"),
  file.path(dest, "figure6_render_QC.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Figure 6 rendered; all eight precomputed omnibus tests preserved.\n")
