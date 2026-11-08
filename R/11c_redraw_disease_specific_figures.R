#!/usr/bin/env Rscript
# Plot-only delivery repair. No model/input changes or fitting are performed.
# Contract: DS1 displays all eight target estimates; DS2 directly tests disease
# specificity. Two single-panel forests, R-only, editable vectors and 600-dpi
# rasters, all intervals retained, nominal versus BH inference distinguished.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(svglite); library(ragg)
})
root <- file.path("outputs", "hrs_elsa_multistate_nextsteps_20260909")
figdir <- file.path(root, "figures", "disease_specific")
resdir <- file.path(root, "results", "disease_specific")
qadir <- file.path(root, "qc", "disease_specific", "final_figure_review_20260911")
dir.create(qadir, recursive = TRUE, showWarnings = FALSE)
archive <- file.path(qadir, "before_redraw")
dir.create(archive, showWarnings = FALSE)
old <- list.files(figdir, pattern = "^figure_ds[12].*\\.(pdf|svg|png|tiff)$", full.names = TRUE)
for (f in old) if (!file.exists(file.path(archive, basename(f)))) {
  stopifnot(file.copy(f, file.path(archive, basename(f))))
}
source("R/helpers/rendered_panel_alignment.R")
local_qa <- normalizePath(".runtime/figure_qa", winslash = "/", mustWork = TRUE)
Sys.setenv(PYTHONPATH = paste(c(local_qa, Sys.getenv("PYTHONPATH")), collapse = .Platform$path.sep))
target <- fread(file.path(resdir, "disease_specific_target_edge_hr_mi20.csv"))
contrast <- fread(file.path(resdir, "disease_specific_heart_vs_stroke_contrasts_mi20.csv"))
stopifnot(nrow(target) == 8L, nrow(contrast) == 4L,
          all(target$m_fitted == 20L), all(contrast$m_fitted == 20L))
target[, row_order := match(cohort, c("HRS", "ELSA")) * 10L + target_order]
target[, display_label := paste0(cohort, "  |  ", c(
  "No disease -> heart disease only", "No disease -> stroke only",
  "HTN only -> HTN + heart disease", "HTN only -> HTN + stroke"
)[target_order])]
target[, display_label := factor(display_label, levels = rev(display_label[order(row_order)]))]
contrast[, row_order := match(cohort, c("HRS", "ELSA")) * 10L +
           match(origin_state_label, c("No disease", "Hypertension only"))]
contrast[, display_label := paste0(cohort, "  |  ", origin_state_label)]
contrast[, display_label := factor(display_label, levels = rev(display_label[order(row_order)]))]
stopifnot(all(is.finite(target$conf_low)), all(is.finite(target$conf_high)),
          min(target$conf_low) > 0.8, max(target$conf_high) < 1.6,
          all(is.finite(contrast$conf_low)), all(is.finite(contrast$conf_high)),
          min(contrast$conf_low) > 0.7, max(contrast$conf_high) < 1.45)
theme_final <- theme_classic(base_size = 7, base_family = "Helvetica") + theme(
  axis.text = element_text(size = 7, colour = "#26323D"),
  axis.title = element_text(size = 7.5),
  plot.title = element_text(size = 9, face = "bold", margin = margin(b = 5)),
  plot.subtitle = element_text(size = 7, margin = margin(b = 7)),
  plot.title.position = "plot", plot.caption.position = "plot",
  plot.caption = element_text(size = 6.5, lineheight = 1.35, hjust = 0,
                              margin = margin(t = 8), colour = "#4B5563"),
  legend.position = "top", legend.title = element_blank(),
  legend.text = element_text(size = 7),
  panel.grid.major.x = element_line(colour = "#E3E7EA", linewidth = 0.25),
  plot.margin = margin(8, 10, 8, 8)
)
p1 <- ggplot(target, aes(hazard_ratio, display_label, colour = target_outcome, shape = target_outcome)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#68737C", linewidth = 0.4) +
  geom_segment(aes(x = conf_low, xend = conf_high, yend = display_label), linewidth = 0.65) +
  geom_point(size = 2.1) +
  scale_x_log10(limits = c(0.8, 1.6), breaks = c(0.8, 1, 1.2, 1.4, 1.6),
                labels = c("0.8", "1.0", "1.2", "1.4", "1.6")) +
  scale_colour_manual(values = c("Heart disease" = "#2166AC", "Stroke" = "#B2182B")) +
  scale_shape_manual(values = c("Heart disease" = 15, "Stroke" = 16)) +
  labs(title = "Disease-specific associations with insomnia burden",
       subtitle = "Adjusted M4-MI20 panel Markov models; per one additional insomnia symptom (0-3)",
       x = "Hazard ratio (95% CI; logarithmic scale)", y = NULL,
       caption = paste(
         "Points and lines: Rubin-pooled HRs and finite-imputation t 95% CIs; target-edge P values are unadjusted.",
         "Independent unit: participant. HRS n=29,029; ELSA n=11,097; 20 imputations per cohort. HTN: hypertension.",
         "Post-hoc secondary extension using frozen primary MI; HRS stroke includes TIA; ELSA W4-W10 excludes death.",
         "Broad heart disease does not distinguish myocardial infarction, angina or heart failure. See accompanying source data.",
         sep = "\n")) + theme_final
p2 <- ggplot(contrast, aes(hr_ratio_heart_vs_stroke, display_label,
                           colour = origin_state_label, shape = origin_state_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#68737C", linewidth = 0.4) +
  geom_segment(aes(x = conf_low, xend = conf_high, yend = display_label), linewidth = 0.65) +
  geom_point(size = 2.1) +
  scale_x_log10(limits = c(0.7, 1.45), breaks = c(0.7, 0.8, 1, 1.2, 1.4),
                labels = c("0.7", "0.8", "1.0", "1.2", "1.4")) +
  scale_colour_manual(values = c("No disease" = "#0072B2", "Hypertension only" = "#D55E00"),
                      breaks = c("No disease", "Hypertension only")) +
  scale_shape_manual(values = c("No disease" = 16, "Hypertension only" = 15),
                     breaks = c("No disease", "Hypertension only")) +
  labs(title = "Direct comparison of heart-disease and stroke associations",
       subtitle = "Within the same origin state: HR ratio = HR for heart disease / HR for stroke",
       x = "Heart-disease / stroke HR ratio (95% CI; logarithmic scale)", y = NULL,
       caption = paste(
         "Points and lines: Rubin-pooled HR ratios and finite-imputation t 95% CIs; full within-model covariance retained.",
         "Independent unit: participant. HRS n=29,029; ELSA n=11,097; 20 imputations per cohort.",
         "BH correction covers two contrasts per cohort; all BH P > 0.32. CIs are pointwise, not multiplicity-adjusted.",
         "Ratios above 1 indicate a larger heart-disease coefficient. A non-significant difference does not imply equivalence.",
         sep = "\n")) + theme_final
export_plot <- function(p, stem, height) {
  dest <- file.path(figdir, stem)
  require_rendered_panel_alignment(p, file.path(qadir, stem), 7.2, height)
  pdf(paste0(dest, ".pdf"), width = 7.2, height = height,
      family = "Helvetica", useDingbats = FALSE)
  print(p); dev.off()
  svglite(paste0(dest, ".svg"), width = 7.2, height = height,
          system_fonts = list(sans = "Arial"))
  print(p + theme(text = element_text(family = "Arial"))); dev.off()
  agg_tiff(paste0(dest, ".tiff"), width = 7.2, height = height, units = "in", res = 600, compression = "lzw")
  print(p + theme(text = element_text(family = "Arial"))); dev.off()
  agg_png(paste0(dest, ".png"), width = 7.2, height = height, units = "in", res = 600)
  print(p + theme(text = element_text(family = "Arial"))); dev.off()
}
export_plot(p1, "figure_ds1_target_transition_forest", 5.0)
export_plot(p2, "figure_ds2_heart_vs_stroke_contrasts", 4.0)
fwrite(target, file.path(figdir, "source_data", "figure_ds1_final_display_source.csv"))
fwrite(contrast, file.path(figdir, "source_data", "figure_ds2_final_display_source.csv"))
writeLines(c("R-only plot repair; all 8 target estimates and all 4 contrasts retained.",
             "Only layout, labels, redundant encoding and tick coverage changed; model tables are unchanged.",
             "Each figure has one plot area: multi-panel alignment is not applicable.",
             "Original exports retained in before_redraw; rerun PDF text/collision audits before delivery."),
           file.path(qadir, "redraw_scope.txt"))
