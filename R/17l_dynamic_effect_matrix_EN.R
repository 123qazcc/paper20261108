#!/usr/bin/env Rscript
# Complete observed estimates only; no refit, repooling, or outcome selection.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("R/17j_dynamic_extra_figure_helpers.R", encoding = "UTF-8")
stopifnot(nrow(target) == 32L, all(target$HR > 0), all(target$HR >= .5 & target$HR <= 2))
matrix_names <- c("Onset vs\npersistent low\n[Primary]", "Persistent high\nvs remission\n[Primary]",
                  "Remission vs\npersistent low\n[Contextual]", "Persistent high\nvs persistent low\n[Contextual]")
matrix_panels <- lapply(c("CORE", "BOTH"), function(spec) {
  z <- copy(target[target$model == spec])
  z[, x := match(comparison, comparisons)]
  z[, y := 5 - match(transition, paths)]
  z[, cell_label := sprintf("%.2f%s\n[%.2f–%.2f]", HR, ifelse(p_BH_family < .05, " *", ""), lower, upper)]
  ggplot(z, aes(x, y, fill = log(HR))) +
    geom_tile(width = .97, height = .94, colour = "white", linewidth = .3) +
    geom_text(aes(label = cell_label), family = "Arial", size = 2.35, lineheight = 1.25, colour = "#1E2931") +
    scale_fill_gradient2(low = "#AAC2D3", mid = "#FAFAFA", high = "#DFA879", midpoint = 0,
      limits = log(c(.5, 2)), breaks = log(c(.5, .75, 1, 1.5, 2)), labels = c("0.5", "0.75", "1", "1.5", "2"),
      name = "HR (log colour scale)") +
    scale_x_continuous(breaks = 1:4, labels = matrix_names, limits = c(.5, 4.5), expand = c(0, 0), position = "top") +
    scale_y_continuous(breaks = 4:1, labels = path_names, limits = c(.5, 4.5), expand = c(0, 0)) +
    labs(title = paste0(spec, " specification"), x = NULL, y = NULL) +
    theme_minimal(base_size = 7, base_family = "Arial") +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 6.5, colour = "#27313B"),
      axis.text.x = element_text(lineheight = 1.1, margin = margin(b = 7)), axis.text.y = element_text(hjust = 1),
      plot.title = element_text(size = 8, face = "bold", margin = margin(b = 9)),
      legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5),
      legend.key.width = unit(18, "mm"), legend.key.height = unit(3, "mm"),
      plot.margin = margin(6, 6, 6, 5))
})
matrix_figure <- wrap_plots(matrix_panels, ncol = 1, guides = "collect") +
  plot_annotation(tag_levels = "a", title = "Sleep-pattern associations across cardiovascular pathways",
    subtitle = "Complete effect matrix | CORE and additional CESD-9 / sleep-duration adjustment | n = 9,397",
    caption = paste("Cells: HR [pointwise 95% CI], Rubin pooled across 20 imputations. * BH-adjusted P < 0.05.",
      "BH families: eight primary and eight contextual contrasts, separately per specification. Colour encodes HR, not significance.",
      "HTN: hypertension. No target disease: no HTN, heart disease or stroke. All 32 estimates are shown.", sep = "\n"),
    theme = publication_theme) &
  theme(legend.position = "bottom", plot.tag = element_text(family = "Arial", face = "bold", size = 9))
export_plot(matrix_figure, matrix_panels, "figure5_effect_matrix", 175)
fwrite(target, file.path(dest, "figure5_source_data.csv"), bom = TRUE)
assert_scientific_unchanged()
cat("Figure 5 rendered; complete 32-row matrix; final QA required.\n")
