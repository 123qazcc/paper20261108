#!/usr/bin/env Rscript
# Presentation-only reuse: remove obsolete report-table references in captions.
options(charls.visual.helpers_only = TRUE)
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
source("R/17i_charls_dynamic_figures_EN.R", encoding = "UTF-8")
font <- "Arial"
# Re-declare the unchanged inherited caption theme for standalone source QA.
publication_theme <- theme(plot.title = element_text(family = "Arial", size = 9, face = "bold"),
  plot.subtitle = element_text(family = "Arial", size = 7),
  plot.caption = element_text(family = "Arial", size = 6.3, hjust = 0, lineheight = 1.15, margin = margin(t = 8)))
out <- "outputs/Three_cohort_sleep_CVD_20261001"
figdir <- file.path(out, "figures")
srcdir <- file.path(out, "source_data")
qadir <- "outputs/three_cohort_integrated_20261001/figure_qa"
for (d in c(figdir, srcdir, qadir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

fig4 <- fig2 + plot_annotation(theme = publication_theme,
  caption = paste("CORE: blue circles; BOTH: orange triangles. Bars: pointwise 95% CIs (20 imputations).",
    "BH correction: eight primary and eight contextual contrasts, separately per specification.",
    "HTN: hypertension. No target disease: no HTN, heart disease or stroke. Patterns use observed W1-W2 sleep.", sep = "\n"))
fig5 <- fig3 + plot_annotation(theme = publication_theme,
  caption = paste("Points: ratios of association HRs, not disease incidence rates. Bars: pointwise 95% CIs (20 imputations).",
    "BH correction: four primary and four contextual direct contrasts per specification. HTN: hypertension.", sep = "\n"))

export_integrated_dynamic <- function(plot, panels, stem, height_mm) {
  path <- file.path(figdir, stem); qa <- file.path(qadir, stem)
  w <- 183 / 25.4; h <- height_mm / 25.4
  require_patchwork_panel_alignment(plot, paste0(qa, ".alignment-layout.json"),
    paste0(qa, ".alignment.json"), width_in = w, height_in = h,
    overlay_svg = paste0(qa, ".alignment.svg"), tolerance_pt = 1.5, gutter_tolerance_pt = 1.5,
    strict = TRUE, audit_script = file.path(skill, "scripts/audit_panel_alignment.py"), python = python)
  cairo_pdf(paste0(path, ".pdf"), width = w, height = h, family = font); print(plot); dev.off()
  svglite(paste0(path, ".svg"), width = w, height = h); print(plot); dev.off()
  render_raster(plot, paste0(path, ".tiff"), w, h, 600, tiff = TRUE)
  render_raster(plot, paste0(path, ".png"), w, h, 300)
  for (j in seq_along(panels)) {
    render_raster(panels[[j]], paste0(qa, ".panel_", letters[j], ".png"), w / 2, h / 2, 300)
  }
  write_json(list(backend = "R", width_mm = 183, height_mm = height_mm,
    tiff_dpi = 600, png_dpi = 300, font = font,
    change = "Removed obsolete table reference from caption only; scientific data and panel objects unchanged"),
    paste0(qa, ".export.json"), pretty = TRUE, auto_unbox = TRUE)
}
export_integrated_dynamic(fig4, target_panels, "Figure4_dynamic_patterns", 175)
export_integrated_dynamic(fig5, direct_panels, "Figure5_dynamic_direct", 145)
fwrite(target, file.path(srcdir, "Figure4_dynamic_patterns.csv"), bom = TRUE)
fwrite(direct, file.path(srcdir, "Figure5_dynamic_direct.csv"), bom = TRUE)
stopifnot(identical(vapply(input_files, sha, ""), input_hashes),
  identical(unname(vapply(names(frozen), sha, "")), unname(frozen)))
write_json(list(input_hashes = as.list(input_hashes), frozen_source_hashes = as.list(frozen),
  changed_scientific_values = FALSE, models_refitted = FALSE,
  figure4_rows = nrow(target), figure5_rows = nrow(direct)),
  file.path(qadir, "dynamic_source_provenance.json"), auto_unbox = TRUE, pretty = TRUE)
