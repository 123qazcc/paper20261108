#!/usr/bin/env Rscript
# Presentation only: no fitting or pooling.
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(jsonlite)
  library(digest); library(svglite); library(ragg); library(grid)
})
setDTthreads(1L)
root <- "outputs/charls_dynamic_sleep_20260920"
dest <- file.path(root, "stage3_delivery/visual_en_v3")
figdir <- file.path(dest, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
review <- fromJSON(file.path(root, "model_qc/dynamic_independent_pool_review.json"))
stopifnot(review$status == "COMPUTATIONAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_REVIEW",
          review$checks_passed == 230L, review$checks_total == 230L)
frozen <- c("R/17c_fit_charls_dynamic_models.R" = "41f8a6924d7ac4e8c07fc0e87bdd7b82ac3962bc1d2e61a1aa953c8847a6279e",
            "R/17d_independently_verify_charls_dynamic_pool.R" = "07a2a6673b3e3250582d91e6d55ec8ea7c68b7ae3b9851856f0425b8b032c2c6",
            "R/17e_build_charls_dynamic_delivery.R" = "ba3b9eafafc40ea32b3a28da211441b6b039c4184d53b92da655f1920344b1c3",
            "analysis/charls_stage3_dynamic_sleep_control.json" = "01c9bd7d4f33bb683b1a03cf72ce54d4c6a872e26f5333895db90d901ee158b4")
sha <- function(p) digest(file = p, algo = "sha256")
stopifnot(identical(unname(vapply(names(frozen), sha, "")), unname(frozen)))
input_files <- file.path(root, "results", c("four_target_dynamic_contrasts_MI20.csv",
  "heart_stroke_dynamic_direct_MI20.csv", "joint_tests_MI20.csv"))
input_hashes <- vapply(input_files, sha, "")
target <- fread(input_files[1]); direct <- fread(input_files[2]); joint <- fread(input_files[3])
stopifnot(nrow(target) == 32L, nrow(direct) == 16L,
  uniqueN(target, by = c("model", "comparison", "transition")) == 32L,
  uniqueN(direct, by = c("model", "comparison", "origin")) == 16L,
  all(is.finite(target$HR)), all(target$lower > 0), all(target$upper >= target$HR),
  all(is.finite(direct$HR)), all(direct$lower > 0), all(direct$upper >= direct$HR))
font <- "Arial"
cols <- c(CORE = "#2C6082", BOTH = "#BF7040")
shapes <- c(CORE = 16, BOTH = 17)
comparisons <- c("onset_vs_persistent_low", "persistent_high_vs_remission",
                 "remission_vs_persistent_low", "persistent_high_vs_persistent_low")
comparison_names <- c("Onset vs persistent low", "Persistent high vs remission",
                      "Remission vs persistent low", "Persistent high vs persistent low")
paths <- c("0-2", "0-4", "1-3", "1-5")
path_names <- c("No target disease\n→ Heart disease only", "No target disease\n→ Stroke only",
                "HTN only\n→ HTN + heart disease", "HTN only\n→ HTN + stroke")
short_paths <- c("No target disease\n→ Heart disease only", "No target disease\n→ Stroke only",
                 "HTN only\n→ HTN + heart disease", "HTN only\n→ HTN + stroke")
theme_pub <- function() theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.text = element_text(size = 6.5, colour = "#27313B"),
    axis.text.y = element_text(family = font),
    axis.title = element_text(size = 7, family = font), axis.line = element_line(linewidth = .3),
    axis.ticks = element_line(linewidth = .3),
    panel.grid.major.x = element_line(colour = "#E8ECEF", linewidth = .25),
    plot.title = element_text(size = 8, family = font, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 6.5, family = font, colour = "#596773", margin = margin(b = 6)),
    plot.tag = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 7), legend.title = element_blank(),
    plot.margin = margin(7, 7, 7, 5))
fmt_p <- function(p) ifelse(p < .001, formatC(p, digits = 2, format = "e"),
                           formatC(p, digits = 4, format = "f"))
skill <- Sys.getenv("NATURE_FIGURE_SKILL_ROOT", "private/render_tools/nature-figure")
python <- Sys.getenv("NATURE_FIGURE_PYTHON", "python")
stopifnot(file.exists(file.path(skill, "scripts/panel_alignment.R")), file.exists(python))
source(file.path(skill, "scripts/panel_alignment.R"))

# Local measurement adapter: same Cairo font/device as final PDF, with resolved
# grid viewport sizes. Strict auditor and 1.5pt tolerances remain unchanged.
write_patchwork_panel_layout <- function(plot, manifest_path, width_in, height_in,
    panel_ids = NULL, row_groups = NULL, column_groups = NULL, exemptions = list()) {
  probe <- tempfile(fileext = ".pdf")
  cairo_pdf(probe, width = width_in, height = height_in, family = font)
  on.exit({dev.off(); unlink(probe)}, add = TRUE)
  g <- patchwork::patchworkGrob(plot)
  grid.newpage(); grid.draw(g); grid.force()
  pr <- g$layout[grepl("^panel(?:-[0-9]+)?$", g$layout$name, perl = TRUE), ]
  pr <- pr[order(pr$t, pr$l), ]
  if (is.null(panel_ids)) panel_ids <- letters[seq_len(nrow(pr))]
  stopifnot(length(panel_ids) == nrow(pr), nrow(pr) >= 2L)
  bounds <- function(row) {
    vp <- sprintf("%s.%d-%d-%d-%d", row$name, row$t, row$l, row$b, row$r)
    seekViewport(vp)
    lo <- deviceLoc(unit(0, "npc"), unit(0, "npc"), valueOnly = TRUE)
    hi <- deviceLoc(unit(1, "npc"), unit(1, "npc"), valueOnly = TRUE)
    upViewport(0)
    c(lo$x, lo$y, hi$x, hi$y) * 72
  }
  panels <- lapply(seq_len(nrow(pr)), function(i) {
    row <- pr[i, ]; b <- bounds(row)
    stopifnot(b[3] > b[1], b[4] > b[2])
    tag <- g$layout[g$layout$name == sub("^panel", "tag", row$name), ]
    ans <- list(id = panel_ids[i], bbox_pt = unname(b), grid_id = "patchwork-grid-1",
      row_start = row$t - 1L, row_stop = row$b, col_start = row$l - 1L, col_stop = row$r)
    if (nrow(tag) == 1L) { tb <- bounds(tag); ans$panel_label_anchor_pt <- unname(c(tb[1], tb[4])) }
    ans
  })
  m <- list(schema_version = 1L, backend = "r-patchwork-cairo-resolved-viewports",
    figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
    panels = panels, exemptions = exemptions)
  if (!is.null(row_groups)) m$row_groups <- lapply(row_groups, as.character)
  if (!is.null(column_groups)) m$column_groups <- lapply(column_groups, as.character)
  write_json(m, manifest_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  invisible(m)
}

primary_panels <- lapply(1:2, function(i) {
  z <- copy(target[model == "CORE" & comparison == comparisons[i]])
  z[, y := 5 - match(transition, paths)]
  z[, numeric_label := sprintf("%.2f (%.2f–%.2f)", HR, lower, upper)]
  z[, bh_label := fmt_p(p_BH_family)]
  ggplot(z, aes(x = HR, y = y)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = .35, colour = "#657583") +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .14,
      linewidth = .6, colour = cols[["CORE"]]) +
    geom_point(size = 2.2, colour = cols[["CORE"]], shape = 16) +
    geom_text(aes(x = 3.8, label = numeric_label), hjust = 0, size = 2.45, family = font) +
    geom_text(aes(x = 13.5, label = bh_label), hjust = 0, size = 2.45, family = font) +
    annotate("text", x = 3.8, y = 4.7, label = "HR (95% CI)", hjust = 0,
      size = 2.45, fontface = "bold", family = font) +
    annotate("text", x = 13.5, y = 4.7, label = "BH P", hjust = 0,
      size = 2.45, fontface = "bold", family = font) +
    scale_x_log10(limits = c(.48, 24), breaks = c(.5, 1, 2, 3), expand = c(0, 0)) +
    scale_y_continuous(limits = c(.5, 5), breaks = 4:1, labels = path_names, expand = c(0, 0)) +
    labs(title = comparison_names[i], subtitle = if (i == 1) "Primary contrast | Low frequency at W1" else "Primary contrast | High frequency at W1",
         x = "Transition hazard ratio (log scale)", y = NULL) + theme_pub()
})
fig1 <- wrap_plots(primary_panels, ncol = 1) +
  plot_annotation(tag_levels = "a", title = "Dynamic restless-sleep patterns and cardiovascular transitions",
    subtitle = "CHARLS | Primary CORE model | W2 landmark cohort: 9,397 participants") &
  theme(plot.title = element_text(family = font), plot.subtitle = element_text(family = font),
        plot.tag = element_text(family = font, size = 9, face = "bold"))

target_panels <- lapply(seq_along(comparisons), function(i) {
  z <- copy(target[comparison == comparisons[i]])
  z[, y := 5 - match(transition, paths) + ifelse(model == "CORE", .13, -.13)]
  ggplot(z, aes(HR, y, colour = model, shape = model)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = .3, colour = "#6B7780") +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .11, linewidth = .48) +
    geom_point(size = 1.9) +
    scale_colour_manual(values = cols) + scale_shape_manual(values = shapes) +
    scale_x_log10(limits = c(.45, 3.6), breaks = c(.5, 1, 2, 3)) +
    scale_y_continuous(limits = c(.5, 4.5), breaks = 4:1, labels = short_paths, expand = c(0, 0)) +
    labs(title = comparison_names[i], subtitle = if (i <= 2) "Primary contrast" else "Contextual contrast",
      x = "HR (95% CI)", y = NULL) + theme_pub()
})
fig2 <- wrap_plots(target_panels, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a", title = "Complete pattern contrasts and adjustment robustness",
    subtitle = "BOTH additionally adjusts for W1 CESD-9 and night-time sleep duration | n = 9,397") &
  theme(legend.position = "top", plot.tag = element_text(family = font, size = 9, face = "bold"),
        plot.title = element_text(family = font), plot.subtitle = element_text(family = font))

direct_panels <- lapply(seq_along(comparisons), function(i) {
  z <- copy(direct[comparison == comparisons[i]])
  z[, y := 3 - origin_order + ifelse(model == "CORE", .11, -.11)]
  ggplot(z, aes(HR, y, colour = model, shape = model)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = .3, colour = "#6B7780") +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .10, linewidth = .48) +
    geom_point(size = 1.9) +
    scale_colour_manual(values = cols) + scale_shape_manual(values = shapes) +
    scale_x_log10(limits = c(.33, 3), breaks = c(.5, 1, 2)) +
    scale_y_continuous(limits = c(.5, 2.5), breaks = 2:1,
      labels = c("No target disease", "HTN only"), expand = c(0, 0)) +
    labs(title = comparison_names[i], subtitle = if (i <= 2) "Primary direct contrast" else "Contextual direct contrast",
      x = "Heart-disease HR / stroke HR (95% CI)", y = "Origin state") + theme_pub()
})
fig3 <- wrap_plots(direct_panels, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a", title = "Direct comparison of heart-disease and stroke associations",
    subtitle = "Primary direct contrasts do not meet BH-adjusted P < 0.05 | n = 9,397") &
  theme(legend.position = "top", plot.tag = element_text(family = font, size = 9, face = "bold"),
        plot.title = element_text(family = font), plot.subtitle = element_text(family = font))

render_raster <- function(plot, file, w, h, dpi, tiff = FALSE) {
  # ragg's Windows path handling requires an ASCII temporary output here.
  # R then copies the completed bitmap unchanged to the Unicode project path.
  tmp <- tempfile("charls-visual-", fileext = if (tiff) ".tiff" else ".png")
  if (tiff) agg_tiff(tmp, width = w, height = h, units = "in", res = dpi, compression = "lzw")
  else agg_png(tmp, width = w, height = h, units = "in", res = dpi)
  print(plot); dev.off()
  stopifnot(file.exists(tmp), file.info(tmp)$size > 0, file.copy(tmp, file, overwrite = TRUE))
  unlink(tmp)
}
export_plot <- function(plot, panels, stem, height_mm) {
  path <- file.path(figdir, stem); w <- 183 / 25.4; h <- height_mm / 25.4
  require_patchwork_panel_alignment(plot, paste0(path, ".alignment-layout.json"),
    paste0(path, ".alignment.json"), width_in = w, height_in = h,
    overlay_svg = paste0(path, ".alignment.svg"), tolerance_pt = 1.5, gutter_tolerance_pt = 1.5,
    strict = TRUE, audit_script = file.path(skill, "scripts/audit_panel_alignment.py"), python = python)
  cairo_pdf(paste0(path, ".pdf"), width = w, height = h, family = font); print(plot); dev.off()
  svglite(paste0(path, ".svg"), width = w, height = h); print(plot); dev.off()
  render_raster(plot, paste0(path, ".png"), w, h, 600)
  render_raster(plot, paste0(path, ".tiff"), w, h, 600, tiff = TRUE)
  render_raster(plot, paste0(path, ".preview.png"), w, h, 300)
  for (j in seq_along(panels)) {
    render_raster(panels[[j]], paste0(path, ".panel_", letters[j], ".png"),
      if (length(panels) == 2L) w else w / 2, h / 2, 300)
  }
  files <- paste0(path, c(".pdf", ".svg", ".png", ".tiff"))
  write_json(list(backend = "R", width_mm = 183, height_mm = height_mm, dpi = 600,
    font = font, files = basename(files), sha256 = as.list(vapply(files, sha, ""))),
    paste0(path, ".export.json"), pretty = TRUE, auto_unbox = TRUE)
}
publication_theme <- theme(plot.title = element_text(family = "Arial", size = 9, face = "bold"),
  plot.subtitle = element_text(family = "Arial", size = 7),
  plot.caption = element_text(family = "Arial", size = 6.3, hjust = 0, lineheight = 1.15, margin = margin(t = 8)))
fig1 <- fig1 + plot_annotation(theme = publication_theme,
  caption = paste("Points and bars: Rubin-pooled HRs and pointwise 95% CIs (20 imputations). BH P: eight primary contrasts.",
    "HTN: hypertension. No target disease: no HTN, heart disease or stroke. High-frequency restless sleep: >=3 days/week.", sep = "\n"))
fig2 <- fig2 + plot_annotation(theme = publication_theme,
  caption = paste("CORE: blue circles; BOTH: orange triangles. Bars: pointwise 95% CIs (20 imputations).",
    "BH correction: eight primary and eight contextual contrasts, separately per specification (Tables 3-4).",
    "HTN: hypertension. No target disease: no HTN, heart disease or stroke. Patterns use observed W1-W2 sleep.", sep = "\n"))
fig3 <- fig3 + plot_annotation(theme = publication_theme,
  caption = paste("Points: ratios of association HRs, not disease incidence rates. Bars: pointwise 95% CIs (20 imputations).",
    "BH correction: four primary and four contextual direct contrasts per specification (Table 6). HTN: hypertension.", sep = "\n"))
if (!isTRUE(getOption("charls.visual.helpers_only", FALSE))) {
export_plot(fig1, primary_panels, "figure1_primary_dynamic", 150)
export_plot(fig2, target_panels, "figure2_full_patterns", 175)
export_plot(fig3, direct_panels, "figure3_direct_heart_stroke", 145)
fwrite(target[model == "CORE" & role == "primary"], file.path(dest, "figure1_source_data.csv"), bom = TRUE)
fwrite(target, file.path(dest, "figure2_source_data.csv"), bom = TRUE)
fwrite(direct, file.path(dest, "figure3_source_data.csv"), bom = TRUE)
fwrite(joint, file.path(dest, "joint_tests_source_data.csv"), bom = TRUE)
stopifnot(identical(vapply(input_files, sha, ""), input_hashes),
          identical(unname(vapply(names(frozen), sha, "")), unname(frozen)))
write_json(list(status = "FIGURES_RENDERED_AWAITING_FINAL_QA", accepted_models = 40L,
  independent_checks_passed = 230L, input_hashes = as.list(input_hashes),
  frozen_source_hashes = as.list(frozen), changed_scientific_values = FALSE,
  failed_dispatcher_guard_cleared = FALSE, automatic_next_stage = FALSE,
  source = "R/17i_charls_dynamic_figures_EN.R", generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  file.path(dest, "render_status.json"), pretty = TRUE, auto_unbox = TRUE)
cat("PRESENTATION_REPAIR_RENDERED: scientific inputs unchanged; final QA still required.\n")
}
