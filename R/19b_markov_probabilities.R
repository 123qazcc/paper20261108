#!/usr/bin/env Rscript
# Existing aggregate results only: no model fitting, pooling, interpolation or normalization.
# Figure contract: quantitative grids. Figure 3 decomposes fixed-profile endpoint
# probabilities by state; Figure 4 quantifies their paired absolute consequences.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid); library(jsonlite)
  library(svglite); library(ragg); library(digest)
})
out <- "outputs/Three_cohort_Markov_results_20261001"
figdir <- file.path(out, "figures"); srcdir <- file.path(out, "source_data")
qadir <- "outputs/three_cohort_markov_revision_20261001/figure_qa"
for (d in c(figdir, srcdir, qadir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
skill <- Sys.getenv("NATURE_FIGURE_SKILL_ROOT", "private/render_tools/nature-figure")
python <- Sys.getenv("NATURE_FIGURE_PYTHON", "python")
source(file.path(skill, "scripts/panel_alignment.R"))
font <- "Arial"
paths <- c(marginal = "outputs/hrs_elsa_multistate_nextsteps_20260909/results/absolute_risk/absolute_quantities_mi20_95ci.csv",
  paired = "outputs/hrs_elsa_multistate_nextsteps_20260909/results/absolute_risk/fixed_score_2_vs_0_paired_differences.csv")
hashes <- vapply(paths, function(p) digest(file = p, algo = "sha256"), "")
marginal_raw <- read.csv(paths[["marginal"]], check.names = FALSE)
paired_raw <- read.csv(paths[["paired"]], check.names = FALSE)
stopifnot(nrow(marginal_raw) == 264L, nrow(paired_raw) == 66L)
score_colours <- c("Exactly 0 symptoms" = "#275F85", "Exactly 2 symptoms" = "#B27036")
score_shapes <- c("Exactly 0 symptoms" = 16, "Exactly 2 symptoms" = 17)
cohort_colours <- c(HRS = "#275F85", ELSA = "#B27036")
cohort_shapes <- c(HRS = 16, ELSA = 17)
theme_pub <- function() theme_classic(base_size = 7, base_family = "Arial") + theme(
  text = element_text(family = font, colour = "#26323A"),
  axis.title = element_text(size = 7), axis.text = element_text(size = 6.5, colour = "#26323A"),
  axis.line = element_line(linewidth = .28, colour = "#6F7D86"),
  axis.ticks = element_line(linewidth = .25, colour = "#6F7D86"),
  axis.ticks.length = unit(1, "mm"),
  panel.grid.major.y = element_line(linewidth = .2, colour = "#E8EDF0"),
  plot.title = element_text(size = 7.5, face = "bold", lineheight = 1.08, margin = margin(b = 3)),
  plot.subtitle = element_text(size = 7, colour = "#52616A", margin = margin(b = 5)),
  plot.tag = element_text(size = 8, face = "bold"),
  legend.title = element_blank(), legend.text = element_text(size = 7),
  legend.key.width = unit(8, "mm"), legend.key.height = unit(3, "mm"),
  legend.margin = margin(0, 0, 0, 0), legend.position = "bottom",
  plot.background = element_blank(), panel.background = element_blank(),
  legend.background = element_blank(), legend.key = element_blank(),
  plot.margin = margin(5, 5, 5, 5))
annotation_theme <- theme(text = element_text(family = "Arial"),
  plot.title = element_text(size = 10, face = "bold", margin = margin(b = 4)),
  plot.subtitle = element_text(size = 7.5, lineheight = 1.2, margin = margin(b = 6)),
  plot.caption = element_text(size = 6.5, lineheight = 1.2, hjust = 0, margin = margin(t = 6)),
  plot.background = element_blank(),
  plot.margin = margin(6, 7, 6, 7))

# Use final Cairo-resolved viewports, preserving the strict 1.5-point gate.
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
    seekViewport(sprintf("%s.%d-%d-%d-%d", row$name, row$t, row$l, row$b, row$r))
    lo <- deviceLoc(unit(0, "npc"), unit(0, "npc"), valueOnly = TRUE)
    hi <- deviceLoc(unit(1, "npc"), unit(1, "npc"), valueOnly = TRUE)
    upViewport(0); c(lo$x, lo$y, hi$x, hi$y) * 72
  }
  panels <- lapply(seq_len(nrow(pr)), function(i) {
    row <- pr[i, ]; b <- bounds(row)
    stopifnot(b[3] > b[1], b[4] > b[2])
    tag <- g$layout[g$layout$name == sub("^panel", "tag", row$name), ]
    ans <- list(id = panel_ids[i], bbox_pt = unname(b), grid_id = "patchwork-grid-1",
      row_start = row$t - 1L, row_stop = row$b, col_start = row$l - 1L, col_stop = row$r)
    if (nrow(tag) == 1L) {tb <- bounds(tag); ans$panel_label_anchor_pt <- unname(c(tb[1], tb[4]))}
    ans
  })
  manifest <- list(schema_version = 1L, backend = "r-patchwork-cairo-resolved-viewports",
    figure = list(width_pt = width_in * 72, height_pt = height_in * 72), panels = panels, exemptions = exemptions)
  if (!is.null(row_groups)) manifest$row_groups <- lapply(row_groups, as.character)
  if (!is.null(column_groups)) manifest$column_groups <- lapply(column_groups, as.character)
  write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  invisible(manifest)
}
export_figure <- function(p, stem, height_mm) {
  w <- 183 / 25.4; h <- height_mm / 25.4; dst <- file.path(figdir, stem); qa <- file.path(qadir, stem)
  require_patchwork_panel_alignment(p, paste0(qa, ".alignment-layout.json"),
    paste0(qa, ".alignment.json"), width_in = w, height_in = h,
    audit_script = file.path(skill, "scripts/audit_panel_alignment.py"), python = python,
    overlay_svg = paste0(qa, ".alignment.svg"), tolerance_pt = 1.5, gutter_tolerance_pt = 1.5, strict = TRUE)
  cairo_pdf(paste0(dst, ".pdf"), width = w, height = h, family = font); print(p); dev.off()
  svglite(paste0(dst, ".svg"), width = w, height = h); print(p); dev.off()
  agg_tiff(paste0(dst, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw")
  print(p); dev.off()
  agg_png(paste0(dst, ".png"), width = w, height = h, units = "in", res = 300); print(p); dev.off()
  write_json(list(backend = "R", width_mm = 183, height_mm = height_mm,
    tiff_dpi = 600, png_dpi = 300, font = font, glyph_floor_pt = 6.5,
    interpolation = FALSE, renormalization = FALSE, refitting = FALSE), paste0(qa, ".export.json"), auto_unbox = TRUE, pretty = TRUE)
}

prob <- subset(marginal_raw, score %in% c(0, 2) & grepl("^endpoint_S[0-4]$", metric))
stopifnot(nrow(prob) == 54L, all(prob$m == 20L), !anyDuplicated(prob[c("cohort", "score", "horizon", "metric")]),
  all(prob$conf_low <= prob$estimate), all(prob$estimate <= prob$conf_high),
  all(prob$conf_low >= 0), all(prob$conf_high <= 1))
prob$series <- factor(ifelse(prob$score == 0, "Exactly 0 symptoms", "Exactly 2 symptoms"), levels = names(score_colours))
prob$display_estimate <- prob$estimate * 100; prob$display_low <- prob$conf_low * 100; prob$display_high <- prob$conf_high * 100
prob$display_unit <- "Percent"
prob$state <- sub("endpoint_", "", prob$metric)
prob$display_x <- prob$horizon + ifelse(prob$score == 0, -.12, .12)
prob$x_offset_note <- "Symmetric 0.12-year visual dodge separates intervals at each source horizon; no interpolation."
state_titles <- c(S0 = "S0 | No target\ndisease", S1 = "S1 | Hypertension\nonly",
  S2 = "S2 | Heart disease\nor stroke only", S3 = "S3 | Cardiovascular\nmultimorbidity", S4 = "S4 | Death\n(absorbing)")
ylimits <- list(S0 = c(0, 100), S1 = c(0, 35), S2 = c(0, 18), S3 = c(0, 11), S4 = c(0, 22))
ybreaks <- list(S0 = c(0, 25, 50, 75, 100), S1 = c(0, 10, 20, 30), S2 = c(0, 5, 10, 15), S3 = c(0, 5, 10), S4 = c(0, 10, 20))
panels3 <- list()
for (cohort_name in c("HRS", "ELSA")) for (s in paste0("S", 0:4)) {
  d <- prob[prob$cohort == cohort_name & prob$state == s, ]
  p <- ggplot(d, aes(x = display_x, y = display_estimate, colour = series, shape = series, group = series)) +
    scale_x_continuous(breaks = c(2, 5, 10), limits = c(1.25, 10.75), expand = c(0, 0)) +
    scale_y_continuous(limits = ylimits[[s]], breaks = ybreaks[[s]], expand = expansion(mult = c(0, .04))) +
    scale_colour_manual(values = score_colours, drop = FALSE) + scale_shape_manual(values = score_shapes, drop = FALSE) +
    labs(title = state_titles[[s]], subtitle = cohort_name, x = "Years", y = if (s == "S0") "Probability (%)" else NULL) + theme_pub()
  if (nrow(d)) {
    stopifnot(nrow(d) == 6L, max(d$display_high) < ylimits[[s]][2])
    p <- p + geom_line(linewidth = .4) +
      geom_errorbar(aes(ymin = display_low, ymax = display_high), width = .35, linewidth = .4) +
      geom_point(size = 1.5, stroke = .3)
  } else {
    stopifnot(cohort_name == "ELSA", s == "S4")
    p <- p + annotate("text", x = 6, y = 12.7,
      label = "Death was not\nrepresented in the\nELSA primary model.\n\nNot estimable", size = 6.5 / ggplot2::.pt,
      family = font, colour = "#52616A", lineheight = 1.2) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid = element_blank())
  }
  panels3[[length(panels3) + 1L]] <- p
}
fig3 <- wrap_plots(panels3, ncol = 5, guides = "collect") +
  plot_annotation(title = "State probabilities under maintained insomnia scores",
    subtitle = "From S0 | Nuisance covariates fixed at fitted means | State-specific probability scales",
    caption = "Points: estimates with pointwise 95% CIs from 20 imputations. Lines connect the 2-, 5- and 10-year estimates only.\nHRS includes death; ELSA excludes death. Fixed-score conditional projections are not intervention effects.",
    tag_levels = "a", theme = annotation_theme) & theme(legend.position = "bottom")
write.csv(prob, file.path(srcdir, "Figure3_Markov_state_probabilities.csv"), row.names = FALSE, fileEncoding = "UTF-8")
export_figure(fig3, "Figure3_Markov_state_probabilities", 151)

dif <- subset(paired_raw, metric %in% c("ever_S3", "endpoint_S3", "S3_free_years"))
stopifnot(nrow(dif) == 18L, all(dif$m == 20L), !anyDuplicated(dif[c("cohort", "horizon", "metric")]),
  all(dif$conf_low <= dif$estimate), all(dif$estimate <= dif$conf_high))
dif$display_multiplier <- ifelse(dif$metric == "S3_free_years", 1, 100)
dif$display_estimate <- dif$estimate * dif$display_multiplier
dif$display_low <- dif$conf_low * dif$display_multiplier; dif$display_high <- dif$conf_high * dif$display_multiplier
dif$display_unit <- ifelse(dif$metric == "S3_free_years", "Years", "Percentage points")
dif$display_x <- dif$horizon + ifelse(dif$cohort == "HRS", -.1, .1)
dif$x_offset_note <- "Symmetric 0.1-year visual dodge separates cohort intervals at source horizons."
delta_panel <- function(d, title, subtitle, years = FALSE) {
  p <- ggplot(d, aes(x = display_x, y = display_estimate, colour = cohort, shape = cohort, group = cohort)) +
    geom_hline(yintercept = 0, colour = "#89969E", linetype = "dashed", linewidth = .35) +
    geom_line(linewidth = .45) + geom_errorbar(aes(ymin = display_low, ymax = display_high), width = .32, linewidth = .5) +
    geom_point(size = 1.9, stroke = .4) +
    scale_colour_manual(values = cohort_colours, limits = c("HRS", "ELSA"), drop = FALSE) +
    scale_shape_manual(values = cohort_shapes, limits = c("HRS", "ELSA"), drop = FALSE) +
    scale_x_continuous(breaks = c(2, 5, 10), limits = c(1.2, 10.8), expand = c(0, 0)) +
    labs(title = title, subtitle = subtitle, x = "Projection horizon (years)",
      y = if (years) "Difference (years)" else "Difference (percentage points)") + theme_pub() +
    theme(plot.title = element_text(size = 8), plot.subtitle = element_text(size = 7))
  if (years) p <- p + scale_y_continuous(limits = c(-.43, .025), breaks = c(-.4, -.3, -.2, -.1, 0), expand = c(0, 0))
  else p <- p + scale_y_continuous(limits = c(-.22, 3.5), breaks = 0:3, expand = c(0, 0))
  p
}
p4a <- delta_panel(subset(dif, metric == "ever_S3"), "Ever reaches cardiovascular multimorbidity", "Cumulative probability within the horizon")
p4b <- delta_panel(subset(dif, metric == "endpoint_S3"), "In cardiovascular multimorbidity at the horizon", "Endpoint probability; HRS requires being alive")
p4c <- delta_panel(subset(dif, metric == "S3_free_years" & cohort == "HRS"), "HRS | Restricted time without multimorbidity", "Expected time alive and in S0-S2", TRUE)
p4d <- delta_panel(subset(dif, metric == "S3_free_years" & cohort == "ELSA"), "ELSA | Restricted time without multimorbidity", "Expected time in S0-S2; mortality not represented", TRUE)
p4c <- p4c + guides(colour = "none", shape = "none")
p4d <- p4d + guides(colour = "none", shape = "none")
fig4 <- wrap_plots(p4a, p4b, p4c, p4d, ncol = 2, guides = "collect") +
  plot_annotation(title = "Paired absolute differences: exactly 2 versus 0 symptoms",
    subtitle = "From S0 | Maintained scores | Other covariates fixed at fitted means",
    caption = "Points and bars: paired differences with pointwise 95% CIs; complete covariance retained within each imputation.\nELSA cumulative and endpoint S3 probabilities coincide by construction. Cross-cohort levels have different mortality scope.",
    tag_levels = "a", theme = annotation_theme) & theme(legend.position = "bottom")
write.csv(dif, file.path(srcdir, "Figure4_Markov_absolute_differences.csv"), row.names = FALSE, fileEncoding = "UTF-8")
export_figure(fig4, "Figure4_Markov_absolute_differences", 153)
stopifnot(identical(hashes, vapply(paths, function(p) digest(file = p, algo = "sha256"), "")))
write_json(list(source_files = as.list(paths), sha256 = as.list(hashes),
  preservation = "Source hashes unchanged after rendering", backend = "R", fitting = FALSE,
  filters = list(Figure3 = "score 0 or 2; all endpoint states; all stored 2/5/10-year horizons", Figure4 = "ever_S3, endpoint_S3, S3_free_years; both cohorts; all stored horizons"),
  row_counts = list(marginal_input = 264L, marginal_displayed = 54L, paired_input = 66L, paired_displayed = 18L),
  transformations = c("Probability multiplied by 100 for percent or percentage-point display", "Years unchanged", "Small horizontal dodges only; no fitted/interpolated values", "No renormalization"),
  figure_contract = list(archetype = "quantitative grid", Figure3 = "State decomposition of fixed-score projections with cohort mortality boundary", Figure4 = "Paired absolute consequences for entry into multimorbidity and restricted time")),
  file.path(qadir, "19b_source_provenance.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message("Figure 3 and Figure 4 exported. PDF glyph and collision QA required before final delivery.")
