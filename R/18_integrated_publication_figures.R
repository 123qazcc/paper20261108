#!/usr/bin/env Rscript
# Presentation only: final aggregate estimates are read without fitting or pooling.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid); library(jsonlite)
  library(svglite); library(ragg); library(digest)
})
out <- "outputs/Three_cohort_sleep_CVD_20261001"
figdir <- file.path(out, "figures")
srcdir <- file.path(out, "source_data")
qadir <- "outputs/three_cohort_integrated_20261001/figure_qa"
for (d in c(figdir, srcdir, qadir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
skill <- Sys.getenv("NATURE_FIGURE_SKILL_ROOT", "private/render_tools/nature-figure")
python <- Sys.getenv("NATURE_FIGURE_PYTHON", "python")
source(file.path(skill, "scripts/panel_alignment.R"))
font <- "Arial"
colours <- c(HRS = "#275F85", ELSA = "#B27036", CHARLS = "#337C6E", CORE = "#617585", BOTH = "#337C6E")
shapes <- c(HRS = 16, ELSA = 17, CHARLS = 15, CORE = 1, BOTH = 16)
txt <- 8 / ggplot2::.pt
small <- 7.5 / ggplot2::.pt
paths <- c(
  main = "outputs/hrs_elsa_multistate_final_20260909/results/primary_multistate_results_exact.csv",
  specific = "outputs/hrs_elsa_multistate_nextsteps_20260909/results/disease_specific/disease_specific_target_edge_hr_mi20.csv",
  direct = "outputs/hrs_elsa_multistate_nextsteps_20260909/results/disease_specific/disease_specific_heart_vs_stroke_contrasts_mi20.csv",
  charls = "outputs/charls_external_extension_20260913/results/four_target_MI20.csv",
  charls_direct = "outputs/charls_external_extension_20260913/results/heart_stroke_contrasts_MI20.csv",
  joint = "outputs/charls_cesd9_duration_20260914/results/four_target_MI20.csv",
  absolute = "outputs/hrs_elsa_multistate_nextsteps_20260909/results/absolute_risk/fixed_score_2_vs_0_paired_differences.csv"
)
hashes <- vapply(paths, function(p) digest(file = p, algo = "sha256"), "")
read_source <- function(key) read.csv(paths[[key]], check.names = FALSE, fileEncoding = "UTF-8")
save_source <- function(d, name) write.csv(d, file.path(srcdir, name), row.names = FALSE, fileEncoding = "UTF-8", na = "")
check_estimates <- function(d, n, positive = TRUE) {
  stopifnot(nrow(d) == n, all(is.finite(as.matrix(d[c("estimate", "lower", "upper")]))),
    all(d$lower <= d$estimate), all(d$estimate <= d$upper), !anyDuplicated(d[c("group", "series")]))
  if (positive) stopifnot(all(d$lower > 0))
  invisible(d)
}
theme_pub <- function() theme_void(base_size = 8, base_family = "Arial") + theme(
  text = element_text(family = font, colour = "#26323A"),
  plot.title = element_text(size = 10, face = "bold", margin = margin(b = 5)),
  plot.subtitle = element_text(size = 8, lineheight = 1.2, colour = "#52616A", margin = margin(b = 8)),
  plot.caption = element_text(size = 7.5, lineheight = 1.2, hjust = 0, margin = margin(t = 8)),
  plot.tag = element_text(size = 9, face = "bold"),
  plot.margin = margin(7, 9, 7, 9), legend.position = "none")

# Forest-table columns use a linear physical layout; estimates are mapped by
# the declared log or identity transform into the dedicated forest region.
forest_table <- function(d, labels, series_order, limits, ticks, title = NULL,
    subtitle = NULL, caption = NULL, axis_title = "Transition intensity ratio (log scale)",
    estimate_header = "Ratio (95% CI)", log_scale = TRUE, reference = 1,
    number_digits = 3, signed = FALSE, group_header = "Transition", series_header = "Cohort") {
  stopifnot(all(d$group %in% names(labels)), all(d$series %in% series_order),
    all(d$lower >= limits[1]), all(d$upper <= limits[2]))
  if (log_scale) stopifnot(all(d$lower > 0), all(limits > 0))
  transform <- if (log_scale) log else identity
  map <- function(x) 75 + 46 * (transform(x) - transform(limits[1])) / diff(transform(limits))
  ng <- length(labels); ns <- length(series_order)
  d$y <- (ng - match(d$group, names(labels))) * (ns + 1) + ns - match(d$series, series_order)
  d$x <- map(d$estimate); d$xl <- map(d$lower); d$xu <- map(d$upper)
  number_format <- if (signed) paste0("%+.", number_digits, "f (%+.", number_digits, "f to %+.", number_digits, "f)") else
    paste0("%.", number_digits, "f (%.", number_digits, if (log_scale) "f–%." else "f to %.", number_digits, "f)")
  d$numeric <- sprintf(number_format, d$estimate, d$lower, d$upper)
  centers <- data.frame(group = names(labels), label = unname(labels),
    y = (ng - seq_len(ng)) * (ns + 1) + (ns - 1) / 2)
  ymax <- max(d$y); header_y <- ymax + 1.35
  tick <- data.frame(x = map(ticks), label = format(ticks, trim = TRUE, scientific = FALSE))
  sep <- if (ng > 1) data.frame(y = (ng - seq_len(ng - 1)) * (ns + 1) - 1) else data.frame(y = numeric())
  p <- ggplot() +
    geom_segment(data = sep, aes(x = 0, xend = 173, y = y, yend = y), colour = "#E5EAEE", linewidth = .25) +
    annotate("segment", x = map(reference), xend = map(reference), y = -.5, yend = ymax + .6,
      colour = "#9AA6AD", linetype = "dashed", linewidth = .35) +
    geom_errorbar(data = d, aes(y = y, xmin = xl, xmax = xu, colour = series),
      orientation = "y", width = .16, linewidth = .5) +
    geom_point(data = d, aes(x = x, y = y, colour = series, shape = series), size = 2.1, stroke = .6) +
    geom_text(data = centers, aes(x = 0, y = y, label = label), hjust = 0, size = txt,
      lineheight = 1.15, family = font) +
    geom_text(data = d, aes(x = 54, y = y, label = series, colour = series), hjust = 0,
      size = txt, family = font) +
    geom_text(data = d, aes(x = 130, y = y, label = numeric), hjust = 0, size = txt,
      family = font) +
    annotate("text", x = c(0, 54, 130), y = header_y,
      label = c(group_header, series_header, estimate_header), hjust = 0, size = txt,
      fontface = "bold", family = font) +
    annotate("segment", x = 0, xend = 173, y = ymax + .78, yend = ymax + .78,
      linewidth = .4, colour = "#89979F") +
    annotate("segment", x = 75, xend = 121, y = -.95, yend = -.95, linewidth = .3) +
    geom_segment(data = tick, aes(x = x, xend = x, y = -.95, yend = -1.10), linewidth = .3) +
    geom_text(data = tick, aes(x = x, y = -1.5, label = label), size = small, family = font) +
    annotate("text", x = 98, y = -2.2, label = axis_title, size = small, family = font) +
    scale_colour_manual(values = colours) + scale_shape_manual(values = shapes) +
    scale_x_continuous(limits = c(-.5, 174), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-2.65, ymax + 1.9), expand = c(0, 0)) +
    labs(title = title, subtitle = subtitle, caption = caption) + theme_pub()
  p
}

# Measure resolved Cairo viewports at final physical size. The upstream strict
# alignment checks and default 1.5-point tolerance are unchanged.
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
    upViewport(0)
    c(lo$x, lo$y, hi$x, hi$y) * 72
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
  w <- 183 / 25.4; h <- height_mm / 25.4
  dst <- file.path(figdir, stem); qa <- file.path(qadir, stem)
  if (inherits(p, "patchwork")) {
    require_patchwork_panel_alignment(p, paste0(qa, ".alignment-layout.json"),
      paste0(qa, ".alignment.json"), width_in = w, height_in = h,
      audit_script = file.path(skill, "scripts/audit_panel_alignment.py"), python = python,
      overlay_svg = paste0(qa, ".alignment.svg"), strict = TRUE)
  } else write_json(list(verdict = "NOT APPLICABLE", applicable = FALSE,
    reason = "Single plot area", panel_count = 1), paste0(qa, ".alignment.json"), auto_unbox = TRUE, pretty = TRUE)
  cairo_pdf(paste0(dst, ".pdf"), width = w, height = h, family = font); print(p); dev.off()
  svglite(paste0(dst, ".svg"), width = w, height = h); print(p); dev.off()
  agg_tiff(paste0(dst, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw")
  print(p); dev.off()
  agg_png(paste0(dst, ".png"), width = w, height = h, units = "in", res = 300)
  print(p); dev.off()
  write_json(list(backend = "R", width_mm = 183, height_mm = height_mm,
    tiff_dpi = 600, png_dpi = 300, font = font, rows_preserved = TRUE),
    paste0(qa, ".export.json"), auto_unbox = TRUE, pretty = TRUE)
  invisible(dst)
}

main_raw <- read_source("main")
main <- subset(main_raw, method == "Markov" & model == "M4_MI20" & confirmatory_family == "nonfatal_four_edge")
main$group <- main$transition; main$series <- main$cohort
main$estimate <- main$hazard_ratio; main$lower <- main$conf_low; main$upper <- main$conf_high
check_estimates(main, 8)
save_source(main, "Figure1_HRS_ELSA_main.csv")
main_labels <- c("S0->S1" = "No target disease\n→ HTN only",
  "S0->S2" = "No target disease\n→ Heart disease or stroke only",
  "S1->S3" = "HTN only\n→ Cardiovascular multimorbidity",
  "S2->S3" = "Heart disease or stroke only\n→ Cardiovascular multimorbidity")
fig1 <- forest_table(main, main_labels, c("HRS", "ELSA"), c(.85, 1.4), c(.9, 1, 1.2, 1.4),
  title = "Insomnia burden and cardiovascular disease accumulation",
  subtitle = "Four nonfatal transitions | Per additional insomnia symptom (0–3 score)",
  caption = "HRS: 29,029 participants; ELSA: 11,097 participants. Points and lines show ratios and 95% CIs.\nHTN, hypertension. Target diseases: hypertension, heart disease and stroke.")
export_figure(fig1, "Figure1_HRS_ELSA_main_transitions", 108)

specific_raw <- read_source("specific"); charls_raw <- read_source("charls")
specific <- data.frame(group = paste0("P", specific_raw$target_order), series = specific_raw$cohort,
  estimate = specific_raw$hazard_ratio, lower = specific_raw$conf_low, upper = specific_raw$conf_high,
  p = specific_raw$p_value, exposure_unit = "One additional insomnia symptom", source_transition = specific_raw$transition)
specific <- rbind(specific, data.frame(group = paste0("P", match(charls_raw$transition, c("0-2", "0-4", "1-3", "1-5"))),
  series = "CHARLS", estimate = charls_raw$HR, lower = charls_raw$lower, upper = charls_raw$upper,
  p = charls_raw$p, exposure_unit = "One higher restless-sleep frequency category", source_transition = charls_raw$transition))
check_estimates(specific, 12)
direct_raw <- read_source("direct"); charls_direct_raw <- read_source("charls_direct")
direct <- data.frame(group = ifelse(direct_raw$origin_state_label == "No disease", "G0", "G1"),
  series = direct_raw$cohort, estimate = direct_raw$hr_ratio_heart_vs_stroke,
  lower = direct_raw$conf_low, upper = direct_raw$conf_high, p = direct_raw$p_value,
  p_BH_two = direct_raw$p_bh_two_contrasts_within_cohort,
  exposure_unit = "One additional insomnia symptom")
direct <- rbind(direct, data.frame(group = ifelse(charls_direct_raw$origin == "No disease", "G0", "G1"),
  series = "CHARLS", estimate = charls_direct_raw$HR, lower = charls_direct_raw$lower, upper = charls_direct_raw$upper,
  p = charls_direct_raw$p, p_BH_two = charls_direct_raw$p_BH_2,
  exposure_unit = "One higher restless-sleep frequency category"))
check_estimates(direct, 6)
save_source(specific, "Figure2a_disease_specific_paths.csv"); save_source(direct, "Figure2b_direct_heart_stroke_ratios.csv")
target_labels <- c(P1 = "No target disease\n→ Heart disease only", P2 = "No target disease\n→ Stroke only",
  P3 = "HTN only\n→ HTN + heart disease", P4 = "HTN only\n→ HTN + stroke")
p2a <- forest_table(specific, target_labels, c("HRS", "ELSA", "CHARLS"), c(.7, 1.6), c(.8, 1, 1.2, 1.6),
  title = "Disease-specific transition associations")
p2b <- forest_table(direct, c(G0 = "No target disease", G1 = "HTN only"), c("HRS", "ELSA", "CHARLS"),
  c(.7, 1.6), c(.8, 1, 1.2, 1.6), title = "Direct heart disease versus stroke comparison",
  axis_title = "Heart-disease ratio / stroke ratio (log scale)", group_header = "Origin state")
fig2 <- (p2a / p2b) + plot_layout(heights = c(1.75, 1)) +
  plot_annotation(title = "Disease-specific associations across three cohorts",
    subtitle = "HRS/ELSA: per additional insomnia symptom; CHARLS: per higher restless-sleep frequency category.\nDifferent exposure units: qualitative cross-cohort comparison only; no three-cohort pooling.",
    caption = "Points and lines show ratios and 95% CIs. In panel b, ratios above 1 favour a larger heart-disease coefficient.\nDirect comparisons use within-model covariance. HRS stroke includes TIA; cohort disease definitions differ.",
    tag_levels = "a", theme = theme_pub()) & theme(plot.tag = element_text(size = 9, face = "bold", family = font))
export_figure(fig2, "Figure2_three_cohort_disease_specific", 210)

joint <- read_source("joint")
joint$group <- paste0("P", match(joint$transition, c("0-2", "0-4", "1-3", "1-5")))
joint$series <- joint$model; joint$estimate <- joint$HR
check_estimates(joint, 8)
save_source(joint, "Figure3_CHARLS_CORE_BOTH.csv")
fig3 <- forest_table(joint, target_labels, c("CORE", "BOTH"), c(.8, 1.25), c(.8, .9, 1, 1.1, 1.25),
  title = "Restless sleep after joint adjustment in CHARLS",
  subtitle = "12,036 participants | Per higher restless-sleep frequency category (0–3)", series_header = "Model",
  caption = "CORE: primary adjustment. BOTH: CORE plus CESD-9 and sleep duration; the same stage-2 imputations.\nPoints and lines show ratios and 95% CIs. Joint adjustment does not establish mediation.")
export_figure(fig3, "Figure3_CHARLS_joint_adjustment", 108)

absolute_raw <- read_source("absolute")
absolute <- subset(absolute_raw, horizon == 10 & metric %in% c("endpoint_S3", "ever_S3", "S3_free_years"))
absolute$group <- absolute$metric; absolute$series <- absolute$cohort
absolute$lower <- absolute$conf_low; absolute$upper <- absolute$conf_high
absolute$source_estimate <- absolute$estimate; absolute$source_lower <- absolute$lower; absolute$source_upper <- absolute$upper
probability <- absolute$metric != "S3_free_years"
absolute[probability, c("estimate", "lower", "upper")] <- absolute[probability, c("estimate", "lower", "upper")] * 100
absolute$display_unit <- ifelse(probability, "Percentage points", "Years")
check_estimates(absolute, 6, positive = FALSE)
save_source(absolute, "Figure6_absolute_scale_differences.csv")
p6a <- forest_table(absolute[probability, ], c(endpoint_S3 = "In S3 at year 10", ever_S3 = "Ever reaches S3\nwithin 10 years"),
  c("HRS", "ELSA"), c(-.1, 3.5), 0:3, title = "Probability differences", log_scale = FALSE,
  reference = 0, number_digits = 2, signed = FALSE, estimate_header = "Difference (95% CI)",
  group_header = "Projected measure", axis_title = "Difference (percentage points)")
p6b <- forest_table(absolute[!probability, ], c(S3_free_years = "Alive and S3-free time\nover 10 years"),
  c("HRS", "ELSA"), c(-.42, .02), c(-.4, -.3, -.2, -.1, 0), title = "Restricted-time differences", log_scale = FALSE,
  reference = 0, number_digits = 3, signed = FALSE, estimate_header = "Difference (95% CI)",
  group_header = "Projected measure", axis_title = "Difference (years)")
fig6 <- (p6a / p6b) + plot_layout(heights = c(1.4, 1)) +
  plot_annotation(title = "Absolute-scale model projections at 10 years",
    subtitle = "From S0 | Maintained insomnia scores of 2 versus 0 | Other covariates fixed at fitted means",
    caption = "HRS includes death; ELSA excludes mortality and is not a competing-death risk projection.\nAlive-and-S3-free time can be reduced by death as well as S3 entry. These projections are not intervention effects.",
    tag_levels = "a", theme = theme_pub()) & theme(plot.tag = element_text(size = 9, face = "bold", family = font))
export_figure(fig6, "Figure6_HRS_ELSA_absolute_scale", 136)

stopifnot(identical(hashes, vapply(paths, function(p) digest(file = p, algo = "sha256"), "")))
write_json(list(source_files = as.list(paths), sha256 = as.list(hashes),
  preservation = "Input hashes unchanged after render", backend = "R", fitting = FALSE,
  filters = list(main_source_rows = nrow(main_raw), main_selected_rows = nrow(main),
    disease_specific_selected_rows = nrow(specific), direct_selected_rows = nrow(direct),
    joint_selected_rows = nrow(joint), absolute_source_rows = nrow(absolute_raw), absolute_selected_rows = nrow(absolute))),
  file.path(qadir, "source_provenance.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message("Export complete; final PDF text, collision, and source audits are required before delivery.")
