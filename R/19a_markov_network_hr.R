#!/usr/bin/env Rscript
# Presentation-only extraction of existing aggregate MI20 Markov results.
if (.Platform$OS.type == "windows") Sys.setlocale("LC_CTYPE", "English_United States.utf8")
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid)
  library(jsonlite); library(svglite); library(ragg); library(digest)
})
out <- "outputs/Three_cohort_Markov_results_20261001"
figdir <- file.path(out, "figures")
srcdir <- file.path(out, "source_data")
qadir <- "outputs/three_cohort_markov_revision_20261001/figure_qa"
for (d in c(figdir, srcdir, qadir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
skill <- "private/render_tools/nature-figure"
python <- "python"
font <- "Arial"
source(file.path(skill, "scripts/panel_alignment.R"))
# Reuse only the resolved-Cairo geometry helper from the existing R style.
helper_expr <- parse("R/18_integrated_publication_figures.R", encoding = "UTF-8")
for (expr in helper_expr) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("write_patchwork_panel_layout"))) eval(expr)
}
txt <- 7.6 / ggplot2::.pt
small <- 7.2 / ggplot2::.pt
blue <- "#275F85"; orange <- "#A46C32"; grey <- "#78828A"; ink <- "#26323A"
theme_pub <- function() theme_void(base_size = 8, base_family = font) +
  theme(text = element_text(family = font, colour = ink),
    plot.title = element_text(size = 9, face = "bold", margin = margin(b = 2)),
    plot.subtitle = element_text(size = 7.6, colour = "#52616A", margin = margin(b = 4)),
    plot.caption = element_text(size = 7.4, lineheight = 1.15, hjust = 0, margin = margin(t = 6)),
    plot.tag = element_text(size = 9, face = "bold", family = font),
    plot.margin = margin(5, 7, 5, 7), legend.position = "none")
export_figure <- function(p, stem, height_mm) {
  w <- 183 / 25.4; h <- height_mm / 25.4
  dst <- file.path(figdir, stem); qa <- file.path(qadir, stem)
  require_patchwork_panel_alignment(p, paste0(qa, ".alignment-layout.json"),
    paste0(qa, ".alignment.json"), width_in = w, height_in = h,
    panel_ids = c("a", "b", "c"), column_groups = list(c("a", "b", "c")),
    audit_script = file.path(skill, "scripts/audit_panel_alignment.py"), python = python,
    overlay_svg = paste0(qa, ".alignment.svg"), tolerance_pt = 1.5, gutter_tolerance_pt = 1.5, strict = TRUE)
  cairo_pdf(paste0(dst, ".pdf"), width = w, height = h, family = font); print(p); dev.off()
  svglite(paste0(dst, ".svg"), width = w, height = h); print(p); dev.off()
  agg_tiff(paste0(dst, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw")
  print(p); dev.off()
  agg_png(paste0(dst, ".png"), width = w, height = h, units = "in", res = 300)
  print(p); dev.off()
  write_json(list(backend = "R", width_mm = 183, height_mm = height_mm,
    tiff_dpi = 600, png_dpi = 300, font = font),
    paste0(qa, ".export.json"), auto_unbox = TRUE, pretty = TRUE)
}
root <- "outputs/hrs_elsa_multistate_final_20260909/results"
source_files <- c(q = file.path(root, "markov_q_intensities.csv"),
  hr = file.path(root, "markov_hazard_ratios.csv"),
  fit = file.path(root, "markov_fit_metrics.csv"),
  wald = file.path(root, "markov_confirmatory_joint_wald.csv"))
hashes <- vapply(source_files, function(p) digest(file = p, algo = "sha256"), "")
read_source <- function(k) read.csv(source_files[[k]], fileEncoding = "UTF-8", check.names = FALSE)
q <- subset(read_source("q"), model_id == "M4_MI20" & covariate_profile == "sample means")
hr <- subset(read_source("hr"), model_id == "M4_MI20" & term == "insomnia_score_primary")
fit <- subset(read_source("fit"), model_id == "M4_MI20")
wald <- subset(read_source("wald"), model_id == "M4_MI20")
ids <- c("hrs_5state_full", "elsa_4state_w4_w10", "elsa_5state_through_2012")
edges <- c("S0->S1", "S0->S2", "S1->S3", "S2->S3", "S0->S4", "S1->S4", "S2->S4", "S3->S4")
for (d in list(q, hr)) {
  stopifnot(nrow(d) == 20, !anyDuplicated(d[c("analysis_id", "transition")]),
    all(d$m_expected == 20), all(d$m_fitted == 20), all(d$pooling_status == "RUBIN_POOLED_COMPLETE"))
  for (i in seq_along(ids)) stopifnot(setequal(d$transition[d$analysis_id == ids[i]],
    if (i == 2) edges[1:4] else edges))
}
stopifnot(all(is.finite(q$intensity)), all(q$conf_low > 0),
  all(q$conf_low <= q$intensity), all(q$intensity <= q$conf_high),
  all(is.finite(hr$hazard_ratio)), all(hr$conf_low > 0),
  all(hr$conf_low <= hr$hazard_ratio), all(hr$hazard_ratio <= hr$conf_high))
q$display_value <- sprintf("%.4f /year", q$intensity)
q$intensity_unit <- "per year"
q$interval_definition <- "95% CI; delta variance within imputation and Rubin pooling of log intensity"
hr$display_value <- sprintf("%.3f (%.3f-%.3f)", hr$hazard_ratio, hr$conf_low, hr$conf_high)
hr$p_scope <- ifelse(grepl("S4", hr$transition), "Raw P; mortality edges are exploratory",
  "Raw P and BH-adjusted P across four nonfatal transitions within analysis")
write.csv(q, file.path(srcdir, "Figure1_Markov_transition_networks.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(hr, file.path(srcdir, "Figure2_Markov_complete_HR.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(wald, file.path(srcdir, "Figure2_Markov_joint_Wald.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")

titles <- c("HRS: five-state primary model",
  "ELSA: four-state primary model",
  "ELSA: mortality available only through 2012")
subtitles <- c("29,029 participants | 8,984 observed deaths | Follow-up through 2020",
  "11,097 participants | Waves 4-10 | Death is not included",
  "9,244 participants | 186 observed deaths | Separate five-state extension")

# Explicit planar routing: progressive disease edges occupy the central diamond;
# mortality edges enter S4 via separate outer corridors. No missing direct edge is drawn.
nodes5 <- data.frame(state = 0:4, x = c(12, 60, 60, 112, 164), y = c(0, 17, -17, 0, 0),
  label = c("S0\nNo disease", "S1\nHTN only", "S2\nOne CVD", "S3\nMultimorbidity", "S4\nDeath"))
routes <- list(
  "S0->S1" = rbind(c(25, 2), c(47, 15)),
  "S0->S2" = rbind(c(25, -2), c(47, -15)),
  "S1->S3" = rbind(c(73, 15), c(99, 2)),
  "S2->S3" = rbind(c(73, -15), c(99, -2)),
  "S0->S4" = rbind(c(12, -9), c(12, -36), c(164, -36), c(164, -9)),
  "S1->S4" = rbind(c(60, 26), c(60, 34), c(164, 34), c(164, 9)),
  "S2->S4" = rbind(c(73, -17), c(128, -17), c(151, -3)),
  "S3->S4" = rbind(c(125, 0), c(151, 0)))
positions <- data.frame(transition = edges, x = c(34, 34, 86, 86, 85, 113, 118, 139),
  y = c(22, -22, 22, 0, -42, 41, -26, 9))
network <- function(i) {
  d <- q[q$analysis_id == ids[i], ]
  n <- if (i == 2) nodes5[1:4, ] else nodes5
  colour <- if (i == 1) blue else orange
  d <- merge(d, positions, by = "transition", sort = FALSE)
  d$edge_colour <- ifelse(d$death_transition == 1, grey, colour)
  path <- do.call(rbind, lapply(d$transition, function(e) {
    xy <- routes[[e]]; data.frame(transition = e, x = xy[, 1], y = xy[, 2],
      death = grepl("S4", e), colour = if (grepl("S4", e)) grey else colour)
  }))
  p <- ggplot() +
    geom_path(data = path[!path$death, ], aes(x, y, group = transition), linewidth = .55,
      colour = colour, lineend = "round", arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
    geom_path(data = path[path$death, ], aes(x, y, group = transition), linewidth = .5,
      colour = grey, linetype = "dashed", lineend = "butt",
      arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
    geom_rect(data = n, aes(xmin = x - 13, xmax = x + 13, ymin = y - 9, ymax = y + 9),
      fill = "#F3F6F7", colour = "#A6B1B7", linewidth = .35) +
    geom_text(data = n, aes(x, y = y + 1, label = label), size = small, family = font, lineheight = 1.02) +
    geom_text(data = d, aes(x, y, label = display_value, colour = edge_colour),
      size = small, family = font) +
    scale_colour_identity() +
    scale_x_continuous(limits = c(-3, 179), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-48, 48), expand = c(0, 0)) +
    labs(title = titles[i], subtitle = subtitles[i]) + theme_pub()
  if (i == 2) p <- p + annotate("text", x = 154, y = 0, label = "Four nonfatal\ntransitions",
    size = txt, colour = "#52616A", family = font, lineheight = 1.25)
  p
}
fig1 <- (network(1) / network(2) / network(3)) +
  plot_annotation(title = "Markov transition intensities",
    subtitle = "Adjusted intensity q at the sample-mean covariate profile | All 20 allowed transitions",
    caption = paste("Solid coloured arrows: nonfatal transitions. Dashed grey arrows: death transitions.",
      "Labels show q per year; these are intensities, not probabilities. No direct S0 to S3 intensity.",
      "HTN, hypertension; CVD, heart disease or stroke. S2 excludes HTN. 95% CIs are in Source Data.",
      sep = "\n"), tag_levels = "a", theme = theme_pub()) &
  theme(plot.tag = element_text(size = 9, face = "bold", family = font))
export_figure(fig1, "Figure1_Markov_transition_networks", 188)

# Identical logarithmic effect mapping and columns in all three panels.
limits <- c(.65, 2.4)
stopifnot(all(limits > 0), all(hr$conf_low > 0), all(hr$hazard_ratio > 0))
map <- function(x) 78 + 50 * (log(x) - log(limits[1])) / diff(log(limits))
labels <- c("S0->S1" = "S0 to S1: No disease to HTN",
  "S0->S2" = "S0 to S2: No disease to one CVD",
  "S1->S3" = "S1 to S3: HTN to multimorbidity",
  "S2->S3" = "S2 to S3: CVD to multimorbidity",
  "S0->S4" = "S0 to S4: No disease to death",
  "S1->S4" = "S1 to S4: HTN to death",
  "S2->S4" = "S2 to S4: CVD to death",
  "S3->S4" = "S3 to S4: Multimorbidity to death")
forest <- function(i) {
  d <- hr[hr$analysis_id == ids[i], ]
  d <- d[match(if (i == 2) edges[1:4] else edges, d$transition), ]
  n <- nrow(d); d$y <- rev(seq_len(n))
  d$x <- map(d$hazard_ratio); d$xl <- map(d$conf_low); d$xu <- map(d$conf_high)
  d$colour <- ifelse(grepl("S4", d$transition), grey, if (i == 1) blue else orange)
  d$label <- labels[d$transition]
  tick <- data.frame(x = map(c(.75, 1, 1.5, 2)), label = c("0.75", "1.0", "1.5", "2.0"))
  p <- ggplot() +
    annotate("segment", x = map(1), xend = map(1), y = .5, yend = n + .4,
      linewidth = .35, linetype = "dashed", colour = "#A6B1B7") +
    geom_errorbar(data = d, aes(y = y, xmin = xl, xmax = xu, colour = colour),
      orientation = "y", width = .18, linewidth = .5) +
    geom_point(data = d, aes(x, y, colour = colour), size = 1.8) +
    geom_text(data = d, aes(x = 0, y = y, label = label), hjust = 0, size = txt, family = font) +
    geom_text(data = d, aes(x = 134, y = y, label = display_value), hjust = 0, size = txt, family = font) +
    annotate("text", x = c(0, 134), y = n + 1.2, label = c("Transition", "HR (95% CI)"),
      hjust = 0, size = txt, fontface = "bold", family = font) +
    annotate("segment", x = 0, xend = 178, y = n + .5, yend = n + .5, colour = "#A6B1B7", linewidth = .35) +
    annotate("segment", x = 78, xend = 128, y = .22, yend = .22, linewidth = .3) +
    geom_segment(data = tick, aes(x = x, xend = x, y = .22, yend = .08), linewidth = .3) +
    geom_text(data = tick, aes(x = x, y = -.3, label = label), size = small, family = font) +
    scale_colour_identity() +
    scale_x_continuous(limits = c(-1, 180), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-.65, n + 1.65), expand = c(0, 0)) +
    labs(title = titles[i], subtitle = subtitles[i]) + theme_pub()
  if (i != 2) p <- p + annotate("segment", x = 0, xend = 178, y = 4.5, yend = 4.5,
    colour = "#E0E5E8", linewidth = .25)
  p
}
fig2 <- (forest(1) / forest(2) / forest(3)) + plot_layout(heights = c(9.85, 5.85, 9.85)) +
  plot_annotation(title = "Insomnia and every allowed Markov transition",
    subtitle = "Hazard ratio per additional insomnia symptom (0-3 score) | Fully adjusted MI20 models",
    caption = paste("Points and lines: pooled HR and 95% CI; horizontal axes use the same logarithmic scale.",
      "Coloured estimates: nonfatal transitions. Grey estimates: exploratory mortality transitions.",
      "No significance symbols are used. Joint tests and multiplicity details are reported in the figure legend.",
      sep = "\n"), tag_levels = "a", theme = theme_pub()) &
  theme(plot.tag = element_text(size = 9, face = "bold", family = font))
export_figure(fig2, "Figure2_Markov_complete_HR", 178)

legend1 <- paste(
"Fig. 1 | Adjusted Markov transition intensities.",
"a, HRS five-state primary analysis (29,029 participants; 8,984 observed deaths).",
"b, ELSA waves 4-10 four-state primary analysis (11,097 participants), which excludes mortality.",
"c, ELSA five-state extension with mortality available only through 2012 (9,244 participants; 186 observed deaths).",
"Each arrow is one allowed direct transition; coloured solid arrows are nonfatal and grey dashed arrows lead to death.",
"Labels report q per year at the fitted model-matrix sample means, including the insomnia score.",
"Intensities are not finite-horizon probabilities. S0 to S3 has no direct intensity; progression can occur through S1 or S2.",
"All 20 estimates use complete 20-imputation Rubin pooling on the log-intensity scale.",
"95% confidence limits and the complete observed event support are provided in Source Data.",
"HTN, hypertension; CVD, heart disease or stroke. S2 denotes one CVD without HTN; S3 denotes cardiovascular multimorbidity.",
"No model was refitted for this figure.")
legend2 <- paste(
"Fig. 2 | Insomnia associations with all allowed Markov transitions.",
"a-c, HRS primary, ELSA nonfatal primary and ELSA mortality-limited extension, respectively; sample sizes are shown in the panels.",
"Points show transition-specific hazard ratios per additional insomnia symptom (0-3 score); horizontal lines show 95% confidence intervals from complete MI20 Rubin pooling.",
"All panels use the same logarithmic HR scale. Coloured estimates represent nonfatal transitions; grey estimates represent exploratory death transitions.",
"Models adjust for age, sex, race, education, wealth, marital status, smoking, drinking, inactivity, BMI, diabetes, lung disease and frequent pain.",
"Insomnia effects vary by transition; nuisance effects are constrained common across transitions.",
"Four-edge nonfatal joint Wald tests use Rubin total covariance and a large-sample chi-square reference (4 df): HRS, chi-square=50.007, P=3.60e-10; ELSA primary, 17.791, P=0.00136; ELSA mortality-limited extension, 15.041, P=0.00462.",
"Source Data include raw P values and BH-adjusted P values within each four-edge nonfatal family. Death-edge P values are unadjusted and exploratory.",
"No significance symbols are used. ELSA mortality follow-up ends in 2012.")
writeLines(c(legend1, "", legend2), file.path(srcdir, "Figure1_2_legends.txt"), useBytes = TRUE)
write_json(list(source_files = as.list(source_files), sha256 = as.list(hashes),
  filters = list(q = "model_id=M4_MI20 and covariate_profile=sample means",
    hr = "model_id=M4_MI20 and term=insomnia_score_primary"),
  displayed_rows = c(network = nrow(q), forest = nrow(hr)),
  figure_contract = list(network = "Quantitative network: primary model structure, external nonfatal comparison, and mortality-window boundary",
    forest = "Quantitative grid: all transition-specific insomnia associations including null and mortality effects"),
  adaptation = "Structural adaptation of existing R forest-table styling and resolved-Cairo alignment helper; new network geometry"),
  file.path(qadir, "markov_figure_provenance.json"), auto_unbox = TRUE, pretty = TRUE)
message("Rendered Figure1 and Figure2; final PDF glyph and collision QA required.")

