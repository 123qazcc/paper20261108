#!/usr/bin/env Rscript
# Figure 4: descriptive observed exposure flow, not disease transition risks.
# Quantitative single-panel alluvial: all four selected-cohort groups retained.
source("R/17j_dynamic_extra_figure_helpers.R", encoding = "UTF-8")

pattern_file <- file.path(root, "stage3_delivery/table2_dynamic_pattern_distribution.csv")
pattern_sha <- sha(pattern_file)
flow <- fread(pattern_file)
expected_order <- c("persistent_low", "onset", "remission", "persistent_high")
expected_counts <- c(4702L, 1369L, 1483L, 1843L)
stopifnot(nrow(flow) == 4L, uniqueN(flow$pattern) == 4L,
          setequal(flow$pattern, expected_order))
flow <- flow[match(expected_order, pattern)]
stopifnot(identical(as.integer(flow$people), expected_counts),
          sum(flow$people) == 9397L,
          all(abs(flow$percent - 100 * flow$people / sum(flow$people)) < 1e-10))
flow[, `:=`(pattern_label = c("Persistent low", "Onset", "Remission", "Persistent high"),
            w1 = c("Low", "Low", "High", "High"),
            w2 = c("Low", "High", "Low", "High"))]
n_total <- sum(flow$people)
gap <- .065 * n_total
high_left <- sum(flow[w1 == "High", people])
high_right <- sum(flow[w2 == "High", people])
# Every ribbon's vertical thickness is its exact participant count at both
# wave nodes. The blank node gap is constant and does not represent people.
flow[, left_low := c(high_left + gap + people[2], high_left + gap,
                    people[4], 0)]
flow[, right_low := c(high_right + gap + people[3], people[4],
                     high_right + gap, 0)]
flow[, `:=`(left_high = left_low + people, right_high = right_low + people)]
stopifnot(all(flow$left_high - flow$left_low == flow$people),
          all(flow$right_high - flow$right_low == flow$people))

curve_x <- seq(0, 1, length.out = 201)
ease <- function(x) x^2 * (3 - 2*x)
ribbons <- rbindlist(lapply(seq_len(nrow(flow)), function(i) {
  z <- flow[i]
  lower <- z$left_low + (z$right_low - z$left_low) * ease(curve_x)
  upper <- z$left_high + (z$right_high - z$left_high) * ease(curve_x)
  data.table(pattern = z$pattern, x = c(curve_x, rev(curve_x)),
             y = c(lower, rev(upper)))
}))
ribbons[, pattern := factor(pattern, levels = c("persistent_low", "persistent_high", "remission", "onset"))]
flow[, `:=`(label_x = .205,
            label_y = (left_low + left_high)/2 + (right_low - left_low)*ease(.205))]
flow[, label := sprintf("%s\n%s (%.1f%%)", pattern_label,
                       format(people, big.mark = ",", trim = TRUE), percent)]
nodes <- data.table(wave = rep(c("W1", "W2"), each = 2),
                    frequency = rep(c("Low", "High"), 2),
                    people = c(n_total-high_left, high_left, n_total-high_right, high_right),
                    xmin = rep(c(-.018, .982), each = 2),
                    xmax = rep(c(.018, 1.018), each = 2),
                    ymin = c(high_left + gap, 0, high_right + gap, 0),
                    ymax = c(n_total + gap, high_left, n_total + gap, high_right))
nodes[, label := sprintf("%s frequency\n%s (%.1f%%)", frequency,
                         format(people, big.mark = ",", trim = TRUE), 100*people/n_total)]
palette_flow <- c(persistent_low = "#B7D0DE", onset = "#DDB477",
                  remission = "#8FB7AC", persistent_high = "#BF7769")
plot4 <- ggplot() +
  geom_polygon(data = ribbons, aes(x, y, group = pattern, fill = pattern),
               colour = "white", linewidth = .25, alpha = .96) +
  geom_rect(data = nodes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#52616D", colour = NA) +
  geom_text(data = flow, aes(label_x, label_y, label = label),
            family = font, colour = "#17242B", size = 2.75, lineheight = 1.2) +
  geom_text(data = nodes[wave == "W1"],
            aes(x = -.055, y = (ymin+ymax)/2, label = label),
            hjust = 1, family = font, colour = "#26343E", size = 2.55, lineheight = 1.22) +
  geom_text(data = nodes[wave == "W2"],
            aes(x = 1.055, y = (ymin+ymax)/2, label = label),
            hjust = 0, family = font, colour = "#26343E", size = 2.55, lineheight = 1.22) +
  annotate("text", x = c(0, 1), y = n_total + gap + .065*n_total,
           label = c("W1", "W2 landmark"), family = font, fontface = "bold", size = 3.0) +
  scale_fill_manual(values = palette_flow, guide = "none") +
  coord_cartesian(xlim = c(-.37, 1.37), ylim = c(-.025*n_total, n_total + gap + .11*n_total),
                  expand = FALSE, clip = "off") +
  labs(title = "Observed changes in restless-sleep frequency",
       subtitle = "CHARLS W1-W2 landmark cohort | 9,397 participants with observed exposure at both waves",
       caption = paste("Ribbon thickness represents participant counts; percentages use the full landmark cohort (n = 9,397).",
         "High frequency: restless sleep on >=3 days in the past week. Low frequency: <3 days.",
         "Descriptive selected-cohort flows, not cardiovascular transitions or population probabilities. No model CIs apply.", sep = "\n")) +
  theme_void(base_family = font, base_size = 7) +
  theme(plot.title = element_text(size = 9, face = "bold", margin = margin(b = 5)),
        plot.subtitle = element_text(size = 7, colour = "#586772", margin = margin(b = 8)),
        plot.caption = element_text(size = 6.4, hjust = 0, lineheight = 1.2, margin = margin(t = 8)),
        plot.margin = margin(9, 10, 9, 10), plot.background = element_rect(fill = "white", colour = NA))
export_single(plot4, "figure4_pattern_flow", 127)
flow[, `:=`(cohort_n = n_total, width_unit = "participants", node_gap_count_units = gap)]
fwrite(flow, file.path(dest, "figure4_source_data.csv"), bom = TRUE)
assert_scientific_unchanged()
stopifnot(sha(pattern_file) == pattern_sha)
write_json(list(status = "RENDERED_REQUIRES_VISUAL_AND_PDF_QA", figure = 4L,
  backend = "R", archetype = "single-panel quantitative alluvial",
  inputs_sha256 = pattern_sha, n = n_total, patterns = nrow(flow), rows_excluded = 0L,
  left_total = sum(nodes[wave == "W1", people]), right_total = sum(nodes[wave == "W2", people]),
  exact_count_widths = TRUE, uncertainty = "Not applicable: descriptive cohort census",
  scientific_inputs_unchanged = TRUE, disease_probability_claim = FALSE,
  source = "R/17k_dynamic_pattern_flow_EN.R"),
  file.path(dest, "qa/figure4_render_checks.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Figure 4 rendered: four observed patterns; counts, inputs and frozen sources verified.\n")
