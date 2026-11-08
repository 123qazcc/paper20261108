#!/usr/bin/env Rscript
# Empirical observed endpoint counts; never uses a fitted Q matrix.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggalluvial); library(patchwork)
  library(digest); library(jsonlite)
})
options(stringsAsFactors = FALSE)
input_root <- file.path("outputs", "hrs_elsa_multistate_final_20260909", "data_multistate")
output_root <- file.path("outputs", "hrs_elsa_multistate_nextsteps_20260909")
result_root <- file.path(output_root, "results", "empirical_flows")
data_root <- file.path(output_root, "data_empirical_flows")
figure_root <- file.path(output_root, "figures", "empirical_flows")
source_root <- file.path(figure_root, "source_data")
qc_root <- file.path(output_root, "qc", "empirical_flows")
report_root <- file.path(output_root, "reports", "empirical_flows")
for (d in c(result_root, data_root, figure_root, source_root, qc_root, report_root)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
script_version <- "2026-09-11-empirical-flows-v2-visible-ribbons"
state_levels <- c("S0", "S1", "S2", "S3", "Death", "Not observed")
state_labels <- c(S0 = "No target disease", S1 = "Hypertension only", S2 = "One CVD without hypertension",
                  S3 = "Cardiovascular multimorbidity", Death = "Recorded death (HRS)", `Not observed` = "Not observed")
state_colours <- c(S0 = "#77A5A5", S1 = "#D6AC64", S2 = "#718FB3", S3 = "#B77671", Death = "#595E65", `Not observed` = "#D5D8DC")
checks <- list(); summaries <- list(); audits <- list(); grids <- list(); composition <- list(); flows <- list(); patterns <- list(); refs <- list()
qc <- function(cohort, check, pass, details = "") {
  checks[[length(checks)+1L]] <<- data.table(cohort = cohort, check = check, pass = isTRUE(pass), details = as.character(details))
  if (!isTRUE(pass)) stop(cohort, ": ", check, " failed: ", details)
}
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

for (cohort_name in c("HRS", "ELSA")) {
  stem <- tolower(cohort_name)
  pp <- file.path(input_root, paste0(stem, "_msm_panel.rds"))
  ap <- file.path(input_root, paste0(stem, "_exposure_panel_augmented.rds"))
  p <- as.data.table(readRDS(pp)); a <- as.data.table(readRDS(ap))
  waves <- if (cohort_name == "HRS") 8:15 else 4:10
  entry_wave_fixed <- min(waves)
  entry <- unique(p[, .(id, entry_wave, analysis_entry_time = exposure_entry_time)])
  qc(cohort_name, "one_entry_per_participant", !anyDuplicated(entry$id))
  qc(cohort_name, "all_final_participants_in_entry_table", nrow(entry) == uniqueN(p$id))
  visits <- a[event_type == "visit" & wave %in% waves]
  visits <- merge(visits, entry, by = "id", all = FALSE)
  visits <- visits[wave >= entry_wave & event_time >= analysis_entry_time - 1e-8]
  qc(cohort_name, "actual_visits_unique_id_wave", !anyDuplicated(visits[, .(id, wave)]))
  visits[, complete_disease_items := complete.cases(hypertension_raw, heart_disease_raw, stroke_raw)]
  visits[, observed_state := fifelse(complete_disease_items & state %in% 0:3, paste0("S", state), "Not observed")]
  visits[, not_observed_reason := fifelse(observed_state == "Not observed", "Visit present; disease items incomplete", NA_character_)]
  cross <- merge(p[event_type == "visit" & wave %in% waves, .(id, wave, frozen_state = state)],
                 visits[, .(id, wave, state)], by = c("id", "wave"), all.x = TRUE)
  qc(cohort_name, "augmented_visit_state_matches_final_panel", all(!is.na(cross$state) & cross$state == cross$frozen_state))
  wave_ref <- visits[, .(reference_time = median(event_time), min_visit_time = min(event_time), max_visit_time = max(event_time),
                         n_actual_visits = .N), by = wave]
  setorder(wave_ref, wave)
  qc(cohort_name, "all_waves_present", identical(wave_ref$wave, as.integer(waves)))
  wave_ref[, `:=`(cohort = cohort_name, reference_rule = "Median actual post-entry visit time; only places missing-visit HRS death")]
  refs[[cohort_name]] <- wave_ref
  deaths <- p[event_type == "death" & state == 4L, .(death_time = min(event_time)), by = id]
  grid <- merge(entry[, join_key := 1L], data.table(wave = waves, join_key = 1L), by = "join_key", allow.cartesian = TRUE)[, join_key := NULL]
  entry[, join_key := NULL]
  grid <- grid[wave >= entry_wave]
  grid <- merge(grid, wave_ref[, .(wave, reference_time)], by = "wave", all.x = TRUE)
  grid <- merge(grid, visits[, .(id, wave, observed_state, not_observed_reason, event_time,
                                complete_disease_items, hypertension_raw, heart_disease_raw, stroke_raw,
                                frozen_history_state = state)], by = c("id", "wave"), all.x = TRUE)
  grid <- merge(grid, deaths, by = "id", all.x = TRUE)
  grid[, actual_visit := as.integer(!is.na(event_time))]
  grid[is.na(observed_state), `:=`(observed_state = "Not observed", not_observed_reason = "No visit in supplied analysis panel")]
  grid[, known_death_before_reference := !is.na(death_time) & death_time <= reference_time]
  grid[, carried_recorded_death := FALSE]
  if (cohort_name == "HRS") grid[actual_visit == 0L & known_death_before_reference, `:=`(observed_state = "Death", not_observed_reason = NA_character_, carried_recorded_death = TRUE)]
  grid[, `:=`(cohort = cohort_name, fixed_entry_cohort = as.integer(entry_wave == entry_wave_fixed),
              state_display = factor(observed_state, levels = state_levels))]
  setorder(grid, id, wave)
  grid[, previous_state := shift(observed_state), by = id]
  grid[, next_state := shift(observed_state, type = "lead"), by = id]
  grid[, next_wave := shift(wave, type = "lead"), by = id]
  grid[, returned_after_not_observed := as.integer(observed_state %in% paste0("S", 0:3) & previous_state == "Not observed")]
  qc(cohort_name, "disease_state_requires_actual_complete_visit", all(grid[observed_state %in% paste0("S", 0:3), actual_visit == 1L & complete_disease_items]))
  qc(cohort_name, "no_living_state_after_recorded_death", grid[previous_state == "Death" & observed_state != "Death", .N] == 0L)
  if (cohort_name == "ELSA") qc(cohort_name, "ELSA_no_death_created", !any(grid$observed_state == "Death"))
  fixed <- grid[fixed_entry_cohort == 1L]
  n_fixed <- uniqueN(fixed$id)
  qc(cohort_name, "fixed_cohort_all_waves_equal_denominator", all(fixed[, .N, by = wave]$N == n_fixed))
  qc(cohort_name, "fixed_cohort_person_wave_complete", nrow(fixed) == n_fixed * length(waves))
  denom <- grid[, .(n_ever_entered = .N, n_actual_visit = sum(actual_visit),
                    n_observed_complete = sum(observed_state %in% paste0("S", 0:3)),
                    n_not_observed = sum(observed_state == "Not observed"),
                    n_no_visit = sum(observed_state == "Not observed" & actual_visit == 0L),
                    n_incomplete_disease = sum(observed_state == "Not observed" & actual_visit == 1L),
                    n_recorded_death = sum(observed_state == "Death"), n_new_entrants = sum(entry_wave == wave),
                    n_known_limited_ELSA_death_in_not_observed = sum(cohort_name == "ELSA" & observed_state == "Not observed" & known_death_before_reference),
                    n_returned = sum(returned_after_not_observed, na.rm = TRUE)), by = wave]
  denom[, cohort := cohort_name]
  qc(cohort_name, "wave_denominator_conservation", all(denom$n_ever_entered == denom$n_observed_complete + denom$n_not_observed + denom$n_recorded_death))
  comp <- grid[observed_state %in% paste0("S", 0:3), .(n = .N), by = .(wave, observed_state)]
  comp <- merge(CJ(wave = waves, observed_state = paste0("S", 0:3)), comp, by = c("wave", "observed_state"), all.x = TRUE)
  comp[is.na(n), n := 0L]
  comp <- merge(comp, denom, by = "wave", all.x = TRUE)
  comp[, observed_proportion := n / n_observed_complete]
  qc(cohort_name, "observed_composition_sums_one", all(abs(comp[, sum(observed_proportion), by = wave]$V1 - 1) < 1e-12))
  fc <- fixed[!is.na(next_wave), .(n = .N), by = .(wave, next_wave, from_state = observed_state, to_state = next_state)]
  fc[, `:=`(cohort = cohort_name, n_fixed = n_fixed, proportion_fixed = n/n_fixed,
            observed_direct_jump = from_state == "S0" & to_state == "S3")]
  qc(cohort_name, "all_adjacent_pair_flows_conserve_fixed_n", all(fc[, sum(n), by = wave]$V1 == n_fixed))
  wide <- dcast(fixed[, .(id, wave, observed_state)], id ~ wave, value.var = "observed_state")
  setnames(wide, as.character(waves), paste0("W", waves))
  pattern_cols <- paste0("W", waves)
  pattern <- wide[, .(n = .N), by = pattern_cols]
  pattern[, trajectory_id := .I]
  pattern[, `:=`(cohort = cohort_name, n_fixed = n_fixed)]
  pattern_long <- melt(pattern, id.vars = c("trajectory_id", "cohort", "n", "n_fixed"), measure.vars = pattern_cols,
                       variable.name = "wave_name", value.name = "observed_state")
  pattern_long[, `:=`(wave = as.integer(sub("W", "", wave_name)), proportion_fixed = n/n_fixed,
                      state_display = factor(observed_state, levels = state_levels))]
  patterns[[cohort_name]] <- pattern_long
  fwrite(wide, file.path(data_root, paste0(stem, "_fixed_entry_state_wide.csv.gz")))
  saveRDS(grid, file.path(data_root, paste0(stem, "_observed_person_wave_grid.rds")), compress = "gzip")
  fwrite(grid, file.path(data_root, paste0(stem, "_observed_person_wave_grid.csv.gz")), na = "")
  fwrite(pattern, file.path(source_root, paste0(stem, "_fixed_entry_trajectory_counts.csv")))
  fwrite(pattern_long, file.path(source_root, paste0(stem, "_fixed_entry_trajectory_lodes.csv.gz")))
  fwrite(denom, file.path(result_root, paste0(stem, "_wave_denominators.csv")))
  fixed_counts <- fixed[, .(n = .N), by = .(wave, observed_state)]
  fixed_counts[, `:=`(cohort = cohort_name, n_fixed = n_fixed, proportion_fixed = n/n_fixed)]
  fwrite(fixed_counts, file.path(source_root, paste0(stem, "_fixed_entry_state_counts.csv")))
  composition[[cohort_name]] <- comp; flows[[cohort_name]] <- fc; grids[[cohort_name]] <- grid
  summaries[[cohort_name]] <- data.table(cohort = cohort_name, final_analysis_participants = nrow(entry),
    fixed_entry_wave = entry_wave_fixed, fixed_entry_n = n_fixed,
    delayed_entry_n = sum(entry$entry_wave != entry_wave_fixed), total_grid_cells = nrow(grid),
    actual_post_entry_visits = sum(grid$actual_visit), incomplete_disease_visits = sum(grid$actual_visit == 1L & grid$observed_state == "Not observed"),
    fixed_returns_n = fixed[returned_after_not_observed == 1L, .N], fixed_returning_people = fixed[returned_after_not_observed == 1L, uniqueN(id)],
    fixed_terminal_not_observed = fixed[wave == max(waves) & observed_state == "Not observed", .N],
    fixed_terminal_recorded_death = fixed[wave == max(waves) & observed_state == "Death", .N],
    fixed_direct_S0_S3_jumps = fixed[observed_state == "S0" & next_state == "S3", .N],
    supplied_death_records = nrow(deaths), death_after_last_wave_reference = sum(deaths$death_time > max(wave_ref$reference_time)),
    n_unique_trajectories = nrow(pattern))
  audits[[cohort_name]] <- data.table(cohort = cohort_name, path_role = c("frozen_analysis_panel", "augmented_visit_panel"),
    sha256 = c(digest(file = pp, algo = "sha256"), digest(file = ap, algo = "sha256")),
    input_rows = c(nrow(p), nrow(a)), input_participants = c(uniqueN(p$id), uniqueN(a$id)))
}
summary_table <- rbindlist(summaries)
fwrite(summary_table, file.path(result_root, "empirical_flow_summary.csv"))
fwrite(rbindlist(refs), file.path(result_root, "wave_reference_times.csv"))
fwrite(rbindlist(audits), file.path(qc_root, "input_fingerprints.csv"))
fwrite(rbindlist(composition), file.path(source_root, "observed_wave_composition.csv"))
fwrite(rbindlist(flows), file.path(source_root, "fixed_entry_adjacent_flow_counts.csv"))
fwrite(rbindlist(checks), file.path(qc_root, "empirical_flow_numerical_checks.csv"))

base_theme <- function() theme_classic(base_size = 7, base_family = "Helvetica") + theme(
  axis.line = element_line(linewidth = .3, colour = "#333333"),
  axis.ticks = element_line(linewidth = .3), axis.text = element_text(size = 6.3, colour = "#333333"),
  axis.title = element_text(size = 7), legend.position = "bottom", legend.text = element_text(size = 6.5),
  legend.title = element_blank(), legend.key.size = grid::unit(3.4, "mm"),
  plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 6.5),
  plot.tag = element_text(size = 9, face = "bold"), plot.margin = margin(5, 5, 5, 5),
  panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = .25))
comp_plots <- list(); flow_plots <- list()
for (cohort_name in c("HRS", "ELSA")) {
  dt <- composition[[cohort_name]]; pl <- patterns[[cohort_name]]
  waves <- sort(unique(dt$wave))
  dt[, state_display := factor(observed_state, levels = state_levels)]
  dt[, wave_factor := factor(wave, levels = waves)]
  pl[, wave_factor := factor(wave, levels = waves)]
  labs <- unique(dt[, .(wave, n_observed_complete, n_not_observed)])
  setorder(labs, wave)
  axis_labs <- paste0("W", labs$wave, "\nn=", fmt(labs$n_observed_complete), "\nN/O=", fmt(labs$n_not_observed))
  comp_plots[[cohort_name]] <- ggplot(dt, aes(x = wave_factor, y = observed_proportion, fill = state_display)) +
    geom_col(width = .68, colour = "white", linewidth = .18) +
    scale_fill_manual(values = state_colours, limits = state_levels[1:4], labels = state_labels[1:4], drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 1, .25), labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, .015))) + coord_cartesian(ylim = c(0, 1)) +
    scale_x_discrete(labels = axis_labs) +
    labs(x = NULL, y = "Share of observed disease states", title = paste0(cohort_name, ": observed state composition"),
         subtitle = "n: complete observed visits; N/O: not observed among those already entered") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + base_theme()
  ss <- summary_table[cohort == cohort_name]
  flow_plots[[cohort_name]] <- ggplot(pl, aes(x = wave_factor, stratum = state_display,
                                       alluvium = trajectory_id, y = proportion_fixed, fill = state_display)) +
    # GeomFlow removes complete.cases(data), including colour; NA silently
    # drops every ribbon. A transparent nonmissing outline preserves counts.
    geom_flow(width = .12, alpha = .6, colour = "transparent", aes.flow = "forward", show.legend = FALSE) +
    geom_stratum(width = .12, colour = "white", linewidth = .1, show.legend = TRUE) +
    scale_fill_manual(values = state_colours, limits = state_levels, labels = state_labels, drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 1, .25), labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, .015))) + coord_cartesian(ylim = c(0, 1)) +
    scale_x_discrete(labels = paste0("W", waves), expand = expansion(mult = c(.025, .025))) +
    labs(x = "Survey wave (spacing does not represent elapsed time)", y = "Share of fixed entry cohort",
         title = paste0(cohort_name, ": fixed W", ss$fixed_entry_wave, " entry cohort, n=", fmt(ss$fixed_entry_n)),
         subtitle = if (cohort_name == "HRS") "Recorded death carried forward; missed visits remain Not observed" else "W4-W10: mortality incompletely ascertained; loss never coded as death") +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE)) + base_theme()
}
for (nm in names(flow_plots)) {
  rendered_flow <- ggplot_build(flow_plots[[nm]])$data[[1L]]
  qc(nm, "flow_layer_nonmissing_and_positive_width", nrow(rendered_flow) > 0L &&
     all(complete.cases(rendered_flow)) && all(rendered_flow$ymax > rendered_flow$ymin))
}
fwrite(rbindlist(checks), file.path(qc_root, "empirical_flow_numerical_checks.csv"))
fig_comp <- comp_plots$HRS / comp_plots$ELSA + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
fig_flow <- flow_plots$HRS / flow_plots$ELSA + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")

write_patchwork_panel_layout_resolved <- function(plot, manifest_path, width_in, height_in,
  panel_ids = NULL, row_groups = NULL, column_groups = NULL, exemptions = list()) {
  grob <- patchwork::patchworkGrob(plot)
  cells <- grob$layout[grepl("^panel(?:-[0-9]+)?$", grob$layout$name, perl = TRUE), , drop = FALSE]
  cells <- cells[order(cells$t, cells$l), , drop = FALSE]
  if (is.null(panel_ids)) panel_ids <- letters[seq_len(nrow(cells))]
  stopifnot(nrow(cells) >= 2, length(panel_ids) == nrow(cells))
  probe <- tempfile(fileext = ".pdf")
  grDevices::pdf(probe, width = width_in, height = height_in, family = "Helvetica", useDingbats = FALSE)
  on.exit({grDevices::dev.off(); unlink(probe)}, add = TRUE)
  grid::grid.newpage(); grid::grid.draw(grob); grid::grid.force()
  resolve <- function(units, total, horizontal) {
    null <- grid::unitType(units) == "null"
    v <- if (horizontal) grid::convertWidth(units, "pt", valueOnly = TRUE) else grid::convertHeight(units, "pt", valueOnly = TRUE)
    if (any(null)) {
      remaining <- total - sum(v[!null]); weights <- as.numeric(units[null])
      stopifnot(remaining > 0, sum(weights) > 0)
      v[null] <- remaining * weights/sum(weights)
    }
    stopifnot(all(is.finite(v)), all(v >= 0)); v
  }
  xe <- c(0, cumsum(resolve(grob$widths, width_in*72, TRUE)))
  ye <- c(0, cumsum(resolve(grob$heights, height_in*72, FALSE)))
  panels <- lapply(seq_len(nrow(cells)), function(i) {
    z <- cells[i, ]
    list(id = panel_ids[i], bbox_pt = unname(c(xe[z$l], height_in*72-ye[z$b+1L], xe[z$r+1L], height_in*72-ye[z$t])),
      grid_id = "empirical-flow-grid", row_start = as.integer(z$t-1L), row_stop = as.integer(z$b),
      col_start = as.integer(z$l-1L), col_stop = as.integer(z$r))
  })
  manifest <- list(schema_version = 1L, backend = "r-patchwork-resolved-null-tracks",
    figure = list(width_pt = width_in*72, height_pt = height_in*72), panels = panels, exemptions = exemptions)
  if (!is.null(row_groups)) manifest$row_groups <- lapply(row_groups, as.character)
  if (!is.null(column_groups)) manifest$column_groups <- lapply(column_groups, as.character)
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  invisible(manifest)
}
export <- function(plot, base, height_mm) {
  width <- 183/25.4; height <- height_mm/25.4
  source("R/helpers/rendered_panel_alignment.R", local = environment())
  require_rendered_panel_alignment(plot, base, width, height)
  grDevices::pdf(paste0(base, ".pdf"), width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  print(plot); dev.off()
  svglite::svglite(paste0(base, ".svg"), width = width, height = height); print(plot); dev.off()
  ragg::agg_tiff(paste0(base, ".tiff"), width = width, height = height, units = "in", res = 600, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(paste0(base, ".png"), width = width, height = height, units = "in", res = 300); print(plot); dev.off()
}
export(fig_comp, file.path(figure_root, "fig_observed_wave_composition"), 180)
export(fig_flow, file.path(figure_root, "fig_fixed_entry_empirical_flows"), 185)

writeLines(c(
  "Fig. | Observed disease-state composition by survey wave.",
  "a, HRS W8-W15. b, ELSA W4-W10. Stacks show the proportions of S0-S3 among actual post-entry visits with complete current hypertension, heart-disease and stroke items. States retain the frozen cumulative diagnosis-history definition at those observed visits. Wave labels give complete observed n and Not observed (N/O) n among all final-analysis participants who have entered by that wave; N/O includes absent visits and incomplete disease items. HRS recorded deaths are excluded from both the observed denominator and N/O after they are placed on the wave grid; ELSA mortality is incompletely ascertained, so N/O may include deaths. Delayed entrants join the observed-composition display after entry. These are unweighted descriptive sample counts, without confidence intervals or model-based transition inference. Full denominators, including known HRS deaths, are provided in Source Data.",
  "", "图：按访视波次展示实际观察者的疾病状态组成。", "每波n是当波三项疾病信息完整的真实访视人数，N/O是已入组但该波未观察到完整状态的人数。延迟入组者自入组起纳入。HRS已记录死亡按波次参考时间归入死亡后不计入观察者或N/O；ELSA的N/O可能含未查明死亡。各波柱形仅描述观察者组成，不代表固定队列患病概率。"
), file.path(figure_root, "fig_observed_wave_composition_caption.md"), useBytes = TRUE)
writeLines(c(
  "Fig. | Empirical adjacent-wave flows in fixed entry cohorts.",
  "a, HRS participants entering the final analysis at W8. b, ELSA participants entering at W4. Every initial participant contributes exactly one endpoint at every displayed wave; ribbon widths are paired participant counts divided by the fixed cohort size. Not observed denotes an absent visit or incomplete current disease items and permits later return to an observed state. At actual complete visits, frozen cumulative disease-history states are used; living disease states are never filled across unobserved waves. HRS death is assigned to a missing-visit wave only when recorded death precedes that wave's median visit date and is then absorbing; a recorded living visit takes precedence at its own wave. ELSA has no displayed death state because ascertainment ends in 2012. Direct observed S0-S3 jumps remain direct, without inventing an intermediate disease order. Ribbons describe adjacent-wave paired flows, not identifiable continuous trajectories between visits. Wave spacing is categorical, not elapsed time; ELSA W10 spans 2021-2023. No fitted transition intensities, uncertainty intervals or significance tests are used.",
  "", "图：固定入组队列的相邻波次经验状态流。", "HRS固定W8、ELSA固定W4入组者，每人在每波均占一个位置。缺失访视或不完整疾病信息单列Not observed，允许其后返回。HRS仅按已有死亡时间携带吸收死亡状态；ELSA不将失访推定为死亡。疾病状态不跨缺失波次填补，同一间隔观察到的多病跳跃保持端点，不假定发病顺序。流带为真实个体相邻波次端点的聚合计数，不是拟合转移概率。"
), file.path(figure_root, "fig_fixed_entry_empirical_flows_caption.md"), useBytes = TRUE)

report <- c("# 经验状态组成与固定队列状态流", "", "本步骤已完成两类互补展示：逐波实际观察者的疾病组成，以及固定入组队列的相邻波次个体状态流。使用冻结分析数据，不重新拟合多状态模型。", "", "## 数据及口径")
for (i in seq_len(nrow(summary_table))) {
  z <- summary_table[i]
  report <- c(report, sprintf("- %s：最终分析参与者%s人；固定W%d入组队列%s人，晚于该波入组%s人。延迟入组者参与逐波观察者组成，自入组前不计入分母；只在固定队列流图中排除。", z$cohort, fmt(z$final_analysis_participants), z$fixed_entry_wave, fmt(z$fixed_entry_n), fmt(z$delayed_entry_n)),
    sprintf("  固定队列共有%s种完整端点模式，缺失后返回%s次（%s人）；相邻已观察访视中S0直接进入S3的端点跳跃%s次。末波Not observed为%s人，已记录死亡为%s人。", fmt(z$n_unique_trajectories), fmt(z$fixed_returns_n), fmt(z$fixed_returning_people), fmt(z$fixed_direct_S0_S3_jumps), fmt(z$fixed_terminal_not_observed), fmt(z$fixed_terminal_recorded_death)))
}
report <- c(report, "", "## 如何读图", "", "组成图的柱形分母是每波三项疾病信息完整的实际访视者，n/N/O直接标在各波下方。它不是固定队列的患病率曲线，不能把后期观察者中较高的多病比例全归为疾病进展。固定队列流图保持人数不变，缺失位置单列，因此可直接看到返回访视、失访和已记录死亡如何影响观察者构成。", "", "流图的带宽来自相邻波次真实参与者配对计数；没有从拟合q矩阵制造流量。相邻端点相同不意味着两次访视之间状态始终不变。S0→S3等跳跃不被强行拆成某种发病先后顺序。", "", "## 缺失、死亡与时间", "", "当波存在真实访视且三项疾病信息完整，才使用冻结的累计疾病史状态。未访视或任一疾病项缺失均记Not observed，不新增疾病LOCF。该图保留原清洗中累计诊断史的定义，并非每波原始自报是否患病的简单快照。", "", "HRS死亡时间只用于在缺少真实访视的波次放置死亡：死亡时间须不晚于该波真实访视时间的中位数，真实活体访视在所属波次优先。之后死亡作为吸收状态携带。波次中位参考日期及最早/最晚访视日期完整导出；这是一种将记录死亡放入离散波次图的规则，不是精确统一日历时点的状态占有率。", "", "ELSA已供应死亡记录仅延续到2012年，因此W4–W10主流图不设置死亡状态，既有早期死亡也留在Not observed并在源数据标记已知有限死亡，失访者可能包含后期死亡。不能比较两队列死亡流或用ELSA的N/O作为纯失访率。", "", "所有可用研究波次均保留。横轴是等间隔波次而非等间隔时间，ELSA W10实际访视跨2021–2023年，与W9之间存在较长间隔。", "", "## 完整源数据及验收", "", "已输出逐人逐波长表和固定队列宽表、去重端点模式及频数、每波观察/缺失/死亡/延迟入组分母、相邻波流量和全部数值QC。固定队列每波人数及每段流量守恒，实际疾病状态必须对应完整真实访视，ELSA没有由缺失创建的死亡。图形由R生成PDF、SVG、600-dpi TIFF及PNG；自动文本、碰撞与面板对齐检查另存QC。", "", "这些是未加权经验样本描述，无p值或置信区间；统计推断仍应来自此前冻结的多状态模型。")
writeLines(report, file.path(report_root, "empirical_state_flow_report_zh.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(qc_root, "sessionInfo.txt"))
cat("12b complete. Numerical checks:", sum(rbindlist(checks)$pass), "/", nrow(rbindlist(checks)), "\n")
print(summary_table)
