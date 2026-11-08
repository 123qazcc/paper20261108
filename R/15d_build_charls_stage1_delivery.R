#!/usr/bin/env Rscript
# CHARLS first research step: four disease paths and two direct comparisons.
# This read-only consumer never fits a model or reads complete-case pilot effects.
# Run: Rscript R/15d_build_charls_stage1_delivery.R preflight or build
# Figure contract: analysis/charls_stage1_figure_contract.md
options(stringsAsFactors = FALSE, warn = 1)
arguments <- commandArgs(trailingOnly = TRUE)
mode <- if (length(arguments) == 1L) arguments else ""
if (!mode %in% c("preflight", "build")) stop("Use preflight or build", call. = FALSE)

required_packages <- c("data.table", "ggplot2", "svglite", "ragg", "jsonlite", "digest")
available_packages <- vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = TRUE)
if (!all(available_packages)) stop("Missing existing R dependencies: ",
  paste(required_packages[!available_packages], collapse = ", "), "; no installation attempted.")
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
setDTthreads(1L)
out <- "outputs/charls_external_extension_20260913"
destination <- file.path(out, "stage1_delivery")
stopifnot(dir.exists(file.path(out, "model_qc")),
          file.exists("analysis/charls_stage1_figure_contract.md"))
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
qc <- file.path(out, "model_qc")
sha <- function(path) digest::digest(file = path, algo = "sha256")
utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
checks <- list()
source_hashes <- list()
figure_stems <- list()
add <- function(check, pass, detail = "", kind = "NUMERICAL") {
  checks[[length(checks) + 1L]] <<- data.table(check = check, pass = isTRUE(pass),
    kind = kind, detail = as.character(detail))
  isTRUE(pass)
}
guard <- function(check, pass, detail = "", kind = "NUMERICAL") {
  if (!add(check, pass, detail, kind)) stop(paste(check, detail, sep = ": "), call. = FALSE)
}
safe_json <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}
write_json <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null")
}
publish <- function(status, detail = "", signature = NULL) {
  result_checks <- rbindlist(checks, fill = TRUE)
  fwrite(result_checks, file.path(destination, if (mode == "preflight") "preflight_checks.csv" else "numerical_checks.csv"))
  status_path <- file.path(destination, "build_status.json")
  previous_status <- safe_json(status_path)
  # Readiness rechecks must not erase an existing build's state or provenance.
  if (mode == "preflight" && !is.null(previous_status) &&
      identical(previous_status$status, "DRAFT_REQUIRES_HUMAN_FIGURE_AND_SCIENTIFIC_REVIEW"))
    status_path <- file.path(destination, "preflight_status.json")
  write_json(list(status = status, mode = mode, generated_utc = utc(), detail = detail,
    source_hashes = source_hashes, signature = signature, figure_stems = figure_stems,
    checks_passed = sum(result_checks$pass), checks_total = nrow(result_checks),
    selected_backend = "r", automatic_final_approval = FALSE,
    next_stage_authorized = FALSE,
    human_figure_review = "PENDING", human_scientific_review = "PENDING",
    script_sha256 = sha("R/15d_build_charls_stage1_delivery.R"),
    figure_contract_sha256 = sha("analysis/charls_stage1_figure_contract.md")),
    status_path)
  message(status, if (nzchar(detail)) paste0(": ", detail) else "")
  invisible(status)
}
same_number <- function(a, b, tolerance = 2e-7) {
  if (length(a) != length(b) || !identical(dim(a), dim(b))) return(FALSE)
  a <- as.numeric(a); b <- as.numeric(b)
  if (!identical(is.na(a), is.na(b)) || !identical(is.infinite(a), is.infinite(b))) return(FALSE)
  infinite <- is.infinite(a)
  if (any(a[infinite] != b[infinite])) return(FALSE)
  finite <- is.finite(a) & is.finite(b)
  all(abs(a[finite] - b[finite]) <= tolerance * pmax(1, abs(a[finite]), abs(b[finite])))
}

bundle_path <- file.path(out, "data/charls_model_bundle.rds")
model_paths <- file.path(out, "models", sprintf("MI_%02d.rds", 1:20))
result_paths <- file.path(out, "results", c("all_edge_MI20.csv", "four_target_MI20.csv",
  "heart_stroke_contrasts_MI20.csv", "joint_tests_MI20.csv", "rubin_matrices_MI20.rds"))
review_path <- file.path(qc, "independent_pool_review.csv")
attestation_path <- file.path(destination, "numerical_review_attestation.json")
required_hash_paths <- c(bundle_path, model_paths, result_paths,
  "R/14b_fit_charls_disease_models.R", "R/14e_independently_verify_charls_pool.R",
  file.path(qc, "model_bundle_review.csv"), "analysis/charls_model_stage2_contract.md")

run <- function() {
  add("existing_R_dependencies", all(available_packages),
    paste(vapply(required_packages, function(p) paste(p, as.character(utils::packageVersion(p))), ""), collapse = "; "), "DEPENDENCY")
  add("msm_version_metadata_available", requireNamespace("msm", quietly = TRUE),
      "Existing version metadata is used only to reconstruct the frozen signature.", "DEPENDENCY")
  # Invalidate any previous successful build before touching source/effect data.
  if (mode == "build") publish("BUILD_RUNNING", "Rechecking all upstream eligibility and source hashes; previous figures are not newly approved.")
  primary <- safe_json(file.path(qc, "background_runner_status.json"))
  review <- safe_json(file.path(qc, "background_review_status.json"))
  readiness <- c(
    add("primary_is_pooled", !is.null(primary) && identical(primary$status, "POOLED_AWAITING_REVIEW_AND_FIGURES"),
      if (is.null(primary)) "Missing/unreadable runner status" else primary$status, "READINESS"),
    add("background_numerical_review_passed", !is.null(review) &&
      identical(review$independent_pool_review, "NUMERICAL_CHECKS_PASSED"),
      if (is.null(review)) "Missing/unreadable review status" else review$independent_pool_review, "READINESS"),
    add("all_20_models_and_pooled_outputs_exist", all(file.exists(required_hash_paths)),
      paste(required_hash_paths[!file.exists(required_hash_paths)], collapse = " | "), "READINESS"),
    add("independent_review_and_attestation_exist", all(file.exists(c(review_path, attestation_path))),
      "No effect data are read until readiness and hash-attestation gates pass.", "READINESS"))
  if (!all(readiness)) {
    publish("NOT_READY", "Main MI20 results are not yet eligible; no pilot effects, medical estimates or figures were read or generated.")
    return(0L)
  }

  attestation <- safe_json(attestation_path)
  guard("attestation_status_schema", !is.null(attestation) &&
    identical(attestation$status, "NUMERICAL_CHECKS_PASSED") &&
    all(c("source_hashes", "signature", "review_csv_sha256", "review_script_sha256", "review_finished_utc") %in% names(attestation)),
    kind = "GATE")
  attested_hashes <- unlist(attestation$source_hashes, use.names = TRUE)
  guard("attestation_contains_all_required_sources", length(attested_hashes) > 0L &&
    !anyDuplicated(names(attested_hashes)) && all(required_hash_paths %in% names(attested_hashes)), kind = "GATE")
  guard("attested_paths_exist", all(file.exists(names(attested_hashes))), kind = "GATE")
  actual_hashes <- vapply(names(attested_hashes), sha, "")
  guard("all_attested_source_hashes_match", identical(unname(actual_hashes), unname(attested_hashes)), kind = "GATE")
  source_hashes <<- as.list(actual_hashes)
  guard("independent_review_hash_matches", identical(sha(review_path), attestation$review_csv_sha256), kind = "GATE")
  guard("independent_review_script_hash_matches", identical(sha("R/14e_independently_verify_charls_pool.R"),
    attestation$review_script_sha256), kind = "GATE")
  attestation_hash <- sha(attestation_path)
  source_hashes[[review_path]] <<- sha(review_path)
  source_hashes[[attestation_path]] <<- attestation_hash
  source_hashes[["R/15d_build_charls_stage1_delivery.R"]] <<- sha("R/15d_build_charls_stage1_delivery.R")
  source_hashes[["analysis/charls_stage1_figure_contract.md"]] <<- sha("analysis/charls_stage1_figure_contract.md")
  independent <- fread(review_path)
  guard("independent_review_columns", all(c("check", "pass", "review_status") %in% names(independent)), kind = "GATE")
  guard("every_independent_check_passed", nrow(independent) > 250L && !anyNA(independent$pass) &&
    all(independent$pass) && all(independent$review_status == "COMPUTATIONAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_REVIEW"), kind = "GATE")
  signature <- digest::digest(list(bundle = sha(bundle_path),
    script = sha("R/14b_fit_charls_disease_models.R"), msm = as.character(utils::packageVersion("msm"))), algo = "sha256")
  guard("attested_signature_current", identical(signature, attestation$signature), kind = "GATE")
  # Preflight intentionally stops before opening any effect table or full model.
  if (mode == "preflight") {
    publish("PREFLIGHT_READY_FOR_BUILD", "Hash-bound independent review is available. No model effects opened and no figures generated.", signature)
    return(0L)
  }

  bundle <- readRDS(bundle_path)
  guard("bundle_approved_with_twenty_imputations", identical(bundle$status, "APPROVED_FOR_PILOT_AND_MI") &&
    all(bundle$preparation_checks$pass) && length(bundle$imputations) == 20L)
  guard("signed_preparation_inputs_current", all(file.exists(names(bundle$input_hashes))) &&
    identical(unname(vapply(names(bundle$input_hashes), sha, "")), unname(bundle$input_hashes)))
  bundle_review <- fread(file.path(qc, "model_bundle_review.csv"))
  guard("bundle_independent_review_passed", nrow(bundle_review) > 0L && all(bundle_review$pass))
  edges <- as.data.table(bundle$edges)
  targets <- match(c("0-2", "0-4", "1-3", "1-5"), edges$transition)
  guard("nine_state_twenty_edge_four_target_schema", nrow(edges) == 20L && !anyNA(targets) &&
    !anyDuplicated(targets) && identical(as.integer(edges$exposure_group[targets]), 1:4))
  n_people <- nrow(bundle$baseline); n_records <- nrow(bundle$panel)
  n_deaths <- sum(bundle$panel$state == 9L)
  for (i in 1:20) {
    saved <- readRDS(model_paths[i])
    guard(sprintf("MI_%02d_accepted_identity_signature", i), identical(saved$status, "ACCEPTED") &&
      identical(saved$signature, signature) && identical(as.integer(saved$imputation), i) &&
      identical(saved$model, sprintf("MI_%02d", i)) && saved$n_people == n_people && saved$n_records == n_records)
    rm(saved); invisible(gc(FALSE))
  }
  pooled <- readRDS(result_paths[5L])
  guard("pooled_signature_twenty_imputations", identical(pooled$signature, signature) &&
    identical(as.integer(pooled$m), 20L) && identical(dim(pooled$Q), c(20L, 20L)))
  for (name in c("Q", "qbar", "U", "B", "T")) guard(paste0("finite_Rubin_", name), all(is.finite(pooled[[name]])))
  guard("Rubin_components_consistent", same_number(pooled$qbar, rowMeans(pooled$Q)) &&
    same_number(pooled$B, stats::cov(t(pooled$Q))) && same_number(pooled$T, pooled$U + 1.05 * pooled$B))
  C <- matrix(0, 2L, nrow(edges))
  C[1L, targets[1:2]] <- c(1, -1); C[2L, targets[3:4]] <- c(1, -1)
  guard("direct_comparison_matrix_and_targets", identical(as.integer(pooled$targets), targets) && same_number(pooled$C, C, 0))
  scalar <- function(L) {
    effect <- as.numeric(L %*% pooled$qbar)
    within <- diag(L %*% pooled$U %*% t(L)); between <- diag(L %*% pooled$B %*% t(L))
    total <- diag(L %*% pooled$T %*% t(L))
    guard("positive_scalar_variances", all(within > 0) && all(total > 0) && all(between >= -1e-12))
    df <- rep(Inf, length(effect)); nonzero <- between > 1e-15
    r <- 1.05 * between[nonzero] / within[nonzero]
    df[nonzero] <- 19 * (1 + 1 / r)^2
    se <- sqrt(total); half_width <- stats::qt(.975, df) * se
    data.table(log_effect = effect, se = se, df = df, HR = exp(effect),
      lower = exp(effect - half_width), upper = exp(effect + half_width),
      p = 2 * stats::pt(abs(effect / se), df, lower.tail = FALSE))
  }
  joint <- function(L) {
    V <- L %*% pooled$T %*% t(L); V <- (V + t(V)) / 2
    eigenvalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    guard("joint_covariance_positive_definite", min(eigenvalues) > 1e-12 * max(eigenvalues) && qr(V)$rank == nrow(L))
    standardized <- forwardsolve(t(chol(V)), as.numeric(L %*% pooled$qbar))
    statistic <- sum(standardized^2)
    data.table(chi_square = statistic, df = nrow(L), p = stats::pchisq(statistic, nrow(L), lower.tail = FALSE))
  }
  all_edges <- fread(result_paths[1L]); target_table <- fread(result_paths[2L])
  contrasts <- fread(result_paths[3L]); joints <- fread(result_paths[4L])
  guard("all_edges_complete_and_ordered", nrow(all_edges) == 20L && identical(all_edges$transition, edges$transition))
  guard("all_four_targets_retained_and_ordered", nrow(target_table) == 4L &&
    identical(target_table$transition, c("0-2", "0-4", "1-3", "1-5")))
  guard("two_direct_comparisons_defined", nrow(contrasts) == 2L &&
    identical(contrasts$origin, c("No disease", "HTN only")) && all(contrasts$contrast == "Heart minus stroke"))
  guard("two_joint_tests_defined", nrow(joints) == 2L && identical(joints$test,
    c("Four target coefficients = 0", "Two heart-stroke coefficient differences = 0")) &&
    all(joints$reference == "Asymptotic chi-square Wald on Rubin total covariance"))
  compare <- function(label, actual, expected) {
    guard(paste0(label, "_schema"), nrow(actual) == nrow(expected) && all(names(expected) %in% names(actual)))
    for (column in names(expected)) guard(paste0(label, "_", column), same_number(actual[[column]], expected[[column]]))
  }
  expected_all <- scalar(diag(nrow(edges)))
  expected_targets <- expected_all[targets]
  expected_targets[, p_BH_4 := p.adjust(p, "BH", n = 4L)]
  expected_contrasts <- scalar(C)
  expected_contrasts[, p_BH_2 := p.adjust(p, "BH", n = 2L)]
  expected_joints <- rbind(joint(diag(nrow(edges))[targets, , drop = FALSE]), joint(C))
  compare("all_edges", all_edges, expected_all)
  compare("four_targets", target_table, expected_targets)
  compare("direct_comparisons", contrasts, expected_contrasts)
  compare("joint_tests", joints, expected_joints)
  guard("BH4_and_global_gate_decisions", same_number(target_table$global_four_edge_p, rep(joints$p[1L], 4L)) &&
    identical(as.logical(target_table$passes_global_and_BH), joints$p[1L] < .05 & target_table$p_BH_4 < .05))
  guard("all_plot_values_positive_with_complete_CI", all(is.finite(c(target_table$HR, target_table$lower,
    target_table$upper, contrasts$HR, contrasts$lower, contrasts$upper))) &&
    all(target_table$lower > 0 & target_table$lower <= target_table$HR & target_table$upper >= target_table$HR) &&
    all(contrasts$lower > 0 & contrasts$lower <= contrasts$HR & contrasts$upper >= contrasts$HR))

  target_zh <- c("无目标疾病→单纯心脏病", "无目标疾病→单纯卒中", "单纯高血压→高血压＋心脏病", "单纯高血压→高血压＋卒中")
  origin_zh <- c("无目标疾病起点", "单纯高血压起点")
  target_table[, `:=`(path_zh = target_zh, estimand = "Conditional transition intensity ratio per one higher W1 restless-sleep frequency category", n_people = n_people,
    n_records = n_records, imputations = 20L, ci_type = "Pointwise 95% finite-MI t interval")]
  contrasts[, `:=`(origin_zh = origin_zh, estimand = "Heart transition HR divided by stroke transition HR at the same origin; not an incidence-rate ratio",
    n_people = n_people, n_records = n_records, imputations = 20L, ci_type = "Pointwise 95% finite-MI t interval")]
  joints[, test_zh := c("四条目标睡眠系数同时为零", "两个起点的心脏病－卒中系数差同时为零")]
  all_edges[, interpretation := "Exploratory supplementary output; non-target coefficients include equality constraints"]
  fwrite(target_table, file.path(destination, "table1_four_target_paths.csv"), bom = TRUE)
  fwrite(contrasts, file.path(destination, "table2_heart_stroke_direct_comparisons.csv"), bom = TRUE)
  fwrite(joints, file.path(destination, "table3_joint_tests.csv"), bom = TRUE)
  fwrite(all_edges, file.path(destination, "supplement_all_edges_exploratory.csv"), bom = TRUE)
  sample_summary <- data.table(n_people = n_people, n_records = n_records,
    accepted_interval_observed_death_records = n_deaths, imputations = 20L,
    exposure = "W1 restless-sleep frequency category 0-3",
    covariates = paste(bundle$core_covariates, collapse = "; "),
    independent_unit = "participant", model = "9-state continuous-time panel Markov")
  fwrite(sample_summary, file.path(destination, "sample_and_model_summary.csv"), bom = TRUE)

  format_p <- function(p) vapply(p, function(x) if (x == 0) "低于双精度可表示范围" else formatC(x, digits = 3, format = "g"), "")
  format_hr <- function(e, l, u) sprintf("%.3f（%.3f–%.3f）", e, l, u)
  path_status <- ifelse(target_table$passes_global_and_BH,
    ifelse(target_table$HR > 1, "全局＋BH(4)门控后支持正向关联", "全局＋BH(4)门控后支持负向关联"),
    ifelse(target_table$p < .05, "仅有名义关联，未通过全局＋BH(4)门控", "当前未检出明确关联；不等于无效应"))
  contrast_status <- ifelse(contrasts$p_BH_2 < .05,
    ifelse(contrasts$HR > 1, "BH(2)后支持心脏病HR较大", "BH(2)后支持心脏病HR较小"),
    "BH(2)后未检出明确差异；不等于等效")
  table_lines <- c("# CHARLS 第一研究环节结果表（待人工复核草稿）", "",
    "## 表 1｜四条目标路径", "", "| 路径 | HR（逐项95% CI） | 双侧P值 | BH(4) | 门控判断 |", "|---|---:|---:|---:|---|",
    paste0("| ", target_zh, " | ", format_hr(target_table$HR, target_table$lower, target_table$upper), " | ",
      format_p(target_table$p), " | ", format_p(target_table$p_BH_4), " | ", path_status, " |"), "",
    paste0("四目标全局渐近 Wald P＝", format_p(joints$p[1L]), "；支持条件为全局P<0.05且对应BH(4)<0.05。"), "",
    "## 表 2｜心脏病与卒中的直接比较", "", "| 起点 | 心脏病HR／卒中HR（逐项95% CI） | 双侧P值 | BH(2) | 判断 |", "|---|---:|---:|---:|---|",
    paste0("| ", origin_zh, " | ", format_hr(contrasts$HR, contrasts$lower, contrasts$upper), " | ",
      format_p(contrasts$p), " | ", format_p(contrasts$p_BH_2), " | ", contrast_status, " |"), "",
    "比值大于1表示同一睡眠频率增量对应的心脏病转移HR较大，不表示心脏病发生率高于卒中；单项比较用BH(2)，不额外用两差异联合检验作为门。", "",
    "## 表 3｜联合检验", "", "| 检验 | 渐近χ² | 自由度 | P |", "|---|---:|---:|---:|",
    paste0("| ", joints$test_zh, " | ", sprintf("%.3f", joints$chi_square), " | ", joints$df, " | ", format_p(joints$p), " |"), "",
    "联合检验采用Rubin总协方差上的渐近χ² Wald，不是有限插补D1检验。标量效应的95% CI和双侧P使用large-complete-sample Rubin有限插补t自由度，区间未经同时覆盖率校正。", "",
    paste0("n＝", n_people, "名参与者，", n_records, "条面板观测；20份插补不是20个独立样本。"),
    "睡眠暴露每升高一个0–3频率类别，不是每增加一天或一个失眠症状。无目标疾病仅指无高血压、心脏病和卒中。",
    "本表只来自完整且独立核验通过的MI20主分析；完整病例运行预检结果未被读取或替代主结果。")
  writeLines(table_lines, file.path(destination, "第一环节_结果表_草稿.md"), useBytes = TRUE)

  plot_dir <- file.path(destination, "figures")
  dir.create(plot_dir, showWarnings = FALSE)
  forest <- function(data, labels, title, x_title) {
    z <- copy(data); z[, y := rev(seq_len(.N))]
    limits_log <- range(log(c(z$lower, z$upper, 1)))
    padding <- max(diff(limits_log) * .13, .06)
    limits <- exp(limits_log + c(-padding, padding))
    breaks <- exp(pretty(log(limits), n = 4L))
    breaks <- sort(unique(c(1, breaks[breaks >= limits[1L] & breaks <= limits[2L]])))
    ggplot(z, aes(x = HR, y = y)) +
      geom_vline(xintercept = 1, linewidth = .35, linetype = "dashed", colour = "#777777") +
      geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .13,
        linewidth = .6, colour = "#246A99") +
      geom_point(shape = 15, size = 2.2, colour = "#246A99") +
      scale_x_log10(limits = limits, breaks = breaks, labels = function(x) formatC(x, digits = 3, format = "g")) +
      scale_y_continuous(breaks = z$y, labels = labels, limits = c(.4, nrow(z) + .6), expand = expansion(mult = 0)) +
      labs(title = title, subtitle = sprintf("CHARLS; n = %s participants; 20 imputations", format(n_people, big.mark = ",")),
        x = x_title, y = NULL,
        caption = "Points and bars: pooled estimate and pointwise 95% CI (finite-MI t).\nExact P values and multiplicity adjustments: companion tables.") +
      theme_classic(base_size = 8, base_family = "sans") +
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_line(linewidth = .35), axis.ticks.x = element_line(linewidth = .35),
        axis.text = element_text(size = 8, colour = "#252525"),
        axis.text.y = element_text(margin = margin(r = 9)),
        axis.title.x = element_text(size = 8, margin = margin(t = 8)),
        plot.title = element_text(size = 10, face = "bold", margin = margin(b = 6)),
        plot.subtitle = element_text(size = 8, margin = margin(b = 8)),
        plot.caption = element_text(size = 7, hjust = 0, lineheight = 1.15, margin = margin(t = 10)),
        plot.title.position = "plot", plot.caption.position = "plot",
        plot.margin = margin(10, 12, 10, 10), legend.position = "none")
  }
  export <- function(plot, stem, height_mm) {
    width_mm = 183
    width_in <- width_mm / 25.4; height_in <- height_mm / 25.4
    grDevices::pdf(paste0(stem, ".pdf"), width = width_in, height = height_in,
      family = "Helvetica", useDingbats = FALSE, onefile = TRUE)
    device_open <- TRUE
    on.exit(if (device_open) grDevices::dev.off(), add = TRUE)
    # Initialise the page before grob layout; layout may initialise grid itself.
    # A second newpage after ggplotGrob would preserve an unwanted blank page.
    grid::grid.newpage()
    grob <- ggplotGrob(plot)
    grid::grid.draw(grob); grid::grid.force()
    cells <- grob$layout[grepl("^panel$", grob$layout$name), , drop = FALSE]
    guard("single_rendered_panel", nrow(cells) == 1L, kind = "FIGURE")
    z <- cells[1L, ]
    grid::seekViewport(paste0(z$name, ".", z$t, "-", z$l, "-", z$b, "-", z$r))
    lo <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"), valueOnly = TRUE)
    hi <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"), valueOnly = TRUE)
    rect <- c(lo$x, lo$y, hi$x, hi$y) * 72
    grid::upViewport(0)
    guard("rendered_panel_geometry_finite_positive", all(is.finite(rect)) && rect[3L] > rect[1L] && rect[4L] > rect[2L], kind = "FIGURE")
    grDevices::dev.off(); device_open <- FALSE
    write_json(list(schema_version = 1L, backend = "r-pdf-rendered-viewports",
      figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
      panels = list(list(id = "single", bbox_pt = rect, grid_id = "single-panel",
        row_start = 0L, row_stop = 1L, col_start = 0L, col_stop = 1L)), exemptions = list()),
      paste0(stem, ".alignment-layout.json"))
    write_json(list(verdict = "NOT APPLICABLE", reason = "One measured plot-area; no comparable multi-panel rectangles.",
      generated_utc = utc(), pdf_sha256 = sha(paste0(stem, ".pdf"))), paste0(stem, ".alignment.json"))
    svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
    device_open <- TRUE; print(plot); grDevices::dev.off(); device_open <- FALSE
    ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in, units = "in", res = 600)
    device_open <- TRUE; print(plot); grDevices::dev.off(); device_open <- FALSE
    ragg::agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in", res = 600,
      compression = "lzw")
    device_open <- TRUE; print(plot); grDevices::dev.off(); device_open <- FALSE
    guard("four_figure_exports_nonempty", all(file.info(paste0(stem, c(".pdf", ".svg", ".png", ".tiff")))$size > 0), kind = "FIGURE")
  }
  stem1 <- file.path(plot_dir, "figure1_four_target_paths")
  stem2 <- file.path(plot_dir, "figure2_heart_stroke_direct_comparisons")
  plot1 <- forest(target_table,
    c("No target disease to\nHeart disease only", "No target disease to\nStroke only",
      "HTN only to\nHTN + heart disease", "HTN only to\nHTN + stroke"),
    "Restless-sleep frequency and four disease transitions",
    "Transition HR per category higher W1 restless-sleep frequency")
  plot2 <- forest(contrasts, c("No target disease\nat origin", "HTN only\nat origin"),
    "Direct comparison of heart and stroke associations", "Heart transition HR / stroke transition HR")
  export(plot1, stem1, 108); export(plot2, stem2, 90)
  fwrite(target_table, paste0(stem1, ".source_data.csv"), bom = TRUE)
  fwrite(contrasts, paste0(stem2, ".source_data.csv"), bom = TRUE)
  figure_stems <<- list(stem1, stem2)
  common_caption <- c(paste0("独立单位为参与者，n＝", n_people, "；", n_records, "条重复面板观测，20份基线协变量插补。"),
    "9状态连续时间面板Markov模型；核心调整年龄、性别、W1原始教育、婚姻、城乡、非住房金融财富、BMI、吸烟、饮酒、糖尿病、慢性肺病和疼痛。混杂变量系数在所有允许边共享。",
    "方点为Rubin合并点估计，横线为逐项95%有限插补t置信区间；竖虚线为无关联/无差异参考值1。区间不是同时置信区间。",
    "暴露为过去一周睡眠不安频率0–3类别每提高一级，不是每天或失眠症状个数。HTN为高血压；无目标疾病指无高血压、心脏病和卒中。",
    "频率类别的线性对数效应是本阶段模型假设，尚未在本阶段检验非线性。蓝色仅用于显示效应，不编码显著性。源数据仅含可公开汇总估计，不含个人原始记录。")
  writeLines(c("# 图 1｜基线睡眠不安频率与四条目标疾病转移（待人工复核）", "",
    "四条目标路径全部显示，不按显著性筛选。HR为给定当前状态及调整协变量的条件转移强度比。",
    common_caption, paste0("标量检验为双侧有限插补t检验；四项采用BH(4)，并要求四目标全局渐近Wald P<0.05。全局P＝", format_p(joints$p[1L]), "。"),
    "精确P、BH(4)、门控判断和数值区间见同行源数据；联合检验用Rubin总协方差的渐近χ²，不是D1。Source data are provided as a Source Data file."),
    paste0(stem1, ".caption_zh.md"), useBytes = TRUE)
  writeLines(c("# 图 2｜同一起点心脏病与卒中睡眠关联的直接比较（待人工复核）", "",
    "横轴为心脏病转移HR除以卒中转移HR，用同模型完整协方差计算后以Rubin规则合并；不是两种疾病发生率之比。",
    common_caption, "比值大于1表示睡眠与心脏病转移的HR较大；小于1表示其较小。两个直接比较以双侧有限插补t检验并采用BH(2)，不额外施加联合检验门控。",
    paste0("两差异联合渐近Wald P＝", format_p(joints$p[2L]), "，作为整体差异检验单列；不是D1。单条关联显著与否不同不构成两效应有差异的证据，未检出差异也不证明等效。"),
    "精确P、BH(2)和数值区间见同行源数据。Source data are provided as a Source Data file."),
    paste0(stem2, ".caption_zh.md"), useBytes = TRUE)

  heart_rows <- c(1L, 3L)
  heart_supported <- sum(target_table$passes_global_and_BH[heart_rows])
  cross_cohort_note <- if (all(target_table$HR[heart_rows] > 1)) {
    if (heart_supported == 2L) "CHARLS两条心脏病路径的点估计方向与既有两库一致，且本库两条均通过既定全局＋BH(4)门控。这可作为不同睡眠测量下的方向一致外部拓展证据，不自动称三库成功复制。"
    else if (heart_supported == 1L) "CHARLS两条心脏病路径的点估计方向与既有两库一致，但仅一条通过本库既定全局＋BH(4)门控；另一条的方向一致不能代替充分的统计证据。不自动称三库成功复制。"
    else "CHARLS两条心脏病路径仅在点估计方向上与既有两库一致，均未通过本库既定全局＋BH(4)门控；当前不足以称为获得明确重复支持，也不证明无关联。"
  } else {
    "CHARLS至少一条心脏病路径的点估计未延续既有两库的正方向。应先审视区间精度、睡眠测量和疾病定义差异，并按既定顺序检查情绪症状、睡眠时长及稳健性；不据此修改纳排或追逐显著性。"
  }
  cross_cohort_difference_note <- if (any(contrasts$p_BH_2 < .05))
    "CHARLS至少一个直接比较通过BH(2)，属于本库新的病种差异线索；HRS/ELSA原有正式差异检验未获支持，因此不能称三库一致的心脏病特异性。"
    else "CHARLS两个直接比较均未通过BH(2)，当前仍不足以声称心脏病与卒中的关联强度不同；各库未检出差异并不证明疾病效应相等。"
  writeLines(c("# CHARLS 第一研究环节复盘（待人工科学及图形审核草稿）", "",
    "状态：DRAFT_REQUIRES_HUMAN_FIGURE_AND_SCIENTIFIC_REVIEW", "",
    "本环节是四条疾病路径及直接比较；与前期代码内部stage2命名不同。下面内容由真实MI20合并结果生成，不能作为自动批准的最终结论。", "",
    "## 本轮回答了什么", "",
    paste0("- ", target_zh, "：HR＝", format_hr(target_table$HR, target_table$lower, target_table$upper),
      "；双侧P＝", format_p(target_table$p), "，BH(4)＝", format_p(target_table$p_BH_4), "。", path_status, "。"), "",
    paste0("- ", origin_zh, "的心脏病HR／卒中HR＝", format_hr(contrasts$HR, contrasts$lower, contrasts$upper),
      "；双侧P＝", format_p(contrasts$p), "，BH(2)＝", format_p(contrasts$p_BH_2), "。", contrast_status, "。"), "",
    paste0("四目标联合χ²(", joints$df[1L], ")＝", sprintf("%.3f", joints$chi_square[1L]), "，P＝", format_p(joints$p[1L]),
      "；两差异联合χ²(", joints$df[2L], ")＝", sprintf("%.3f", joints$chi_square[2L]), "，P＝", format_p(joints$p[2L]), "。"), "",
    "这些是统计关联及关联强度比较，不能解释为睡眠造成疾病、某条路径是机制或具有治疗效果。目标路径逐项支持与两种疾病效应差异是不同问题。", "",
    "## 与HRS／ELSA既有结果的关系", "",
    "既有HRS／ELSA疾病拆分中，无目标疾病→单纯心脏病及单纯高血压→高血压＋心脏病两条路径均获得正向名义支持；卒中相关支持不一致。正式心脏病与卒中差异检验未获支持，不能由单条路径显著性不同推断心脏病特异性。",
    cross_cohort_note, cross_cohort_difference_note,
    "以上为既有结论与本轮真实结果的定性对照；未新增跨库交互检验、异质性检验或三库荟萃。", "",
    "## 推断范围与仍未完成的核验", "",
    paste0("- n＝", n_people, "名参与者、", n_records, "条面板记录；纳入", n_deaths, "条被接受的区间观察死亡记录。重复访谈和20份插补都不视为新增独立参与者；未作全国患病率或全国绝对风险推断。"),
    "- 9状态包含8种存活疾病组合和死亡，12条单病增加边及8条死亡边。病史阳性前瞻性保持；部分未知疾病以兼容状态集合处理。疾病访谈间多病跳跃由潜在路径积分，不人为指定先后或精确发病时刻。",
    "- 死亡采用已接受登记访谈月份月末的已死亡上界，obstype=1；不当成精确死亡日，也不宣称全员完整死亡链接。两条争议末次访谈仍需后续敏感性检查，当前保留既有审计口径，不因结果方向自行改判。",
    "- 优化、Hessian和独立数值重建通过不证明Markov/时间齐性、无未测混杂、失访可忽略或MICE与完整多状态似然严格同源。自由参数Hessian为已保存拟合的数值导数，本轮未重新独立数值微分原始似然。",
    "- 标量采用large-complete-sample Rubin有限插补t自由度，而非Barnard–Rubin。联合检验是Rubin总协方差的渐近卡方Wald，不是有限插补D1。CI为逐项区间；BH(4)+全局门与BH(2)两个检验族分别报告。",
    "- CHARLS单题睡眠不安频率不等同于HRS/ELSA三症状评分；这是获知两库结果后的外部拓展，不称预注册验证、不直接三库荟萃。心脏病为宽泛综合问项，不能据此解释为已分别测量的冠心病、心梗或心衰亚型。",
    "- 所有四条目标均保留。全20边输出（含四条目标）另存为探索性附表，非目标边受既定共享参数约束；不据其名义P筛选新主要结论。不重启为S0→S1寻找显著性的纳排重试。", "",
    "## 交付前必须人工复核", "",
    "1. 核对两图全部点和95% CI与源表一致、对数坐标与参考线无歧义；最终尺寸检查文字、误差线和裁切。",
    "2. 阅读源代码、PDF字号与碰撞审计；任何FAIL、NOT AUDITABLE或WARN须显式记录并修复/解释，不能由脚本自动批准。",
    "3. 将关联结论、病种直接比较和不确定性分别表述；尤其不以一显著一不显著声称心脏病强于卒中。", "",
    "## 下一步如何受本轮指导（尚未执行）", "",
    "下一环节仍按约定优先CESD-9／睡眠时长，不因本轮显著性改变纳排或删除卒中路径。其作用是检查情绪症状和睡眠时长能否解释部分关联，而非通过调参制造支持。应在匹配的分析样本和经复核的插补规格中比较调整前后。",
    "如果区间较宽，下一轮须优先报告精度和样本损失，不能把未显著当无关联；如果存在病种差异证据，也须在后续调整中检验稳健性，不能直接升格为机制。",
    "本阶段仅报告静态模型；CESD-9、睡眠时长、动态模式及机制不在本阶段估计。"),
    file.path(destination, "第一环节_结果复盘_草稿.md"), useBytes = TRUE)
  # A hash-bound attestation must remain current throughout table/figure build.
  guard("source_hashes_unchanged_through_build", identical(unname(vapply(names(attested_hashes), sha, "")), unname(attested_hashes)), kind = "GATE")
  guard("independent_review_unchanged_through_build", identical(sha(review_path), attestation$review_csv_sha256), kind = "GATE")
  guard("numerical_attestation_unchanged_through_build", identical(sha(attestation_path), attestation_hash), kind = "GATE")
  figure_inventory <- rbindlist(lapply(figure_stems, function(stem) {
    files <- paste0(stem, c(".pdf", ".svg", ".png", ".tiff", ".source_data.csv", ".caption_zh.md", ".alignment-layout.json", ".alignment.json"))
    data.table(path = files, sha256 = vapply(files, sha, ""), size_bytes = file.info(files)$size)
  }))
  fwrite(figure_inventory, file.path(destination, "figure_inventory.csv"))
  writeLines(c("# 第一环节绘图审核待办", "",
    "当前仅完成R原生导出与源数据数值核对，不宣称PDF文字、碰撞或人工科学审核通过。", "",
    "对每份最终 PDF 运行源代码验证、audit_pdf_text和audit_figure_collisions，保存全部JSON及exit code。", "",
    "源代码预检的EDITABLE-TEXT规则只识别Cairo PDF设备，因此可能对原生Helvetica PDF给出FAIL。必须保留该原始FAIL，不用无效代码欺骗检查；需通过最终PDF的真实可选文本、字形和几何实证另行人工复核。本构建不自动豁免该FAIL或宣布源预检全通过。", "",
    "单面板对齐为NOT APPLICABLE，实际R PDF绘图区坐标已保存。导出PNG须人工逐图view_image，在最终物理尺寸下检查；非R工具只用于PDF几何/文字审核，不重绘科学图。", "",
    "原始FAIL阻止交付批准；WARN须人工解释；NOT AUDITABLE不等于通过。最终人工审核不得由脚本自动更改为PASS。"),
    file.path(destination, "figure_review_pending.md"), useBytes = TRUE)
  publish("DRAFT_REQUIRES_HUMAN_FIGURE_AND_SCIENTIFIC_REVIEW",
    "MI20-only tables, two R forest plots and Chinese recap generated; no next research stage started.", signature)
  0L
}

exit_code <- tryCatch(run(), error = function(e) {
  add("delivery_execution", FALSE, conditionMessage(e), "ERROR")
  publish("BLOCKED_REQUIRES_REVIEW", conditionMessage(e))
  1L
})
quit(save = "no", status = exit_code)
