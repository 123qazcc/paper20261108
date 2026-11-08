#!/usr/bin/env Rscript
# Post-estimation only: full delta-method covariance + Rubin MI20 uncertainty.
suppressPackageStartupMessages({library(data.table); library(msm); library(ggplot2)})
for (pkg in c("expm", "numDeriv", "digest", "svglite", "ragg")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg)
}
version <- "2026-09-11-v2-frechet"
root <- normalizePath(Sys.getenv("GA_PROJECT_ROOT", "."), winslash = "/", mustWork = TRUE)
primary <- file.path(root, "outputs", "hrs_elsa_multistate_final_20260909")
out <- file.path(root, "outputs", "hrs_elsa_multistate_nextsteps_20260909")
resdir <- file.path(out, "results", "absolute_risk")
figdir <- file.path(out, "figures", "absolute_risk")
qcdir <- file.path(out, "qc", "absolute_risk")
repdir <- file.path(out, "reports", "absolute_risk")
cachedir <- file.path(resdir, "component_cache")
for (d in c(resdir, figdir, qcdir, repdir, cachedir, file.path(figdir, "source_data"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE)
scores <- 0:3
horizons <- c(2, 5, 10)
configs <- list(
  HRS = list(prefix = "hrs_5state_full", directory = "hrs_mi20_final", nstate = 5L),
  ELSA = list(prefix = "elsa_4state_w4_w10", directory = "elsa_w10tri_final", nstate = 4L)
)

edges_of <- function(fit) {
  a <- fit$qmodel$imatrix
  rbindlist(lapply(seq_len(nrow(a)), function(i) data.table(from = i, to = which(a[i, ] != 0))))
}
q_from <- function(theta, score, mean_score, edges, nstate) {
  k <- nrow(edges)
  q <- matrix(0, nstate, nstate)
  q[cbind(edges$from, edges$to)] <- exp(theta[seq_len(k)] + theta[k + seq_len(k)] * (score - mean_score))
  diag(q) <- -rowSums(q)
  q
}
metric_catalog <- function(nstate) {
  rbindlist(lapply(scores, function(s) rbindlist(lapply(horizons, function(t) {
    data.table(score = s, horizon = t,
      metric = c(paste0("endpoint_S", 0:(nstate - 1)), "ever_S3", paste0("years_S", 0:(nstate - 1)), "S3_free_years"),
      quantity = c(rep("endpoint probability", nstate), "S3 cumulative incidence", rep("restricted state years", nstate), "restricted S3-free years"),
      bound = c(rep(1, nstate + 1), rep(t, nstate + 1)))
  }))))
}
evaluate <- function(theta, mean_score, edges, nstate) {
  unlist(lapply(scores, function(s) {
    q <- q_from(theta, s, mean_score, edges, nstate)
    q_s3 <- q; q_s3[4L, ] <- 0
    augmented <- rbind(cbind(q, diag(nstate)), matrix(0, nstate, 2L * nstate))
    unlist(lapply(horizons, function(t) {
      p <- as.matrix(expm::expm(q * t))[1L, ]
      incidence <- as.matrix(expm::expm(q_s3 * t))[1L, 4L]
      years <- as.matrix(expm::expm(augmented * t))[1L, nstate + seq_len(nstate)]
      c(p, incidence, years, sum(years[1:3]))
    }), use.names = FALSE)
  }), use.names = FALSE)
}
frechet_jacobian <- function(theta, mean_score, edges, nstate) {
  k <- nrow(edges)
  do.call(rbind, lapply(scores, function(s) {
    q <- q_from(theta, s, mean_score, edges, nstate)
    q_s3 <- q; q_s3[4L, ] <- 0
    augmented <- rbind(cbind(q, diag(nstate)), matrix(0, nstate, 2L * nstate))
    do.call(rbind, lapply(horizons, function(t) {
      vapply(seq_len(2L * k), function(j) {
        edge_j <- (j - 1L) %% k + 1L
        derivative_q <- matrix(0, nstate, nstate)
        a <- edges$from[edge_j]; b <- edges$to[edge_j]
        rate_derivative <- q[a, b] * if (j <= k) 1 else (s - mean_score)
        derivative_q[a, b] <- rate_derivative; derivative_q[a, a] <- -rate_derivative
        derivative_s3 <- derivative_q; derivative_s3[4L, ] <- 0
        derivative_augmented <- matrix(0, 2L * nstate, 2L * nstate)
        derivative_augmented[seq_len(nstate), seq_len(nstate)] <- derivative_q
        p <- expm::expmFrechet(q * t, derivative_q * t, expm = FALSE)$Lexpm[1L, ]
        incidence <- expm::expmFrechet(q_s3 * t, derivative_s3 * t, expm = FALSE)$Lexpm[1L, 4L]
        years <- expm::expmFrechet(augmented * t, derivative_augmented * t, expm = FALSE)$Lexpm[1L, nstate + seq_len(nstate)]
        c(p, incidence, years, sum(years[1:3]))
      }, numeric(2L * nstate + 2L))
    }))
  }))
}
fixed_step_jacobian <- function(fun, theta, step) {
  vapply(seq_along(theta), function(j) {
    upper <- lower <- theta; upper[j] <- upper[j] + step; lower[j] <- lower[j] - step
    (fun(upper) - fun(lower)) / (2 * step)
  }, numeric(length(fun(theta))))
}
mean_profile <- function(fit, score) {
  profile <- as.list(fit$qcmodel$covmeans)
  profile$insomnia_score_primary <- score
  profile
}
extract_component <- function(path, cohort, config) {
  hash <- digest::digest(file = path, algo = "sha256")
  imp <- as.integer(sub(".*__imp([0-9]+)\\.rds$", "\\1", basename(path)))
  cachefile <- file.path(cachedir, sprintf("%s_imp%02d.rds", cohort, imp))
  signature <- list(version = version, sha256 = hash, scores = scores, horizons = horizons)
  if (file.exists(cachefile)) {
    old <- readRDS(cachefile)
    if (identical(old$signature, signature)) return(old)
  }
  wrapper <- readRDS(path); fit <- wrapper$fit
  stopifnot(inherits(fit, "msm"), fit$opt$convergence == 0, isTRUE(fit$foundse), wrapper$imputation == imp)
  edge <- edges_of(fit); k <- nrow(edge); ns <- config$nstate
  coef_matrices <- coef(fit)
  j <- match("insomnia_score_primary", attr(coef_matrices, "covlabels"))
  stopifnot(!is.na(j), nrow(fit$qmodel$imatrix) == ns)
  index <- c(seq_len(k), k + (j - 1L) * k + seq_len(k))
  theta <- as.numeric(fit$paramdata$params[index])
  sigma <- fit$covmat[index, index, drop = FALSE]
  stopifnot(all(is.finite(theta)), all(is.finite(sigma)), min(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values) > 0)
  coefficient_error <- max(abs(theta[k + seq_len(k)] - coef_matrices$insomnia_score_primary[cbind(edge$from, edge$to)]))
  mean_score <- unname(fit$qcmodel$covmeans["insomnia_score_primary"])
  q_error <- max(vapply(scores, function(s) max(abs(q_from(theta, s, mean_score, edge, ns) - as.matrix(msm::qmatrix.msm(fit, covariates = mean_profile(fit, s), ci = "none")))), numeric(1)))
  catalog <- metric_catalog(ns)
  fun <- function(z) evaluate(z, mean_score, edge, ns)
  point <- fun(theta)
  jac <- frechet_jacobian(theta, mean_score, edge, ns)
  # Fixed steps avoid cancellation from very small relative perturbations of
  # near-zero regression coefficients. Neither numerical grid is used for CIs.
  jac_check <- fixed_step_jacobian(fun, theta, 1e-4)
  jac_check_half <- fixed_step_jacobian(fun, theta, 5e-5)
  derivative_error <- max(abs(jac - jac_check), abs(jac - jac_check_half))
  within_raw <- jac %*% sigma %*% t(jac)
  within_raw <- (within_raw + t(within_raw)) / 2
  stopifnot(length(point) == nrow(catalog), all(point > 0), all(point < catalog$bound), all(diag(within_raw) >= 0), all(is.finite(within_raw)))
  catalog[, point := point]
  p_sum_error <- max(abs(catalog[grepl("^endpoint_", metric), .(s = sum(point)), by = .(score, horizon)]$s - 1))
  years_sum <- catalog[grepl("^years_S", metric), .(s = sum(point)), by = .(score, horizon)]
  years_sum_error <- max(abs(years_sum$s - years_sum$horizon))
  reference_error <- 0
  for (s in scores) for (t in horizons) {
    expected <- as.numeric(msm::pmatrix.msm(fit, t = t, covariates = mean_profile(fit, s), ci = "none")[1L, ])
    observed <- catalog[score == s & horizon == t & grepl("^endpoint_", metric), point]
    reference_error <- max(reference_error, max(abs(expected - observed)))
  }
  # The block-exponential integral is independently checked against msm once per fit.
  los_ref <- as.numeric(msm::totlos.msm(fit, start = 1L, tot = 10, covariates = mean_profile(fit, 2), ci = "none"))
  los_error <- max(abs(los_ref - catalog[score == 2 & horizon == 10 & grepl("^years_S", metric), point]))
  stopifnot(coefficient_error < 1e-9, q_error < 1e-9, reference_error < 1e-8, los_error < 1e-6, derivative_error < 1e-6, p_sum_error < 1e-8, years_sum_error < 1e-7)
  scale <- 1 / (point * (1 - point / catalog$bound))
  transformed <- qlogis(point / catalog$bound)
  within_logit <- within_raw * outer(scale, scale)
  audit <- data.table(cohort = cohort, imputation = imp, model_sha256 = hash, coefficient_error = coefficient_error, q_error = q_error,
    msm_probability_error = reference_error, msm_expected_years_error = los_error, derivative_grid_error = derivative_error,
    probability_sum_error = p_sum_error, state_years_sum_error = years_sum_error, covariance_min_eigen = min(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values), status = "PASS")
  compact <- list(signature = signature, cohort = cohort, imputation = imp, catalog = catalog,
    transformed = transformed, within_logit = within_logit, within_raw = within_raw, audit = audit)
  saveRDS(compact, cachefile, compress = "gzip")
  message("Derived uncertainty: ", cohort, " ", imp, "/20")
  compact
}
rubin_scalar <- function(estimates, variances) {
  m <- length(estimates); stopifnot(m == 20, all(is.finite(estimates)), all(is.finite(variances)), all(variances >= -1e-12))
  u <- mean(pmax(variances, 0)); b <- var(estimates); total <- u + (1 + 1 / m) * b
  r <- if (u > 0) (1 + 1 / m) * b / u else 0
  df <- if (r < 1e-14) Inf else (m - 1) * (1 + 1 / r)^2
  se <- sqrt(total); critical <- qt(0.975, df)
  data.table(estimate = mean(estimates), standard_error = se, conf_low = mean(estimates) - critical * se,
    conf_high = mean(estimates) + critical * se, within_variance = u, between_variance = b, total_variance = total,
    degrees_freedom = df, m = m)
}

all_components <- list(); marginal <- list(); differences <- list(); audits <- list()
for (cohort in names(configs)) {
  config <- configs[[cohort]]
  paths <- sort(list.files(file.path(primary, "models", config$directory), pattern = paste0("^", config$prefix, "__M4_MI20__imp[0-9]{2}\\.rds$"), full.names = TRUE))
  stopifnot(length(paths) == 20)
  components <- lapply(paths, extract_component, cohort = cohort, config = config)
  stopifnot(identical(vapply(components, `[[`, integer(1), "imputation"), 1:20))
  catalog <- copy(components[[1]]$catalog)
  for (z in components) stopifnot(identical(z$catalog[, !"point"], catalog[, !"point"]))
  for (i in seq_len(nrow(catalog))) {
    pooled <- rubin_scalar(vapply(components, function(z) z$transformed[i], numeric(1)), vapply(components, function(z) z$within_logit[i, i], numeric(1)))
    pooled[, `:=`(estimate = plogis(estimate) * catalog$bound[i], conf_low = plogis(conf_low) * catalog$bound[i], conf_high = plogis(conf_high) * catalog$bound[i],
      standard_error_scale = "logit(quantity / bound)", arithmetic_mean_point = mean(vapply(components, function(z) z$catalog$point[i], numeric(1))))]
    marginal[[length(marginal) + 1L]] <- cbind(data.table(cohort = cohort), catalog[i, !"point"], pooled)
  }
  for (i in which(catalog$score == 2)) {
    j <- which(catalog$score == 0 & catalog$horizon == catalog$horizon[i] & catalog$metric == catalog$metric[i])
    stopifnot(length(j) == 1)
    est <- vapply(components, function(z) z$catalog$point[i] - z$catalog$point[j], numeric(1))
    vari <- vapply(components, function(z) z$within_raw[i, i] + z$within_raw[j, j] - 2 * z$within_raw[i, j], numeric(1))
    differences[[length(differences) + 1L]] <- cbind(data.table(cohort = cohort, comparison = "Exactly 2 versus exactly 0 symptoms, maintained during projection"), catalog[i, .(horizon, metric, quantity)], rubin_scalar(est, vari))
  }
  all_components[[cohort]] <- rbindlist(lapply(components, function(z) cbind(data.table(cohort = cohort, imputation = z$imputation), z$catalog,
    data.table(logit_estimate = z$transformed, within_logit_variance = diag(z$within_logit), within_raw_variance = diag(z$within_raw)))))
  audits[[cohort]] <- rbindlist(lapply(components, `[[`, "audit"))
}
marginal <- rbindlist(marginal); differences <- rbindlist(differences); audits <- rbindlist(audits)
marginal[, `:=`(ci_definition = "Pointwise 95% confidence interval; within-MI delta method, logit-scale Rubin MI20 pooling", profile = "Nuisance covariates fixed at each fitted model's own means", mortality_scope = fifelse(cohort == "HRS", "Death represented", "Nonfatal model; death not represented"))]
differences[, `:=`(ci_definition = "Pointwise paired 95% confidence interval; within-MI full covariance and identity-scale Rubin MI20 pooling", profile = "Fixed scores; nuisance covariates at fitted means; not a causal effect")]
fwrite(marginal, file.path(resdir, "absolute_quantities_mi20_95ci.csv"))
fwrite(differences, file.path(resdir, "fixed_score_2_vs_0_paired_differences.csv"))
fwrite(rbindlist(all_components), file.path(resdir, "absolute_quantity_mi_components.csv.gz"))
fwrite(audits, file.path(qcdir, "parameter_and_prediction_validation.csv"))
frozen <- fread(file.path(primary, "results", "markov_transition_probability_matrices.csv"))
old <- frozen[model_id == "M4_MI20" & analysis_id %in% c("hrs_5state_full", "elsa_4state_w4_w10") & from_state == 0 & grepl("^insomnia_[0-3]$", profile), .(cohort, score = as.integer(sub("insomnia_", "", profile)), horizon = time_years, metric = paste0("endpoint_S", to_state), old_point = probability)]
stopifnot(nrow(old) == 108L, !anyDuplicated(old[, .(cohort, score, horizon, metric)]))
reproduction <- merge(marginal[grepl("^endpoint_", metric), .(cohort, score, horizon, metric, arithmetic_mean_point)], old, by = c("cohort", "score", "horizon", "metric"))
reproduction[, absolute_error := abs(arithmetic_mean_point - old_point)]
stopifnot(nrow(reproduction) == 108L, max(reproduction$absolute_error) < 1e-8)
fwrite(reproduction, file.path(qcdir, "frozen_probability_point_reproduction.csv"))

palette <- c(`0` = "#286DA8", `2` = "#CB6840")
common_theme <- theme_classic(base_family = "Helvetica", base_size = 7) + theme(
  plot.title = element_text(size = 8.5, face = "bold"), plot.subtitle = element_text(size = 6.5),
  axis.text = element_text(size = 6.5), axis.title = element_text(size = 7), strip.text = element_text(size = 7, face = "bold"),
  legend.position = "bottom", legend.text = element_text(size = 6.5), plot.caption = element_text(size = 6, hjust = 0, lineheight = 1.35),
  plot.margin = margin(5, 7, 5, 5))
export <- function(plot, stem, height) {
  stopifnot(height > 0)
  source(file.path(root, "R", "helpers", "rendered_panel_alignment.R"), local = environment())
  require_rendered_panel_alignment(plot, stem, 180 / 25.4, height / 25.4)
  grDevices::pdf(paste0(stem, ".pdf"), width = 180 / 25.4, height = height / 25.4, family = "Helvetica", useDingbats = FALSE)
  print(plot); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = 180 / 25.4, height = height / 25.4); print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = 180, height = height, units = "mm", res = 600, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = 180, height = height, units = "mm", res = 600); print(plot); dev.off()
}
plotdata <- marginal[score %in% c(0, 2) & metric %in% c("endpoint_S3", "ever_S3")]
plotdata[, label := factor(metric, levels = c("ever_S3", "endpoint_S3"), labels = c("Reached multimorbidity by the horizon", "In multimorbidity at the horizon"))]
plotdata[, score_label := factor(score)]
plotdata[, cohort := factor(cohort, levels = c("HRS", "ELSA"))]
p <- ggplot(plotdata, aes(horizon, estimate * 100, colour = score_label, group = score_label)) +
  geom_line(linewidth = .45) + geom_errorbar(aes(ymin = conf_low * 100, ymax = conf_high * 100), width = .22, linewidth = .4) + geom_point(size = 1.7) +
  facet_grid(label ~ cohort) + scale_x_continuous(breaks = horizons) + scale_colour_manual(values = palette, labels = c("Exactly 0 symptoms", "Exactly 2 symptoms"), name = NULL) +
  labs(title = "Absolute multimorbidity projections with full MI uncertainty", subtitle = "Starting in S0; fixed insomnia score and nuisance covariates at fitted means",
    x = "Projection horizon (years)", y = "Model-based probability (%)",
    caption = "Pointwise 95% CIs: delta method plus Rubin MI20 pooling. Fixed-score projections are not causal effects.\nHRS includes death; ELSA does not. In ELSA, S3 is absorbing, so both rows coincide.") + common_theme
export(p, file.path(figdir, "fig_absolute_multimorbidity_95ci"), 140)
fwrite(plotdata, file.path(figdir, "source_data", "fig_absolute_multimorbidity_95ci.csv"))
year_data <- marginal[score %in% c(0, 2) & metric == "S3_free_years"]
year_data[, `:=`(score_label = factor(score), cohort_label = fifelse(cohort == "HRS", "HRS: alive and free of multimorbidity", "ELSA: death is not represented"))]
year_data[, cohort_label := factor(cohort_label, levels = c("HRS: alive and free of multimorbidity", "ELSA: death is not represented"))]
p_year <- ggplot(year_data, aes(horizon, estimate, colour = score_label, group = score_label)) + geom_line(linewidth = .45) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = .2, linewidth = .4) + geom_point(size = 1.7) +
  facet_wrap(~cohort_label, nrow = 1) + scale_x_continuous(breaks = horizons) + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  scale_colour_manual(values = palette, labels = c("Exactly 0 symptoms", "Exactly 2 symptoms"), name = NULL) +
  labs(title = "Restricted expected time without cardiovascular multimorbidity", subtitle = "Starting in S0; sum of expected time in S0-S2 up to each horizon", x = "Projection horizon (years)", y = "Restricted expected years",
    caption = "Pointwise 95% CIs include within- and between-imputation uncertainty. These are not lifetime life-expectancy estimates.\nELSA's nonfatal state model cannot estimate years both alive and disease-free; cross-cohort levels are not comparable.") + common_theme
export(p_year, file.path(figdir, "fig_restricted_multimorbidity_free_years_95ci"), 102)
fwrite(year_data, file.path(figdir, "source_data", "fig_restricted_multimorbidity_free_years_95ci.csv"))

# Paired contrasts make small changes interpretable without comparing
# overlapping marginal error bars or dropping covariance between profiles.
paired_plot <- differences[metric %in% c("ever_S3", "S3_free_years")]
paired_plot[, `:=`(display_scale = fifelse(metric == "ever_S3", 100, 1),
  cohort = factor(cohort, levels = c("HRS", "ELSA")),
  metric_label = factor(metric, levels = c("ever_S3", "S3_free_years"),
    labels = c("S3 cumulative probability (percentage points)", "Restricted S3-free time (years)")))]
p_paired <- ggplot(paired_plot, aes(horizon, estimate * display_scale, colour = cohort, shape = cohort)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .35, colour = "#777777") +
  geom_errorbar(aes(ymin = conf_low * display_scale, ymax = conf_high * display_scale),
    position = position_dodge(width = .35), width = .18, linewidth = .4) +
  geom_point(position = position_dodge(width = .35), size = 1.9) +
  facet_wrap(~metric_label, scales = "free_y", nrow = 1) + scale_x_continuous(breaks = horizons) +
  scale_colour_manual(values = c(HRS = "#286DA8", ELSA = "#CB6840"), name = NULL) + scale_shape_discrete(name = NULL) +
  labs(title = "Paired absolute differences: fixed score 2 versus score 0", subtitle = "Within-model covariance retained; S0 origin and covariates at fitted means",
    x = "Projection horizon (years)", y = "Difference (fixed score 2 minus score 0)",
    caption = "Pointwise 95% CIs: identity-scale Rubin MI20 pooling; not multiplicity-adjusted. Fixed-score differences are not causal effects.\nHRS represents death; ELSA does not. ELSA's restricted time is not time alive and disease-free.") + common_theme
export(p_paired, file.path(figdir, "fig_paired_absolute_differences_95ci"), 105)
fwrite(paired_plot, file.path(figdir, "source_data", "fig_paired_absolute_differences_95ci.csv"))

ten <- differences[horizon == 10 & metric %in% c("ever_S3", "endpoint_S3", "S3_free_years")]
lines <- c("# 绝对风险和受限状态时间的不确定性", "", paste0("生成时间：", Sys.time()), "",
  "已复用HRS和ELSA各20个冻结主M4-MI模型，通过矩阵指数的Frechet导数，在每个插补中传播基线强度及失眠系数的完整协方差，再按Rubin规则合并插补内与插补间不确定性。未重新拟合主模型。",
  "", "## 估计量与边界", "",
  "从S0起始，在其他模型矩阵协变量保持各模型样本均值、失眠分数持续固定为0、1、2或3的条件下，估计2、5、10年结局状态概率、到达S3的累积概率和各状态受限期望年数。这些是恒定强度模型的条件性投影；固定2分不等于实际>=2分组，也不代表干预的因果效应。",
  "边际置信区间在logit(数量/上界)尺度上作delta-Rubin合并，以保持合法范围。2分相对0分的差值使用同一模型完整协方差及原始尺度配对合并。区间是点对点95%置信区间，不是同时置信带或个体预测区间；样本均值剖面作为固定条件，不含其抽样不确定性。",
  "HRS在S3时点的概率仅计当时仍活着的S3参与者；曾经到达S3的累积概率将S3在派生矩阵中设为吸收态，允许到达S3后死亡者继续计入。ELSA四状态不含死亡，S3本身吸收，因此两种概率相同，不能据此验证死亡或比较健康寿命。",
  "受限无S3时间仅积分至指定时限：HRS为活着且未进入S3的期望年数；ELSA未考虑死亡，不能称为活着且健康的年数。",
  "", "## 10年固定2分相对0分的配对差值", "")
for (i in seq_len(nrow(ten))) {
  z <- ten[i]; is_year <- z$metric == "S3_free_years"; scale <- if (is_year) 1 else 100
  lines <- c(lines, sprintf("- %s，%s：%.3f（95%% CI %.3f至%.3f）%s。", z$cohort, z$metric, z$estimate * scale, z$conf_low * scale, z$conf_high * scale, if (is_year) "年" else "个百分点"))
}
lines <- c(lines, "", "## 验证", "", "40/40模型参数映射、强度矩阵、概率矩阵、受限年数积分及导数步长交叉检查通过。108个既有概率点估计在算术平均尺度精确复现；新的有界区间表以logit池化的中心值呈现，另保留算术平均列以便核对。完整系数协方差、配对不确定性和每个插补的派生结果可复算。", "", "脚本：R/12a_absolute_risk_uncertainty.R。全部分数/时限/状态和配对差值在results/absolute_risk中；图件及源数据在figures/absolute_risk中。")
writeLines(enc2utf8(lines), file.path(repdir, "absolute_risk_uncertainty_report_zh.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(qcdir, "session_info.txt"))
message("Absolute-risk MI20 uncertainty completed")
