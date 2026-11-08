suppressPackageStartupMessages({
  library(data.table)
  library(haven)
  library(survey)
})

options(datatable.print.nrows = 50L, survey.lonely.psu = "adjust")

stage_root <- Sys.getenv("GA_STAGE_ROOT", unset = "private/stage")
input_root <- Sys.getenv(
  "GA_CLEAN_INPUT_ROOT",
  unset = file.path(stage_root, "output", "data_clean")
)
output_root <- Sys.getenv(
  "GA_ANALYSIS_OUTPUT_ROOT",
  unset = file.path(stage_root, "analysis_output")
)
derived_root <- file.path(output_root, "data_derived")
qc_root <- file.path(output_root, "qc")
dir.create(derived_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)

locf_down <- function(x) {
  if (!length(x)) return(x)
  idx <- which(!is.na(x))
  if (!length(idx)) return(x)
  first <- idx[1L]
  if (first < length(x)) {
    for (i in seq.int(first + 1L, length(x))) {
      if (is.na(x[i])) x[i] <- x[i - 1L]
    }
  }
  x
}

valid_positive_weight <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x) | x <= 0] <- NA_real_
  x
}

read_elsa_cross_sectional_weights <- function(root) {
  pieces <- list()
  for (w in 4:10) {
    path <- file.path(root, "elsa", sprintf("w%d.dta", w))
    if (!file.exists(path)) next
    meta <- read_dta(path, n_max = 0)
    nm <- names(meta)
    id_idx <- which(tolower(nm) %in% c("idauniq", "idauniqc"))[1L]
    wt_idx <- which(tolower(nm) == sprintf("w%dxwgt", w))[1L]
    if (is.na(id_idx) || is.na(wt_idx)) next
    selected <- read_dta(
      path,
      col_select = tidyselect::all_of(c(nm[id_idx], nm[wt_idx]))
    )
    selected <- as.data.table(selected)
    setnames(selected, names(selected), c("id", "sampling_weight"))
    selected[, `:=`(
      id = as.numeric(id),
      wave = as.integer(w),
      sampling_weight = as.numeric(sampling_weight)
    )]
    pieces[[length(pieces) + 1L]] <- selected[, .(id, wave, sampling_weight)]
  }
  if (!length(pieces)) {
    return(data.table(id = numeric(), wave = integer(), sampling_weight = numeric()))
  }
  unique(rbindlist(pieces, use.names = TRUE), by = c("id", "wave"))
}

build_exposure_panel <- function(dt, cohort_name) {
  x <- as.data.table(copy(dt))
  setorder(x, id, event_time, dead, wave, na.last = TRUE)

  x[, `:=`(
    insomnia_score_primary = as.integer(insomnia_score),
    insomnia_onset_primary = as.integer(insomnia_onset),
    insomnia_maintenance_primary = as.integer(insomnia_maintenance),
    insomnia_nonrestorative_primary = as.integer(insomnia_nonrestorative),
    insomnia_assessment_complete = as.integer(!is.na(insomnia_score_observed))
  )]
  x[, insomnia_burden_cat := fifelse(
    is.na(insomnia_score_primary), NA_character_,
    fifelse(insomnia_score_primary == 0L, "0 symptoms",
            fifelse(insomnia_score_primary == 1L, "1 symptom", ">=2 symptoms"))
  )]
  x[, insomnia_burden_code := fifelse(
    is.na(insomnia_score_primary), NA_integer_,
    pmin(insomnia_score_primary, 2L)
  )]
  x[, insomnia_high_primary := fifelse(
    is.na(insomnia_score_primary), NA_integer_,
    as.integer(insomnia_score_primary >= 2L)
  )]

  x[, insomnia_assessment_number := cumsum(!is.na(insomnia_score_observed)), by = id]
  x[, insomnia_score_cumavg := {
    valid <- !is.na(insomnia_score_observed)
    n_valid <- cumsum(valid)
    running_sum <- cumsum(fcoalesce(as.numeric(insomnia_score_observed), 0))
    fifelse(n_valid > 0L, running_sum / n_valid, NA_real_)
  }, by = id]
  x[, insomnia_score_lag1_row := shift(insomnia_score_primary), by = id]

  x[, `:=`(
    insomnia_score_prev_valid = {
      out <- rep(NA_real_, .N)
      idx <- which(!is.na(insomnia_score_observed))
      if (length(idx) >= 2L) out[idx[-1L]] <- insomnia_score_observed[idx[-length(idx)]]
      locf_down(out)
    },
    insomnia_prev_valid_wave = {
      out <- rep(NA_integer_, .N)
      idx <- which(!is.na(insomnia_score_observed))
      if (length(idx) >= 2L) out[idx[-1L]] <- wave[idx[-length(idx)]]
      locf_down(out)
    },
    insomnia_dynamic_pattern = {
      out <- rep(NA_character_, .N)
      idx <- which(!is.na(insomnia_score_observed))
      if (length(idx)) out[idx[1L]] <- "no prior assessment"
      if (length(idx) >= 2L) {
        current_high <- as.integer(insomnia_score_observed[idx] >= 2L)
        previous_high <- shift(current_high)
        out[idx[-1L]] <- fifelse(
          previous_high[-1L] == 0L & current_high[-1L] == 0L, "persistent low",
          fifelse(previous_high[-1L] == 0L & current_high[-1L] == 1L, "incident high",
                  fifelse(previous_high[-1L] == 1L & current_high[-1L] == 0L, "remission",
                          "persistent high"))
        )
      }
      locf_down(out)
    },
    insomnia_dynamic_source_wave = {
      out <- fifelse(!is.na(insomnia_score_observed), as.integer(wave), NA_integer_)
      locf_down(out)
    }
  ), by = id]
  x[, insomnia_assessment_gap_waves := fifelse(
    !is.na(insomnia_prev_valid_wave) & !is.na(insomnia_dynamic_source_wave),
    insomnia_dynamic_source_wave - insomnia_prev_valid_wave,
    NA_integer_
  )]
  x[, insomnia_dynamic_code := fcase(
    insomnia_dynamic_pattern == "persistent low", 0L,
    insomnia_dynamic_pattern == "incident high", 1L,
    insomnia_dynamic_pattern == "remission", 2L,
    insomnia_dynamic_pattern == "persistent high", 3L,
    default = NA_integer_
  )]

  if (cohort_name == "HRS") {
    x[, `:=`(
      insomnia_early_waking_observed = as.integer(insomnia_early_waking),
      insomnia_early_waking_current = locf_down(as.integer(insomnia_early_waking)),
      insomnia_early_waking_source_wave = locf_down(
        fifelse(!is.na(insomnia_early_waking), as.integer(wave), NA_integer_)
      ),
      sleep_medication_observed = as.integer(sleep_medication),
      sleep_medication_current = locf_down(as.integer(sleep_medication)),
      sleep_medication_source_wave = locf_down(
        fifelse(!is.na(sleep_medication), as.integer(wave), NA_integer_)
      )
    ), by = id]
    x[, insomnia_score_extended4 := fifelse(
      !is.na(insomnia_score_primary) & !is.na(insomnia_early_waking_current),
      insomnia_score_primary + insomnia_early_waking_current,
      NA_integer_
    )]
  } else {
    x[, `:=`(
      insomnia_early_waking_observed = NA_integer_,
      insomnia_early_waking_current = NA_integer_,
      insomnia_early_waking_source_wave = NA_integer_,
      insomnia_score_extended4 = NA_integer_,
      sleep_medication_observed = NA_integer_,
      sleep_medication_current = NA_integer_,
      sleep_medication_source_wave = NA_integer_
    )]
  }

  x[, exposure_entry_time := {
    ok <- event_type == "visit" & state %in% 0:3 & !is.na(insomnia_score_primary)
    if (any(ok)) min(event_time[ok]) else NA_real_
  }, by = id]
  x[, n_state_records_from_exposure_entry := {
    entry <- exposure_entry_time[1L]
    if (is.na(entry)) 0L else sum(!is.na(state) & event_time >= entry)
  }, by = id]
  x[, analysis_sample := as.integer(
    !is.na(exposure_entry_time) & n_state_records_from_exposure_entry >= 2L
  )]
  x[, analysis_time_since_entry := fifelse(
    analysis_sample == 1L & event_time >= exposure_entry_time,
    event_time - exposure_entry_time,
    NA_real_
  )]
  x[, analysis_baseline := as.integer(
    analysis_sample == 1L & event_type == "visit" &
      abs(event_time - exposure_entry_time) < 1e-8
  )]
  x[]
}

format_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

weighted_mean_sd <- function(v, w) {
  ok <- !is.na(v) & !is.na(w) & is.finite(v) & is.finite(w) & w > 0
  if (!any(ok)) return("")
  v <- as.numeric(v[ok])
  w <- as.numeric(w[ok])
  m <- sum(w * v) / sum(w)
  s <- sqrt(sum(w * (v - m)^2) / sum(w))
  sprintf("%.1f (%.1f)", m, s)
}

weighted_n_pct <- function(v, w, level) {
  ok <- !is.na(v) & !is.na(w) & is.finite(w) & w > 0
  if (!any(ok)) return("")
  n <- sum(v[ok] == level, na.rm = TRUE)
  pct <- 100 * sum(w[ok] * (v[ok] == level), na.rm = TRUE) / sum(w[ok])
  sprintf("%s (%.1f%%)", format(n, big.mark = ",", scientific = FALSE), pct)
}

weighted_missing <- function(v, w) {
  ok_w <- !is.na(w) & is.finite(w) & w > 0
  if (!any(ok_w)) return("")
  n <- sum(is.na(v))
  pct <- 100 * sum(w[ok_w] * is.na(v[ok_w])) / sum(w[ok_w])
  sprintf("%s (%.1f%%)", format(n, big.mark = ",", scientific = FALSE), pct)
}

weighted_p <- function(v, g, w, type) {
  ok <- !is.na(v) & !is.na(g) & !is.na(w) & is.finite(w) & w > 0
  if (sum(ok) < 10L || uniqueN(g[ok]) < 2L || uniqueN(v[ok]) < 2L) return(NA_real_)
  analyzed_v <- if (type == "continuous") v[ok] else factor(v[ok])
  d <- data.frame(.v = analyzed_v, .g = factor(g[ok]), .w = w[ok])
  design <- tryCatch(
    svydesign(ids = ~1, weights = ~.w, data = d),
    error = function(e) NULL
  )
  if (is.null(design)) return(NA_real_)
  tryCatch({
    if (type == "continuous") {
      fit <- svyglm(.v ~ .g, design = design)
      as.numeric(regTermTest(fit, ~.g)$p)
    } else {
      suppressWarnings(as.numeric(
        svychisq(~.v + .g, design = design, statistic = "F")$p.value
      ))
    }
  }, error = function(e) NA_real_)
}

make_row <- function(section, variable, level, values, p_value = "") {
  data.table(
    section = section,
    variable = variable,
    level = level,
    Overall = values[[1L]],
    group_1 = values[[2L]],
    group_2 = values[[3L]],
    group_3 = values[[4L]],
    group_4 = if (length(values) >= 5L) values[[5L]] else NA_character_,
    `P value` = p_value
  )
}

build_baseline_table <- function(b, group_var, group_levels, group_headers, cohort_name) {
  stopifnot(length(group_levels) == length(group_headers))
  b <- copy(b)
  b[, table_group := get(group_var)]
  subsets <- c(list(b), lapply(group_levels, function(g) b[table_group == g]))
  subset_values <- function(fun) vapply(subsets, fun, character(1))

  out <- list()
  add <- function(row) out[[length(out) + 1L]] <<- row
  add(make_row(
    "Population", "Participants, No.", "",
    subset_values(function(z) format(nrow(z), big.mark = ",", scientific = FALSE)), ""
  ))
  add(make_row(
    "Population", "Positive sampling weight, No.", "",
    subset_values(function(z) format(sum(!is.na(z$analysis_weight)), big.mark = ",", scientific = FALSE)), ""
  ))

  specs <- list(
    list("Demographic", "Age, years", "", "continuous", "age", NA),
    list("Demographic", "Female", "", "binary", "female", 1L),
    list("Demographic", "Race/ethnicity", if (cohort_name == "HRS") "Non-Hispanic White" else "White", "categorical", "race_group", if (cohort_name == "HRS") "Non-Hispanic White" else "White"),
    list("Demographic", "Race/ethnicity", if (cohort_name == "HRS") "Non-Hispanic Black" else "Non-White", "categorical", "race_group", if (cohort_name == "HRS") "Non-Hispanic Black" else "Non-White"),
    list("Demographic", "Race/ethnicity", "Hispanic", "categorical", "race_group", "Hispanic", cohort_name == "HRS"),
    list("Demographic", "Race/ethnicity", "Other", "categorical", "race_group", "Other", cohort_name == "HRS"),
    list("Demographic", "Comparable education level", "Below upper secondary", "categorical", "education_level", 1L),
    list("Demographic", "Comparable education level", "Upper secondary/vocational", "categorical", "education_level", 2L),
    list("Demographic", "Comparable education level", "Tertiary", "categorical", "education_level", 3L),
    list("Demographic", "Married/partnered", "", "binary", "married_partnered", 1L),
    list("Lifestyle", "Current smoker", "", "binary", "current_smoker", 1L),
    list("Lifestyle", "Current drinker", "", "binary", "current_drinker", 1L),
    list("Lifestyle", "Physically inactive", "", "binary", "physically_inactive", 1L),
    list("Lifestyle", "Body mass index, kg/m2", "", "continuous", "bmi", NA),
    list("Lifestyle", "Obesity (BMI >=30 kg/m2)", "", "binary", "obesity", 1L),
    list("Health", "Diabetes", "", "binary", "diabetes", 1L),
    list("Health", "Chronic lung disease", "", "binary", "lung_disease", 1L),
    list("Health", "Frequent pain", "", "binary", "frequent_pain", 1L),
    list("Health", "CES-D summary score", "", "continuous", "depression_score", NA),
    list("Sleep", "Insomnia symptom score (0-3)", "", "continuous", "insomnia_score_primary", NA),
    list("Sleep", "Difficulty initiating sleep", "", "binary", "insomnia_onset_primary", 1L),
    list("Sleep", "Night waking", "", "binary", "insomnia_maintenance_primary", 1L),
    list("Sleep", "Nonrestorative sleep", "", "binary", "insomnia_nonrestorative_primary", 1L),
    list("Sleep", "Sleep medication use", "", "binary", "sleep_medication_current", 1L, cohort_name == "HRS"),
    list("Baseline state", "S0: No target disease", "", "categorical", "state", 0L),
    list("Baseline state", "S1: Hypertension only", "", "categorical", "state", 1L),
    list("Baseline state", "S2: One CVD without hypertension", "", "categorical", "state", 2L),
    list("Baseline state", "S3: Multimorbidity", "", "categorical", "state", 3L),
    list("Missing data", "Missing comparable education level", "", "missing", "education_level", NA),
    list("Missing data", "Missing BMI", "", "missing", "bmi", NA),
    list("Missing data", "Missing CES-D score", "", "missing", "depression_score", NA)
  )

  for (sp in specs) {
    include <- if (length(sp) >= 7L) isTRUE(sp[[7L]]) else TRUE
    if (!include) next
    section <- sp[[1L]]; variable <- sp[[2L]]; level_label <- sp[[3L]]
    type <- sp[[4L]]; var <- sp[[5L]]; level <- sp[[6L]]
    values <- if (type == "continuous") {
      subset_values(function(z) weighted_mean_sd(z[[var]], z$analysis_weight))
    } else if (type == "missing") {
      subset_values(function(z) weighted_missing(z[[var]], z$analysis_weight))
    } else {
      subset_values(function(z) weighted_n_pct(z[[var]], z$analysis_weight, level))
    }
    p_type <- if (type == "continuous") "continuous" else "categorical"
    p <- if (type == "missing") {
      weighted_p(as.integer(is.na(b[[var]])), b$table_group, b$analysis_weight, "categorical")
    } else {
      weighted_p(b[[var]], b$table_group, b$analysis_weight, p_type)
    }
    if (type == "categorical") {
      first_level <- !any(vapply(out, function(z) {
        nrow(z) && z$section[1L] == section && z$variable[1L] == variable
      }, logical(1)))
      p_text <- if (first_level) format_p(p) else ""
    } else {
      p_text <- format_p(p)
    }
    add(make_row(section, variable, level_label, values, p_text))
  }

  result <- rbindlist(out, use.names = TRUE, fill = TRUE)
  if (length(group_headers) == 3L) {
    setnames(result, c("group_1", "group_2", "group_3"), group_headers)
    result[, group_4 := NULL]
  } else {
    setnames(result, c("group_1", "group_2", "group_3", "group_4"), group_headers)
  }
  result[, cohort := cohort_name]
  setcolorder(result, c("cohort", "section", "variable", "level", "Overall", group_headers, "P value"))
  result[]
}

prepare_baseline <- function(panel, cohort_name) {
  b <- panel[analysis_baseline == 1L]
  setorder(b, id, event_time)
  b <- b[, .SD[1L], by = id]
  b[, female := fifelse(!is.na(sex_male), as.integer(sex_male == 0L), NA_integer_)]
  b[, married_partnered := fifelse(
    !is.na(marital_status), as.integer(marital_status %in% c(1L, 2L, 3L)), NA_integer_
  )]
  b[, obesity := fifelse(!is.na(bmi), as.integer(bmi >= 30), NA_integer_)]
  b[, education_level := fifelse(raeducl %in% 1:3, as.integer(raeducl), NA_integer_)]
  if (cohort_name == "HRS") {
    b[, race_group := fcase(
      hispanic == 1L, "Hispanic",
      hispanic == 0L & race_code == 1L, "Non-Hispanic White",
      hispanic == 0L & race_code == 2L, "Non-Hispanic Black",
      hispanic == 0L & race_code == 3L, "Other",
      default = NA_character_
    )]
    b[, physically_inactive := fifelse(
      !is.na(vigorous_activity) & !is.na(moderate_activity),
      as.integer(vigorous_activity == 0L & moderate_activity == 0L),
      NA_integer_
    )]
  } else {
    b[, race_group := fifelse(
      is.na(race_white), NA_character_, fifelse(race_white == 1L, "White", "Non-White")
    )]
    b[, physically_inactive := fifelse(
      !is.na(vigorous_activity) & !is.na(moderate_activity),
      as.integer(vigorous_activity == 5L & moderate_activity == 5L),
      NA_integer_
    )]
  }
  b[, analysis_weight := valid_positive_weight(sampling_weight)]
  b[]
}

hrs <- readRDS(file.path(input_root, "hrs_main_analysis.rds"))
elsa <- readRDS(file.path(input_root, "elsa_main_analysis.rds"))

elsa_weights <- read_elsa_cross_sectional_weights(stage_root)
elsa <- as.data.table(elsa)
if ("sampling_weight" %in% names(elsa)) elsa[, sampling_weight := NULL]
elsa <- merge(elsa, elsa_weights, by = c("id", "wave"), all.x = TRUE, sort = FALSE)

hrs_panel <- build_exposure_panel(hrs, "HRS")
elsa_panel <- build_exposure_panel(elsa, "ELSA")

hrs_baseline <- prepare_baseline(hrs_panel, "HRS")
elsa_baseline <- prepare_baseline(elsa_panel, "ELSA")

table_burden <- rbindlist(list(
  build_baseline_table(
    hrs_baseline, "insomnia_burden_code", c(0L, 1L, 2L),
    c("0 symptoms", "1 symptom", ">=2 symptoms"), "HRS"
  ),
  build_baseline_table(
    elsa_baseline, "insomnia_burden_code", c(0L, 1L, 2L),
    c("0 symptoms", "1 symptom", ">=2 symptoms"), "ELSA"
  )
), use.names = TRUE, fill = TRUE)

table_state <- rbindlist(list(
  build_baseline_table(
    hrs_baseline, "state", c(0L, 1L, 2L, 3L),
    c("S0", "S1", "S2", "S3"), "HRS"
  ),
  build_baseline_table(
    elsa_baseline, "state", c(0L, 1L, 2L, 3L),
    c("S0", "S1", "S2", "S3"), "ELSA"
  )
), use.names = TRUE, fill = TRUE)

exposure_dictionary <- data.table::rbindlist(list(
  data.table(
    hierarchy = c("Primary", "Primary", "Primary", "Primary", "Primary"),
    cohort = c("HRS + ELSA", "HRS + ELSA", "HRS + ELSA", "HRS + ELSA", "HRS + ELSA"),
    variable = c(
      "insomnia_onset_primary", "insomnia_maintenance_primary",
      "insomnia_nonrestorative_primary", "insomnia_score_primary", "insomnia_burden_code"
    ),
    label = c(
      "Difficulty initiating sleep", "Night waking", "Nonrestorative sleep",
      "Primary insomnia symptom score", "Categorical insomnia burden"
    ),
    construction = c(
      "HRS: most of the time; ELSA: >=3 times/week",
      "HRS: most of the time; ELSA: >=3 times/week",
      "HRS: rarely/never feels rested; ELSA: >=3 times/week",
      "Sum of three harmonized binary symptoms (0-3)",
      "0 symptoms / 1 symptom / >=2 symptoms; high burden is >=2"
    ),
    timing = c(rep("Latest valid assessment at interval start", 4L), "Derived from current primary score"),
    missingness = c(
      rep("Observed values retained; latest valid value carried forward", 3L),
      "Complete three-symptom score; latest complete score carried forward",
      "Missing if primary score is missing"
    )
  ),
  data.table(
    hierarchy = c("Secondary", "Secondary", "Secondary", "Secondary"),
    cohort = c("HRS + ELSA", "HRS", "HRS", "HRS + ELSA"),
    variable = c(
      "component symptoms", "insomnia_score_extended4", "sleep_medication_current",
      "insomnia_score_cumavg"
    ),
    label = c(
      "Three symptom domains analyzed separately", "Extended insomnia score (0-4)",
      "Sleep medication use", "Cumulative mean insomnia score"
    ),
    construction = c(
      "Binary onset, maintenance, and nonrestorative symptoms",
      "Primary 0-3 score plus early waking",
      "Current/latest valid medication indicator",
      "Running mean over complete insomnia assessments"
    ),
    timing = c(
      "Interval-start time-varying", "Interval-start time-varying",
      "Interval-start covariate/effect modifier", "Updated at each complete assessment and carried forward"
    ),
    missingness = c(
      "Uses harmonized current values", "Unavailable in ELSA",
      "Unavailable in ELSA", "Missing before first complete assessment"
    )
  ),
  data.table(
    hierarchy = c("Dynamic", "Dynamic", "Dynamic"),
    cohort = c("HRS + ELSA", "HRS + ELSA", "HRS + ELSA"),
    variable = c(
      "insomnia_dynamic_pattern", "insomnia_score_prev_valid", "insomnia_score_lag1_row"
    ),
    label = c(
      "Adjacent valid-assessment pattern", "Previous valid assessment score",
      "One-record lag of current primary score"
    ),
    construction = c(
      "Persistent low / incident high / remission / persistent high using >=2 cutoff",
      "Score from the immediately preceding complete assessment",
      "Primary score at the preceding ordered state record"
    ),
    timing = c(
      "Updated at complete assessments; latest pattern carried between assessments",
      "Updated only when a new complete assessment occurs", "Lagged one ordered panel record"
    ),
    missingness = c(
      "No prior assessment until two complete assessments exist",
      "Missing until the second complete assessment", "Missing at first panel record"
    )
  )
), use.names = TRUE, fill = TRUE)

qc <- rbindlist(lapply(list(HRS = hrs_panel, ELSA = elsa_panel), function(x) {
  cohort_name <- unique(x$cohort)[1L]
  b <- if (cohort_name == "HRS") hrs_baseline else elsa_baseline
  data.table(
    cohort = cohort_name,
    metric = c(
      "Exposure panel rows", "Exposure panel participants", "Final analysis participants",
      "Final analysis baseline score 0", "Final analysis baseline score 1",
      "Final analysis baseline score >=2", "Positive baseline sampling weight",
      "Participants with >=2 complete sleep assessments", "Rows with dynamic pattern after two assessments",
      "Baseline S0", "Baseline S1", "Baseline S2", "Baseline S3"
    ),
    value = c(
      nrow(x), uniqueN(x$id), nrow(b),
      sum(b$insomnia_burden_code == 0L), sum(b$insomnia_burden_code == 1L),
      sum(b$insomnia_burden_code == 2L), sum(!is.na(b$analysis_weight)),
      x[, uniqueN(id[insomnia_assessment_number >= 2L])],
      sum(x$insomnia_dynamic_pattern %chin% c("persistent low", "incident high", "remission", "persistent high"), na.rm = TRUE),
      sum(b$state == 0L), sum(b$state == 1L), sum(b$state == 2L), sum(b$state == 3L)
    )
  )
}), use.names = TRUE)

notes <- data.table(
  topic = c(
    "Analysis baseline", "Primary exposure timing", "ELSA design missingness", "Dynamic pattern",
    "Weights", "Comparable education", "HRS activity definition", "ELSA activity definition", "Depression covariate",
    "ELSA mortality limitation"
  ),
  note = c(
    paste(
      "First age-eligible visit with a valid disease state and complete primary insomnia score, with at least one subsequent state/death record.",
      "ELSA uses the documented broad heart phenotype and backcasts labelled Wave-10 component disputes; unidentified coexisting other-heart status remains missing with lower/upper sensitivity bounds."
    ),
    "The latest valid insomnia assessment is assigned at the interval start; no backward filling is used.",
    "ELSA primary sleep symptoms are scheduled at waves 4, 6, and 8. Design-missing waves use the latest valid value and remain flagged separately from item missingness.",
    "Patterns compare adjacent complete insomnia assessments using high burden >=2: persistent low, incident high, remission, persistent high.",
    "HRS uses the harmonized respondent cross-sectional weight; ELSA uses the wave-specific raw cross-sectional weight (waves 4-9). Nonpositive/missing weights are excluded from weighted summaries.",
    "Baseline tables use the three-level harmonized raeducl variable: below upper secondary, upper secondary/vocational, and tertiary. ELSA raedyrs_e is not treated as directly comparable schooling years.",
    "Physically inactive: both vigorous and moderate activity coded never (0).",
    "Physically inactive: both vigorous and moderate activity coded hardly ever/never (5).",
    "Current workbook reports the harmonized CES-D summary. Item-level removal of a sleep item must be completed before the final adjusted model if required by the protocol.",
    "The supplied Harmonized ELSA EOL A2 death-year file covers deaths only through 2012; later mortality must be supplemented before the final death-transition analysis."
  )
)

saveRDS(hrs_panel, file.path(derived_root, "hrs_exposure_panel.rds"), compress = "xz")
saveRDS(elsa_panel, file.path(derived_root, "elsa_exposure_panel.rds"), compress = "xz")
fwrite(hrs_panel, file.path(derived_root, "hrs_exposure_panel.csv.gz"), na = "")
fwrite(elsa_panel, file.path(derived_root, "elsa_exposure_panel.csv.gz"), na = "")
saveRDS(hrs_baseline, file.path(derived_root, "hrs_baseline.rds"), compress = "xz")
saveRDS(elsa_baseline, file.path(derived_root, "elsa_baseline.rds"), compress = "xz")
fwrite(hrs_baseline, file.path(derived_root, "hrs_baseline.csv.gz"), na = "")
fwrite(elsa_baseline, file.path(derived_root, "elsa_baseline.csv.gz"), na = "")
fwrite(table_burden, file.path(derived_root, "table1_by_insomnia_burden.csv"), na = "")
fwrite(table_state, file.path(derived_root, "table1_by_baseline_state.csv"), na = "")
fwrite(exposure_dictionary, file.path(derived_root, "exposure_dictionary.csv"), na = "")
fwrite(qc, file.path(qc_root, "exposure_baseline_qc.csv"), na = "")
fwrite(notes, file.path(qc_root, "exposure_baseline_notes.csv"), na = "")

cat("Exposure construction and baseline summaries completed.\n")
print(qc)
