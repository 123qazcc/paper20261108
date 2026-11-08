#!/usr/bin/env Rscript

# Stage 1 cleaning for the HRS x ELSA insomnia/CVD multi-state study and
# the CHARLS x ELSA sleep-duration/BP mechanism module.
#
# The GlobalAgingR 2.0.0 harmonized cohorts are the longitudinal spine.
# Original ELSA core files supply sleep items that are absent from the package
# dataset. Source data are never modified.

options(stringsAsFactors = FALSE, width = 120)

required_packages <- c("GlobalAgingR", "data.table", "dplyr", "haven", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(haven)
  library(tidyr)
  library(GlobalAgingR)
})

stage_root <- Sys.getenv("GA_STAGE_ROOT", "private/stage")
output_root <- Sys.getenv("GA_OUTPUT_ROOT", file.path(stage_root, "output"))
data_output <- file.path(output_root, "data_clean")
qc_output <- file.path(output_root, "qc")
log_output <- file.path(output_root, "logs")
dir.create(data_output, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_output, recursive = TRUE, showWarnings = FALSE)
dir.create(log_output, recursive = TRUE, showWarnings = FALSE)

run_started <- Sys.time()
package_version <- as.character(utils::packageVersion("GlobalAgingR"))
if (package_version != "2.0.0") {
  warning("This pipeline was validated with GlobalAgingR 2.0.0; found ", package_version)
}

message("GlobalAgingR version: ", package_version)
message("Stage root: ", normalizePath(stage_root, winslash = "/", mustWork = TRUE))

raw_numeric <- function(x) {
  if (is.null(x)) return(numeric())
  y <- unclass(x)
  if (is.factor(y)) y <- as.character(y)
  suppressWarnings(as.numeric(y))
}

to_numeric_frame <- function(x, variables) {
  variables <- intersect(variables, names(x))
  out <- as.data.frame(x[, variables, drop = FALSE])
  out[] <- lapply(out, raw_numeric)
  out
}

valid_range <- function(x, lower, upper) {
  x <- raw_numeric(x)
  x[!is.finite(x) | x < lower | x > upper] <- NA_real_
  x
}

binary01 <- function(x) {
  x <- raw_numeric(x)
  out <- rep(NA_integer_, length(x))
  out[x == 0] <- 0L
  out[x == 1] <- 1L
  out
}

combine_binary_any <- function(...) {
  fields <- cbind(...)
  if (!ncol(fields)) return(integer())
  out <- rep(NA_integer_, nrow(fields))
  any_yes <- rowSums(fields == 1L, na.rm = TRUE) > 0L
  all_known_no <- rowSums(!is.na(fields)) == ncol(fields) &
    rowSums(fields == 0L, na.rm = TRUE) == ncol(fields)
  out[any_yes] <- 1L
  out[!any_yes & all_known_no] <- 0L
  out
}

yes1_no2 <- function(x) {
  x <- raw_numeric(x)
  out <- rep(NA_integer_, length(x))
  out[x == 1] <- 1L
  out[x == 2] <- 0L
  out
}

# Wave 10 uses feed-forward checks for diagnoses recorded at an earlier wave.
# A respondent may answer the current HEEverXX item "No" and then confirm at
# HENvrXX=1 that the carried-forward diagnosis is correct.  Such a record is
# still an ever-diagnosis.  Explicit rejection/misdiagnosis codes remain No;
# don't-know/refusal remains unknown.
ever_with_feedforward <- function(ever, prior_confirmation) {
  ever <- raw_numeric(ever)
  prior_confirmation <- raw_numeric(prior_confirmation)
  out <- rep(NA_integer_, length(ever))
  out[ever == 1 | prior_confirmation == 1] <- 1L
  # Only labelled questionnaire codes are interpreted.  Unlabelled raw codes
  # 85/86 occur in the deposit but are not documented as negative responses;
  # retaining them as unknown avoids silently inventing a disease-free status.
  rejected_or_not_routed <- prior_confirmation %in% c(-2, -1, 2, 3, 4)
  out[ever == 2 & rejected_or_not_routed] <- 0L
  out
}

is_explicit_dispute <- function(prior_confirmation) {
  raw_numeric(prior_confirmation) %in% c(2, 3, 4)
}

cum_diagnosis <- function(x) {
  x <- as.integer(x)
  out <- rep(NA_integer_, length(x))
  known <- FALSE
  ever <- FALSE
  for (i in seq_along(x)) {
    if (!is.na(x[i])) {
      known <- TRUE
      if (x[i] == 1L) ever <- TRUE
    }
    if (ever) out[i] <- 1L else if (known) out[i] <- 0L
  }
  out
}

derive_cvd_state <- function(hypertension, heart_disease, stroke, dead = 0L) {
  h <- as.integer(hypertension)
  c <- as.integer(heart_disease)
  s <- as.integer(stroke)
  d <- as.integer(dead)
  yes_count <- rowSums(cbind(h, c, s) == 1L, na.rm = TRUE)
  known_count <- rowSums(!is.na(cbind(h, c, s)))
  out <- rep(NA_integer_, length(h))
  out[d == 1L] <- 4L
  alive <- d != 1L
  out[alive & known_count == 3L & yes_count == 0L] <- 0L
  out[alive & h == 1L & c == 0L & s == 0L] <- 1L
  # S2 means exactly one of heart disease or stroke with confirmed absence of
  # hypertension, so all three component diagnoses must be observed.
  out[alive & known_count == 3L & h == 0L & yes_count == 1L & (c == 1L | s == 1L)] <- 2L
  out[alive & yes_count >= 2L] <- 3L
  out
}

add_cumulative_diseases_and_state <- function(dt) {
  setorder(dt, id, event_time, wave)
  dt[, hypertension := cum_diagnosis(hypertension_raw), by = id]
  dt[, heart_disease := cum_diagnosis(heart_disease_raw), by = id]
  dt[, stroke := cum_diagnosis(stroke_raw), by = id]
  dt[, diabetes := cum_diagnosis(diabetes_raw), by = id]
  dt[, lung_disease := cum_diagnosis(lung_disease_raw), by = id]
  dt[, state := derive_cvd_state(hypertension, heart_disease, stroke, dead)]
  dt[]
}

first_nonmissing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x)) x[1] else NA_real_
}

append_death_rows <- function(visits, deaths, study_end_year) {
  visits <- as.data.table(copy(visits))
  deaths <- as.data.table(copy(deaths))
  if (!nrow(deaths)) return(visits)
  deaths <- deaths[!is.na(death_year) & death_year <= study_end_year]
  deaths <- deaths[, .(
    death_year = first_nonmissing(death_year),
    death_source = first_nonmissing_character(death_source)
  ), by = id]

  last_visits <- visits[event_type == "visit" & !is.na(event_time), .(
    last_visit_time = max(event_time),
    last_visit_age = first_nonmissing(age[which.max(event_time)])
  ), by = id]
  deaths <- merge(deaths, last_visits, by = "id", all = FALSE)
  deaths[, death_time_midyear := death_year + 0.5]
  deaths[, death_time_adjusted := as.integer(death_time_midyear <= last_visit_time)]
  deaths[, event_time := pmax(death_time_midyear, last_visit_time + 1 / 365.25)]
  deaths[, age := ifelse(!is.na(last_visit_age), last_visit_age + event_time - last_visit_time, NA_real_)]
  deaths[!is.na(age) & (age < 18 | age > 110), age := NA_real_]

  # NA row selection preserves every source column's storage type.
  death_rows <- visits[rep(NA_integer_, nrow(deaths))]
  death_rows[, cohort := visits$cohort[1]]
  death_rows[, id := deaths$id]
  death_rows[, wave := 999L]
  death_rows[, interview_year := floor(deaths$event_time)]
  death_rows[, interview_month := NA_real_]
  death_rows[, event_time := deaths$event_time]
  death_rows[, event_type := "death"]
  death_rows[, dead := 1L]
  death_rows[, state := 4L]
  death_rows[, death_year := deaths$death_year]
  death_rows[, death_time_midyear := deaths$death_time_midyear]
  death_rows[, death_time_adjusted := deaths$death_time_adjusted]
  death_rows[, death_source := deaths$death_source]
  death_rows[, age := deaths$age]

  out <- rbindlist(list(visits, death_rows), use.names = TRUE, fill = TRUE)
  setorder(out, id, event_time, wave)
  out
}

first_nonmissing_character <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[1] else NA_character_
}

locf_down <- function(x) {
  if (!length(x)) return(x)
  last <- NA
  for (i in seq_along(x)) {
    if (!is.na(x[i])) last <- x[i] else if (!is.na(last)) x[i] <- last
  }
  x
}

make_main_analysis <- function(dt) {
  x <- as.data.table(copy(dt))
  setorder(x, id, event_time, wave)
  x[, eligible_visit := event_type == "visit" & !is.na(age) & age >= 50 & state %in% 0:3]
  entry <- x[eligible_visit == TRUE, .(entry_time = min(event_time)), by = id]
  x <- merge(x, entry, by = "id", all = FALSE)
  x <- x[event_time >= entry_time]
  x[, time_since_entry := event_time - entry_time]
  keep_ids <- x[, .(
    n_valid_states = sum(!is.na(state)),
    n_living_visits = sum(event_type == "visit" & state %in% 0:3),
    n_insomnia = sum(!is.na(insomnia_score))
  ), by = id][n_valid_states >= 2L & n_living_visits >= 1L & n_insomnia >= 1L, id]
  x <- x[id %in% keep_ids]
  x[, eligible_visit := NULL]
  setorder(x, id, event_time, wave)
  x[]
}

read_package_data <- function(item, object_name = item) {
  env <- new.env(parent = baseenv())
  utils::data(list = item, package = "GlobalAgingR", envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("GlobalAgingR data object not found: ", object_name, " (item ", item, ")")
  }
  get(object_name, envir = env, inherits = FALSE)
}

message("Cleaning HRS...")
hrs_package <- read_package_data("hrs")
hrs_filtered <- GlobalAgingR::CH_filter(hrs_package, wave %in% 8:15)
hrs_vars <- c(
  "hhidpn", "wave", "iwendm", "iwendy", "ragey_e", "ragender", "raracem", "rahispan",
  "raedyrs", "raeducl", "mstath", "wtresp", "proxy", "fallslp", "wakent", "wakeup",
  "rested", "rxslp", "hibpe", "hearte", "stroke", "diabe", "lunge", "painfr", "cesd",
  "smokev", "smoken", "drink", "vgactx", "mdactx", "ltactx", "bmi", "radyear",
  "radtimtdth", "shlt"
)
hrs <- as.data.table(to_numeric_frame(hrs_filtered, hrs_vars))
hrs[, cohort := "HRS"]
hrs[, id := hhidpn]
hrs[, wave := as.integer(wave)]
hrs[, interview_year := valid_range(iwendy, 1990, 2030)]
hrs[is.na(interview_year), interview_year := 1990 + 2 * wave]
hrs[, interview_month := valid_range(iwendm, 1, 12)]
hrs[, event_time := interview_year + ifelse(is.na(interview_month), 0.5, (interview_month - 0.5) / 12)]
hrs[, event_type := "visit"]
hrs[, dead := 0L]
hrs[, age := valid_range(ragey_e, 18, 110)]
hrs[, sex_male := binary01(ragender)]
hrs[, race_code := valid_range(raracem, 1, 3)]
hrs[, hispanic := binary01(rahispan)]
hrs[, education_years := valid_range(raedyrs, 0, 30)]
hrs[, marital_status := mstath]
hrs[, sampling_weight := valid_range(wtresp, 0, Inf)]
hrs[, proxy_interview := binary01(proxy)]
hrs[, self_rated_health := valid_range(shlt, 1, 5)]
hrs[, insomnia_onset := fifelse(fallslp %in% 1:3, as.integer(fallslp == 1), NA_integer_)]
hrs[, insomnia_maintenance := fifelse(wakent %in% 1:3, as.integer(wakent == 1), NA_integer_)]
hrs[, insomnia_early_waking := fifelse(wakeup %in% 1:3, as.integer(wakeup == 1), NA_integer_)]
hrs[, insomnia_nonrestorative := fifelse(rested %in% 1:3, as.integer(rested == 3), NA_integer_)]
hrs[, insomnia_n_observed := rowSums(!is.na(cbind(
  insomnia_onset, insomnia_maintenance, insomnia_nonrestorative
)))]
hrs[, insomnia_score_observed := fifelse(
  insomnia_n_observed == 3L,
  insomnia_onset + insomnia_maintenance + insomnia_nonrestorative,
  NA_integer_
)]
hrs[, insomnia_score := insomnia_score_observed]
hrs[, insomnia_source_wave := fifelse(!is.na(insomnia_score_observed), wave, NA_integer_)]
setorder(hrs, id, event_time, wave)
hrs[, c("insomnia_onset", "insomnia_maintenance", "insomnia_nonrestorative", "insomnia_score", "insomnia_source_wave") :=
      lapply(.SD, locf_down),
    by = id,
    .SDcols = c("insomnia_onset", "insomnia_maintenance", "insomnia_nonrestorative", "insomnia_score", "insomnia_source_wave")]
hrs[, insomnia_locf_reason := fifelse(
  !is.na(insomnia_score_observed), "observed",
  fifelse(is.na(insomnia_score), "missing", "latest_valid_locf")
)]
hrs[, insomnia_high := fifelse(!is.na(insomnia_score), as.integer(insomnia_score >= 2L), NA_integer_)]
hrs[, sleep_medication := binary01(rxslp)]
hrs[, hypertension_raw := binary01(hibpe)]
hrs[, heart_disease_raw := binary01(hearte)]
hrs[, stroke_raw := binary01(stroke)]
hrs[, diabetes_raw := binary01(diabe)]
hrs[, lung_disease_raw := binary01(lunge)]
hrs[, frequent_pain := binary01(painfr)]
hrs[, depression_score := valid_range(cesd, 0, 8)]
hrs[, ever_smoked := binary01(smokev)]
hrs[, current_smoker := binary01(smoken)]
hrs[, current_drinker := binary01(drink)]
hrs[, vigorous_activity := valid_range(vgactx, 0, 4)]
hrs[, moderate_activity := valid_range(mdactx, 0, 4)]
hrs[, light_activity := valid_range(ltactx, 0, 4)]
hrs[, bmi := valid_range(bmi, 10, 80)]
hrs[, death_year := valid_range(radyear, 1900, 2030)]
hrs[, death_time_midyear := NA_real_]
hrs[, death_time_adjusted := NA_integer_]
hrs[, death_source := NA_character_]
hrs <- add_cumulative_diseases_and_state(hrs)
hrs_deaths <- hrs[!is.na(death_year), .(
  death_year = first_nonmissing(death_year), death_source = "GlobalAgingR HRS radyear"
), by = id]
hrs <- append_death_rows(hrs, hrs_deaths, study_end_year = 2020)

message("Cleaning ELSA...")
elsa_package <- read_package_data("elsa")
elsa_filtered <- GlobalAgingR::CH_filter(elsa_package, wave %in% 4:9)
elsa_vars <- c(
  "idauniqc", "wave", "iwindm", "iwindy", "agey", "ragender", "raracem", "raeduc_e",
  "raedyrs_e", "raeducl", "mstath", "hibpe", "hearte", "stroke", "diabe", "lunge",
  "angine", "hrtatte", "conhrtfe", "hrtmre", "hrtrhme",
  "painfr", "cesd", "smokev", "smoken", "drink", "vgactx_e", "mdactx_e", "ltactx_e",
  "mbmi", "rxhibp", "systo", "diasto", "radyear", "shlt"
)
elsa <- as.data.table(to_numeric_frame(elsa_filtered, elsa_vars))
elsa[, cohort := "ELSA"]
elsa[, id := idauniqc]
elsa[, wave := as.integer(wave)]
elsa[, interview_year := valid_range(iwindy, 1990, 2030)]
elsa_year_map <- c(`4` = 2008, `5` = 2010, `6` = 2012, `7` = 2014, `8` = 2016, `9` = 2018)
elsa[is.na(interview_year), interview_year := unname(elsa_year_map[as.character(wave)])]
elsa[, interview_month := valid_range(iwindm, 1, 12)]
elsa[, event_time := interview_year + ifelse(is.na(interview_month), 0.5, (interview_month - 0.5) / 12)]
elsa[, event_type := "visit"]
elsa[, dead := 0L]
elsa[, age := valid_range(agey, 18, 110)]
elsa[, sex_male := binary01(ragender)]
elsa[, race_white := binary01(raracem)]
elsa[, education_years := valid_range(raedyrs_e, 0, 30)]
elsa[, education_code := raeduc_e]
elsa[, marital_status := mstath]
elsa[, self_rated_health := valid_range(shlt, 1, 5)]

# Read the exact Harmonized ELSA G3 source underlying GlobalAgingR.  The
# package remains the longitudinal spine, but its Wave-5 `angine` column is
# entirely empty and a handful of wave rows lose all disease fields during
# packaging.  The source-wide variables restore those values by id/wave and
# make the documented broad heart phenotype independently reproducible.
elsa_harmonized_file <- file.path(stage_root, "elsa_harmonized", "h_elsa_g3.dta")
if (!file.exists(elsa_harmonized_file)) {
  stop("Staged Harmonized ELSA G3 source missing: ", elsa_harmonized_file)
}
elsa_disease_bases <- c(
  "hibpe", "hearte", "stroke", "angine", "hrtatte", "conhrtfe", "hrtmre", "hrtrhme"
)
elsa_harmonized_vars <- c(
  "idauniqc",
  unlist(lapply(4:9, function(w) paste0("r", w, elsa_disease_bases)), use.names = FALSE)
)
elsa_harmonized_wide <- haven::read_dta(
  elsa_harmonized_file,
  col_select = dplyr::all_of(elsa_harmonized_vars)
)
elsa_harmonized_wide <- as.data.table(elsa_harmonized_wide)
elsa_harmonized_panel <- rbindlist(lapply(4:9, function(w) {
  out <- data.table(
    id = raw_numeric(elsa_harmonized_wide[["idauniqc"]]),
    wave = as.integer(w)
  )
  for (variable in elsa_disease_bases) {
    out[[paste0("source_", variable)]] <- binary01(
      elsa_harmonized_wide[[paste0("r", w, variable)]]
    )
  }
  out
}))
if (elsa_harmonized_panel[, anyDuplicated(paste(id, wave))]) {
  stop("Harmonized ELSA disease panel is not unique by id/wave")
}
elsa <- merge(elsa, elsa_harmonized_panel, by = c("id", "wave"), all.x = TRUE, sort = FALSE)

elsa_package_source_audit <- rbindlist(lapply(elsa_disease_bases, function(variable) {
  rbindlist(lapply(4:9, function(w) {
    package_value <- binary01(elsa[wave == w][[variable]])
    source_value <- elsa[wave == w][[paste0("source_", variable)]]
    same <- (is.na(package_value) & is.na(source_value)) |
      (!is.na(package_value) & !is.na(source_value) & package_value == source_value)
    data.table(
      wave = w,
      variable = variable,
      package_nonmissing = sum(!is.na(package_value)),
      source_nonmissing = sum(!is.na(source_value)),
      package_source_mismatch = sum(!same)
    )
  }))
}))
fwrite(
  elsa_package_source_audit,
  file.path(qc_output, "elsa_package_harmonized_disease_audit.csv"),
  na = ""
)

elsa[, hypertension_raw := fcoalesce(source_hibpe, binary01(hibpe))]
elsa[, stroke_raw := fcoalesce(source_stroke, binary01(stroke))]
elsa[, heart_aggregate_source_raw := fcoalesce(source_hearte, binary01(hearte))]
elsa[, heart_angina_raw := fcoalesce(source_angine, binary01(angine))]
elsa[, heart_mi_raw := fcoalesce(source_hrtatte, binary01(hrtatte))]
elsa[, heart_failure_raw := fcoalesce(source_conhrtfe, binary01(conhrtfe))]
elsa[, heart_murmur_raw := fcoalesce(source_hrtmre, binary01(hrtmre))]
elsa[, heart_rhythm_raw := fcoalesce(source_hrtrhme, binary01(hrtrhme))]

# Harmonized ELSA has no separate distributed variable for "other heart
# trouble".  Aggregate-positive rows with no named positive component identify
# the residual other/unclassified component.  When a named component is
# positive, that component already establishes the broad phenotype.
elsa_named_heart <- cbind(
  elsa$heart_angina_raw, elsa$heart_mi_raw, elsa$heart_failure_raw,
  elsa$heart_murmur_raw, elsa$heart_rhythm_raw
)
elsa_named_heart_yes <- rowSums(elsa_named_heart == 1L, na.rm = TRUE) > 0L
elsa[, heart_other_ambiguous_raw := as.integer(
  heart_aggregate_source_raw == 1L & elsa_named_heart_yes
)]
elsa[, heart_other_raw := fcase(
  heart_aggregate_source_raw == 1L & !elsa_named_heart_yes, 1L,
  heart_aggregate_source_raw == 0L, 0L,
  # When the broad aggregate is positive and a named component already
  # explains that positivity, coexisting "other heart trouble" is not
  # identifiable from the distributed Harmonized fields.  Preserve that as
  # unknown; coding it No would turn an identification limit into absence.
  default = NA_integer_
)]
elsa[, heart_disease_raw := combine_binary_any(
  heart_angina_raw, heart_mi_raw, heart_failure_raw,
  heart_murmur_raw, heart_rhythm_raw, heart_other_raw
)]
elsa[, diabetes_raw := binary01(diabe)]
elsa[, lung_disease_raw := binary01(lunge)]
elsa[, frequent_pain := binary01(painfr)]
elsa[, depression_score := valid_range(cesd, 0, 8)]
elsa[, ever_smoked := binary01(smokev)]
elsa[, current_smoker := binary01(smoken)]
elsa[, current_drinker := binary01(drink)]
elsa[, vigorous_activity := valid_range(vgactx_e, 2, 5)]
elsa[, moderate_activity := valid_range(mdactx_e, 2, 5)]
elsa[, light_activity := valid_range(ltactx_e, 2, 5)]
elsa[, bmi := valid_range(mbmi, 10, 80)]
elsa[, antihypertensive_medication := binary01(rxhibp)]
elsa[, systolic_bp_raw := systo]
elsa[, diastolic_bp_raw := diasto]
elsa[, systolic_bp := valid_range(systo, 70, 260)]
elsa[, diastolic_bp := valid_range(diasto, 40, 150)]
elsa[!is.na(systolic_bp) & !is.na(diastolic_bp) & systolic_bp <= diastolic_bp,
     c("systolic_bp", "diastolic_bp") := .(NA_real_, NA_real_)]
elsa[, pulse_pressure := systolic_bp - diastolic_bp]
elsa[, death_year := valid_range(radyear, 1900, 2030)]
elsa[, death_time_midyear := NA_real_]
elsa[, death_time_adjusted := NA_integer_]
elsa[, death_source := NA_character_]

elsa_sleep_files <- file.path(stage_root, "elsa", paste0("w", c(4, 6, 8), ".dta"))
if (!all(file.exists(elsa_sleep_files))) {
  stop("ELSA staged source file(s) missing: ", paste(elsa_sleep_files[!file.exists(elsa_sleep_files)], collapse = ", "))
}
elsa_sleep <- rbindlist(lapply(c(4L, 6L, 8L), function(w) {
  x <- haven::read_dta(
    file.path(stage_root, "elsa", paste0("w", w, ".dta")),
    col_select = dplyr::any_of(c("idauniq", "heslpa", "heslpb", "heslpd", "heslpe", "heslpf"))
  )
  x <- as.data.table(to_numeric_frame(x, names(x)))
  setnames(x, tolower(names(x)))
  x[, wave := w]
  x[, id := idauniq]
  x[, insomnia_onset_observed := fifelse(heslpa %in% 1:4, as.integer(heslpa == 4), NA_integer_)]
  x[, insomnia_maintenance_observed := fifelse(heslpb %in% 1:4, as.integer(heslpb == 4), NA_integer_)]
  x[, insomnia_nonrestorative_observed := fifelse(heslpd %in% 1:4, as.integer(heslpd == 4), NA_integer_)]
  x[, insomnia_n_observed := rowSums(!is.na(cbind(
    insomnia_onset_observed, insomnia_maintenance_observed, insomnia_nonrestorative_observed
  )))]
  x[, insomnia_score_observed := fifelse(
    insomnia_n_observed == 3L,
    insomnia_onset_observed + insomnia_maintenance_observed + insomnia_nonrestorative_observed,
    NA_integer_
  )]
  x[, sleep_hours_observed := valid_range(heslpe, 1, 18)]
  x[, sleep_quality_observed := fifelse(heslpf %in% 1:4, as.integer(heslpf), NA_integer_)]
  x[, .(
    id, wave, insomnia_onset_observed, insomnia_maintenance_observed,
    insomnia_nonrestorative_observed, insomnia_n_observed, insomnia_score_observed,
    sleep_hours_observed, sleep_quality_observed
  )]
}))
elsa <- merge(elsa, elsa_sleep, by = c("id", "wave"), all.x = TRUE)
elsa[, sleep_hours_centered7 := sleep_hours_observed - 7]

elsa_w10_file <- file.path(stage_root, "elsa", "w10.dta")
if (!file.exists(elsa_w10_file)) stop("ELSA Wave 10 staged source file missing: ", elsa_w10_file)
w10_vars <- c(
  "idauniq", "iintdatm", "iintdaty", "indager", "indsex",
  "HEEverBP", "HENvrBP", "HEEverAN", "HENvrAN", "HEEverHF", "HENvrHF",
  "HEEverHM", "HENvrHM", "HEEverAH", "HENvrAH", "HEEverOH", "HENvrOH",
  "HEEverST", "HENvrST",
  "HEAgeMI", "HEOthHA", "HEMedHA", "HEEverDI", "HEEverCL"
)
w10_raw <- haven::read_dta(elsa_w10_file, col_select = dplyr::any_of(w10_vars))
w10 <- as.data.table(to_numeric_frame(w10_raw, names(w10_raw)))
setnames(w10, tolower(names(w10)))
w10[, cohort := "ELSA"]
w10[, id := idauniq]
w10[, wave := 10L]
w10[, interview_year := valid_range(iintdaty, 1990, 2030)]
w10[, interview_month := valid_range(iintdatm, 1, 12)]
w10[, event_time := interview_year + ifelse(is.na(interview_month), 0.5, (interview_month - 0.5) / 12)]
w10[, event_type := "visit"]
w10[, dead := 0L]
w10[, age := valid_range(indager, 18, 110)]
w10[indager == -7, age := 90]
w10[, sex_male := fifelse(indsex %in% 1:2, as.integer(indsex == 1), NA_integer_)]
w10[, hypertension_raw := ever_with_feedforward(heeverbp, henvrbp)]
w10[, stroke_raw := ever_with_feedforward(heeverst, henvrst)]
# Harmonized ELSA defines doctor-diagnosed heart problems broadly: angina,
# myocardial infarction/heart attack, congestive heart failure, heart murmur,
# abnormal heart rhythm, or other heart trouble.  Direct HEEver responses are
# combined with HENvr=1 confirmations of diagnoses carried forward from prior
# waves.  Wave 10 does not archive a
# direct HEEverMI field.  Instead, HEAgeMI, HEOthHA and HEMedHA are routed only
# to respondents whose interview establishes a prior heart attack.  Therefore
# any substantive response to one of those routed fields (including a
# don't-know/refused answer about age, recurrence or treatment) is affirmative
# evidence of heart-attack history; -1/-2 denote non-applicable routing, not a
# negative answer about treatment alone.
w10[, heart_angina_raw := ever_with_feedforward(heeveran, henvran)]
w10[, heart_failure_raw := ever_with_feedforward(heeverhf, henvrhf)]
w10[, heart_murmur_raw := ever_with_feedforward(heeverhm, henvrhm)]
w10[, heart_rhythm_raw := ever_with_feedforward(heeverah, henvrah)]
w10[, heart_other_raw := ever_with_feedforward(heeveroh, henvroh)]
w10[, heart_other_ambiguous_raw := 0L]
w10_mi_route <- cbind(
  raw_numeric(w10$heagemi), raw_numeric(w10$heothha), raw_numeric(w10$hemedha)
)
w10_mi_route_observed <- apply(w10_mi_route, 1L, function(value) {
  any(!is.na(value) & !(value %in% c(-1, -2)))
})
w10_mi_route_not_applicable <- apply(w10_mi_route, 1L, function(value) {
  all(!is.na(value) & value %in% c(-1, -2))
})
w10[, heart_mi_raw := fcase(
  w10_mi_route_observed, 1L,
  w10_mi_route_not_applicable, 0L,
  default = NA_integer_
)]
w10[, heart_aggregate_source_raw := NA_integer_]
w10[, heart_disease_raw := combine_binary_any(
  heart_angina_raw, heart_mi_raw, heart_failure_raw,
  heart_murmur_raw, heart_rhythm_raw, heart_other_raw
)]
w10[, diabetes_raw := yes1_no2(heeverdi)]
w10[, lung_disease_raw := yes1_no2(heevercl)]
w10[, death_year := NA_real_]
w10[, death_time_midyear := NA_real_]
w10[, death_time_adjusted := NA_integer_]
w10[, death_source := NA_character_]
w10[, insomnia_onset_observed := NA_integer_]
w10[, insomnia_maintenance_observed := NA_integer_]
w10[, insomnia_nonrestorative_observed := NA_integer_]
w10[, insomnia_n_observed := NA_integer_]
w10[, insomnia_score_observed := NA_integer_]
w10[, sleep_hours_observed := NA_real_]
w10[, sleep_quality_observed := NA_integer_]

elsa <- rbindlist(list(elsa, w10), use.names = TRUE, fill = TRUE)

# ELSA differs from HRS in permitting respondents to dispute a diagnosis that
# was carried forward from an earlier interview.  Harmonized ELSA backcasts
# such disputes to the earlier waves; Wave 10 is newer than the G3 harmonized
# file, so perform the same component-specific backcast here.  Only prior
# affirmative reports are changed; missing/non-interview records stay missing.
elsa_w10_dispute_specs <- list(
  hypertension = list(raw = "hypertension_raw", henvr = "henvrbp"),
  angina = list(raw = "heart_angina_raw", henvr = "henvran"),
  heart_failure = list(raw = "heart_failure_raw", henvr = "henvrhf"),
  heart_murmur = list(raw = "heart_murmur_raw", henvr = "henvrhm"),
  abnormal_rhythm = list(raw = "heart_rhythm_raw", henvr = "henvrah"),
  other_heart = list(raw = "heart_other_raw", henvr = "henvroh"),
  stroke = list(raw = "stroke_raw", henvr = "henvrst")
)
elsa_w10_dispute_audit <- rbindlist(lapply(names(elsa_w10_dispute_specs), function(component) {
  spec <- elsa_w10_dispute_specs[[component]]
  disputed_ids <- unique(w10[is_explicit_dispute(get(spec$henvr)), id])
  prior_selector <- elsa$wave < 10L & elsa$id %in% disputed_ids
  prior_positive <- prior_selector & elsa[[spec$raw]] == 1L
  prior_unknown_other_resolved <- component == "other_heart" & prior_selector &
    is.na(elsa[[spec$raw]]) & elsa$heart_aggregate_source_raw == 1L
  out <- data.table(
    component = component,
    henvr_variable = spec$henvr,
    disputed_ids = length(disputed_ids),
    prior_positive_rows_backcast_to_no = sum(prior_positive, na.rm = TRUE),
    prior_positive_ids_backcast_to_no = uniqueN(elsa$id[which(prior_positive %in% TRUE)]),
    prior_ambiguous_other_rows_resolved_to_no = sum(prior_unknown_other_resolved, na.rm = TRUE)
  )
  if (any(prior_positive, na.rm = TRUE)) {
    set(elsa, i = which(prior_positive %in% TRUE), j = spec$raw, value = 0L)
  }
  if (any(prior_unknown_other_resolved, na.rm = TRUE)) {
    set(elsa, i = which(prior_unknown_other_resolved %in% TRUE), j = spec$raw, value = 0L)
    set(elsa, i = which(prior_unknown_other_resolved %in% TRUE), j = "heart_other_ambiguous_raw", value = 0L)
  }
  out
}))
elsa_w10_unlabelled_audit <- rbindlist(lapply(
  c("henvrbp", "henvran", "henvrhf", "henvrhm", "henvrah", "henvroh", "henvrst"),
  function(variable) data.table(
    variable = variable,
    code_85 = sum(raw_numeric(w10[[variable]]) == 85, na.rm = TRUE),
    code_86 = sum(raw_numeric(w10[[variable]]) == 86, na.rm = TRUE)
  )
))
fwrite(
  elsa_w10_dispute_audit,
  file.path(qc_output, "elsa_w10_dispute_backcast_audit.csv"),
  na = ""
)
fwrite(
  elsa_w10_unlabelled_audit,
  file.path(qc_output, "elsa_w10_unlabelled_henvr_codes.csv"),
  na = ""
)

# Rebuild the broad heart field after component-specific dispute backcasting.
# The archived Wave-10 file has no HEEverMI/HENvrMI pair, so historical MI
# disputes cannot be backcast and are retained as an explicit source limit.
elsa[, heart_disease_raw := combine_binary_any(
  heart_angina_raw, heart_mi_raw, heart_failure_raw,
  heart_murmur_raw, heart_rhythm_raw, heart_other_raw
)]
elsa[, `:=`(
  heart_other_lower_raw = fifelse(
    heart_other_ambiguous_raw == 1L & is.na(heart_other_raw), 0L, heart_other_raw
  ),
  heart_other_upper_raw = fifelse(
    heart_other_ambiguous_raw == 1L & is.na(heart_other_raw), 1L, heart_other_raw
  )
)]
elsa[, heart_disease_raw_lower := combine_binary_any(
  heart_angina_raw, heart_mi_raw, heart_failure_raw,
  heart_murmur_raw, heart_rhythm_raw, heart_other_lower_raw
)]
elsa[, heart_disease_raw_upper := combine_binary_any(
  heart_angina_raw, heart_mi_raw, heart_failure_raw,
  heart_murmur_raw, heart_rhythm_raw, heart_other_upper_raw
)]
setorder(elsa, id, event_time, wave)
elsa[, sleep_measurement_wave := as.integer(wave %in% c(4L, 6L, 8L))]
elsa[, insomnia_source_wave := fifelse(!is.na(insomnia_score_observed), wave, NA_integer_)]
elsa[, `:=`(
  insomnia_onset = insomnia_onset_observed,
  insomnia_maintenance = insomnia_maintenance_observed,
  insomnia_nonrestorative = insomnia_nonrestorative_observed,
  insomnia_score = insomnia_score_observed
)]
elsa[, c("insomnia_onset", "insomnia_maintenance", "insomnia_nonrestorative", "insomnia_score", "insomnia_source_wave") :=
       lapply(.SD, locf_down),
     by = id,
     .SDcols = c("insomnia_onset", "insomnia_maintenance", "insomnia_nonrestorative", "insomnia_score", "insomnia_source_wave")]
elsa[, insomnia_locf_reason := fifelse(
  !is.na(insomnia_score_observed), "observed",
  fifelse(is.na(insomnia_score), "missing",
          fifelse(sleep_measurement_wave == 1L, "item_missing_locf", "design_locf"))
)]
elsa[, insomnia_design_imputed := as.integer(insomnia_locf_reason == "design_locf")]
elsa[, insomnia_general_missing := as.integer(sleep_measurement_wave == 1L & is.na(insomnia_score_observed))]
elsa[, insomnia_high := fifelse(!is.na(insomnia_score), as.integer(insomnia_score >= 2L), NA_integer_)]
elsa <- add_cumulative_diseases_and_state(elsa)
elsa[, heart_disease_lower := cum_diagnosis(heart_disease_raw_lower), by = id]
elsa[, heart_disease_upper := cum_diagnosis(heart_disease_raw_upper), by = id]
elsa[, state_heart_lower := derive_cvd_state(hypertension, heart_disease_lower, stroke, dead)]
elsa[, state_heart_upper := derive_cvd_state(hypertension, heart_disease_upper, stroke, dead)]

elsa_eol_file <- file.path(stage_root, "elsa_eol", "h_elsa_eol_a2.dta")
elsa_deaths <- data.table(id = numeric(), death_year = numeric(), death_source = character())
if (file.exists(elsa_eol_file)) {
  eol <- haven::read_dta(elsa_eol_file, col_select = dplyr::any_of(c("idauniqc", "raxyear")))
  eol <- as.data.table(to_numeric_frame(eol, names(eol)))
  elsa_deaths <- eol[, .(
    id = idauniqc,
    death_year = valid_range(raxyear, 1900, 2030),
    death_source = "Harmonized ELSA EOL A2 raxyear (coverage through 2012)"
  )]
}
elsa <- append_death_rows(elsa, elsa_deaths, study_end_year = 2023)
elsa[dead == 1L, `:=`(state_heart_lower = 4L, state_heart_upper = 4L)]

message("Cleaning CHARLS...")
charls_package <- read_package_data("CHARLS", object_name = "charls")
charls_filtered <- GlobalAgingR::CH_filter(charls_package, wave %in% 1:5)
charls_vars <- c(
  "ID", "wave", "iwm", "iwy", "age", "ragender", "raeduc_c", "raeducl", "marry",
  "hcoresd", "hibpe", "hearte", "stroke", "diabe", "lunge", "sleep_night", "systo",
  "diasto", "rxhibp", "srh", "smokev", "smoken", "drinkl", "vgactx_c", "mdactx_c",
  "ltactx_c", "bmi", "cesd10", "radyear", "iwstat"
)
charls <- as.data.table(to_numeric_frame(charls_filtered, charls_vars))
charls[, cohort := "CHARLS"]
charls[, id := ID]
charls[, wave := as.integer(wave)]
charls[, interview_year := valid_range(iwy, 1990, 2030)]
charls_year_map <- c(`1` = 2011, `2` = 2013, `3` = 2015, `4` = 2018, `5` = 2020)
charls[is.na(interview_year), interview_year := unname(charls_year_map[as.character(wave)])]
charls[, interview_month := valid_range(iwm, 1, 12)]
charls[, event_time := interview_year + ifelse(is.na(interview_month), 0.5, (interview_month - 0.5) / 12)]
charls[, event_type := "visit"]
charls[, dead := 0L]
charls[, age := valid_range(age, 18, 110)]
charls[, sex_male := binary01(ragender)]
charls[, education_code := raeduc_c]
charls[, marital_status := marry]
charls[, coresident_child := binary01(hcoresd)]
charls[, self_rated_health := valid_range(srh, 1, 5)]
charls[, hypertension_raw := binary01(hibpe)]
charls[, heart_disease_raw := binary01(hearte)]
# GlobalAgingR 2.0.0 retains the original CHARLS Wave 4 no=2 coding for stroke.
charls[, stroke_raw := fifelse(stroke == 1, 1L, fifelse(stroke %in% c(0, 2), 0L, NA_integer_))]
charls[, diabetes_raw := binary01(diabe)]
charls[, lung_disease_raw := binary01(lunge)]
charls[, sleep_hours_observed := valid_range(sleep_night, 1, 18)]
charls[, sleep_hours_centered7 := sleep_hours_observed - 7]
charls[, systolic_bp_raw := systo]
charls[, diastolic_bp_raw := diasto]
charls[, systolic_bp := valid_range(systo, 70, 260)]
charls[, diastolic_bp := valid_range(diasto, 40, 150)]
charls[!is.na(systolic_bp) & !is.na(diastolic_bp) & systolic_bp <= diastolic_bp,
       c("systolic_bp", "diastolic_bp") := .(NA_real_, NA_real_)]
charls[, pulse_pressure := systolic_bp - diastolic_bp]
charls[, antihypertensive_medication := binary01(rxhibp)]
charls[, ever_smoked := binary01(smokev)]
charls[, current_smoker := binary01(smoken)]
charls[, current_drinker := binary01(drinkl)]
charls[, vigorous_activity_days := valid_range(vgactx_c, 0, 7)]
charls[, moderate_activity_days := valid_range(mdactx_c, 0, 7)]
charls[, light_activity_days := valid_range(ltactx_c, 0, 7)]
charls[, bmi := valid_range(bmi, 10, 80)]
charls[, depression_score := valid_range(cesd10, 0, 30)]
charls[, insomnia_score_observed := NA_integer_]
charls[, insomnia_score := NA_integer_]
charls[, death_year := valid_range(radyear, 1900, 2030)]
charls[, death_time_midyear := NA_real_]
charls[, death_time_adjusted := NA_integer_]
charls[, death_source := NA_character_]
charls <- add_cumulative_diseases_and_state(charls)
charls_deaths <- charls[!is.na(death_year), .(
  death_year = first_nonmissing(death_year), death_source = "GlobalAgingR CHARLS radyear"
), by = id]
charls <- append_death_rows(charls, charls_deaths, study_end_year = 2020)

hrs_main <- make_main_analysis(hrs)
elsa_main <- make_main_analysis(elsa)

make_mechanism_sets <- function(dt, waves, cohort_name) {
  vars <- c(
    "cohort", "id", "wave", "interview_year", "interview_month", "event_time", "age",
    "sex_male", "education_years", "education_code", "marital_status", "bmi",
    "sleep_hours_observed", "sleep_hours_centered7", "sleep_quality_observed",
    "systolic_bp", "diastolic_bp", "pulse_pressure", "antihypertensive_medication",
    "current_smoker", "current_drinker", "depression_score", "vigorous_activity",
    "moderate_activity", "light_activity", "vigorous_activity_days",
    "moderate_activity_days", "light_activity_days"
  )
  x <- as.data.table(copy(dt))[event_type == "visit" & wave %in% waves]
  x <- x[, intersect(vars, names(x)), with = FALSE]
  x <- x[!is.na(age) & age >= 50]
  eligible_ids <- x[, .(
    n_waves = uniqueN(wave),
    any_sleep = any(!is.na(sleep_hours_observed)),
    any_sbp = any(!is.na(systolic_bp))
  ), by = id][n_waves >= 2L & any_sleep & any_sbp, id]
  fiml <- x[id %in% eligible_ids]
  setorder(fiml, id, wave)
  complete_ids <- fiml[, .(
    all_waves = all(waves %in% wave),
    complete_sleep_sbp = all(!is.na(sleep_hours_observed[match(waves, wave)])) &&
      all(!is.na(systolic_bp[match(waves, wave)]))
  ), by = id][all_waves & complete_sleep_sbp, id]
  complete_input <- as.data.frame(fiml[id %in% complete_ids])
  if (nrow(complete_input)) {
    # GlobalAgingR::CH_extract_patient_waves() coerces patient_ids to character.
    # Match that type explicitly so its internal anti_join remains type-safe.
    complete_input$ID <- as.character(complete_input$id)
    complete3 <- suppressWarnings(GlobalAgingR::CH_extract_patient_waves(
      complete_input,
      patient_ids = unique(as.character(complete_input$ID)),
      wave_numbers = waves,
      id_col = "ID",
      wave_col = "wave"
    ))
    complete3$ID <- NULL
    complete3 <- as.data.table(complete3)
  } else {
    complete3 <- as.data.table(complete_input)
  }
  message(cohort_name, " mechanism FIML IDs: ", uniqueN(fiml$id),
          "; complete three-wave IDs: ", uniqueN(complete3$id))
  list(fiml = fiml, complete3 = complete3)
}

charls_mechanism <- make_mechanism_sets(charls, c(1L, 2L, 3L), "CHARLS")
elsa_mechanism <- make_mechanism_sets(elsa, c(4L, 6L, 8L), "ELSA")

save_dataset <- function(x, stem, write_csv = TRUE) {
  x <- as.data.table(copy(x))
  setorder(x, id, event_time, na.last = TRUE)
  saveRDS(x, file.path(data_output, paste0(stem, ".rds")), compress = "gzip")
  if (write_csv) {
    data.table::fwrite(x, file.path(data_output, paste0(stem, ".csv.gz")), compress = "gzip", na = "")
  }
}

save_dataset(hrs, "hrs_clean")
save_dataset(elsa, "elsa_clean")
save_dataset(charls, "charls_clean")
save_dataset(hrs_main, "hrs_main_analysis")
save_dataset(elsa_main, "elsa_main_analysis")
save_dataset(charls_mechanism$fiml, "charls_mechanism_fiml")
save_dataset(charls_mechanism$complete3, "charls_mechanism_complete3")
save_dataset(elsa_mechanism$fiml, "elsa_mechanism_fiml")
save_dataset(elsa_mechanism$complete3, "elsa_mechanism_complete3")

wave_summary <- function(dt) {
  x <- as.data.table(copy(dt))[event_type == "visit"]
  x[, .(
    n_rows = .N,
    n_ids = uniqueN(id),
    median_age = median(age, na.rm = TRUE),
    n_state_observed = sum(state %in% 0:3, na.rm = TRUE),
    n_insomnia_observed = sum(!is.na(insomnia_score_observed)),
    mean_insomnia_score = mean(insomnia_score_observed, na.rm = TRUE),
    n_sleep_hours = if ("sleep_hours_observed" %in% names(x)) sum(!is.na(sleep_hours_observed)) else 0L,
    n_systolic_bp = if ("systolic_bp" %in% names(x)) sum(!is.na(systolic_bp)) else 0L,
    hypertension_prevalence = mean(hypertension, na.rm = TRUE),
    heart_disease_prevalence = mean(heart_disease, na.rm = TRUE),
    stroke_prevalence = mean(stroke, na.rm = TRUE),
    state_s0 = sum(state == 0L, na.rm = TRUE),
    state_s1 = sum(state == 1L, na.rm = TRUE),
    state_s2 = sum(state == 2L, na.rm = TRUE),
    state_s3 = sum(state == 3L, na.rm = TRUE)
  ), by = .(cohort, wave, interview_year)]
}

transition_summary <- function(dt) {
  x <- as.data.table(copy(dt))[!is.na(state)]
  setorder(x, id, event_time, wave)
  x[, next_state := shift(state, type = "lead"), by = id]
  x[, next_time := shift(event_time, type = "lead"), by = id]
  x[!is.na(next_state), .(
    n_intervals = .N,
    n_ids = uniqueN(id),
    median_interval_years = median(next_time - event_time, na.rm = TRUE)
  ), by = .(cohort, from_state = state, to_state = next_state)]
}

quality_checks <- function(dt) {
  x <- as.data.table(copy(dt))
  setorder(x, id, event_time, wave)
  live <- x[event_type == "visit"]
  living_states <- live[!is.na(state), .(states = list(state)), by = id]
  nonmonotone <- living_states[, sum(vapply(states, function(z) any(diff(z) < 0), logical(1)))]
  data.table(
    cohort = x$cohort[1],
    check = c(
      "duplicate_id_wave_visits", "duplicate_id_event_time", "age_outside_18_110",
      "state_outside_0_4", "living_state_reversals", "death_not_last",
      "insomnia_score_outside_0_3", "clean_sbp_outside_70_260",
      "clean_dbp_outside_40_150", "clean_sbp_not_above_dbp"
    ),
    value = c(
      live[, sum(duplicated(paste(id, wave)))],
      x[, sum(duplicated(paste(id, sprintf("%.8f", event_time))))],
      x[, sum(!is.na(age) & (age < 18 | age > 110))],
      x[, sum(!is.na(state) & !state %in% 0:4)],
      nonmonotone,
      x[, sum(dead == 1L & seq_len(.N) != .N), by = id][, sum(V1)],
      x[, sum(!is.na(insomnia_score) & !insomnia_score %in% 0:3)],
      if ("systolic_bp" %in% names(x)) x[, sum(!is.na(systolic_bp) & (systolic_bp < 70 | systolic_bp > 260))] else 0,
      if ("diastolic_bp" %in% names(x)) x[, sum(!is.na(diastolic_bp) & (diastolic_bp < 40 | diastolic_bp > 150))] else 0,
      if (all(c("systolic_bp", "diastolic_bp") %in% names(x))) x[, sum(!is.na(systolic_bp) & !is.na(diastolic_bp) & systolic_bp <= diastolic_bp)] else 0
    )
  )
}

wave_qc <- rbindlist(lapply(list(hrs, elsa, charls), wave_summary), fill = TRUE)
transition_qc <- rbindlist(lapply(list(hrs_main, elsa_main), transition_summary), fill = TRUE)
checks_qc <- rbindlist(lapply(list(hrs, elsa, charls), quality_checks), fill = TRUE)
data.table::fwrite(wave_qc, file.path(qc_output, "cohort_wave_summary.csv"), na = "")
data.table::fwrite(transition_qc, file.path(qc_output, "observed_transition_summary.csv"), na = "")
data.table::fwrite(checks_qc, file.path(qc_output, "quality_checks.csv"), na = "")

source_manifest <- data.table(
  cohort = c("HRS", "ELSA", "ELSA", "ELSA", "ELSA", "ELSA", "CHARLS"),
  source = c(
    "GlobalAgingR::hrs", "GlobalAgingR::elsa",
    "Harmonized_ELSA_G3/h_elsa_g3.dta",
    "data_raw/2_ELSA_UK/wave4, wave6, wave8 core data",
    "data_raw/2_ELSA_UK/wave10/wave_10_elsa_data_eul_v1.dta",
    "hamonized data/Harmonized ELSA EOL A2",
    "GlobalAgingR::CHARLS (object charls)"
  ),
  role = c(
    "Harmonized longitudinal spine, insomnia, CVD and covariates",
    "Harmonized W4-W9 longitudinal spine, CVD, BP and covariates",
    paste(
      "Authoritative W4-W9 six-component heart phenotype and aggregate flag;",
      "also supplies W5 angina values absent from GlobalAgingR::elsa"
    ),
    "Original heslpa/heslpb/heslpd/heslpe/heslpf sleep items",
    paste(
      "Wave 10 visit extension; heart problems combine angina, heart attack,",
      "heart failure, heart murmur, abnormal rhythm and other heart trouble"
    ),
    "Death year supplement; available deaths end in 2012",
    "Harmonized W1-W5 longitudinal spine, sleep duration, BP, CVD and covariates"
  ),
  limitations = c(
    "Death date is represented by reported death-year midpoint",
    "Package longitudinal file ends at W9",
    paste(
      "A positive aggregate flag cannot reveal coexisting residual other heart trouble",
      "when a named component is already positive; primary=unknown, lower/upper bounds retained"
    ),
    "Sleep symptom batteries available at W4/W6/W8 only",
    paste(
      "HENvr component codes 2/3/4 backcast prior positives and 85/86 remain unknown;",
      "no archived HEEverMI/HENvrMI, so MI routing is recoverable but MI dispute backcast is not"
    ),
    "ELSA deaths after 2012 are not available in supplied files",
    "Wave 4 stroke no=2 is corrected to 0"
  )
)
data.table::fwrite(source_manifest, file.path(qc_output, "source_manifest.csv"), na = "")

variable_dictionary <- data.table(
  variable = c(
    "id", "wave", "event_time", "time_since_entry", "state", "hypertension",
    "heart_disease", "heart_disease_lower", "heart_disease_upper",
    "state_heart_lower", "state_heart_upper",
    "heart_angina_raw", "heart_mi_raw", "heart_failure_raw", "heart_murmur_raw",
    "heart_rhythm_raw", "heart_other_raw", "heart_other_ambiguous_raw",
    "stroke", "insomnia_onset", "insomnia_maintenance",
    "insomnia_nonrestorative", "insomnia_score", "insomnia_score_observed",
    "insomnia_source_wave", "insomnia_locf_reason", "sleep_hours_observed",
    "sleep_hours_centered7", "systolic_bp", "diastolic_bp", "pulse_pressure",
    "death_year", "death_time_midyear", "death_time_adjusted"
  ),
  definition = c(
    "Cohort-specific stable person identifier",
    "Original survey wave; 999 denotes an appended death row",
    "Decimal interview/death time in calendar years",
    "Years from first eligible age>=50 state observation",
    "0 no CVD; 1 hypertension only; 2 one heart/stroke class without hypertension; 3 >=2 CVD classes; 4 death",
    "Cumulative doctor-diagnosed hypertension history",
    "Cumulative doctor-diagnosed heart disease history",
    "ELSA cumulative heart history under the lower-bound residual-other rule",
    "ELSA cumulative heart history under the upper-bound residual-other rule",
    "Living/death state recomputed with heart_disease_lower",
    "Living/death state recomputed with heart_disease_upper",
    "Visit-level angina component; ELSA only",
    "Visit-level myocardial-infarction/heart-attack component; ELSA only",
    "Visit-level congestive-heart-failure component; ELSA only",
    "Visit-level heart-murmur component; ELSA only",
    "Visit-level abnormal-heart-rhythm component; ELSA only",
    "Visit-level identifiable residual other-heart-trouble component; ELSA only",
    "ELSA flag for aggregate-positive visits where coexisting other trouble is unidentifiable",
    "Cumulative doctor-diagnosed stroke history",
    "Binary difficulty initiating sleep",
    "Binary waking repeatedly during the night",
    "Binary non-restorative sleep",
    "Sum of three harmonized insomnia symptoms, 0-3",
    "Insomnia score at a wave where all three source items were observed",
    "ELSA wave supplying the carried-forward insomnia score",
    "ELSA exposure provenance: observed/design_locf/item_missing_locf/missing",
    "Reported weekday/night sleep duration in hours; valid range 1-18",
    "Sleep hours minus 7; mechanism interpretation center",
    "Clean mean systolic BP, 70-260 mmHg and greater than DBP",
    "Clean mean diastolic BP, 40-150 mmHg and lower than SBP",
    "Systolic minus diastolic BP",
    "Reported calendar year of death",
    "Death year plus 0.5 for the primary interval approximation",
    "1 if death midpoint was moved just after the last observed visit to preserve ordering"
  ),
  coding = c(
    "numeric", "integer", "numeric", "numeric", "0-4", "0/1", "0/1",
    "0/1", "0/1", "0-4", "0-4", "0/1/NA", "0/1/NA", "0/1/NA",
    "0/1/NA", "0/1/NA", "0/1/NA", "0/1", "0/1",
    "0/1", "0/1", "0/1", "0-3", "0-3", "integer", "text", "hours", "hours",
    "mmHg", "mmHg", "mmHg", "year", "decimal year", "0/1"
  )
)
data.table::fwrite(variable_dictionary, file.path(qc_output, "variable_dictionary.csv"), na = "")

dataset_counts <- rbindlist(list(
  data.table(dataset = "hrs_clean", rows = nrow(hrs), ids = uniqueN(hrs$id)),
  data.table(dataset = "elsa_clean", rows = nrow(elsa), ids = uniqueN(elsa$id)),
  data.table(dataset = "charls_clean", rows = nrow(charls), ids = uniqueN(charls$id)),
  data.table(dataset = "hrs_main_analysis", rows = nrow(hrs_main), ids = uniqueN(hrs_main$id)),
  data.table(dataset = "elsa_main_analysis", rows = nrow(elsa_main), ids = uniqueN(elsa_main$id)),
  data.table(dataset = "charls_mechanism_fiml", rows = nrow(charls_mechanism$fiml), ids = uniqueN(charls_mechanism$fiml$id)),
  data.table(dataset = "charls_mechanism_complete3", rows = nrow(charls_mechanism$complete3), ids = uniqueN(charls_mechanism$complete3$id)),
  data.table(dataset = "elsa_mechanism_fiml", rows = nrow(elsa_mechanism$fiml), ids = uniqueN(elsa_mechanism$fiml$id)),
  data.table(dataset = "elsa_mechanism_complete3", rows = nrow(elsa_mechanism$complete3), ids = uniqueN(elsa_mechanism$complete3$id))
))
data.table::fwrite(dataset_counts, file.path(qc_output, "dataset_counts.csv"), na = "")

summary_lines <- c(
  "# Stage 1 cleaning summary",
  "",
  paste0("Run time: ", format(run_started, "%Y-%m-%d %H:%M:%S %Z")),
  paste0("GlobalAgingR: ", package_version),
  "",
  "## Delivered datasets",
  "",
  paste0("- ", dataset_counts$dataset, ": ", format(dataset_counts$rows, big.mark = ","),
         " rows, ", format(dataset_counts$ids, big.mark = ","), " people"),
  "",
  "## Frozen cleaning decisions",
  "",
  "- HRS: waves 8-15 (2006-2020); three-item insomnia score uses most-of-time onset, most-of-time night waking, and rarely/never rested.",
  "- ELSA: waves 4-10; positive insomnia items are three-or-more-times/week at waves 4/6/8; the latest complete score is carried forward with provenance flags.",
  "- ELSA heart disease is reconstructed from angina, myocardial infarction/heart attack, congestive heart failure, heart murmur, abnormal rhythm, and other heart trouble. W4-W9 use the Harmonized ELSA G3 source (including W5 angina missing from the package); W10 component disputes 2/3/4 backcast prior positives, while 85/86 remain unknown.",
  "- When a positive historical aggregate cannot identify coexisting residual other heart trouble, the primary phenotype preserves unknown; lower- and upper-bound heart histories and states are retained for protocol sensitivity analyses.",
  "- The archived W10 file has no HEEverMI/HENvrMI pair: substantive myocardial-infarction routing confirms MI, but an MI-specific historical dispute cannot be backcast.",
  "- CHARLS: waves 1-5; sleep duration is restricted to 1-18 hours and BP to SBP 70-260 / DBP 40-150 mmHg with SBP>DBP.",
  "- CVD diagnoses are cumulative within person. The four living states and absorbing death state follow the supplied study plan.",
  "- Main analysis files enter follow-up at the first age>=50 valid state and require >=2 state observations plus >=1 insomnia score.",
  "",
  "## Known source limitation",
  "",
  "The supplied Harmonized ELSA EOL A2 file contains death years only through 2012. ELSA mortality after 2012 must be refreshed before the final multi-state death analysis; this is recorded in source_manifest.csv and is not silently treated as complete."
)
writeLines(summary_lines, file.path(qc_output, "cleaning_summary.md"), useBytes = TRUE)

run_metadata <- data.table(
  run_started = format(run_started, "%Y-%m-%d %H:%M:%S %Z"),
  run_finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  R_version = R.version.string,
  GlobalAgingR_version = package_version,
  stage_root = normalizePath(stage_root, winslash = "/", mustWork = TRUE)
)
data.table::fwrite(run_metadata, file.path(log_output, "run_metadata.csv"), na = "")
capture.output(sessionInfo(), file = file.path(log_output, "session_info.txt"))

message("Cleaning completed. Output: ", normalizePath(output_root, winslash = "/", mustWork = TRUE))
