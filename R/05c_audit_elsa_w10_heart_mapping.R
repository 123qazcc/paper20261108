#!/usr/bin/env Rscript

# Compact reporting layer for the current ELSA heart-phenotype audit.
#
# The independent reconstruction itself lives in R/05b_verify_disease_states.R:
# it reads the Harmonized ELSA G3 source and the Wave 10 interviewer file,
# rebuilds all six components, applies component-specific W10 dispute backcasts,
# and compares the result with the cleaned and analysis datasets. Keeping one
# authoritative reconstruction prevents this report from drifting back to the
# obsolete direct-HEEver-only rule.

suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE, width = 160)

stage_root <- Sys.getenv("GA_STAGE_ROOT", unset = "private/stage")
output_root <- Sys.getenv(
  "GA_MULTISTATE_OUTPUT_ROOT",
  unset = file.path(stage_root, "multistate_output")
)
multistate_qc <- file.path(stage_root, "multistate_output", "qc")
clean_qc <- file.path(stage_root, "output", "qc")
report_qc <- file.path(output_root, "qc")
dir.create(report_qc, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  verification = file.path(multistate_qc, "disease_state_verification.csv"),
  backcast = file.path(multistate_qc, "elsa_harmonized_component_backcast_verification.csv"),
  dispute = file.path(clean_qc, "elsa_w10_dispute_backcast_audit.csv"),
  unlabelled = file.path(clean_qc, "elsa_w10_unlabelled_henvr_codes.csv")
)
missing_paths <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing_paths)) {
  stop(
    "Current ELSA heart audit inputs are missing. Run R/05b_verify_disease_states.R first: ",
    paste(missing_paths, collapse = "; ")
  )
}

verification <- fread(paths$verification)
backcast <- fread(paths$backcast)
dispute <- fread(paths$dispute)
unlabelled <- fread(paths$unlabelled)

elsa_heart_checks <- verification[
  cohort == "ELSA" &
    (grepl("HEART", check_id) | source %chin% c(
      "w10_raw_heart_mapping", "elsa_harmonized_source_backcast"
    ))
]
failed_checks <- elsa_heart_checks[status != "PASS" | violations != 0L, .N]

component_order <- c("angina", "heart_failure", "heart_murmur", "heart_rhythm", "other_heart")
metadata <- data.table(
  component = c(component_order, "myocardial_infarction", "broad_heart_disease"),
  wave10_rule = c(
    rep("HEEver plus HENvr: 1 confirms; 2/3/4 disputes and backcasts; 85/86 is unknown", 5L),
    "Any substantive HEAgeMI/HEOthHA/HEMedHA route confirms; all -1/-2 is no route; no MI dispute field is archived",
    "OR of angina, MI/heart attack, CHF, murmur, abnormal rhythm, and identifiable other heart trouble"
  ),
  historical_rule = c(
    rep("Harmonized ELSA G3 component history with component-specific W10 backcast", 5L),
    "Harmonized MI history retained because HENvrMI is unavailable",
    "Primary tri-state residual-other rule; lower and upper bounds retained"
  )
)

metric <- rbindlist(list(
  data.table(
    audit_family = "independent_verification",
    metric = c("ELSA heart-related checks", "failed checks", "backcast comparison rows", "backcast mismatches"),
    n = c(
      nrow(elsa_heart_checks), failed_checks, nrow(backcast),
      if ("mismatch" %in% names(backcast)) sum(backcast$mismatch, na.rm = TRUE) else
        sum(backcast$status != "PASS", na.rm = TRUE)
    )
  ),
  data.table(
    audit_family = "wave10_dispute_backcast",
    metric = c("audit rows", "component-participant disputes", "prior positive component rows set to no"),
    n = c(
      nrow(dispute),
      if ("disputed_ids" %in% names(dispute)) sum(dispute$disputed_ids, na.rm = TRUE) else NA_real_,
      if ("prior_positive_rows_backcast_to_no" %in% names(dispute)) sum(dispute$prior_positive_rows_backcast_to_no, na.rm = TRUE) else NA_real_
    )
  ),
  data.table(
    audit_family = "wave10_unlabelled_codes",
    metric = c("85/86 audit rows", "85/86 observations retained as unknown"),
    n = c(
      nrow(unlabelled),
      if (all(c("code_85", "code_86") %in% names(unlabelled)))
        sum(unlabelled$code_85, unlabelled$code_86, na.rm = TRUE) else NA_real_
    )
  )
), use.names = TRUE, fill = TRUE)

metric[, status := ifelse(failed_checks == 0L, "PASS", "FAIL")]
fwrite(metadata, file.path(report_qc, "elsa_w10_heart_variable_metadata.csv"), na = "")
fwrite(metric, file.path(report_qc, "elsa_w10_heart_mapping_summary.csv"), na = "")

print(metric)
if (failed_checks > 0L) stop("ELSA heart-disease verification contains failed checks")
message("Current ELSA Wave 10 heart audit report completed.")
