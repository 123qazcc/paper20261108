#!/usr/bin/env Rscript
# Read-only audit of published stage-3 CSV exports. No fitting, pooling,
# changes to existing fits, or figure generation; outputs are new review evidence.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table); library(jsonlite); library(digest)})
root <- "outputs/charls_dynamic_sleep_20260920"
dest <- file.path(root, "stage3_delivery")
review <- fromJSON(file.path(root, "model_qc/dynamic_independent_pool_review.json"))
stopifnot(review$status == "COMPUTATIONAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_REVIEW",
          review$checks_total == 230L, review$checks_passed == 230L)
checks <- list()
add <- function(label, value) {
  checks[[length(checks) + 1L]] <<- data.table(check = label, pass = isTRUE(value))
  stopifnot(isTRUE(value))
}
same <- function(x, y) isTRUE(all.equal(as.data.frame(x), as.data.frame(y),
                                      tolerance = 1e-12, check.attributes = TRUE))
targets <- fread(file.path(root, "results/four_target_dynamic_contrasts_MI20.csv"))
direct <- fread(file.path(root, "results/heart_stroke_dynamic_direct_MI20.csv"))
exports <- list(
  table1_dynamic_landmark_flow = fread(file.path(root, "tables/dynamic_landmark_flow.csv")),
  table2_dynamic_pattern_distribution = fread(file.path(root, "tables/dynamic_pattern_distribution.csv")),
  table3_primary_dynamic_four_paths_CORE_BOTH = targets[role == "primary"],
  table4_contextual_dynamic_four_paths_CORE_BOTH = targets[role == "contextual"],
  table5_four_level_pattern_omnibus_CORE_BOTH = fread(file.path(root, "results/four_target_pattern_omnibus_MI20.csv")),
  table6_heart_stroke_direct_dynamic_CORE_BOTH = direct,
  table7_joint_tests_CORE_BOTH = fread(file.path(root, "results/joint_tests_MI20.csv")),
  table8_CORE_BOTH_change_descriptive = fread(file.path(root, "results/matched_CORE_BOTH_change_descriptive.csv")),
  table9_observed_known_pair_support = fread(file.path(root, "tables/four_target_observed_known_pairs_by_pattern.csv"))
)
files <- file.path(dest, paste0(names(exports), ".csv"))
for (i in seq_along(exports)) {
  add(paste0(names(exports)[i], "_exists"), file.exists(files[i]))
  add(paste0(names(exports)[i], "_matches_source"), same(exports[[i]], fread(files[i])))
}
add("unique_32_target_keys", uniqueN(targets, by = c("model", "comparison", "transition")) == 32L && nrow(targets) == 32L)
add("unique_16_direct_keys", uniqueN(direct, by = c("model", "comparison", "origin")) == 16L && nrow(direct) == 16L)
add("target_BH_families_8", all(targets[, .N, by = .(model, role)]$N == 8L))
add("direct_BH_families_4", all(direct[, .N, by = .(model, role)]$N == 4L))
add("target_BH_recomputed", all(targets[, abs(p.adjust(p, "BH") - p_BH_family) < 1e-12, by = .(model, role)]$V1))
add("direct_BH_recomputed", all(direct[, abs(p.adjust(p, "BH") - p_BH_family) < 1e-12, by = .(model, role)]$V1))
add("pattern_sum_is_9397", sum(exports[[2L]]$people) == 9397L)
add("all_16_pattern_path_cells_present", nrow(exports[[9L]]) == 16L && all(exports[[9L]]$observed_known_pairs > 0L))
add("primary_supported_keys_match_CORE_BOTH", identical(
  targets[model == "CORE" & role == "primary" & p_BH_family < .05, paste(comparison, transition)],
  targets[model == "BOTH" & role == "primary" & p_BH_family < .05, paste(comparison, transition)]))
add("no_primary_direct_BH_rejection", all(direct[role == "primary", p_BH_family] >= .05))
tab <- rbindlist(checks)
fwrite(tab, file.path(dest, "table_export_review_20260929.csv"), bom = TRUE)
write_json(list(status = "EXPORTED_TABLES_MATCH_VERIFIED_NUMERIC_RESULTS",
  checks_total = nrow(tab), checks_passed = sum(tab$pass), accepted_models = 40L,
  upstream_independent_checks = 230L, tables = basename(files),
  table_sha256 = as.list(setNames(vapply(files, function(p) digest(file = p, algo = "sha256"), ""), basename(files))),
  visual_review_complete = FALSE, scheduler_terminal_cleared = FALSE,
  automatic_next_stage = FALSE), file.path(dest, "table_export_review_20260929.json"),
  pretty = TRUE, auto_unbox = TRUE)
cat("TABLE_EXPORT_AUDIT_PASS:", nrow(tab), "/", nrow(tab), "checks. No figure or final-delivery approval.\n")
