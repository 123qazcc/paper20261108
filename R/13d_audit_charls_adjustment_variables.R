#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table); library(haven); library(digest)})
out <- "outputs/charls_external_extension_20260913/adjustment_audit"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
env <- new.env(); data("CHARLS", package = "GlobalAgingR", envir = env)
p <- as.data.table(env$charls)
p[, id := sprintf("%012.0f", as.numeric(ID))]
b <- as.data.table(readRDS("outputs/charls_external_extension_20260913/data/charls_w1_descriptive_baseline.rds"))
demo_path <- "private/stage/charls_phase2_root/w1_demographic.dta"
stopifnot(digest(file = demo_path, algo = "sha256") == "46d94759c9596b8198237c9460af2dc0550d8ef93812575826cc33279271cac0")
raw <- as.data.table(read_dta(demo_path, col_select = c("ID", "bd001")))
raw[, id := paste0(substr(ID, 1, 9), "0", substr(ID, 10, 11))]
raw[, original_education := as.numeric(bd001)]
raw[, education4_raw := fcase(original_education %in% 1:3, 1L, original_education == 4, 2L,
  original_education == 5, 3L, original_education %in% 6:11, 4L, default = NA_integer_)]
raw[, education3_raw := fcase(original_education %in% 1:5, 1L, original_education %in% 6:7, 2L,
  original_education %in% 8:11, 3L, default = NA_integer_)]
j <- merge(p[wave == 1, .(id, education4_package = as.numeric(raeduc_c), education3_package = as.numeric(raeducl))],
  raw[, .(id, original_education, education4_raw, education3_raw)], by = "id", all.x = TRUE)
fwrite(j[, .N, by = .(original_education, education4_package, education3_package)], file.path(out, "education_crosswalk_counts.csv"))
eq_na <- function(a,b) all((is.na(a)&is.na(b)) | (!is.na(a)&!is.na(b)&a==b))
comparison <- data.table(comparison=c("raw_W1_education4_vs_harmonized", "raw_W1_education3_vs_harmonized"),
  exact_match=c(eq_na(j$education4_raw,j$education4_package),eq_na(j$education3_raw,j$education3_package)),
  value_differences=c(sum(j$education4_raw!=j$education4_package,na.rm=TRUE),sum(j$education3_raw!=j$education3_package,na.rm=TRUE)),
  missingness_differences=c(sum(xor(is.na(j$education4_raw),is.na(j$education4_package))),sum(xor(is.na(j$education3_raw),is.na(j$education3_package)))))
fwrite(comparison,file.path(out,"raw_harmonized_education_comparison.csv"))
checks <- data.table(check = c("raw_education4_valid_range", "raw_education3_valid_range", "raw_unique_id"),
  pass = c(all(j$education4_raw %in% 1:4 | is.na(j$education4_raw)),all(j$education3_raw %in% 1:3 | is.na(j$education3_raw)),!anyDuplicated(raw$id)))
fwrite(checks,file.path(out,"adjustment_checks.csv"))
vars <- c("ID", "id", "raeduc_c", "raeducl", "hatotfa", "hhcperc", "hrural", "nation", "marry",
  "phy_acta", "phy_actb", "vgact_c", "mdact_c", "ltact_c", paste0("da042s",1:15))
z <- p[wave == 1, ..vars]
z[, setdiff(vars,c("ID","id")) := lapply(.SD, as.numeric), .SDcols = setdiff(vars,c("ID","id"))]
pain_matrix <- as.matrix(z[,paste0("da042s",1:15),with=FALSE])
z[, any_pain := fifelse(rowSums(pain_matrix == 1,na.rm=TRUE)>0,1L,
  fifelse(rowSums(pain_matrix==0,na.rm=TRUE)==15,0L,NA_integer_))]
z[raw,on="id",`:=`(education4_W1=i.education4_raw,education3_W1=i.education3_raw)]
saveRDS(z,file.path(out,"all_W1_adjustment_supplement.rds"))
z <- z[id %in% b$id]
available <- rbindlist(lapply(c(setdiff(vars,c("ID","id")),"any_pain","education4_W1","education3_W1"),function(v) data.table(variable=v,
  n=nrow(z),n_observed=sum(!is.na(z[[v]])),n_missing=sum(is.na(z[[v]])),missing_pct=100*mean(is.na(z[[v]])))))
fwrite(available,file.path(out,"baseline_adjustment_missingness.csv"))
saveRDS(z,file.path(out,"baseline_adjustment_supplement.rds"))
writeLines(c("# CHARLS baseline adjustment variable audit", "",
  "Raw W1 bd001 supports four categories: less than primary (1-3); primary (4); middle school (5); high school or above (6-11).",
  "The harmonized education fields do not match this raw W1 answer person by person. Reasons may include cross-wave reconciliation; that process was not established by this audit.",
  "Use education4_W1 derived only from original W1 answers for this extension, preserving missingness. Do not treat the phase-1 literal package label as an independently validated classification.",
  "Physical activity is unavailable for most W1 participants in the package; retain the extent and examine sampling/routing before adjustment. Do not label missing activity as inactive.",
  "hatotfa is labelled non-housing financial wealth, not total household wealth; negative values represent a different scale than missingness and are preserved.",
  "No association model or p-value selection was performed."),file.path(out,"adjustment_audit.md"))
print(checks);print(available)
if(!all(checks$pass))stop("Raw education derivation checks failed")
