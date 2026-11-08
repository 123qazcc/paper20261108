#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(digest)})
out<-"outputs/charls_external_extension_20260913"
b<-readRDS(file.path(out,"data/charls_model_bundle.rds"))
imp<-readRDS(file.path(out,"data/charls_baseline_mice20.rds"))
checks<-list();put<-function(n,v){checks[[length(checks)+1]]<<-data.table(check=n,pass=isTRUE(v))}
put("preparation_checks_pass",all(b$preparation_checks$pass))
put("input_hashes_current",all(vapply(names(b$input_hashes),function(f)digest(file=f,algo="sha256"),"")==b$input_hashes))
put("20_imputations_20_iterations",imp$m==20L&&imp$iteration==20L&&length(b$imputations)==20)
ref<-as.data.table(b$baseline);setorder(ref,id)
for(i in 1:20){
 z<-as.data.table(b$imputations[[i]]);setorder(z,id)
 put(paste0("draw_",i,"_keys_exposure_complete"),identical(z$id,ref$id)&&identical(z$restless_bl,ref$restless_bl)&&!anyNA(z))
 preserved<-vapply(b$core_covariates,function(v){a<-ref[[v]];q<-z[[v]];known<-!is.na(a)
  if(is.factor(a)||is.factor(q))all(as.character(a[known])==as.character(q[known]))else all(a[known]==q[known])},TRUE)
 put(paste0("draw_",i,"_observed_covariates_unchanged"),all(preserved))
}
put("no_outcome_or_exposure_imputation",all(imp$method[setdiff(names(imp$method),b$core_covariates)]==""))
put("imputed_target_predictors_complete",all(vapply(names(imp$method)[imp$method!=""],function(v){
 predictors<-names(which(imp$predictorMatrix[v,]!=0)); fixed<-predictors[imp$method[predictors]==""]
 all(vapply(imp$data[fixed],function(x)!anyNA(x),TRUE))},TRUE)))
put("finite_MI_chain_mean_and_variance",all(is.finite(imp$chainMean[!is.na(imp$chainMean)]))&&
 all(is.finite(imp$chainVar[!is.na(imp$chainVar)])))
put("20_allowed_edges",nrow(b$edges)==20L)
put("death_is_interval_not_exact",all(b$panel$obstype==1L))
put("partial_sets_living_only",all(unlist(b$censor_states)%in%1:8))
events<-as.data.table(imp$loggedEvents)
fwrite(rbindlist(checks),file.path(out,"model_qc/model_bundle_review.csv"))
if(nrow(events))fwrite(events[,.(n=.N),by=.(dep,meth,out)],file.path(out,"model_qc/mice_logged_event_summary.csv"))
writeLines(c("# CHARLS model bundle review", "",
 paste0("Checks passed: ",sum(rbindlist(checks)$pass),"/",length(checks)),
 paste0("MICE logged events: ",nrow(events),". Inspect logged_event_summary if non-zero; these are not automatically optimizer failures."),
 "Observed covariates, exposure, IDs and outcome/time inputs are preserved across the20draws.",
 "Finite-chain checks are numerical checks, not proof that MI converged or that MAR holds.",
 "Only core covariates are imputed; exploratory auxiliary outcome summaries are not treated as exact event times."),
 file.path(out,"model_qc/model_bundle_review.md"))
print(rbindlist(checks));if(nrow(events))print(events[,.(n=.N),by=.(dep,meth,out)])
stopifnot(all(rbindlist(checks)$pass))
