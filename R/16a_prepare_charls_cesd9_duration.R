#!/usr/bin/env Rscript
# Research stage 2 only: shared fresh MI, no disease-transition model fitting.
options(stringsAsFactors=FALSE, contrasts=c("contr.treatment","contr.treatment"), warn=1)
suppressPackageStartupMessages({library(data.table); library(mice); library(digest); library(jsonlite)})
setDTthreads(1L)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L, args %in% c("prepare","resume"))
resume <- args=="resume"
out <- "outputs/charls_cesd9_duration_20260914"
old <- "outputs/charls_external_extension_20260913"
for (folder in c("data","model_qc","tables","models","results","reports","logs"))
  dir.create(file.path(out,folder),recursive=TRUE,showWarnings=FALSE)
bundle_path <- file.path(out,"data/charls_stage2_model_bundle.rds")
mids_path <- file.path(out,"data/charls_stage2_mice20.rds")
if (file.exists(bundle_path) || (!resume && file.exists(mids_path)))
  stop("Existing preparation retained. Review or reuse it; never overwrite by rerunning preparation.")
if(resume && !file.exists(mids_path))stop("Resume requires the saved completed MICE object.")
checks <- list()
sha <- function(p) digest(file=p,algo="sha256")
utc <- function() format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")
putcheck <- function(name,pass,detail="") {
  checks[[length(checks)+1L]] <<- data.table(check=name,pass=isTRUE(pass),detail=as.character(detail))
  if (!isTRUE(pass)) stop(paste(name,detail,sep=": "),call.=FALSE)
}
json <- function(x,p) write_json(x,p,pretty=TRUE,auto_unbox=TRUE,digits=NA,null="null")
state <- function(status,detail="") {
  json(list(status=status,started_or_updated_utc=utc(),pid=Sys.getpid(),stage=2L,
    detail=detail,model_fits_started=FALSE,next_stage_authorized=FALSE),
    file.path(out,"model_qc/preparation_status.json"))
}
save_new <- function(x,p) {
  if(file.exists(p))stop("Refusing to overwrite ",p)
  tmp<-tempfile("stage2_",tmpdir=dirname(p),fileext=".rds")
  saveRDS(x,tmp,compress="gzip")
  if(!file.rename(tmp,p))stop("Cannot atomically publish ",p,"; temporary file retained: ",tmp)
}
same_observed <- function(a,b) {
  known<-!is.na(a)
  identical(is.na(a[known]),is.na(b[known])) &&
    identical(as.character(a[known]),as.character(b[known]))
}
run <- function() {
  state("AUDITING_BASELINE_INPUTS")
  control_file<-"analysis/charls_stage2_cesd9_duration_control.json"
  contract_file<-"analysis/charls_stage2_cesd9_duration_contract.md"
  old_bundle<-file.path(old,"data/charls_model_bundle.rds")
  old_mids<-file.path(old,"data/charls_baseline_mice20.rds")
  observed_file<-file.path(old,"data/charls_observed_exposures_tristate_panel.rds")
  raw_file<-"private/stage/charls_field_audit/w1_health.dta"
  teaching_file<-"private/source/reference_scripts/msm_example.R"
  sources<-c(control_file,contract_file,old_bundle,old_mids,observed_file,raw_file,
    file.path(old,"stage1_delivery/delivery_review.json"),
    "R/13a_prepare_charls_external_extension.R","R/14a_prepare_charls_model_bundle.R",
    "R/14b_fit_charls_disease_models.R","R/16a_prepare_charls_cesd9_duration.R",teaching_file)
  if(resume)sources<-c(sources,mids_path,file.path(out,"data/charls_stage2_preimputation.rds"),
    file.path(out,"model_qc/preparation_recovery_20260914/16a_original_before_resume.R"))
  putcheck("all_sources_exist",all(file.exists(sources)))
  source_hashes<-setNames(vapply(sources,sha,""),sources)
  control<-read_json(control_file,simplifyVector=TRUE)
  putcheck("stage2_scope_check",isTRUE(control$stage_enabled) &&
    control$research_stage==2 && isTRUE(control$stop_after_stage) &&
    identical(control$automatic_next_stage,FALSE) &&
    identical(control$allow_later_stages,FALSE) &&
    identical(control$models,c("CORE","BOTH")) && control$required_model_fits==40L)
  prior_review<-read_json(file.path(old,"stage1_delivery/delivery_review.json"),simplifyVector=TRUE)
  putcheck("stage1_review_package_exists",identical(prior_review$status,"STAGE1_REVIEW_PACKAGE_READY") &&
    prior_review$completed_imputations==20L && identical(prior_review$signature,control$primary_stage1_signature))
  for(p in names(prior_review$source_hashes))
    putcheck(paste0("stage1_source_preserved:",p),identical(sha(p),prior_review$source_hashes[[p]]))
  bun<-readRDS(old_bundle); old_imp<-readRDS(old_mids)
  base<-as.data.table(copy(bun$baseline)); panel<-copy(bun$panel)
  core<-bun$core_covariates
  putcheck("fixed_sample_and_panel",nrow(base)==12036L && nrow(panel)==53786L &&
    !anyDuplicated(base$id) && setequal(base$id,panel$id) && all(panel$obstype==1L))
  putcheck("approved_old_bundle",identical(bun$status,"APPROVED_FOR_PILOT_AND_MI") &&
    all(bun$preparation_checks$pass))
  for(p in names(bun$input_hashes))putcheck(paste0("old_bundle_input:",p),identical(sha(p),bun$input_hashes[[p]]))
  putcheck("mids_original_rows_match",nrow(old_imp$data)==nrow(base) &&
    identical(as.numeric(old_imp$data$restless_bl),base$restless_bl))
  for(v in core)putcheck(paste0("old_mids_original_core:",v),
    identical(is.na(old_imp$data[[v]]),is.na(base[[v]])) && same_observed(base[[v]],old_imp$data[[v]]))
  d<-as.data.table(readRDS(observed_file))
  w<-d[wave==1][match(base$id,id)]
  putcheck("W1_identity",identical(w$id,base$id) && all(w$wave==1L) &&
    identical(w$restless_score_observed,base$restless_bl))
  vars<-c("botherl","mindtsl","depresl","effortl","fhopel","fearll","sleeprl","whappyl","flonel","goingl")
  raw_fields<-sprintf("dc%03d",9:18)
  putcheck("raw_W1_hash",identical(sha(raw_file),"87bae0f7e75899ba3c844985aceaed42dd9e60483edbebc9c402d3cb866a810e"))
  raw<-as.data.table(haven::read_dta(raw_file,col_select=c("ID",raw_fields,"da049")))
  raw[,id:=as.character(ID)]
  raw[,id:=paste0(substr(id,1,9),"0",substr(id,10,11))]
  putcheck("raw_id_unique",!anyDuplicated(raw$id) && all(base$id %in% raw$id))
  raw<-raw[match(base$id,id)]
  valid4<-function(x){z<-as.numeric(x);z[!z %in% 1:4]<-NA_real_;z}
  for(j in seq_along(vars)) {
    a<-valid4(raw[[raw_fields[j]]]); b<-valid4(w[[vars[j]]])
    putcheck(paste0("raw_question_matches_package:",vars[j]),identical(a,b))
  }
  keep<-setdiff(vars,"sleeprl")
  item_names<-paste0("cesd_",keep)
  raw_items<-as.data.frame(lapply(w[,..keep],valid4))
  for(v in keep)raw_items[[v]]<-if(v %in% c("fhopel","whappyl"))4-raw_items[[v]] else raw_items[[v]]-1
  names(raw_items)<-item_names
  sums<-rowSums(as.matrix(raw_items))
  putcheck("CESD9_observed_sum_reproduced",identical(as.numeric(sums),as.numeric(w$cesd9_no_sleep)))
  putcheck("CESD9_missing_count",sum(is.na(sums))==597L && sum(!is.na(sums))==11439L)
  night<-as.numeric(raw$da049)
  night[!is.finite(night)|night<1|night>18]<-NA_real_
  putcheck("night_raw_matches_existing_observed",identical(night,as.numeric(w$sleep_hours_observed)))
  putcheck("night_missing_count",sum(is.na(night))==114L)
  putcheck("joint_complete_count",sum(!is.na(sums)&!is.na(night))==11362L)
  knots<-unname(quantile(night,c(1/3,2/3),na.rm=TRUE,type=7))
  boundary<-range(night,na.rm=TRUE)
  putcheck("frozen_spline_knots",identical(knots,c(6,7)) && identical(boundary,c(1,13)))
  derive<-function(z) {
    total<-rowSums(as.matrix(z[,item_names,drop=FALSE]))
    ns<-splines::ns(z$night_hours,knots=knots,Boundary.knots=boundary,intercept=FALSE)
    ans<-as.data.table(z[,c("restless_bl",core),drop=FALSE])
    for(v in setdiff(core,"education4"))if(is.factor(ans[[v]]))ans[,(v):=as.numeric(as.character(get(v)))]
    ans[,education4:=factor(education4,levels=1:4)]
    ans[,cesd9_5:=as.numeric(total)/5]
    for(j in 1:3)ans[,(paste0("night_ns",j)):=as.numeric(ns[,j])]
    ans[,id:=base$id]
    setcolorder(ans,c("id","restless_bl",core,"cesd9_5",paste0("night_ns",1:3)))
    ans
  }
  mi_data<-old_imp$data
  for(v in item_names)mi_data[[v]]<-raw_items[[v]]
  mi_data$night_hours<-night
  extra<-c(item_names,"night_hours")
  targets<-c(core,extra)
  methods<-setNames(rep("",ncol(mi_data)),names(mi_data))
  methods[core]<-old_imp$method[core]
  for(v in extra)if(anyNA(mi_data[[v]]))methods[v]<-"pmm"
  pred<-matrix(1L,ncol(mi_data),ncol(mi_data),dimnames=list(names(mi_data),names(mi_data)))
  diag(pred)<-0L;pred[methods=="",]<-0L
  fixed<-setdiff(names(mi_data),targets)
  putcheck("fixed_fields_have_no_missing",!anyNA(mi_data[fixed]))
  putcheck("no_CESD_total_in_predictor_matrix",!any(c("cesd9_5","cesd9_total",paste0("night_ns",1:3)) %in% names(mi_data)))
  miss<-rbindlist(lapply(names(mi_data),function(v)data.table(variable=v,n=nrow(mi_data),
    missing=sum(is.na(mi_data[[v]])),method=methods[[v]],
    role=if(v %in% extra)"added_adjustment_item_or_hours" else if(v %in% core)"core_covariate" else "fixed_exposure_or_outcome_summary")))
  fwrite(miss,file.path(out,"tables/imputation_missingness.csv"),bom=TRUE)
  fwrite(data.table(n_people=12036L,n_records=53786L,CESD9_observed=sum(!is.na(sums)),
    night_observed=sum(!is.na(night)),both_observed=sum(!is.na(sums)&!is.na(night)),
    retained_for_joint_MI=12036L),file.path(out,"tables/stage2_sample_coverage.csv"),bom=TRUE)
  if(!resume)save_new(list(data=mi_data,ids=base$id,raw_item_names=item_names,
    method=methods,predictor_matrix=pred,spline_spec=list(knots=knots,boundary=boundary,df=3L,quantile_type=7L)),
    file.path(out,"data/charls_stage2_preimputation.rds"))
  mi_warnings<-character()
  if(resume) {
    state("RESUMING_COMPLETED_MICE_REVIEW","No MICE refit; correcting data.table memory-reference comparison only.")
    origin<-readRDS(file.path(out,"data/charls_stage2_preimputation.rds"))
    imp<-readRDS(mids_path)
    putcheck("resume_original_data_identical",identical(origin$data,mi_data) && identical(imp$data,mi_data))
    putcheck("resume_methods_predictors_identical",identical(origin$method,methods) &&
      identical(origin$predictor_matrix,pred) && identical(imp$method,methods) &&
      identical(dimnames(imp$predictorMatrix),dimnames(pred)) &&
      identical(dim(imp$predictorMatrix),dim(pred)) &&
      identical(as.numeric(imp$predictorMatrix),as.numeric(pred)))
    old_runtime<-read_json(file.path(out,"model_qc/imputation_runtime_review.json"),simplifyVector=TRUE)
    mi_warnings<-unlist(old_runtime$warnings,use.names=FALSE)
    message("Reusing completed joint MICE20; no imputation rerun.")
  } else {
    state("MICE20_RUNNING","Fresh joint original data; no reuse of completed stage1 draws; no model fitting.")
    message("MICE20 started: n=",nrow(mi_data),"; targets=",sum(methods!=""),"; total iterations=20.")
    imp<-withCallingHandlers(mice(mi_data,m=20,maxit=20,method=methods,predictorMatrix=pred,
      seed=20260914,printFlag=FALSE),warning=function(w){mi_warnings<<-c(mi_warnings,conditionMessage(w));invokeRestart("muffleWarning")})
    save_new(imp,mids_path)
  }
  events<-as.data.table(imp$loggedEvents)
  if(!ncol(events))events<-data.table(it=integer(),im=integer(),dep=character(),meth=character(),out=character())
  fwrite(events,file.path(out,"model_qc/mice_logged_events.csv"),bom=TRUE)
  json(list(warnings=as.list(unique(mi_warnings)),logged_events=nrow(events),
    m=imp$m,iteration=imp$iteration,methods=as.list(imp$method)),
    file.path(out,"model_qc/imputation_runtime_review.json"))
  putcheck("mice_no_logged_events",nrow(events)==0L,"Nonzero events require review before models.")
  putcheck("mice_no_runtime_warnings",length(mi_warnings)==0L,paste(unique(mi_warnings),collapse="; "))
  putcheck("20_draws_20_iterations",imp$m==20L&&imp$iteration==20L)
  draws<-vector("list",20L)
  specs<-list(CORE=core,CESD9=c(core,"cesd9_5"),
    DURATION=c(core,paste0("night_ns",1:3)),BOTH=c(core,"cesd9_5",paste0("night_ns",1:3)))
  specs<-specs[control$models]
  for(i in 1:20) {
    z<-complete(imp,i)
    putcheck(paste0("draw",i,"_complete"),!anyNA(z))
    putcheck(paste0("draw",i,"_all_observed_preserved"),all(vapply(names(mi_data),function(v)same_observed(mi_data[[v]],z[[v]]),TRUE)))
    putcheck(paste0("draw",i,"_fixed_fields_preserved"),all(vapply(fixed,function(v)identical(mi_data[[v]],z[[v]]),TRUE)))
    putcheck(paste0("draw",i,"_nine_item_support"),all(as.matrix(z[,item_names]) %in% 0:3))
    for(v in extra)putcheck(paste0("draw",i,"_donor_support:",v),
      all(z[[v]][is.na(mi_data[[v]])] %in% mi_data[[v]][!is.na(mi_data[[v]])]))
    draws[[i]]<-derive(z)
    putcheck(paste0("draw",i,"_CESD9_range"),all(draws[[i]]$cesd9_5>=0 & draws[[i]]$cesd9_5<=27/5))
    for(s in names(specs)){
      mm<-model.matrix(reformulate(c("restless_bl",specs[[s]])),data=draws[[i]])
      putcheck(paste0("draw",i,"_",s,"_design_full_rank"),nrow(mm)==12036L && qr(mm)$rank==ncol(mm) && all(is.finite(mm)))
    }
  }
  observed_base<-derive(mi_data)
  putcheck("baseline_CESD_missing_preserved",identical(is.na(observed_base$cesd9_5),is.na(sums)))
  putcheck("panel_unchanged",identical(as.data.frame(panel),as.data.frame(bun$panel)),
    "Compare all columns, types, values and row order; data.table internal memory references are not data.")
  warm_paths<-file.path(old,"models",sprintf("MI_%02d.rds",1:20))
  putcheck("all_stage1_warm_files_exist",all(file.exists(warm_paths)))
  warm_hashes<-setNames(vapply(warm_paths,sha,""),warm_paths)
  for(p in names(source_hashes))putcheck(paste0("input_unchanged:",p),identical(sha(p),source_hashes[[p]]))
  save_new(list(method=imp$method,predictor_matrix=imp$predictorMatrix,
    chainMean=imp$chainMean,chainVar=imp$chainVar,warnings=unique(mi_warnings)),
    file.path(out,"model_qc/mice_diagnostics.rds"))
  preparation_checks<-rbindlist(checks)
  bundle<-list(status="APPROVED_FOR_STAGE2_MODELS",research_stage=2L,panel=panel,
    baseline=observed_base,imputations=draws,edges=copy(bun$edges),q_initial=bun$q_initial,
    core_covariates=core,model_covariates=specs,censor_codes=bun$censor_codes,censor_states=bun$censor_states,
    original_mi_data=mi_data,raw_item_names=item_names,imputation_target_variables=targets,
    spline_spec=list(knots=knots,boundary=boundary,df=3L,quantile_type=7L),
    input_hashes=source_hashes,warm_start_hashes=warm_hashes,
    primary_stage1_signature=control$primary_stage1_signature,preparation_checks=preparation_checks,
    estimand="W1 restless-sleep frequency HR after CESD9 and/or natural-spline nightly duration adjustment",
    m=20L,seed=20260914L,models=names(specs))
  save_new(bundle,bundle_path)
  fwrite(preparation_checks,file.path(out,"model_qc/preparation_checks.csv"),bom=TRUE)
  state("PREPARED_AWAITING_INDEPENDENT_REVIEW","CORE and BOTH share completed MI20; 40 fits specified; no fit has started.")
  message("PREPARED_AWAITING_INDEPENDENT_REVIEW: all20 shared draws ready; n=12036; no models fitted.")
  0L
}
exit_code<-tryCatch(run(),error=function(e){
  if(length(checks))fwrite(rbindlist(checks),file.path(out,"model_qc/preparation_checks.csv"),bom=TRUE)
  state("FAILED_PREPARATION_REVIEW_REQUIRED",conditionMessage(e))
  message("FAILED_PREPARATION_REVIEW_REQUIRED: ",conditionMessage(e));1L
})
quit(save="no",status=exit_code)
