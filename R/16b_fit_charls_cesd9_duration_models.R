#!/usr/bin/env Rscript
# Stage 2: CORE/BOTH with a common frozen set of 20 imputations.
# Reuses the exact signed stage-1 Rubin pooling function, without sourcing its CLI.
options(stringsAsFactors=FALSE,contrasts=c("contr.treatment","contr.treatment"))
suppressPackageStartupMessages({library(data.table);library(msm);library(digest);library(jsonlite)})
setDTthreads(1L)
args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=1L,args[1]%in%c("validate","component","pool"))
root<-"outputs/charls_cesd9_duration_20260914"
old<-"outputs/charls_external_extension_20260913"
script<-"R/16b_fit_charls_cesd9_duration_models.R"
bundle_path<-file.path(root,"data/charls_stage2_model_bundle.rds")
mids_path<-file.path(root,"data/charls_stage2_mice20.rds")
sha<-function(p)digest(file=p,algo="sha256")
control<-read_json("analysis/charls_stage2_cesd9_duration_control.json",simplifyVector=TRUE)
stopifnot(control$research_stage==2L,isTRUE(control$stage_enabled),isTRUE(control$stop_after_stage),
  identical(control$automatic_next_stage,FALSE),identical(control$allow_later_stages,FALSE),
  identical(control$models,c("CORE","BOTH")),control$required_model_fits==40L,
  identical(control$model_directory,"private/model_store/CHARLS_CESD9_duration_20260914/models"))
model_root<-control$model_directory
bun<-readRDS(bundle_path)
stopifnot(identical(bun$status,"APPROVED_FOR_STAGE2_MODELS"),all(bun$preparation_checks$pass),
  identical(bun$models,c("CORE","BOTH")),length(bun$imputations)==20L)
review_path<-file.path(root,"model_qc/model_bundle_review.json")
review<-read_json(review_path,simplifyVector=TRUE)
stopifnot(identical(review$status,"PREPARATION_NUMERICAL_CHECKS_PASS"),
  review$checks_total==review$checks_passed,review$checks_total>0,
  identical(review$bundle_sha256,sha(bundle_path)),identical(review$mids_sha256,sha(mids_path)),
  identical(review$review_script_sha256,sha("R/16c_verify_charls_stage2_bundle.R")),
  identical(review$review_csv_sha256,sha(file.path(root,"model_qc/model_bundle_review.csv"))))
stopifnot(all(fread(file.path(root,"model_qc/model_bundle_review.csv"))$pass))
stopifnot(all(vapply(names(bun$input_hashes),sha,"")==bun$input_hashes))
edges<-as.data.table(bun$edges)
target<-match(c("0-2","0-4","1-3","1-5"),edges$transition)
expo<-"restless_bl"
stopifnot(nrow(edges)==20L,!anyNA(target),length(unique(target))==4L,
  nrow(bun$baseline)==12036L,nrow(bun$panel)==53786L)
bundle_sha<-sha(bundle_path);script_sha<-sha(script)
signature_for<-function(spec)digest(list(bundle=bundle_sha,script=script_sha,
  msm=as.character(packageVersion("msm")),spec=spec),algo="sha256")
save_atomic<-function(x,path) {
  if(file.exists(path))stop("Existing result retained, no overwrite: ",path)
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp<-tempfile("charls_stage2_",tmpdir=dirname(path),fileext=".rds")
  saveRDS(x,tmp,compress="gzip")
  if(!file.rename(tmp,path))stop("Publish failed; complete temporary object retained: ",tmp)
}
make_sample<-function(spec,imp) {
  covs<-bun$model_covariates[[spec]]
  b<-as.data.table(copy(bun$imputations[[imp]]))
  stopifnot(nrow(b)==12036L,!anyDuplicated(b$id),identical(b$id,bun$baseline$id))
  z<-merge(bun$panel,b[,c("id",expo,covs),with=FALSE],by="id",sort=FALSE)
  setorder(z,id,time);z[,subject_id:=.GRP,by=id]
  stopifnot(nrow(z)==53786L,uniqueN(z$id)==12036L,
    !anyNA(z[,c("id","time","state",expo,covs),with=FALSE]))
  as.data.frame(z)
}
design<-function(spec,z) {
  form<-reformulate(c(expo,bun$model_covariates[[spec]]))
  mm<-model.matrix(form,z)
  stopifnot(all(is.finite(mm)),qr(mm)$rank==ncol(mm))
  labels<-setdiff(colnames(mm),"(Intercept)")
  con<-setNames(lapply(labels,function(v)if(v==expo)edges$exposure_group else rep(1L,nrow(edges))),labels)
  list(formula=form,constraints=con,labels=labels)
}
warm_source<-function(imp) {
  p<-file.path(old,"models",sprintf("MI_%02d.rds",imp))
  stopifnot(p%in%names(bun$warm_start_hashes),identical(sha(p),bun$warm_start_hashes[[p]]))
  w<-readRDS(p)
  stopifnot(identical(w$status,"ACCEPTED"),identical(w$signature,bun$primary_stage1_signature),
    as.integer(w$imputation)==imp,identical(w$model,sprintf("MI_%02d",imp)),
    w$n_people==12036L,w$n_records==53786L)
  w
}
inspect_fit<-function(fit) {
  stopifnot(inherits(fit,"msm"),identical(as.integer(fit$opt$convergence),0L),isTRUE(fit$foundse),
    is.finite(as.numeric(logLik(fit))),is.matrix(fit$covmat),all(is.finite(fit$covmat)))
  eigenvalues<-eigen((fit$opt$hessian+t(fit$opt$hessian))/2,symmetric=TRUE,only.values=TRUE)$values
  stopifnot(all(is.finite(eigenvalues)),min(eigenvalues)>1e-8)
  k<-nrow(edges);stopifnot(as.integer(fit$qmodel$npars)==k)
  col<-match(expo,fit$qcmodel$covlabels);stopifnot(!is.na(col))
  idx<-k+(col-1L)*k+seq_len(k)
  theta<-as.numeric(fit$paramdata$params[idx]);V<-fit$covmat[idx,idx,drop=FALSE]
  stopifnot(all(is.finite(theta)),max(abs(V-t(V)))<1e-8,all(diag(V)>=0))
  V<-(V+t(V))/2
  ev<-eigen(V,symmetric=TRUE,only.values=TRUE)$values
  stopifnot(min(ev)>=-1e-8*max(1,max(abs(ev))))
  cm<-coef(fit)[[expo]]
  stopifnot(isTRUE(all.equal(theta,vapply(seq_len(k),function(j)
    cm[edges$from_index[j],edges$to_index[j]],0.0),tolerance=1e-10)))
  list(block=list(theta=theta,covariance=V),mineig=min(eigenvalues))
}
verify_saved<-function(obj,spec,imp,z,des) {
  stopifnot(identical(obj$status,"ACCEPTED"),identical(obj$signature,signature_for(spec)),
    identical(obj$spec,spec),identical(obj$model,sprintf("MI_%02d",imp)),
    as.integer(obj$imputation)==imp,obj$n_people==12036L,obj$n_records==53786L,
    identical(obj$constraints,des$constraints),identical(obj$model_input_hash,digest(z,algo="sha256")),
    identical(obj$draw_hash,digest(as.data.frame(bun$imputations[[imp]]),algo="sha256")),
    identical(obj$warm_source_sha256,unname(bun$warm_start_hashes[file.path(old,"models",sprintf("MI_%02d.rds",imp))])))
  audited<-inspect_fit(obj$fit)
  stopifnot(isTRUE(all.equal(audited$block,obj$block,tolerance=1e-12)),
    identical(as.character(obj$fit$qcmodel$covlabels),des$labels))
  invisible(TRUE)
}
fit_one<-function(spec,imp) {
  label<-sprintf("%s_MI_%02d",spec,imp)
  path<-file.path(model_root,spec,sprintf("MI_%02d.rds",imp))
  z<-make_sample(spec,imp);des<-design(spec,z)
  if(file.exists(path)) {
    obj<-readRDS(path);verify_saved(obj,spec,imp,z,des)
    message("Verified accepted cache: ",label);return(invisible(obj))
  }
  warm<-warm_source(imp);q<-warm$q_at_means;diag(q)<-0
  starts<-setNames(lapply(des$labels,function(v) {
    if(v%in%names(warm$warm_covariates))warm$warm_covariates[[v]] else rep(0,nrow(edges))
  }),des$labels)
  stopifnot(all(vapply(starts,function(x)length(x)==nrow(edges)&&all(is.finite(x)),TRUE)),
    all(bun$core_covariates[bun$core_covariates!="education4"]%in%names(warm$warm_covariates)))
  begin<-Sys.time();warnings<-character();error<-NULL;fit<-NULL;audited<-NULL
  status_path<-file.path(root,"model_qc",paste0(label,"_status.csv"))
  fwrite(data.table(model=label,status="RUNNING",pid=Sys.getpid(),started=as.character(begin),
    signature=signature_for(spec),model_file=path),status_path)
  fitargs<-list(formula=state~time,subject=quote(subject_id),data=z,qmatrix=q,
    covariates=des$formula,constraint=des$constraints,covinits=starts,obstype=z$obstype,
    center=TRUE,hessian=TRUE,control=list(fnscale=10000,maxit=1500,reltol=1e-9,trace=0))
  if(length(bun$censor_codes)) {
    fitargs$censor<-bun$censor_codes;fitargs$censor.states<-bun$censor_states
  }
  fit<-tryCatch(withCallingHandlers(do.call(msm::msm,fitargs),warning=function(w) {
    warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  }),error=function(e){error<<-conditionMessage(e);NULL})
  if(!is.null(fit))audited<-tryCatch(inspect_fit(fit),error=function(e){error<<-conditionMessage(e);NULL})
  accepted<-!is.null(audited)
  obj<-list(status=if(accepted)"ACCEPTED"else"REJECTED_REVIEW_REQUIRED",signature=signature_for(spec),
    spec=spec,model=sprintf("MI_%02d",imp),imputation=imp,fit=fit,
    block=if(accepted)audited$block else NULL,constraints=des$constraints,
    model_input_hash=digest(z,algo="sha256"),draw_hash=digest(as.data.frame(bun$imputations[[imp]]),algo="sha256"),
    warm_source_sha256=unname(bun$warm_start_hashes[file.path(old,"models",sprintf("MI_%02d.rds",imp))]),
    bundle_sha256=bundle_sha,script_sha256=script_sha,review_sha256=sha(review_path),
    warnings=unique(warnings),error=error,n_people=uniqueN(z$id),n_records=nrow(z),
    hessian_min_eigenvalue=if(accepted)audited$mineig else NA_real_,
    elapsed_seconds=as.numeric(difftime(Sys.time(),begin,units="secs")))
  if(accepted)verify_saved(obj,spec,imp,z,des)
  save_atomic(obj,path)
  fwrite(data.table(model=label,status=obj$status,pid=Sys.getpid(),n_people=obj$n_people,n_records=obj$n_records,
    elapsed_seconds=obj$elapsed_seconds,convergence=if(is.null(fit))NA_integer_ else fit$opt$convergence,
    found_se=if(is.null(fit))FALSE else fit$foundse,hessian_min_eigenvalue=obj$hessian_min_eigenvalue,
    warnings=paste(unique(warnings),collapse="; "),error=if(is.null(error))""else error,
    signature=obj$signature,model_file=path,model_sha256=sha(path)),status_path)
  message(label,": ",obj$status,"; seconds=",round(obj$elapsed_seconds,1))
  if(!accepted)stop("Rejected fit retained for review; no automatic omission/retry: ",label)
  invisible(obj)
}
pool_all<-function() {
  # Load only this function definition from the frozen, signed original source.
  # No top-level execution, pilot, old-results write, or model refit occurs.
  exprs<-parse("R/14b_fit_charls_disease_models.R")
  hits<-Filter(function(e)is.call(e)&&identical(e[[1]],as.name("<-"))&&
    identical(e[[2]],as.name("pool_models")),as.list(exprs))
  stopifnot(length(hits)==1L)
  env<-new.env(parent=environment())
  eval(hits[[1]],env)
  # All 40 must pass before any per-spec pooling begins.
  manifest<-list()
  for(spec in control$models)for(i in 1:20) {
    p<-file.path(model_root,spec,sprintf("MI_%02d.rds",i))
    stopifnot(file.exists(p));z<-make_sample(spec,i);obj<-readRDS(p)
    verify_saved(obj,spec,i,z,design(spec,z))
    manifest[[length(manifest)+1L]]<-data.table(spec=spec,imputation=i,path=p,sha256=sha(p))
    rm(obj,z)
  }
  for(spec in control$models) {
    env$out<-file.path(root,"pooled",spec);dir.create(file.path(env$out,"results"),recursive=TRUE,showWarnings=FALSE)
    dir.create(file.path(env$out,"model_qc"),recursive=TRUE,showWarnings=FALSE)
    env$model_dir<-file.path(model_root,spec);env$signature<-signature_for(spec)
    if(file.exists(file.path(env$out,"results/rubin_matrices_MI20.rds")))
      stop("Existing pooled result retained; explicit review required before repooling")
    env$pool_models()
  }
  for(name in c("all_edge_MI20","four_target_MI20","heart_stroke_contrasts_MI20","joint_tests_MI20")) {
    combined<-rbindlist(lapply(control$models,function(spec) {
      d<-fread(file.path(root,"pooled",spec,"results",paste0(name,".csv")));d[,model:=spec];d
    }),fill=TRUE)
    fwrite(combined,file.path(root,"results",paste0(name,".csv")),bom=TRUE)
  }
  a<-fread(file.path(root,"pooled/CORE/results/four_target_MI20.csv"))
  b<-fread(file.path(root,"pooled/BOTH/results/four_target_MI20.csv"))
  stopifnot(identical(a$transition,b$transition))
  fwrite(data.table(transition=a$transition,CORE_logHR=a$log_effect,BOTH_logHR=b$log_effect,
    change_logHR=b$log_effect-a$log_effect,CORE_HR=a$HR,BOTH_HR=b$HR,
    interpretation="Descriptive correlated-model change; no mediation proportion or change P/CI"),
    file.path(root,"results/matched_adjustment_change_descriptive.csv"),bom=TRUE)
  fwrite(rbindlist(manifest),file.path(root,"model_qc/accepted_40_model_manifest.csv"),bom=TRUE)
  write_json(list(status="POOLED_AWAITING_INDEPENDENT_REVIEW_AND_FIGURES",models=control$models,
    accepted_fits=40L,m=20L,automatic_next_stage=FALSE,primary_model="BOTH"),
    file.path(root,"model_qc/pooled_status.json"),pretty=TRUE,auto_unbox=TRUE)
}
if(args[1]=="validate") {
  stopifnot(length(args)==1L)
  for(spec in control$models)for(i in 1:20) {
    z<-make_sample(spec,i);d<-design(spec,z)
    stopifnot(length(d$constraints[[expo]])==20L)
  }
  for(i in 1:20)warm_source(i)
  write_json(list(status="FIT_PREFLIGHT_PASS",required_fits=40L,model_directory=model_root,
    bundle_sha256=bundle_sha,script_sha256=script_sha,no_models_fitted=TRUE),
    file.path(root,"model_qc/fit_preflight.json"),pretty=TRUE,auto_unbox=TRUE)
  message("FIT_PREFLIGHT_PASS: 40 designs; unchanged shared sample; frozen warm starts; D-drive output.")
} else if(args[1]=="component") {
  stopifnot(length(args)==3L,args[2]%in%control$models)
  imp<-as.integer(args[3]);stopifnot(!is.na(imp),imp%in%1:20)
  fit_one(args[2],imp)
} else {stopifnot(length(args)==1L);pool_all()}
