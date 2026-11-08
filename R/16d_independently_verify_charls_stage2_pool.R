#!/usr/bin/env Rscript
# Independent audit of the completed stage-2 CORE/BOTH pool.
# Reads saved fits and published tables only; never refits or edits model objects.
options(stringsAsFactors=FALSE, contrasts=c("contr.treatment","contr.treatment"))
suppressPackageStartupMessages({library(data.table);library(msm);library(digest);library(jsonlite)})
setDTthreads(1L)
root <- "outputs/charls_cesd9_duration_20260914"
model_root <- "private/model_store/CHARLS_CESD9_duration_20260914/models"
qc <- file.path(root,"model_qc")
bundle_path <- file.path(root,"data/charls_stage2_model_bundle.rds")
fit_script <- "R/16b_fit_charls_cesd9_duration_models.R"
sha <- function(p) digest(file=p,algo="sha256")
same_number <- function(a,b,tolerance=1e-8) {
  if(length(a)!=length(b)||!identical(dim(a),dim(b)))return(FALSE)
  a<-as.numeric(a);b<-as.numeric(b)
  if(!identical(is.na(a),is.na(b)))return(FALSE)
  keep<-!is.na(a);a<-a[keep];b<-b[keep]
  if(!identical(is.infinite(a),is.infinite(b)))return(FALSE)
  inf<-is.infinite(a);if(any(a[inf]!=b[inf]))return(FALSE)
  a<-a[!inf];b<-b[!inf]
  all(abs(a-b)<=tolerance*pmax(1,abs(a),abs(b)))
}
stop_if <- function(x,msg){if(!isTRUE(x))stop(msg,call.=FALSE)}
extract_fit <- function(fit,edges,exposure="restless_bl") {
  stop_if(inherits(fit,"msm"),"saved fit is not an msm object")
  p<-fit$paramdata
  stop_if(length(p$fixedpars)==0L&&length(p$auxpars)==0L&&all(p$constr>0),"unsupported constraint schema")
  stop_if(identical(as.integer(fit$opt$convergence),0L)&&isTRUE(fit$foundse),"fit gate failed")
  stop_if(is.finite(as.numeric(logLik(fit))),"non-finite likelihood")
  h<-as.matrix(fit$opt$hessian)
  stop_if(all(is.finite(h))&&same_number(h,t(h),1e-10),"invalid Hessian")
  h<-(h+t(h))/2
  ev<-eigen(h,symmetric=TRUE,only.values=TRUE)$values
  stop_if(min(ev)>1e-8,"Hessian not positive definite")
  stop_if(nrow(h)==length(p$optpars)&&nrow(h)==length(fit$opt$par),"Hessian dimensions mismatch")
  gi<-match(p$constr,p$constr[p$optpars])
  stop_if(!anyNA(gi)&&!anyDuplicated(p$constr[p$optpars]),"constraint map failed")
  expand<-matrix(0,length(p$constr),length(p$optpars));expand[cbind(seq_along(gi),gi)]<-1
  theta_all<-drop(expand%*%as.numeric(fit$opt$par))
  cov_all<-expand%*%(2*chol2inv(chol(h)))%*%t(expand)
  stop_if(same_number(theta_all,as.numeric(p$params),1e-9),"expanded coefficients disagree")
  stop_if(same_number(cov_all,fit$covmat,2e-7),"expanded covariance disagrees")
  q<-fit$qmodel$imatrix;cells<-which(q!=0,arr.ind=TRUE);cells<-cells[order(cells[,1L],cells[,2L]),,drop=FALSE]
  stop_if(nrow(cells)==nrow(edges)&&identical(as.integer(cells[,1L]),as.integer(edges$from_index))&&
    identical(as.integer(cells[,2L]),as.integer(edges$to_index)),"edge ordering differs")
  labs<-fit$qcmodel$covlabels;col<-match(exposure,labs);stop_if(!is.na(col),"exposure not found")
  qidx<-which(p$plabs=="qcov")
  stop_if(length(qidx)==nrow(edges)*length(labs),"unexpected qcov layout")
  idx<-matrix(qidx,nrow=nrow(edges))[,col]
  coefficient_matrix<-coef(fit)[[exposure]]
  th<-as.numeric(theta_all[idx]);V<-cov_all[idx,idx,drop=FALSE]
  stop_if(same_number(th,as.numeric(coefficient_matrix[cells]),1e-10),"coef matrix disagrees")
  stop_if(all(is.finite(V))&&all(diag(V)>=0),"exposure covariance invalid")
  V<-(V+t(V))/2;ve<-eigen(V,symmetric=TRUE,only.values=TRUE)$values
  stop_if(min(ve)>=-1e-8*max(1,abs(ve)),"exposure covariance not PSD")
  list(theta=th,covariance=V,hessian_min=min(ev),likelihood=as.numeric(logLik(fit)))
}
pool <- function(theta,covs) {
  m<-ncol(theta);k<-nrow(theta);qbar<-rowMeans(theta)
  B<-tcrossprod(sweep(theta,1L,qbar,"-"))/(m-1L)
  U<-Reduce(`+`,lapply(covs,function(x)x/m));T<-U+(m+1)/m*B
  list(m=m,Q=theta,qbar=qbar,U=U,B=B,T=T)
}
scalar <- function(pl,L) {
  est<-as.numeric(L%*%pl$qbar);u<-diag(L%*%pl$U%*%t(L));b<-diag(L%*%pl$B%*%t(L));tt<-diag(L%*%pl$T%*%t(L))
  stop_if(all(u>0)&&all(tt>0)&&all(b>=-1e-12),"invalid Rubin scalar variance")
  df<-rep(Inf,length(est));nz<-b>1e-15;r<-(pl$m+1)/pl$m*b[nz]/u[nz];df[nz]<-(pl$m-1)*((r+1)/r)^2
  se<-sqrt(tt);hw<-qt(.975,df)*se
  data.table(log_effect=est,se=se,df=df,HR=exp(est),lower=exp(est-hw),upper=exp(est+hw),p=2*pt(abs(est/se),df,lower.tail=FALSE))
}
joint <- function(pl,L) {
  V<-L%*%pl$T%*%t(L);V<-(V+t(V))/2;ev<-eigen(V,symmetric=TRUE,only.values=TRUE)$values
  stop_if(all(is.finite(ev))&&min(ev)>1e-12*max(ev)&&qr(V)$rank==nrow(L),"invalid joint covariance")
  a<-as.numeric(L%*%pl$qbar);data.table(chi_square=drop(t(a)%*%solve(V,a)),df=nrow(L),p=pchisq(drop(t(a)%*%solve(V,a)),nrow(L),lower.tail=FALSE))
}
checks<-list();add<-function(check,pass,detail=""){checks[[length(checks)+1L]]<<-data.table(check=check,pass=isTRUE(pass),detail=detail);if(!isTRUE(pass))stop(paste(check,detail,sep=": "),call.=FALSE)}
run<-function(){
  bun<-readRDS(bundle_path);control<-fromJSON("analysis/charls_stage2_cesd9_duration_control.json",simplifyVector=TRUE)
  add("stage2_control_fixed",identical(as.character(control$models),c("CORE","BOTH"))&&control$required_model_fits==40L)
  add("bundle_status",identical(bun$status,"APPROVED_FOR_STAGE2_MODELS")&&all(bun$preparation_checks$pass))
  add("bundle_review_credential",file.exists(file.path(qc,"model_bundle_review.json"))&&
    identical(fromJSON(file.path(qc,"model_bundle_review.json"),simplifyVector=TRUE)$status,"PREPARATION_NUMERICAL_CHECKS_PASS"))
  edges<-as.data.table(bun$edges);add("frozen_9state_20edge_contract",nrow(edges)==20L)
  targets<-match(c("0-2","0-4","1-3","1-5"),edges$transition);add("target_order",!anyNA(targets))
  review<-list();all_pooled<-list()
  for(spec in control$models){
    sig<-digest(list(bundle=sha(bundle_path),script=sha(fit_script),msm=as.character(packageVersion("msm")),spec=spec),algo="sha256")
    th<-matrix(NA_real_,20L,20L);covs<-vector("list",20L)
    for(i in 1:20){
      p<-file.path(model_root,spec,sprintf("MI_%02d.rds",i));add(p,file.exists(p));obj<-readRDS(p)
      add(paste0(spec,"_",i,"_identity"),identical(obj$status,"ACCEPTED")&&identical(obj$signature,sig)&&as.integer(obj$imputation)==i)
      rec<-extract_fit(obj$fit,edges);add(paste0(spec,"_",i,"_saved_block"),same_number(rec$theta,obj$block$theta,1e-9)&&same_number(rec$covariance,obj$block$covariance,2e-7))
      add(paste0(spec,"_",i,"_sample"),obj$n_people==12036L&&obj$n_records==53786L)
      th[,i]<-rec$theta;covs[[i]]<-rec$covariance
    }
    pl<-pool(th,covs);all_pooled[[spec]]<-pl
    published<-readRDS(file.path(root,"pooled",spec,"results/rubin_matrices_MI20.rds"))
    for(f in c("Q","qbar","U","B","T"))add(paste0(spec,"_Rubin_",f),same_number(pl[[f]],published[[f]],2e-7))
    all_tab<-fread(file.path(root,"pooled",spec,"results/all_edge_MI20.csv"));target_tab<-fread(file.path(root,"pooled",spec,"results/four_target_MI20.csv"));contrast_tab<-fread(file.path(root,"pooled",spec,"results/heart_stroke_contrasts_MI20.csv"));joint_tab<-fread(file.path(root,"pooled",spec,"results/joint_tests_MI20.csv"))
    ex_all<-scalar(pl,diag(nrow(edges)));ex_target<-ex_all[targets];ex_target[,p_BH_4:=p.adjust(p,"BH",n=4L)];C<-matrix(0,2,nrow(edges));C[1,targets[1:2]]<-c(1,-1);C[2,targets[3:4]]<-c(1,-1);ex_con<-scalar(pl,C);ex_con[,p_BH_2:=p.adjust(p,"BH",n=2L)];ex_joint<-rbind(joint(pl,diag(nrow(edges))[targets,,drop=FALSE]),joint(pl,C))
    for(j in names(ex_all))add(paste0(spec,"_all_",j),same_number(all_tab[[j]],ex_all[[j]],2e-7))
    for(j in names(ex_target))add(paste0(spec,"_target_",j),same_number(target_tab[[j]],ex_target[[j]],2e-7))
    for(j in names(ex_con))add(paste0(spec,"_contrast_",j),same_number(contrast_tab[[j]],ex_con[[j]],2e-7))
    for(j in names(ex_joint))add(paste0(spec,"_joint_",j),same_number(joint_tab[[j]],ex_joint[[j]],2e-7))
    review[[spec]]<-data.table(model=spec,m=20L,targets_global_p=ex_joint$p[1],contrast_joint_p=ex_joint$p[2])
  }
  fwrite(rbindlist(checks),file.path(qc,"stage2_independent_pool_review.csv"),bom=TRUE)
  fwrite(rbindlist(review),file.path(qc,"stage2_independent_pool_summary.csv"),bom=TRUE)
  write_json(list(status="COMPUTATIONAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_REVIEW",checks_total=length(checks),checks_passed=sum(vapply(checks,function(x)x$pass,TRUE)),models=as.character(control$models),accepted_fits=40L,independent_reviewer_script_sha256=sha("R/16d_independently_verify_charls_stage2_pool.R"),bundle_sha256=sha(bundle_path)),file.path(qc,"stage2_independent_pool_review.json"),pretty=TRUE,auto_unbox=TRUE)
  writeLines(c("# CHARLS stage-2 independent pool review","",paste0("Status: COMPUTATIONAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_REVIEW"),paste0("Checks passed: ",sum(vapply(checks,function(x)x$pass,TRUE)),"/",length(checks)),"","This review reconstructed the constrained exposure coefficients and Hessian covariance from all 40 saved fits, independently recomputed Rubin U/B/T, scalar intervals, BH corrections, direct heart-stroke contrasts and joint Wald tests, and compared them with the published CORE/BOTH tables.","","This is a computational verification, not proof of MAR, MICE convergence, Markov/time-homogeneity, no unmeasured confounding, or causal mediation."),file.path(qc,"stage2_independent_pool_review.md"),useBytes=TRUE)
  message("COMPUTATIONAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_REVIEW: ",sum(vapply(checks,function(x)x$pass,TRUE)),"/",length(checks))
}
tryCatch(run(),error=function(e){fwrite(rbindlist(checks),file.path(qc,"stage2_independent_pool_review.csv"),bom=TRUE);message("FAILED_REVIEW_REQUIRED: ",conditionMessage(e));quit(save="no",status=1L)})
