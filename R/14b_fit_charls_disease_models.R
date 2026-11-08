#!/usr/bin/env Rscript
# Continuous-time panel likelihood; execution requires a signed preparation bundle.
# Extends GlobalAgingR teaching msm workflow with structured Q, constraints and MI.
options(stringsAsFactors=FALSE, contrasts=c("contr.treatment","contr.treatment"))
suppressPackageStartupMessages({library(data.table);library(msm);library(digest)})
setDTthreads(1)
args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args))args[1] else "pilot"
component <- if(length(args)>1)as.integer(args[2]) else NA_integer_
out <- "outputs/charls_external_extension_20260913"
model_dir <- file.path(out,"models")
for (sub in c("models","model_qc","results","logs"))dir.create(file.path(out,sub),recursive=TRUE,showWarnings=FALSE)
bundle_path <- file.path(out,"data/charls_model_bundle.rds")
stopifnot(file.exists(bundle_path))
bun <- readRDS(bundle_path)
stopifnot(identical(bun$status,"APPROVED_FOR_PILOT_AND_MI"),all(bun$preparation_checks$pass))
bundle_review<-fread(file.path(out,"model_qc/model_bundle_review.csv"))
stopifnot(all(bundle_review$pass))
sha <- function(f)digest(file=f,algo="sha256")
stopifnot(all(vapply(names(bun$input_hashes),sha,"")==bun$input_hashes))
signature <- digest(list(bundle=sha(bundle_path),script=sha("R/14b_fit_charls_disease_models.R"),
                         msm=as.character(packageVersion("msm"))),algo="sha256")
edges <- as.data.table(bun$edges)
target <- match(c("0-2","0-4","1-3","1-5"),edges$transition)
stopifnot(!anyNA(target),length(unique(target))==4)
expo <- "restless_bl"
covs <- bun$core_covariates
form <- reformulate(c(expo,covs))
save_atomic <- function(x,path) {
  tmp <- tempfile("charls_",tmpdir=dirname(path),fileext=".rds")
  saveRDS(x,tmp,compress="gzip")
  if(!file.copy(tmp,path,overwrite=TRUE))stop("Cannot publish ",path)
  unlink(tmp)
}
make_sample <- function(imp=0L) {
  b <- if(imp==0L)copy(bun$baseline) else as.data.table(bun$imputations[[imp]])
  if(imp==0L)b<-b[complete.cases(b[,c("id",expo,covs),with=FALSE])]
  z<-merge(bun$panel,b[,c("id",expo,covs),with=FALSE],by="id",sort=FALSE)
  stopifnot(!anyNA(z[,c("id","time","state",expo,covs),with=FALSE]))
  setorder(z,id,time)
  z[,subject_id:=.GRP,by=id]
  as.data.frame(z)
}
fit_one <- function(imp) {
  label <- if(imp==0L)"pilot_CC" else sprintf("MI_%02d",imp)
  path <- file.path(model_dir,paste0(label,".rds"))
  if(file.exists(path)) {
    saved<-readRDS(path)
    if(!identical(saved$signature,signature))stop("Existing fit signature differs: ",label)
    if(identical(saved$status,"ACCEPTED")){message("Accepted cached ",label);return(invisible(saved))}
    stop("Prior rejected fit retained; review before explicit retry: ",label)
  }
  z<-make_sample(imp)
  mm<-model.matrix(form,z)
  covlabels<-setdiff(colnames(mm),"(Intercept)")
  constraint<-setNames(lapply(covlabels,function(x)if(x==expo)edges$exposure_group else rep(1L,nrow(edges))),covlabels)
  q<-bun$q_initial
  starts<-NULL
  pilot_path<-file.path(model_dir,"pilot_CC.rds")
  if(imp>0L&&file.exists(pilot_path)) {
    pilot<-readRDS(pilot_path)
    stopifnot(identical(pilot$status,"ACCEPTED"),identical(pilot$signature,signature))
    q<-pilot$q_at_means
    diag(q)<-0
    starts<-pilot$warm_covariates[covlabels]
  }
  begin<-Sys.time(); warnings<-character(); error<-NULL
  fwrite(data.table(model=label,status="RUNNING",pid=Sys.getpid(),started=as.character(begin),signature=signature),
    file.path(out,"model_qc",paste0(label,"_status.csv")))
  fitargs<-list(formula=state~time,subject=quote(subject_id),data=z,qmatrix=q,
    covariates=form,constraint=constraint,covinits=starts,obstype=z$obstype,center=TRUE,
    hessian=TRUE,control=list(fnscale=10000,maxit=1500,reltol=1e-9,trace=0))
  if(length(bun$censor_codes)){
    fitargs$censor<-bun$censor_codes;fitargs$censor.states<-bun$censor_states
  }
  fit<-tryCatch(withCallingHandlers(do.call(msm::msm,fitargs),warning=function(w){
    warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  }),error=function(e){error<<-conditionMessage(e);NULL})
  accepted<-FALSE; block<-NULL; qhat<-NULL; warm<-NULL; mineig<-NA_real_; covariance_ok<-FALSE
  if(!is.null(fit)) {
    eig<-tryCatch(eigen((fit$opt$hessian+t(fit$opt$hessian))/2,symmetric=TRUE,only.values=TRUE)$values,error=function(e)NA_real_)
    mineig<-min(eig)
    covariance_ok<-is.matrix(fit$covmat)&&all(is.finite(fit$covmat))
    accepted<-identical(as.integer(fit$opt$convergence),0L)&&isTRUE(fit$foundse)&&covariance_ok&&
      is.finite(mineig)&&mineig>1e-8&&is.finite(as.numeric(logLik(fit)))
    if(accepted) {
      extracted<-tryCatch({
        k<-nrow(edges);stopifnot(as.integer(fit$qmodel$npars)==k)
        col<-match(expo,fit$qcmodel$covlabels);stopifnot(!is.na(col))
        idx<-k+(col-1L)*k+seq_len(k)
        theta<-as.numeric(fit$paramdata$params[idx]);V<-fit$covmat[idx,idx,drop=FALSE]
        stopifnot(max(abs(V-t(V)))<1e-8,all(diag(V)>=0))
        V<-(V+t(V))/2
        ev<-eigen(V,symmetric=TRUE,only.values=TRUE)$values
        stopifnot(min(ev)>=-1e-8*max(1,max(abs(ev))))
        cm<-coef(fit)[[expo]]
        stopifnot(isTRUE(all.equal(theta,vapply(seq_len(k),function(j)cm[edges$from_index[j],edges$to_index[j]],0.0),tolerance=1e-10)),
          all(is.finite(theta)),all(is.finite(V)))
        qhat<-msm::qmatrix.msm(fit,covariates="mean",ci="none")
        warm<-setNames(lapply(seq_along(fit$qcmodel$covlabels),function(j)
          as.numeric(fit$paramdata$params[k+(j-1)*k+seq_len(k)])),fit$qcmodel$covlabels)
        list(theta=theta,covariance=V)
      },error=function(e){error<<-paste("extraction:",conditionMessage(e));NULL})
      if(is.null(extracted))accepted<-FALSE else block<-extracted
    }
  }
  status<-if(accepted)"ACCEPTED" else "REJECTED_REVIEW_REQUIRED"
  obj<-list(signature=signature,status=status,model=label,imputation=imp,fit=fit,
    block=block,q_at_means=qhat,warm_covariates=warm,constraints=constraint,
    warnings=unique(warnings),error=error,elapsed_seconds=as.numeric(difftime(Sys.time(),begin,units="secs")),
    n_people=length(unique(z$id)),n_records=nrow(z),hessian_min_eigenvalue=mineig)
  save_atomic(obj,path)
  fwrite(data.table(model=label,status=status,pid=Sys.getpid(),n_people=obj$n_people,n_records=obj$n_records,
    elapsed_seconds=obj$elapsed_seconds,convergence=if(is.null(fit))NA_integer_ else fit$opt$convergence,
    found_se=if(is.null(fit))FALSE else fit$foundse,covariance_finite=covariance_ok,hessian_min_eigenvalue=mineig,
    warnings=paste(unique(warnings),collapse="; "),error=if(is.null(error))"" else error,signature=signature),
    file.path(out,"model_qc",paste0(label,"_status.csv")))
  message(label,": ",status,"; seconds=",round(obj$elapsed_seconds,1))
  if(!accepted)stop("Fit not accepted: ",label)
  invisible(obj)
}
pool_models <- function() {
  paths<-file.path(model_dir,sprintf("MI_%02d.rds",1:20))
  stopifnot(all(file.exists(paths)))
  fits<-lapply(paths,readRDS)
  stopifnot(all(vapply(fits,function(f)identical(f$status,"ACCEPTED")&&identical(f$signature,signature),TRUE)))
  stopifnot(identical(vapply(fits,function(f)as.integer(f$imputation),1L),1:20),
    length(unique(vapply(fits,function(f)f$n_people,1L)))==1L,
    all(vapply(fits,function(f)length(f$block$theta)==nrow(edges)&&
      identical(dim(f$block$covariance),c(nrow(edges),nrow(edges)))&&
      all(is.finite(f$block$covariance))&&all(is.finite(f$block$theta)),TRUE)))
  Q<-sapply(fits,function(f)f$block$theta)
  m<-ncol(Q);qbar<-rowMeans(Q);U<-Reduce(`+`,lapply(fits,function(f)f$block$covariance))/m
  B<-cov(t(Q));T<-U+(1+1/m)*B
  scalar<-function(L) {
    est<-drop(L%*%qbar);u<-diag(L%*%U%*%t(L));b<-diag(L%*%B%*%t(L));se<-sqrt(diag(L%*%T%*%t(L)))
    r<-(1+1/m)*b/u;df<-ifelse(b>1e-15,(m-1)*(1+1/r)^2,Inf)
    data.table(log_effect=est,se=se,df=df,HR=exp(est),lower=exp(est-qt(.975,df)*se),
      upper=exp(est+qt(.975,df)*se),p=2*pt(-abs(est/se),df))
  }
  all<-cbind(edges,scalar(diag(nrow(edges))))
  target_table<-all[target];target_table[,p_BH_4:=p.adjust(p,"BH")]
  C<-matrix(0,2,nrow(edges));C[1,target[1:2]]<-c(1,-1);C[2,target[3:4]]<-c(1,-1)
  contrast<-cbind(data.table(origin=c("No disease","HTN only"),contrast="Heart minus stroke"),scalar(C))
  contrast[,p_BH_2:=p.adjust(p,"BH")]
  joint<-function(L,label){v<-L%*%T%*%t(L);v<-(v+t(v))/2
    eigenvalues<-eigen(v,symmetric=TRUE,only.values=TRUE)$values
    stopifnot(all(is.finite(eigenvalues)),min(eigenvalues)>1e-12*max(eigenvalues),qr(v)$rank==nrow(v))
    a<-L%*%qbar;stat<-drop(t(a)%*%solve(v,a));
    data.table(test=label,chi_square=stat,df=nrow(L),p=pchisq(stat,nrow(L),lower.tail=FALSE),
      reference="Asymptotic chi-square Wald on Rubin total covariance")}
  tests<-rbind(joint(diag(nrow(edges))[target,,drop=FALSE],"Four target coefficients = 0"),
               joint(C,"Two heart-stroke coefficient differences = 0"))
  target_table[,`:=`(global_four_edge_p=tests$p[1],passes_global_and_BH=(tests$p[1]<.05 & p_BH_4<.05))]
  fwrite(all,file.path(out,"results/all_edge_MI20.csv"));fwrite(target_table,file.path(out,"results/four_target_MI20.csv"))
  fwrite(contrast,file.path(out,"results/heart_stroke_contrasts_MI20.csv"));fwrite(tests,file.path(out,"results/joint_tests_MI20.csv"))
  save_atomic(list(signature=signature,m=m,Q=Q,U=U,B=B,T=T,qbar=qbar,C=C,targets=target),file.path(out,"results/rubin_matrices_MI20.rds"))
  fwrite(data.table(status="MI20_POOLED_AWAITING_INDEPENDENT_REVIEW_AND_FIGURES",m=20L,signature=signature),
    file.path(out,"model_qc/delivery_status.csv"))
  print(target_table);print(contrast);print(tests)
}
if(mode=="pilot")fit_one(0L) else if(mode=="component"){
  stopifnot(component%in%1:20);fit_one(component)
}else if(mode=="pool")pool_models() else stop("Use pilot, component <1..20>, or pool")
