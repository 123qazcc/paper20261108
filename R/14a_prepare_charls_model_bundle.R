#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,contrasts=c("contr.treatment","contr.treatment"))
suppressPackageStartupMessages({library(data.table);library(mice);library(digest)})
setDTthreads(2)
out<-"outputs/charls_external_extension_20260913"
for(s in c("data","model_qc","tables","logs"))dir.create(file.path(out,s),recursive=TRUE,showWarnings=FALSE)
rawdir<-file.path(out,"raw_disease_audit");mortdir<-file.path(out,"mortality_audit")
rawgate<-readRDS(file.path(rawdir,"disease_routing_gate.rds"))
stopifnot(identical(rawgate$status,"AUDITED_RAW_HISTORY_READY_WITH_LIMITATIONS"),isTRUE(rawgate$all_checks_pass))
# Mortality audit must have its completed acceptance checklist before preparation.
mort_checks<-fread(file.path(mortdir,"mortality_checks.csv"))
stopifnot(all(mort_checks$pass))
files<-c(file.path(out,"data/charls_observed_exposures_tristate_panel.rds"),
 file.path(rawdir,c("charls_raw_disease_components_audited.rds","charls_raw_baseline_covariates_audited.rds","disease_routing_gate.rds")),
 file.path(mortdir,c("charls_death_candidates.rds","charls_raw_tracking_status.rds","mortality_checks.csv")),
 file.path(out,"adjustment_audit/all_W1_adjustment_supplement.rds"),
 "analysis/charls_model_stage2_contract.md","R/13b_audit_charls_disease_routing.R","R/13c_audit_charls_mortality.R",
 "R/13d_audit_charls_adjustment_variables.R","R/14a_prepare_charls_model_bundle.R")
stopifnot(all(file.exists(files)))
hashes<-setNames(vapply(files,function(f)digest(file=f,algo="sha256"),""),files)
d<-as.data.table(readRDS(files[1]));raw<-as.data.table(readRDS(files[2]))
rawcov<-as.data.table(readRDS(files[3]));death<-as.data.table(readRDS(file.path(mortdir,"charls_death_candidates.rds")))
track<-as.data.table(readRDS(file.path(mortdir,"charls_raw_tracking_status.rds")))
extra<-as.data.table(readRDS(file.path(out,"adjustment_audit/all_W1_adjustment_supplement.rds")))
d<-merge(d,raw,by=c("id","wave"),all.x=TRUE,sort=FALSE)
# Exclude death-wave health artifacts before accumulating any reported history.
d<-d[!track[died==1,.(id,wave)],on=.(id,wave)]
setorder(d,id,event_time)
everpos<-function(x){positive<-cumsum(!is.na(x)&x==1)>0;ifelse(positive,1L,ifelse(!is.na(x)&x==0,0L,NA_integer_))}
for(v in c("hypertension","heart","stroke"))d[,(paste0(v,"_hist")):=everpos(get(paste0(v,"_raw_audited"))),by=id]
d[, bits:=hypertension_hist+2L*heart_hist+4L*stroke_hist]
baseline<-copy(d[wave==1 & !is.na(age) & age>=50 & !is.na(restless_score_observed) & !is.na(bits)])
baseline<-merge(baseline,rawcov[,setdiff(names(rawcov),"wave"),with=FALSE],by="id",all.x=TRUE,sort=FALSE)
baseline[extra,on="id",education4_raw:=i.education4_W1]
baseline[,`:=`(restless_bl=as.numeric(restless_score_observed),age5=(age-65)/5,female=1-sex_male,
 education4=factor(education4_raw,levels=1:4),married=fifelse(marital_status %in% 1:3,1,
  fifelse(marital_status %in% c(4,5,7,8),0,NA_real_)),rural=as.numeric(hrural),
 wealth_asinh=asinh(hatotfa/10000),bmi5=(bmi-25)/5,smoker=as.numeric(current_smoker_raw_audited),
 drinker=as.numeric(current_drinker_raw_audited),diabetes=as.numeric(diabetes_raw_audited),
 lung=as.numeric(lung_raw_audited),pain=as.numeric(pain_raw_audited),entry_time=event_time)]
baseline[!is.finite(wealth_asinh),wealth_asinh:=NA_real_]
core<-c("age5","female","education4","married","rural","wealth_asinh","bmi5","smoker","drinker","diabetes","lung","pain")
saveRDS(baseline,file.path(out,"data/charls_raw_W1_baseline_stage2.rds"))

# Every genuine living visit contributes a compatible state set. Unknown components
# are not diagnosed as negative and do not require estimating misclassification.
panel<-merge(d,baseline[,.(id,entry_time)],by="id",sort=FALSE)
panel[,time:=event_time-entry_time]
panel[,state:=as.integer(bits)+1L]
panel[,pattern:=9L*fifelse(is.na(hypertension_hist),2L,as.integer(hypertension_hist))+
 3L*fifelse(is.na(heart_hist),2L,as.integer(heart_hist))+fifelse(is.na(stroke_hist),2L,as.integer(stroke_hist))]
panel[is.na(state),state:=100L+pattern]
states<-data.table(bits=0:7,hypertension=bitwAnd(0:7,1L),heart=as.integer(bitwAnd(0:7,2L)>0),
 stroke=as.integer(bitwAnd(0:7,4L)>0))
codes<-sort(unique(panel[state>=100,state]))
censor_states<-lapply(codes,function(code){r<-panel[state==code][1];
 states[(is.na(r$hypertension_hist)|hypertension==r$hypertension_hist)&
  (is.na(r$heart_hist)|heart==r$heart_hist)&(is.na(r$stroke_hist)|stroke==r$stroke_hist),bits+1L]})
panel<-panel[,.(id,wave,time,state,obstype=1L,record_type="living_visit")]
dd<-merge(death[usable==TRUE],baseline[,.(id,entry_time)],by="id")
death_rows<-dd[,.(id,wave=first_death_wave,time=model_observation_time-entry_time,state=9L,obstype=1L,record_type="interval_death")]
panel<-rbindlist(list(panel,death_rows),use.names=TRUE)
setorder(panel,id,time)
eligible<-panel[,.(n=.N,positive_followup=any(time>0)),by=id][n>=2 & positive_followup,id]
panel<-panel[id %in% eligible];baseline<-baseline[id %in% eligible]
setorder(baseline,id);setorder(panel,id,time)
checks<-list()
addcheck<-function(label,value){checks[[length(checks)+1]]<<-data.table(check=label,pass=isTRUE(value))}
addcheck("unique_subject_time",!anyDuplicated(panel[,.(id,time)]))
addcheck("no_negative_times",all(is.finite(panel$time)&panel$time>=0))
addcheck("every_entry_time_zero",all(panel[,min(time),by=id]$V1==0))
addcheck("baseline_IDs_equal_panel_IDs",setequal(baseline$id,panel$id))
addcheck("valid_baseline_exposure",all(baseline$restless_bl %in% 0:3))
addcheck("compatible_partial_state_sets",all(lengths(censor_states)>1)&all(unlist(censor_states)%in%1:8))
endcheck<-panel[,.(death_n=sum(state==9),last_state=last(state)),by=id]
addcheck("death_once_and_last",all(endcheck$death_n<=1 & (endcheck$death_n==0 | endcheck$last_state==9)))

# Fixed lattice: 12 single-disease increments and 8 death edges.
edges<-CJ(from=0:7,to=0:8)[(to==8)|(bitwAnd(from,to)==from & (to-from)%in%c(1L,2L,4L))]
edges[,`:=`(from_index=from+1L,to_index=to+1L,transition=paste(from,to,sep="-"))]
targets<-c("0-2","0-4","1-3","1-5")
edges[,exposure_group:=fcase(transition==targets[1],1L,transition==targets[2],2L,
 transition==targets[3],3L,transition==targets[4],4L,to==8,8L,to-from==1,5L,to-from==2,6L,to-from==4,7L)]
setorder(edges,from,to)
q<-matrix(0,9,9);q[cbind(edges$from_index,edges$to_index)]<-ifelse(edges$to==8,.025,.025)
addcheck("20_edges_and_four_free_targets",nrow(edges)==20 & all(targets %in% edges$transition) & all(edges$exposure_group[match(targets,edges$transition)]==1:4))

# MI uses baseline predictors and disease-specific observed-history summaries.
history<-merge(d[id %in% eligible],baseline[,.(id,entry_time)],by="id",sort=FALSE)
aux<-history[,{
 ti<-event_time-entry_time
 first<-function(x){ii<-which(!is.na(x)&x==1);if(length(ii))min(ti[ii])else max(ti)}
 .(htn0=hypertension_hist[which.min(ti)],heart0=heart_hist[which.min(ti)],stroke0=stroke_hist[which.min(ti)],
   ever_htn=as.numeric(any(hypertension_hist==1,na.rm=TRUE)),ever_heart=as.numeric(any(heart_hist==1,na.rm=TRUE)),
   ever_stroke=as.numeric(any(stroke_hist==1,na.rm=TRUE)),first_htn_observed=first(hypertension_hist),
   first_heart_observed=first(heart_hist),first_stroke_observed=first(stroke_hist),last_living_time=max(ti),n_visits=.N)
},by=id]
aux[,recorded_death:=as.numeric(id %in% dd$id)]
aux[dd,on="id",death_or_last_time:=i.model_observation_time-i.entry_time]
aux[is.na(death_or_last_time),death_or_last_time:=last_living_time]
mi_data<-merge(baseline[,c("id","restless_bl",core),with=FALSE],aux,by="id",sort=TRUE)
ids<-mi_data$id;mi_data[,id:=NULL]
addcheck("auxiliary_variables_complete",!anyNA(mi_data[,setdiff(names(mi_data),core),with=FALSE]))
miss<-rbindlist(lapply(names(mi_data),function(v)data.table(variable=v,n=nrow(mi_data),missing=sum(is.na(mi_data[[v]])),
  role=if(v %in% core)"baseline_covariate" else "fixed_exposure_or_observed_outcome_auxiliary")))
fwrite(miss,file.path(out,"model_qc/imputation_missingness.csv"))
method<-setNames(rep("",ncol(mi_data)),names(mi_data))
for(v in core)if(anyNA(mi_data[[v]]))method[v]<-if(v=="education4")"polyreg" else if(v %in%c("female","married","rural","smoker","drinker","diabetes","lung","pain"))"logreg" else "pmm"
for(v in names(method)[method=="logreg"])mi_data[,(v):=factor(get(v),levels=0:1)]
pred<-matrix(1L,ncol(mi_data),ncol(mi_data),dimnames=list(names(mi_data),names(mi_data)));diag(pred)<-0L
pred[method=="",]<-0L
prepchecks<-rbindlist(checks);fwrite(prepchecks,file.path(out,"model_qc/preparation_checks.csv"));stopifnot(all(prepchecks$pass))
fwrite(data.table(status="IMPUTATION_RUNNING",n_people=nrow(baseline),n_records=nrow(panel),n_deaths=nrow(death_rows[id %in% eligible]),pid=Sys.getpid()),
 file.path(out,"model_qc/preparation_status.csv"))
message("Preparing MICE20: people=",nrow(baseline),"; records=",nrow(panel),"; deaths=",sum(panel$state==9))
imp<-mice(as.data.frame(mi_data),m=20,maxit=20,method=method,predictorMatrix=pred,seed=20260913,printFlag=FALSE)
draws<-lapply(1:20,function(i){z<-as.data.table(complete(imp,i));z[,id:=ids]
 for(v in setdiff(core,"education4"))if(is.factor(z[[v]]))z[,(v):=as.numeric(as.character(get(v)))]
 z[,education4:=factor(education4,levels=1:4)];z[,c("id","restless_bl",core),with=FALSE]})
addcheck("20_complete_covariate_draws",length(draws)==20 & all(vapply(draws,function(z)!anyNA(z),TRUE)))
addcheck("MI_preserves_ID_and_exposure",all(vapply(draws,function(z)identical(z$id,ids)&identical(z$restless_bl,mi_data$restless_bl),TRUE)))
saveRDS(imp,file.path(out,"data/charls_baseline_mice20.rds"))
fwrite(as.data.table(imp$loggedEvents),file.path(out,"model_qc/mice_logged_events.csv"))
saveRDS(list(method=imp$method,predictor_matrix=imp$predictorMatrix,chainMean=imp$chainMean,chainVar=imp$chainVar),
 file.path(out,"model_qc/mice_diagnostics.rds"))
checks<-rbindlist(checks);stopifnot(all(checks$pass))
bun<-list(status="APPROVED_FOR_PILOT_AND_MI",panel=panel,baseline=baseline[,c("id","restless_bl",core),with=FALSE],
 imputations=draws,edges=edges,q_initial=q,core_covariates=core,censor_codes=codes,censor_states=censor_states,
 input_hashes=hashes,preparation_checks=checks,estimand="W1 restless-sleep frequency and prospective reported disease-history transitions; interval-observed death")
saveRDS(bun,file.path(out,"data/charls_model_bundle.rds"))
fwrite(edges,file.path(out,"tables/model_allowed_edges.csv"))
fwrite(checks,file.path(out,"model_qc/preparation_checks.csv"))
fwrite(data.table(status="READY_FOR_PILOT",n_people=nrow(baseline),n_records=nrow(panel),n_deaths=sum(panel$state==9),
 n_partial=sum(panel$state>=100),m=20L),file.path(out,"model_qc/preparation_status.csv"))
message("Model bundle ready; no association estimates produced by preparation.")
