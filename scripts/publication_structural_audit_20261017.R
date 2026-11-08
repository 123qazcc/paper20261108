#!/usr/bin/env Rscript
Sys.setlocale('LC_CTYPE','English_United States.utf8')
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(digest)})
out<-'outputs/publication_methods_audit_20261017';F<-'outputs/hrs_elsa_multistate_final_20260909'
C1<-'outputs/charls_external_extension_20260913';C3<-'outputs/charls_dynamic_sleep_20260920'
wr<-function(z,n)fwrite(z,file.path(out,paste0(n,'.csv')),bom=TRUE,na='NA')
# Parse only pure functions; never source the cleaning pipeline.
expr<-parse('R/01_clean_three_cohorts.R');en<-new.env(parent=globalenv())
for(e in expr)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&as.character(e[[2]])[1]%in%c('derive_cvd_state','cum_diagnosis'))eval(e,en)
states<-CJ(h=c(0L,1L,NA_integer_),heart=c(0L,1L,NA_integer_),stroke=c(0L,1L,NA_integer_))
states[,state:=en$derive_cvd_state(h,heart,stroke)];wr(states,'HRS_ELSA_all_27_state_classifications')
cc<-list()
for(c in c('hrs','elsa')){x<-as.data.table(readRDS(file.path(F,'data_multistate',paste0(c,'_msm_panel.rds'))));rec<-en$derive_cvd_state(x$hypertension,x$heart_disease,x$stroke,x$dead);cc[[c]]<-data.table(cohort=c,records=nrow(x),mismatches=sum((is.na(rec)!=is.na(x$state))|(!is.na(rec)&!is.na(x$state)&rec!=x$state)))}
wr(rbindlist(cc),'state_reconstruction_checks')
# Structural cause of logged events in the actual MIDS design.
ls<-list()
for(a in c('hrs_5state_full','elsa_5state_through_2012')){m<-readRDS(file.path(F,'data_multistate/mi_baseline',paste0(a,'__M4_MI20__mids.rds')));d<-m$data;ls[[a]]<-data.table(analysis=a,n=nrow(d),baseline_death_count=sum(as.character(d$state_msm)=='5'),death_equals_laststate5=all(d$death_observed==as.integer(as.character(d$last_state_msm)=='5')),logged_iterations=paste(range(m$loggedEvents$it),collapse='-'),logged_chains=uniqueN(m$loggedEvents$im),target_variables=uniqueN(m$loggedEvents$dep))}
wr(rbindlist(ls),'logged_events_structural_cause_checks')
b1<-readRDS(file.path(C1,'data/charls_model_bundle.rds'));b3<-readRDS(file.path(C3,'data/charls_dynamic_model_bundle.rds'))
base<-as.data.table(copy(b1$baseline));panel<-as.data.table(b1$panel);ids<-panel[wave==2 & record_type=='living_visit',unique(id)]
base<-base[id%in%ids];base[,included:=id%in%b3$baseline$id]
sel<-list()
for(v in setdiff(names(base),c('id','included'))){x<-base[[v]];catvar<-is.factor(x)||v%in%c('female','married','rural','smoker','drinker','diabetes','lung','pain');for(lv in if(catvar)sort(unique(as.character(x[!is.na(x)])))else'continuous'){
val<-if(catvar)as.numeric(as.character(x)==lv)else as.numeric(x);aa<-val[base$included];bb<-val[!base$included];ma<-mean(aa,na.rm=TRUE);mb<-mean(bb,na.rm=TRUE);den<-if(catvar)sqrt((ma*(1-ma)+mb*(1-mb))/2)else sqrt((var(aa,na.rm=TRUE)+var(bb,na.rm=TRUE))/2);sel[[paste(v,lv)]]<-data.table(variable=v,level=lv,n_included=sum(base$included),n_excluded=sum(!base$included),included_mean_or_fraction=ma,excluded_mean_or_fraction=mb,included_missing=sum(is.na(aa)),excluded_missing=sum(is.na(bb)),SMD=if(den>0)(ma-mb)/den else NA_real_)}}
wr(rbindlist(sel),'dynamic_selection_among_W2_living_SMD')
p<-as.data.table(b3$panel);w3<-p[wave==3 & record_type=='living_visit' & state%in%1:8,.(id,entry=time,entry_state=state)]
post<-merge(p,w3,by='id');post<-post[time>=entry];eligible<-post[,.(n=.N,positive=any(time>entry)),by=id][n>=2&positive,id]
wr(w3[id%in%eligible,.(n=.N),by=entry_state],'proposed_W3_landmark_feasibility_only')
wr(data.table(candidate='W3 re-landmark using frozen W1 W2 pattern',known_W3_living=nrow(w3),with_subsequent_record=length(eligible),status='DESIGN FEASIBILITY ONLY - NO EFFECT MODEL FITTED'),'proposed_sensitivity_feasibility')
codes<-c('R/01_clean_three_cohorts.R','R/03_build_harmonized_exposures_baseline.R','R/05_prepare_multistate_data.R','R/07_fit_markov_models.R','R/07b_protocol_sensitivities.R')
wr(rbindlist(lapply(codes,function(p){q<-file.path(F,'code',p);data.table(script=p,sha256=digest(file=p,algo='sha256'),frozen_copy_exists=file.exists(q),matches_frozen=if(file.exists(q))identical(digest(file=p,algo='sha256'),digest(file=q,algo='sha256'))else NA)})),'frozen_code_hash_checks')
wr(data.table(package=c('R','data.table','msm','mice','posterior'),audit_runtime_version=c(as.character(getRversion()),vapply(c('data.table','msm','mice'),function(p)as.character(packageVersion(p)),''),as.character(packageVersion('posterior',lib.loc='.runtime/charls_diagnostics_R')))),'audit_runtime_versions')
print(rbindlist(cc));print(rbindlist(ls));print(rbindlist(sel)[abs(SMD)>.1]);print(length(eligible))
