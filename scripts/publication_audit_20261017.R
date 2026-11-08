#!/usr/bin/env Rscript
# Read-only audit of frozen inputs. No imputation, model fitting or original-file mutation.
options(stringsAsFactors=FALSE, digits=17)
Sys.setlocale('LC_CTYPE','English_United States.utf8')
if (dir.exists('.runtime/charls_diagnostics_R')) .libPaths(c(normalizePath('.runtime/charls_diagnostics_R'), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(jsonlite); library(digest)})
setDTthreads(1L)
out <- 'outputs/publication_methods_audit_20261017'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
F <- 'outputs/hrs_elsa_multistate_final_20260909'
C1 <- 'outputs/charls_external_extension_20260913'
C2 <- 'outputs/charls_cesd9_duration_20260914'
C3 <- 'outputs/charls_dynamic_sleep_20260920'
inputs <- character()
rd <- function(p) {inputs <<- unique(c(inputs,p)); readRDS(p)}
wr <- function(z,n) fwrite(z,file.path(out,paste0(n,'.csv')),bom=TRUE,na='NA')
qs <- function(x) as.list(setNames(as.numeric(quantile(x,c(0,.25,.5,.75,1),type=7,na.rm=TRUE)),c('min','q1','median','q3','max')))
fu <- list(); per <- list(); waves <- list(); carries <- list(); carrydist<-list()
for (a in c('hrs_5state_full','elsa_4state_w4_w10','elsa_5state_through_2012')) {
  c <- if(grepl('^hrs',a)) 'hrs' else 'elsa'
  x <- as.data.table(rd(file.path(F,'data_multistate',paste0(c,'_msm_panel.rds'))))
  x[,state_msm:=as.integer(state)+1L]; x[,time_msm:=event_time-exposure_entry_time]
  if(a=='elsa_4state_w4_w10') x<-x[event_type=='visit' & wave %in% 4:10 & state_msm %in% 1:4]
  else if(a=='elsa_5state_through_2012') x<-x[event_time<=2012.999999 & state_msm %in% 1:5]
  else x<-x[state_msm %in% 1:5]
  x<-x[!is.na(id)&is.finite(time_msm)&time_msm>=0]; x<-x[id %in% x[,.N,by=id][N>=2,id]]
  setorder(x,id,time_msm,state_msm)
  z<-x[,.(followup=max(event_time)-min(event_time),n_records=.N,death=any(state_msm==5)),by=id]
  expect<-c(hrs_5state_full=29029,elsa_4state_w4_w10=11097,elsa_5state_through_2012=9244)[a]
  stopifnot(nrow(z)==expect,all(z$followup>0))
  fu[[a]]<-cbind(data.table(analysis=a,n=nrow(z),records=nrow(x),deaths=sum(z$death)),as.data.table(qs(z$followup)))
  per[[a]]<-z[,.(analysis=a,followup,n_records,death)] # deliberately no person identifiers exported
  w<-x[event_type=='visit',.(visits=.N,complete_observed=sum(!is.na(insomnia_score_observed)),carried=sum(is.na(insomnia_score_observed)&!is.na(insomnia_score_primary))),by=wave]
  w[,analysis:=a]; waves[[a]]<-w
  src<-as.data.table(rd(paste0('private/stage/output/data_clean/',c,'_clean.rds')))
  src<-unique(src[event_type=='visit'&!is.na(insomnia_score_observed),.(id,insomnia_source_wave=wave,source_time=event_time)])
  m<-merge(x[event_type=='visit'],src,by=c('id','insomnia_source_wave'),all.x=TRUE)
  stopifnot(nrow(m)==nrow(x[event_type=='visit']))
  stopifnot(!anyNA(m$source_time))
  m[,carry_years:=event_time-source_time]
  cc<-m[is.na(insomnia_score_observed)&!is.na(insomnia_score_primary)]
  carries[[a]]<-cbind(data.table(analysis=a,all_visits=nrow(m),carried_visits=nrow(cc),people_carried=uniqueN(cc$id),source_time_unresolved=sum(is.na(cc$source_time)),max_wave_gap=max(cc$wave-cc$insomnia_source_wave)),as.data.table(qs(cc$carry_years)))
  carrydist[[a]]<-cc[,.(analysis=a,wave,source_wave=insomnia_source_wave,carry_years)]
}
b1<-rd(file.path(C1,'data/charls_model_bundle.rds'))
b3<-rd(file.path(C3,'data/charls_dynamic_model_bundle.rds'))
for(a in c('charls_static','charls_dynamic')) {
  b<-if(a=='charls_static') b1 else b3
  x<-as.data.table(b$panel); z<-x[,.(followup=max(time)-min(time),n_records=.N,death=any(state==9)),by=id]
  stopifnot(nrow(z)==if(a=='charls_static')12036 else 9397,all(z$followup>0))
  fu[[a]]<-cbind(data.table(analysis=a,n=nrow(z),records=nrow(x),deaths=sum(z$death)),as.data.table(qs(z$followup)))
  per[[a]]<-z[,.(analysis=a,followup,n_records,death)]
}
wr(rbindlist(fu),'followup_summary');wr(rbindlist(per),'followup_distribution_no_identifiers')
wr(rbindlist(waves),'sleep_measurement_by_wave');wr(rbindlist(carries),'sleep_carry_summary');wr(rbindlist(carrydist),'sleep_carry_distribution_no_identifiers')
paths<-c(HRS=file.path(F,'data_multistate/mi_baseline/hrs_5state_full__M4_MI20__mids.rds'),
 ELSA_main=file.path(F,'data_multistate/mi_baseline/elsa_4state_w4_w10__M4_MI20__mids.rds'),
 ELSA_death=file.path(F,'data_multistate/mi_baseline/elsa_5state_through_2012__M4_MI20__mids.rds'),
 CHARLS_static=file.path(C1,'data/charls_baseline_mice20.rds'),
 CHARLS_duration=file.path(C2,'data/charls_stage2_mice20.rds'),
 CHARLS_dynamic=file.path(C3,'data/charls_dynamic_mice20.rds'))
miss<-list(); traces<-list(); diag<-list(); events<-list(); summaries<-list();preds<-list();support<-list()
for (a in names(paths)) {
  imp<-rd(paths[[a]])
  miss[[a]]<-rbindlist(lapply(names(imp$data),function(v)data.table(analysis=a,variable=v,n=nrow(imp$data),missing=sum(is.na(imp$data[[v]])),missing_pct=100*mean(is.na(imp$data[[v]])),method=imp$method[[v]],levels_or_range=if(is.factor(imp$data[[v]]))paste(levels(imp$data[[v]]),collapse=' | ')else paste(range(imp$data[[v]],na.rm=TRUE),collapse=' | '))))
  le<-as.data.table(imp$loggedEvents); if(nrow(le)) {le[,analysis:=a];events[[a]]<-le}
  summaries[[a]]<-data.table(analysis=a,m=imp$m,iterations=imp$iteration,seed=imp$seed,mice_version=as.character(imp$version),logged_events=nrow(le),chain_summary_caution='Short imputation chains; chain summary diagnostics are screening, not proof of convergence or MAR')
  pm<-as.data.table(as.table(imp$predictorMatrix));setnames(pm,c('target','predictor','included'));pm[,analysis:=a];preds[[a]]<-pm
  vars<-names(imp$method)[nzchar(imp$method)]
  for(v in vars) {
    native<-imp$data[[v]]; filled<-imp$imp[[v]]
    obs<-native[!is.na(native)]; valid<-if(is.factor(native))all(as.character(unlist(filled))%in%levels(native)) else all(is.finite(as.matrix(filled))) && all(as.matrix(filled)>=min(obs)&as.matrix(filled)<=max(obs))
    support[[paste(a,v)]]<-data.table(analysis=a,variable=v,method=imp$method[[v]],n_imputed_per_chain=nrow(filled),all_20_complete=!anyNA(filled),within_observed_support=valid)
    for(metric in c('mean','variance')) {
      mat<-if(metric=='mean')imp$chainMean[v,,] else imp$chainVar[v,,]
      dt<-CJ(iteration=seq_len(nrow(mat)),chain=seq_len(ncol(mat)));dt[,value:=mat[cbind(iteration,chain)]];dt[,`:=`(analysis=a,variable=v,metric=metric)];traces[[paste(a,v,metric)]]<-dt
      finite<-all(is.finite(mat)); wh<-if(finite)mean(apply(mat,2,var))else NA_real_
      rh<-if(finite && is.finite(wh)&&wh>0) sqrt(((nrow(mat)-1)/nrow(mat)*wh+var(colMeans(mat)))/wh) else NA_real_
      rr<-if(finite && requireNamespace('posterior',quietly=TRUE))tryCatch(posterior::rhat(mat),error=function(e)NA_real_)else NA_real_
      ac<-if(finite)mean(apply(mat,2,function(y) if(sd(y)>0)cor(y[-length(y)],y[-1])else NA_real_),na.rm=TRUE)else NA_real_
      diag[[paste(a,v,metric)]]<-data.table(analysis=a,variable=v,metric=metric,iterations=nrow(mat),chains=ncol(mat),all_finite=finite,classic_Rhat=rh,rank_split_Rhat=rr,mean_lag1_ac=ac,first_half_mean=mean(mat[seq_len(floor(nrow(mat)/2)),],na.rm=TRUE),second_half_mean=mean(mat[(floor(nrow(mat)/2)+1):nrow(mat),],na.rm=TRUE),reason=if(!finite)'Undefined trace statistic, e.g. variance with one missing observation'else'')
    }
  }
}
wr(rbindlist(miss),'mi_missingness_all_variables');wr(rbindlist(summaries),'mi_settings');wr(rbindlist(traces),'mi_chain_traces');wr(rbindlist(diag),'mi_chain_diagnostics');wr(rbindlist(events,fill=TRUE),'mi_logged_events_full');wr(rbindlist(preds),'mi_predictor_matrix_long');wr(rbindlist(support),'mi_support_checks')
le<-rbindlist(events,fill=TRUE);wr(le[,.N,by=.(analysis,dep,meth,out)],'mi_logged_event_counts')
# New descriptive selection audit uses original incomplete W1 covariates, not filled MI values.
base<-as.data.table(copy(b1$baseline));base[,included:=as.character(id)%in%as.character(b3$baseline$id)]
pp<-as.data.table(b1$panel);setorder(pp,id,time);blstate<-pp[,.SD[1],by=id][,.(id,baseline_state=state)]
base<-merge(base,blstate,by='id'); sel<-list()
for(v in setdiff(names(base),c('id','included'))) {
  x<-base[[v]]; catvar<-is.factor(x)||v%in%c('female','married','rural','smoker','drinker','diabetes','lung','pain','baseline_state')
  for(lv in if(catvar)sort(unique(as.character(x[!is.na(x)])))else 'continuous') {
    vals<-if(catvar)as.numeric(as.character(x)==lv)else as.numeric(x)
    aa<-vals[base$included];bb<-vals[!base$included];ma<-mean(aa,na.rm=TRUE);mb<-mean(bb,na.rm=TRUE)
    den<-if(catvar)sqrt((ma*(1-ma)+mb*(1-mb))/2)else sqrt((var(aa,na.rm=TRUE)+var(bb,na.rm=TRUE))/2)
    sel[[paste(v,lv)]]<-data.table(variable=v,level=lv,n_included=sum(base$included),n_excluded=sum(!base$included),included_missing=sum(is.na(aa)),excluded_missing=sum(is.na(bb)),included_mean_or_fraction=ma,excluded_mean_or_fraction=mb,SMD=if(den>0)(ma-mb)/den else NA_real_)
  }
}
wr(rbindlist(sel),'dynamic_selection_baseline_SMD')
# Exposure vs W1-W2 disease change: descriptive temporal overlap, not prospective effect.
lm<-as.data.table(b3$landmark)
join<-merge(lm,blstate,by='id');join[,disease_changed_by_W2:=baseline_state!=w2_state]
wr(join[,.(n=.N,disease_changed_by_W2=sum(disease_changed_by_W2),pct_disease_changed=100*mean(disease_changed_by_W2)),by=pattern],'dynamic_exposure_window_disease_change')
wr(as.data.table(b1$edges),'charls_allowed_edges');wr(data.table(code=b3$censor_codes,compatible_states=vapply(b3$censor_states,paste,collapse='|',FUN.VALUE='')),'dynamic_censor_state_sets')
write_json(list(created=format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'),scope='Read-only audit; no fitted scientific models changed',input_hashes=setNames(vapply(inputs,function(p)digest(file=p,algo='sha256'),''),inputs),R=R.version.string,posterior=if(requireNamespace('posterior',quietly=TRUE))as.character(packageVersion('posterior'))else NA_character_),file.path(out,'provenance.json'),pretty=TRUE,auto_unbox=TRUE)
print(rbindlist(fu)); print(rbindlist(carries));print(rbindlist(summaries));cat('AUDIT_COMPLETE\n')
