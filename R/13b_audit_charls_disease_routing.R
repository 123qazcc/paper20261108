#!/usr/bin/env Rscript
# Prospective raw-answer disease audit; no exposure-outcome association models.
# Does not edit phase 13a, its baseline, or any HRS/ELSA artifact.
options(stringsAsFactors=FALSE, width=140)
suppressPackageStartupMessages({library(data.table);library(haven);library(digest)})
root <- normalizePath(getwd(),winslash='/',mustWork=TRUE)
out <- file.path(root,'outputs/charls_external_extension_20260913/raw_disease_audit')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source_root <- Sys.getenv('CHARLS_SOURCE_ROOT','private/source')
raw_root <- file.path(source_root,'data_raw/1_CHARLS_China')
stage <- 'private/stage/charls_field_audit'
sha <- function(p) digest(file=p,algo='sha256')
put <- function(x,n) fwrite(x,file.path(out,paste0(n,'.csv')),bom=TRUE,na='')
checks<-list()
check<-function(n,ok,detail='')checks[[length(checks)+1L]]<<-data.table(check=n,pass=isTRUE(ok),detail=as.character(detail))
protected_paths<-c('R/13a_prepare_charls_external_extension.R','data_clean/charls_clean.rds','data_clean/hrs_clean.rds','data_clean/elsa_clean.rds',
 'R/11c_fit_disease_specific_models.R','outputs/charls_external_extension_20260913/data/charls_w1_descriptive_baseline.rds')
protected<-data.table(path=protected_paths,before=vapply(protected_paths,sha,''))
files<-data.table(wave=1:5,year=c(2011,2013,2015,2018,2020),
 original=file.path(raw_root,c('2011charls/health_status_and_functioning.dta','2013charls/Health_Status_and_Functioning.dta',
 '2015charls/Health_Status_and_Functioning.dta','2018charls/Health_Status_and_Functioning.dta','2020charls/Health_Status_and_Functioning.dta')),
 staged=file.path(stage,paste0('w',1:5,'_health.dta')))
files[,original_sha256:=vapply(original,sha,'')]
files[,staged_sha256:=vapply(staged,sha,'')]
for(w in 1:5)check(paste0('original_stage_hash_W',w),files$original_sha256[w]==files$staged_sha256[w])
stopifnot(all(vapply(checks,function(x)x$pass,TRUE)))
put(files,'raw_source_manifest')
e<-new.env(parent=baseenv());data('CHARLS',package='GlobalAgingR',envir=e)
pkg<-as.data.table(GlobalAgingR::CH_filter(e$charls,wave %in% 1:5))
pkg[,id:=sprintf('%012.0f',as.numeric(ID))]
check('package_version_2_0_0',as.character(packageVersion('GlobalAgingR'))=='2.0.0')
check('package_unique_keys',!anyDuplicated(pkg[,.(id,wave)]))
tri12<-function(x)fifelse(x==1,1L,fifelse(x==2,0L,NA_integer_))
flag<-function(x)replace(as.logical(x),is.na(x),FALSE)

# A raw negative is never propagated into a later missing assessment. For W2/W3,
# an interval-only negative needs a confirmed negative prior record. W3 entrants
# with only a since-last-interview question retain UNKNOWN, as required by the
# Harmonized CHARLS D manual (PDF pp128-129, printed pp127-128).
map_one_wave<-function(a,w,self_possible=FALSE){
 a<-copy(a);n<-nrow(a)
 a[,`:=`(audited=NA_integer_,provenance='unknown_no_defensible_current_route',
          ambiguity=FALSE,denial=FALSE,contradiction=FALSE)]
 assign_where<-function(cond,val,label){
  z<-which(flag(cond));if(length(z))a[z,`:=`(audited=as.integer(val),provenance=label)]
 }
 if(w==1L){assign_where(a$direct==1,1L,'direct_doctor_yes');assign_where(a$direct==2,0L,'direct_doctor_no')}
 if(w %in% 2:3){
  # Direct lifetime question is asked to new interviewees.
  assign_where(a$direct==1,1L,'direct_doctor_yes')
  assign_where(a$direct==2,0L,'direct_doctor_no')
  # Record agreement is not itself "yes": interpretation needs prior status.
  assign_where(a$preload==1 & a$confirm==1 & is.na(a$incident),1L,'doctor_preload_yes_confirmed')
  # Negative preloads at W3 are often blank. Only use actual previous-wave raw
  # history together with current confirmation AND the current interval answer.
  known_previous_negative<-flag(a$preload==2) |
   flag(is.na(a$preload) & a$prior_wave==w-1L & a$prior_ever==0L)
  previous_any_positive<-flag(a$preload==1) | flag(a$selfpreload==1)
  assign_where(known_previous_negative & !previous_any_positive & a$confirm==1 & a$incident==2,
               0L,'confirmed_previous_no_and_interval_doctor_no')
  assign_where(previous_any_positive & a$confirm==2 & a$incident==2,
               0L,'previous_positive_denied_and_interval_doctor_no')
  # A physician diagnosis since the preceding interview establishes ever yes
  # irrespective of missing earlier history, and is not backdated.
  assign_where(a$incident==1,1L,'interval_doctor_yes')
  if(!self_possible){
   assign_where(known_previous_negative & a$confirm==2 & is.na(a$incident),
                1L,'previous_doctor_no_disputed_to_yes_current_report')
  } else {
   # For HTN, the confirmation loop combines doctor and self-known history.
   # Do not promote a disputed combined negative to a doctor diagnosis without
   # a direct/current doctor-positive answer.
   z<-which(flag(known_previous_negative & !previous_any_positive & a$confirm==2 & is.na(a$incident) & is.na(a$direct)))
   if(length(z))a[z,`:=`(ambiguity=TRUE,provenance='unknown_disputed_combined_doctor_self_history')]
  }
  a[,denial:=flag(previous_any_positive & confirm==2)]
  # Impossible concurrent direct and repeat-route responses or opposing explicit
  # answers remain unknown rather than letting a priority order conceal conflict.
  a[,contradiction:=flag(direct %in% 1:2 & incident %in% 1:2 & direct!=incident) |
    flag(preload==1 & confirm==1 & incident==2)]
  a[contradiction==TRUE,`:=`(audited=NA_integer_,ambiguity=TRUE,provenance='unknown_internally_conflicting_raw_routes')]
  a[is.na(audited) & incident==2 & !ambiguity,
    `:=`(ambiguity=TRUE,provenance='unknown_interval_no_without_identifiable_previous_negative')]
 }
 if(w==4L){
  assign_where(a$preload==1 & a$confirm %in% 1:3,1L,'doctor_preload_confirmed_better_worse_same')
  assign_where(a$preload==1 & a$onset %in% 1:2 & is.na(a$confirm),1L,'doctor_preload_and_diagnosis_time_answer')
  assign_where(a$direct==1,1L,'direct_doctor_yes')
  assign_where(a$direct==2,0L,'direct_doctor_no')
  # W4 99 literally means NEVER had disease, unlike W5 99 (not at prior visit).
  a[,denial:=flag(confirm==99 | onset==99 | treatment_never==99 | control==99)]
  assign_where(a$denial & !(a$direct==1) & !(a$confirm %in% 1:3),0L,'explicit_never_disease_denial')
  # R's NA boolean handling must not suppress a denial when direct was skipped.
  assign_where(a$denial & !flag(a$direct==1) & !(a$confirm %in% 1:3),0L,'explicit_never_disease_denial')
  a[,contradiction:=denial & (flag(direct==1) | confirm %in% 1:3)]
  a[contradiction==TRUE,`:=`(audited=NA_integer_,ambiguity=TRUE,provenance='unknown_positive_and_never_disease_conflict')]
  if(self_possible)a[is.na(audited) & combinedpreload==1 & is.na(preload) & confirm %in% 1:3,
       `:=`(ambiguity=TRUE,provenance='unknown_self_known_not_doctor_confirmed')]
 }
 if(w==5L){
  assign_where(a$preload==1 & a$confirm %in% 1:3,1L,'doctor_preload_confirmed_better_worse_same')
  assign_where(a$direct==1,1L,'direct_doctor_yes')
  assign_where(a$direct==2,0L,'direct_doctor_no')
  a[,denial:=flag(confirm==99)]
  # W5 "not at previous interview" is followed by DA003. Never use it alone
  # as a current absence or erase an earlier positive.
  a[confirm==99 & is.na(direct),`:=`(audited=NA_integer_,ambiguity=TRUE,
      provenance='unknown_previous_denial_missing_current_doctor_answer')]
  a[,contradiction:=flag(confirm %in% 1:3 & direct==2)]
  a[contradiction==TRUE,`:=`(audited=NA_integer_,ambiguity=TRUE,provenance='unknown_internally_conflicting_raw_routes')]
 }
 a
}

# Synthetic safety tests, not empirical effect-size checks.
blank<-function(n=1L)data.table(id=as.character(seq_len(n)),direct=NA_real_,preload=NA_real_,selfpreload=NA_real_,
 confirm=NA_real_,incident=NA_real_,onset=NA_real_,combinedpreload=NA_real_,treatment_never=NA_real_,control=NA_real_,
 prior_wave=NA_integer_,prior_ever=NA_integer_)
a<-blank();a[,`:=`(confirm=1,incident=2)]
check('W3_interval_only_new_entrant_stays_unknown',is.na(map_one_wave(a,3)$audited))
a[,`:=`(prior_wave=2L,prior_ever=0L)]
check('W3_confirmed_no_plus_previous_observed_no_maps_no',map_one_wave(a,3)$audited==0L)
a[,`:=`(prior_wave=1L)]
check('W3_nonadjacent_previous_no_not_filled',is.na(map_one_wave(a,3)$audited))
a<-blank();a[,`:=`(preload=1,confirm=2,incident=2)]
check('W3_explicit_denial_retained_as_raw_no',map_one_wave(a,3)$audited==0L)
a<-blank();a[,`:=`(confirm=99)]
check('W5_prior_denial_alone_not_current_no',is.na(map_one_wave(a,5)$audited))
check('W4_never_disease_maps_raw_no',map_one_wave(a,4)$audited==0L)
a[,`:=`(direct=1)]
check('W4_positive_never_conflict_unknown',is.na(map_one_wave(a,4)$audited))
check('W5_prior_denial_then_new_doctor_yes_is_yes',map_one_wave(a,5)$audited==1L)

longs<-list();dictionaries<-list();prior<-list();baseline_cov<-NULL
for(w in 1:5){
 raw<-as.data.table(haven::read_dta(files$staged[w]))
 raw[,id:=as.character(ID)]
 if(w==1L)raw[,id:=paste0(substr(id,1,9),'0',substr(id,10,11))]
 check(paste0('raw_unique_ID_W',w),!anyDuplicated(raw$id))
 check(paste0('raw_12digit_ID_W',w),all(grepl('^[0-9]{12}$',raw$id)))
 check(paste0('all_raw_health_IDs_in_package_W',w),all(raw$id %in% pkg[wave==w,id]))
 if(w==1L){
  baseline_cov<-raw[,.(id,wave=1L,diabetes_raw_audited=tri12(as.numeric(da007_3_)),
   lung_raw_audited=tri12(as.numeric(da007_5_)),
   current_smoker_raw_audited=fifelse(as.numeric(da059)==2,0L,
    fifelse(as.numeric(da059)==1 & as.numeric(da061)==1,1L,
     fifelse(as.numeric(da059)==1 & as.numeric(da061)==2,0L,NA_integer_))),
   current_drinker_raw_audited=fifelse(as.numeric(da067) %in% 1:2,1L,
    fifelse(as.numeric(da067)==3,0L,NA_integer_)),pain_raw_audited=tri12(as.numeric(da041)))]
 }
 for(k in 1:3){
  item<-c(1L,7L,8L)[k];dis<-c('hypertension','heart','stroke')[k];pv<-c('hibpe','hearte','stroke')[k]
  fields<-if(w==1L)c(direct=sprintf('da007_%d_',item)) else if(w %in% 2:3)
   c(direct=sprintf('da007_%d_',item),preload=sprintf('zda007_%d_',item),selfpreload=sprintf('zda008_%d_',item),
     confirm=sprintf('da007_w2_1_%d_',item),incident=sprintf('da007_w2_2_%d_',item)) else if(w==4L)
   c(direct=sprintf('da007_%d_',item),preload=sprintf('zdiagnosed_%d_',item),combinedpreload=sprintf('zdisease_%d_',item),
     confirm=sprintf('da010_w2_2_%d_',item),onset=sprintf('da009_%d_',item),
     treatment_never=sprintf('da010_w4_%d__s99',item),control=sprintf('da010_w2_1_%d_',item)) else
   c(direct=sprintf('da003_%d_',item),preload=sprintf('zdisease_%d_',item),confirm=sprintf('da002_%d_',item))
  a<-blank(nrow(raw));a[,id:=raw$id]
  for(f in names(fields)){
   v<-fields[[f]]
   if(v %in% names(raw)){
    a[,(f):=as.numeric(raw[[v]])]
    dictionaries[[length(dictionaries)+1L]]<-data.table(wave=w,disease=dis,role=f,variable=v,
      label=paste(attr(raw[[v]],'label'),collapse=';'),
      values=paste(names(attr(raw[[v]],'labels')),attr(raw[[v]],'labels'),collapse=';'))
   }
  }
  if(w>1L){
   ix<-match(a$id,prior[[dis]]$id)
   a[,`:=`(prior_wave=prior[[dis]]$prior_wave[ix],prior_ever=prior[[dis]]$prior_ever[ix])]
  }
  a<-map_one_wave(a,w,self_possible=k==1L)
  a[,`:=`(wave=w,disease=dis,raw_no_after_previous_positive=flag(audited==0L & prior_ever==1L))]
  p<-pkg[wave==w,.(id,package_value=as.numeric(get(pv)))]
  p[,package_tri:=fifelse(package_value==1,1L,fifelse(package_value %in% c(0,2),0L,NA_integer_))]
  a<-merge(a,p,by='id',all=TRUE,sort=FALSE)
  a[,`:=`(wave=w,disease=dis)]
  a[is.na(provenance),`:=`(provenance='unknown_no_raw_health_module_record',ambiguity=FALSE,
     denial=FALSE,contradiction=FALSE,raw_no_after_previous_positive=FALSE)]
  if(w>1L){
   z<-which(a$provenance=='unknown_no_raw_health_module_record')
   if(length(z)){
    ix<-match(a$id[z],prior[[dis]]$id)
    a[z,`:=`(prior_wave=prior[[dis]]$prior_wave[ix],prior_ever=prior[[dis]]$prior_ever[ix])]
   }
  }
  check(paste0('mapped_tri_values_',dis,'_W',w),all(is.na(a$audited)|a$audited %in% 0:1))
  # Prospective history used only for interpreting the next wave's confirmation.
  nextprior<-a[,.(id,prior_wave=w,prior_ever=fifelse(prior_ever==1L | audited==1L,1L,audited))]
  nextprior[is.na(prior_ever) & a$audited==1L,prior_ever:=1L]
  # fifelse(NA | FALSE) must not accidentally discard an explicit raw negative.
  nextprior[,prior_ever:=fifelse(flag(a$prior_ever==1L)|flag(a$audited==1L),1L,a$audited)]
  if(w>1L){old<-prior[[dis]][!id %in% nextprior$id];nextprior<-rbindlist(list(nextprior,old))}
  prior[[dis]]<-nextprior
  longs[[length(longs)+1L]]<-a
 }
}
long<-rbindlist(longs,fill=TRUE);setorder(long,id,wave,disease)
check('long_unique_ID_wave_disease',!anyDuplicated(long[,.(id,wave,disease)]))
stopifnot(!anyDuplicated(long[,.(id,wave,disease)]))
wide<-dcast(long,id+wave~disease,value.var=c('audited','provenance','ambiguity','denial','contradiction','raw_no_after_previous_positive'))
for(dis in c('hypertension','heart','stroke')){
 setnames(wide,paste0('audited_',dis),paste0(dis,'_raw_audited'))
 for(p in c('provenance','ambiguity','denial','contradiction','raw_no_after_previous_positive'))setnames(wide,paste0(p,'_',dis),paste0(dis,'_',p))
}
check('audited_panel_key_set_equals_package',setequal(paste(wide$id,wide$wave),paste(pkg$id,pkg$wave)))
check('audited_panel_unique_keys',!anyDuplicated(wide[,.(id,wave)]))
check('baseline_original_direct_answers_only',all(long[wave==1,!is.na(audited)]==long[wave==1,direct %in% 1:2]))
basecompare<-merge(baseline_cov,pkg[wave==1,.(id,diabetes=as.numeric(diabe),lung=as.numeric(lunge),
  smoker=as.numeric(smoken),drinker=as.numeric(drinkl))],by='id')
covcompare<-rbindlist(lapply(c(diabetes='diabetes_raw_audited',lung='lung_raw_audited',smoker='current_smoker_raw_audited',drinker='current_drinker_raw_audited'),
 function(v){j<-which(c(diabetes='diabetes_raw_audited',lung='lung_raw_audited',smoker='current_smoker_raw_audited',drinker='current_drinker_raw_audited')==v)
  p<-c('diabetes','lung','smoker','drinker')[j];data.table(variable=v,n=nrow(basecompare),raw_missing=sum(is.na(basecompare[[v]])),
   package_missing=sum(is.na(basecompare[[p]])),known_discordance=sum(basecompare[[v]]!=basecompare[[p]],na.rm=TRUE),
   missing_discordance=sum(xor(is.na(basecompare[[v]]),is.na(basecompare[[p]]))))}))
put(covcompare,'baseline_covariate_raw_package_comparison')
summ<-long[,.(records=.N,package_missing=sum(is.na(package_tri)),audited_missing=sum(is.na(audited)),
 known_discordance=sum(audited!=package_tri,na.rm=TRUE),package_no_raw_yes=sum(package_tri==0 & audited==1,na.rm=TRUE),
 package_yes_raw_no=sum(package_tri==1 & audited==0,na.rm=TRUE),filled_package_missing=sum(is.na(package_tri)&!is.na(audited)),
 package_known_now_unknown=sum(!is.na(package_tri)&is.na(audited)),denials=sum(denial),ambiguous=sum(ambiguity),
 contradictions=sum(contradiction),raw_no_after_previous_positive=sum(raw_no_after_previous_positive)),by=.(wave,disease)]
setorder(summ,wave,disease);put(summ,'raw_package_disease_comparison')
put(long[,.N,by=.(wave,disease,provenance,audited)],'mapped_provenance_counts')
put(long[,.N,by=.(wave,disease,direct,preload,selfpreload,confirm,incident,onset,combinedpreload,
 treatment_never,control,prior_wave,prior_ever,audited,provenance,package_tri)],'per_ID_comparison_aggregated_patterns')
put(rbindlist(dictionaries),'raw_variable_dictionary')
saveRDS(wide,file.path(out,'charls_raw_disease_components_audited.rds'))
saveRDS(long,file.path(out,'charls_raw_disease_routing_per_ID.rds'))
saveRDS(baseline_cov,file.path(out,'charls_raw_baseline_covariates_audited.rds'))

# Official local source documents and exact provenance, with no disease subtyping
# inferred from abbreviated labels. These page references are PDF page numbers.
docs<-data.table(document=c('W1 questionnaire','W2 questionnaire','W3 questionnaire','W4 questionnaire','W5 questionnaire','W5 routing flowchart','Harmonized CHARLS D manual'),
 path=c(file.path(raw_root,c('2011charls/CHARLS_wave1_问卷.pdf','2013charls/CHARLS_Wave2_Main_Questionnaire.pdf',
 '2015charls/CHARLS_2015_Questionnaire.pdf','2018charls/CHARLS_2018_Household_Questionnaire.pdf',
 '2020charls/CHARLS_2020_Questionnaire_Chinese.pdf','2020charls/CHARLS_2020_Flowchart_Chinese.pdf')),
 file.path(source_root,'harmonized/CHARLS/H_CHARLS_D_Data/Harmonized_CHARLS_D.pdf')),
 relevant_pdf_pages=c('60-61;67;71-72','75-76','60-61','71-75','32','28','128-129'),
 reviewed=c(rep('local full relevant text; W3 page61 and W5 flowchart page28 also visually reviewed',6),'local construction text and dispute labels'))
docs[,sha256:=vapply(path,sha,'')];put(docs,'questionnaire_source_manifest')
protected[,after:=vapply(path,sha,'')];put(protected,'protected_file_hashes')
check('all_protected_files_unchanged',all(protected$before==protected$after))
qa<-rbindlist(checks);put(qa,'disease_routing_QA')
status<-list(schema='charls-prospective-raw-routing-v1',status=if(all(qa$pass))'AUDITED_RAW_HISTORY_READY_WITH_LIMITATIONS' else 'FAILED',
 timestamp=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),all_checks_pass=all(qa$pass),checks=nrow(qa),
 raw_panel_rows=nrow(wide),baseline_changes_required=TRUE,
 interpretation='Prospective self-reported physician-diagnosis history; never reclassify prior waves using future disputes; apply ever-positive-only downstream.',
 limitations=c('W3 life-history entrants with interval-only negatives remain unknown.',
  'HTN record confirmation combines physician/self-known conditions; ambiguous routes remain unknown.',
  'Raw diagnosis denial is recorded, not proof earlier diagnosis was clinically false.',
  'This gate does not audit mortality coverage or permit causal/association claims.'),
 script_sha256=sha(file.path(root,'R/13b_audit_charls_disease_routing.R')))
saveRDS(status,file.path(out,'disease_routing_gate.rds'))
saveRDS(sessionInfo(),file.path(out,'session_info.rds'))
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
report<-c('# CHARLS 原始疾病路由审计（第二阶段）','',
 '## 结论与采用边界','',
 '三种疾病可构建独立、前瞻性原始回答历史，但不能仅给包内 W3 缺失补值。GlobalAgingR 的既往疾病字段与原始同期回答存在实质差异；本地 Harmonized CHARLS D 手册明确指出，会根据后续争议修改以前波次。第一阶段基线表必须保留为旧版本，并在主分析采用原始历史后重新生成。',
 '原始历史是曾报告医生诊断，不等同经过病历核实的临床真值。后续否认保留标志，不回改过去波次；下游应使用曾阳性仅正向保留规则。不得把先前阴性自动填进缺失访谈。','',
 '## 各波路由','',
 '- W1：DA007_1/7/8 医生诊断直接回答；1=是、2=否，其他未知。高血压不纳入单独自知患病 DA008。',
 '- W2/W3：先区分新访直接问题、预载诊断确认与自上次访问以来新诊断。DA007_W2_1=1 仅表示记录正确，不等于患病；需结合预载值。肯定的新诊断为阳性。否定的区间回答只有在本次确认前次阴性且有可识别的前次阴性，或明确否认前次阳性时，才为原始阴性。',
 '- W3：手册指出部分生活史/2015新入组者只问了自上次访问以来的患病，缺乏可靠曾经患病史。本脚本对没有可识别前次阴性的区间否定保留未知；不将 W3 的空预载直接视为无病。',
 '- W2/W3 高血压：确认环节合并医生诊断与自知患病；没有明确医生诊断支持的反向争议保持未知。',
 '- W4：医生诊断直接回答，或医生诊断预载加本波好转/恶化/相同确认。问卷 99 指从未患病，记录否认；与同波明确阳性矛盾时保持未知。',
 '- W5：DA003 医生诊断，或医生诊断预载加 DA002 的好转/恶化/相同确认。这里 DA002=99 只指上次访问未患病，之后必须查看 DA003；不能单凭 99 判当前无病。','',
 '## 疾病名称','',
 '心脏病为心肌梗死、冠心病、心绞痛、充血性心力衰竭或其他心脏问题的合并问题，不能拆成这些单独疾病。W4 原始变量标签的 Heart Attack 是缩写，完整问卷仍是广义心脏病。卒中为中风，W4 明确含脑梗和脑出血，但没有独立亚型回答。','',
 '## 原始与包数据逐人对照（各波疾病记录，不是独立人数相加）','',
 '|波次|疾病|包缺失|原始审计后缺失|双方已知但不同|恢复包缺失|包已知改为未知|原始否认|同波冲突|',
 '|---|---|---:|---:|---:|---:|---:|---:|---:|',
 apply(as.data.frame(summ[,.(wave,disease,package_missing,audited_missing,known_discordance,filled_package_missing,package_known_now_unknown,denials,contradictions)]),1,
  function(r)paste0('|',paste(r,collapse='|'),'|')),'',
 '## 输出与复核','',
 '- `charls_raw_disease_components_audited.rds`：id（12位字符）、wave、三疾病 `*_raw_audited`，以及逐病来源、歧义、否认、冲突、既往阳性后原始阴性标志。',
 '- `charls_raw_disease_routing_per_ID.rds`：内部逐ID原始问项、前次历史、包值及映射结果，供独立复核；汇报仅提供聚合计数。',
 '- `charls_raw_baseline_covariates_audited.rds`：W1 直接医生诊断糖尿病、肺病，当前吸烟、过去一年饮酒和经常疼痛；不利用未来争议。',
 '- `disease_routing_gate.rds`：仅为疾病映射就绪门，不代表死亡覆盖或关联模型已通过。',
 sprintf('- 自动质控：%d/%d 通过。三库清洗、第一阶段原件和冻结 HRS/ELSA 疾病拆分代码散列不变。',sum(qa$pass),nrow(qa)),
 '- 采用 nature-statistics 的测量/分析/推断分层：本审计无 HR、q 或 p 值，不作关联、因果或疾病特异强弱结论。','',
 '## 本地来源','',
 paste0('- ',docs$document,'；PDF 页 ',docs$relevant_pdf_pages,'；`',docs$path,'`。'))
writeLines(report,file.path(out,'CHARLS_疾病原始路由审计说明.md'),useBytes=TRUE)
file.copy(file.path(root,'R/13b_audit_charls_disease_routing.R'),file.path(out,'13b_audit_charls_disease_routing.R'),overwrite=TRUE)
print(summ);print(status[c('status','checks','raw_panel_rows')])
if(!all(qa$pass))stop('Raw disease routing audit failed; inspect QA before use.')
