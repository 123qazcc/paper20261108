#!/usr/bin/env Rscript
# Read-only study inputs; descriptive support, not event-time or intensity estimation.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(data.table);library(digest)})
setDTthreads(1L)
out<-"outputs/charls_external_extension_20260913"
dest<-file.path(out,"observational_support")
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
path<-file.path(out,"data/charls_model_bundle.rds")
bundle<-readRDS(path)
stopifnot(identical(bundle$status,"APPROVED_FOR_PILOT_AND_MI"),all(bundle$preparation_checks$pass))
checks<-list()
check<-function(name,pass){checks[[length(checks)+1L]]<<-data.table(check=name,pass=isTRUE(pass))}
check("signed_inputs_unchanged",all(vapply(names(bundle$input_hashes),function(f)digest(file=f,algo="sha256"),"")==bundle$input_hashes))
p<-copy(as.data.table(bundle$panel));setorder(p,id,time)
b<-copy(as.data.table(bundle$baseline));setorder(b,id)
check("unique_ID_and_observation_keys",!anyDuplicated(b$id)&&!anyDuplicated(p[,.(id,time)]))
check("baseline_panel_IDs_match",setequal(p$id,b$id))
check("all_panel_observations_are_type1",all(p$obstype==1L))
p[,visit_index:=seq_len(.N),by=id]
cc_ids<-b[complete.cases(b[,c("restless_bl",bundle$core_covariates),with=FALSE]),id]
labels<-c("无三种疾病","仅高血压","仅心脏病","高血压+心脏病","仅卒中",
  "高血压+卒中","心脏病+卒中","三病共存","已登记死亡")
dictionary<-data.table(model_state=1:9,bits=c(0:7,NA_integer_),label=labels)
fwrite(dictionary,file.path(dest,"state_dictionary.csv"))
compatible<-function(s)if(s%in%1:9)as.integer(s)else bundle$censor_states[[match(s,bundle$censor_codes)]]
reachable<-function(from,to)any(vapply(compatible(from),function(a)
  any(vapply(compatible(to),function(z)if(z==9L)TRUE else if(a==9L)FALSE else
    bitwAnd(a-1L,z-1L)==a-1L,TRUE)),TRUE))
pairs<-function(panel,bridge_partial=FALSE){
  z<-copy(panel)
  if(bridge_partial)z<-z[state%in%1:9]
  z[,`:=`(next_state=shift(state,type="lead"),next_time=shift(time,type="lead"),
    next_wave=shift(wave,type="lead"),next_index=shift(visit_index,type="lead")),by=id]
  z<-z[!is.na(next_state)]
  z[,`:=`(elapsed_years=next_time-time,partial_records_between=next_index-visit_index-1L,
    nominal_waves_between=pmax(0L,next_wave-wave-1L),
    exact_endpoints=state%in%1:9 & next_state%in%1:9)]
  z[,pair_class:=fcase(!exact_endpoints,"Partial-state endpoint(s)",next_state==9L,"Death registration interval",
    state==next_state,"Same fully observed living state",default="Living-state change")]
  z[exact_endpoints & next_state<=8L & state!=next_state,
    pair_class:=ifelse(bitwAnd(state-1L,next_state-1L)==state-1L & (next_state-state)%in%c(1L,2L,4L),
      "Single-disease increase endpoints","Multiple-disease increase endpoints")]
  z[,transition:=ifelse(exact_endpoints,paste(state-1L,next_state-1L,sep="-"),NA_character_)]
  z
}
# Small isolated assertions: a partial observation is not a healthy state;
# deleting it for a secondary descriptive pair cannot erase the gap flag.
toy<-data.table(id=rep(c("test1","test2"),each=3),time=rep(0:2,2),wave=rep(1:3,2),
  state=c(1L,101L,3L,2L,2L,6L),visit_index=rep(1:3,2))
aa<-pairs(toy);bb<-pairs(toy,TRUE)
check("pair_helper_preserves_partial_and_bridge_flags",nrow(aa)==4L&&sum(!aa$exact_endpoints)==2L&&
  nrow(bb)==3L&&bb[id=="test1",partial_records_between]==1L&&bb[id=="test1",transition]=="0-2")
rm(toy,aa,bb)
adj<-pairs(p);bridged<-pairs(p,TRUE)
check("adjacent_pair_count",nrow(adj)==nrow(p)-nrow(b))
check("strictly_positive_interval_lengths",all(adj$elapsed_years>0)&&all(bridged$elapsed_years>0))
check("adjacent_pairs_do_not_skip_recorded_visits",all(adj$partial_records_between==0L))
keys<-unique(adj[,.(state,next_state)])
keys[,compatible_path:=mapply(reachable,state,next_state)]
check("all_observations_have_compatible_forward_path",all(keys$compatible_path))
check("no_interval_starts_after_death",all(adj$state!=9L))
check("all_death_records_have_previous_observation",sum(adj$next_state==9L)==sum(p$state==9L))
check("fully_observed_death_not_assumed_exact_transition",all(adj[next_state==9L,pair_class]%in%
  c("Death registration interval","Partial-state endpoint(s)")))

make_summaries<-function(z,sample_name,pairing_name){
  allowed<-copy(as.data.table(bundle$edges))
  actual<-z[exact_endpoints==TRUE,.(observed_endpoint_pairs=.N,people=uniqueN(id),
    pairs_bridging_partial=sum(partial_records_between>0L),pairs_spanning_missing_waves=sum(nominal_waves_between>0L)),by=transition]
  allowed<-merge(allowed,actual,by="transition",all.x=TRUE,sort=FALSE)
  for(v in c("observed_endpoint_pairs","people","pairs_bridging_partial","pairs_spanning_missing_waves"))
    set(allowed,which(is.na(allowed[[v]])),v,0L)
  setorder(allowed,from_index,to_index)
  allowed[,`:=`(cohort=sample_name,pair_definition=pairing_name,
    interpretation=ifelse(to_index==9L,"Last seen living state to death registration; latent diseases may intervene",
      "Observed single-disease endpoint pattern; onset time unknown"))]
  full<-z[,.(pairs=.N,people=uniqueN(id)),by=pair_class]
  full[,`:=`(cohort=sample_name,pair_definition=pairing_name)]
  matrix<-z[exact_endpoints==TRUE,.(pairs=.N,people=uniqueN(id)),by=.(state,next_state)]
  matrix[,`:=`(cohort=sample_name,pair_definition=pairing_name)]
  list(edges=allowed,classes=full,matrix=matrix)
}
summaries<-list();j<-0L
for(sample_name in c("Primary_MI_cohort","Complete_case_pilot_cohort"))for(definition in c("All_adjacent_observations","Consecutive_fully_observed_endpoints")){
  z<-if(definition=="All_adjacent_observations")adj else bridged
  if(sample_name=="Complete_case_pilot_cohort")z<-z[id%in%cc_ids]
  j<-j+1L;summaries[[j]]<-make_summaries(z,sample_name,definition)
}
all_edges<-rbindlist(lapply(summaries,`[[`,"edges"))
classes<-rbindlist(lapply(summaries,`[[`,"classes"))
all_matrix<-rbindlist(lapply(summaries,`[[`,"matrix"))
fwrite(all_edges,file.path(dest,"allowed_edge_endpoint_support.csv"))
fwrite(classes,file.path(dest,"all_interval_classes.csv"))
fwrite(all_matrix,file.path(dest,"fully_observed_endpoint_matrix_long.csv"))
fwrite(keys,file.path(dest,"state_pair_path_compatibility.csv"))
targets<-c("0-2","0-4","1-3","1-5")
target_names<-c("无病→仅心脏病","无病→仅卒中","仅高血压→高血压+心脏病","仅高血压→高血压+卒中")
target_table<-data.table(transition=targets)
target_table[,target:=target_names]
for(sample_name in c("Primary_MI_cohort","Complete_case_pilot_cohort"))for(definition in c("All_adjacent_observations","Consecutive_fully_observed_endpoints")){
  z<-all_edges[cohort==sample_name & pair_definition==definition]
  newname<-paste(sample_name,definition,sep="__")
  target_table[,(newname):=z$observed_endpoint_pairs[match(transition,z$transition)]]
}
target_table[,origin_model_state:=c(1L,1L,2L,2L)]
target_table[,people_with_observed_origin_and_followup:=vapply(origin_model_state,function(s)uniqueN(adj[state==s,id]),1L)]
check("four_targets_retained_even_if_zero",nrow(target_table)==4L&&!anyNA(target_table))
fwrite(target_table,file.path(dest,"four_target_observation_support.csv"))
stratified<-rbindlist(lapply(seq_along(targets),function(k)rbindlist(lapply(0:3,function(level){
  ids_at_level<-b[restless_bl==level,id]
  z<-adj[id%in%ids_at_level]
  hit<-z[exact_endpoints==TRUE & transition==targets[k]]
  data.table(transition=targets[k],target=target_names[k],baseline_sleep_frequency_score=level,
    cohort_people=length(ids_at_level),people_observed_at_origin_with_followup=uniqueN(z[state==c(1L,1L,2L,2L)[k],id]),
    observed_endpoint_pairs=nrow(hit),people_with_endpoint_pattern=uniqueN(hit$id))
}))))
fwrite(stratified,file.path(dest,"four_target_support_by_baseline_sleep.csv"))
strata_totals<-stratified[,.(pairs=sum(observed_endpoint_pairs)),by=transition]
check("all_four_exposure_categories_retained",nrow(stratified)==16L&&all(stratified[,uniqueN(baseline_sleep_frequency_score),by=transition]$V1==4L))
check("exposure_strata_recover_full_target_counts",identical(strata_totals$pairs,
  target_table$Primary_MI_cohort__All_adjacent_observations))
check("no_repeated_target_pattern_per_person",all(stratified$observed_endpoint_pairs==stratified$people_with_endpoint_pattern))

fu<-p[,.(last_observation_time=max(time),observations=.N,partial_observations=sum(state>=100L),
  has_registered_death=any(state==9L)),by=id]
baseline_n<-nrow(readRDS(file.path(out,"data/charls_raw_W1_baseline_stage2.rds")))
flow<-data.table(stage=c("Raw_W1_descriptive_baseline","Primary_model_cohort","Complete_case_pilot",
  "Primary_people_with_partial_state_observations","Primary_people_with_death_registration"),
  people=c(baseline_n,nrow(b),length(cc_ids),
    sum(fu$partial_observations>0L),sum(fu$has_registered_death)))
fwrite(flow,file.path(dest,"cohort_counts.csv"))
followup<-data.table(people=nrow(b),records=nrow(p),adjacent_intervals=nrow(adj),
  observation_years=sum(fu$last_observation_time),median_followup=median(fu$last_observation_time),
  q25_followup=as.numeric(quantile(fu$last_observation_time,.25)),q75_followup=as.numeric(quantile(fu$last_observation_time,.75)),
  max_followup=max(fu$last_observation_time))
fwrite(followup,file.path(dest,"followup_summary.csv"))
wave<-p[,.(observations=.N,people=uniqueN(id),fully_observed_living=sum(state%in%1:8),
  partially_observed_living=sum(state>=100),death_registration=sum(state==9)),by=wave]
fwrite(wave,file.path(dest,"wave_observation_coverage.csv"))
state_coverage<-dictionary[,.(model_state,label)]
state_coverage[,`:=`(baseline_people=vapply(model_state,function(s)p[time==0 & state==s,uniqueN(id)],1L),
  people_ever_fully_observed=vapply(model_state,function(s)p[state==s,uniqueN(id)],1L),
  observed_origin_with_followup=vapply(model_state,function(s)adj[state==s,uniqueN(id)],1L))]
fwrite(state_coverage,file.path(dest,"state_observation_support.csv"))
check("baseline_state_counts_sum_to_cohort",sum(state_coverage$baseline_people)==nrow(b))
check("wave_coverage_sums_to_panel",sum(wave$observations)==nrow(p)&&
  all(wave$fully_observed_living+wave$partially_observed_living+wave$death_registration==wave$observations))
qa<-rbindlist(checks);fwrite(qa,file.path(dest,"support_checks.csv"));stopifnot(all(qa$pass))
fwrite(data.table(input=c("signed_model_bundle","support_script"),sha256=vapply(c(path,
  "R/15a_charls_observed_transition_support.R"),function(f)digest(file=f,algo="sha256"),"")),file.path(dest,"input_hashes.csv"))

lines<-c("# CHARLS 原始病史主队列：四条目标转移的观察支持", "",
  "本报告不拟合模型、不重插补、不改纳排、不估计q或HR。它更新原先包口径的历史计数，不能把两套口径拼在一起。", "",
  sprintf("主模型 %s 人；完整病例预检 %s 人；%s 条观测组成 %s 个相邻观察区间。",
    format(nrow(b),big.mark=","),format(length(cc_ids),big.mark=","),format(nrow(p),big.mark=","),format(nrow(adj),big.mark=",")),
  sprintf("累计观察窗口跨度 %.1f 人×年；个体跨度中位数 %.2f 年（四分位数 %.2f–%.2f 年）。死亡登记时点只是死亡已发生的上界，因此该数值不是准确存活人年或各状态的真实在险时间。",
    followup$observation_years,followup$median_followup,followup$q25_followup,followup$q75_followup),"",
  "| 目标转移 | 严格相邻观测、两端完整 | 允许跨部分状态后的完整两端 | 完整病例预检：严格相邻 |",
  "|---|---:|---:|---:|")
for(i in 1:4)lines<-c(lines,sprintf("| %s | %d | %d | %d |",target_names[i],
  target_table$Primary_MI_cohort__All_adjacent_observations[i],
  target_table$Primary_MI_cohort__Consecutive_fully_observed_endpoints[i],
  target_table$Complete_case_pilot_cohort__All_adjacent_observations[i]))
lines<-c(lines,"", "## 如何解释", "",
  "严格相邻指相邻的实际模型观测，不保证是相邻调查波；缺访波数另表记录。跨部分状态的计数仅作描述，没有从正式模型删除部分状态。",
  "数字是观察到的端点组合数，不是连续时间模型估计的事件数、年转移强度或风险比。没有假定访视之间准确何时发病。",
  "跨波多病增加仍作为多病跳跃保留；不按希望得到的病种顺序分配。死亡区间按上一次所见状态列示，区间内可能先新增疾病，因此不把它当作该状态直接死亡的确认事件。",
  "部分状态不能视为健康，缺访不能视为存活无病。所有观察对的前向兼容路径另行核查。",
  "两端完整和桥接计数是同一批参与者的不同摘要，不能相加。各状态观察人数也可能重复，不能相加为总人数。",
  "另保存四档基线睡眠频率下的全部16个目标转移单元及对应起点观察人数；未换算粗率、未检验或筛选组别。",
  "完整病例与主队列的差别由已冻结协变量完整性决定；不按目标边计数或显著性调整样本。计数较少提示可能精度不足，但不决定删除目标边。",
  "没有进行组间检验；不以事件较多或较少推断睡眠关联强弱，更不能据此宣布心脏病效应大于卒中。正式差异需等待同模型协方差的直接对比。", "",
  sprintf("核查：%d/%d项通过。HRS/ELSA文件及在算的CHARLS模型不改动。",sum(qa$pass),nrow(qa)))
writeLines(lines,file.path(dest,"CHARLS_目标转移观察支持.md"),useBytes=TRUE)
print(target_table);print(flow);print(classes[cohort=="Primary_MI_cohort" & pair_definition=="All_adjacent_observations"])
message("Observed support audit complete; no association or transition-intensity estimates produced.")
