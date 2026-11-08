#!/usr/bin/env Rscript
# Targeted source recheck only. No event recoding or model changes.
suppressPackageStartupMessages({library(data.table);library(haven);library(digest)})
setDTthreads(1L)
out<-"outputs/charls_external_extension_20260913"
dest<-file.path(out,"secondary_readiness")
loc<-readRDS(file.path(dest,"internal/living_status_conflict_locators.rds"))
loc<-as.data.table(loc)
stopifnot(nrow(loc)==2L,all(loc$wave==4L))
rawroot<-"private/source/data_raw/1_CHARLS_China"
paths<-c(sample=file.path(rawroot,"2018charls/Sample_Infor.dta"),
  cognition=file.path(rawroot,"2018charls/Cognition.dta"),
  exit=file.path(rawroot,"2020charls/Exit_Module.dta"))
stopifnot(all(file.exists(paths)))
read_raw<-function(path,columns){
  tmp<-tempfile(fileext=".dta")
  stopifnot(file.copy(path,tmp),digest(file=path,algo="sha256")==digest(file=tmp,algo="sha256"))
  on.exit(unlink(tmp),add=TRUE)
  z<-as.data.table(haven::read_dta(tmp,col_select=all_of(columns)))
  z[,id:=if(is.numeric(ID))sprintf("%012.0f",as.numeric(ID))else trimws(as.character(ID))]
  z[id%in%loc$id]
}
s<-read_raw(paths["sample"],c("ID","died","iyear","imonth"))
cog<-read_raw(paths["cognition"],c("ID","dc015"))
e<-read_raw(paths["exit"],c("ID","exb001_1","exb001_2","exb002"))
z<-merge(loc,s[,.(id,source_died=as.numeric(died),source_iyear=as.numeric(iyear),source_imonth=as.numeric(imonth))],by="id",all.x=TRUE)
z<-merge(z,cog[,.(id,source_sleep_code=as.numeric(dc015))],by="id",all.x=TRUE)
z<-merge(z,e[,.(id,source_exit_year=as.numeric(exb001_1),source_exit_month=as.numeric(exb001_2),source_calendar=as.numeric(exb002))],by="id",all.x=TRUE)
z[,source_sleep_score:=ifelse(source_sleep_code%in%1:4,source_sleep_code-1,NA_real_)]
bundle<-readRDS(file.path(out,"data/charls_model_bundle.rds"))
p<-as.data.table(bundle$panel)
z[,no_later_model_observation:=vapply(id,function(person){a<-p[id==person];max(a$wave)==4L},TRUE)]
checks<-data.table(check=c("two_source_rows_per_module","raw_tracking_alive_flag_and_month",
  "death_dates_precede_W4_interview","sleep_value_and_missingness_match_raw",
  "both_death_endpoints_already_quarantined","both_disputed_visits_are_terminal","signed_inputs_unchanged"),
  pass=c(nrow(s)==2L&&nrow(cog)==2L&&nrow(e)==2L,
    all(z$source_died==0 & z$source_iyear==2018 & z$source_imonth==7),
    all(z$source_exit_year+z$source_exit_month/12 < z$source_iyear+(z$source_imonth-1)/12),
    identical(z$source_sleep_score,as.numeric(z$restless_score_observed)),
    all(!z$death_usable & grepl("solar_death_before_alive",z$conflicts)),
    all(z$no_later_model_observation),
    all(vapply(names(bundle$input_hashes),function(f)digest(file=f,algo="sha256"),"")==bundle$input_hashes)))
fwrite(checks,file.path(dest,"independent_terminal_source_checks.csv"))
stopifnot(all(checks$pass))
saveRDS(z,file.path(dest,"internal/independent_terminal_source_review.rds"))
fwrite(data.table(source=c(names(paths),"review_script"),sha256=vapply(c(paths,
  "R/15c_verify_charls_disputed_terminal_visits.R"),function(f)digest(file=f,algo="sha256"),"")),
  file.path(dest,"independent_terminal_source_hashes.csv"))
writeLines(c("# 两条争议末次访谈：独立原始字段复核", "",
  "重新读取原始W4追踪、认知和W5退出文件；源文件散列及两条内部定位记录另存。公开说明不披露参与者ID。",
  "7项字段一致性复核通过：两条W4原始追踪均标存活并有访谈年月，W5所填死亡年月却早于该访谈；睡眠频率的1条有效值及1条缺失与原始dc015一致。",
  "两个死亡终点均已被既定死亡审计隔离；这两条W4行均为末次模型观察，不作任何后续区间起点。因此争议睡眠值没有进入已枚举时间更新区间的起点暴露，但作为两个结局端点仍有真实性限制。",
  "字段复现不证明哪个来源反映真实生存状况；本次没有依据单一来源改判死亡、删除人或追溯修改疾病史。",
  "REVIEW_REQUIRED继续保留。W2动态有后续信息的样本人数不受这两条行是否计入影响；此人数稳定性不表示拟合结果一定不受影响。",
  "正式结果交付时应明确这一小范围来源冲突，并按原死亡审计已经提出的整人排除/较早截尾敏感性检查其影响。任何变体另建分支，不覆盖当前主模型。"),
  file.path(dest,"CHARLS_争议末次访谈独立复核.md"),useBytes=TRUE)
print(checks)
