#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
out<-"outputs/charls_external_extension_20260913"
b<-as.data.table(readRDS(file.path(out,"data/charls_raw_W1_baseline_stage2.rds")))
groups<-c("总体","＜1天/周","1–2天/周","3–4天/周","5–7天/周")
continuous<-c(age="年龄（岁）",bmi="BMI",sleep_hours_observed="夜间睡眠时长（小时）",cesd9_no_sleep="CESD-9（不含睡眠题）",
 wealth_asinh="非住房金融财富：asinh(value/10000)")
binary<-c(female="女性",smoker="当前吸烟（原始W1）",drinker="过去一年饮酒（原始W1）",diabetes="糖尿病（原始W1）",
 lung="慢性肺病（原始W1）",pain="经常身体疼痛（原始W1）",married="已婚/同居",rural="农村居住")
state_names<-c("无高血压/心脏病/卒中","单纯高血压","单纯心脏病","高血压＋心脏病","单纯卒中","高血压＋卒中","心脏病＋卒中","三病共存")
education<-c("未完成小学/私塾","小学","初中","高中及以上")
rows<-list()
for(g in 0:4){
 z<-if(g==0)b else b[restless_bl==g-1]
 add<-function(key,label,summary,available,missing){rows[[length(rows)+1]]<<-data.table(variable=key,label=label,group=groups[g+1],
   summary=summary,n_group=nrow(z),n_available=available,n_missing=missing)}
 add("n","人数",as.character(nrow(z)),nrow(z),0L)
 for(v in names(continuous)){x<-z[[v]];ok<-x[!is.na(x)];add(v,continuous[[v]],sprintf("%.2f (%.2f)",mean(ok),sd(ok)),length(ok),sum(is.na(x)))}
 for(v in names(binary)){x<-z[[v]];den<-sum(!is.na(x));n<-sum(x==1,na.rm=TRUE)
  add(v,binary[[v]],sprintf("%d/%d (%.1f%%)",n,den,100*n/den),den,sum(is.na(x)))}
 for(s in 0:7)add(paste0("state",s),state_names[s+1],sprintf("%d (%.1f%%)",sum(z$bits==s),100*mean(z$bits==s)),nrow(z),0)
 for(e in 1:4){x<-as.integer(z$education4);n<-sum(x==e,na.rm=TRUE);den<-sum(!is.na(x))
  add(paste0("education",e),paste0("教育：",education[e]),sprintf("%d/%d (%.1f%%)",n,den,100*n/den),den,sum(is.na(x)))}
}
t<-rbindlist(rows);fwrite(t,file.path(out,"tables/baseline_raw_W1_stage2_with_denominators.csv"),bom=TRUE)
w<-dcast(t,variable+label~group,value.var="summary");w<-w[match(unique(t$variable),variable)];setcolorder(w,c("variable","label",groups))
fwrite(w,file.path(out,"tables/baseline_raw_W1_stage2.csv"),bom=TRUE)
md<-c("# CHARLS 原始W1口径基线表（第二阶段更新）","",
 "N=12,319；这是W1描述队列，不按未来回访筛选。拟合队列另要求后续有信息观察，当前12,036人。",
 "",
 "疾病使用当时原始回答，不利用后续争议回改基线。糖尿病、肺病、吸烟、饮酒、疼痛亦采用原始W1回答；教育来自W1 bd001。",
 "第一阶段按包数据生成的基线表保留为历史核查版，疾病/教育/协变量分布不能与此版本混用。",
 "",
 "连续变量为均值（标准差），分类变量注明非缺失分母。未加调查权重，不作基线p值检验；缺失分母见CSV长表。",
 "",paste0("| 变量 | ",paste(groups,collapse=" | ")," |"),"|---|---|---|---|---|---|")
for(i in seq_len(nrow(w)))md<-c(md,paste0("| ",w$label[i]," | ",paste(unlist(w[i,..groups]),collapse=" | ")," |"))
writeLines(md,file.path(out,"reports/CHARLS_原始W1口径更新基线表.md"),useBytes=TRUE)
old<-as.data.table(readRDS(file.path(out,"data/charls_w1_descriptive_baseline.rds")))
j<-merge(old[,.(id,package_state=disease_state_bits)],b[,.(id,raw_state=bits)],by="id",all=TRUE)
fwrite(j[,.N,by=.(package_state,raw_state)],file.path(out,"tables/baseline_package_to_raw_state_crosswalk.csv"))
cat("Raw baseline table updated. N=",nrow(b),"; state differences=",sum(j$package_state!=j$raw_state,na.rm=TRUE),"\n")
