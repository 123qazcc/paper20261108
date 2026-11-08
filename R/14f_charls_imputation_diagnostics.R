#!/usr/bin/env Rscript
# Diagnostic graphics only. No fitting, re-imputation or participant filtering.
options(stringsAsFactors=FALSE)
diagnostic_lib<-".runtime/charls_diagnostics_R"
if(dir.exists(diagnostic_lib)).libPaths(c(normalizePath(diagnostic_lib),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(jsonlite);library(digest)})
setDTthreads(1)
args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args))args[1] else "missingness"
stopifnot(mode %in% c("missingness","complete"))
out <- "outputs/charls_external_extension_20260913"
figdir <- file.path(out,"figures/imputation_qc")
srcdir <- file.path(out,"figure_source_data/imputation_qc")
qcdir <- file.path(out,"model_qc/imputation_diagnostics")
for(p in c(figdir,srcdir,qcdir))dir.create(p,recursive=TRUE,showWarnings=FALSE)
source("R/helpers/charls_cairo_panel_alignment.R")
local_qa<-normalizePath(".runtime/figure_qa",winslash="/",mustWork=TRUE)
Sys.setenv(PYTHONPATH=paste(c(local_qa,Sys.getenv("PYTHONPATH")),collapse=.Platform$path.sep))
python <- "python"
skill_scripts <- "private/render_tools/nature-figure/scripts"
need <- c("ggplot2","jsonlite","svglite","ragg","mice")
stopifnot(all(vapply(need,requireNamespace,TRUE,quietly=TRUE)))
theme_set(theme_classic(base_size=7.5,base_family="Arial")+
  theme(axis.text=element_text(size=7,colour="#303030"),
    axis.title=element_text(size=7.5),axis.line=element_line(linewidth=.3),
    axis.ticks=element_line(linewidth=.3),strip.background=element_blank(),
    strip.text=element_text(size=7.5,face="bold",hjust=0),
    plot.title=element_text(size=9,face="bold"),plot.subtitle=element_text(size=7.5),
    plot.caption=element_text(size=7,hjust=0),panel.spacing=grid::unit(4,"mm"),
    plot.margin=margin(7,9,7,7),legend.position="none"))
run_audit <- function(script,arguments,logfile,blocking=TRUE) {
  output<-suppressWarnings(system2(python,c(shQuote(file.path(skill_scripts,script)),arguments),stdout=TRUE,stderr=TRUE))
  status<-attr(output,"status");if(is.null(status))status<-0L
  writeLines(output,logfile)
  if(status!=0L&&blocking)stop("Figure QA did not pass; inspect ",logfile)
  invisible(status)
}
export_plot <- function(p,name,height_mm) {
  stem<-file.path(figdir,name);auditstem<-file.path(qcdir,name);width_mm<-183;w<-width_mm/25.4;h<-height_mm/25.4
  require_charls_panel_alignment(p,auditstem,w,h,family="Arial")
  svglite::svglite(paste0(stem,".svg"),width=w,height=h);print(p);dev.off()
  grDevices::cairo_pdf(paste0(stem,".pdf"),width=w,height=h,family="Arial")
  print(p);dev.off()
  ragg::agg_tiff(paste0(stem,".tiff"),width=w,height=h,units="in",res=600,compression="lzw")
  print(p);dev.off()
  ragg::agg_png(paste0(stem,".png"),width=w,height=h,units="in",res=600)
  print(p);dev.off()
  text_status<-run_audit("audit_pdf_text.py",c(shQuote(paste0(stem,".pdf")),"--min-pt","5","--json"),
    paste0(auditstem,".text-audit.json"),blocking=FALSE)
  run_audit("audit_figure_collisions.py",c(shQuote(paste0(stem,".pdf")),"--json-out",
    shQuote(paste0(auditstem,".collision-audit.json")),"--overlay-pdf",shQuote(paste0(auditstem,".collision-audit.pdf"))),
    paste0(auditstem,".collision-audit.log"),blocking=FALSE)
  # Raw skill checks are preserved, not relabelled as pass. The adapted check
  # accounts for Cairo matrices and separately positioned lines in one paint operation.
  status<-system2(python,c(shQuote("scripts/audit_charls_cairo_geometry.py"),shQuote(paste0(stem,".pdf")),
    "--json-out",shQuote(paste0(auditstem,".cairo-rendered-audit.json"))))
  if(status!=0L)stop("Effective rendered geometry did not pass: ",auditstem)
  invisible(stem)
}
labels <- c(age5="Age",female="Sex",education4="Education (four categories)",married="Marital status",
  rural="Rural residence",wealth_asinh="Non-housing financial wealth",bmi5="Body mass index",
  smoker="Current smoking",drinker="Alcohol use",diabetes="Diabetes",lung="Chronic lung disease",pain="Physical pain")
miss <- fread(file.path(out,"model_qc/imputation_missingness.csv"))[role=="baseline_covariate"]
stopifnot(nrow(miss)==12L,setequal(miss$variable,names(labels)),length(unique(miss$n))==1L,
  all(miss$missing>=0 & miss$missing<=miss$n))
miss[,percent:=100*missing/n]
setorder(miss,-missing,variable)
miss[,display:=factor(labels[variable],levels=rev(labels[variable]))]
miss[,label:=sprintf("%s (%.2f%%)",format(missing,big.mark=",",trim=TRUE),percent)]
fwrite(miss,file.path(srcdir,"missingness_all_core_variables.csv"))
p<-ggplot(miss,aes(percent,display))+
  geom_segment(aes(x=0,xend=percent,yend=display),linewidth=.55,colour="#B9C7D0")+
  geom_point(size=2,colour="#27658C")+
  geom_text(aes(x=percent+.45,label=label),hjust=0,size=7/ggplot2::.pt,family="Arial")+
  scale_x_continuous(limits=c(0,max(miss$percent)+5.5),breaks=seq(0,20,5),expand=expansion(mult=c(.01,.01)))+
  labs(x="Participants with missing baseline values (%)",y=NULL,
    title="CHARLS: missingness before covariate imputation",
    subtitle=sprintf("Primary model cohort: %s participants; all 12 core covariates shown",format(miss$n[1],big.mark=",")),
    caption="Observed missing counts and percentages; no confidence intervals or hypothesis tests.")
export_plot(p,"charls_core_covariate_missingness",116)
if(mode=="missingness"){
  message("Missingness figure ready. Full MICE diagnostics intentionally wait for completed imputation.")
  quit(save="no",status=0)
}

# The complete mode is read-only with respect to all analysis inputs.
imp_path<-file.path(out,"data/charls_baseline_mice20.rds")
bun_path<-file.path(out,"data/charls_model_bundle.rds")
if(!all(file.exists(c(imp_path,bun_path))))stop("NOT READY: completed mids and model bundle are required.")
imp<-readRDS(imp_path);bun<-readRDS(bun_path)
stopifnot(imp$m==20L,imp$iteration==20L,nrow(imp$data)==miss$n[1],all(bun$preparation_checks$pass))
variables<-names(imp$method)[imp$method!=""]
stopifnot(setequal(variables,miss[missing>0,variable]))
binary<-intersect(variables,c("female","married","rural","smoker","drinker","diabetes","lung","pain"))
for(v in binary)stopifnot(identical(levels(imp$data[[v]]),c("0","1")))
trace<-rbindlist(lapply(variables,function(v){
  mn<-imp$chainMean[v,,,drop=TRUE];va<-imp$chainVar[v,,,drop=TRUE]
  stopifnot(identical(dim(mn),c(20L,20L)),identical(dim(va),c(20L,20L)))
  z<-CJ(iteration=1:20,chain=1:20)
  z[,`:=`(variable=v,mean_native=mn[cbind(iteration,chain)],variance_native=va[cbind(iteration,chain)])]
  # mice's factor summaries use 1-based integer level indices.
  z[,mean_display:=mean_native-if(v %in% binary)1 else 0]
  z[,sd_display:=sqrt(variance_native)]
  z
}))
fwrite(trace,file.path(srcdir,"chain_trace_all_iterations_all_chains.csv"))
trace_names<-paste0(letters[seq_along(variables)],"  ",labels[variables],
  ifelse(variables=="education4"," [code 1-4]",ifelse(variables%in%binary," [proportion]"," [model scale]")))
trace[,panel:=factor(variable,levels=variables,labels=trace_names)]
for(metric in c("mean_display","sd_display")){
  shown<-copy(trace)
  if(metric=="sd_display")shown[,panel:=factor(variable,levels=variables,
    labels=gsub("[proportion]","[binary 0/1]",trace_names,fixed=TRUE))]
  z<-shown[is.finite(get(metric))]
  absent<-shown[,.(available=any(is.finite(get(metric)))),by=panel][available==FALSE]
  if(nrow(absent))absent[,`:=`(iteration=10,value=0,label="Undefined: only one missing value")]
  p<-ggplot(z,aes(iteration,.data[[metric]],group=chain))+
    geom_line(colour="#27658C",alpha=.45,linewidth=.25)+
    facet_wrap(~panel,ncol=2,scales="free_y",drop=FALSE)+
    scale_x_continuous(breaks=c(1,5,10,15,20),limits=c(1,20),expand=expansion(mult=.02))+
    labs(x="MICE iteration",y=if(metric=="mean_display")"Mean of imputed values" else "SD of imputed values",
      title=paste("CHARLS: imputation-chain",if(metric=="mean_display")"means" else "standard deviations"),
      subtitle="All 20 chains and all 20 iterations; only originally missing cells contribute",
      caption="Trace summaries are diagnostics, not confidence intervals or proof of convergence.\nEducation codes are categorical summaries, not a continuous exposure model.")
  if(nrow(absent))p<-p+geom_text(data=absent,aes(x=iteration,y=value,label=label),inherit.aes=FALSE,
    size=6.5/ggplot2::.pt,family="Arial")
  if(metric=="sd_display")p<-p+scale_y_continuous(limits=function(x)c(max(0,x[1]),max(.1,x[2])))
  export_plot(p,paste0("charls_mice20_",metric),208)
}

# Export every observed/imputed distribution; do not pool draws as extra people.
summaries<-list();frequencies<-list();checks<-list()
addcheck<-function(n,v)checks[[length(checks)+1L]]<<-data.table(check=n,pass=isTRUE(v))
num<-function(x)if(is.factor(x))as.numeric(as.character(x))else as.numeric(x)
for(v in variables){
  obs<-num(imp$data[[v]][!is.na(imp$data[[v]])]);mat<-as.matrix(imp$imp[[v]])
  for(i in 0:20){
    x<-if(i==0L)obs else num(mat[,i]);q<-quantile(x,c(0,.25,.5,.75,1),names=FALSE)
    summaries[[length(summaries)+1L]]<-data.table(variable=v,imputation=i,
      group=if(i==0)"Observed cells" else "Imputed cells",n=length(x),mean=mean(x),
      sd=sd(x),min=q[1],q25=q[2],median=q[3],q75=q[4],max=q[5])
    if(v %in% c(binary,"education4")){
      valid<-if(v=="education4")1:4 else 0:1
      addcheck(paste(v,i,"categorical_support",sep="_"),all(x%in%valid))
      frequencies[[length(frequencies)+1L]]<-data.table(variable=v,imputation=i,
        value=valid,n=as.integer(table(factor(x,levels=valid))),denominator=length(x))
    }else if(i>0L){
      addcheck(paste(v,i,"PMM_observed_donor_support",sep="_"),all(x%in%obs))
    }
  }
}
fwrite(rbindlist(summaries),file.path(qcdir,"observed_imputed_distribution_summary.csv"))
fwrite(rbindlist(frequencies),file.path(qcdir,"observed_imputed_category_counts.csv"))
addcheck("all_trace_means_finite",all(is.finite(trace$mean_display)))
addcheck("nonnegative_finite_variances_where_defined",all(trace$variance_native[is.finite(trace$variance_native)]>=0))
fwrite(rbindlist(checks),file.path(qcdir,"support_checks.csv"))
warnings<-character()
diagnostics<-rbindlist(lapply(c("mean","sd"),function(parameter){
  diagnostic_parameter<-parameter
  ans<-withCallingHandlers(mice::convergence(imp,diagnostic="ac",parameter=parameter),warning=function(w){
    warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  ans<-as.data.table(ans);ans<-ans[vrb %in% variables];ans[,parameter:=diagnostic_parameter]
  # The installed mice PSRF branch requires absent rstan. Use the documented
  # standalone posterior rank-normalized/folded split-Rhat implementation.
  stopifnot(requireNamespace("posterior",quietly=TRUE))
  ans[,rhat_posterior:=vapply(seq_len(.N),function(i){
    a<-if(diagnostic_parameter=="mean")imp$chainMean[as.character(vrb[i]),,,drop=TRUE]else sqrt(imp$chainVar[as.character(vrb[i]),,,drop=TRUE])
    x<-a[seq_len(.it[i]),,drop=FALSE]
    if(nrow(x)<4L)return(NA_real_)
    posterior::rhat(x)
  },0.0)]
  ans[,rhat_method:="posterior rank-normalized and folded split-Rhat; not mice::convergence(all)"]
  ans
}),fill=TRUE)
fwrite(diagnostics,file.path(qcdir,"mice_convergence_all_iterations.csv"))
fwrite(data.table(package=c("mice","posterior"),version=vapply(c("mice","posterior"),function(p)as.character(packageVersion(p)),"")),
  file.path(qcdir,"diagnostic_package_versions.csv"))
fwrite(data.table(warning=unique(warnings)),file.path(qcdir,"mice_convergence_warnings.csv"))
fwrite(data.table(input=c("completed_mids","model_bundle","diagnostics_script"),
  sha256=vapply(c(imp_path,bun_path,"R/14f_charls_imputation_diagnostics.R"),function(f)digest(file=f,algo="sha256"),"")),
  file.path(qcdir,"input_hashes.csv"))
writeLines(c("# CHARLS 多重插补诊断：待人工审阅", "",
  paste0("支持集与数值检查：",sum(rbindlist(checks)$pass),"/",length(checks)," 通过。"),
  "20 份插补、20 轮轨迹完整保存；参与者是独立分析单位。连续变量图使用模型中的缩放单位。",
  "教育为多类别插补；均值只是类别编码摘要。二元均值转换为0/1比例。",
  "一个缺失值的链内标准差无定义；极少缺失或恒定轨迹的自相关/PSRF可能为NA或不稳定。",
  "当前 mice 的 PSRF 依赖未安装的 rstan；自相关仍由 mice 计算，Rhat 改用项目隔离库中的 posterior::rhat。",
  "缺失与观测分布不要求一致；支持集通过不证明条件模型、MAR或多状态拟合有效。",
  "请联合查看轨迹、MICE loggedEvents、均值/标准差诊断的迭代趋势与模型收敛。",
  "本脚本不根据单一PSRF阈值宣布收敛，不修改插补或模型数据。",
  "方法参考：https://amices.org/mice/reference/convergence.html ; https://mc-stan.org/posterior/reference/rhat.html"),file.path(qcdir,"diagnostic_review_pending.md"))
stopifnot(all(rbindlist(checks)$pass))
message("Diagnostic artifacts written; substantive trace and warning review is still required.")
