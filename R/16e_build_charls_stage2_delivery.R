#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,warn=1)
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L,args[[1]]=="build")
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(svglite);library(ragg);library(jsonlite);library(digest)})
setDTthreads(1L)
root<-"outputs/charls_cesd9_duration_20260914";dest<-file.path(root,"stage2_delivery");figdir<-file.path(dest,"figures")
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
sha<-function(p)digest(file=p,algo="sha256")
review<-fromJSON(file.path(root,"model_qc/stage2_independent_pool_review.json"),simplifyVector=TRUE)
ps<-fromJSON(file.path(root,"model_qc/pooled_status.json"),simplifyVector=TRUE)
stopifnot(identical(ps$status,"POOLED_AWAITING_INDEPENDENT_REVIEW_AND_FIGURES"),ps$accepted_fits==40L,
  identical(review$status,"COMPUTATIONAL_CHECKS_PASS_REQUIRES_SUBSTANTIVE_REVIEW"),review$checks_passed==review$checks_total)
four<-fread(file.path(root,"results/four_target_MI20.csv"));contrast<-fread(file.path(root,"results/heart_stroke_contrasts_MI20.csv"));joint<-fread(file.path(root,"results/joint_tests_MI20.csv"))
stopifnot(nrow(four)==8L,nrow(contrast)==4L,nrow(joint)==4L,all(is.finite(four$HR)),all(four$lower>0),all(four$upper>=four$HR),all(contrast$lower>0))
write_csv<-function(x,p){fwrite(x,p,bom=TRUE);stopifnot(file.exists(p),file.info(p)$size>0)}
write_csv(four,file.path(dest,"table1_four_target_paths_CORE_vs_BOTH.csv"));write_csv(contrast,file.path(dest,"table2_heart_stroke_direct_CORE_vs_BOTH.csv"));write_csv(joint,file.path(dest,"table3_joint_tests_CORE_vs_BOTH.csv"))
path_labels<-c("0–2: no disease → heart disease only","0–4: no disease → stroke only","1–3: hypertension only → hypertension + heart disease","1–5: hypertension only → hypertension + stroke")
four[,path_label:=factor(path_labels[match(transition,c("0-2","0-4","1-3","1-5"))],levels=rev(path_labels))]
four[,model_label:=factor(model,levels=c("CORE","BOTH"),labels=c("CORE","BOTH: CESD-9 + duration"))]
four[,display:=factor(paste0(as.character(path_label),"  |  ",as.character(model_label)),levels=rev(unique(paste0(as.character(path_label),"  |  ",as.character(model_label)))))]
contrast[,origin_label:=factor(ifelse(origin=="No disease","No target disease at origin","Hypertension only at origin"),levels=rev(c("No target disease at origin","Hypertension only at origin")))]
contrast[,model_label:=factor(model,levels=c("CORE","BOTH"),labels=c("CORE","BOTH: CESD-9 + duration"))]
contrast[,display:=factor(paste0(as.character(origin_label),"  |  ",as.character(model_label)),levels=rev(unique(paste0(as.character(origin_label),"  |  ",as.character(model_label)))))]
cols<-c("CORE"="#4D5966","BOTH: CESD-9 + duration"="#B44B3E")
theme_pub<-theme_classic(base_family="Arial",base_size=8)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=7.2,colour="#1F2933",margin=margin(r=6)),axis.text.x=element_text(size=7.2,colour="#1F2933"),axis.title.x=element_text(size=8,margin=margin(t=7)),plot.title=element_text(size=10,face="bold",hjust=0),plot.subtitle=element_text(size=7.5,colour="#4B5563",hjust=0,margin=margin(b=6)),plot.caption=element_text(size=6.5,colour="#4B5563",hjust=0,lineheight=1.15,margin=margin(t=8)),legend.position="top",legend.title=element_blank(),legend.text=element_text(size=7.2),panel.grid.major.x=element_line(colour="#E5E7EB",linewidth=.28),panel.grid.minor=element_blank(),plot.margin=margin(8,10,8,8))
forest_plot<-function(d,title,subtitle,xlabel,caption){
  lim<-range(log(c(d$lower,d$upper,1)));pad<-max(diff(lim)*.12,.06);lims<-exp(lim+c(-pad,pad));br<-sort(unique(c(1,exp(pretty(log(lims),n=5L)))));br<-br[br>=lims[1]&br<=lims[2]]
  ggplot(d,aes(x=HR,y=display,colour=model_label))+geom_vline(xintercept=1,linetype="dashed",linewidth=.35,colour="#6B7280")+geom_segment(aes(x=lower,xend=upper,yend=display),linewidth=.65)+geom_point(shape=15,size=2.25)+scale_colour_manual(values=cols,drop=FALSE)+scale_x_log10(limits=lims,breaks=br,labels=function(x)formatC(x,digits=3,format="g"))+scale_y_discrete(drop=FALSE)+labs(title=title,subtitle=subtitle,x=xlabel,y=NULL,caption=caption)+theme_pub
}
sub<-"CHARLS; n=12,036; 53,786 panel records; 20 MI datasets/model"
cap<-paste(strwrap("Squares/bars: Rubin-pooled HR and pointwise 95% finite-MI t interval. Exposure: one higher W1 restless-sleep category (0–3). CORE=core covariates; BOTH=CORE plus CESD-9 and sleep-duration spline. Same panel and 20 completed datasets; adjustment change is descriptive, not mediation.",width=68),collapse="\n")
p1<-forest_plot(four,"Four target transitions: CORE vs BOTH",sub,"Transition HR per one category higher restless-sleep frequency",cap)
p2<-forest_plot(contrast,"Heart versus stroke: CORE vs BOTH",sub,"Heart-transition HR / stroke-transition HR",paste(strwrap("Ratio uses complete within-model covariance; it is not an incidence-rate ratio. A ratio above 1 means the sleep association is larger for the heart transition at the same origin. BH(2) values are in the companion table.",width=80),collapse="\n"))
export_one<-function(p,stem,wmm=183,hmm=112){wi<-wmm/25.4;hi<-hmm/25.4;grDevices::cairo_pdf(paste0(stem,".pdf"),width=wi,height=hi,family="Arial");print(p);dev.off();svglite::svglite(paste0(stem,".svg"),width=wi,height=hi);print(p);dev.off();ragg::agg_png(paste0(stem,".png"),width=wi,height=hi,units="in",res=600);print(p);dev.off();ragg::agg_tiff(paste0(stem,".tiff"),width=wi,height=hi,units="in",res=600,compression="lzw");print(p);dev.off();stopifnot(all(file.info(paste0(stem,c(".pdf",".svg",".png",".tiff")))$size>0));write_json(list(verdict="NOT APPLICABLE",reason="Single plot area; no comparable multi-panel rectangles.",backend="R ggplot2",figure=basename(stem),pdf_sha256=sha(paste0(stem,".pdf"))),paste0(stem,".alignment.json"),pretty=TRUE,auto_unbox=TRUE)}
s1<-file.path(figdir,"figure1_four_targets_CORE_vs_BOTH");s2<-file.path(figdir,"figure2_heart_stroke_CORE_vs_BOTH");export_one(p1,s1,183,112);export_one(p2,s2,183,88)
write_csv(four,file.path(figdir,"figure1_four_targets_CORE_vs_BOTH.source_data.csv"));write_csv(contrast,file.path(figdir,"figure2_heart_stroke_CORE_vs_BOTH.source_data.csv"))
fmtp<-function(x)formatC(x,digits=4,format="g")
writeLines(c("# Figure 1 | Four target transition associations","","All four prespecified transitions are shown for CORE and BOTH. Squares and bars are Rubin-pooled HRs and pointwise 95% finite-MI t intervals. This is an association analysis; adjustment changes are descriptive, not mediation estimates.","",paste0("Four-target global Wald P: CORE = ",fmtp(joint[model=="CORE"&test=="Four target coefficients = 0",p]),"; BOTH = ",fmtp(joint[model=="BOTH"&test=="Four target coefficients = 0",p]),".")),file.path(figdir,"figure1_four_targets_CORE_vs_BOTH.caption.md"),useBytes=TRUE)
writeLines(c("# Figure 2 | Direct heart-versus-stroke comparison","","The plotted quantity is the heart-transition HR divided by the stroke-transition HR at the same origin. It is not an incidence-rate ratio. BH(2)-adjusted values and the joint two-contrast Wald test are reported in the companion table. A difference in significance between two individual transitions is not itself evidence of a disease-specific difference."),file.path(figdir,"figure2_heart_stroke_CORE_vs_BOTH.caption.md"),useBytes=TRUE)
write_json(list(status="DRAFT_REQUIRES_HUMAN_FIGURE_AND_SCIENTIFIC_REVIEW",backend="R",independent_review=review,figures=c(s1,s2),automatic_final_approval=FALSE,human_figure_review="PENDING",human_scientific_review="PENDING"),file.path(dest,"delivery_status.json"),pretty=TRUE,auto_unbox=TRUE)
writeLines(c("# CHARLS stage 2 delivery review","","Status: DRAFT_REQUIRES_HUMAN_FIGURE_AND_SCIENTIFIC_REVIEW","","40/40 CORE and BOTH fits were accepted; 227/227 independent computational checks passed. Tables and two R figures were generated from pooled outputs.","","The numerical audit does not prove MAR, MICE convergence, Markov/time-homogeneity, no unmeasured confounding or causal mediation. Inspect final PDFs/PNGs at final size and run source/PDF/collision QA before manuscript use."),file.path(dest,"delivery_review.md"),useBytes=TRUE)
message("STAGE2_DELIVERY_DRAFT_READY: two figures; three tables; 227/227 numerical checks")
