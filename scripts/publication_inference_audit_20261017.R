#!/usr/bin/env Rscript
# Independently re-extract/re-pool saved fits and audit existing correction families.
# No fitting, no replacement of historical output.
Sys.setlocale('LC_CTYPE','English_United States.utf8')
options(stringsAsFactors=FALSE,digits=17)
suppressPackageStartupMessages({library(data.table);library(msm);library(jsonlite);library(digest)})
setDTthreads(1L)
out<-'outputs/publication_methods_audit_20261017';dir.create(out,recursive=TRUE,showWarnings=FALSE)
F<-'outputs/hrs_elsa_multistate_final_20260909';N<-'outputs/hrs_elsa_multistate_nextsteps_20260909'
C1<-'outputs/charls_external_extension_20260913';C2<-'outputs/charls_cesd9_duration_20260914';C3<-'outputs/charls_dynamic_sleep_20260920'
wr<-function(z,n){z<-copy(as.data.table(z));for(v in names(z))if(is.double(z[[v]]))z[,(v):=ifelse(is.na(get(v)),'NA',sprintf('%.17g',get(v)))];fwrite(z,file.path(out,paste0(n,'.csv')),bom=TRUE,na='NA')}
checks<-list(); add<-function(name,err,tol=1e-10,n=1L) checks[[length(checks)+1]]<<-data.table(check=name,n=n,max_abs_error=err,tolerance=tol,pass=is.finite(err)&&err<=tol)
compare<-function(name,a,b,tol=1e-10){stopifnot(length(a)==length(b));ok<-is.finite(a)&is.finite(b);add(name,if(any(ok))max(abs(a[ok]-b[ok]))else Inf,tol,sum(ok));add(paste0(name,' NA pattern'),sum(is.na(a)!=is.na(b)),0,length(a))}
calc<-function(q,U,B,T,L){e<-drop(L%*%q);u<-diag(L%*%U%*%t(L));b<-diag(L%*%B%*%t(L));se<-sqrt(diag(L%*%T%*%t(L)));r<-1.05*b/u;df<-ifelse(b>1e-15,19*(1+1/r)^2,Inf);data.table(log_effect=e,se=se,df=df,HR=exp(e),lower=exp(e-qt(.975,df)*se),upper=exp(e+qt(.975,df)*se),p=2*pt(-abs(e/se),df))}
jcalc<-function(q,T,L){z<-drop(L%*%q);V<-L%*%T%*%t(L);v<-drop(t(z)%*%solve(V,z));data.table(statistic=v,df=nrow(L),p=pchisq(v,nrow(L),lower.tail=FALSE))}
scalarcheck<-function(a,b,label,map=c(HR='HR',lower='lower',upper='upper',p='p')) for(n in names(map)) compare(paste(label,n),a[[n]],b[[map[[n]]]])
bh<-function(d,pcol,bcol,groups,label,subset=NULL){d<-copy(d);if(!is.null(subset))d<-d[subset];d<-d[is.finite(get(bcol))]; if(!nrow(d))return(invisible(NULL));d[,p_recalculated_BH:=p.adjust(get(pcol),'BH'),by=groups];compare(paste(label,'BH'),d$p_recalculated_BH,d[[bcol]]);d[,family_audit:=label];bhrows[[length(bhrows)+1]]<<-d[,c(groups,pcol,bcol,'p_recalculated_BH','family_audit'),with=FALSE]}
bhrows<-list();mainout<-list();mainj<-list();mainmat<-list();modelmeta<-list()
for(a in c('hrs_5state_full','elsa_4state_w4_w10','elsa_5state_through_2012')) {
  folder<-file.path(F,'models',if(grepl('hrs',a))'hrs_mi20_final'else'elsa_w10tri_final')
  paths<-file.path(folder,sprintf('%s__M4_MI20__imp%02d.rds',a,1:20));stopifnot(all(file.exists(paths)))
  bb<-list(); vv<-list()
  for(i in 1:20){obj<-readRDS(paths[i]);fit<-obj$fit;stopifnot(inherits(fit,'msm'),fit$opt$convergence==0,isTRUE(fit$foundse))
    im<-fit$qmodel$imatrix;ix<-which(t(im)!=0,arr.ind=TRUE);edge<-data.table(from=ix[,2]-1L,to=ix[,1]-1L);edge[,transition:=paste0('S',from,'->S',to)]
    k<-nrow(edge);cf<-coef(fit);ci<-match('insomnia_score_primary',attr(cf,'covlabels'));idx<-k+(ci-1L)*k+seq_len(k)
    bb[[i]]<-as.numeric(fit$paramdata$params[idx]);vv[[i]]<-fit$covmat[idx,idx,drop=FALSE]
    cf2<-vapply(seq_len(k),function(j)cf[['insomnia_score_primary']][edge$from[j]+1,edge$to[j]+1],1.)
    add(paste(a,i,'coefficient mapping'),max(abs(bb[[i]]-cf2)))
    if(i==1){modelmeta[[a]]<-list(edges=edge,constraints=obj$constraints,call=paste(deparse(fit$call),collapse=' '),n_free=length(fit$opt$par));wr(edge,paste0(a,'_allowed_edges'))}
    rm(obj,fit); if(i%%5==0)cat(a,'read',i,'/20\n')
  }
  Q<-do.call(cbind,bb);q<-rowMeans(Q);U<-Reduce('+',vv)/20;B<-cov(t(Q));T<-U+1.05*B
  tab<-cbind(data.table(analysis=a),edge,calc(q,U,B,T,diag(k)));target<-which(edge$to!=4L);j<-jcalc(q,T,diag(k)[target,,drop=FALSE]);j[,analysis:=a];mainj[[a]]<-j
  tab[,p_BH:=NA_real_];tab[target,p_BH:=p.adjust(p,'BH')];tab[,gate_pass:=j$p<.05]
  stored<-fread(file.path(F,'results/markov_hazard_ratios.csv'))[analysis_id==a & model_id=='M4_MI20' & term=='insomnia_score_primary'];stored<-stored[match(edge$transition,transition)]
  stopifnot(nrow(stored)==k);scalarcheck(tab,stored,a,c(HR='hazard_ratio',lower='conf_low',upper='conf_high',p='p_value'))
  # Mortality extension is not included in the primary confirmatory family.
  if(a!='elsa_5state_through_2012')compare(paste(a,'stored primary BH'),tab[target,p_BH],stored[target,p_fdr_bh_confirmatory])
  else tab[,p_BH:=NA_real_]
  mainout[[a]]<-tab;mainmat[[a]]<-list(qbar=q,U=U,B=B,T=T,edges=edge)
}
wr(rbindlist(mainout),'main_markov_fullprecision');wr(rbindlist(mainj),'main_markov_joint_fullprecision');saveRDS(mainmat,file.path(out,'main_independently_repooled_matrices.rds'))
write_json(modelmeta,file.path(out,'saved_model_structure.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
# Disease-specific target and direct contrasts from saved unrounded covariance matrices.
spl<-readRDS(file.path(N,'results/disease_specific/disease_specific_target_rubin_matrices.rds'))
splout<-list();spldirect<-list()
for(a in names(spl)){s<-spl[[a]];d<-cbind(data.table(cohort=a,transition=s$transition_order),calc(s$qbar,s$within,s$between,s$total,diag(4)));old<-fread(file.path(N,'results/disease_specific/disease_specific_target_edge_hr_mi20.csv'))[cohort==a];old<-old[match(s$transition_order,transition)];scalarcheck(d,old,paste('split',a),c(HR='hazard_ratio',lower='conf_low',upper='conf_high',p='p_value'));d[,p_BH:=NA_real_];splout[[a]]<-d
  d<-cbind(data.table(cohort=a,origin=c('No disease','HTN only')),calc(s$qbar,s$within,s$between,s$total,s$contrast_matrix));d[,p_BH:=p.adjust(p,'BH')];old<-fread(file.path(N,'results/disease_specific/disease_specific_heart_vs_stroke_contrasts_mi20.csv'))[cohort==a];scalarcheck(d,old,paste('split direct',a),c(HR='hr_ratio_heart_vs_stroke',lower='conf_low',upper='conf_high',p='p_value'));compare(paste('split',a,'BH'),d$p_BH,old$p_bh_two_contrasts_within_cohort);spldirect[[a]]<-d}
wr(rbindlist(splout),'disease_specific_fullprecision');wr(rbindlist(spldirect),'disease_specific_direct_fullprecision')
charout<-list();chardirect<-list();charjoint<-list()
sources<-c(CHARLS_static=file.path(C1,'results'),CHARLS_stage2_CORE=file.path(C2,'pooled/CORE/results'),CHARLS_stage2_BOTH=file.path(C2,'pooled/BOTH/results'))
for(a in names(sources)){root<-sources[[a]];s<-readRDS(file.path(root,'rubin_matrices_MI20.rds'));L<-diag(length(s$qbar))[s$targets,,drop=FALSE];d<-calc(s$qbar,s$U,s$B,s$T,L);old<-fread(file.path(root,'four_target_MI20.csv'));scalarcheck(d,old,a);d[,`:=`(analysis=a,transition=old$transition,p_BH=p.adjust(p,'BH'))];compare(paste(a,'BH target'),d$p_BH,old$p_BH_4);charout[[a]]<-d
  d<-calc(s$qbar,s$U,s$B,s$T,s$C);old<-fread(file.path(root,'heart_stroke_contrasts_MI20.csv'));scalarcheck(d,old,paste(a,'direct'));d[,`:=`(analysis=a,origin=old$origin,p_BH=p.adjust(p,'BH'))];compare(paste(a,'BH direct'),d$p_BH,old$p_BH_2);chardirect[[a]]<-d
  j<-rbind(cbind(data.table(analysis=a,test='four target coefficients zero'),jcalc(s$qbar,s$T,L)),cbind(data.table(analysis=a,test='two heart minus stroke differences zero'),jcalc(s$qbar,s$T,s$C)));old<-fread(file.path(root,'joint_tests_MI20.csv'));compare(paste(a,'joint P'),j$p,old$p);charjoint[[a]]<-j}
wr(rbindlist(charout),'charls_static_stage2_fullprecision');wr(rbindlist(chardirect),'charls_static_stage2_direct_fullprecision');wr(rbindlist(charjoint),'charls_static_stage2_joint_fullprecision')
dyn<-list();dyndir<-list();dynj<-list();dynom<-list()
for(a in c('CORE','BOTH')){root<-file.path(C3,'pooled',a,'results');s<-readRDS(file.path(root,'rubin_matrices_MI20.rds'))
  d<-cbind(as.data.table(s$target_contrast_meta),calc(s$qbar,s$U,s$B,s$T,s$target_contrast_matrix));d[,p_BH:=p.adjust(p,'BH'),by=role];d[,model:=a];old<-fread(file.path(root,'four_target_dynamic_contrasts_MI20.csv'));scalarcheck(d,old,paste('dynamic',a));compare(paste('dynamic',a,'BH target'),d$p_BH,old$p_BH_family);dyn[[a]]<-d
  d<-cbind(as.data.table(s$direct_contrast_meta),calc(s$qbar,s$U,s$B,s$T,s$direct_contrast_matrix));d[,p_BH:=p.adjust(p,'BH'),by=role];d[,model:=a];old<-fread(file.path(root,'heart_stroke_dynamic_direct_MI20.csv'));scalarcheck(d,old,paste('dynamic direct',a));compare(paste('dynamic',a,'BH direct'),d$p_BH,old$p_BH_family);dyndir[[a]]<-d
  target<-unique(s$target_contrast_meta$target_index);k<-nrow(s$edges);Lall<-diag(3*k)[unlist(lapply(0:2,function(i)target+i*k)),,drop=FALSE]
  Ls<-list(s$target_contrast_matrix[s$target_contrast_meta$role=='primary',,drop=FALSE],s$target_contrast_matrix[s$target_contrast_meta$role=='contextual',,drop=FALSE],Lall,s$direct_contrast_matrix[s$direct_contrast_meta$role=='primary',,drop=FALSE],s$direct_contrast_matrix[s$direct_contrast_meta$role=='contextual',,drop=FALSE]);old<-fread(file.path(root,'joint_tests_MI20.csv'));j<-rbindlist(lapply(Ls,function(L)jcalc(s$qbar,s$T,L)));compare(paste('dynamic',a,'joint P'),j$p,old$p);j[,`:=`(model=a,test=old$test)];dynj[[a]]<-j
  j<-rbindlist(lapply(target,function(t)jcalc(s$qbar,s$T,diag(3*k)[t+(0:2)*k,,drop=FALSE])));old<-fread(file.path(root,'four_target_pattern_omnibus_MI20.csv'));j[,`:=`(model=a,transition=old$transition,p_BH=p.adjust(p,'BH'))];compare(paste('dynamic',a,'omnibus P'),j$p,old$p);compare(paste('dynamic',a,'omnibus BH'),j$p_BH,old$p_BH_4);dynom[[a]]<-j}
wr(rbindlist(dyn),'dynamic_fullprecision');wr(rbindlist(dyndir),'dynamic_direct_fullprecision');wr(rbindlist(dynj),'dynamic_joint_fullprecision');wr(rbindlist(dynom),'dynamic_omnibus_fullprecision')
# Additional existing families; corrections are checked, not extended to new families.
ds<-file.path(N,'results/symptom_domains');d<-fread(file.path(ds,'symptom_domain_pairwise_contrasts.csv'));bh(d,'p_value','p_bh',c('analysis_id','family'),'symptom pairwise')
fr<-file.path(N,'results/formal_inference');d<-fread(file.path(fr,'onset_progression_contrasts.csv'));bh(d,'p_value','p_fdr_bh_path_specific','cohort','onset progression')
d<-fread(file.path(fr,'hrs_death_origin_contrasts.csv'));bh(d,'p_value','p_fdr_bh',character(),'HRS death origins')
d<-fread(file.path(fr,'cross_cohort_heterogeneity.csv'));bh(d,'p_value_difference','p_fdr_bh_difference',character(),'cross cohort differences')
d<-fread(file.path(F,'results/protocol_sensitivity_estimates.csv'));bh(d,'p_value','p_fdr_bh_sensitivity',c('cohort','scenario'),'protocol sensitivities')
wr(rbindlist(bhrows,fill=TRUE),'additional_BH_recalculation')
cox<-fread(file.path(F,'results/cox_transition_estimates.csv'));sx<-cox[term_role=='exposure'&(scenario%in%c('observed_sleep','landmark2','interval_le4','lag_available','previous_valid_available')|weighted==1)]
wr(sx,'executed_Cox_sensitivity_raw_outputs');wr(unique(sx[,.(cohort,analysis_id,scenario,basis,transition,weighted,m_used,run_status,status)]),'executed_Cox_sensitivity_inventory')
reg<-fread(file.path(F,'results/protocol_sensitivity_model_registry.csv'));wr(reg[,.N,by=.(cohort,scenario,status)],'protocol_sensitivity_completion')
completed<-reg[status=='completed',spec_id];files<-list.files('private/stage/multistate_output/models/protocol_sensitivities',pattern='\\.rds$',recursive=TRUE,full.names=TRUE);wr(data.table(spec_id=completed,saved_model_found=vapply(completed,function(s)any(grepl(s,basename(files),fixed=TRUE)),TRUE)),'protocol_saved_model_presence')
wr(rbindlist(checks),'inference_numeric_checks');print(rbindlist(checks)[pass==FALSE]);cat('Checks:',length(checks),'Failures:',sum(!rbindlist(checks)$pass),'\n')
