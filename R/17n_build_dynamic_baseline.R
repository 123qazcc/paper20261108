#!/usr/bin/env Rscript
# Descriptive baseline only. Reads frozen stage-3 inputs; does not fit or impute.
options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(data.table); library(jsonlite); library(digest)})
setDTthreads(1L)
stage <- 'outputs/charls_dynamic_sleep_20260920'
out <- file.path(stage, 'baseline_20260930')
delivery <- 'private/model_store/CHARLS_dynamic_sleep_20260920/baseline_20260930'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
dir.create(delivery, recursive=TRUE, showWarnings=FALSE)
inputs <- c(file.path(stage,'data/charls_dynamic_model_bundle.rds'),
            file.path(stage,'data/charls_dynamic_preimputation.rds'),
            'outputs/charls_external_extension_20260913/data/charls_raw_W1_baseline_stage2.rds',
            'R/17c_fit_charls_dynamic_models.R', 'analysis/charls_stage3_dynamic_sleep_contract.md')
sha <- function(p) digest(file=p, algo='sha256')
hashes <- setNames(vapply(inputs,sha,''),inputs)
checks <- list()
check <- function(label, ok, detail='') {
  checks[[length(checks)+1L]] <<- data.table(check=label,pass=isTRUE(ok),detail=detail)
  if (!isTRUE(ok)) stop(label, ': ', detail)
}
num <- function(x) if (is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)
bun <- readRDS(inputs[1]); pre <- readRDS(inputs[2])
b <- as.data.table(copy(bun$baseline)); m <- pre$data; lm <- as.data.table(copy(bun$landmark))
check('frozen_cohort',nrow(b)==9397L && !anyDuplicated(b$id) && bun$status=='APPROVED_FOR_DYNAMIC_MODELS')
check('id_order_and_original_data',identical(as.character(b$id),as.character(pre$ids)) && identical(m,bun$original_mi_data))
lm <- lm[match(b$id,id)]
check('landmark_id_order',identical(as.character(b$id),as.character(lm$id)))
levels_pattern <- c('persistent_low','onset','remission','persistent_high')
labels_pattern <- c('Persistent low','Onset','Remission','Persistent high')
counts <- c(4702L,1369L,1483L,1843L)
check('four_frozen_pattern_counts',identical(as.integer(table(factor(b$pattern,levels=levels_pattern))),counts))
check('pattern_consistency',identical(as.character(m$pattern),as.character(b$pattern)) && identical(as.character(lm$pattern),as.character(b$pattern)))
states <- c('No target disease','HTN only','Heart disease only','HTN + heart disease',
            'Stroke only','HTN + stroke','Heart disease + stroke','HTN + heart disease + stroke')
z <- data.frame(pattern=factor(b$pattern,levels=levels_pattern,labels=labels_pattern),
                age=5*b$age5+65, female=b$female, education4=num(b$education4),
                married=b$married,rural=b$rural,wealth_asinh=b$wealth_asinh,
                bmi=5*b$bmi5+25,smoker=b$smoker,drinker=b$drinker,
                diabetes=b$diabetes,lung=b$lung,pain=b$pain,
                cesd9=5*b$cesd9_5,night_hours=num(m$night_hours),
                w1_state=num(m$htn0)+2*num(m$heart0)+4*num(m$stroke0),
                w2_state=num(lm$w2_state)-1L)
check('state_totals_W1',identical(as.integer(table(factor(z$w1_state,levels=0:7))),c(6167L,1835L,580L,598L,69L,95L,10L,43L)))
check('state_totals_W2',identical(as.integer(table(factor(z$w2_state,levels=0:7))),c(5466L,2011L,718L,852L,94L,151L,24L,81L)))
check('complete_states_and_exposure',!anyNA(z[,c('pattern','w1_state','w2_state')]))
check('CESD9_is_complete_item_sum',isTRUE(all.equal(z$cesd9,rowSums(m[,pre$raw_item_names]),check.attributes=FALSE)))
# Independent provenance cross-check against the audited raw W1 baseline, by ID.
raw <- as.data.table(readRDS(inputs[3])); raw <- raw[match(b$id,id)]
map <- c(age='age',female='female',education4='education4',married='married',rural='rural',
         wealth_asinh='wealth_asinh',bmi='bmi',smoker='smoker',drinker='drinker',diabetes='diabetes',
         lung='lung',pain='pain',cesd9='cesd9_no_sleep',night_hours='sleep_hours_observed',w1_state='bits')
check('raw_W1_ids',identical(as.character(raw$id),as.character(b$id)))
for (v in names(map)) check(paste0('raw_W1_match_',v),isTRUE(all.equal(num(z[[v]]),num(raw[[map[v]]]),check.attributes=FALSE)))
labels <- c(age='Age, years',female='Female',education4='Education',married='Married or cohabiting',
  rural='Rural residence',wealth_asinh='Non-housing financial wealth, asinh(value/10000)',
  bmi='Body mass index, kg/m2',smoker='Current smoking',drinker='Alcohol use in the past year',
  diabetes='Diabetes',lung='Chronic lung disease',pain='Frequent bodily pain',
  cesd9='CESD-9 score excluding sleep, 0-27',night_hours='Night-time sleep duration, hours',
  w1_state='W1 disease state',w2_state='W2 landmark disease state')
binary <- c('female','married','rural','smoker','drinker','diabetes','lung','pain')
continuous <- c('age','wealth_asinh','bmi','cesd9','night_hours')
edu <- c('Less than primary / traditional school','Primary school','Middle school','High school or above')
for (v in binary) z[[v]] <- factor(z[[v]],levels=0:1,labels=c('No','Yes'))
z$education4 <- factor(z$education4,levels=1:4,labels=edu)
for (v in c('w1_state','w2_state')) z[[v]] <- factor(z[[v]],levels=0:7,labels=states)
for (v in names(labels)) attr(z[[v]],'label') <- labels[[v]]
vars <- c('age','female','education4','married','rural','wealth_asinh','bmi','smoker','drinker','diabetes','lung','pain','cesd9','night_hours')
allvars <- c(vars,'w1_state','w2_state')
check('known_missing_counts',identical(as.integer(vapply(z[vars],function(x) sum(is.na(x)),0L)),
                                    c(0L,0L,3L,0L,0L,945L,1199L,2L,0L,55L,16L,4L,437L,82L)))

# Real package invocation, not a reproduction labelled as GlobalAgingR.
native <- list()
native_log <- capture.output({
  native$w1 <- GlobalAgingR::generate_summary_table(
    as.formula(paste('~',paste(c(vars,'w1_state'),collapse=' + '),'| pattern')),data=z)
  native$w2 <- GlobalAgingR::generate_summary_table(~ w2_state | pattern,data=z)
})
writeLines(native_log,file.path(out,'globalaging_invocation.log'),useBytes=TRUE)
writeLines(c('GlobalAgingR::generate_summary_table',
             deparse(getExportedValue('GlobalAgingR','generate_summary_table'))),
           file.path(out,'globalaging_function_source.txt'))
for (name in names(native)) {
  nt <- as.data.frame(native[[name]])
  fwrite(nt,file.path(out,paste0('GlobalAgingR_native_',name,'.csv')),bom=TRUE)
  writeLines(c('<!doctype html><html><head><meta charset="utf-8"></head><body>',
    '<p>Unmodified GlobalAgingR table. Categorical percentages use total group N, including missing values in the denominator. Continuous second rows show Median [Min, Max], not IQR.</p>',
    as.character(native[[name]]),'</body></html>'),file.path(out,paste0('GlobalAgingR_native_',name,'.html')),useBytes=TRUE)
}

# Structured exact summaries support compact publication tables and independent checks.
groups <- c('Overall',labels_pattern)
summary <- list()
for (g in groups) {
  d <- if (g=='Overall') z else z[z$pattern==g,,drop=FALSE]
  for (v in allvars) {
    x <- d[[v]]; valid <- x[!is.na(x)]; n <- nrow(d); nvalid <- length(valid)
    if (v %in% continuous) {
      qs <- unname(quantile(valid,c(.25,.5,.75),type=7))
      summary[[length(summary)+1L]] <- data.table(variable=v,label=labels[[v]],wave='W1',group=g,
        kind='continuous',level='',n_group=n,n_available=nvalid,n_missing=n-nvalid,
        count=NA_integer_,percent_available=NA_real_,mean=mean(valid),sd=sd(valid),
        q1=qs[1],median=qs[2],q3=qs[3],min=min(valid),max=max(valid))
    } else for (lev in levels(x)) {
      count <- sum(valid==lev)
      summary[[length(summary)+1L]] <- data.table(variable=v,label=labels[[v]],wave=if(v=='w2_state') 'W2' else 'W1',
        group=g,kind='categorical',level=lev,n_group=n,n_available=nvalid,n_missing=n-nvalid,
        count=count,percent_available=100*count/nvalid,mean=NA_real_,sd=NA_real_,q1=NA_real_,
        median=NA_real_,q3=NA_real_,min=NA_real_,max=NA_real_)
    }
  }
}
s <- rbindlist(summary)
check('categorical_valid_totals',all(s[kind=='categorical',.(ok=sum(count)==unique(n_available)),by=.(variable,group)]$ok))
check('available_plus_missing_equals_group',all(s$n_available+s$n_missing==s$n_group))
check('overall_equals_group_counts',all(s[group!='Overall',.(total=sum(n_group)),by=.(variable,level)]$total==9397))
fwrite(s,file.path(out,'baseline_summary_numeric_long.csv'),bom=TRUE)
saveRDS(z,file.path(out,'observed_values_for_independent_audit.rds'))

getrow <- function(v,g,lev='') s[variable==v & group==g & level==lev]
disp <- function(v,g,lev='') {
  r <- getrow(v,g,lev); stopifnot(nrow(r)==1L)
  if (v %in% c('age','bmi')) sprintf('%.2f (%.2f)',r$mean,r$sd)
  else if (v %in% continuous) sprintf('%.2f [%.2f, %.2f]',r$median,r$q1,r$q3)
  else sprintf('%s (%.1f)',format(r$count,big.mark=',',trim=TRUE),r$percent_available)
}
make_row <- function(v,label,lev='',section='W1 characteristics') {
  vals <- setNames(vapply(groups,function(g) disp(v,g,lev),''),groups)
  c(list(variable=v,label=label,section=section,level=lev),as.list(vals))
}
main <- list()
for (v in vars) {
  suffix <- if (v %in% c('age','bmi')) ', mean (SD)' else if (v %in% continuous) ', median [Q1, Q3]' else ', n (%)'
  if (v=='education4') for(lev in edu) main[[length(main)+1L]] <- make_row(v,paste('Education:',lev),lev)
  else main[[length(main)+1L]] <- make_row(v,paste0(labels[[v]],suffix),if(v %in% binary)'Yes' else '')
}
main <- rbindlist(main)
disease <- rbindlist(lapply(c('w1_state','w2_state'),function(v) rbindlist(lapply(states,function(lev)
  make_row(v,lev,lev,if(v=='w1_state')'W1 disease state' else 'W2 landmark disease state')))))
missing <- unique(s[,.(variable,label,wave,group,n_group,n_available,n_missing)])
missing[,missing_percent:=100*n_missing/n_group]
missing_display <- dcast(missing,variable+label+wave~group,value.var='n_missing')
for (g in groups) {
  dat <- missing[group==g]
  missing_display[, (g):=vapply(variable,function(v) {
    r<-dat[variable==v];sprintf('%s (%.2f)',format(r$n_missing,big.mark=',',trim=TRUE),r$missing_percent)
  },'')]
}
missing_display <- missing_display[match(allvars,variable)]

# Descriptive imbalance only; no baseline P values and no covariate selection.
tt <- tableone::CreateTableOne(vars=allvars,strata='pattern',data=z,
  factorVars=setdiff(allvars,continuous),includeNA=FALSE,test=FALSE,smd=TRUE,addOverall=TRUE)
sm <- tableone::ExtractSmd(tt)
pairs <- combn(labels_pattern,2L)
paircols <- setdiff(colnames(sm),'average')
check('six_SMD_pairs',length(paircols)==6L)
sm <- sm[match(allvars,rownames(sm)),paircols,drop=FALSE]
check('finite_SMD',all(is.finite(sm)))
smd <- as.data.table(sm);setnames(smd,paste(pairs[1,],pairs[2,],sep=' vs '))
smd[,`:=`(variable=allvars,label=unname(labels[allvars]))]
smd[,max_pairwise_SMD:=apply(sm,1L,max)]
setcolorder(smd,c('variable','label',paste(pairs[1,],pairs[2,],sep=' vs '),'max_pairwise_SMD'))
outputs <- list(Table1_W1_characteristics=main,TableS1_W1_W2_disease_states=disease,
                TableS2_missingness=missing_display,TableS3_pairwise_SMD=smd,
                missingness_numeric=missing,baseline_summary_numeric_long=s)
for (name in names(outputs)) fwrite(outputs[[name]],file.path(out,paste0(name,'.csv')),bom=TRUE)
payload <- list(n=9397L,groups=groups,counts=c(9397L,counts),
  pattern_keys=levels_pattern,main=main,disease=disease,missing=missing_display,
  smd=smd,numeric=s,missing_numeric=missing,
  metadata=list(software=R.version.string,GlobalAgingR=as.character(packageVersion('GlobalAgingR')),
    table1=as.character(packageVersion('table1')),tableone=as.character(packageVersion('tableone')),
    summary='Observed, unweighted, one participant per row; no imputed values used',
    continuous='Age and BMI: mean (SD). Wealth, CESD-9 and sleep duration: median [Q1, Q3], R type-7 quantiles',
    categorical='n (%) of non-missing values within each column. Missing percentages use total column N',
    native_difference='GlobalAgingR native output uses total column N for categorical percentages and Median [Min, Max]; compact outputs are explicitly reformatted from the same observed data',
    smd='Absolute standardized differences for six group pairs. Continuous SMD uses mean/SD even when the table displays median/IQR. Multilevel categorical SMD is generalized; not level-specific. No P values',
    sampling='Selected W2 landmark cohort; not nationally representative estimates',
    source_files=basename(inputs[1:3])))
write_json(payload,file.path(out,'baseline_payload.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA,na='null')
check('frozen_sources_unchanged',identical(unname(vapply(inputs,sha,'')),unname(hashes)))
fwrite(rbindlist(checks),file.path(out,'baseline_build_checks.csv'),bom=TRUE)
write_json(list(status='BASELINE_DESCRIPTIVE_BUILD_PASS',n=9397L,patterns=counts,checks_passed=length(checks),
  input_hashes=as.list(hashes),model_refits=0L,imputations_run=0L,new_exclusions=0L,
  GlobalAgingR_invoked=TRUE,output_directory=out,delivery_directory=delivery),
  file.path(out,'baseline_build_review.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('BASELINE_DESCRIPTIVE_BUILD_PASS; n=9397; groups=4702/1369/1483/1843; no model refits\n')
