#!/usr/bin/env Rscript
# Independent baseline audit: base-R statistics from frozen source RDS.
# Does not call the builder, GlobalAgingR, tableone, model fitting, or imputation.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(jsonlite); library(digest)})
out <- 'outputs/charls_dynamic_sleep_20260920/baseline_20260930'
inputs <- c('outputs/charls_dynamic_sleep_20260920/data/charls_dynamic_model_bundle.rds',
            'outputs/charls_dynamic_sleep_20260920/data/charls_dynamic_preimputation.rds',
            'R/17c_fit_charls_dynamic_models.R')
hash <- function(p) digest(file=p,algo='sha256')
before <- setNames(vapply(inputs,hash,''),inputs)
review <- list()
check <- function(name,ok,detail='') {
  review[[length(review)+1L]] <<- data.frame(check=name,pass=isTRUE(ok),detail=as.character(detail))
  if(!isTRUE(ok)) stop(name,': ',detail)
}
near <- function(a,b,tol=1e-10) {
  a <- unname(a); b <- unname(b)
  if(length(a)!=length(b)||!identical(is.na(a),is.na(b))) return(FALSE)
  k <- !is.na(a)
  all(abs(a[k]-b[k])<=tol*pmax(1,abs(a[k]),abs(b[k])))
}
to_num <- function(x) as.numeric(as.character(x))
readcsv <- function(name) read.csv(file.path(out,paste0(name,'.csv')),check.names=FALSE,fileEncoding='UTF-8-BOM')
writecsv <- function(x,name) write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE,na='',fileEncoding='UTF-8')
run <- function() {
  bun <- readRDS(inputs[1]); pre <- readRDS(inputs[2])
  base <- as.data.frame(bun$baseline); entry <- as.data.frame(bun$landmark); predata <- pre$data
  check('unique_9397_people',nrow(base)==9397L && !anyDuplicated(base$id))
  check('source_ID_alignment',identical(base$id,pre$ids) && identical(base$id,entry$id))
  lev <- c('Persistent low','Onset','Remission','Persistent high')
  expected_key <- c('persistent_low','onset','remission','persistent_high')[1+2*as.integer(entry$w1_score>=2)+as.integer(entry$w2_score>=2)]
  check('independently_rederived_patterns',identical(expected_key,as.character(base$pattern)))
  group <- lev[1+2*as.integer(entry$w1_score>=2)+as.integer(entry$w2_score>=2)]
  check('group_N',identical(as.integer(table(factor(group,levels=lev))),c(4702L,1369L,1483L,1843L)))
  z <- data.frame(age=65+5*base$age5,female=base$female,education4=to_num(base$education4),
    married=base$married,rural=base$rural,wealth_asinh=base$wealth_asinh,bmi=25+5*base$bmi5,
    smoker=base$smoker,drinker=base$drinker,diabetes=base$diabetes,lung=base$lung,pain=base$pain,
    cesd9=rowSums(predata[,pre$raw_item_names,drop=FALSE]),night_hours=predata$night_hours,
    w1_state=predata$htn0+2*predata$heart0+4*predata$stroke0,w2_state=entry$w2_state-1L)
  check('CESD_derived_from_all9_observed_items',near(z$cesd9,5*base$cesd9_5))
  continuous <- c('age','wealth_asinh','bmi','cesd9','night_hours')
  bin <- c('female','married','rural','smoker','drinker','diabetes','lung','pain')
  state_labels <- c('No target disease','HTN only','Heart disease only','HTN + heart disease','Stroke only','HTN + stroke','Heart disease + stroke','HTN + heart disease + stroke')
  names_by_code <- list(education4=c('Less than primary / traditional school','Primary school','Middle school','High school or above'),w1_state=state_labels,w2_state=state_labels)
  for(v in bin) names_by_code[[v]] <- c('No','Yes')
  value_codes <- lapply(names_by_code,function(a)seq_along(a)-1L)
  value_codes$education4 <- 1:4
  sample_var <- function(x) {u<-x[!is.na(x)]; sum((u-sum(u)/length(u))^2)/(length(u)-1)}
  quant <- function(x,p) {x<-sort(x[!is.na(x)]); h<-(length(x)-1)*p+1; j<-floor(h); f<-h-j; if(j>=length(x))tail(x,1) else x[j]+f*(x[j+1]-x[j])}
  s <- readcsv('baseline_summary_numeric_long')
  check('numeric_keys_unique',!anyDuplicated(s[,c('variable','group','level')]))
  check('expected_numeric_rows',nrow(s)==205L)
  numeric_checks <- list()
  for(i in seq_len(nrow(s))) {
    r <- s[i,]; ids <- if(r$group=='Overall')seq_len(nrow(z)) else which(group==r$group)
    x <- z[[r$variable]][ids]; valid <- x[!is.na(x)]
    check(paste0('denominator_',i),r$n_group==length(ids)&&r$n_available==length(valid)&&r$n_missing==sum(is.na(x)))
    check(paste0('wave_',i),r$wave==if(r$variable=='w2_state')'W2' else 'W1')
    if(r$variable %in% continuous) {
      ex <- c(mean=sum(valid)/length(valid),sd=sqrt(sample_var(valid)),q1=quant(valid,.25),median=quant(valid,.5),q3=quant(valid,.75),min=min(valid),max=max(valid))
      obs <- as.numeric(r[,names(ex)])
      check(paste0('continuous_exact_',i),near(obs,ex))
      numeric_checks[[length(numeric_checks)+1L]] <- data.frame(variable=r$variable,group=r$group,level='',field=names(ex),expected=ex,observed=obs,abs_error=abs(obs-ex))
    } else {
      index <- match(r$level,names_by_code[[r$variable]])
      check(paste0('category_label_',i),!is.na(index))
      n <- length(which(valid==value_codes[[r$variable]][index]))
      ex <- c(count=n,percent_available=100*n/length(valid)); obs <- as.numeric(r[,names(ex)])
      check(paste0('categorical_exact_',i),near(obs,ex))
      numeric_checks[[length(numeric_checks)+1L]] <- data.frame(variable=r$variable,group=r$group,level=r$level,field=names(ex),expected=ex,observed=obs,abs_error=abs(obs-ex))
    }
  }
  for(file in c('Table1_W1_characteristics','TableS1_W1_W2_disease_states')) {
    d <- readcsv(file)
    for(i in seq_len(nrow(d)))for(g in c('Overall',lev)) {
      r <- s[s$variable==d$variable[i]&s$group==g&s$level==d$level[i],]
      check(paste(file,'display_source',i,g),nrow(r)==1L)
      expected <- if(d$variable[i] %in% c('age','bmi'))sprintf('%.2f (%.2f)',r$mean,r$sd) else if(d$variable[i] %in% continuous)sprintf('%.2f [%.2f, %.2f]',r$median,r$q1,r$q3) else sprintf('%s (%.1f)',format(r$count,big.mark=',',trim=TRUE),r$percent_available)
      check(paste(file,'display',i,g),identical(d[[g]][i],expected))
    }
  }
  miss <- readcsv('missingness_numeric'); display <- readcsv('TableS2_missingness')
  for(i in seq_len(nrow(miss))) {
    r<-miss[i,]; ids<-if(r$group=='Overall')seq_len(nrow(z)) else which(group==r$group)
    n<-sum(is.na(z[[r$variable]][ids])); pct<-100*n/length(ids)
    check(paste0('missing_exact_',i),r$n_missing==n&&r$n_available==length(ids)-n&&r$n_group==length(ids)&&near(r$missing_percent,pct))
    expected<-sprintf('%s (%.2f)',format(n,big.mark=',',trim=TRUE),pct)
    check(paste0('missing_display_',i),identical(display[[r$group]][match(r$variable,display$variable)],expected))
  }
  # Independent Cohen-type absolute SMD; categorical covariance is that of a
  # multinomial indicator vector, averaged between groups without sample-N weights.
  smd <- readcsv('TableS3_pairwise_SMD'); pair <- combn(lev,2L); pair_names<-apply(pair,2,paste,collapse=' vs ')
  check('SMD_pair_names',identical(names(smd)[3:8],pair_names))
  smd_check <- list()
  for(v in names(z))for(j in 1:6) {
    x<-z[[v]][group==pair[1,j]]; y<-z[[v]][group==pair[2,j]]
    x<-x[!is.na(x)];y<-y[!is.na(y)]
    if(v %in% continuous) estimate<-abs(mean(x)-mean(y))/sqrt((sample_var(x)+sample_var(y))/2)
    else {
      codes<-value_codes[[v]]
      px<-vapply(codes,function(k)mean(x==k),0.0);py<-vapply(codes,function(k)mean(y==k),0.0)
      px<-px[-1];py<-py[-1]
      cx<-diag(px,nrow=length(px))-outer(px,px);cy<-diag(py,nrow=length(py))-outer(py,py)
      delta<-px-py; pooled<-(cx+cy)/2
      estimate<-sqrt(drop(t(delta)%*%solve(pooled,delta)))
    }
    obs<-smd[[pair_names[j]]][match(v,smd$variable)]
    check(paste('SMD_exact',v,j),near(obs,estimate))
    smd_check[[length(smd_check)+1L]]<-data.frame(variable=v,pair=pair_names[j],independent_SMD=estimate,tableone_SMD=obs,abs_error=abs(obs-estimate))
  }
  for(i in seq_len(nrow(smd)))check(paste0('max_SMD_',i),near(smd$max_pairwise_SMD[i],max(as.numeric(smd[i,pair_names]))))
  payload<-fromJSON(file.path(out,'baseline_payload.json'))
  check('payload_N',payload$n==9397L&&identical(payload$counts,c(9397L,4702L,1369L,1483L,1843L)))
  check('native_difference_disclosed',grepl('total column N',payload$metadata$native_difference,fixed=TRUE)&&grepl('Median [Min, Max]',payload$metadata$native_difference,fixed=TRUE))
  check('SMD_mean_SD_disclosed',grepl('mean/SD',payload$metadata$smd,fixed=TRUE)&&grepl('not level-specific',payload$metadata$smd,fixed=TRUE))
  check('observed_unweighted_disclosed',grepl('Observed, unweighted',payload$metadata$summary,fixed=TRUE))
  # Test native categorical percentages against total N, distinguishing them
  # from the compact nonmissing denominator. W1 native row 1 is the N row.
  native<-readcsv('GlobalAgingR_native_w1')
  labels<-trimws(native[[1]])
  for(v in c(bin,'education4','w1_state')) {
    label <- unique(s$label[s$variable==v]); start<-which(labels==label)
    check(paste0('native_heading_',v),length(start)==1L)
    for(k in seq_along(names_by_code[[v]]))for(g in lev) {
      txt <- native[[g]][start+k]; r<-s[s$variable==v&s$group==g&s$level==names_by_code[[v]][k],]
      number<-as.numeric(gsub(',','',sub(' .*','',txt)))
      percentage<-as.numeric(sub('.*\\(([0-9.]+)%\\).*','\\1',txt))
      check(paste('native_count_pct',v,k,g),number==r$count && abs(percentage-100*r$count/r$n_group)<=.05000001)
    }
  }
  check('frozen_sources_unchanged',identical(unname(vapply(inputs,hash,'')),unname(before)))
  writecsv(do.call(rbind,numeric_checks),'independent_numeric_review')
  writecsv(do.call(rbind,smd_check),'independent_SMD_review')
  writecsv(do.call(rbind,review),'independent_baseline_review_checks')
  top <- do.call(rbind,lapply(c('Persistent low vs Onset','Remission vs Persistent high'),function(p){a<-smd_check[vapply(smd_check,function(d)d$pair==p,TRUE)];a<-do.call(rbind,a);a[order(-a$independent_SMD),][1:5,]}))
  writecsv(top,'primary_contrast_top5_SMD')
  write_json(list(status='INDEPENDENT_BASELINE_REVIEW_PASS',checks_total=length(review),checks_passed=sum(vapply(review,function(x)x$pass,TRUE)),
    n=9397L,numeric_rows=nrow(s),numeric_values=nrow(do.call(rbind,numeric_checks)),pairwise_SMD_values=length(smd_check),
    independent_methods='Base-R source reconstruction; sample variance; explicitly computed type-7 quantiles; categorical counts; multinomial covariance SMD. No use of builder functions, GlobalAgingR or tableone for recomputation.',
    source_hashes=as.list(before),model_fits=0L,imputations_run=0L,new_exclusions=0L,reviewer_script='R/17o_verify_dynamic_baseline.R',
    reviewed_utc=format(Sys.time(),'%Y-%m-%dT%H:%M:%SZ',tz='UTC')),file.path(out,'independent_baseline_review.json'),pretty=TRUE,auto_unbox=TRUE,digits=NA)
  print(top);cat('PASS checks=',length(review),' numeric_rows=',nrow(s),' SMD=',length(smd_check),'\n')
}
tryCatch(run(),error=function(e){
  if(length(review))writecsv(do.call(rbind,review),'independent_baseline_review_checks')
  write_json(list(status='INDEPENDENT_BASELINE_REVIEW_FAILED',error=conditionMessage(e),checks_total=length(review),checks_passed=sum(vapply(review,function(x)x$pass,TRUE))),file.path(out,'independent_baseline_review.json'),pretty=TRUE,auto_unbox=TRUE)
  stop(e)
})
