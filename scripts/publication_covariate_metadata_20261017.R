Sys.setlocale('LC_CTYPE','English_United States.utf8')
suppressPackageStartupMessages(library(data.table))
out <- 'outputs/publication_methods_audit_20261017'
meta <- list(); counts <- list()
for(cohort in c('hrs','elsa')) {
  e <- new.env(); data(list=cohort,package='GlobalAgingR',envir=e)
  d <- get(cohort,envir=e)
  for(v in intersect(c('raeducl','mstath','vgactx','mdactx','vgactx_e','mdactx_e','smoken','drink','diabe','lunge','painfr'),names(d))) {
    labs <- attr(d[[v]],'labels')
    meta[[paste(cohort,v)]] <- data.table(cohort=cohort,variable=v,variable_label=attr(d[[v]],'label'),value=as.numeric(labs),value_label=names(labs))
  }
  p <- as.data.table(readRDS(file.path('outputs/hrs_elsa_multistate_final_20260909/data_multistate',paste0(cohort,'_msm_panel.rds'))))
  b <- p[analysis_baseline==1]
  counts[[cohort]] <- b[,.(n=.N),by=.(marital_status,married_bl)]
  counts[[cohort]][,cohort:=cohort]
  rm(d,e,p,b);gc()
}
fwrite(rbindlist(meta),file.path(out,'package_covariate_value_labels.csv'),bom=TRUE)
fwrite(rbindlist(counts),file.path(out,'baseline_marital_coding_counts.csv'),bom=TRUE)
print(rbindlist(counts))
