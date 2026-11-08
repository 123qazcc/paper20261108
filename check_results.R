# Run from a fresh unpacked repository: Rscript --vanilla check_results.R
# Offline review of saved aggregate statistics; no participant data or packages.
options(stringsAsFactors=FALSE)
if (!dir.exists("reference")) stop("Run this command from the unpacked repository root.")
tab <- function(name, n) {
  x <- read.csv(file.path("reference", paste0(name, ".csv")), check.names=FALSE,fileEncoding="UTF-8-BOM")
  stopifnot(nrow(x)==n)
  x
}
same <- function(x,y,tol=1e-10) {
  stopifnot(length(x)>0, length(x)==length(y), all(is.finite(x)), all(is.finite(y)),
            all(abs(x-y) <= tol*pmax(1,abs(y))))
}
scalar <- function(x) {
  stopifnot(all(x$se>0), all(x$df>0))
  delta <- qt(.975, x$df)*x$se
  same(exp(x$log_effect), x$HR)
  same(exp(x$log_effect-delta), x$lower)
  same(exp(x$log_effect+delta), x$upper)
  same(2*pt(-abs(x$log_effect/x$se), x$df), x$p)
}
bh <- function(x,g,n) {
  stopifnot(length(g)==nrow(x),!anyNA(g))
  for (i in split(seq_len(nrow(x)),g,drop=TRUE)) {
    stopifnot(length(i)==n)
    same(p.adjust(x$p[i], "BH"),x$p_BH[i])
  }
}
load_function <- function(file,name) {
  ast <- parse(file,encoding="UTF-8")
  sel <- Filter(function(e) is.call(e) && identical(e[[1]],as.name("<-")) &&
                        identical(e[[2]],as.name(name)) && is.call(e[[3]]) &&
                        identical(e[[3]][[1]],as.name("function")),as.list(ast))
  stopifnot(length(sel)==1)
  env <- new.env(parent=globalenv()); eval(sel[[1]],env); env[[name]]
}
state <- load_function("R/01_clean_three_cohorts.R","derive_cvd_state")
truth <- tab("HRS_ELSA_all_27_state_classifications",27)
stopifnot(identical(as.integer(truth$state), state(truth$h,truth$heart,truth$stroke)),
          all(state(truth$h,truth$heart,truth$stroke,rep(1L,27))==4L))
pool <- load_function("R/17d_independently_verify_charls_dynamic_pool.R","pool")
toy <- pool(matrix(c(1,0,2,2,3,4),2),rep(list(diag(2)),3))
same(toy$qbar,c(2,2)); same(as.numeric(toy$T), as.numeric(diag(2)+(4/3)*matrix(c(1,2,2,4),2)))
main <- tab("main_markov_fullprecision",20); scalar(main)
m <- subset(main,to!=4); stopifnot(nrow(m)==12,all(is.na(main$p_BH[main$to==4])))
main_family <- subset(m,analysis!="elsa_5state_through_2012")
bh(main_family,main_family$analysis,4)
stopifnot(all(is.na(m$p_BH[m$analysis=="elsa_5state_through_2012"])))
ds <- tab("disease_specific_fullprecision",8); scalar(ds); stopifnot(all(is.na(ds$p_BH)))
dd <- tab("disease_specific_direct_fullprecision",4); scalar(dd); bh(dd,dd$cohort,2)
s <- tab("charls_static_stage2_fullprecision",12); scalar(s); bh(s,s$analysis,4)
sdirect <- tab("charls_static_stage2_direct_fullprecision",6); scalar(sdirect); bh(sdirect,sdirect$analysis,2)
d <- tab("dynamic_fullprecision",32); scalar(d); bh(d,interaction(d$model,d$role),8)
direct <- tab("dynamic_direct_fullprecision",16); scalar(direct); bh(direct,interaction(direct$model,direct$role),4)
for (spec in c("CORE","BOTH")) {
  primary <- subset(d,model==spec & role=="primary" & p_BH < .05)
  stopifnot(setequal(paste(primary$comparison,primary$transition),
    c("onset_vs_persistent_low 0-4","persistent_high_vs_remission 0-2","persistent_high_vs_remission 1-3")),
    !any(direct$p_BH[direct$model==spec & direct$role=="primary"] < .05))
}
for (nm in c(main_markov_joint_fullprecision=3,charls_static_stage2_joint_fullprecision=6,
             dynamic_joint_fullprecision=10,dynamic_omnibus_fullprecision=8) |> names()) {
  n <- switch(nm,main_markov_joint_fullprecision=3,charls_static_stage2_joint_fullprecision=6,
              dynamic_joint_fullprecision=10,dynamic_omnibus_fullprecision=8)
  x <- tab(nm,n); same(pchisq(x$statistic,x$df,lower.tail=FALSE),x$p)
  if (nm=="dynamic_omnibus_fullprecision") bh(x,x$model,4)
}
fg <- function(file,n) {
  x <- read.csv(file.path("figure_data",paste0(file,".csv")),check.names=FALSE,fileEncoding="UTF-8-BOM")
  stopifnot(nrow(x)==n); x
}
match_values <- function(a,ka,b,kb,cols_a,cols_b) {
  key <- function(x,k) do.call(paste,c(x[k],sep="|"))
  ka <- key(a,ka); kb <- key(b,kb)
  stopifnot(!anyDuplicated(ka),!anyDuplicated(kb),setequal(ka,kb))
  j <- match(ka,kb)
  for (i in seq_along(cols_a)) same(a[[cols_a[i]]],b[[cols_b[i]]][j])
}
f2 <- fg("Figure2_Markov_complete_HR",20)
f2nf <- subset(f2,to_state!=4)
f2nf$p <- f2nf$p_value; f2nf$p_BH <- f2nf$p_fdr_bh_confirmatory
bh(f2nf,f2nf$analysis_id,4)
match_values(main,c("analysis","transition"),f2,c("analysis_id","transition"),
             c("HR","lower","upper","p"),c("hazard_ratio","conf_low","conf_high","p_value"))
for (i in 7:8) {
  direct_test <- i==8
  f <- if (direct_test) fg("Figure8_dynamic_direct",16) else fg("Figure7_dynamic_patterns",32)
  by <- c("model","comparison",if(direct_test) "origin" else "transition")
  match_values(if(direct_test) direct else d,by,f,by,c("HR","lower","upper","p_BH"),c("HR","lower","upper","p_BH_family"))
}
cat("PASS: three-cohort HRs, pointwise CIs, P values, family-specific BH, joint-test P values, and key figure tables.\n")
cat("Source inspection and aggregate checks only; no individual-level model was refitted.\n")
