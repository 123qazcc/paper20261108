#!/usr/bin/env Rscript
# ASCII bootstrap: explicit UTF-8 for English publication figures and symbols.
if (.Platform$OS.type == "windows") {
  stopifnot(nzchar(Sys.setlocale("LC_CTYPE", "English_United States.utf8")))
}
source("R/17i_charls_dynamic_figures_EN.R", encoding = "UTF-8")
