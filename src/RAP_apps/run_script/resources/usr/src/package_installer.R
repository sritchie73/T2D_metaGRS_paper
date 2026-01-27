#!/bin/Rscript

# N.b. needs to be run on the run_script workstation to make sure compiled 
# packages are installed correctly
#
# An interactive run_script workstation can be spun up by running 
# dx run run_script -icmd='while true; do sleep 10; done' -iscript='' --ssh -y
#
install.packages(c(
 "data.table", "R.utils", "foreach", "doMC", "ggplot2", "BiocManager", "remotes",
 "bit64", "Rcpp", "RcppArmadillo", "RcppEigen", "tidyverse", "devtools", "qrng",
 "MendelianRandomization","susieR", "epiR", "coloc", "incidence", "prevalence",
 "mvtnorm", "outbreaks", "odbc", "docopt", "patchwork", "cowplot", "openxlsx",
 "RNOmni", "pROC", "ggforce", "ggh4x", "ggpp", "ggrastr", "ggstance", "ggthemes", 
 "palettetown", "caret", "NetRep", "ukbnmr", "readstata13", "nricens", "hexbin",
 "microbenchmark", "flashClust", "bigstatsr", "bigsnpr", "Rmpfr", "lme4", "lemon",
 "txtplot", "munsell", "medflex", "olsrr", "expss"
))

BiocManager::install(c(
  "snpStats", "WGCNA"
))

remotes::install_github('erocoar/gghalves')

pkgdir <- sprintf("R/x86_64-pc-linux-gnu-library/%s.%s/", R.version$major, gsub("\\.[0-9]", "", R.version$minor))

system(sprintf("tar -czvf Rpackages.tar.gz -C %s .", pkgdir, pkgdir))

system("dx upload Rpackages.tar.gz --destination $DX_PROJECT_CONTEXT_ID --brief")
