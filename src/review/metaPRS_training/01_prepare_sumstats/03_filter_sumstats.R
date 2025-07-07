library(data.table)
source("src/functions/filter_sumstats.R")

system("mkdir -p data/filtered_sumstats/filtered_gwas", wait=TRUE)

# For each GWAS, filter to set of HapMap variants above, orienting to common effect allele,
# and output in LDpred2 format.
gwas_list <- fread("data/gwas_summary_stats/curated_gwas_table.txt", na.strings="-")
gwas_list[, Samples := as.integer(gsub(",", "", Samples))]
gwas_list[, Cases := as.integer(gsub(",", "", Cases))]
gwas_list[, Controls := as.integer(gsub(",", "", Controls))]

######################################################### 
# Filter the 14 case-control FinnGenn summary statistics
#########################################################
for (this_prs in gwas_list[PRS %like% "FinnGen" & !is.na(Cases), PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".*\n", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/FinnGen/finngen_R12_%s.gz", this_pheno_fname))
  gwas_ss <- gwas_ss[,.(chr=`#chrom`, pos_b38=pos, EA=alt, OA=ref, beta, beta_se=sebeta, neg_log10_p=mlogp, EAF=af_alt)]
  gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_samples=this_gwas$Samples, total_cases=this_gwas$Cases, total_controls=this_gwas$Controls)
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

##############################################################
# Filter the 1 quantitative trait FinnGenn summary statistics
##############################################################
for (this_prs in gwas_list[PRS %like% "FinnGen" & is.na(Cases), PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".*\n", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/FinnGen/finngen_R12_%s.gz", this_pheno_fname))
  gwas_ss <- gwas_ss[,.(chr=`#chrom`, pos_b38=pos, EA=alt, OA=ref, beta, beta_se=sebeta, neg_log10_p=mlogp, EAF=af_alt)]
  gwas_ss <- filter_sumstats(gwas_ss, type="continuous", total_samples=this_gwas$Samples)
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

##############################################################
# Filter the 8 pQTL CKB summary statistics
##############################################################
for (this_prs in gwas_list[`GWAS catalog accession or other download source` %like% 'pheweb.ckbiobank.org', PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".* ", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/CKB_pQTLs/phenocode-%s.tsv.gz", this_pheno_fname))
  gwas_ss <- gwas_ss[,.(chr=chrom, pos_b38=pos, EA=alt, OA=ref, beta, beta_se=sebeta, neg_log10_p=-log10(pval), EAF=af)]
  gwas_ss <- filter_sumstats(gwas_ss, type="continuous", total_samples=this_gwas$Samples)
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
} 


