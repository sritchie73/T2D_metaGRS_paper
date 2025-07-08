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

##############################################################
# Filter the 12 GIANT consortium summary statistics
##############################################################
for (this_prs in gwas_list[`GWAS catalog accession or other download source` %like% "GIANT", PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".*\n", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/GIANT/%s", this_pheno_fname), skip=8, header=FALSE) # Header row corrupted

  # Sumstats files are a mix of those with just rsID, and those that also have chromosome and position on GRCh36
  if (ncol(gwas_ss) == 8) { 
    setnames(gwas_ss, c("MarkerName", "Allele1", "Allele2", "FreqAllele1HapMapCEU", "b", "se", "p", "N"))
		gwas_ss <- gwas_ss[,.(rsid=MarkerName, EA=Allele1, OA=Allele2, beta=b, beta_se=se, samples=N, EAF=FreqAllele1HapMapCEU, neg_log10_p=-log10(p))]
  } else if (ncol(gwas_ss) == 10) {
    setnames(gwas_ss, c("MarkerName", "Chr", "Pos", "Allele1", "Allele2", "FreqAllele1HapMapCEU", "b", "se", "p", "N"))
		gwas_ss <- gwas_ss[,.(chr=Chr, pos_b36=Pos, EA=Allele1, OA=Allele2, beta=b, beta_se=se, samples=N, EAF=FreqAllele1HapMapCEU, neg_log10_p=-log10(p))]
  } else {
    stop("unrecognised number of columns in GIANT summary statistics")
  }

  gwas_ss <- filter_sumstats(gwas_ss, type="continuous", total_samples=this_gwas$Samples)
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
} 

#####################################################################################
# Filter GWAS summary statistics for CAD published by Koyama et al. 2020 on Figshare
#####################################################################################

gwas_ss <- fread("data/gwas_summary_stats/CAD_BBJ/BBJCAD_2020.sumstats.gz")
gwas_ss <- gwas_ss[, .(chr=CHR, pos_b37=POS, EA=ALT, OA=REF, EAF=AAF, beta=BETA, beta_se=SE, neg_log10_p=-log10(P), samples=N)]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=25892, total_controls=142336)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/filtered_gwas/CAD_BBJ.txt.gz")

##############################################################
# Filter the 9 G&H case-control ExWAS
##############################################################
for (this_prs in gwas_list[PRS %like% "GH_Exome" & !is.na(Cases), PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".*\n", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/GenesAndHealth/ExWAS/2024_02_05_%s_GNH_singlevariantExWAS_%s.regenie.gz", this_pheno_fname, this_pheno_fname))
  gwas_ss <- gwas_ss[,.(chr=CHROM, pos_b38=GENPOS, EA=ALLELE1, OA=ALLELE0, beta=BETA, beta_se=SE, neg_log10_p=LOG10P, EAF=A1FREQ, samples=N)]
  gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=this_gwas$Cases, total_controls=this_gwas$Controls)
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

##############################################################
# Filter the 18 G&H quantitative trait ExWAS
##############################################################
for (this_prs in gwas_list[PRS %like% "GH_Exome" & is.na(Cases), PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".*\n", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/GenesAndHealth/ExWAS/2024_05_08_%s.residual_GNH_singlevariantExWAS_%s.residual.regenie.gz", this_pheno_fname, this_pheno_fname))
  gwas_ss <- gwas_ss[,.(chr=CHROM, pos_b38=GENPOS, EA=ALLELE1, OA=ALLELE0, beta=BETA, beta_se=SE, neg_log10_p=LOG10P, EAF=A1FREQ, samples=N)]
  gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

##############################################################
# Filter the 11 G&H case-control PheWAS
##############################################################
for (this_prs in gwas_list[`PubMed ID` == 31504546 & !is.na(Cases), PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".*\n", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/GenesAndHealth/PheWAS/2025_05_23_%s_singlevariant51kGSA-TOPMEDr3-GWAS_%s.regenie.gz", this_pheno_fname, this_pheno_fname))
  gwas_ss <- gwas_ss[,.(chr=CHROM, pos_b38=GENPOS, EA=ALLELE1, OA=ALLELE0, beta=BETA, beta_se=SE, neg_log10_p=LOG10P, EAF=A1FREQ, samples=N)]
  gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=this_gwas$Cases, total_controls=this_gwas$Controls)
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}


##############################################################
# Filter the 4 G&H quantitative trait PheWAS
##############################################################
for (this_prs in gwas_list[`PubMed ID` == 31504546 & is.na(Cases), PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".*\n", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/GenesAndHealth/PheWAS/2025_05_13_%s_singlevariant51kGSA-TOPMEDr3-GWAS_%s.regenie.gz", this_pheno_fname, this_pheno_fname))
  gwas_ss <- gwas_ss[,.(chr=CHROM, pos_b38=GENPOS, EA=ALLELE1, OA=ALLELE0, beta=BETA, beta_se=SE, neg_log10_p=LOG10P, EAF=A1FREQ, samples=N)]
  gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}
