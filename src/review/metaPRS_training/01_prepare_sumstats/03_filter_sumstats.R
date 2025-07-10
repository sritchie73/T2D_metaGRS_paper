library(data.table)
source("src/functions/filter_sumstats.R")

system("mkdir -p data/filtered_sumstats/filtered_gwas", wait=TRUE)

# For each GWAS, filter to set of HapMap variants above, orienting to common effect allele,
# and output in LDpred2 format.
gwas_list <- fread("data/gwas_summary_stats/curated_gwas_table.txt", na.strings="-")
gwas_list[, Samples := as.integer(gsub(",", "", Samples))]
gwas_list[, Cases := as.integer(gsub(",", "", Cases))]
gwas_list[, Controls := as.integer(gsub(",", "", Controls))]

##
## Start with GWASs that are either not in the GWAS Catalog, or have not been harmonized
## by the GWAS Catalog into their standard format
##

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
# Filter the 11 qunatitative CKB summary statistics
##############################################################
for (this_prs in gwas_list[`GWAS catalog accession or other download source` %like% 'pheweb.ckbiobank.org' & is.na(Cases), PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".* ", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/CKB/phenocode-%s.tsv.gz", this_pheno_fname))
  if ("chrom" %in% names(gwas_ss)) {
		gwas_ss <- gwas_ss[,.(chr=chrom, pos_b38=pos, EA=alt, OA=ref, beta, beta_se=sebeta, neg_log10_p=-log10(pval), EAF=af)]
  } else {
    # SBP and PP have different column headings to other 9 traits for some reason
    gwas_ss <- gwas_ss[,.(chr=chromosome, pos_b38=base_pair_location, EA=effect_allele, OA=other_allele, beta, beta_se=standard_error,
                          neg_log10_p=-log10(p_value), EAF=effect_allele_frequency)]
  }
  gwas_ss <- filter_sumstats(gwas_ss, type="continuous", total_samples=this_gwas$Samples)
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
} 

##############################################################
# Filter the 7 case-control GWAS CKB summary statistics
##############################################################
for (this_prs in gwas_list[`GWAS catalog accession or other download source` %like% 'pheweb.ckbiobank.org' & !is.na(Cases), PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".* ", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/CKB/phenocode-%s.tsv.gz", this_pheno_fname))
  gwas_ss <- gwas_ss[,.(chr=chrom, pos_b38=pos, EA=alt, OA=ref, beta, beta_se=sebeta, neg_log10_p=-log10(pval), EAF=af)]
  gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_samples=this_gwas$Samples, total_cases=this_gwas$Cases, total_controls=this_gwas$Controls)
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

#####################################################################################
# Filter GWAS summary statistics for T2D available through DIAGRAM website
#####################################################################################
gwas_ss <- fread("data/gwas_summary_stats/DIAGRAM/Mahajan.NatGenet2018b.T2D-noUKBB.European.zip")
gwas_ss <- gwas_ss[,.(chr=Chr, pos_b37=Pos, EA, OA=NEA, beta=Beta, beta_se=SE, EAF, neg_log10_p=-log10(Pvalue))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=55005, total_controls=400308)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/filtered_gwas/T2D_2018_no_UKB.txt.gz")

##
## Moving on to GWASs with summary stats available through the GWAS Catalog, but which
## haven't been harmonized by them into their standard format
##

#####################################################################################
# Filter 18 GWAS summary statistics bundled into GCST008673 on the GWAS Catalog
#####################################################################################
for (this_prs in gwas_list[`PubMed ID` == 31367044, PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".* ", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/GWAS_Catalog/GCST008673/LockeAE_prePMID_%s_sex-combined.gz", this_pheno_fname))
  gwas_ss[, c("REF", "ALT") := tstrsplit(gsub("_.*", "", gsub(".*:[0-9]*_", "", MARKER_ID)), "/")]
  gwas_ss[, EAF := AC / (NS*2)]
  gwas_ss <- gwas_ss[,.(chr=CHROM, pos_b37=BEG, EA=ALT, OA=REF, beta=BETA, beta_se=SEBETA, neg_log10_p=-log10(PVALUE), EAF, samples=NS)] 
  gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

#####################################################################################
# Filter Psychiatric Genomics Consortium GWAS bip2021_noUKBB
#####################################################################################
gwas_ss <- fread("data/gwas_summary_stats/PGC/bip2021_noUKBB/daner_bip_pgc3_nm_noukbiobank.gz", fill=TRUE)
gwas_ss[is.na(BP), names(gwas_ss)[1:12] := tstrsplit(CHR, " ")] # 88 malformed rows??
gwas_ss[, samples := Nca + Nco]
gwas_ss[, EAF := (FRQ_A_40463 * Nca + FRQ_U_313436 * Nco)/(samples)] # Frequency reported for cases and controls separately
gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=BP, EA=A1, OA=A2, beta=log(OR), beta_se=SE, neg_log10_p=-log10(P), EAF, samples, cases=Nca, controls=Nco)]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=40463, total_controls=313436, total_samples=353899)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/filtered_gwas/Bipolar_PCG_2021_no_UKB.txt.gz")

#####################################################################################
# Filter Psychiatric Genomics Consortium GWAS bip2019
#####################################################################################
gwas_ss <- fread("data/gwas_summary_stats/PGC/daner_PGC_BIP32b_mds7a_0416a.gz")
gwas_ss[, samples := Nca + Nco]
gwas_ss[, EAF := (FRQ_A_20352 * Nca + FRQ_U_31358 * Nco)/(samples)] # Frequency reported for cases and controls separately
gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=BP, EA=A1, OA=A2, beta=log(OR), beta_se=SE, neg_log10_p=-log10(P), EAF, samples, cases=Nca, controls=Nco)]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=20352, total_controls=31358, total_samples=51710)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/filtered_gwas/Bipolar_PCG_2019.txt.gz")

#####################################################################################
# Filter Psychiatric Genomics Consortium GWAS mdd2025
#####################################################################################
gwas_ss <- fread("data/gwas_summary_stats/PGC/mdd2025/daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz")
gwas_ss[, samples := Nca + Nco]
gwas_ss[, EAF := (FRQ_A_357636 * Nca + FRQ_U_1281936 * Nco)/(samples)] # Frequency reported for cases and controls separately
gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=BP, EA=A1, OA=A2, beta=log(OR), beta_se=SE, neg_log10_p=-log10(P), EAF, samples, cases=Nca, controls=Nco)]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=357636, total_controls=1281936, total_samples=1639572)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/filtered_gwas/Depression_PCG_2025_no_UKB.txt.gz")

#####################################################################################
# Filter 5 Psychiatric Genomics Consortium GWASs in study mdd2023diverse
#####################################################################################
for (this_prs in gwas_list[PRS %like% "Depression_PCG_2023", PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".*\n", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/PGC/mdd2023diverse/%s.gz", this_pheno_fname))
  gwas_ss[, EA := toupper(EA)]
  gwas_ss[, NEA := toupper(NEA)]
  gwas_ss[, samples := Ncase + Ncontrol]
  gwas_ss <- gwas_ss[,.(chr=Chromosome, pos_b37=Position, EA, OA=NEA, beta=logOR, beta_se=SE, neg_log10_p=-log10(P), EAF, cases=Ncase, controls=Ncontrol, samples)]
  gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=this_gwas$Cases, total_controls=this_gwas$Controls, total_samples=this_gwas$Samples)
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

#####################################################################################
# Filter 5 Psychiatric Genomics Consortium GWASs in study scz2022
#####################################################################################
for (this_prs in gwas_list[PRS %like% "Schizophrenia_PCG_2022", PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  this_pheno_fname <- this_gwas[, gsub(".*\n", "", `GWAS catalog accession or other download source`)]
  gwas_ss <- fread(sprintf("data/gwas_summary_stats/PGC/scz2022/%s", this_pheno_fname))

  # Per-SNP case and control counts missing in African American and Latino ancestry GWASs
  if (!("NCAS" %in% names(gwas_ss))) {
    gwas_ss[, NCAS := this_gwas$Cases]
    gwas_ss[, NCON := this_gwas$Controls]
  }

	gwas_ss[, samples := NCAS + NCON]
	gwas_ss[, EAF := (FCAS * NCAS + FCON * NCON)/(samples)] # Frequency reported for cases and controls separately
  gwas_ss <- gwas_ss[,.(chr=CHROM, pos_b37=POS, EA=A1, OA=A2, beta=BETA, beta_se=SE, neg_log10_p=-log10(PVAL), EAF, cases=NCAS, controls=NCON, samples)]
  gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=this_gwas$Cases, total_controls=this_gwas$Controls, total_samples=this_gwas$Samples)
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}


#####################################################################################
# Filter Psychiatric Genomics Consortium GWAS scz2018clozuk
#####################################################################################
gwas_ss <- fread("data/gwas_summary_stats/PGC/scz2018clozuk/CLOZUK_PGC2noclo.METAL.assoc.dosage.fix.gz")
gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=BP, EA=A1, OA=A2, beta=log(OR), beta_se=SE, neg_log10_p=-log10(P))] # No EAF column
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=40675, total_controls=64643)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/filtered_gwas/Schizophrenia_PCG_2018.txt.gz")

#####################################################################################
# Filter the 16 UGR GWASs (PMID: 31675503)
#####################################################################################
for (this_prs in gwas_list[`PubMed ID` == 31675503, PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  gcstid <- this_gwas[,`GWAS catalog accession or other download source`]
  gwas_fname <- list.files(pattern="*.txt.gz", path=sprintf("data/gwas_summary_stats/GWAS_Catalog/PMID_31675503/%s/", gcstid), full.names=TRUE)
  gwas_ss <- fread(gwas_fname, fill=TRUE)

  gwas_ss[, c("chr", "pos_b37", "EA", "OA") := tstrsplit(snpid, ":")]
  gwas_ss[, pos_b37 := as.integer(pos_b37)]

  # Harmonize meta-analysis study specific columns
  if (this_prs == "HbA1c_UGR") gwas_ss[, no_DCC := 7526 - no_uganda] # column missing data for some reason

  if ("af_uganda" %in% names(gwas_ss)) {
		gwas_ss[is.na(af_uganda), c("af_uganda", "no_uganda") := 0]
  } else {
    gwas_ss[, c("af_uganda", "no_uganda") := 0]
  }

  if ("af_DCC" %in% names(gwas_ss)) {
		gwas_ss[is.na(af_DCC), c("af_DCC", "no_DCC") := 0]
  } else {
    gwas_ss[, c("af_DCC", "no_DCC") := 0]
  }

  if ("af_DDS" %in% names(gwas_ss)) {
		gwas_ss[is.na(af_DDS), c("af_DSS", "no_DSS") := 0]
  } else {
    gwas_ss[, c("af_DDS", "no_DSS") := 0]
  }

  if ("af_AADM" %in% names(gwas_ss)) {
		gwas_ss[is.na(af_AADM), c("af_AADM", "no_AADM") := 0]
  } else {
    gwas_ss[, c("af_AADM", "no_AADM") := 0]
  }

  # Compute overall allele frequency across meta-analysis
  gwas_ss[, samples := no_uganda + no_DCC + no_DDS + no_AADM]
  gwas_ss[, EAF := (af_uganda*no_uganda + af_DCC*no_DCC + af_DDS*no_DDS + af_AADM*no_AADM)/samples]

  gwas_ss <- gwas_ss[, .(chr, pos_b37, EA, OA, beta=beta_re, beta_se=se_re, neg_log10_p=-log10(pval_re2), EAF, samples)]
  gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

#####################################################################################
# Filter the 3 Biobank Japan Smoking Status GWAS (PMID: 31089300)
#####################################################################################
for (this_prs in gwas_list[`PubMed ID` == 31089300 & !is.na(Cases), PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  gcstid <- this_gwas[,`GWAS catalog accession or other download source`]
  gwas_fname <- list.files(pattern="*_autosome_.*.txt.gz", path=sprintf("data/gwas_summary_stats/GWAS_Catalog/PMID_31089300/%s/", gcstid), full.names=TRUE)
  gwas_ss <- fread(gwas_fname, fill=TRUE)
  gwas_ss[, BETA := bolt_lmm_fix(BETA, this_gwas$Cases, this_gwas$Controls)] # Conversion of BOLT-LMM linear model beta to log odds ratio 
  gwas_ss[, SE := bolt_lmm_fix(SE, this_gwas$Cases, this_gwas$Controls)]
  gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=POS, EA=A1, OA=A2, beta=BETA, beta_se=SE, EAF=A1Frq, neg_log10_p=-log10(P))]
  gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=this_gwas$Cases, total_controls=this_gwas$Controls)
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

#####################################################################################
# Filter the 3 Biobank Japan Smoking Severity GWAS (PMID: 31089300)
#####################################################################################
for (this_prs in gwas_list[`PubMed ID` == 31089300 & is.na(Cases), PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  gcstid <- this_gwas[,`GWAS catalog accession or other download source`]
  gwas_fname <- list.files(pattern="*_autosome_.*.txt.gz", path=sprintf("data/gwas_summary_stats/GWAS_Catalog/PMID_31089300/%s/", gcstid), full.names=TRUE)
  gwas_ss <- fread(gwas_fname, fill=TRUE)
  gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=POS, EA=A1, OA=A2, beta=BETA, beta_se=as.numeric(SE), EAF=A1Frq, neg_log10_p=-log10(P))] # ignore warning message about NA coercion, one SNP has a SE of '.' for CigsPerDayFemale_BBJ
  gwas_ss <- filter_sumstats(gwas_ss, type="continuous", total_samples=this_gwas$Samples)
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

#####################################################################################
# Filter the 4 ancestry specific BMI ExWAS (PMID: 29273807)
#####################################################################################
for (this_prs in gwas_list[`PubMed ID` == 29273807, PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  gcstid <- this_gwas[,`GWAS catalog accession or other download source`]
  gwas_fname <- list.files(pattern="*.txt.gz", path=sprintf("data/gwas_summary_stats/GWAS_Catalog/PMID_29273807/%s/", gcstid), full.names=TRUE)
  gwas_ss <- fread(gwas_fname, fill=TRUE, na.strings=c("NA", "-"))
  
  # Process weird MAF format
  setnames(gwas_ss, names(gwas_ss)[names(gwas_ss) %like% "MAF" & !(names(gwas_ss) %like% "ExAC")], "MAF_pairs")
  pairs <- strsplit(gwas_ss$MAF_pairs, ",")
  pairs <- sapply(pairs, unique)
  pairs <- sapply(pairs, function(x) {  if(length(x) > 1) { x[!(x %like% "-")] } else { x } })
  gwas_ss[, MAF_pair := unlist(pairs)]
  gwas_ss[, c("MA", "MAF") := tstrsplit(MAF_pair, ":", type.convert=TRUE)]
  gwas_ss[MA == ALT, EAF := MAF]
  gwas_ss[is.na(EAF) & MA == flip_strand(ALT), EAF := MAF]
  gwas_ss[is.na(EAF) & MA == REF, EAF := 1 - MAF]
  gwas_ss[is.na(EAF) & MA == flip_strand(REF), EAF := 1 - MAF]

  gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=POS, EA=ALT, OA=REF, beta, beta_se=se, neg_log10_p=-log10(Pvalue), EAF)]
  gwas_ss <- filter_sumstats(gwas_ss, type="continuous", total_samples=this_gwas$Samples)
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

#####################################################################################
# Filter the 4 multi-ancestry Lipid ExWAS (PMID: 29083408)
#####################################################################################
for (this_prs in gwas_list[`PubMed ID` == 29083408, PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  gcstid <- this_gwas[,`GWAS catalog accession or other download source`]
  gwas_fname <- list.files(pattern="*.txt", path=sprintf("data/gwas_summary_stats/GWAS_Catalog/PMID_29083408/%s/", gcstid), full.names=TRUE)
  gwas_fname <- gwas_fname[!(gwas_fname %like% "md5sum.txt")] 
  gwas_ss <- fread(gwas_fname)

  gwas_ss[, c("chr", "pos_b37") := tstrsplit(SNP_hg19, ":", type.convert=TRUE)]
  setnames(gwas_ss, names(gwas_ss)[names(gwas_ss) %like% "Freq.A1"], "EAF")

  gwas_ss <- gwas_ss[,.(chr, pos_b37, EA=A1, OA=A2, beta=Beta, beta_se=SE, neg_log10_p=-log10(`P-value`), EAF, samples=N)]
  gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

#####################################################################################
# Filter the 4 GERA Lipid GWAS (PMID: 29507422)
#####################################################################################
for (this_prs in gwas_list[`PubMed ID` == 29507422, PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  gcstid <- this_gwas[,`GWAS catalog accession or other download source`]
  gwas_fname <- list.files(pattern="*.tsv.gz", path=sprintf("data/gwas_summary_stats/GWAS_Catalog/PMID_29507422/%s/", gcstid), full.names=TRUE)
  gwas_ss <- fread(gwas_fname)
  setnames(gwas_ss, c("Allele 1", "Allele 2", "Effect allele (EA)", "Effect allele frequency (EAF)", "Sample size", "Estimate Effect", "P value"),
           c("a1", "a2", "EA", "EAF", "samples", "beta", "Pval"))
  
  gwas_ss[EA == a1, OA := a2]
  gwas_ss[EA == a2, OA := a1]

  gwas_ss <- gwas_ss[,.(chr=chromosome, pos_b37=position, EA, OA, beta, neg_log10_p=-log10(Pval), EAF, samples)]
  gwas_ss[, beta_se := lm_se(beta, neg_log10_p, samples, 13)] # 13 covariates: age^3, sex, BMI^3, and 10 PCs 

  gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

#####################################################################################
# Filter the 3 GERA Blood Pressure GWAS (PMID: 27841878)
#####################################################################################
for (this_prs in gwas_list[`PubMed ID` == 27841878, PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  gcstid <- this_gwas[,`GWAS catalog accession or other download source`]
  gwas_fname <- list.files(pattern="*.tsv.gz", path=sprintf("data/gwas_summary_stats/GWAS_Catalog/PMID_27841878/%s/", gcstid), full.names=TRUE)
  gwas_ss <- fread(gwas_fname)
  setnames(gwas_ss, c("Allele 1", "Allele 2", "Effect allele (EA)", "Effect allele frequency (EAF)", "Sample size", "Estimate Effect", "P value"),
           c("a1", "a2", "EA", "EAF", "samples", "beta", "Pval"))
  
  gwas_ss[EA == a1, OA := a2]
  gwas_ss[EA == a2, OA := a1]

  gwas_ss <- gwas_ss[,.(chr=chromosome, pos_b37=position, EA, OA, beta, neg_log10_p=-log10(Pval), EAF, samples)]
  gwas_ss[, beta_se := lm_se(beta, neg_log10_p, samples, 14)] # 14 covariates: age, age^2, sex, BMI, and 10 PCs 

  gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

#####################################################################################
# Filter the 2017 DIAGRAM T2D GWAS (PMID: 28566273)
#####################################################################################
gwas_ss <- fread("data/gwas_summary_stats/GWAS_Catalog/PMID_28566273/METAANALYSIS_DIAGRAM_SE1.txt")
gwas_ss[, c("chr", "pos_b37") := tstrsplit(`Chr:Position`, ":", type.convert=TRUE)]
gwas_ss <- gwas_ss[,.(chr, pos_b37, EA=Allele1, OA=Allele2, beta=Effect, beta_se=StdErr,
                      samples=TotalSampleSize, neg_log10_p=-log10(`P-value`))] # No EAF column
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=26676, total_controls=132532, total_samples=159208)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/filtered_gwas/T2D_2017.txt.gz")

#####################################################################################
# Filter the 2019 African T2D GWAS (PMID: 31049640)
#####################################################################################
gwas_ss <- fread("data/gwas_summary_stats/GWAS_Catalog/PMID_31049640/ChenJ_31049640")
gwas_ss[, c("chr", "pos_b37") := tstrsplit(MarkerName, ":", keep=1:2, type.convert=TRUE)]
gwas_ss <- gwas_ss[,.(chr, pos_b37, EA=toupper(Allele1), OA=toupper(Allele2), beta=Effect, beta_se=StdErr, n_eff=Weight,
                      EAF=Freq1, neg_log10_p=-log10(`P-value`))]
gwas_ss <- filter_sumstats(gwas_ss)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/filtered_gwas/T2D_Africa.txt.gz")

#####################################################################################
# Filter the 2019 BBJ T2D GWAS (PMID: 30718926)
#####################################################################################
gwas_ss <- fread("data/gwas_summary_stats/GWAS_Catalog/PMID_30718926/BBJ_BetaBased1.MAF_001.AtLeast2studies.AllChr.txt.gz")
gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=POS, EA=ALT, OA=REF, beta=BETA, beta_se=SE, samples=N,
                      EAF=Frq, neg_log10_p=-log10(P))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=6967, total_controls=49670)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/filtered_gwas/T2D_Japan.txt.gz")

#####################################################################################
# Filter the 2018 PCOS GWAS (PMID: 30566500)
#####################################################################################
gwas_ss <- fread("data/gwas_summary_stats/GWAS_Catalog/PMID_30566500/PCOS_top_10000_SNPS.txt") # includes 23andMe, hence only top 10k SNPs
gwas_ss[, c("chr", "pos_b37") := tstrsplit(MarkerName, ":", keep=1:2, type.convert=TRUE)]
gwas_ss <- gwas_ss[, .(chr, pos_b37, EA=toupper(Allele1), OA=toupper(Allele2), EAF=Freq1, beta=Effect, beta_se=StdErr, samples=TotalSampleSize, neg_log10_p=-log10(Pvalue))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=10074, total_controls=103164)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/filtered_gwas/PCOS_2018.txt.gz")

#####################################################################################
# Filter the 2016 Educational Attainment GWAS (PMID: 27225129)
#####################################################################################
gwas_ss <- fread("data/gwas_summary_stats/GWAS_Catalog/PMID_27225129/Okbay_27225129-EduYears_Discovery_5000.txt")
gwas_ss <- gwas_ss[, .(chr=CHR, pos_b37=POS, EA=A1, OA=A2, EAF, beta=Beta, beta_se=SE, neg_log10_p=-log10(Pval))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous", total_samples=280007)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/filtered_gwas/EduYears_EUR_noUKB.txt.gz")

##
## Now onto GWASs that have harmonized summary statistics in the GWAS Catalog
##

#####################################################################################
# Filter the 58 MVP quantitative trait GWASs (PMID: 39024449)
#####################################################################################
for (this_prs in gwas_list[`PubMed ID` == 39024449 & is.na(Cases), PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  gcstid <- this_gwas[,`GWAS catalog accession or other download source`]
  gwas_fname <- list.files(pattern=gcstid, path="data/gwas_summary_stats/GWAS_Catalog/Harmonized/", full.names=TRUE)
  gwas_ss <- fread(gwas_fname)

  gwas_ss <- gwas_ss[, .(chr=chromosome, pos_b38=base_pair_location, EA=effect_allele, OA=other_allele, beta, beta_se=standard_error,
                         neg_log10_p=-log10(p_value), EAF=effect_allele_frequency, samples=n)]

  gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

#####################################################################################
# Filter the 38 MVP case/control GWASs (PMID: 39024449)
#####################################################################################
for (this_prs in gwas_list[`PubMed ID` == 39024449 & !is.na(Cases), PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  gcstid <- this_gwas[,`GWAS catalog accession or other download source`]
  gwas_fname <- list.files(pattern=gcstid, path="data/gwas_summary_stats/GWAS_Catalog/Harmonized/", full.names=TRUE)
  gwas_ss <- fread(gwas_fname)
  
  gwas_ss <- gwas_ss[, .(chr=chromosome, pos_b38=base_pair_location, EA=effect_allele, OA=other_allele, 
                         beta=log(odds_ratio), beta_se=se_from_ci(log(odds_ratio), log(ci_lower), log(ci_upper)),
                         neg_log10_p=-log10(p_value), EAF=effect_allele_frequency, samples=n, cases=num_cases, controls=num_controls)]

  gwas_ss <- filter_sumstats(gwas_ss, type="case/control")
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}

#####################################################################################
# Filter the 29 BBJ+TWB meta-analyses quantitative trait GWASs (PMID: 38116116)
#####################################################################################
for (this_prs in gwas_list[`PubMed ID` == 38116116, PRS]) {
  this_gwas <- gwas_list[PRS == this_prs]
  gcstid <- this_gwas[,`GWAS catalog accession or other download source`]
  gwas_fname <- list.files(pattern=gcstid, path="data/gwas_summary_stats/GWAS_Catalog/Harmonized/", full.names=TRUE)
  gwas_ss <- fread(gwas_fname)

  gwas_ss <- gwas_ss[, .(chr=chromosome, pos_b38=base_pair_location, EA=effect_allele, OA=other_allele, beta, beta_se=standard_error,
                         neg_log10_p=-log10(p_value), EAF=effect_allele_frequency, samples=n)]

  gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
  fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file=sprintf("data/filtered_sumstats/filtered_gwas/%s.txt.gz", this_prs))
}




