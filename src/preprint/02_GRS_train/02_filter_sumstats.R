library(data.table)
library(foreach)
source("src/functions/flip_strand.R")

# Create directories for each GWAS
for (dd in list.dirs("data/gwas_summary_stats", full.names=FALSE, recursive=FALSE)) {
  system(sprintf("mkdir -p data/filtered_sumstats/%s", dd), wait=TRUE)
}

# For each GWAS, filter to set of HapMap variants above, orienting to common effect allele,
# and output in LDpred2 format.

# Note n_eff for case control studies is computed as 4/(1/cases + 1/controls) and sample size 
# for continuous traits. Case/control numbers are obtained directly from the study/GWAS Catalog,
# then also weighted by total sample size (if provided and varying by SNP in the summary stats)
# or total number of contributing studies (if provided for each SNP in the absence of per-SNP 
# sample size). 
n_eff <- function(cases, controls, per_snp_N=1, total_N=1) {
  cases <- cases / total_N * per_snp_N
  controls <- controls / total_N * per_snp_N
  4 / (1/cases + 1/controls)
}

# Some GWAS are missing the standard error, for additive effects in linear models we can 
# back-calculate the T-statistic using the sample size and P-value
lm_se <- function(beta, neg_log10_p, samples, n_covar=0) {
  df <- samples - 2 - n_covar # in practice for any sufficiently large N for GWAS, subtracting right number of covariates is irrelevant
  t_stat <- qt(10^(-neg_log10_p)/2, df = df)
  abs(beta)/abs(t_stat)
}

# Correction of BOLT-LMM beta and se for case %
# https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html#x1-470008
bolt_lmm_fix <- function(x, cases, controls) {
  u <- cases / (cases + controls)
  x  / (u * (1-u))
}

# Stop bad fread when running out of space in /tmp
fread <- function(...) {
  data.table::fread(..., tmpdir="tmp")
}

# Filter and format sumstats orienting to UKB alleles
filter_sumstats <- function(gwas_ss, type, total_samples, total_cases, total_controls) {
  # Get chromosome and position if missing
  if ("rsid" %in% names(gwas_ss)) {
    gwas_ss[varset, on = .(rsid=rsid_1000G), c("chr", "pos_b37") := .(i.chr, i.pos_b37)]
  }

  # Filter to autosomal chromosomes
  if (typeof(gwas_ss$chr) != "integer") {
		gwas_ss <- gwas_ss[chr %in% 1:22]
		gwas_ss[, chr := as.integer(chr)]
  }
 
  # Get positions on build 37 if necessary
	if ("pos_b38" %in% names(gwas_ss)) {
		gwas_ss[varset, on = .(chr, pos_b38), pos_b37 := i.pos_b37]
		gwas_ss <- gwas_ss[!is.na(pos_b37)]
	} else if ("pos_b36" %in% names(gwas_ss)) {
		gwas_ss[varset, on = .(chr, pos_b36), pos_b37 := i.pos_b37]
		gwas_ss <- gwas_ss[!is.na(pos_b37)]
	} 

  # Add in EAF column if missing to simplify downstream code
  if (!("EAF" %in% names(gwas_ss))) {
    gwas_ss[, EAF := NA]
  }

  # Require finite and non-missing weights and standard errors
  gwas_ss <- gwas_ss[is.finite(beta)]
  gwas_ss <- gwas_ss[is.finite(beta_se)]

  # Filter to variants in the candidate variant set by chromosome and position
  gwas_ss <- gwas_ss[varset[,.(chr, pos_b37)], on = .(chr, pos_b37), nomatch=0]

  # Add in rsID from the candidate varset
  gwas_ss[varset, on = .(chr, pos_b37), rsid := i.rsid_1000G]

	# Make sure we can also match by allele
	gwas_ss[, allele_match := FALSE]
  gwas_ss[varset, on = .(chr, pos_b37, EA=effect_allele, OA=other_allele), c("allele_match", "oriented", "flipped") := .(TRUE, TRUE, FALSE)]
  gwas_ss[varset, on = .(chr, pos_b37, EA=other_allele, OA=effect_allele), c("allele_match", "oriented", "flipped") := .(TRUE, FALSE, FALSE)]

	# Some SNPs may be on the opposite strand in the GWAS, which we detect and fix here:
  gwas_ss[!(allele_match), c("EA", "OA", "flipped") := .(flip_strand(EA), flip_strand(OA), TRUE)]
  gwas_ss[varset, on = .(chr, pos_b37, EA=effect_allele, OA=other_allele), c("allele_match", "oriented") := .(TRUE, TRUE)]
  gwas_ss[varset, on = .(chr, pos_b37, EA=other_allele, OA=effect_allele), c("allele_match", "oriented") := .(TRUE, FALSE)]

	# Remove any variants which matched on chromosome and position (or rsid) but not on alleles even after checking for strand mismatch
	gwas_ss <- gwas_ss[!gwas_ss[!(allele_match)], on = .(chr, pos_b37)] # anti-join in case multi-allelic variants with match for some alleles

	# Drop variants that were multi-allelic in the gwas
	mult <- gwas_ss[,.N,by=.(chr, pos_b37)][N > 1]
	gwas_ss <- gwas_ss[!mult, on =.(chr, pos_b37)]

	# For strand ambiguous alleles (A/T or G/C SNPs) we will assume same strand orientation as (majority) rest of sumstats 
  pct_flipped <- gwas_ss[EA != flip_strand(OA), sum(flipped)/.N]
  if (pct_flipped < 0.5) {
    gwas_ss[(flipped) & EA == flip_strand(OA), c("EA", "OA", "EAF", "beta", "oriented", "flipped") := .(OA, EA, 1 - EAF, -beta, !oriented, FALSE)]
  } else {
    gwas_ss[!(flipped) & EA == flip_strand(OA), c("EA", "OA", "EAF", "beta", "oriented", "flipped") := .(OA, EA, 1 - EAF, -beta, !oriented, TRUE)]
  }

  # Now fix orientation of effect alleles to match the candidate varset
  gwas_ss[!(oriented), c("EA", "OA", "EAF", "beta", "oriented") := .(OA, EA, 1-EAF, -beta, TRUE)]

  # If the GWAS has EAF information, we can also cross-check strand ambiguous SNPs with EAF in UKB
  bad <- gwas_ss[EA == flip_strand(OA) & EAF > 0.42 & EAF < 0.58]
  gwas_ss <- gwas_ss[!bad, on = .(chr, pos_b37)]
    
	# Add in UKB allele frequency
	gwas_ss[varset, on=.(chr, pos_b37), EAF_UKB := i.EAF_UKB]


  # Check cases where EAF does not agree with UKB EAF for strand ambiguous alleles
  # In these cases we need to (1) flip the strand of the alleles, and (2) reverse the orientation
  # of effect/other alleles. This (1) basically means we don't have to change the 'EA' and 'OA'
  # columns because the two operations cancel out, but (2) we have to flip the 'beta' and 'EAF'
	gwas_ss[EA == flip_strand(OA) & ( (EAF < 0.5 & EAF_UKB > 0.5) | (EAF > 0.5 & EAF_UKB < 0.5) ),
					c("beta", "EAF", "flipped") := .(-beta, 1-EAF, !flipped)]

  # Compute effective sample size if not provided
  if (!missing(type) && type == "continuous") {
    if ("samples" %in% names(gwas_ss)) { 
      gwas_ss[, n_eff := samples]
    } else {
      gwas_ss[, n_eff := total_samples]
    } 
  } else if (!missing(type) && type == "case/control") {
    if (!("cases" %in% names(gwas_ss))) {
      gwas_ss[, cases := total_cases]
    }
    if (!("controls" %in% names(gwas_ss))) {
      gwas_ss[, controls := total_controls]
    }
    if ("samples" %in% names(gwas_ss)) { 
      if (missing(total_samples)) {
        total_samples <- gwas_ss[, max(samples)]
      }
      gwas_ss[, n_eff := n_eff(cases, controls, samples, total_samples)]
    } else {
      gwas_ss[, n_eff := n_eff(cases, controls)]
    }
  } 

  # Extract columns of interest and return
  gwas_ss[, .(chr, rsid, pos=pos_b37, a1=EA, a0=OA, a1freq=EAF, beta, beta_se, neg_log10_p, n_eff)]
}

# Load candidate variant set
varset <- fread("data/filtered_sumstats/filtered_oriented_SNPs.txt")

####################
# Bipolar GWAS
####################
gwas_ss <- fread("data/gwas_summary_stats/Bipolar/daner_PGC_BIP32b_mds7a_0416a.gz")
gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=BP, EA=A1, OA=A2, 
                      EAF=(FRQ_A_20352*20352 + FRQ_U_31358*31358)/(20352 + 31358), 
                      beta=log(OR), beta_se=SE, neg_log10_p=-log10(P),
                      cases=Nca, controls=Nco)]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/Bipolar/gwas_filtered_oriented_SNPs.txt.gz")

##############
# BMI GWAS
##############
gwas_ss <- fread("data/gwas_summary_stats/BMI/25673413-GCST002783-EFO_0004340.h.tsv.gz")
gwas_ss <- gwas_ss[, .(chr=hm_chrom, pos_b38=hm_pos, EA=hm_effect_allele, OA=hm_other_allele, 
                       EAF=eaf_ref, beta=hm_beta, beta_se=standard_error, neg_log10_p=-log10(p_value),
                       samples=n)]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/BMI/gwas_filtered_oriented_SNPs.txt.gz")

#############
# CAD GWAS
#############
# Here the effect sample size per allele isn't provided, but the number of studies contributing to the meta analysis 
# (for each variant) is, so we weight n_eff by the number of contributing studies instead of number of contributing 
# samples
gwas_ss <- fread("data/gwas_summary_stats/CAD/26343387-GCST003116-EFO_0000378-build37.f.tsv.gz")
gwas_ss <- gwas_ss[, .(chr=chromosome, pos_b37=base_pair_location, EA=effect_allele, OA=other_allele, 
                       EAF=effect_allele_frequency, beta, beta_se=standard_error, neg_log10_p=-log10(p_value),
                       samples=nstudy)]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=61289, total_controls=126310)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/CAD/gwas_filtered_oriented_SNPs.txt.gz")

####################
# CAD Japan GWAS
####################
gwas_ss <- fread("data/gwas_summary_stats/CAD_Japan/BBJCAD_2020.sumstats.gz")
gwas_ss <- gwas_ss[, .(chr=CHR, pos_b37=POS, EA=ALT, OA=REF, EAF=AAF, beta=BETA, beta_se=SE, neg_log10_p=-log10(P), samples=N)]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=25982, total_controls=142336)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/CAD_Japan/gwas_filtered_oriented_SNPs.txt.gz")

############
# DBP GWAS
############
gwas_ss <- fread("data/gwas_summary_stats/DBP/gera-dbp.tsv.gz")
gwas_ss <- gwas_ss[, .(chr=chromosome, pos_b37=position, EA=`Effect allele (EA)`, A1=`Allele 1`, A2=`Allele 2`, 
                       EAF=`Effect allele frequency (EAF)`, samples=`Sample size`, beta=`Estimate Effect`, 
                       neg_log10_p=-log10(`P value`))]

# Determine non-effect allele
gwas_ss[, OA := ifelse(EA == A1, A2, A1)]

# Compute standard error since not provided 
gwas_ss[, beta_se := lm_se(beta, neg_log10_p, samples, 14)] # covariates were age + age^2 + BMI + sex + 10 PCs

# Filter
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")

# Write out
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/DBP/gwas_filtered_oriented_SNPs.txt.gz")

############
# SBP GWAS
############
# Note same study as above - so same GWAS summary stats format and issues
gwas_ss <- fread("data/gwas_summary_stats/SBP/gera-sbp.tsv.gz")
gwas_ss <- gwas_ss[, .(chr=chromosome, pos_b37=position, EA=`Effect allele (EA)`, A1=`Allele 1`, A2=`Allele 2`, 
                       EAF=`Effect allele frequency (EAF)`, samples=`Sample size`, beta=`Estimate Effect`,
                       neg_log10_p=-log10(`P value`))]
gwas_ss[, OA := ifelse(EA == A1, A2, A1)]
gwas_ss[, beta_se := lm_se(beta, neg_log10_p, samples, 14)] # covariates were age + age^2 + BMI + sex + 10 PCs
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/SBP/gwas_filtered_oriented_SNPs.txt.gz")

################
# EduYears GWAS
################
gwas_ss <- fread("data/gwas_summary_stats/EduYears/Okbay_27225129-EduYears_Main.txt.gz")
gwas_ss <- gwas_ss[, .(chr=CHR, pos_b37=POS, EA=A1, OA=A2, EAF, beta=Beta, beta_se=SE, neg_log10_p=-log10(Pval))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous", total_samples=405072)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/EduYears/gwas_filtered_oriented_SNPs.txt.gz")

#########################
# Fasting Glucose GWAS
#########################
gwas_ss <- fread("data/gwas_summary_stats/FastingGlucose/GCST90002232_buildGRCh37.tsv.gz")
gwas_ss <- gwas_ss[,.(chr=chromosome, pos_b37=base_pair_location, EA=effect_allele, OA=other_allele, 
                      EAF=effect_allele_frequency, beta, beta_se=standard_error, 
                      neg_log10_p=-log10(as.numeric(p_value)), samples=sample_size)]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/FastingGlucose/gwas_filtered_oriented_SNPs.txt.gz")

#########################
# Fasting Insulin GWAS
#########################
gwas_ss <- fread("data/gwas_summary_stats/FastingInsulin/GCST90002238_buildGRCh37.tsv.gz")
gwas_ss <- gwas_ss[,.(chr=chromosome, pos_b37=base_pair_location, EA=effect_allele, OA=other_allele, 
                      EAF=effect_allele_frequency, beta, beta_se=standard_error, 
                      neg_log10_p=-log10(p_value), samples=sample_size)]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/FastingInsulin/gwas_filtered_oriented_SNPs.txt.gz")

#########################
# HBA1C GWAS
#########################
gwas_ss <- fread("data/gwas_summary_stats/HBA1C/GCST90002244_buildGRCh37.tsv.gz")
gwas_ss <- gwas_ss[,.(chr=chromosome, pos_b37=base_pair_location, EA=effect_allele, OA=other_allele, 
                      EAF=effect_allele_frequency, beta, beta_se=standard_error, 
                      neg_log10_p=-log10(p_value), samples=sample_size)]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/HBA1C/gwas_filtered_oriented_SNPs.txt.gz")

#########################
# HDL GWAS
#########################
gwas_ss <- fread("data/gwas_summary_stats/HDL/24097068-GCST002223-EFO_0004612-build37.f.tsv.gz")
gwas_ss <- gwas_ss[,.(chr=as.integer(gsub("chr", "", chromosome)), pos_b37=base_pair_location, 
                      EA=toupper(A1), OA=toupper(A2), beta, beta_se=standard_error, samples=N,
                      EAF=eaf_ref, neg_log10_p=-log10(as.numeric(`P-value`)))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/HDL/gwas_filtered_oriented_SNPs.txt.gz")

#########################
# LDL GWAS
#########################
gwas_ss <- fread("data/gwas_summary_stats/LDL/24097068-GCST002222-EFO_0004611-build37.f.tsv.gz")
gwas_ss <- gwas_ss[,.(chr=as.integer(gsub("chr", "", chromosome)), pos_b37=base_pair_location, 
                      EA=toupper(effect_allele), OA=toupper(other_allele), beta, beta_se=standard_error, 
                      EAF=eaf_ref, neg_log10_p=-log10(as.numeric(p_value)), samples=n)]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/LDL/gwas_filtered_oriented_SNPs.txt.gz")

#########################
# Total Cholesterol GWAS
#########################
gwas_ss <- fread("data/gwas_summary_stats/TChol/24097068-GCST002221-EFO_0004574-build37.f.tsv.gz")
gwas_ss <- gwas_ss[,.(rsid=variant_id, EA=toupper(effect_allele), OA=toupper(other_allele), samples=n,
                      beta, beta_se=standard_error, EAF=eaf_ref, neg_log10_p=-log10(p_value))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/TChol/gwas_filtered_oriented_SNPs.txt.gz")

#########################
# Triglycerides GWAS
#########################
gwas_ss <- fread("data/gwas_summary_stats/Triglycerides/24097068-GCST002216-EFO_0004530-build37.f.tsv.gz")
gwas_ss <- gwas_ss[,.(rsid=variant_id, EA=toupper(effect_allele), OA=toupper(other_allele), 
                      beta, beta_se=standard_error, samples=n, EAF=eaf_ref, neg_log10_p=-log10(p_value))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/Triglycerides/gwas_filtered_oriented_SNPs.txt.gz")

##########################
# Hip Circumference GWAS
##########################
gwas_ss <- fread("data/gwas_summary_stats/HipCircumference/GIANT_2015_HIP_COMBINED_EUR.txt.gz")
gwas_ss <- gwas_ss[,.(rsid=MarkerName, EA=Allele1, OA=Allele2, beta=b, beta_se=se, samples=N,
                      EAF=FreqAllele1HapMapCEU, neg_log10_p=-log10(p))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/HipCircumference/gwas_filtered_oriented_SNPs.txt.gz")

############################################
# Hip Circumference adjusting for BMI GWAS
############################################
gwas_ss <- fread("data/gwas_summary_stats/HipCircumferenceAdjBMI/GIANT_2015_HIPadjBMI_COMBINED_EUR.txt.gz")
gwas_ss <- gwas_ss[,.(rsid=MarkerName, EA=Allele1, OA=Allele2, beta=b, beta_se=se, samples=N,
                      EAF=FreqAllele1HapMapCEU, neg_log10_p=-log10(p))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/HipCircumferenceAdjBMI/gwas_filtered_oriented_SNPs.txt.gz")

##########################
# Waist Circumference GWAS
##########################
gwas_ss <- fread("data/gwas_summary_stats/WaistCircumference/GIANT_2015_WC_COMBINED_EUR.txt.gz")
gwas_ss <- gwas_ss[,.(rsid=MarkerName, EA=Allele1, OA=Allele2, beta=b, beta_se=se, samples=N,
                      EAF=FreqAllele1HapMapCEU, neg_log10_p=-log10(p))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/WaistCircumference/gwas_filtered_oriented_SNPs.txt.gz")

############################################
# Waist Circumference adjusting for BMI GWAS
############################################
gwas_ss <- fread("data/gwas_summary_stats/WaistCircumferenceAdjBMI/GIANT_2015_WCadjBMI_COMBINED_EUR.txt.gz")
gwas_ss <- gwas_ss[,.(rsid=MarkerName, EA=Allele1, OA=Allele2, beta=b, beta_se=se, samples=N,
                      EAF=FreqAllele1HapMapCEU, neg_log10_p=-log10(p))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/WaistCircumferenceAdjBMI/gwas_filtered_oriented_SNPs.txt.gz")

##########################
# Waist-Hip ratio GWAS
##########################
gwas_ss <- fread("data/gwas_summary_stats/WHR/GIANT_2015_WHR_COMBINED_EUR.txt.gz")
gwas_ss <- gwas_ss[,.(rsid=MarkerName, EA=Allele1, OA=Allele2, beta=b, beta_se=se, samples=N,
                      EAF=FreqAllele1HapMapCEU, neg_log10_p=-log10(p))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/WHR/gwas_filtered_oriented_SNPs.txt.gz")

############################################
# Waist-Hip ratio adjusting for BMI GWAS
############################################
gwas_ss <- fread("data/gwas_summary_stats/WHRadjBMI/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz")
gwas_ss <- gwas_ss[,.(rsid=MarkerName, EA=Allele1, OA=Allele2, beta=b, beta_se=se, samples=N,
                      EAF=FreqAllele1HapMapCEU, neg_log10_p=-log10(p))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/WHRadjBMI/gwas_filtered_oriented_SNPs.txt.gz")

#########################
# HOMA-B GWAS
#########################
gwas_ss <- fread("data/gwas_summary_stats/HOMAB/20081858-GCST005180-EFO_0004469-build37.f.tsv.gz")
gwas_ss <- gwas_ss[, .(rsid=variant_id, EA=toupper(effect_allele), OA=toupper(other_allele), beta, beta_se=standard_error,
                       EAF=effect_allele_frequency, neg_log10_p=-log10(p_value))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous", total_samples=36466)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/HOMAB/gwas_filtered_oriented_SNPs.txt.gz")

#########################
# HOMA-IR GWAS
#########################
gwas_ss <- fread("data/gwas_summary_stats/HOMAIR/20081858-GCST005179-EFO_0004501-build37.f.tsv.gz")
gwas_ss <- gwas_ss[, .(rsid=variant_id, EA=toupper(effect_allele), OA=toupper(other_allele), beta, beta_se=standard_error,
                       EAF=effect_allele_frequency, neg_log10_p=-log10(p_value))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous", total_samples=37037)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/HOMAIR/gwas_filtered_oriented_SNPs.txt.gz")

###############
# Leptin GWAS
###############
gwas_ss <- fread("data/gwas_summary_stats/Leptin/GCST90007310_buildGRCh37.tsv.gz")
gwas_ss <- gwas_ss[beta != ".", .(chr=chromosome, pos_b37=base_pair_location, EA=effect_allele, OA=other_allele, beta=as.numeric(beta), 
                                  beta_se=standard_error, samples=sample_size, EAF=effect_allele_frequency, 
                                  neg_log10_p=-log10(as.numeric(p_value)))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/Leptin/gwas_filtered_oriented_SNPs.txt.gz")

#################################
# Leptin adjusting for BMI GWAS
#################################
gwas_ss <- fread("data/gwas_summary_stats/LeptinAdjBMI/GCST90007322_buildGRCh37.tsv.gz")
gwas_ss <- gwas_ss[beta != ".", .(chr=chromosome, pos_b37=base_pair_location, EA=effect_allele, OA=other_allele, beta=as.numeric(beta), 
                                  beta_se=standard_error, samples=sample_size, EAF=effect_allele_frequency, 
                                  neg_log10_p=-log10(as.numeric(p_value)))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/LeptinAdjBMI/gwas_filtered_oriented_SNPs.txt.gz")

#####################
# Schizophrenia GWAS
#####################
gwas_ss <- fread("data/gwas_summary_stats/Schizophrenia/CLOZUK_PGC2noclo.METAL.assoc.dosage.fix.gz")
gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=BP, EA=A1, OA=A2, beta=log(OR), beta_se=SE, neg_log10_p=-log10(P))] # No EAF column
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=40675, total_controls=64643)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/Schizophrenia/gwas_filtered_oriented_SNPs.txt.gz")

########################
# Smoking Status GWAS
########################
gwas_ss <- fread("data/gwas_summary_stats/SmokingStatus/All_2018_EVERSMK_BBJ_autosome_Pcorrected.txt.gz")
gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=POS, EA=A1, OA=A2, beta=bolt_lmm_fix(BETA, 83830, 81626), beta_se=bolt_lmm_fix(SE, 83830, 81626),
                      EAF=A1Frq, neg_log10_p=-log10(P))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=83830, total_controls=81626)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/SmokingStatus/gwas_filtered_oriented_SNPs.txt.gz")

###############################
# Smoking Status (Males) GWAS
###############################
gwas_ss <- fread("data/gwas_summary_stats/SmokingStatusMale/Male_2018_EVERSMK_BBJ_autosome_Pcorrected.txt.gz")
gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=POS, EA=A1, OA=A2, beta=bolt_lmm_fix(BETA, 67773, 21905), beta_se=bolt_lmm_fix(SE, 67773, 21905),
                      EAF=A1Frq, neg_log10_p=-log10(P))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=67773, total_controls=21905)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/SmokingStatusMale/gwas_filtered_oriented_SNPs.txt.gz")

###############################
# Smoking Status (Females) GWAS
###############################
gwas_ss <- fread("data/gwas_summary_stats/SmokingStatusFemale/Female_2018_EVERSMK_BBJ_autosome_Pcorrected.txt.gz")
gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=POS, EA=A1, OA=A2, beta=bolt_lmm_fix(BETA, 16077, 59721), beta_se=bolt_lmm_fix(SE, 16077, 59721),
                      EAF=A1Frq, neg_log10_p=-log10(P))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=16077, total_controls=59721)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/SmokingStatusFemale/gwas_filtered_oriented_SNPs.txt.gz")

###################################################
# Smoking Status - Cigarettes smoked per day GWAS
###################################################
gwas_ss <- fread("data/gwas_summary_stats/SmokingCigsPerDay/All_2018_CPD_BBJ_autosome_Pcorrected.txt.gz")
gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=POS, EA=A1, OA=A2, beta=BETA, beta_se=SE, EAF=A1Frq, neg_log10_p=-log10(P))]
gwas_ss <- filter_sumstats(gwas_ss, type="continuous", total_samples=72655)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/SmokingCigsPerDay/gwas_filtered_oriented_SNPs.txt.gz")

########################
# Stroke (Any) GWAS
########################
gwas_ss <- fread("data/gwas_summary_stats/StrokeAS/MEGASTROKE.1.AS.EUR.out.gz")
gwas_ss <- gwas_ss[,.(rsid=MarkerName, EA=toupper(Allele1), OA=toupper(Allele2), 
                      beta=Effect, beta_se=StdErr, EAF=Freq1, 
                      neg_log10_p=-log10(`P-value`))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=40585, total_controls=406111)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/StrokeAS/gwas_filtered_oriented_SNPs.txt.gz")

###########################
# Stroke (Ischaemic) GWAS
###########################
gwas_ss <- fread("data/gwas_summary_stats/StrokeIS/MEGASTROKE.2.AIS.EUR.out.gz")
gwas_ss <- gwas_ss[,.(rsid=MarkerName, EA=toupper(Allele1), OA=toupper(Allele2), 
                      beta=Effect, beta_se=StdErr, EAF=Freq1, 
                      neg_log10_p=-log10(`P-value`))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=34217, total_controls=406111)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/StrokeIS/gwas_filtered_oriented_SNPs.txt.gz")

##############################
# Stroke (Large Artery) GWAS
##############################
gwas_ss <- fread("data/gwas_summary_stats/StrokeLAS/MEGASTROKE.3.LAS.EUR.out.gz")
gwas_ss <- gwas_ss[,.(rsid=MarkerName, EA=toupper(Allele1), OA=toupper(Allele2), 
                      beta=Effect, beta_se=StdErr, EAF=Freq1, 
                      neg_log10_p=-log10(`P-value`))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=4373, total_controls=406111)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/StrokeLAS/gwas_filtered_oriented_SNPs.txt.gz")

##############################
# Stroke (Cardioembolic) GWAS
##############################
gwas_ss <- fread("data/gwas_summary_stats/StrokeCES/MEGASTROKE.4.CES.EUR.out.gz")
gwas_ss <- gwas_ss[,.(rsid=MarkerName, EA=toupper(Allele1), OA=toupper(Allele2), 
                      beta=Effect, beta_se=StdErr, EAF=Freq1, 
                      neg_log10_p=-log10(`P-value`))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=7193, total_controls=406111)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/StrokeCES/gwas_filtered_oriented_SNPs.txt.gz")

#############################
# Stroke (Small Vessel) GWAS
#############################
gwas_ss <- fread("data/gwas_summary_stats/StrokeSVS/MEGASTROKE.5.SVS.EUR.out.gz")
gwas_ss <- gwas_ss[,.(rsid=MarkerName, EA=toupper(Allele1), OA=toupper(Allele2), 
                      beta=Effect, beta_se=StdErr, EAF=Freq1, 
                      neg_log10_p=-log10(`P-value`))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=5386, total_controls=406111)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/StrokeSVS/gwas_filtered_oriented_SNPs.txt.gz")

################################
# T2D 2017 GWAS 
################################
gwas_ss <- fread("data/gwas_summary_stats/T2D_2017/METAANALYSIS_DIAGRAM_SE1.txt")
gwas_ss[, c("chr", "pos_b37") := tstrsplit(`Chr:Position`, ":", type.convert=TRUE)]
gwas_ss <- gwas_ss[,.(chr, pos_b37, EA=Allele1, OA=Allele2, beta=Effect, beta_se=StdErr,
                      samples=TotalSampleSize, neg_log10_p=-log10(`P-value`))] # No EAF column
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=26676, total_controls=132532, total_samples=159208)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/T2D_2017/gwas_filtered_oriented_SNPs.txt.gz")

######################
# T2D 2018 GWAS
######################
gwas_ss <- fread("data/gwas_summary_stats/T2D_2018_no_UKB/Mahajan.NatGenet2018b.T2D-noUKBB.European.gz")
gwas_ss <- gwas_ss[,.(chr=Chr, pos_b37=Pos, EA, OA=NEA, beta=Beta, beta_se=SE, EAF, neg_log10_p=-log10(Pvalue))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=55005, total_controls=400308)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/T2D_2018_no_UKB/gwas_filtered_oriented_SNPs.txt.gz")

##################
# T2D Africa GWAS
##################
gwas_ss <- fread("data/gwas_summary_stats/T2D_Africa/ChenJ_31049640.gz")
gwas_ss[, c("chr", "pos_b37") := tstrsplit(MarkerName, ":", keep=1:2, type.convert=TRUE)]
gwas_ss <- gwas_ss[,.(chr, pos_b37, EA=toupper(Allele1), OA=toupper(Allele2), beta=Effect, beta_se=StdErr, n_eff=Weight,
                      EAF=Freq1, neg_log10_p=-log10(`P-value`))]
gwas_ss <- filter_sumstats(gwas_ss)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/T2D_Africa/gwas_filtered_oriented_SNPs.txt.gz")

#####################
# T2D East Asia GWAS
#####################
gwas_ss <- fread("data/gwas_summary_stats/T2D_EastAsia/SpracklenCN_prePMID_T2D_ALL_Primary.txt.gz")
gwas_ss <- gwas_ss[,.(chr=Chr, pos_b37=Pos, EA=toupper(EA), OA=toupper(NEA), beta=Beta, beta_se=SE, n_eff=Neff,
                      EAF, neg_log10_p=-log10(P))]
gwas_ss <- filter_sumstats(gwas_ss)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/T2D_EastAsia/gwas_filtered_oriented_SNPs.txt.gz")

#######################
# T2D Exome GWAS
#######################
gwas_ss <- fread("data/gwas_summary_stats/T2D_Exome/29632382-GCST007515-EFO_0001360-build37.f.tsv.gz")
gwas_ss <- gwas_ss[,.(chr=chromosome, pos_b37=base_pair_location, EA=effect_allele, OA=other_allele, beta, beta_se=standard_error, n_eff=Neff,
                      EAF=effect_allele_frequency, neg_log10_p=-log10(p_value))]
gwas_ss <- filter_sumstats(gwas_ss)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/T2D_Exome/gwas_filtered_oriented_SNPs.txt.gz")

####################################
# T2D Exome GWAS adjusting for BMI
####################################
gwas_ss <- fread("data/gwas_summary_stats/T2D_Exome_Adj_BMI/29632382-GCST007516-EFO_0001360-build37.f.tsv.gz")
gwas_ss <- gwas_ss[,.(chr=chromosome, pos_b37=base_pair_location, EA=effect_allele, OA=other_allele, beta, beta_se=standard_error, n_eff=Neff,
                      EAF=effect_allele_frequency, neg_log10_p=-log10(p_value))]
gwas_ss <- filter_sumstats(gwas_ss)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/T2D_Exome_Adj_BMI/gwas_filtered_oriented_SNPs.txt.gz")

##############################
# T2D FinnGen R6 GWAS
##############################
gwas_ss <- fread("data/gwas_summary_stats/T2D_FinnGen_v6/summary_stats_finngen_R6_E4_DM2.gz")
gwas_ss <- gwas_ss[,.(chr=`#chrom`, pos_b38=pos, EA=alt, OA=ref, beta, beta_se=sebeta, 
                      cases=n_hom_cases+n_het_cases, controls=n_hom_controls+n_het_controls,
                      neg_log10_p=mlogp, EAF=af_alt)] 
gwas_ss <- filter_sumstats(gwas_ss, type="case/control")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/T2D_FinnGen_v6/gwas_filtered_oriented_SNPs.txt.gz")

#################
# T2D GERA GWAS
#################
gwas_ss <- fread("data/gwas_summary_stats/T2D_GERA/Results_DIA2_GERA_chr_1_to_22.txt.gz")
gwas_ss <- gwas_ss[,.(chr, pos_b37=position, EA=alleleB, OA=alleleA, beta=frequentist_add_beta_1, beta_se=frequentist_add_se_1,
                      neg_log10_p=-log10(frequentist_add_pvalue), EAF=freq_effect_allele)]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=6967, total_controls=49670)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/T2D_GERA/gwas_filtered_oriented_SNPs.txt.gz")

#################
# T2D Japan GWAS
#################
gwas_ss <- fread("data/gwas_summary_stats/T2D_Japan/BBJ_BetaBased1.MAF_001.AtLeast2studies.AllChr.txt.gz")
gwas_ss <- gwas_ss[,.(chr=CHR, pos_b37=POS, EA=ALT, OA=REF, beta=BETA, beta_se=SE, samples=N,
                      EAF=Frq, neg_log10_p=-log10(P))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control", total_cases=6967, total_controls=49670)
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/T2D_Japan/gwas_filtered_oriented_SNPs.txt.gz")

################
# T2D PAGE GWAS
################
gwas_ss <- fread("data/gwas_summary_stats/T2D_PAGE/31217584-GCST008048-EFO_0001360-build37.f.tsv.gz")
gwas_ss <- gwas_ss[,.(chr=chromosome, pos_b37=base_pair_location, EA=effect_allele, OA=other_allele, 
                      beta, beta_se=standard_error, cases=n_cas, controls=n-n_cas,
                      EAF=effect_allele_frequency, neg_log10_p=-log10(p_value))]
gwas_ss <- filter_sumstats(gwas_ss, type="case/control")
fwrite(gwas_ss, sep="\t", quote=FALSE, compress="gzip", file="data/filtered_sumstats/T2D_PAGE/gwas_filtered_oriented_SNPs.txt.gz")

