library(data.table)
library(foreach)
source("src/functions/flip_strand.R")
source("src/functions/eaf_to_maf.R")

# Load candidate variant set
varset <- fread("data/filtered_sumstats/candidate_varset.txt")

# Load variant information extracted from All of Us Research Platform by Carles
#
# ACAF Threshold subset of All of Us:
# - Variants that are frequent in the AoU computed ancestry subpopulations.
# - The cutoff we use is population-specific allele frequency (AF) > 1% OR population-specific allele count (AC) > 100, in any computed ancestry subpopulations
#
# https://support.researchallofus.org/hc/en-us/articles/14929793660948-Smaller-Callsets-for-Analyzing-Short-Read-WGS-SNP-Indel-Data-with-Hail-MT-VCF-and-PLINK
aou <- foreach(chrIdx = 1:22, .combine=rbind) %do% {
  fread(sprintf("data/All_of_Us/candidate_bim_files/candidate_varset_chr%s.bim", chrIdx))
}
setnames(aou, c("chr", "varid", "cm", "pos", "alt", "ref"))
aou[, chr := as.integer(gsub("chr", "", chr))]

# Filter to variants present in All of Us by chr and position (29,303 variants not present in ACAF subset)
varset <- varset[unique(aou[,.(chr, pos)]), on = .(chr, pos_b38=pos), nomatch=0]

# Make sure we can also match by allele to All of Us, orienting effect alleles to the ALT allele in All of Us
varset[, allele_match := FALSE]
varset[aou, on = .(chr, pos_b38=pos, effect_allele=alt, other_allele=ref), c("allele_match", "oriented", "flipped") := .(TRUE, TRUE, FALSE)]
varset[aou, on = .(chr, pos_b38=pos, effect_allele=ref, other_allele=alt), c("allele_match", "oriented", "flipped") := .(TRUE, FALSE, FALSE)]

# Some SNPs may be on the opposite strand in 1000G compared to All of Us, which we detect and fix here (n=1,237 SNPs)
varset[!(allele_match), c("effect_allele", "other_allele", "flipped") := .(flip_strand(effect_allele), flip_strand(other_allele), TRUE)]
varset[aou, on = .(chr, pos_b38=pos, effect_allele=alt, other_allele=ref), c("allele_match", "oriented") := .(TRUE, TRUE)]
varset[aou, on = .(chr, pos_b38=pos, effect_allele=ref, other_allele=alt), c("allele_match", "oriented") := .(TRUE, FALSE)]

# Remove any variants which matched on chromosome and position 
# but not on alleles even after checking for strand mismatch 
# (N=88 SNPs)
varset <- varset[(allele_match)]

# Orient effect alleles to the ALT allele in All of Us - LDpred2 loads this allele as the
# dosage, so we ultimately want to align gwas summary stat effect alleles to this
varset[!(oriented),
  c("effect_allele", "other_allele",  "EAS_EAF_1000G", "AMR_EAF_1000G", "AFR_EAF_1000G", "EUR_EAF_1000G", "SAS_EAF_1000G", "oriented") :=
  .(other_allele, effect_allele, 1-EAS_EAF_1000G, 1-AMR_EAF_1000G, 1-AFR_EAF_1000G, 1-EUR_EAF_1000G, 1-SAS_EAF_1000G, TRUE)]

# Add in SNP ID from All of Us
varset[aou, on = .(chr, pos_b38=pos, effect_allele=alt, other_allele=ref), AoU_varID := varid]

# Add in allele frequencies from ACAF callset
aou_af <- fread("data/All_of_Us/candidate_varset_AF.txt.gz")
aou_af[,AoU_varID := paste0("chr", gsub("-", ":", vid))]
varset[aou_af, on = .(AoU_varID), 
  c("rsid_AoU", "AoU_varID2", "EUR_EAF_AoU", "AFR_EAF_AoU", "AMR_EAF_AoU", "EAS_EAF_AoU", "SAS_EAF_AoU", "MID_EAF_AoU", "OTH_EAF_AoU") :=
  .(dbsnp_rsid, vid, gvs_eur_af, gvs_afr_af, gvs_amr_af, gvs_eas_af, gvs_sas_af, gvs_mid_af, gvs_oth_af)]

# SNPs with allele counts < 20 have their frequency set to NA for that ancestry in All of Us. To simplify filtering, set these to 0 in our
# varset table
varset[is.na(EUR_EAF_AoU), EUR_EAF_AoU := 0]
varset[is.na(AFR_EAF_AoU), AFR_EAF_AoU := 0]
varset[is.na(AMR_EAF_AoU), AMR_EAF_AoU := 0]
varset[is.na(EAS_EAF_AoU), EAS_EAF_AoU := 0]
varset[is.na(SAS_EAF_AoU), SAS_EAF_AoU := 0]
varset[is.na(MID_EAF_AoU), MID_EAF_AoU := 0]
varset[is.na(OTH_EAF_AoU), OTH_EAF_AoU := 0]
stop()

# Filter to variants that have >1% frequency in any All of Us ancestry - this filter is applied later on by the
# per-ancestry LDpred2 quality control when calculating the LD matrix (assuming 10K samples)
varset <- varset[
  maf(EUR_EAF_AoU) >= 0.01 |   # 1,295,104 SNPs pass, 319,973 below 1% frequency
  maf(AFR_EAF_AoU) >= 0.01 |   # 1,485,124 SNPs pass, 129,953 below 1% frequency
  maf(AMR_EAF_AoU) >= 0.01 |   # 1,406,229 SNPs pass, 208,848 below 1% frequency
  maf(EAS_EAF_AoU) >= 0.01 |   # 1,204,683 SNPs pass, 410,394 below 1% frequency
  maf(SAS_EAF_AoU) >= 0.01 |   # 1,313,545 SNPs pass, 301,532 below 1% frequency
  maf(MID_EAF_AoU) >= 0.01 |   # 1,354,677 SNPs pass, 260,400 below 1% frequency
  maf(OTH_EAF_AoU) >= 0.01     # 1,465,460 SNPs pass, 149,617 below 1% frequency
  # 1,537,448 SNPs pass overall, 77,629 below 1% frequency in all ancestries
]

# Filter any strand ambiguous SNPs that will be difficult to match in All of Us based on MAF
# Reminder - we have already applied similar filtering based on allele frequencies in 1000 Genomes populations
maf_cutoff <- 0.42
varset[, ambig := effect_allele == flip_strand(other_allele)]
varset <- varset[!(ambig) | (
  maf(EUR_EAF_AoU) <= maf_cutoff & # 578 SNPs not matchable in All of Us EUR population
  maf(AFR_EAF_AoU) <= maf_cutoff & # 6,828 SNPs not matchable in All of Us AFR population
  maf(AMR_EAF_AoU) <= maf_cutoff & # 528 SNPs not matchable in All of Us AMR population
  maf(EAS_EAF_AoU) <= maf_cutoff & # 758 SNPs not matchable in All of Us EAS population
  maf(SAS_EAF_AoU) <= maf_cutoff & # 392 SNPs not matchable in All of Us SAS population
  maf(MID_EAF_AoU) <= maf_cutoff & # 2,171 SNPs not matchable in All of Us MID population
  maf(OTH_EAF_AoU) <= maf_cutoff   # 1,013 SNPs not matchable in All of Us OTH population
  # 9,981 SNPs total excluded
)]

# Filter out any strand ambiguous SNPs where the allele frequency is inconsistent (in 
# terms of being above or below 50%) across ancestries, which make strand alignment
# challenging in All of Us (which we don't do on an ancestry-specific basis).
ambig <- varset[(ambig)]
ambig <- melt(ambig, id.vars="AoU_varID", measure.vars=c("EUR_EAF_AoU", "AFR_EAF_AoU", "AMR_EAF_AoU", "EAS_EAF_AoU", "SAS_EAF_AoU", "MID_EAF_AoU", "OTH_EAF_AoU"), variable.name="ancestry", value.name="EAF")
ambig[, ancestry := gsub("_.*", "", ancestry)]
ambig <- ambig[,.(consistent = all(EAF < 0.5) | all(EAF > 0.5)), by=.(AoU_varID)]
bad <- ambig[!(consistent)] # 3,601 SNPs excluded
varset <- varset[!bad, on = .(AoU_varID)]

# Reorganise columns
varset <- varset[,.(chr, pos_b36, pos_b37, pos_b38, rsid_HapMap3, rsid_1000G, rsid_AoU, AoU_varID, AoU_varID2, effect_allele, other_allele, 
  EAF_1000G, EUR_EAF_1000G, AFR_EAF_1000G, AMR_EAF_1000G, EAS_EAF_1000G, SAS_EAF_1000G, 
  EUR_EAF_AoU, AFR_EAF_AoU, AMR_EAF_AoU, EAS_EAF_AoU, SAS_EAF_AoU, MID_EAF_AoU, OTH_EAF_AoU,
  ASW_HapMap3, CEU_HapMap3, CHB_HapMap3, CHD_HapMap3, GIH_HapMap3, JPT_HapMap3, 
  LWK_HapMap3, MEX_HapMap3, MKK_HapMap3, TSI_HapMap3, YRI_HapMap3)]

# Write out
fwrite(varset[order(pos_b38)][order(chr)], sep="\t", quote=FALSE, "data/filtered_sumstats/filtered_oriented_SNPs.txt")

