library(data.table)
library(foreach)
source("src/functions/flip_strand.R")

# Set up parallelisation if on compute node
source("src/functions/par_setup.R")

# Create output directory if needed
system("mkdir -p data/filtered_sumstats", wait=TRUE)

# Make temporary working directory if needed
system("mkdir -p tmp", wait=TRUE)

# Load in HapMap3 variant set - LDpred2 recommends using this variant set
vars_hapmap <- lapply(list.files("data/HapMap3/", pattern="*.map$", full.names=TRUE), fread)
names(vars_hapmap) <- list.files("data/HapMap3/", pattern="*.map$")
names(vars_hapmap) <- gsub(".qc.poly.recode.map", "", names(vars_hapmap))
names(vars_hapmap) <- gsub("hapmap3_r1_b36_fwd.", "", names(vars_hapmap))
vars_hapmap <- rbindlist(vars_hapmap, idcol="pop")
setnames(vars_hapmap, c("V1", "V2", "V3", "V4"), c("chr", "rsid", "cm", "pos"))
vars_hapmap[, cm := TRUE]
vars_hapmap <- dcast(vars_hapmap, chr + pos + rsid ~ pop, value.var="cm", fill=FALSE)
vars_hapmap <- vars_hapmap[chr %in% 1:22]
fwrite(vars_hapmap, sep="\t", quote=FALSE, "tmp/hapmap3_autosomal_variants_b36.txt")

# Lift over to hg19/GRCh37, in which UKB is imputed to and stored
vars_hapmap[, liftOver_pos := sprintf("chr%s:%s-%s", chr, pos, pos)]
fwrite(vars_hapmap[,.(liftOver_pos)], col.names=FALSE, quote=FALSE, file="tmp/hapmap3_autosome_b36.pos")

cmd <- "liftOver -positions"
cmd <- paste(cmd, "tmp/hapmap3_autosome_b36.pos")
cmd <- paste(cmd, "data/liftOver/hg18ToHg19.over.chain.gz")
cmd <- paste(cmd, "tmp/hapmap3_autosome_b37.pos")
cmd <- paste(cmd, "tmp/hapmap3_autosome_b36_b37_unmapped.txt")
system(cmd, wait=TRUE)

unmapped <- fread("tmp/hapmap3_autosome_b36_b37_unmapped.txt", header=FALSE)
unmapped <- cbind(unmapped[seq(2, .N, by=2)], unmapped[seq(1, .N, by=2)])
setnames(unmapped, c("liftOver_pos", "reason")) # All '#Deleted in new'
vars_hapmap <- vars_hapmap[!unmapped, on = .(liftOver_pos)]

b37 <- fread("tmp/hapmap3_autosome_b37.pos", header=FALSE)
b37[, c("chr_b37", "pos_b37", "pos_b37.2") := tstrsplit(V1, ":|-")]
b37 <- b37[,.(chr_b37, pos_b37)]
b37[, chr_b37 := as.integer(gsub("chr", "", chr_b37))]
b37[, pos_b37 := as.integer(pos_b37)]
vars_hapmap <- cbind(b37, vars_hapmap)

# About 1.55 million variants
vars_hapmap <- vars_hapmap[, .(chr, rsid_b36=rsid, pos_b36=pos, pos_b37, 
                               ASW, CEU, CHB, CHD, GIH, JPT, LWK, MEX, MKK, TSI, YRI)]

# Some of the GWAS we are building GRS from are exclusively EXOME GWAS, so it also makes 
# sense to try and include common exome variants even if those were not part of HapMap3.
ex1 <- fread("data/gwas_summary_stats/Leptin/GCST90007310_buildGRCh37.tsv.gz")
ex2 <- fread("data/gwas_summary_stats/T2D_Exome/29632382-GCST007515-EFO_0001360-build37.f.tsv.gz")

exome <- unique(rbind(
  ex1[, .(chr=chromosome, pos=base_pair_location)],
  ex2[chromosome %in% 1:22, .(chr=chromosome, pos=base_pair_location)]
))

# Remove any variants already in HapMap3 set (leaves around 208K variants)
exome <- exome[!vars_hapmap, on = .(chr, pos=pos_b37)]

# Next, load in the 1000 genomes data, so we can get the alleles at each position.
# Importantly, we want to orient all GWAS to have the same effect allele so we can
# later derive a metaGRS. 
vars_1kg <- foreach(this_chr = 1:22, .combine=rbind) %do% {
  fread(sprintf("data/1000G/plink_format/pgen/chr%s.pvar", this_chr), skip="#CHROM")
}
vars_1kg <- rbind(
  vars_1kg[vars_hapmap[,.(chr, pos_b37)], on = .(`#CHROM`=chr, POS=pos_b37), nomatch=0],
  vars_1kg[exome, on = .(`#CHROM`=chr, POS=pos), nomatch=0]
)
gc()

# Remove any duplicate chromosome positions (N=168)
mult <- vars_1kg[,.N,by=.(`#CHROM`, POS)][N > 1]
vars_1kg <- vars_1kg[!mult, on = .(`#CHROM`, POS)]

# Filter to bi-allelic SNPs (N=11288 - note most of these are multi-allelic ALT, only 219 are indels)
vars_1kg <- vars_1kg[nchar(REF) == 1 & nchar(ALT) == 1]

# Split out info column
vars_1kg[, c("AC", "AF", "AN", "NS", "DP", "EAS_AF", "AMR_AF", "AFR_AF", "EUR_AF", "SAS_AF", "AA", "VT", "EX_TARGET") := tstrsplit(INFO, ";")]
vars_1kg[, INFO := NULL]

# Extract alelle frequency information
vars_1kg[, AF := as.numeric(gsub("AF=", "", AF))]
vars_1kg[, EAS_AF := as.numeric(gsub("EAS_AF=", "", EAS_AF))]
vars_1kg[, AMR_AF := as.numeric(gsub("AMR_AF=", "", AMR_AF))]
vars_1kg[, AFR_AF := as.numeric(gsub("AFR_AF=", "", AFR_AF))]
vars_1kg[, EUR_AF := as.numeric(gsub("EUR_AF=", "", EUR_AF))]
vars_1kg[, SAS_AF := as.numeric(gsub("SAS_AF=", "", SAS_AF))]

# Filter to variants present in UKB data (N=96 in 1KG not present in UKB by chr and pos)
ukb_vars <- foreach(this_chr = 1:22, .combine=rbind) %do% {
  fread(sprintf("data/UKB/genetics/imputed_bed/ukb_imp_v3_dedup_chr%s.bim", this_chr), header=FALSE)
}
setnames(ukb_vars, c("chr", "rsid", "cM", "pos", "A1", "A2"))
vars_1kg <- vars_1kg[unique(ukb_vars[,.(chr, pos)]), on = .(`#CHROM`=chr, POS=pos), nomatch=0]

# Get info on allele frequency in UKB
for (this_chr in 1:22) {
  ukb_maf <- fread(sprintf("data/UKB/genetics/reference_files/ukb_impv3_chr%s_snpstats.txt", this_chr), skip="alternate_ids")
  ukb_vars[ukb_maf, on = .(chr=chromosome, pos=position, A1=minor_allele, A2=major_allele),
           c("A1freq", "INFO") := .(minor_allele_frequency, impute_info)]
  ukb_vars[ukb_maf, on = .(chr=chromosome, pos=position, A2=minor_allele, A1=major_allele),
           c("A1freq", "INFO") := .(1-minor_allele_frequency, impute_info)]
}
# Note warnings arising above are fread discarding footer line when reading files 

# Make sure we can also match by allele to UK Biobank
vars_1kg[, allele_match := FALSE]
vars_1kg[ukb_vars, on = .(`#CHROM`=chr, POS=pos, ALT=A1, REF=A2), c("allele_match", "oriented", "flipped") := .(TRUE, TRUE, FALSE)]
vars_1kg[ukb_vars, on = .(`#CHROM`=chr, POS=pos, ALT=A2, REF=A1), c("allele_match", "oriented", "flipped") := .(TRUE, FALSE, FALSE)]

# Some SNPs may be on the opposite strand in UK Biobank, which we detect and fix here:
vars_1kg[!(allele_match), c("ALT", "REF", "flipped") := .(flip_strand(ALT), flip_strand(REF), TRUE)]
vars_1kg[ukb_vars, on = .(`#CHROM`=chr, POS=pos, ALT=A1, REF=A2), c("allele_match", "oriented") := .(TRUE, TRUE)]
vars_1kg[ukb_vars, on = .(`#CHROM`=chr, POS=pos, ALT=A2, REF=A1), c("allele_match", "oriented") := .(TRUE, FALSE)]

# Remove any variants which matched on chromosome and position 
# but not on alleles even after checking for strand mismatch 
# (N=2 SNPs)
vars_1kg <- vars_1kg[(allele_match)]

# For strand ambiguous alleles (A/T or G/C SNPs) for now 
# we will assume same strand orientation as (majority) of 1000g
pct_flipped <- vars_1kg[ALT != flip_strand(REF), sum(flipped)/.N]
if (pct_flipped < 0.5) {
  # pct_flipped should be 0, so code executed is here:
	vars_1kg[(flipped) & ALT == flip_strand(REF), c("ALT", "REF", "oriented", "flipped") := .(flip_strand(ALT), flip_strand(REF), !oriented, FALSE)]
} else {
	vars_1kg[!(flipped) & ALT == flip_strand(REF), c("ALT", "REF", "oriented", "flipped") := .(flip_strand(ALT), flip_strand(REF), !oriented, TRUE)]
}

# Orient allele we want to use as the effect allele (ALT) to first occurring allele 
# in the .bim format UKB files (A1) - LDpred loads this allele (A1) as the dosage, so
# we ultimately want to align gwas summary stat effect alleles to this allele anyway.
vars_1kg[!(oriented), 
         c("ALT", "REF", "AF", "EAS_AF", "AMR_AF", "AFR_AF", "EUR_AF", "SAS_AF", "oriented") := 
         .(REF, ALT, 1-AF, 1-EAS_AF, 1-AMR_AF, 1-AFR_AF, 1-EUR_AF, 1-SAS_AF, TRUE)]

# Add in UKB allele frequency and INFO score
vars_1kg[ukb_vars, on = .(`#CHROM`=chr, POS=pos, ALT=A1, REF=A2), c("UKB_AF", "UKB_INFO") := .(A1freq, INFO)]

# Remove palindromic SNPs which will be difficult to match between datasets, even with known allele-frequencies 
maf_cutoff <- 0.42
palindromic_snps <- vars_1kg[REF == flip_strand(ALT)] # N = 144,205; 8.6% SNPs
bad <- palindromic_snps[(
  (UKB_AF > maf_cutoff & UKB_AF < 1 - maf_cutoff) |  # 15,148 SNPs not matchable to UKB
  (EUR_AF > maf_cutoff & EUR_AF < 1 - maf_cutoff)    # 15,157 SNPs not matchable to 1000G EUR population
  # 16,955 SNPs total; 11.8% of palindromic SNPs, 1% overall, with MAF between 42% and 50%
  # in UKB, 1000G, or any 1000G European population group. 
)]

## bad <- palindromic_snps[(
##   (UKB_AF > maf_cutoff & UKB_AF < 1 - maf_cutoff) |  # 15,148 SNPs not matchable to UKB
##   (AF > maf_cutoff & AF < 1 - maf_cutoff)         |  # 15,571 SNPs not matchable to 1000G (all ancestries combined)
##   (AFR_AF > maf_cutoff & AFR_AF < 1 - maf_cutoff) |  # 14,037 SNPs not matchable to 1000G AFR population
##   (AMR_AF > maf_cutoff & AMR_AF < 1 - maf_cutoff) |  # 15,270 SNPs not matchable to 1000G AMR population
##   (EAS_AF > maf_cutoff & EAS_AF < 1 - maf_cutoff) |  # 14,244 SNPs not matchable to 1000G EAS population
##   (EUR_AF > maf_cutoff & EUR_AF < 1 - maf_cutoff) |  # 15,157 SNPs not matchable to 1000G EUR population
##   (SAS_AF > maf_cutoff & SAS_AF < 1 - maf_cutoff)    # 15,267 SNPs not matchable to 1000G SAS population
##   # 44,740 SNPs total; 31% of palindromic SNPs, 2.7% overall, with MAF between 42% and 50%
##   # in UKB, 1000G, or any 1000G population group. 
## )]

vars_1kg <- vars_1kg[!bad, on = .(`#CHROM`, POS)]

# Drop 5 palindromic SNPs that appear to have swapped strand orientation sometime between 1000G and UKB
# (despite all other candidate SNPs being oriented on the same strand)
vars_1kg <- vars_1kg[REF != flip_strand(ALT) | (EUR_AF < 0.5 & UKB_AF < 0.5) | (EUR_AF > 0.5 & UKB_AF > 0.5)]

## # Uncomment this if you checking against all 1000G populations instead of just Europeans above
## vars_1kg <- vars_1kg[REF != flip_strand(ALT) | (
##        ((AF < 0.5 & UKB_AF < 0.5) | (AF > 0.5 & UKB_AF > 0.5)) &
##        ((AFR_AF < 0.5 & UKB_AF < 0.5) | (AFR_AF > 0.5 & UKB_AF > 0.5)) &
##        ((AMR_AF < 0.5 & UKB_AF < 0.5) | (AMR_AF > 0.5 & UKB_AF > 0.5)) &
##        ((EAS_AF < 0.5 & UKB_AF < 0.5) | (EAS_AF > 0.5 & UKB_AF > 0.5)) &
##        ((EUR_AF < 0.5 & UKB_AF < 0.5) | (EUR_AF > 0.5 & UKB_AF > 0.5)) &
##        ((SAS_AF < 0.5 & UKB_AF < 0.5) | (SAS_AF > 0.5 & UKB_AF > 0.5))
## )]

# Get GRCh38 positions in case we need to liftOver any GWAS
vars_1kg[, liftOver_pos := sprintf("chr%s:%s-%s", `#CHROM`, POS, POS)]
fwrite(vars_1kg[,.(liftOver_pos)], col.names=FALSE, quote=FALSE, file="tmp/1000G_autosome_b37.pos")

cmd <- "liftOver -positions"
cmd <- paste(cmd, "tmp/1000G_autosome_b37.pos")
cmd <- paste(cmd, "data/liftOver/hg19ToHg38.over.chain.gz")
cmd <- paste(cmd, "tmp/1000G_autosome_b38.pos")
cmd <- paste(cmd, "tmp/1000G_autosome_b37_b38_unmapped.txt")
system(cmd, wait=TRUE)

unmapped <- fread("tmp/1000G_autosome_b37_b38_unmapped.txt", header=FALSE)
unmapped <- cbind(unmapped[seq(2, .N, by=2)], unmapped[seq(1, .N, by=2)])
setnames(unmapped, c("liftOver_pos", "reason")) # All '#Deleted in new'
vars_1kg <- vars_1kg[!unmapped, on = .(liftOver_pos)] # Drop variants deleted in hg38

b38 <- fread("tmp/1000G_autosome_b38.pos", header=FALSE)
b38[, c("chr_b38", "pos_b38", "pos_b38.2") := tstrsplit(V1, ":|-")]
b38 <- b38[,.(chr_b38, pos_b38)]
b38[, chr_b38 := as.integer(gsub("chr", "", chr_b38))] # note some now on alternate contigs, becoming 'NA' here (with warning), and removed later
b38[, pos_b38 := as.integer(pos_b38)]
vars_1kg <- cbind(b38, vars_1kg)
vars_1kg <- vars_1kg[!is.na(chr_b38)] # drop (N=200) SNPs on alternate contigs on b38
vars_1kg <- vars_1kg[chr_b38 == `#CHROM`] # drop 1 SNP moved from chr 19 to chr 7 on b38

# Collate filtered variant information
varset <- merge(vars_1kg, vars_hapmap, by.x=c("#CHROM", "POS"), by.y=c("chr", "pos_b37"), all.x=TRUE)
varset <- varset[, .(chr=`#CHROM`, pos_b36, pos_b37=POS, pos_b38, rsid_HapMap3=rsid_b36, rsid_1000G=ID, 
                     effect_allele=ALT, other_allele=REF, INFO_UKB=UKB_INFO, EAF_UKB=UKB_AF, 
                     EAF_1000G=AF, EUR_EAF_1000G=EUR_AF, AFR_EAF_1000G=EUR_AF, AMR_EAF_1000G=AMR_AF,
                     EAS_EAF_1000G=EAS_AF, SAS_EAF_1000G=SAS_AF, ASW_HapMap3=as.logical(ASW), 
                     CEU_HapMap3=as.logical(CEU), CHB_HapMap3=as.logical(CHB), CHD_HapMap3=as.logical(CHD),
                     GIH_HapMap3=as.logical(GIH), JPT_HapMap3=as.logical(JPT),
                     LWK_HapMap3=as.logical(LWK), MEX_HapMap3=as.logical(MEX), MKK_HapMap3=as.logical(MKK),
                     TSI_HapMap3=as.logical(TSI), YRI_HapMap3=as.logical(YRI))]
varset[is.na(pos_b36), c("ASW_HapMap3", "CEU_HapMap3", "CHB_HapMap3", "CHD_HapMap3",
	"GIH_HapMap3", "JPT_HapMap3", "LWK_HapMap3", "MEX_HapMap3",
	"MKK_HapMap3", "TSI_HapMap3", "YRI_HapMap3") := FALSE]

# For exome variants not in HapMap3, also get positions on build36 in case they are needed
varset[is.na(pos_b36), liftOver_pos := sprintf("chr%s:%s-%s", chr, pos_b37, pos_b37)]
fwrite(varset[,.(liftOver_pos)], col.names=FALSE, quote=FALSE, file="tmp/1000G_autosome_b37.pos")

cmd <- "liftOver -positions"
cmd <- paste(cmd, "tmp/1000G_autosome_b37.pos")
cmd <- paste(cmd, "data/liftOver/hg19ToHg18.over.chain.gz")
cmd <- paste(cmd, "tmp/1000G_autosome_b36.pos")
cmd <- paste(cmd, "tmp/1000G_autosome_b37_b36_unmapped.txt")
system(cmd, wait=TRUE)

unmapped <- fread("tmp/1000G_autosome_b37_b36_unmapped.txt", header=FALSE)
unmapped <- cbind(unmapped[seq(2, .N, by=2)], unmapped[seq(1, .N, by=2)])
setnames(unmapped, c("liftOver_pos", "reason")) # All '#Deleted in new'
varset <- varset[unmapped, on = .(liftOver_pos), liftOver_pos := NA] # Leave in variants that were not in b36, but added in b37

b36 <- fread("tmp/1000G_autosome_b36.pos", header=FALSE)
b36[, c("chr_b36", "pos_b36", "pos_b36.2") := tstrsplit(V1, ":|-")]
b36 <- b36[,.(chr_b36, pos_b36)]
b36[, chr_b36 := as.integer(gsub("chr", "", chr_b36))] # note some now on alternate contigs, becoming 'NA' here (with warning), and removed later
b36[, pos_b36 := as.integer(pos_b36)]

varset[!is.na(liftOver_pos), pos_b36 := ifelse(chr == b36$chr, b36$pos_b36, NA)] # exclude 1 variant that moved chromosomes
varset[, liftOver_pos := NULL]

# Load in UK Biobank SNP QC information
ukb_snp_qc <- fread("data/UKB/genetics/reference_files/ukb_snp_qc.txt")

varset[, UKB_directly_genotyped := FALSE]
varset[ukb_snp_qc, on = .(chr=chromosome, pos_b37=position), UKB_directly_genotyped := TRUE]
varset[ukb_snp_qc, on = .(chr=chromosome, pos_b37=position), UKB_used_for_HetMiss := as.logical(in_HetMiss)]
varset[ukb_snp_qc, on = .(chr=chromosome, pos_b37=position), UKB_used_for_PCA := as.logical(in_PCA)]
varset[ukb_snp_qc, on = .(chr=chromosome, pos_b37=position), UKB_used_for_phasing := as.logical(in_Phasing_Input)]
varset[ukb_snp_qc, on = .(chr=chromosome, pos_b37=position), UKB_used_for_kinship := as.logical(in_Relatedness)]
varset[ukb_snp_qc, on = .(chr=chromosome, pos_b37=position), UKB_BiLEVE_chip := ifelse(array == 1, FALSE, TRUE)]
varset[ukb_snp_qc, on = .(chr=chromosome, pos_b37=position), UKB_Axiom_chip := ifelse(array == 0, FALSE, TRUE)]

varset[ukb_snp_qc, on = .(chr=chromosome, pos_b37=position), UKB_BiLEVE_pct_batch_pass_QC :=
  (UKBiLEVEAX_b1_qc + UKBiLEVEAX_b2_qc + UKBiLEVEAX_b3_qc + UKBiLEVEAX_b4_qc + UKBiLEVEAX_b5_qc + 
   UKBiLEVEAX_b6_qc + UKBiLEVEAX_b7_qc + UKBiLEVEAX_b8_qc + UKBiLEVEAX_b9_qc + UKBiLEVEAX_b10_qc + 
   UKBiLEVEAX_b11_qc)/11*100]
varset[!(UKB_BiLEVE_chip), UKB_BiLEVE_pct_batch_pass_QC := NA]

varset[ukb_snp_qc, on = .(chr=chromosome, pos_b37=position), UKB_Axiom_pct_batch_pass_QC :=
  (Batch_b001_qc + Batch_b002_qc + Batch_b003_qc + Batch_b004_qc + Batch_b005_qc + Batch_b006_qc + 
   Batch_b007_qc + Batch_b008_qc + Batch_b009_qc + Batch_b010_qc + Batch_b011_qc + Batch_b012_qc + 
   Batch_b013_qc + Batch_b014_qc + Batch_b015_qc + Batch_b016_qc + Batch_b017_qc + Batch_b018_qc + 
   Batch_b019_qc + Batch_b020_qc + Batch_b021_qc + Batch_b022_qc + Batch_b023_qc + Batch_b024_qc + 
   Batch_b025_qc + Batch_b026_qc + Batch_b027_qc + Batch_b028_qc + Batch_b029_qc + Batch_b030_qc + 
   Batch_b031_qc + Batch_b032_qc + Batch_b033_qc + Batch_b034_qc + Batch_b035_qc + Batch_b036_qc + 
   Batch_b037_qc + Batch_b038_qc + Batch_b039_qc + Batch_b040_qc + Batch_b041_qc + Batch_b042_qc + 
   Batch_b043_qc + Batch_b044_qc + Batch_b045_qc + Batch_b046_qc + Batch_b047_qc + Batch_b048_qc + 
   Batch_b049_qc + Batch_b050_qc + Batch_b051_qc + Batch_b052_qc + Batch_b053_qc + Batch_b054_qc + 
   Batch_b055_qc + Batch_b056_qc + Batch_b057_qc + Batch_b058_qc + Batch_b059_qc + Batch_b060_qc + 
   Batch_b061_qc + Batch_b062_qc + Batch_b063_qc + Batch_b064_qc + Batch_b065_qc + Batch_b066_qc + 
   Batch_b067_qc + Batch_b068_qc + Batch_b069_qc + Batch_b070_qc + Batch_b071_qc + Batch_b072_qc + 
   Batch_b073_qc + Batch_b074_qc + Batch_b075_qc + Batch_b076_qc + Batch_b077_qc + Batch_b078_qc + 
   Batch_b079_qc + Batch_b080_qc + Batch_b081_qc + Batch_b082_qc + Batch_b083_qc + Batch_b084_qc + 
   Batch_b085_qc + Batch_b086_qc + Batch_b087_qc + Batch_b088_qc + Batch_b089_qc + Batch_b090_qc + 
   Batch_b091_qc + Batch_b092_qc + Batch_b093_qc + Batch_b094_qc + Batch_b095_qc)/95*100]
varset[!(UKB_Axiom_chip), UKB_Axiom_pct_batch_pass_QC := NA]

# Clean up temporary files
system("rm liftOver_* tmp/1000G_* tmp/hapmap3_*", wait=TRUE)

# Write out:
fwrite(varset[order(pos_b37)][order(chr)], sep="\t", quote=FALSE, file="data/filtered_sumstats/filtered_oriented_SNPs.txt")

