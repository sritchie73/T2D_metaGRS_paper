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

# Lift over to hg19/GRCh37 which most GWASs and datasets (still) use
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

# Get GRCh38 positions, which are used in All of Us, and some newer GWASs
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

dups <- vars_1kg[,.N,by=.(chr_b38, pos_b38)][N > 1]
vars_1kg <- vars_1kg[!dups, on=.(chr_b38, pos_b38)] # drop 2 SNPs that have been merged into one on b38

# Collate filtered variant information
varset <- merge(vars_1kg, vars_hapmap, by.x=c("#CHROM", "POS"), by.y=c("chr", "pos_b37"), all.x=TRUE)
varset <- varset[, .(chr=`#CHROM`, pos_b36, pos_b37=POS, pos_b38, rsid_HapMap3=rsid_b36, rsid_1000G=ID,
                     effect_allele=ALT, other_allele=REF, 
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

# Remove palindromic SNPs which will be difficult to match between datasets, even with known allele-frequencies 
maf_cutoff <- 0.42
palindromic_snps <- varset[effect_allele == flip_strand(other_allele)] # N = 144,133; 8.6% SNPs

bad <- palindromic_snps[(
  (AFR_EAF_1000G > maf_cutoff & AFR_EAF_1000G < 1 - maf_cutoff) |  # 15,147 SNPs not matchable to 1000G AFR population
  (AMR_EAF_1000G > maf_cutoff & AMR_EAF_1000G < 1 - maf_cutoff) |  # 15,258 SNPs not matchable to 1000G AMR population
  (EAS_EAF_1000G > maf_cutoff & EAS_EAF_1000G < 1 - maf_cutoff) |  # 14,234 SNPs not matchable to 1000G EAS population
  (EUR_EAF_1000G > maf_cutoff & EUR_EAF_1000G < 1 - maf_cutoff) |  # 15,147 SNPs not matchable to 1000G EUR population
  (SAS_EAF_1000G > maf_cutoff & SAS_EAF_1000G < 1 - maf_cutoff)    # 15,259 SNPs not matchable to 1000G SAS population
  # 35,689 SNPs total; 25% of palindromic SNPs, 2.1% overall, with MAF between 42% and 50% in any 1000G continental population
)]

varset <- varset[!bad, on = .(chr, pos_b37)]

# Clean up temporary files
system("rm liftOver_* tmp/1000G_* tmp/hapmap3_*", wait=TRUE)

# Write out:
fwrite(varset[order(pos_b38)][order(chr)], sep="\t", quote=FALSE, file="data/filtered_sumstats/candidate_varset.txt")


