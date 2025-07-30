library(data.table)
library(foreach)
source("src/functions/flip_strand.R")

# Set up parallelisation if on compute node
source("src/functions/par_setup.R")

# Create output directory if needed
system("mkdir -p data/filtered_sumstats", wait=TRUE)

# Make temporary working directory if needed
system("mkdir -p tmp", wait=TRUE)

# Load in HapMap3+ variant set recommended by LDpred2
# https://privefl.github.io/bigsnpr/articles/LDpred2.html
vars_hm3plus <- readRDS("data/LDpred2/map_hm3_plus.rds")
setDT(vars_hm3plus)

# Load in 1000 Genomes data so we can get allele frequencies for reference purposes
vars_1kg <- foreach(this_chr = 1:22, .combine=rbind) %do% {
  fread(sprintf("data/1000G/plink_format/pgen/chr%s.pvar", this_chr), skip="#CHROM")
}

# Filter to variants in HapMap3+ variant set
vars_1kg <- vars_1kg[vars_hm3plus[,.(chr, pos)], on = .(`#CHROM`=chr, POS=pos), nomatch=0]

# Split out info column
columns <- strsplit(vars_1kg$INFO, ";")
kvs <- lapply(columns, strsplit, "=")
cols <- lapply(kvs, function(row) { sapply(row, function(kv) { structure(names=kv[1], kv[2]) }) })
cols <- rbindlist(lapply(cols, as.list), fill=TRUE)
vars_1kg <- cbind(vars_1kg, cols)
vars_1kg[, INFO := NULL]

# Split out multi-allelic sites in 1KG to multiple rows for allele matching to HM3+
mult <- vars_1kg[grepl(",", ALT)]
vars_1kg <- vars_1kg[!grepl(",", ALT)]

split <- foreach(rn = mult[,.I], .combine=rbind) %do% {
  this_row <- mult[rn]
  this_row[, .(`#CHROM`, POS, ID=strsplit(ID, ";")[[1]], REF, ALT=strsplit(ALT, ",")[[1]], QUAL, FILTER,
    AC=strsplit(AC, ",")[[1]], AF=strsplit(AF, ",")[[1]], AN, NS, DP, EAS_AF=strsplit(EAS_AF, ",")[[1]],
    AMR_AF=strsplit(AMR_AF, ",")[[1]], AFR_AF=strsplit(AFR_AF, ",")[[1]], EUR_AF=strsplit(EUR_AF, ",")[[1]],
    SAS_AF=strsplit(SAS_AF, ",")[[1]], AA, VT, EX_TARGET, MULTI_ALLELIC, CIEND, CIPOS,
    CS, END, SVTYPE, MEINFO, SVLEN, TSD, MC, IMPRECISE)]
}
split[nchar(REF) == 1 & nchar(ALT) == 1, VT := "SNP"]
split[nchar(ALT) > 1 & nchar(REF) == 1, VT := "INDEL"]
vars_1kg <- rbind(vars_1kg, split)
vars_1kg <- vars_1kg[order(POS)][order(`#CHROM`)]

# Convert EAFs to numeric
vars_1kg[, AF := as.numeric(AF)]
vars_1kg[, EAS_AF := as.numeric(EAS_AF)]
vars_1kg[, AMR_AF := as.numeric(AMR_AF)]
vars_1kg[, AFR_AF := as.numeric(AFR_AF)]
vars_1kg[, EUR_AF := as.numeric(EUR_AF)]
vars_1kg[, SAS_AF := as.numeric(SAS_AF)]

# Save extended 1KG info
fwrite(vars_1kg, "data/filtered_sumstats/candidate_varset_1KG_info.txt")

# Add in rsids and EAFs to vars_hm3plus
# N.b. matching done to strand flipped alleles first so that strand ambiguous matches (if present)
# assume orientation in HM3+ is the same as 1000G
vars_hm3plus[, c("a1", "a0") := .(a1=flip_strand(a1), a0=flip_strand(a0))]
vars_hm3plus[vars_1kg, on = .(chr=`#CHROM`, pos=POS, a1=ALT, a0=REF), 
  c("rsid_1000G", "EAF_1000G", "EAF_AMR_1000G", "EAF_AFR_1000G", "EAF_EAS_1000G", "EAF_EUR_1000G", "EAF_SAS_1000G") := 
  .(ID, AF, AMR_AF, AFR_AF, EAS_AF, EUR_AF, SAS_AF)]
vars_hm3plus[vars_1kg, on = .(chr=`#CHROM`, pos=POS, a1=REF, a0=ALT), 
  c("rsid_1000G", "EAF_1000G", "EAF_AMR_1000G", "EAF_AFR_1000G", "EAF_EAS_1000G", "EAF_EUR_1000G", "EAF_SAS_1000G") := 
  .(ID, 1-AF, 1-AMR_AF, 1-AFR_AF, 1-EAS_AF, 1-EUR_AF, 1-SAS_AF)]

vars_hm3plus[, c("a1", "a0") := .(a1=flip_strand(a1), a0=flip_strand(a0))]
vars_hm3plus[vars_1kg, on = .(chr=`#CHROM`, pos=POS, a1=ALT, a0=REF), 
  c("rsid_1000G", "EAF_1000G", "EAF_AMR_1000G", "EAF_AFR_1000G", "EAF_EAS_1000G", "EAF_EUR_1000G", "EAF_SAS_1000G") := 
  .(ID, AF, AMR_AF, AFR_AF, EAS_AF, EUR_AF, SAS_AF)]
vars_hm3plus[vars_1kg, on = .(chr=`#CHROM`, pos=POS, a1=REF, a0=ALT), 
  c("rsid_1000G", "EAF_1000G", "EAF_AMR_1000G", "EAF_AFR_1000G", "EAF_EAS_1000G", "EAF_EUR_1000G", "EAF_SAS_1000G") := 
  .(ID, 1-AF, 1-AMR_AF, 1-AFR_AF, 1-EAS_AF, 1-EUR_AF, 1-SAS_AF)]

# Reorganise and rename columns for compatability with downstream scripts
vars_hm3plus <- vars_hm3plus[, .(chr, pos_b36=pos_hg18, pos_b37=pos, pos_b38=pos_hg38, rsid_HapMap3=rsid, rsid_1000G,
                                 effect_allele=a1, other_allele=a0, EAF_UKBB=af_UKBB, EAF_1000G, EAF_AFR_1000G, 
                                 EAF_AMR_1000G, EAF_EAS_1000G, EAF_EUR_1000G, EAF_SAS_1000G, ld, block_id)]

# Write out:
fwrite(vars_hm3plus[order(pos_b38)][order(chr)], sep="\t", quote=FALSE, file="data/filtered_sumstats/candidate_varset.txt")

