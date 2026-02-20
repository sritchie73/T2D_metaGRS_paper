library(data.table)
library(foreach)
library(R.utils)

# Pull down pvar and MAF+INFO information from persistent storage to the local
# cloud workstation
for (this_chr in c(1:22, "X")) {
  system(sprintf("dx download 'common/Imputed Genotypes/TopMed (GRCh38)/ukb21007_c%s_b0_v1.pvar'", this_chr))
  system(sprintf("dx download 'Bulk/Imputation/Imputation from genotype (TOPmed)/helper_files/ukb21007_c%s_b0_v1.sites.vcf.gz'", this_chr))
}

# Add rsIDs from helper files to pvar files (all have ID == '.' otherwise)
for(this_chr in c(1:22, "X")) { 
  pvar <- fread(sprintf("ukb21007_c%s_b0_v1.pvar", this_chr))
  pvar[, `#CHROM` := as.character(`#CHROM`)]

  info <- fread(sprintf("ukb21007_c%s_b0_v1.sites.vcf.gz", this_chr))  
  info[, `#CHROM` := gsub("chr", "", `#CHROM`)]
  pvar[info, on = .(`#CHROM`, POS, REF, ALT), ID := i.ID]
  stopifnot(pvar[ID == ".", .N] == 0) # Check to make sure we don't need to swap alleles or strand
  
  fwrite(pvar, sep="\t", quote=FALSE, file=sprintf("ukb21007_c%s_b0_v1.pvar", this_chr))
}

# Double check to make sure we don't have duplicate SNPs in the pvar files
# (this was previously an issue with a handful of SNPs in the HRCUK10K 
# imputation on CSD3)
dups <- foreach(this_chr = c(1:22, "X"), .combine=rbind) %do% {
  gc()
  # Load in variant information
  pvar <- fread(sprintf("ukb21007_c%s_b0_v1.pvar", this_chr))
  pvar[, `#CHROM` := as.character(`#CHROM`)]
  
  # Find duplicate SNPs we need to strip out
  dups_by_rsid <- pvar[, .N, by=.(ID, REF, ALT)][N > 1]
  dups_by_pos <- pvar[,.N,by=.(POS, REF, ALT)][N > 1]
  dups <- rbind(
    pvar[dups_by_rsid, on = .(ID, REF, ALT)],
    pvar[dups_by_pos, on = .(POS, REF, ALT)]
  )
  dups[, N := NULL]
  dups <- unique(dups)
  return(dups)
}
stopifnot(nrow(dups) == 0) 

# Curate the variant information from the VCF helper files to match what we have
# for the HRCUK10K imputation
for (this_chr in c(1:22, "X")) {
  pvar <- fread(sprintf("ukb21007_c%s_b0_v1.pvar", this_chr))
  pvar[, `#CHROM` := as.character(`#CHROM`)]
  
  info <- fread(sprintf("ukb21007_c%s_b0_v1.sites.vcf.gz", this_chr))  
  info[, `#CHROM` := gsub("chr", "", `#CHROM`)]
  
  info[INFO %like% ".*;.*;.*;.*", c("AF", "R2", "ER2", "TYPE") := tstrsplit(INFO, ";")] # Directly genotyped SNPs
  info[!(INFO %like% ".*;.*;.*;.*"), c("AF", "R2", "TYPE") := tstrsplit(INFO, ";")] # Imputed SNPs
  info[, AF := as.numeric(gsub(".*=", "", AF))]
  info[, R2 := as.numeric(gsub(".*=", "", R2))]
  
  snpstats <- pvar[,.(
    alternate_id=sprintf("%s:%s_%s_%s", `#CHROM`, POS, REF, ALT),
    rsid=ID, chromosome=`#CHROM`, position=POS, ref=REF, alt=ALT
  )]
  
  snpstats[info, on = .(chromosome=`#CHROM`, position=POS, ref=REF, alt=ALT),
    c("eaf", "info", "imputed") := .(AF, R2, TYPE == "IMPUTED")]
  
  snpstats[eaf > 0.5, c("minor_allele", "maf") := .(ref, 1 - eaf)]
  snpstats[eaf <= 0.5, c("minor_allele", "maf") := .(alt, eaf)]

  snpstats <- snpstats[,.(alternate_id, rsid, chromosome, position, ref, alt,
                          maf, minor_allele, info, imputed)]
  
  fwrite(snpstats, sep="\t", quote=FALSE, file=sprintf("ukb21007_c%s_b0_v1_varstats.txt", this_chr))
}

# Clean up intermediate files on persistent storage
system("dx mv 'common/Imputed Genotypes/TopMed (GRCh38)/*.log' trash/")
system("dx mv 'common/Imputed Genotypes/TopMed (GRCh38)/*.pvar' trash/")

# Upload files
system("dx upload *_varstats.txt *.pvar --destination 'common/Imputed Genotypes/TopMed (GRCh38)/'")
