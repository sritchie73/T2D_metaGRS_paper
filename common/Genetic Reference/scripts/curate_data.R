library(data.table)

# For the per-SNP QC information and KING relatedness file, we can simply make
# copies in the relevant folder. Unfortunately it does not appear to be possible
# to make symbolic link or hard link in the project storage so that we can make
# these files more findable, and we have to do a download and upload as 'dx cp'
# only works between projects.
system("dx download 'Bulk/Genotype Results/Genotype calls/ukb_snp_qc.txt'")
system("dx upload ukb_snp_qc.txt --destination 'common/Genetic Reference/'")

system("dx download 'Bulk/Genotype Results/Genotype calls/ukb_rel.dat'")
system("dx upload ukb_rel.dat --destination 'common/Genetic Reference/kinship_relatness.txt'")

# For the per sample information in ukb_sqc_v2.txt we need to additionally give
# this headers, and while we're at it, we might as well split it out into 
# distinct files (e.g. so that the PCs are separate) - we also need to map the
# genetic IDs to the application-specific sample IDs, as the second column in
# ukb_sqc_v2.txt does not map to our eids - we'll use the IIDs in the fam files
# instead. If in doubt, check the genetic_id and eid column pair against those
# same columns in common/application_id_map.csv
system("dx download 'Bulk/Genotype Results/Genotype calls/ukb_sqc_v2.txt'")
system("dx download 'Bulk/Genotype Results/Genotype calls/ukb22418_c22_b0_v2.fam'")
fam <- fread("ukb22418_c22_b0_v2.fam")
sqc <- fread("ukb_sqc_v2.txt")
stopifnot(sqc[,.N] == fam[,.N]) # if not true, need to use different mapping file
setnames(sqc, c("genetic_id", "eid", "genotyping.array", "Batch", "Plate.Name", 
                "Well", "Cluster.CR", "dQC", "Internal.Pico..ng.uL.", "Submitted.Gender", 
                "Inferred.Gender", "X.intensity", "Y.intensity", "Submitted.Plate.Name", 
                "Submitted.Well", "sample.qc.missing.rate", "heterozygosity", 
                "heterozygosity.pc.corrected", "het.missing.outliers", 
                "putative.sex.chromosome.aneuploidy", "in.kinship.table", 
                "excluded.from.kinship.inference", "excess.relatives", 
                "in.white.British.ancestry.subset", "used.in.pca.calculation",
                paste0("PC", 1:40), "in.Phasing.Input.chr1_22", "in.Phasing.Input.chrX", 
                "in.Phasing.Input.chrXY"))
sqc[, eid := fam$V2] # replace application-unknown eids with our application's eids
sqc <- sqc[eid > 0] # drop withdrawn samples

# Save and upload PCs
pcs <- sqc[,c("eid", paste0("PC", 1:40))]
fwrite(pcs, file="pcs.csv")
system("dx upload pcs.csv --destination 'common/Genetic Reference/'")

# Extract and upload sex inferred from genotypes
sex <- sqc[,.(eid, genetic_sex=ifelse(Inferred.Gender == "M", "Male", "Female"))]
fwrite(sex, file="genetic_sex.csv")
system("dx upload genetic_sex.csv --destination 'common/Genetic Reference/'")

# Extract UKB's old White British genetic ancestry call - keep this aside to be
# combined with the new genomic ancestry assignments (Field 30079)
wb <- sqc[,.(eid, white_british=as.logical(in.white.British.ancestry.subset))]

# Subset and write out the rest of the sample qc information
sqc[, c(paste0("PC", 1:40), "Submitted.Gender", "Inferred.Gender", "in.white.British.ancestry.subset") := NULL]
sqc <- sqc[,c(2,1,3:ncol(sqc)), with=FALSE] # Reorder so eid is first column
fwrite(sqc, "extended_sample_qc_information.csv")
system("dx upload extended_sample_qc_information.csv --destination 'common/Genetic Reference/'")

# Pull down genetic ancestry and blood type information extracted by Table Exporter
system("mkdir -p raw_data", wait=TRUE)
system("dx download 'common/Genetic Reference/raw_data/data.csv' -o raw_data/genetic_reference.csv", wait=TRUE)

# Load in raw data and curated field information
raw <- fread("raw_data/genetic_reference.csv")
info <- fread("genetic_reference/field_information.csv")

# Convert field ids to variable names
setnames(raw, paste0("p", as.character(info$field.id)), info$var)

# Write out blood type information
blood_type <- raw[,.(eid, blood_type)]
fwrite(blood_type, "blood_type_haplotypes.csv")
system("dx upload blood_type_haplotypes.csv --destination 'common/Genetic Reference/'")

# Combine genomic ancestry with old white british ancestry assignments
ancestry <- raw[,.(eid, ancestry)]
ancestry[, ancestry := fcase(
  ancestry == 1L, "AFR",
  ancestry == 2L, "AMR",
  ancestry == 3L, "CSA",
  ancestry == 4L, "EAS",
  ancestry == 5L, "EUR",
  ancestry == 6L, "MID"
)]
ancestry <- merge(ancestry, wb, by="eid", all=TRUE)
fwrite(ancestry, "ancestry.csv")
system("dx upload ancestry.csv --destination 'common/Genetic Reference/'")

# Send raw data to deletion folder to reduce storage costs - this needs to be 
# done in two steps, since we can't move two folders with the same name to the
# same location, so here we rename with a random number and then move to trash
rn <- as.integer(Sys.time())
system(sprintf("dx mv 'common/Genetic Reference/raw_data/' 'common/Genetic Reference/raw_data_%s'", rn), wait=TRUE)
system(sprintf("dx mv 'common/Genetic Reference/raw_data_%s/' trash/", rn), wait=TRUE) 
