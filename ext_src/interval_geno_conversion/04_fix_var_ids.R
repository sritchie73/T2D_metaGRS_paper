library(data.table)

out_dir = "/rds/project/jmmh2/rds-jmmh2-post_qc_data/interval/imputed/uk10k_1000g_b37/imputed/plink_format/pgen"

chr_id = Sys.getenv("SLURM_ARRAY_TASK_ID")

# Remove extra row identifier tacked onto the variant IDs
pvar = fread(sprintf("%s/impute_dedup_%s_interval.pvar", out_dir, chr_id))
pvar[, ID := gsub("\\.[0-9]*$", "", ID)]

# For variants without rsIDs, give them IDs based on chromosome position and alleles:
pvar[ID == ".", ID := paste(`#CHROM`, POS, ALT, REF, sep=":")]

# write out updated information
fwrite(pvar, sep="\t", quote=FALSE, file=sprintf("%s/impute_dedup_%s_interval.pvar", out_dir, chr_id))

