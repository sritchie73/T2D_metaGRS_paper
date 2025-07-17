library(data.table)
library(foreach)

# Run for each ancestry separately
ancestry <- "EUR"

# Get list of GWASs
gwas_list <- fread('data/gwas_summary_stats/curated_gwas_table.txt')

# collate results
res <- foreach(this_gwas=gwas_list$PRS, .combine=rbind) %do% {
  qc_res <- readlines(sprintf("output/ldpred2/train/%s/%s/ldpred2_gwasqc_fail_rate.txt", ancestry, this_gwas))[1]
  data.table(GWAS=this_gwas, SNP_qc=qc_res)
}

# Order by fail rate
res[, pct_failed := as.numeric(gsub(".*\\(", "", gsub("%\\).*", "", SNP_qc)))]
res <- res[order(-pct_failed)]

# Write out 
res[, pct_failed := NULL]
fwrite(res, sep="\t", quote=FALSE, file=sprintf("output/ldpred2/train/%s/aggregated_ldpred2_gwasqc_fail_rate.txt", ancestry))

