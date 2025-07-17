library(data.table)
library(foreach)

# Load each ancestry-specific metaPRS, and combine into a single multi-ancestry metaPRS
var_weights <- foreach(ancestry=c("EUR", "AFR", "AMR", "EAS", "SAS", "OTH"), .combine=rbind) %do% {
  anc_weights <- fread(sprintf("output/metaPRS/train/%s/T2D_%s_metaPRS.txt.gz", ancestry, ancestry))
  cbind("ancestry"=ancestry, anc_weights)
}

# Give each ancestry equal (1/6th) weight, n.b. 'mean()' explicitly not used here, as not all
# variants may contribute to all ancestries - treating these as 0s is most consistent with 
# other elements of metaPRS training
var_weights <- var_weights[,.(weight=sum(weight)/6), by=.(AoU_varID, chr, pos, effect_allele, other_allele)]

# write out
fwrite(var_weights, sep="\t", quote=FALSE, compress="gzip", file="output/metaPRS/train/T2D_multiancestry_metaPRS.txt.gz")

