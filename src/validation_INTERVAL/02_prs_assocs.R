library(data.table)
library(survival)
library(foreach)
source("src/functions/glm_test.R")
source("src/functions/factor.R") # set reference group (first level) without specifying all levels
source("src/functions/factor_by_size.R") # set reference group (first level) as largest group

# Load in the phenotype data
pheno <- fread("data/INTERVAL/curated_pheno.tsv")

# Drop people with prevalent diabetes (mix of T1D and T2D)
pheno <- pheno[!(prevalent_diabetes)]

# Load in all PRS to test
pgs <- fread("data/INTERVAL/PRS_levels/all_PRSs.sscore.gz")
setnames(pgs, gsub("_hmPOS_GRCh37", "", names(pgs)))

# Add to phenotype data
pheno <- pheno[pgs, on = .(IID), nomatch=0]

# Get PCs
pcs <- fread("data/INTERVAL/genetic_reference/annot_INT_50PCs_pcs.txt")
setnames(pcs, "ID", "IID")
setnames(pcs, gsub("_", "", names(pcs)))
pcs <- pcs[,.(IID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, 
              PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20)]

# Add to phenotype data
pheno <- pheno[pcs, on = .(IID), nomatch=0]

# Process risk factors and covariates prior to modelling - we want all
# estimates to be per SD change in variable, and for factors, using the
# largest group or lowest risk as reference 
pheno[, age := scale(age)]
pheno[, sex := ifelse(sex == 2, "Female", "Male")]
pheno[, sex := factor(sex, reference="Female")]

# For the HuertaChagoya2023, we need to combine their three scores as they describe at
# https://www.pgscatalog.org/score/PGS003443/
# https://www.pgscatalog.org/score/PGS003444/
# https://www.pgscatalog.org/score/PGS003445/
pheno[, HuertaChagoya2023 := scale(
  scale(PGS003443)*0.531117 +
    scale(PGS003444)*0.5690198 +
    scale(PGS003445)*0.1465538
)]

# Adjust pgs for 20 PCs
pgs_list <- c(intersect(names(pgs)[-1], names(pheno)), "HuertaChagoya2023")
for (this_pgs in pgs_list) {
  setnames(pheno, this_pgs, "this_pgs")
  pheno[, this_pgs := scale(lm(scale(this_pgs) ~ PC1 + PC2 + PC3 + PC4 + 
                                 PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + 
                                 PC16 + PC17 + PC18 + PC19 + PC20)$residuals)]
  setnames(pheno, "this_pgs", this_pgs)
}

mf <- "incident_diabetes ~ %s + age + sex"
assocs <- foreach(this_pgs = pgs_list, .combine=rbind) %do% {
  res <- glm.test(sprintf(mf, this_pgs), "incident_diabetes", pheno)
  res[, PRS := this_pgs]
}

# Rename coefficients for readability
assocs[!(coefficient %in% c("age", "sexMale", "(Intercept)")), coefficient := "PRS"]
assocs[coefficient == "age", coefficient := "Age"]
assocs[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]

# Write out
fwrite(assocs, sep="\t", quote=FALSE, file="output/INTERVAL_tests/incident_diabetes_associations.txt")
