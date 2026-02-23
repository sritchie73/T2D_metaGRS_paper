# Install missing packages on the RAP
install.packages("doMC")
install.packages("pROC")
install.packages("R.utils")
remotes::install_github("sritchie73/dxutils")

# Load libraries
library(data.table)
library(survival)
library(foreach)
library(dxutils)
library(doMC)

# Set up parallel computation
registerDoMC(8)

# Download and load additional scripts needed
dx_download("src/functions/glm_test.R", "src/functions/")
dx_download("src/functions/factor.R", "src/functions/")
dx_download("src/functions/factor_by_size.R", "src/functions/")
source("src/functions/glm_test.R")
source("src/functions/factor.R") # set reference group (first level) without specifying all levels
source("src/functions/factor_by_size.R") # set reference group (first level) as largest group

# Download and load in the phenotype data
dx_download("data/ukb_t2d_pheno.tsv", "input_data/")
pheno <- fread("input_data/ukb_t2d_pheno.tsv")
pheno <- pheno[visit_index == 0]

# Download and load in all PRS to test
dx_download("PRS_levels/all_PRSs.sscore.gz", "input_data/")
pgs <- fread("input_data/all_PRSs.sscore.gz")
setnames(pgs, gsub("_hmPOS_GRCh38", "", names(pgs)))

# Add to phenotype data
pheno <- pheno[pgs, on = .(eid=IID), nomatch=0]

# Download and add genotype PCs
dx_download("common/Genetic Reference/pcs.csv", "input_data/")
pcs <- fread("input_data/pcs.csv")
pcs <- pcs[,.(eid, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, 
              PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20)]
pheno <- pheno[pcs, on = .(eid), nomatch=0]

# Create "OTH" ancestry for people who did not cluster into any genetic ancestry
pheno[genetic_ancestry == "", genetic_ancestry := "OTH"]

# Filter to unrelated individuals. When choosing from pairs of related samples, 
# prioritise keeping T2D cases over non-cases.
dx_download("common/Genetic Reference/kinship_relatness.txt", "input_data/")
kinship <- fread("input_data/kinship_relatness.txt")
kinship <- kinship[Kinship > 0.0884] # cutoff from KING manual for first-degree relatives http://people.virginia.edu/~wc9c/KING/manual.html
kinship <- kinship[ID1 %in% pheno$eid & ID2 %in% pheno$eid]
kinship[pheno, on = .(ID1=eid), ID1_t2d := i.t2d_case]
kinship[pheno, on = .(ID2=eid), ID2_t2d := i.t2d_case]
kinship[, to_drop := ifelse(ID1_t2d & !ID2_t2d, "ID1", "ID2")]
kinship[to_drop == "ID1", to_drop_eid := ID1]
kinship[to_drop == "ID2", to_drop_eid := ID2]
pheno <- pheno[!kinship, on = .(eid=to_drop_eid)]

# Loop through all ancestries to assess associations
assocs <- foreach(this_ancestry = pheno[,unique(genetic_ancestry)], .combine=rbind) %do% {
  cat(sprintf("Testing associations in %s samples...\n", this_ancestry))
  # Filter to this ancestry
  this_pheno <- pheno[genetic_ancestry == this_ancestry]
  
  # Process risk factors and covariates prior to modelling - we want all
  # estimates to be per SD change in variable, and for factors, using the
  # largest group or lowest risk as reference 
  this_pheno[, age := scale(age)]
  this_pheno[, sex := factor(sex, reference="Female")]
  this_pheno[, assessment_centre := factor_by_size(assessment_centre)]
  
  # For the HuertaChagoya2023, we need to combine their three scores as they describe at
  # https://www.pgscatalog.org/score/PGS003443/
  # https://www.pgscatalog.org/score/PGS003444/
  # https://www.pgscatalog.org/score/PGS003445/
  this_pheno[, HuertaChagoya2023 := scale(
    scale(PGS003443)*0.531117 +
      scale(PGS003444)*0.5690198 +
      scale(PGS003445)*0.1465538
  )]
  
  # Adjust pgs for 20 PCs
  pgs_list <- c(intersect(names(pgs)[-1], names(pheno)), "HuertaChagoya2023")
  for (this_pgs in pgs_list) {
    setnames(this_pheno, this_pgs, "this_pgs")
    this_pheno[, this_pgs := scale(lm(scale(this_pgs) ~ PC1 + PC2 + PC3 + PC4 + 
                                        PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + 
                                        PC16 + PC17 + PC18 + PC19 + PC20)$residuals)]
    setnames(this_pheno, "this_pgs", this_pgs)
  }
  
  mf <- "t2d_case ~ %s + age + sex + assessment_centre"
  foreach(this_pgs = pgs_list, .combine=rbind) %dopar% {
    cat(sprintf("Testing PRS %s of %s...\n", which(pgs_list == this_pgs), length(pgs_list)))
    res <- suppressMessages(glm.test(sprintf(mf, this_pgs), "t2d_case", this_pheno))
    cbind("genetic_ancestry"=this_ancestry, "PRS"=this_pgs, res)
  }
}

# Rename coefficients for readability
assocs[!(coefficient %in% c("age", "sexMale", "(Intercept)")), coefficient := "PRS"]
assocs[coefficient == "age", coefficient := "Age"]
assocs[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]

# Write out and upload
fwrite(assocs, sep="\t", quote=FALSE, file="t2d_assocs.tsv")
dx_upload("t2d_assocs.tsv", "output/")
