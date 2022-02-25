library(data.table)

# Load pre-curated data from the QDiabetes folder
pheno <- fread("data/UKB/QDiabetes/output/qdiabetes.txt")

# Rank and classifying participants by T2D information - useful for downstream sanity checks
# and kinship filtering
pheno[, T2D_information := fcase(
  (incident_type_2_diabetes) & !(fasting_undiagnosed_diabetes) & !(fasting_prediabetes), "Incident T2D, clear control definition",
  (type_2_diabetes), "Prevalent T2D, clear control definition",
  !(incident_type_2_diabetes) & !(type_2_diabetes) & !(fasting_undiagnosed_diabetes) & !(fasting_prediabetes), "No prevalent or incident T2D, clear control definition",
  (incident_type_2_diabetes) & (fasting_undiagnosed_diabetes), "Incident T2D, but maybe undiagnosed diabetes",
  (incident_type_2_diabetes) & (fasting_prediabetes), "Incident T2D, but maybe pre-diabetes",
  !(type_2_diabetes) & !(incident_type_2_diabetes) & (fasting_undiagnosed_diabetes), "No prevalent or incident T2D, but maybe undiagnosed diabetes",
  !(type_2_diabetes) & !(incident_type_2_diabetes) & (fasting_prediabetes), "No prevalent or incident T2D, but maybe pre-diabetes",
  !is.na(type_2_diabetes) & is.na(incident_type_2_diabetes) & !(fasting_undiagnosed_diabetes) & !(fasting_prediabetes), "Undetermined incident T2D, but clear control definition",
  !is.na(type_2_diabetes) & is.na(incident_type_2_diabetes) & (fasting_undiagnosed_diabetes), "Undetermined incident T2D, maybe undiagnosed diabetes",
  !is.na(type_2_diabetes) & is.na(incident_type_2_diabetes) & (fasting_prediabetes), "Undetermined incident T2D, maybe pre-diabetes",
  is.na(type_2_diabetes) & !(fasting_undiagnosed_diabetes) & !(fasting_prediabetes), "Undetermined T2D status, but HbA1c and fasting glucose normal",
  is.na(type_2_diabetes) & (fasting_undiagnosed_diabetes), "Undetermined T2D status, maybe undiagnosed diabetes",
  is.na(type_2_diabetes) & (fasting_prediabetes), "Undetermined T2D status, maybe pre-diabetes"
)]

pheno[, T2D_information_rank := fcase(
  (incident_type_2_diabetes) & !(fasting_undiagnosed_diabetes) & !(fasting_prediabetes), 1,
  (type_2_diabetes), 2,
  !(incident_type_2_diabetes) & !(type_2_diabetes) & !(fasting_undiagnosed_diabetes) & !(fasting_prediabetes), 3,
  (incident_type_2_diabetes) & (fasting_undiagnosed_diabetes), 4,
  (incident_type_2_diabetes) & (fasting_prediabetes), 5,
  !(type_2_diabetes) & !(incident_type_2_diabetes) & (fasting_undiagnosed_diabetes), 6,
  !(type_2_diabetes) & !(incident_type_2_diabetes) & (fasting_prediabetes), 7,
  !is.na(type_2_diabetes) & is.na(incident_type_2_diabetes) & !(fasting_undiagnosed_diabetes) & !(fasting_prediabetes), 8,
  !is.na(type_2_diabetes) & is.na(incident_type_2_diabetes) & (fasting_undiagnosed_diabetes), 9,
  !is.na(type_2_diabetes) & is.na(incident_type_2_diabetes) & (fasting_prediabetes), 10,
  is.na(type_2_diabetes) & !(fasting_undiagnosed_diabetes) & !(fasting_prediabetes), 11,
  is.na(type_2_diabetes) & (fasting_undiagnosed_diabetes), 12,
  is.na(type_2_diabetes) & (fasting_prediabetes), 13
)]


# Drop latest sample withdrawals
withdrawals <- fread("data/UKB/latest_withdrawals/output/latest_withdrawals.txt")
pheno <- pheno[!withdrawals, on = .(eid)]

# Add in information about genetic sex and ancestry groupings
setnames(pheno, "ethnicity", "ethnicity_qdiabetes_groups")
anthro <- fread("data/UKB/anthropometrics/output/anthropometrics.txt", na.strings=c("NA", ""))
pheno[anthro, on = .(eid, visit_index), ethnicity_group := i.ethnicity_group]
pheno[anthro, on = .(eid, visit_index), ethnicity_subgroup := i.ethnicity_subgroup]
pheno[anthro, on = .(eid, visit_index), imputed_genotypes := ifelse(is.na(i.genetic_white_british), FALSE, TRUE)]
pheno[anthro, on = .(eid, visit_index), genetic_white_british := i.genetic_white_british]
pheno[anthro, on = .(eid, visit_index), genetic_sex := i.genetic_sex]

# Add in genotyping chip and first 20 PCs
genref <- fread("data/UKB/genetics/reference_files/ukb_sqc_v2.txt")
genref <- genref[, .(eid, chip=genotyping.array,
                     PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, PC11,
                     PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20)]
pheno <- merge(pheno, genref, by="eid", all.x=TRUE)

# Load genetic relatedness information, and flag samples with first degree relationships that we may
# want to remove downstream. When choosing from pairs of related samples, prioritise keeping T2D cases
# over non-cases.
kinship <- fread("data/UKB/genetics/reference_files/kinship_relatedness.txt")
kinship <- kinship[Kinship > 0.0884] # cutoff from KING manual http://people.virginia.edu/~wc9c/KING/manual.html
kinship <- kinship[ID1 %in% pheno$eid & ID2 %in% pheno$eid]
kinship[pheno, on = .(ID1=eid), ID1_rank := i.T2D_information_rank]
kinship[pheno, on = .(ID2=eid), ID2_rank := i.T2D_information_rank]
kinship[, to_drop := ifelse(ID1_rank > ID2_rank, "ID1", "ID2")]
kinship[to_drop == "ID1", to_drop_eid := ID1]
kinship[to_drop == "ID2", to_drop_eid := ID2]
pheno[, first_degree_relative_to_drop := eid %in% kinship$to_drop_eid]
pheno[!(imputed_genotypes), first_degree_relative_to_drop := NA]

# Flag training samples we will use for running ldpred2
ldpred2 <- fread("data/UKB/sample_splits/Stroke_metaGRS_training_samples.txt")
pheno[, ldpred2_samples := eid %in% ldpred2$eid_7439 & (genetic_white_british) & !(first_degree_relative_to_drop)]

# Flag UKB phase 1 samples, i.e. those used to train the Khera 2018
ukb_ph1 <- fread("data/UKB/sample_splits/UKB_phase1_interim_release.txt", skip=5)
pheno[, ukb_phase1_samples := eid %in% ukb_ph1$V1]

# Split data into training and test sets
pheno[, metaGRS_train_samples := (ukb_phase1_samples) & !(ldpred2_samples) & (genetic_white_british) & !(first_degree_relative_to_drop)]
pheno[, metaGRS_test_samples := !(ukb_phase1_samples) & !(ldpred2_samples) & (genetic_white_british) & !(first_degree_relative_to_drop)]

# Write out 
fwrite(pheno, sep="\t", quote=FALSE, file="data/UKB/collated_curated_data.txt")
