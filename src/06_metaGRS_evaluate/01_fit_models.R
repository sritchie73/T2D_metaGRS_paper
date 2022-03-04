library(data.table)
library(survival)
library(foreach)
library(RNOmni)
source("src/functions/glm_test.R")
source("src/functions/cox_test.R")
source("src/functions/factor.R") # set reference group (first level) without specifying all levels
source("src/functions/factor_by_size.R") # set reference group (first level) as largest group
source("src/functions/logit.R") # like log transform, but for percentages

# Make output directory
system("mkdir -p output/UKB_tests", wait=TRUE)

# Load phenotype data
pheno <- fread("data/UKB/collated_curated_data.txt", tmpdir="tmp", na.strings=c("", "NA"))
pheno <- pheno[(metaGRS_test_samples)]

# Load metaGRS and existing T2D PGS
pgs <- fread("data/UKB/T2D_PGS/collated_scores.sscore.gz")
metaGRS <- fread("data/UKB/T2D_metaGRS_levels.txt")
pgs <- metaGRS[pgs, on = .(eid=IID)]

# Add to phenotype data
pheno <- pheno[pgs, on = .(eid), nomatch=0]

# Adjust pgs for 20 PCs
for (this_pgs in names(pgs)[-1]) {
  setnames(pheno, this_pgs, "this_pgs")
  pheno[, this_pgs := scale(lm(scale(this_pgs) ~ PC1 + PC2 + PC3 + PC4 + 
    PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + 
    PC16 + PC17 + PC18 + PC19 + PC20)$residuals)]
  setnames(pheno, "this_pgs", this_pgs)
}

# Ye2021_PGS001357 has effect allele and other allele incorrectly reported to PGS Catalog;
# as-is the score has strong inverse correlation with T2D status, so we need to flip the 
# sign.
pheno[, Ye2021_PGS001357 := -Ye2021_PGS001357]

# Process risk factors and covariates prior to modelling - we want all
# estimates to be per SD change in variable, and for factors, using the
# largest group or lowest risk as reference 
pheno[, age := scale(age)]
pheno[, sex := factor(sex, reference="Female")]
pheno[, bmi := scale(log(bmi))]
pheno[!is.na(townsend), townsend := RankNorm(townsend)]
pheno[, smoking_status := factor(smoking_status, reference="Non")]
pheno[, family_history_diabetes := factor(family_history_diabetes, reference=FALSE)]
pheno[, history_cvd := factor(history_cvd, reference=FALSE)]
pheno[, history_gestational_diabetes := factor(history_gestational_diabetes, reference=FALSE)]
pheno[, history_pcos := factor(history_pcos, reference=FALSE)]
pheno[, history_learning_difficulties := factor(history_learning_difficulties, reference=FALSE)]
pheno[, history_bipolar_schizophrenia := factor(history_bipolar_schizophrenia, reference=FALSE)]
pheno[, hypertension_medication := factor(hypertension_medication, reference=FALSE)]
pheno[, lipid_lowering_medication := factor(lipid_lowering_medication, reference=FALSE)]
pheno[, systematic_corticosteroids := factor(systematic_corticosteroids, reference=FALSE)]
pheno[, atypical_antipsychotics := factor(atypical_antipsychotics, reference=FALSE)]
pheno[, hba1c := scale(log(hba1c))]
pheno[, non_fasting_glucose := scale(log(non_fasting_glucose))]
pheno[, fasting_glucose := scale(log(fasting_glucose))]
pheno[, QDiabetes2013 := scale(logit(QDiabetes2013/100))]
pheno[, QDiabetes2018A := scale(logit(QDiabetes2018A/100))]
pheno[, QDiabetes2018B_fasting := scale(logit(QDiabetes2018B_fasting/100))]
pheno[, QDiabetes2018B_non_fasting := scale(logit(QDiabetes2018B_non_fasting/100))]
pheno[, QDiabetes2018C := scale(logit(QDiabetes2018C/100))]
pheno[, assessment_centre := factor_by_size(assessment_centre)]
pheno[, censor_hospital_nation := factor_by_size(censor_hospital_nation)]
pheno[is.na(earliest_hospital_nation) & !(type_2_diabetes), earliest_hospital_nation := assessment_nation] # N=5 Prevalent T2D cases from self-report who've withdrawn consent for hospital record linkage
pheno[, earliest_hospital_nation := factor_by_size(earliest_hospital_nation)]

# Drop participants free of T2D at baseline who have withdrawn consent for hospital record
# linkage and people with uncertain diabetes status
pheno <- pheno[!is.na(type_2_diabetes) & !is.na(earliest_hospital_nation)]

# First examine prevalent T2D
prev <- rbind(idcol="model",
  "age"=glm.test(type_2_diabetes ~ age, "type_2_diabetes", pheno),
  "sex"=glm.test(type_2_diabetes ~ sex, "type_2_diabetes", pheno),
  "follow_diff"=glm.test(type_2_diabetes ~ assessment_centre + earliest_hospital_nation, "type_2_diabetes", pheno),
  "age + sex"=glm.test(type_2_diabetes ~ age + sex, "type_2_diabetes", pheno),
  "age + follow_diff"=glm.test(type_2_diabetes ~ age + assessment_centre + earliest_hospital_nation, "type_2_diabetes", pheno),
  "sex + follow_diff"=glm.test(type_2_diabetes ~ sex + assessment_centre + earliest_hospital_nation, "type_2_diabetes", pheno),
  "age + sex + follow_diff"=glm.test(type_2_diabetes ~ age + sex + assessment_centre + earliest_hospital_nation, "type_2_diabetes", pheno)
)

mf <- "type_2_diabetes ~ %s + age + sex + assessment_centre + earliest_hospital_nation"
prs_prev <- foreach(this_pgs = names(pgs)[-1], .combine=rbind) %do% {
  res <- glm.test(sprintf(mf, this_pgs), "type_2_diabetes", pheno)
  res[, model := this_pgs]
} 

prev <- rbind(idcol="model_type", "reference"=prev, "pgs"=prs_prev)

# Rename coefficients for readability
prev[coefficient == "age", coefficient := "Age"]
prev[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]
prev[coefficient %like% "earliest_hospital_nation", coefficient := 
  sprintf("Earliest hospital nation: %s vs. England", gsub("earliest_hospital_nation", "", coefficient))]
prev[coefficient %like% "assessment_centre", coefficient :=
  sprintf("Assessment centre: %s vs. %s", gsub("(assessment_centre)|( \\(.*\\))", "", coefficient),
  levels(pheno$assessment_centre)[1])]

# Write out
fwrite(prev, sep="\t", quote=FALSE, file="output/UKB_tests/prevalent_T2D_associations.txt")

# Now do incident T2D
y <- "Surv(incident_censor_years, incident_type_2_diabetes)"
follow_diff <- "assessment_centre + earliest_hospital_nation + censor_hospital_nation"
inci <- rbind(idcol="model",
  "age"=cox.test(paste(y, "~ age"), "incident_type_2_diabetes", pheno),
  "sex"=cox.test(paste(y, "~ sex"), "incident_type_2_diabetes", pheno),
  "follow_diff"=cox.test(paste(y, "~", follow_diff), "incident_type_2_diabetes", pheno),
  "age + sex"=cox.test(paste(y, "~ age + sex"), "incident_type_2_diabetes", pheno),
  "age + follow_diff"=cox.test(paste(y, "~ age +", follow_diff), "incident_type_2_diabetes", pheno),
  "sex + follow_diff"=cox.test(paste(y, "~ sex +", follow_diff), "incident_type_2_diabetes", pheno),
  "age + sex + follow_diff"=cox.test(paste(y, "~ age + sex +", follow_diff), "incident_type_2_diabetes", pheno),
  "strata(sex) + age"=cox.test(paste(y, "~ age + strata(sex)"), "incident_type_2_diabetes", pheno),
  "strata(sex) + follow_diff"=cox.test(paste(y, "~ strata(sex) +", follow_diff), "incident_type_2_diabetes", pheno),
  "strata(sex) + age + follow_diff"=cox.test(paste(y, "~ age + strata(sex) +", follow_diff), "incident_type_2_diabetes", pheno)
)

# Treating sex as a factor variable instead of stratifying for different baseline hazards between males and females
# consistently gives better C-index, so will use that going forward
rf <- c("bmi", "smoking_status", "townsend", "family_history_diabetes", "history_cvd", "history_gestational_diabetes", 
        "history_pcos", "history_learning_difficulties", "history_bipolar_schizophrenia", "hypertension_medication", 
        "lipid_lowering_medication", "systematic_corticosteroids", "atypical_antipsychotics", "hba1c", "fasting_glucose", 
        "non_fasting_glucose")
mf <- paste(y, "~ %s + age + sex +", follow_diff)
rf_inci <- foreach(this_rf = rf, .combine=rbind) %do% {
  res <- cox.test(sprintf(mf, this_rf), "incident_type_2_diabetes", pheno)
  res[, model := this_rf]
} 

qdiab <- c("QDiabetes2018A", "QDiabetes2018B_fasting", "QDiabetes2018B_non_fasting", "QDiabetes2018C", "QDiabetes2013")
qd_inci <- foreach(this_rs = qdiab, .combine=rbind) %do% {
  res <- cox.test(sprintf(mf, this_rs), "incident_type_2_diabetes", pheno)
  res[, model := this_rs]
}

qd_inci2 <- rbind(idcol = "model",
  "QDiabetes2018B_fast_gte_3"=cox.test(sprintf(mf, "QDiabetes2018B_non_fasting"), "incident_type_2_diabetes", pheno[fasting_time >= 3]),
  "QDiabetes2018B_fast_gte_8"=cox.test(sprintf("%s ~ %s + age + sex", y, "QDiabetes2018B_non_fasting"), "incident_type_2_diabetes", pheno[fasting_time >= 8])
)
qd_inci <- rbind(qd_inci, qd_inci2)

pgs_inci <- foreach(this_pgs = names(pgs)[-1], .combine=rbind) %do% {
  res <- cox.test(sprintf(mf, this_pgs), "incident_type_2_diabetes", pheno)
  res[, model := this_pgs]
}

inci <- rbind(idcol="model_type", "reference"=inci, "risk factor"=rf_inci, "QDiabetes"=qd_inci, "PGS"=pgs_inci)

# Rename coefficients for readability
inci[coefficient == "age", coefficient := "Age"]
inci[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]
inci[coefficient %like% "earliest_hospital_nation", coefficient := 
  sprintf("Earliest hospital nation: %s vs. England", gsub("earliest_hospital_nation", "", coefficient))]
inci[coefficient %like% "assessment_centre", coefficient :=
  sprintf("Assessment centre: %s vs. %s", gsub("(assessment_centre)|( \\(.*\\))", "", coefficient),
  levels(pheno$assessment_centre)[1])]
inci[coefficient %like% "censor_hospital_nation", coefficient := 
  sprintf("Censor hospital nation: %s vs. England", gsub("censor_hospital_nation", "", coefficient))]
inci[coefficient %like% "smoking_status", coefficient :=
  sprintf("Smoking status: %s vs. Non-smoker", gsub("smoking_status", "", coefficient))]
inci[coefficient == "bmi", coefficient := "Body mass index"]
inci[coefficient == "townsend", coefficient := "Townsend deprivation index"]
inci[coefficient %like% "QDiabetes", coefficient := gsub("QDiabetes", "QDiabetes ", coefficient)]
inci[coefficient %like% "QDiabetes", coefficient := gsub("2018", "2018 model ", coefficient)]
inci[coefficient %like% "QDiabetes", coefficient := gsub("B_fasting", "B (fasting glucose)", coefficient)]
inci[coefficient %like% "QDiabetes", coefficient := gsub("B_non_fasting", "B (non-fasting glucose)", coefficient)]
inci[coefficient == "fasting_glucose", coefficient := "Glucose (fasting)"]
inci[coefficient == "non_fasting_glucose", coefficient := "Glucose (non-fasting)"]
inci[coefficient == "hba1c", coefficient := "Glycated haemoglobin (HbA1c)"]
inci[coefficient == "family_history_diabetesTRUE", coefficient := "Family history of diabetes"]
inci[coefficient == "history_cvdTRUE", coefficient := "History of cardiovascular disease"]
inci[coefficient == "history_gestational_diabetesTRUE", coefficient := "History of gestational diabetes"]
inci[coefficient == "history_pcosTRUE", coefficient := "History of polycystic ovary syndrome"]
inci[coefficient == "history_learning_difficultiesTRUE", coefficient := "History of learning difficulties"]
inci[coefficient == "history_bipolar_schizophreniaTRUE", coefficient := "History of bipolar or schizophrenia disorders"]
inci[coefficient == "hypertension_medicationTRUE", coefficient := "Taking hypertension medication"]
inci[coefficient == "lipid_lowering_medicationTRUE", coefficient := "Taking lipid lowering medication"]
inci[coefficient == "systematic_corticosteroidsTRUE", coefficient := "Taking systematic corticosteroids"]
inci[coefficient == "atypical_antipsychoticsTRUE", coefficient := "Taking 2nd generation atypical antipsychotics"]
inci[model == "QDiabetes2018B_fast_gte_3" & coefficient %like% "QDiabetes", coefficient := "QDiabetes 2018 model B (fasting ≥ 3hours)"]
inci[model == "QDiabetes2018B_fast_gte_8" & coefficient %like% "QDiabetes", coefficient := "QDiabetes 2018 model B (fasting ≥ 8hours)"]

# Write out
fwrite(inci, sep="\t", quote=FALSE, file="output/UKB_tests/incident_T2D_associations.txt")

# Do it all again, but this time filtering to consistent dataset (QDiabetes2018C - the strongest QDiabetes model)
pheno <- pheno[!is.na(QDiabetes2018C)]

inci <- rbind(idcol="model",
  "age"=cox.test(paste(y, "~ age"), "incident_type_2_diabetes", pheno),
  "sex"=cox.test(paste(y, "~ sex"), "incident_type_2_diabetes", pheno),
  "follow_diff"=cox.test(paste(y, "~", follow_diff), "incident_type_2_diabetes", pheno),
  "age + sex"=cox.test(paste(y, "~ age + sex"), "incident_type_2_diabetes", pheno),
  "age + follow_diff"=cox.test(paste(y, "~ age +", follow_diff), "incident_type_2_diabetes", pheno),
  "sex + follow_diff"=cox.test(paste(y, "~ sex +", follow_diff), "incident_type_2_diabetes", pheno),
  "age + sex + follow_diff"=cox.test(paste(y, "~ age + sex +", follow_diff), "incident_type_2_diabetes", pheno),
  "strata(sex) + age"=cox.test(paste(y, "~ age + strata(sex)"), "incident_type_2_diabetes", pheno),
  "strata(sex) + follow_diff"=cox.test(paste(y, "~ strata(sex) +", follow_diff), "incident_type_2_diabetes", pheno),
  "strata(sex) + age + follow_diff"=cox.test(paste(y, "~ age + strata(sex) +", follow_diff), "incident_type_2_diabetes", pheno)
)

rf_inci <- foreach(this_rf = rf, .combine=rbind) %do% {
  res <- cox.test(sprintf(mf, this_rf), "incident_type_2_diabetes", pheno)
  res[, model := this_rf]
} 

qd_inci <- foreach(this_rs = qdiab, .combine=rbind) %do% {
  res <- cox.test(sprintf(mf, this_rs), "incident_type_2_diabetes", pheno)
  res[, model := this_rs]
}

qd_inci2 <- rbind(idcol = "model",
  "QDiabetes2018B_fast_gte_3"=cox.test(sprintf(mf, "QDiabetes2018B_non_fasting"), "incident_type_2_diabetes", pheno[fasting_time >= 3]),
  "QDiabetes2018B_fast_gte_8"=cox.test(sprintf("%s ~ %s + age + sex", y, "QDiabetes2018B_non_fasting"), "incident_type_2_diabetes", pheno[fasting_time >= 8])
)
qd_inci <- rbind(qd_inci, qd_inci2)

pgs_inci <- foreach(this_pgs = names(pgs)[-1], .combine=rbind) %do% {
  res <- cox.test(sprintf(mf, this_pgs), "incident_type_2_diabetes", pheno)
  res[, model := this_pgs]
}

inci <- rbind(idcol="model_type", "reference"=inci, "risk factor"=rf_inci, "QDiabetes"=qd_inci, "PGS"=pgs_inci)

inci[coefficient == "age", coefficient := "Age"]
inci[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]
inci[coefficient %like% "earliest_hospital_nation", coefficient := 
  sprintf("Earliest hospital nation: %s vs. England", gsub("earliest_hospital_nation", "", coefficient))]
inci[coefficient %like% "assessment_centre", coefficient :=
  sprintf("Assessment centre: %s vs. %s", gsub("(assessment_centre)|( \\(.*\\))", "", coefficient),
  levels(pheno$assessment_centre)[1])]
inci[coefficient %like% "censor_hospital_nation", coefficient := 
  sprintf("Censor hospital nation: %s vs. England", gsub("censor_hospital_nation", "", coefficient))]
inci[coefficient %like% "smoking_status", coefficient :=
  sprintf("Smoking status: %s vs. Non-smoker", gsub("smoking_status", "", coefficient))]
inci[coefficient == "bmi", coefficient := "Body mass index"]
inci[coefficient == "townsend", coefficient := "Townsend deprivation index"]
inci[coefficient %like% "QDiabetes", coefficient := gsub("QDiabetes", "QDiabetes ", coefficient)]
inci[coefficient %like% "QDiabetes", coefficient := gsub("2018", "2018 model ", coefficient)]
inci[coefficient %like% "QDiabetes", coefficient := gsub("B_fasting", "B (fasting glucose)", coefficient)]
inci[coefficient %like% "QDiabetes", coefficient := gsub("B_non_fasting", "B (non-fasting glucose)", coefficient)]
inci[coefficient == "fasting_glucose", coefficient := "Glucose (fasting)"]
inci[coefficient == "non_fasting_glucose", coefficient := "Glucose (non-fasting)"]
inci[coefficient == "hba1c", coefficient := "Glycated haemoglobin (HbA1c)"]
inci[coefficient == "family_history_diabetesTRUE", coefficient := "Family history of diabetes"]
inci[coefficient == "history_cvdTRUE", coefficient := "History of cardiovascular disease"]
inci[coefficient == "history_gestational_diabetesTRUE", coefficient := "History of gestational diabetes"]
inci[coefficient == "history_pcosTRUE", coefficient := "History of polycystic ovary syndrome"]
inci[coefficient == "history_learning_difficultiesTRUE", coefficient := "History of learning difficulties"]
inci[coefficient == "history_bipolar_schizophreniaTRUE", coefficient := "History of bipolar or schizophrenia disorders"]
inci[coefficient == "hypertension_medicationTRUE", coefficient := "Taking hypertension medication"]
inci[coefficient == "lipid_lowering_medicationTRUE", coefficient := "Taking lipid lowering medication"]
inci[coefficient == "systematic_corticosteroidsTRUE", coefficient := "Taking systematic corticosteroids"]
inci[coefficient == "atypical_antipsychoticsTRUE", coefficient := "Taking 2nd generation atypical antipsychotics"]
inci[model == "QDiabetes2018B_fast_gte_3" & coefficient %like% "QDiabetes", coefficient := "QDiabetes 2018 model B (fasting ≥ 3hours)"]
inci[model == "QDiabetes2018B_fast_gte_8" & coefficient %like% "QDiabetes", coefficient := "QDiabetes 2018 model B (fasting ≥ 8hours)"]

# Write out
fwrite(inci, sep="\t", quote=FALSE, file="output/UKB_tests/incident_T2D_associations_QDiabetes2018C_subset.txt")

# Compare how PGS add to QDiabetes 2018 model C
mf <- paste(y, "~ %s + QDiabetes2018C + age + sex +", follow_diff)
multi_inci <- foreach(this_pgs = names(pgs)[-1], .combine=rbind) %do% {
  res <- cox.test(sprintf(mf, this_pgs), "incident_type_2_diabetes", pheno)
  res[, model := this_pgs]
}

# And to HbA1c and BMI
mf <- paste(y, "~ %s + hba1c + bmi + age + sex +", follow_diff)
multi_inci2 <- foreach(this_pgs = names(pgs)[-1], .combine=rbind) %do% {
  res <- cox.test(sprintf(mf, this_pgs), "incident_type_2_diabetes", pheno)
  res[, model := this_pgs]
}

inci <- rbind(idcol="model_type", "PGS + QDiabetes2018C"=multi_inci, "PGS + BMI + HbA1c"=multi_inci2)

inci[coefficient == "age", coefficient := "Age"]
inci[coefficient == "sexMale", coefficient := "Sex: Male vs. Female"]
inci[coefficient %like% "earliest_hospital_nation", coefficient :=
  sprintf("Earliest hospital nation: %s vs. England", gsub("earliest_hospital_nation", "", coefficient))]
inci[coefficient %like% "assessment_centre", coefficient :=
  sprintf("Assessment centre: %s vs. %s", gsub("(assessment_centre)|( \\(.*\\))", "", coefficient),
  levels(pheno$assessment_centre)[1])]
inci[coefficient %like% "censor_hospital_nation", coefficient :=
  sprintf("Censor hospital nation: %s vs. England", gsub("censor_hospital_nation", "", coefficient))]
inci[coefficient %like% "smoking_status", coefficient :=
  sprintf("Smoking status: %s vs. Non-smoker", gsub("smoking_status", "", coefficient))]
inci[coefficient == "bmi", coefficient := "Body mass index"]
inci[coefficient %like% "QDiabetes", coefficient := gsub("QDiabetes", "QDiabetes ", coefficient)]
inci[coefficient %like% "QDiabetes", coefficient := gsub("2018", "2018 model ", coefficient)]
inci[coefficient == "hba1c", coefficient := "Glycated haemoglobin (HbA1c)"]

fwrite(inci, sep="\t", quote=FALSE, file="output/UKB_tests/incident_T2D_multivariate_associations.txt")

