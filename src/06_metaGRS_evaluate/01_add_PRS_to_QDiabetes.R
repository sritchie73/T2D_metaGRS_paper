library(data.table)
library(survival)
library(nricens)
library(ggplot2)
library(scales)

# Make output directory
system("mkdir -p output/UKB_tests", wait=TRUE)

# Load phenotype data 
pheno <- fread("data/UKB/collated_curated_data.txt", tmpdir="tmp", na.strings=c("", "NA"))
pheno <- pheno[(metaGRS_test_samples)]

# Load metaGRS 
metaGRS <- fread("data/UKB/T2D_metaGRS_levels.txt")
pheno <- pheno[metaGRS, on = .(eid), nomatch=0]

# Drop participants without QDiabetes data 
pheno <- pheno[!is.na(QDiabetes2018C)]

# Drop participants with uncertain/unknowable diabetes status
pheno <- pheno[!is.na(incident_type_2_diabetes)]

# add in age tertiles for later QDiabetes stratification analysis
pheno[!is.na(QDiabetes2018C), age_tertile := cut(age, breaks=quantile(age, probs=seq(0, 1, length=4)), include.lowest=TRUE)]

# Adjust PRS for 20 PCs
pheno[, T2D_metaGRS := scale(lm(scale(T2D_metaGRS) ~ PC1 + PC2 + PC3 + PC4 +
	PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 +
	PC16 + PC17 + PC18 + PC19 + PC20)$residuals)]

# Filter to columns we actually need
pheno <- pheno[,.(eid, sex, age, age_tertile, incident_type_2_diabetes, incident_censor_years,
  QDiabetes2018A, QDiabetes2018B_non_fasting, QDiabetes2018C, T2D_metaGRS)]

# Convert QDiabets risk scores from 0-100% to 0-1
pheno[, QDiabetes2018A := QDiabetes2018A/100]
pheno[, QDiabetes2018B_non_fasting := QDiabetes2018B_non_fasting/100]
pheno[, QDiabetes2018C := QDiabetes2018C/100]

# Function to convert QDiabetes absolute risks back to their original linear predictors
absrisk_to_lp <- function(QDRabsrisk, sex, QDRmodel) {
  stopifnot(QDRmodel %in% c("A", "B", "C"))
  stopifnot(all(sort(unique(sex)) == c("Female", "Male")))

  # First convert absrisk back to sex-specific linear predictor
  # Baseline hazards (here the base of the inner log function) are directly from the QDiabetes R package source code
  QDRLP <- rep(NA_real_, length(QDRabsrisk))
  if (QDRmodel == "A") {
    QDRLP[sex == "Male"] <- log(log(base=0.978732228279114, 1 - QDRabsrisk[sex == "Male"]))
    QDRLP[sex == "Female"] <- log(log(base=0.98622727394104, 1 - QDRabsrisk[sex == "Female"]))
  } else if (QDRmodel == "B") {
    QDRLP[sex == "Male"] <- log(log(base=0.985019445419312, 1 - QDRabsrisk[sex == "Male"]))
    QDRLP[sex == "Female"] <- log(log(base=0.990905702114105, 1 - QDRabsrisk[sex == "Female"]))
  } else if (QDRmodel == "C") {
    QDRLP[sex == "Male"] <- log(log(base=0.981181740760803, 1 - QDRabsrisk[sex == "Male"]))
    QDRLP[sex == "Female"] <- log(log(base=0.988788545131683, 1 - QDRabsrisk[sex == "Female"]))
  }

  return(QDRLP)
}

# Function to add T2D metaPRS to QDiabetes risk scores following naive method of Hageman et al.
# https://academic.oup.com/eurjpc/article/30/15/1705/7188647?login=true
add_risk <- function(QDRabsrisk, PRS, sex, time, event, QDRmodel="C") {
  stopifnot(QDRmodel %in% c("A", "B", "C"))
  stopifnot(all(sort(unique(sex)) == c("Female", "Male")))

  # First convert absrisk back to sex-specific linear predictor
  # Baseline hazards (here the base of the inner log function) are directly from the QDiabetes R package source code
  QDRLP <- absrisk_to_lp(QDRabsrisk, sex, QDRmodel)

  # Obtain sex-specific log Hazard Ratios for the PRS with QDiabetes linear predictor as the offset term
  male_logHR <- coxph(Surv(time[sex == "Male"], event[sex == "Male"]) ~ offset(QDRLP[sex == "Male"]) + scale(PRS[sex == "Male"]))$coefficient[1]
  female_logHR <- coxph(Surv(time[sex == "Female"], event[sex == "Female"]) ~ offset(QDRLP[sex == "Female"]) + scale(PRS[sex == "Female"]))$coefficient[1]

  # Calculate modified absolute risk
  new_risk <- rep(NA_real_, length(QDRabsrisk))
  new_risk[sex == "Male"] <- (1 - (1-QDRabsrisk[sex == "Male"])^exp(male_logHR * scale(PRS[sex == "Male"])))
  new_risk[sex == "Female"] <- (1 - (1-QDRabsrisk[sex == "Female"])^exp(female_logHR * scale(PRS[sex == "Female"])))

  # Return results
  attributes(new_risk) <- list("SubdistributionHR"=c("Male"=as.numeric(exp(male_logHR)), "Female"=as.numeric(exp(female_logHR))))
  return(new_risk)
}

# Add the T2D PRS to each QDiabetes risk model:
pheno[,QDiabA_plus_PRS := add_risk(QDiabetes2018A, T2D_metaGRS, sex, incident_censor_years, incident_type_2_diabetes, QDRmodel="A")]
pheno[!is.na(QDiabetes2018B_non_fasting), QDiabB_plus_PRS := add_risk(QDiabetes2018B_non_fasting, T2D_metaGRS, sex, incident_censor_years, incident_type_2_diabetes, QDRmodel="B")]
pheno[, QDiabC_plus_PRS := add_risk(QDiabetes2018C, T2D_metaGRS, sex, incident_censor_years, incident_type_2_diabetes, QDRmodel="C")]

# Write out new dataset
fwrite(pheno, sep="\t", quote=FALSE, file="output/UKB_tests/QDiabetes_plus_metaPRS_absrisks.txt")

# Curate subdistribution hazard ratio information
prs_SHR <- rbind(idcol="QDiabetes_model",
  "model A"=data.table(sex=c("Male", "Female"), HR=attr(pheno$QDiabA_plus_PRS, "SubdistributionHR")), 
  "model B"=data.table(sex=c("Male", "Female"), HR=attr(pheno$QDiabB_plus_PRS, "SubdistributionHR")), 
  "model C"=data.table(sex=c("Male", "Female"), HR=attr(pheno$QDiabC_plus_PRS, "SubdistributionHR"))
)
fwrite(prs_SHR, sep="\t", quote=FALSE, file="output/UKB_tests/QDiabetes_plus_metaPRS_subdistributionHRs.txt")

# Create table with QDiabetes and QDiabetes + metaPRS linear predictors (i.e. for downstream C-index calculation)
lps <- pheno[,.(eid, sex, age_tertile, incident_type_2_diabetes, incident_censor_years,
  QDiabA_LP=absrisk_to_lp(QDiabetes2018A, sex, "A"), QDiabA_PRS_LP=absrisk_to_lp(QDiabA_plus_PRS, sex, "A"),
  QDiabB_LP=absrisk_to_lp(QDiabetes2018B_non_fasting, sex, "B"), QDiabB_PRS_LP=absrisk_to_lp(QDiabB_plus_PRS, sex, "B"),
  QDiabC_LP=absrisk_to_lp(QDiabetes2018C, sex, "C"), QDiabC_PRS_LP=absrisk_to_lp(QDiabC_plus_PRS, sex, "C")
)]

fwrite(lps, sep="\t", quote=FALSE, file="output/UKB_tests/QDiabetes_plus_metaPRS_linear_predictors.txt")

# Sanity check empirical cumulative distributions of risk
eCDF_check <- melt(pheno, id.vars=c("eid", "sex", "age", "age_tertile", "incident_type_2_diabetes", "incident_censor_years"), value.name="absrisk", variable.name="model")
eCDF_check <- eCDF_check[model != "T2D_metaGRS"] # drop column - too lazy to write out full measure.vars
eCDF_check <- eCDF_check[!is.na(absrisk)]
eCDF_check[, QDiabetes_model := fcase(
  model %in% c("QDiabetes2018A", "QDiabA_plus_PRS"), "QDiabetes model A",
  model %in% c("QDiabetes2018B_non_fasting", "QDiabB_plus_PRS"), "QDiabetes model B",
  model %in% c("QDiabetes2018C", "QDiabC_plus_PRS"), "QDiabetes model C"
)]
eCDF_check[, model := ifelse(model %like% "plus_PRS", "QDiabetes + metaPRS", "QDiabetes")]
eCDF_check[, group := ifelse(incident_type_2_diabetes, "Case", "Control")]

# Plot probability of exceeding X% risk, stratified by case/control status
eCDF_check[, eCDF := ecdf(absrisk)(absrisk), by=.(QDiabetes_model, model, incident_type_2_diabetes)]
g <- ggplot(eCDF_check) +
  aes(x=absrisk, y=1-eCDF, color=model) +
  facet_grid(group ~ QDiabetes_model) +
  geom_line(linewidth=0.6) +
  scale_color_manual(values=c("QDiabetes"="#08306b", "QDiabetes + metaPRS"="#bd0026")) +
  geom_vline(xintercept=0.056, color="#fc4e2a", linetype=2, linewidth=0.6) +
  scale_x_continuous("Predicted 10-year T2D risk", labels=scales::percent, limits=c(0,0.3), oob=oob_keep, expand=expansion(mult=0, add=c(0.01, 0.03))) +
  scale_y_continuous("Probability of 10-year T2D risk exceeding X%") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )
ggsave(g, width=7.2, height=4.2, file="output/UKB_tests/QDiabetes_plus_metaPRS_absrisk_distribution_check.png")

# Plot probability of exceeding X% risk, stratified by case/control status and sex
eCDF_check[, eCDF := ecdf(absrisk)(absrisk), by=.(QDiabetes_model, model, incident_type_2_diabetes, sex)]
g <- ggplot(eCDF_check) +
  aes(x=absrisk, y=1-eCDF, color=model) +
  facet_grid(sex + group ~ QDiabetes_model) +
  geom_line(linewidth=0.6) +
  scale_color_manual(values=c("QDiabetes"="#08306b", "QDiabetes + metaPRS"="#bd0026")) +
  geom_vline(xintercept=0.056, color="#fc4e2a", linetype=2, linewidth=0.6) +
  scale_x_continuous("Predicted 10-year T2D risk", labels=scales::percent, limits=c(0,0.3), oob=oob_keep, expand=expansion(mult=0, add=c(0.01, 0.03))) +
  scale_y_continuous("Probability of 10-year T2D risk exceeding X%") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )
ggsave(g, width=7.2, height=7.2, file="output/UKB_tests/QDiabetes_plus_metaPRS_absrisk_distribution_check_by_sex.png")

# Plot probability of exceeding X% risk, stratified by case/control status and age-tertile
eCDF_check[, eCDF := ecdf(absrisk)(absrisk), by=.(QDiabetes_model, model, incident_type_2_diabetes, age_tertile)]
g <- ggplot(eCDF_check) +
  aes(x=absrisk, y=1-eCDF, color=model) +
  facet_grid(age_tertile + group ~ QDiabetes_model) +
  geom_line(linewidth=0.6) +
  scale_color_manual(values=c("QDiabetes"="#08306b", "QDiabetes + metaPRS"="#bd0026")) +
  geom_vline(xintercept=0.056, color="#fc4e2a", linetype=2, linewidth=0.6) +
  scale_x_continuous("Predicted 10-year T2D risk", labels=scales::percent, limits=c(0,0.3), oob=oob_keep, expand=expansion(mult=0, add=c(0.01, 0.03))) +
  scale_y_continuous("Probability of 10-year T2D risk exceeding X%") +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )
ggsave(g, width=7.2, height=7.2, file="output/UKB_tests/QDiabetes_plus_metaPRS_absrisk_distribution_check_by_age_tertile.png")

