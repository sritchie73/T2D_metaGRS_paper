library(data.table)
library(nricens)

# Do a two-stage screening comparison following NICE guidelines and NHS health check best practices;
# using QDiabetes model A (with risk threshold of 5.6%) to prioritise people for follow-up measurement
# of HbA1c or fasting glucose.
#
# After screening people with QDiabetes model A with or without PRS, explore different thresholds of
# risk using QDiabetes model B or model C
pheno <- fread("output/UKB_tests/QDiabetes_plus_metaPRS_absrisks.txt")

pheno[,QDiabA_screen_then_QDiabB := QDiabetes2018A]
pheno[QDiabA_screen_then_QDiabB > 0.056, QDiabA_screen_then_QDiabB := QDiabetes2018B_non_fasting]

pheno[,QDiabA_screen_then_QDiabC := QDiabetes2018A]
pheno[QDiabA_screen_then_QDiabC > 0.056, QDiabA_screen_then_QDiabC := QDiabetes2018C]

pheno[,QDiabA_screen_then_QDiabB_plus_PRS := QDiabA_plus_PRS]
pheno[QDiabA_screen_then_QDiabB_plus_PRS > 0.056, QDiabA_screen_then_QDiabB_plus_PRS := QDiabB_plus_PRS]

pheno[,QDiabA_screen_then_QDiabC_plus_PRS := QDiabA_plus_PRS]
pheno[QDiabA_screen_then_QDiabC_plus_PRS > 0.056, QDiabA_screen_then_QDiabC_plus_PRS := QDiabC_plus_PRS]

pheno[, c("QDiabetes2018B_non_fasting", "QDiabetes2018C", "T2D_metaGRS", "QDiabB_plus_PRS", "QDiabC_plus_PRS") := NULL]

fwrite(pheno, sep="\t", quote=FALSE, file="output/UKB_tests/screening_absrisks.txt")

# Reorganise for NRI analysis
pheno <- rbind(idcol="QDiabetes_model",
  "QDiabetes model A"=pheno[,.(eid, sex, age, age_tertile, incident_type_2_diabetes, incident_censor_years, without_PRS=QDiabetes2018A, with_PRS=QDiabA_plus_PRS)],
  "QDiabetes model B"=pheno[,.(eid, sex, age, age_tertile, incident_type_2_diabetes, incident_censor_years, without_PRS=QDiabA_screen_then_QDiabB, with_PRS=QDiabA_screen_then_QDiabB_plus_PRS)],
  "QDiabetes model C"=pheno[,.(eid, sex, age, age_tertile, incident_type_2_diabetes, incident_censor_years, without_PRS=QDiabA_screen_then_QDiabC, with_PRS=QDiabA_screen_then_QDiabC_plus_PRS)]
)
pheno <- pheno[!is.na(without_PRS) & !is.na(with_PRS)]

# Wrapper function for categorical NRI test
nri.test <- function(data, qdiab_model, risk_threshold) {
  # Extract predicted 10-year risk from QDiabetes and QDiabetes + metaPRS
  comp_risk <- data[QDiabetes_model == qdiab_model]

  # Run NRI analysis
  NRI <- nricens(event = comp_risk$incident_type_2_diabetes, time = comp_risk$incident_censor_years,
                 p.std = comp_risk$without_PRS, p.new = comp_risk$with_PRS,
                 updown = "category", cut = risk_threshold, t0 = 10, niter = 1000)
  system("rm -f Rplots.pdf")

  # Add in sample size and case numbers
  NRI$n <- comp_risk[,.N]
  NRI$nevent <- comp_risk[, sum(incident_type_2_diabetes)]

  return(NRI)
}

nri_lists <- list(
  "5.6%"=list(
     "QDiabetes model A"=nri.test(pheno, "QDiabetes model A", 0.056), # NRI for initial screening step
     "QDiabetes model B"=nri.test(pheno, "QDiabetes model B", 0.056),
     "QDiabetes model C"=nri.test(pheno, "QDiabetes model C", 0.056)
  ),
  "10%"=list(
     "QDiabetes model B"=nri.test(pheno, "QDiabetes model B", 0.1),
     "QDiabetes model C"=nri.test(pheno, "QDiabetes model C", 0.1)
  ),
  "13.2%"=list(
     "QDiabetes model B"=nri.test(pheno, "QDiabetes model B", 0.132)
  ),
  "14.3%"=list(
     "QDiabetes model C"=nri.test(pheno, "QDiabetes model C", 0.143)
  ),
  "15%"=list(
     "QDiabetes model B"=nri.test(pheno, "QDiabetes model B", 0.15),
     "QDiabetes model C"=nri.test(pheno, "QDiabetes model C", 0.15)
  )
)
saveRDS(nri_lists, file="output/UKB_tests/screening_nri.rds")

# Extract tables of estimates
nri_estimates <- rbindlist(idcol="risk_threshold", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="QDiabetes_model", fill=TRUE, lapply(l1, function(l2) {
    cbind(samples=l2$n, cases=l2$nevent, as.data.table(l2$nri, keep.rownames="metric"))
  }))
}))

# Compute bootstrap standard errors
nri_bsse <- rbindlist(idcol="risk_threshold", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="QDiabetes_model", fill=TRUE, lapply(l1, function(l2) {
    se_vec <- apply(l2$bootstrapsample, 2, sd)
    data.table(metric=names(se_vec), SE=se_vec)
  }))
}))
nri_estimates[nri_bsse, on = .(risk_threshold, QDiabetes_model, metric), SE := i.SE]

# Compute 95% CI and P-value from BSSE 
nri_estimates[, L95 := Estimate - qnorm(1-(0.05/2))*SE]
nri_estimates[, U95 := Estimate + qnorm(1-(0.05/2))*SE]
nri_estimates[, Pval := pmin(1, pnorm(abs(Estimate/SE), lower.tail=FALSE)*2)]

# Extract tables of reclassifications for categorical nris
reclassified_cases <- rbindlist(idcol="risk_threshold", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="QDiabetes_model", fill=TRUE, lapply(l1, function(l2) {
    as.data.table(l2[["rtab.case"]])
  }))
}))

reclassified <- rbindlist(idcol="risk_threshold", fill=TRUE, lapply(nri_lists, function(l1) {
  rbindlist(idcol="QDiabetes_model", fill=TRUE, lapply(l1, function(l2) {
    as.data.table(l2[["rtab"]])
  }))
}))

# Merge
setnames(reclassified, "N", "All")
reclassified[reclassified_cases, on = .(risk_threshold, QDiabetes_model, Standard, New), Cases := N]
setnames(reclassified, "Standard", "QDiabetes")
setnames(reclassified, "New", "QDiabetes_plus_PRS")

# Add in total sample size and total cases to reclassified table
reclassified[nri_estimates, on = .(risk_threshold, QDiabetes_model), Total_Samples := i.samples]
reclassified[nri_estimates, on = .(risk_threshold, QDiabetes_model), Total_Cases := i.cases]

fwrite(nri_estimates, sep="\t", quote=FALSE, file="output/UKB_tests/screening_nri_estimates.txt")
fwrite(reclassified, sep="\t", quote=FALSE, file="output/UKB_tests/screening_nri_reclassified.txt")

