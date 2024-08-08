library(data.table)
library(boot)
options(boot.parallel="multicore")
options(boot.ncpus=10) # Takes < 1 minute

# Make output directory
system("mkdir -p output/UKB_tests", wait=TRUE)

# Load relevant phenotype data
absrisks <- fread("output/UKB_tests/QDiabetes_plus_metaPRS_absrisks.txt")
absrisks <- absrisks[,.(eid, incident_type_2_diabetes, incident_censor_years, QDiabetes2018A, QDiabA_plus_PRS)] 

# Bootstrap function to return multiple statistics:
#  - Number of participants at >5.6% risk with QDiabetes model A
#  - Number of participants at >5.6% risk with QDiabetes model A + metaPRS
#  - Number of future T2D cases at >5.6% risk with QDiabetes model A
#  - Number of future T2D cases at >5.6% risk with QDiabetes model A + metaPRS
#  - Number needed to screen (NNS) per T2D event with QDiabetes model A
#  - NNS with QDiabetes model A + metaPRS
#  - delta-NNS 
boot.fun <- function(dt) {
  stats <- c(
    dt[QDiabetes2018A > 0.056, .N], 
    dt[QDiabA_plus_PRS > 0.056, .N],
    dt[(incident_type_2_diabetes) & QDiabetes2018A > 0.056, .N],
    dt[(incident_type_2_diabetes) & QDiabA_plus_PRS > 0.056, .N]
  )
  stats <- c(stats, stats[1]/stats[3], stats[2]/stats[4])
  stats <- c(stats, stats[6]-stats[5])
  return(stats)
}

# Run bootstrap analysis 
surv_cols_idx <- match(c("incident_censor_years", "incident_type_2_diabetes"), names(absrisks))
boot_res <- censboot(absrisks, boot.fun, 1000, index=surv_cols_idx)
saveRDS(boot_res, file="output/UKB_tests/QDiabA_screening_bootstraps.rds")

# Curate results
boot_stats <- data.table(
  metric=c("QDiabA_high_risk", "QDiabA_plus_PRS_high_risk", "QDiabA_high_risk_cases", "QDiabA_plus_PRS_high_risk_cases", "QDiabA_NNS", "QDiabA_plus_PRS_NNS", "delta_NNS"),
  estimate=boot_res$t0, boot_SE = apply(boot_res$t, 2, sd)
)
boot_stats[, boot_L95 := estimate - qnorm(1-(0.05/2))*boot_SE]
boot_stats[, boot_U95 := estimate + qnorm(1-(0.05/2))*boot_SE]
boot_stats[, boot_pval := pmin(1, pnorm(abs(estimate/boot_SE), lower.tail=FALSE)*2)]
fwrite(boot_stats, sep="\t", quote=FALSE, file="output/UKB_tests/QDiabA_screening_bootstraps.txt")



