library(data.table)
library(boot)
options(boot.parallel="multicore")
options(boot.ncpus=10) # Takes < 1 minute

# Make output directory
system("mkdir -p output/UKB_tests", wait=TRUE)

# Load relevant phenotype data
absrisks <- fread("output/UKB_tests/screening_absrisks.txt")
absrisks <- absrisks[,.(eid, incident_type_2_diabetes, incident_censor_years, QDiabA_screen_then_QDiabC, QDiabA_screen_then_QDiabC_plus_PRS)] 

# Bootstrap function to return multiple statistics:
#  - Number of participants at >10% risk with QDiabetes C after prioritisation with QDiabetes A at risk >5.6% 
#  - Number of participants at >10% risk with QDiabetes C + metaPRS after prioritisation with QDiabetes A + metaPRS at risk >5.6% 
#  - Number of future T2D cases at >10% risk with QDiabetes C after prioritisation with QDiabetes A at risk >5.6% 
#  - Number of future T2D cases at >10% risk with QDiabetes C + metaPRS after prioritisation with QDiabetes A + metaPRS at risk >5.6% 
#  - Number needed to treat (NNT) per T2D event at >10% risk with QDiabetes C after prioritisation with QDiabetes A at risk >5.6%
#  - NNT per T2D event at >10% risk with QDiabetes C + metaPRS after prioritisation with QDiabetes A + metaPRS at risk >5.6% 
#  - delta-NNT
boot.fun <- function(dt) {
  stats <- c(
    dt[QDiabA_screen_then_QDiabC > 0.143, .N], 
    dt[QDiabA_screen_then_QDiabC_plus_PRS > 0.143, .N],
    dt[(incident_type_2_diabetes) & QDiabA_screen_then_QDiabC > 0.143, .N],
    dt[(incident_type_2_diabetes) & QDiabA_screen_then_QDiabC_plus_PRS > 0.143, .N]
  )
  stats <- c(stats, stats[1]/stats[3], stats[2]/stats[4])
  stats <- c(stats, stats[6]-stats[5])
  return(stats)
}

# Run bootstrap analysis 
surv_cols_idx <- match(c("incident_censor_years", "incident_type_2_diabetes"), names(absrisks))
boot_res <- censboot(absrisks, boot.fun, 1000, index=surv_cols_idx)
saveRDS(boot_res, file="output/UKB_tests/QDiabA_then_C_screening_bootstraps.rds")

# Curate results
boot_stats <- data.table(
  metric=c("QDiabA_then_C_high_risk", "QDiabA_then_C_plus_PRS_high_risk", "QDiabA_then_C_high_risk_cases", "QDiabA_then_C_plus_PRS_high_risk_cases", "QDiabA_then_C_NNT", "QDiabA_then_C_plus_PRS_NNT", "delta_NNT"),
  estimate=boot_res$t0, boot_SE = apply(boot_res$t, 2, sd)
)
boot_stats[, boot_L95 := estimate - qnorm(1-(0.05/2))*boot_SE]
boot_stats[, boot_U95 := estimate + qnorm(1-(0.05/2))*boot_SE]
boot_stats[, boot_pval := pmin(1, pnorm(abs(estimate/boot_SE), lower.tail=FALSE)*2)]
fwrite(boot_stats, sep="\t", quote=FALSE, file="output/UKB_tests/QDiabA_then_C_screening_bootstraps.txt")



