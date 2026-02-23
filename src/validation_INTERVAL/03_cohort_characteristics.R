library(data.table)

# Load pre-curated phenotype data
pheno <- fread("data/INTERVAL/curated_pheno.tsv")

# Drop people with prevalent diabetes (mix of T1D and T2D)
pheno <- pheno[!(prevalent_diabetes)]

# Load in extended phneotype information so we can get smoking status and BMI
all_pheno <- fread("data/INTERVAL/P1074/INTERVALdata_02FEB2024.csv")

# Compute BMI and remove outliers
all_pheno[, bmi := wt_bl/ht_bl^2]
all_pheno[wt_bl == 777, bmi := NA_real_] # bad coding
all_pheno[ht_bl < 1.47, bmi := NA_real_] # clinical cutoff for dwarfism
all_pheno[ht_bl > 2.1, bmi := NA_real_] # clinical cutoff for gigantism
all_pheno[wt_bl < 50 | wt_bl > 160, bmi := NA_real_] # NHS restrictions for weight

# Add BMI to the main table
pheno[all_pheno, on = .(identifier), bmi := i.bmi]

# Curate smoking status
pheno[all_pheno, on = .(identifier), smoker := fcase(
  is.na(smCurr_bl), FALSE, # Answered "No" to ever smoked
  smCurr_bl == 1, TRUE,
  smCurr_bl == 2, FALSE,
  smCurr_bl == 999, NA  # Do not know/prefer not to answer
)]

# Curate sex
pheno[, sex := ifelse(sex == 2, "Female", "Male")]

# Compute cohort characteristics
cc <- data.table(
  "sample_size_cases"=sprintf("%s (%s/%s)",
    format(pheno[(incident_diabetes), .N], big.mark=","),
    format(pheno[(incident_diabetes) & sex == "Male", .N], big.mark=","),
    format(pheno[(incident_diabetes) & sex == "Female", .N], big.mark=",")
  ),
  "sample_size_controls"=sprintf("%s (%s/%s)",
    format(pheno[!(incident_diabetes), .N], big.mark=","),
    format(pheno[!(incident_diabetes) & sex == "Male", .N], big.mark=","),
    format(pheno[!(incident_diabetes) & sex == "Female", .N], big.mark=",")
  ),
  "age_cases"=sprintf("%s (%s)",
    pheno[(incident_diabetes), round(mean(age), digits=1)],
    pheno[(incident_diabetes), round(sd(age), digits=1)]
  ),
  "age_controls"=sprintf("%s (%s)",
    pheno[!(incident_diabetes), round(mean(age), digits=1)],
    pheno[!(incident_diabetes), round(sd(age), digits=1)]
  ),
  "bmi_cases"=sprintf("%s (%s)",
    pheno[(incident_diabetes) & !is.na(bmi), round(mean(bmi), digits=1)],
    pheno[(incident_diabetes) & !is.na(bmi), round(sd(bmi), digits=1)]
  ),
  "bmi_controls"=sprintf("%s (%s)",
    pheno[!(incident_diabetes) & !is.na(bmi), round(mean(bmi), digits=1)],
    pheno[!(incident_diabetes) & !is.na(bmi), round(sd(bmi), digits=1)]
  ),
  "smokers_cases"=sprintf("%s (%s/%s)",
    format(pheno[(incident_diabetes) & (smoker), .N], big.mark=","),
    format(pheno[(incident_diabetes) & (smoker) & sex == "Male", .N], big.mark=","),
    format(pheno[(incident_diabetes) & (smoker) & sex == "Female", .N], big.mark=",")
  ),
  "smokers_controls"=sprintf("%s (%s/%s)",
    format(pheno[!(incident_diabetes) & (smoker), .N], big.mark=","),
    format(pheno[!(incident_diabetes) & (smoker) & sex == "Male", .N], big.mark=","),
    format(pheno[!(incident_diabetes) & (smoker) & sex == "Female", .N], big.mark=",")
  )
)
fwrite(cc, sep="\t", quote=FALSE, file="output/INTERVAL_tests/case_control_cohort_characteristics.txt")

# Compute cohort characteristics for 10 year T2D risk prediction
# Truncate follow-up to 10-years
pheno[!(incident_diabetes), incident_censor_years := pmin(10, incident_censor_years)]
pheno[(incident_diabetes) & incident_censor_years > 10, c("incident_diabetes", "incident_censor_years") := .(FALSE, 10)]
cc <- data.table(
  "sample_size_cases"=sprintf("%s (%s/%s)",
    format(pheno[(incident_diabetes), .N], big.mark=","),
    format(pheno[(incident_diabetes) & sex == "Male", .N], big.mark=","),
    format(pheno[(incident_diabetes) & sex == "Female", .N], big.mark=",")
  ),
  "sample_size_controls"=sprintf("%s (%s/%s)",
    format(pheno[!(incident_diabetes), .N], big.mark=","),
    format(pheno[!(incident_diabetes) & sex == "Male", .N], big.mark=","),
    format(pheno[!(incident_diabetes) & sex == "Female", .N], big.mark=",")
  ),
  "followup_cases"=sprintf("%s (%s)",
    pheno[(incident_diabetes), round(mean(incident_censor_years), digits=1)],
    pheno[(incident_diabetes), round(sd(incident_censor_years), digits=1)]
  ),
  "followup_controls"=sprintf("%s (%s)",
    pheno[!(incident_diabetes), round(mean(incident_censor_years), digits=1)],
    pheno[!(incident_diabetes), round(sd(incident_censor_years), digits=1)]
  ),
  "age_cases"=sprintf("%s (%s)",
    pheno[(incident_diabetes), round(mean(age), digits=1)],
    pheno[(incident_diabetes), round(sd(age), digits=1)]
  ),
  "age_controls"=sprintf("%s (%s)",
    pheno[!(incident_diabetes), round(mean(age), digits=1)],
    pheno[!(incident_diabetes), round(sd(age), digits=1)]
  ),
  "bmi_cases"=sprintf("%s (%s)",
    pheno[(incident_diabetes) & !is.na(bmi), round(mean(bmi), digits=1)],
    pheno[(incident_diabetes) & !is.na(bmi), round(sd(bmi), digits=1)]
  ),
  "bmi_controls"=sprintf("%s (%s)",
    pheno[!(incident_diabetes) & !is.na(bmi), round(mean(bmi), digits=1)],
    pheno[!(incident_diabetes) & !is.na(bmi), round(sd(bmi), digits=1)]
  ),
  "smokers_cases"=sprintf("%s (%s/%s)",
    format(pheno[(incident_diabetes) & (smoker), .N], big.mark=","),
    format(pheno[(incident_diabetes) & (smoker) & sex == "Male", .N], big.mark=","),
    format(pheno[(incident_diabetes) & (smoker) & sex == "Female", .N], big.mark=",")
  ),
  "smokers_controls"=sprintf("%s (%s/%s)",
    format(pheno[!(incident_diabetes) & (smoker), .N], big.mark=","),
    format(pheno[!(incident_diabetes) & (smoker) & sex == "Male", .N], big.mark=","),
    format(pheno[!(incident_diabetes) & (smoker) & sex == "Female", .N], big.mark=",")
  )
)
fwrite(cc, sep="\t", quote=FALSE, file="output/INTERVAL_tests/incident_t2d_cohort_characteristics.txt")
