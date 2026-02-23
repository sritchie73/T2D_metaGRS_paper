# Install missing packages on the RAP
remotes::install_github("sritchie73/dxutils")

# Load libraries
library(data.table)
library(foreach)
library(dxutils)

# Download and load in the phenotype data
dx_download("data/ukb_t2d_pheno.tsv", "input_data/")
pheno <- fread("input_data/ukb_t2d_pheno.tsv")
pheno <- pheno[visit_index == 0]

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

# Download and add information about bmi
dx_download("common/Anthropometrics/anthropometrics.csv", "input_data/")
anthro <- fread("input_data/anthropometrics.csv")
pheno[anthro[visit_index == 0], on = .(eid), bmi := i.bmi]

# Download and add information about smoking status
dx_download("common/Smoking/smoking.csv", "input_data/")
smoking <- fread("input_data/smoking.csv")
pheno[smoking[visit_index == 0], on = .(eid), smoker := i.current_smoker]

# Compute cohort characteristics in each genetic ancestry
ancestries <- pheno[,.N,by=genetic_ancestry][order(-N)]
cc <- foreach(this_ancestry = ancestries$genetic_ancestry, .combine=rbind) %do% {
  this_pheno <- pheno[genetic_ancestry == this_ancestry]
  data.table(
    "genetic_ancestry"=this_ancestry,
    "sample_size_cases"=sprintf("%s (%s/%s)",
      format(this_pheno[(t2d_case), .N], big.mark=","),
      format(this_pheno[(t2d_case) & sex == "Male", .N], big.mark=","),
      format(this_pheno[(t2d_case) & sex == "Female", .N], big.mark=",")
    ),
    "sample_size_controls"=sprintf("%s (%s/%s)",
      format(this_pheno[!(t2d_case), .N], big.mark=","),
      format(this_pheno[!(t2d_case) & sex == "Male", .N], big.mark=","),
      format(this_pheno[!(t2d_case) & sex == "Female", .N], big.mark=",")
    ),
    "age_cases"=sprintf("%s (%s)",
      this_pheno[(t2d_case), round(mean(age), digits=1)],
      this_pheno[(t2d_case), round(sd(age), digits=1)]
    ),
    "age_controls"=sprintf("%s (%s)",
      this_pheno[!(t2d_case), round(mean(age), digits=1)],
      this_pheno[!(t2d_case), round(sd(age), digits=1)]
    ),
    "bmi_cases"=sprintf("%s (%s)",
      this_pheno[(t2d_case) & !is.na(bmi), round(mean(bmi), digits=1)],
      this_pheno[(t2d_case) & !is.na(bmi), round(sd(bmi), digits=1)]
    ),
    "bmi_controls"=sprintf("%s (%s)",
      this_pheno[!(t2d_case) & !is.na(bmi), round(mean(bmi), digits=1)],
      this_pheno[!(t2d_case) & !is.na(bmi), round(sd(bmi), digits=1)]
    ),
    "smokers_cases"=sprintf("%s (%s/%s)",
      format(this_pheno[(t2d_case) & (smoker), .N], big.mark=","),
      format(this_pheno[(t2d_case) & (smoker) & sex == "Male", .N], big.mark=","),
      format(this_pheno[(t2d_case) & (smoker) & sex == "Female", .N], big.mark=",")
    ),
    "smokers_controls"=sprintf("%s (%s/%s)",
      format(this_pheno[!(t2d_case) & (smoker), .N], big.mark=","),
      format(this_pheno[!(t2d_case) & (smoker) & sex == "Male", .N], big.mark=","),
      format(this_pheno[!(t2d_case) & (smoker) & sex == "Female", .N], big.mark=",")
    )
  )
}
fwrite(cc, sep="\t", quote=FALSE, file="case_control_cohort_characteristics.txt")
dx_upload("case_control_cohort_characteristics.txt", "output/")

# Filter to samples we would use for incident T2D analyses
pheno <- pheno[!is.na(inci_t2d_followup)]

# Truncate follow-up to 10-years
pheno[!(t2d_case), inci_t2d_followup := pmin(10, inci_t2d_followup)]
pheno[(t2d_case) & inci_t2d_followup > 10, c("t2d_case", "inci_t2d_followup") := .(FALSE, 10)]

# Compute cohort characteristics
ancestries <- pheno[,.N,by=genetic_ancestry][order(-N)]
cc <- foreach(this_ancestry = ancestries$genetic_ancestry, .combine=rbind) %do% {
  this_pheno <- pheno[genetic_ancestry == this_ancestry]
  data.table(
    "genetic_ancestry"=this_ancestry,
    "sample_size_cases"=sprintf("%s (%s/%s)",
      format(this_pheno[(t2d_case), .N], big.mark=","),
      format(this_pheno[(t2d_case) & sex == "Male", .N], big.mark=","),
      format(this_pheno[(t2d_case) & sex == "Female", .N], big.mark=",")
    ),
    "sample_size_controls"=sprintf("%s (%s/%s)",
      format(this_pheno[!(t2d_case), .N], big.mark=","),
      format(this_pheno[!(t2d_case) & sex == "Male", .N], big.mark=","),
      format(this_pheno[!(t2d_case) & sex == "Female", .N], big.mark=",")
    ),
    "followup_cases"=sprintf("%s (%s)",
      this_pheno[(t2d_case), round(mean(inci_t2d_followup), digits=1)],
      this_pheno[(t2d_case), round(sd(inci_t2d_followup), digits=1)]
    ),
    "followup_controls"=sprintf("%s (%s)",
      this_pheno[!(t2d_case), round(mean(inci_t2d_followup), digits=1)],
      this_pheno[!(t2d_case), round(sd(inci_t2d_followup), digits=1)]
    ),
    "age_cases"=sprintf("%s (%s)",
      this_pheno[(t2d_case), round(mean(age), digits=1)],
      this_pheno[(t2d_case), round(sd(age), digits=1)]
    ),
    "age_controls"=sprintf("%s (%s)",
      this_pheno[!(t2d_case), round(mean(age), digits=1)],
      this_pheno[!(t2d_case), round(sd(age), digits=1)]
    ),
    "bmi_cases"=sprintf("%s (%s)",
      this_pheno[(t2d_case) & !is.na(bmi), round(mean(bmi), digits=1)],
      this_pheno[(t2d_case) & !is.na(bmi), round(sd(bmi), digits=1)]
    ),
    "bmi_controls"=sprintf("%s (%s)",
      this_pheno[!(t2d_case) & !is.na(bmi), round(mean(bmi), digits=1)],
      this_pheno[!(t2d_case) & !is.na(bmi), round(sd(bmi), digits=1)]
    ),
    "smokers_cases"=sprintf("%s (%s/%s)",
      format(this_pheno[(t2d_case) & (smoker), .N], big.mark=","),
      format(this_pheno[(t2d_case) & (smoker) & sex == "Male", .N], big.mark=","),
      format(this_pheno[(t2d_case) & (smoker) & sex == "Female", .N], big.mark=",")
    ),
    "smokers_controls"=sprintf("%s (%s/%s)",
      format(this_pheno[!(t2d_case) & (smoker), .N], big.mark=","),
      format(this_pheno[!(t2d_case) & (smoker) & sex == "Male", .N], big.mark=","),
      format(this_pheno[!(t2d_case) & (smoker) & sex == "Female", .N], big.mark=",")
    )
  )
}
fwrite(cc, sep="\t", quote=FALSE, file="incident_t2d_cohort_characteristics.txt")
dx_upload("incident_t2d_cohort_characteristics.txt", "output/")
