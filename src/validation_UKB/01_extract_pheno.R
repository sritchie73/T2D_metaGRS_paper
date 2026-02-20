library(data.table)
library(dxutils) # remotes::install_github("sritchie73/dxutils")

# Pull down supplemental functions
dx_download("src/functions/calendar_time.R", "src/functions/")
source("src/functions/calendar_time.R")

# Pull down basic participant demographic data and filter to columns of interest
dx_download("common/Demographics/demographics.csv", "input_data/")
pheno <- fread("input_data/demographics.csv")
pheno <- pheno[,.(eid, visit_index, assessment_date, assessment_centre, age=age_decimal, sex)]

# Add in information about genetic ancestry
dx_download("common/Genetic Reference/ancestry.csv", "input_data/")
pheno[genetic_ancestry, on = .(eid), genetic_ancestry := ancestry]

# Pull down the different sources of T2D cases - those identified through GP
# records (available for ~1/3rd of participants) - along with those identified
# from self-reported medical history, medication use, and as contributing causes
# to hospitalisations or deaths (adjudicated with Eastwood et al. 2016 
# algorithms)
dx_download("common/Endpoints/T2D_GP/events_and_followup.csv", "input_data/primary_care_t2d.csv")
dx_download("common/Endpoints/Diabetes (Eastwood Algorithm)/prevalent_diabetes.csv", "input_data/eastwood_prevalent_t2d.csv")
dx_download("common/Endpoints/Diabetes (Eastwood Algorithm)/incident_diabetes.csv", "input_data/eastwood_incident_t2d.csv")

# Pull down the information about maximum and minimum follow-up times in each data source
dx_download("common/Follow-up/follow_up.csv", "input_data/")

# Build combined table of all possible t2d case sources
t2d <- pheno[,.(eid, visit_index, assessment_date, age)]

eastwood_prev <- fread("input_data/eastwood_prevalent_t2d.csv")
t2d[eastwood_prev, on = .(eid, visit_index), 
    c("eastwood_prev_adjudicated_diabetes", "eastwood_prev_followup") :=
    .(adjudicated_diabetes, age - age_of_onset)]

eastwood_inci <- fread("input_data/eastwood_incident_t2d.csv")
t2d[eastwood_inci, on = .(eid, visit_index), 
  c("eastwood_inci_adjudicated_diabetes", "eastwood_inci_followup", "eastwood_inci_adjudication_note") :=
  .(adjudicated_diabetes, interpolated_follow_years, adjudication_note)]

t2d_gp <- fread("input_data/primary_care_t2d.csv")
t2d[t2d_gp[reason_or_id != "Assessment after maximum censor date"], on = .(eid, visit_index),
  c("gp_t2d", "gp_followup") := .(event, follow_up)]

follow <- fread("input_data/follow_up.csv")
t2d[follow, on = .(eid), c("min_gp_followup", "max_gp_followup", "max_hospital_followup") := 
  .(years_between(assessment_date, min_gp_censor), years_between(assessment_date, max_gp_censor),
    years_between(assessment_date, max_hospital_censor))]

# Get combined column of all possible T2D cases
t2d[, t2d_case := 
  eastwood_prev_adjudicated_diabetes %like% "type 2 diabetes" |
  eastwood_inci_adjudicated_diabetes %like% "type 2 diabetes" |
  (
    eastwood_inci_adjudicated_diabetes == "Prevalent diabetes" & 
    eastwood_inci_adjudication_note %like% "type 2 diabetes records" & 
    !(eastwood_inci_adjudication_note %like% " 0 type 2 diabetes records")
  ) |
  (!is.na(gp_t2d) & (gp_t2d))
]

# For time to incident diabetes, preferentially use the first date in the GP
# records, if that doesn't exist, use the midpoint date inferred from the HES
# records
t2d[(gp_t2d) & gp_followup > 0, inci_t2d_followup := gp_followup]
t2d[!(gp_t2d) | is.na(gp_t2d), inci_t2d_followup := eastwood_inci_followup]
t2d[is.na(inci_t2d_followup) & !(t2d_case) & eastwood_inci_adjudicated_diabetes == "No evidence of diabetes" & max_hospital_followup > 0, inci_t2d_followup := max_hospital_followup]

# Add in to main phenotype table
# NA in inci_t2d_followup column indicates: (1) prevalent non-T2D diabetes, 
# (2) withdrawn consent for EHR linkage, or (3) maximum follow-up in hospital 
# records is prior to assessment visit
pheno[t2d, on = .(eid, visit_index), c("t2d_case", "inci_t2d_followup") := .(i.t2d_case, i.inci_t2d_followup)]

# Save and upload
fwrite(pheno, sep="\t", quote=FALSE, file="ukb_t2d_pheno.tsv")
dx_upload("ukb_t2d_pheno.tsv")
