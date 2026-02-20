library(data.table)
library(foreach)

# Check if we have the necessary python libraries installed for dx extract_dataset
if (system("python3 -c 'import pandas'", ignore.stderr=TRUE)) {
  system("pip install pandas", wait=TRUE)
}

# Extract dataset meta-data to local cloud instance, if it hasn't already been downloaded
if (length(list.files(pattern="data_dictionary")) == 0) {
  project_files <- system("dx ls", intern=TRUE)
  dataset_file <- project_files[project_files %like% "^app" & project_files %like% ".dataset$"]
  system(sprintf("dx extract_dataset %s -ddd", dataset_file), wait=TRUE) # Can ignore warning about pandas version
}

# Auto-detect data-dictionary file then load
data_dict_file <- list.files(pattern="data_dictionary")
data_dict <- fread(data_dict_file)

# Curate biomarker info
biomarker_info <- rbind(use.names=FALSE,
  data.table(var = "tchol", biomarker = "Cholesterol", units="mmol/l", sample_type = "Serum", supplier = "Beckman Coulter",
             instrumentation = "AU5800", analysis_method = "CHO-POD", group = "Cardiovascular"),
  data.table("ldl", "Direct Low Density Lipoprotein", "mmol/l", "Serum", "Beckman Coulter", "AU5800", "Enzymatic selective protection", "Cardiovascular"),
  data.table("hdl", "HDL-Cholesterol", "mmol/l", "Serum", "Beckman Coulter", "AU5800", "Enzyme immunoinhibition", "Cardiovascular"),
  data.table("trig", "Triglyceride", "mmol/l", "Serum", "Beckman Coulter", "AU5800", "GPO-POD", "Cardiovascular"),
  data.table("apoa1", "Apolipoprotein A", "g/l", "Serum", "Beckman Coulter", "AU5800", "Immunoturbidimetric", "Cardiovascular"),
  data.table("apob", "Apolipoprotein B", "g/l", "Serum", "Beckman Coulter", "AU5800", "Immunoturbidimetric", "Cardiovascular"),
  data.table("crp", "C-reactive Protein", "mg/l", "Serum", "Beckman Coulter", "AU5800", "Immunoturbidimetric - high sensitivity", "Cardiovascular"),
  data.table("lpa", "Lipoprotein (a)", "mg/dl", "Serum", "Randox", "AU5800", "Immunoturbidimetric", "Cardiovascular"),
  data.table("vitd25", "Vitamin D", "nmol/l", "Serum", "DiaSorin Ltd.", "LIASON XL", "CLIA", "Bone and joint"),
  data.table("rheuf", "Rheumatoid factor", "iu/ml", "Serum", "Beckman Coulter", "AU5800", "Immunoturbidimetric", "Bone and joint"),
  data.table("alp", "Alkaline Phosphatase", "iu/l", "Serum", "Beckman Coulter", "AU5800", "AMP (IFCC)", "Bone and joint"),
  data.table("calcium", "Calcium", "mmol/l", "Serum", "Beckman Coulter", "AU5800", "Arsenazo III", "Bone and joint"),
  data.table("shbg", "SHBG", "nmol/l", "Serum", "Beckman Coulter", "Unicel DxI 800", "Two step standwich immunoassay", "Cancer"),
  data.table("testos", "Testosterone", "nmol/l", "Serum", "Beckman Coulter", "Unicel DxI 800", "One step competitive", "Cancer"),
  data.table("oest", "Oestradiol", "pmol/l", "Serum", "Beckman Coulter", "Unicel DxI 800", "Two step competitive", "Cancer"),
  data.table("igf1", "IGF-1", "nmol/l", "Serum", "DiaSorin Ltd.", "LIASON XL", "CLIA", "Cancer"),
  data.table("hba1c", "HbA1c", "mmol/mol", "RBC", "Bio-Rad", "VARIANT II Turbo", "HPLC", "Diabetes"),
  data.table("glucose", "Glucose", "mmol/l", "Serum", "Beckman Coulter", "AU5800", "Hexokinase", "Diabetes"),
  data.table("cyst", "Cystatin C", "mg/l", "Serum", "Siemens", "ADVIA 1800", "Latex enhanced immunoturbidimetric", "Renal"),
  data.table("creat", "Creatinine", "µmol/l", "Serum", "Beckman Coulter", "AU5800", "Enzymatic",  "Renal"),
  data.table("protein", "Total protein", "g/l", "Serum", "Beckman Coulter", "AU5800", "Biuret", "Renal"),
  data.table("urea", "Urea", "mmol/l", "Serum", "Beckman Coulter", "AU5800", "GLDH, kinetic", "Renal"),
  data.table("phos", "Phosphate", "mmol/l", "Serum", "Beckman Coulter", "AU5800", "Phophomolybdate complex", "Renal"),
  data.table("uric", "Urate", "µmol/l", "Serum", "Beckman Coulter", "AU5800", "Uricase PAP", "Renal"),
  data.table("uriacc", "Creatinine (enzymatic)", "µmol/l", "Urine", "Beckman Coulter", "AU5400", "Enzymatic", "Renal"),
  data.table("urianac", "Sodium", "mmol/l", "Urine", "Beckman Coulter", "AU5400", "ISE", "Renal"),
  data.table("uriamac", "Microalbumin", "mg/l", "Urine", "Randox", "AU5400", "Immunoturbidimetric", "Renal"),
  data.table("uriakc", "Potassium", "mmol/l", "Urine", "Beckman Coulter", "AU5400", "ISE", "Renal"),
  data.table("alb", "Albumin", "g/l", "Serum", "Beckman Coulter", "AU5800", "BCG", "Liver"),
  data.table("dbili", "Direct Bilirubin", "µmol/l", "Serum", "Beckman Coulter", "AU5800", "DPD", "Liver"),
  data.table("tbili", "Total Bilirubin", "µmol/l", "Serum", "Beckman Coulter", "AU5800", "Photometric colour", "Liver"),
  data.table("ggt", "Gamma Glutamyltransferase", "iu/l", "Serum", "Beckman Coulter", "AU5800", "IFCC", "Liver"),
  data.table("alt", "Alanine aminotransferase", "iu/l", "Serum", "Beckman Coulter", "AU5800", "IFCC", "Liver"),
  data.table("asp", "Aspartate aminotransferase", "iu/l", "Serum", "Beckman Coulter", "AU5800", "IFCC", "Liver"),
  data.table("nonhdl", "Non HDL Cholesterol", "mmol/l", "Serum", "Beckman Coulter", "AU5800", "derived: tchol - hdl", "Cardiovascular"),
  data.table("apobapoa1", "Apolipoprotein B / Apolipoprotein A", "ratio", "Serum", "Beckman Coulter", "AU5800", "derived: apob / apoa1", "Cardiovascular"),
  data.table("hba1c_pct", "HbA1c", "%", "RBC", "Bio-Rad", "VARIANT II Turbo", "derived: hba1c / 10.929 + 2.15", "Diabetes")
)

# Add UKB Field IDs
biomarker_info[, UKB.Field.ID := fcase(
  var == "alt",  30620,
  var == "alb", 30600,
  var == "alp", 30610,
  var == "apoa1", 30630,
  var == "apob", 30640,
  var == "asp", 30650,
  var == "crp", 30710,
  var == "calcium", 30680,
  var == "tchol", 30690,
  var == "creat", 30700,
  var == "cyst", 30720,
  var == "dbili", 30660,
  var == "ggt", 30730,
  var == "glucose", 30740,
  var == "hba1c", 30750,
  var == "hdl", 30760,
  var == "igf1", 30770,
  var == "ldl", 30780,
  var == "lpa", 30790,
  var == "oest", 30800,
  var == "phos", 30810,
  var == "rheuf", 30820,
  var == "shbg", 30830,
  var == "testos", 30850,
  var == "tbili", 30840,
  var == "protein", 30860,
  var == "trig", 30870,
  var == "uric", 30880,
  var == "urea", 30670,
  var == "vitd25", 30890,
  var == "uriacc", 30510,
  var == "uriamac", 30500,
  var == "uriakc", 30520,
  var == "urianac", 30530
)]

# Add Field IDs for biomarker missigness
biomarker_info[, Reportability.Field.ID := fcase(
  var == "alt", 30626,
  var == "alb", 30606,
  var == "alp", 30616,
  var == "apoa1", 30636,
  var == "apob", 30646,
  var == "asp", 30656,
  var == "crp", 30716,
  var == "calcium", 30686,
  var == "tchol", 30696,
  var == "creat", 30706,
  var == "cyst", 30726,
  var == "dbili", 30666,
  var == "ggt", 30736,
  var == "glucose", 30746,
  var == "hba1c", 30756,
  var == "hdl", 30766,
  var == "igf1", 30776,
  var == "ldl", 30786,
  var == "lpa", 30796,
  var == "oest", 30806,
  var == "phos", 30816,
  var == "rheuf", 30826,
  var == "shbg", 30836,
  var == "testos", 30856,
  var == "tbili", 30846,
  var == "protein", 30866,
  var == "trig", 30876,
  var == "uric", 30886,
  var == "urea", 30676,
  var == "vitd25", 30896,
  var == "uriacc", 30515,
  var == "uriamac", 30505,
  var == "uriakc", 30525,
  var == "urianac", 30535
)]

# Curate sample QC information also
sample_info <- rbind(use.names=FALSE,
  data.table(var = "eid", description = "Participant identifier for application 7439 (adiposity project)", UKB.Field.ID = NA_integer_),
  data.table("visit_index", "Timepoint/visit for sample measurement of that participant. 0 = baseline assessment, 1 = first repeat assessment", NA),
  data.table("no_blood_sample_reason", "Reason blood sampling not attempted", 20050),
  data.table("no_urine_sample_reason", "Reason urine sampling not attempted", 20072),
  data.table("no_blood_sample", "TRUE where \"Reason blood sampling not attempted\" field (UKB field ID: 20050) is not empty, FALSE otherwise", NA),
  data.table("no_urine_sample", "TRUE where \"Reason urine sampling not attempted\" field (UKB field ID: 20072) is not empty, FALSE otherwise", NA),
  data.table("fasting_time", "Time since last food/drink (hours)", 74)
)

# Save table of information
system("mkdir -p biomarkers")
fwrite(biomarker_info, file="biomarkers/biomarker_information.csv")
fwrite(sample_info, file="biomarkers/sample_information.csv")

# Ge set of field IDs to extract
fields <- c(
  na.omit(biomarker_info$UKB.Field.ID),
  na.omit(biomarker_info$Reportability.Field.ID),
  na.omit(sample_info$UKB.Field.ID)
)

# Extract list of fields for Table Exporter
fields <- foreach(this_field_id = fields, .combine=rbind) %do% {
  data_dict[name %like% sprintf("^p%s(_|$)", this_field_id)]
}

# Add in participant ID
fields <- rbind(data_dict[entity == "participant" & name == "eid"], fields)

# Write out field IDs for each distinct entity (only participant in this case,
# but generalizes later if we want to extract different entity types)
system("mkdir -p field_lists/")
for (entity_type in unique(fields$entity)) {
  fwrite(fields[entity == entity_type, .(name)], quote=FALSE, col.names=FALSE,
    file=sprintf("field_lists/biomarkers_fields_%s_entity.txt", entity_type))
}

# Upload field list to persistent storage
system("dx upload field_lists/biomarkers_fields_* --destination 'common/Biomarkers/scripts/'", wait=TRUE)