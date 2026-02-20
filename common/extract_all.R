library(data.table)
library(foreach)
library(jsonlite)

# Needs at least a mem1_hdd1_v2_x72 instance

################################################################################
# Functions to help with data extractions
################################################################################

# source() a file saved in project storage on the RAP
# temporarily downloads the file so it can be run locally
rap_source <- function(path) {
  system(sprintf("dx download -f '%s'", path))
  source(basename(path))
  system(sprintf("rm %s", basename(path)))
}

# As above, but for bash scripts
rap_sh <- function(path) {
  system(sprintf("dx download -f '%s'", path))
  system(sprintf("chmod +x %s", basename(path)))
  system(sprintf("./%s", basename(path)))
  system(sprintf("rm %s", basename(path)))
}

# Function to run the curate_fields.R script in the corresponding folder and 
# launch Table Exporter to extract the relevant fields. Returns a data.t  able 
# containing the input folder name, the name of the field file (i.e. as there 
# may be multiple), the job id of the scheduled job and placeholder values for
# the job state and completion status for downstream polling.
re_extract_raw <- function(folder, ...) {
  message(sprintf("Curating fields for re-extraction of 'common/%s/'...", folder))
  # Remove old field IDs
  system(sprintf("dx mv 'common/%s/scripts/*.txt' trash/", folder))
  
  # Re-extract field IDs
  rap_source(sprintf("common/%s/scripts/curate_fields.R", folder))
  
  # Make raw_data directory
  system(sprintf("dx mkdir -p 'common/%s/raw_data'", folder))
  
  # Detect field files to extract
  fields_files <- system(sprintf("dx ls 'common/%s/scripts/*.txt'", folder), intern=TRUE)
  
  # Run Table Exporter on each field file
  message("  Launching Table Exporter...")
  foreach(ff=fields_files, .combine=rbind) %do% {
    # Detect entity type
    entity <- gsub("(.*_fields_)|(_entity.txt)", "", ff)
    
    # Build Table Exporter command
    cmd <- "dx run table-exporter"
    cmd <- paste(cmd, sprintf("-idataset_or_cohort_or_dashboard=%s", dataset_file))
    cmd <- paste(cmd, "-icoding_option=RAW")
    cmd <- paste(cmd, sprintf("-ientity=%s", entity))
    cmd <- paste(cmd, sprintf("-ifield_names_file_txt='common/%s/scripts/%s'", folder, ff))
    cmd <- paste(cmd, "--brief --yes --ignore-reuse")
    cmd <- paste(cmd, sprintf("--destination='common/%s/raw_data/'", folder))
    
    # Some need bigger instance types
    big_entities <- c("olink_instance_0", "hesin", "hesin_oper", "hesin_diag", 
      "gp_clinical", "gp_registrations", "gp_scripts")
    if (entity %in% big_entities || ff == "nmr_fields_participant_entity.txt") {
      cmd <- paste(cmd, "--instance-type='mem1_ssd1_v2_x16'")
    } else {
      cmd <- paste(cmd, "--instance-type='mem1_ssd1_v2_x8'")
    }
    
    # Some need prefixes other than 'data'
    if (entity == "olink_instance_0") cmd <- paste(cmd, "-ioutput='olink_instance_0'")
    if (entity == "olink_instance_2") cmd <- paste(cmd, "-ioutput='olink_instance_2'")
    if (entity == "olink_instance_3") cmd <- paste(cmd, "-ioutput='olink_instance_3'")
    if (ff == "proteomics_fields_participant_entity.txt") cmd <- paste(cmd, "-ioutput='sample_metadata'")
    
    prefix_as_entity <- c("death", "death_cause", "hesin", "hesin_oper", 
      "hesin_diag", "gp_clinical", "gp_registrations", "gp_scripts")
    if (entity %in% prefix_as_entity) cmd <- paste(cmd, sprintf("-ioutput='%s'", entity))
    
    # Launch Table Exporter
    jid <- system(cmd, intern=TRUE)
    
    # Return information and job id
    data.table("folder"=folder, "field_file"=ff, job_id=jid, state="", completed=FALSE)
  }
}

# Get the state of a job on DNA nexus
get_job_state <- function(jid) {
  fromJSON(system(sprintf("dx find jobs --id %s --json", jid), intern=TRUE))$state
}

# Remove curated data from a folder
remove_curated_data <- function(folder) {
  # Get old files to remove
  old_files <- system(sprintf("dx ls '%s/' --obj", folder), intern=TRUE)
  old_files <- old_files[!(old_files %like% "README")]
  
  # Remove old files
  for (ff in old_files) {
    system(sprintf("dx mv '%s/%s' trash/", folder, ff)) 
  }
}

# Function to run the curate_data.R script in the corresponding folder
curate_new <- function(folder) {
  message(sprintf("Curating data for 'common/%s/'...", folder))
  
  # Remove old files
  remove_curated_data(sprintf("common/%s", folder))
  if (folder == "NMR Metabolomics") remove_curated_data("common/NMR Metabolomics/qc_information")
  
  # Run data curation script
  rap_source(sprintf("common/%s/scripts/curate_data.R", folder))
  
  if (folder == "Medications") {
    rap_source("common/Medications/scripts/curate_medication_classes.R")
  }
}

# Function to run curate_endpoint.R on a file that 
curate_endpoint <- function(folder) {
  # Find the name of the the curated endpoint file
  out <- system(sprintf("dx ls '%s' --obj", folder), intern=TRUE)
  out <- out[out != "endpoint_definition.txt"]
  
  if (length(out) > 1) {
    stop(sprintf("Could not determine output file in '%s': too many files", folder))
  } else if (length(out) == 0) {
    message("endpoint not previous curated, setting output file to 'events_and_followup.csv'")
    out <- "events_and_followup.csv"
  }
  
  # Remove it
  system(sprintf("dx mv '%s/%s' trash/", folder, out))
  
  # Pull down curate_endpoint.R script
  if (!file.exists("curate_endpoint.R")) system("dx download common/Endpoints/curate_endpoint.R")
  
  # Run the endpoint curation script
  system(sprintf("Rscript curate_endpoint.R --def-file '%s/endpoint_definition.txt' --output '%s/%s' --verbose", 
    folder, folder, out))
}

# Function to remove a folder to trash while handling name collisions
delete_folder <- function(folder) {
  rn <- as.integer(Sys.time())
  system(sprintf("dx mv '%s' '%s_%s'", folder, folder, rn))
  system(sprintf("dx mv '%s_%s' trash/", folder, rn))
}

################################################################################
# Do the data extractions
################################################################################

# Detect dataset file for Table Exporter
dataset_file <- system("dx ls *.dataset", intern=TRUE)

# Run through all folders and run table exporter
job_info <- rbind(
  re_extract_raw("Demographics"),
  re_extract_raw("Anthropometrics"),
  re_extract_raw("Genetic Reference"),
  re_extract_raw("Blood Pressure"),
  re_extract_raw("Smoking"),
  re_extract_raw("Body Composition"),
  re_extract_raw("Biomarkers"),
  re_extract_raw("Family History"),
  re_extract_raw("Medical History"),
  re_extract_raw("Medications"),
  re_extract_raw("Proteomics"),
  re_extract_raw("NMR Metabolomics"),
  re_extract_raw("Deaths"),
  re_extract_raw("Hospital Records"),
  re_extract_raw("Cancer Register"),
  re_extract_raw("Primary Care"),
  re_extract_raw("Follow-up"),
  re_extract_raw("Endpoints/UKB curated")
)

# Monitor jobs until finished, poll every minute, 1 hour timeout
timeout <- 60
counter <- 0
while (any(!(job_info$completed)) || counter > timeout) {
  Sys.sleep(60)
  counter <- counter + 1
  print(sprintf("Checking Table Exporter job status... (%s/%s)", counter, timeout))
  
  job_info[!(completed), state := get_job_state(job_id), by=.(folder, field_file)]
  job_info[state %in% c("done", "terminated", "failed"), completed := TRUE]
} 

# Make sure that they all finished successfully before continuing
if (any(!(job_info$completed))) {
  fwrite(job_info[!(completed)], file="incomplete_extraction_jobs.csv")
  system("dx upload incomplete_extraction_jobs.csv --destination common/")
  stop("Some jobs did not complete! See common/incomplete_extraction_jobs.csv")
}

# Run the data curation scripts - order has been carefully set here, some of
# these depend on earlier curations in the sequence!
curate_new("Demographics")
curate_new("Anthropometrics")
curate_new("Genetic Reference")
curate_new("Blood Pressure")
curate_new("Smoking")
curate_new("Body Composition")
curate_new("Biomarkers")
curate_new("Family History")
curate_new("Medical History")
curate_new("Medications")
curate_new("Proteomics")
curate_new("NMR Metabolomics")
curate_new("Deaths")
curate_new("Hospital Records")
curate_new("Cancer Register")
curate_new("Primary Care")
curate_new("Follow-up")
curate_new("Endpoints/UKB curated")

# Re-run Eastwood et al. 2016 diabetes adjudication algorithms
remove_curated_data("common/Endpoints/Diabetes (Eastwood Algorithm)")
rap_source("common/Endpoints/Diabetes (Eastwood Algorithm)/scripts/01_adjudicate_prevalent_diabetes.R")
rap_source("common/Endpoints/Diabetes (Eastwood Algorithm)/scripts/02_adjudicate_incident_diabetes.R")

# Re-extract all curated endpoints
cmd <- "dx find data --path 'common/Endpoints/' --name endpoint_definition.txt --json"
endpoint_folders <- fromJSON(system(cmd, intern=TRUE))$describe$folder
for (dd in endpoint_folders) curate_endpoint(dd)

# Re-extract phecodes
# Warning, this launches ~3000 jobs that will take 68 compute days costing ~£420
# There is also a hard limit on 100 concurrent jobs per user, so you will also 
# be locked out of running any interactive sessions or new jobs until these all
# complete (takes about a day).
stop("User attention check! the next lines of code are compute heavy and expensive (see code comment)")
delete_folder("common/Endpoints/Phecodes/PhecodeX")
rap_source('common/Endpoints/Phecodes/scripts/PhecodeX/setup_phecode_batches.R')
rap_sh('common/Endpoints/Phecodes/scripts/PhecodeX/launch_phecode_extraction_jobs.sh') #
