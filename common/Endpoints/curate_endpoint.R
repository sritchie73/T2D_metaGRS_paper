#!/usr/bin/env Rscript
########################################################################################
# Load R package dependencies
########################################################################################
options(repos=c(CRAN="https://cloud.r-project.org/"))

# Load libraries, install if missing
load_pkg <- function(pkg_name) {
  if (suppressMessages(suppressWarnings(!require(pkg_name, character.only=TRUE)))) {
    install.packages(pkg_name)
    suppressMessages(library(pkg_name, character.only=TRUE))
  }
}

load_pkg("docopt")
load_pkg("data.table")
load_pkg("foreach")
load_pkg("bit64")
load_pkg("lubridate")

########################################################################################
# Custom function definitions
########################################################################################

# Function that only prints when args[["--verbose"]] is TRUE
vmessage <- function(...) {
  if (args[["--verbose"]]) message(...)
}

# Function for converting ICD and OPCS codes and code ranges to individual codes
# in UKB format. E.g.:
# 
#  I24.2 becomes I242
#  I20-I25 becomes I20, I21, I22, I23, I24, I25
#
parse_codelist <- function(x) {
  if (length(x) == 0) return(character(0))
  codes <- strsplit(split=", ?", x)[[1]]
  chapters <- codes[grepl("-", codes)]
  codes <- codes[!grepl("-", codes)]
  if (length(chapters) > 0) {
    chapters <- strsplit(split=" ?- ?", chapters)
    chapters <- data.table(start=sapply(chapters, `[`, 1), end=sapply(chapters, `[`, 2))
    chapters[, start := gsub("\\.", "", start)]
    chapters[, end := gsub("\\.", "", end)]
    chapters[, chapter_letter := substr(start, 0, 1)]
    chapters[!(tolower(chapter_letter) %in% letters), chapter_letter := ""]
    chapters[, start := gsub("[a-z]", "", tolower(start))]
    chapters[, end := gsub("[a-z]", "", tolower(end))]
    chapters[, group := .I]
    chapters[nchar(start) == 1, start := paste0(start, "0")]
    chapters[nchar(end) == 1, end := paste0(end, "0")]
    chapter_codes <- chapters[, .(codes=paste0(chapter_letter, sprintf("%02d", seq(start, end)))), by = group][, codes]
    return(sort(c(codes, chapter_codes)))
  } else {
    return(codes)
  }
}

# Format an ICD/OPCS code for regex matching to UKB records
format_code <- function(code) {
  code <- paste0("^", gsub("\\.", "", code))
}

# Function for parsing read code-value pairs and self-report field-code pairs
parse_code_value_pairs <- function(x) {
  if (length(x) == 0) return(character(0))
  codes <- strsplit(split=", ?", x)[[1]]
  key_value_pairs <- strsplit(split="( ?= ?)|( ?!= ?)|( ?> ?)|( ?>= ?)|( ?< ?)|( ?<= ?)", codes)
  operators <- fcase(
    codes %like% ">=", ">=",
    codes %like% "<=", "<=",
    codes %like% "!=", "!=",
    codes %like% "<", "<",
    codes %like% ">", ">",
    codes %like% "=", "=",
    default=NA
  )
  list(
    "codes"=sapply(key_value_pairs, `[`, 1),
    "value_comparison_operators"=operators,
    "values"=sapply(key_value_pairs, `[`, 2)
  )
}

# Function to load a dataset from the RAP
load_from_rap <- function(rap_path) {
  if (args[["--container-mode"]]) {
    fread(file=sprintf("/mnt/project/%s", rap_path), na.strings=c("", "NA"))
  } else {
    if (!dir.exists(sprintf("%s/input_data/", work_dir))) {
      system(sprintf("mkdir -p %s/input_data", work_dir), ignore.stdout=!args[["--verbose"]])
    }
    fname <- basename(rap_path) 
    if (!file.exists(sprintf("%s/input_data/%s", work_dir, fname))) {
      system(sprintf("dx download '%s' -o '%s/input_data/%s'", rap_path, work_dir, fname), ignore.stdout=!args[["--verbose"]])
    }
    fread(file=sprintf("%s/input_data/%s", work_dir, fname), na.strings=c("", "NA"))
  }
}

# Function to compute years between two dates
years_between <- function(d1, d2) {
  as.period(interval(as.Date(d1), as.Date(d2)), unit="years") / years(1)
}

# Get a date after adding a specific number of whole years
add_years <- function(d1, follow) {
  as.IDate(as.Date(d1) + years(follow))
}

# Add a day to a date - used for offsetting min/max censor dates so we can always
# use the < or > operator when finding events rather than needing a <= or >= 
# comparison in the first matching instance.
add_day <- function(d1) {
  as.IDate(as.Date(d1) + days(1))
}

# Get the midpoint between two dates
midpoint <- function(d1, d2) {
  m <- as.IDate(length(d1))
  not_na <- !is.na(d1) & !is.na(d2)
  m[not_na] <- as.IDate(date_decimal((decimal_date(d2[not_na]) - decimal_date(d1[not_na]))/2 + decimal_date(d1[not_na])))
  return(m)
}

########################################################################################
# Set up program options
########################################################################################

"Curate incidence and/or prevalence for a custom endpoint

Given a set of ICD codes (primary care, hospitalisations, and deaths), OPCS 
codes (operations), Read codes (primary care), and UK Biobank self reported 
medical history field codes, along with related options in the provided endpoint 
definition file, curates prevalent and/or incident disease events for each UK 
Biobank participant at each assessment visit. See README.txt for more details.

Usage:
  curate_endpoint.R [--verbose] [--local-input] [--local-output] [--container-mode] [--work-dir <folder>] --def-file <file> --output <file>
  curate_endpoint.R -h | --help

Options:
  -h --help               Show this screen.
  --def-file <file>       File containing endpoint definition, see template.txt 
                          for example, and README.txt for details on options.
  --output <file>         Location to save the output file on persistent project
                          storage.
  --work-dir <folder>     Location on local cloud workstation used to download
                          and store local copies of the required underlying 
                          data, a local copy of the endpoint definition file, 
                          and a local copy of the output file before uploading
                          to project storage. Default: curate_endpoints_workdir/
  --local-input           Look for the endpoint definition file on the local 
                          cloud workstation instead of downloading from 
                          persistent storage.
  --local-output          Output file is saved to the local cloud workstation at 
                          the given location and not uploaded to persistent 
                          storage.
  --container-mode        Set this option when curate_endpoint.R does not have
                          access to project storage, e.g. when running on a
                          swiss-army-knife compute workstations launched via 
                          dx run. This tells curate_endpoint.R to look for the
                          data locally under /mnt/project, the endpoint 
                          definition file in the local home directory (e.g. 
                          uploaded by -iin in dx run swiss-army-knife) and save 
                          the output to local the home directory for upload to
                          the location set by dx run --output argument.
  --verbose               Prints status updates.
  " -> doc

# Parse arguments 
args <- docopt(doc, strict=TRUE)

vmessage("Input arguments read in by docopt:")
vmessage(show(args))

# Check for invalid combinations of arguments when using --container-mode
if (args[["--container-mode"]]) {
  if (args[["--local-input"]]) {
    stop("--local-input not valid with --container-mode")
  }
  if (args[["--local-output"]]) {
    stop("--local-input not valid with --container-mode")
  }
  if (dirname(args[["--def-file"]]) != ".") {
    stop("folder paths not accepted in --def-file when using --container-mode")
  }
  if (dirname(args[["--output"]]) != ".") {
    stop("folder paths not accepted in --def-file when using --container-mode")
  }
  if (!is.null(args[["--work-dir"]])) {
    stop("--work-dir cannot be set when --container-mode")
  }
} 

# Set/create temporary working directory
if (args[["--container-mode"]]) {
  work_dir <- "."
} else if (is.null(args[["--work-dir"]])) {
  work_dir <- "curate_endpoints_workdir"
} else {
  work_dir <- args[["--work-dir"]]
}
system(sprintf("mkdir -p '%s'", work_dir))

########################################################################################
# Load and check input def-file
########################################################################################

if (args[["--local-input"]] || args[["--container-mode"]]) {
  def_file <- args[["--def-file"]]
} else {
  check <- system(sprintf("dx download -f '%s' -o '%s'", args[["--def-file"]], work_dir))
  if (check > 0) {
    stop(sprintf("%s not found in current project storage", args[["--def-file"]]))
  }
  def_file <- system(sprintf("dx ls '%s'", args[["--def-file"]]), intern=TRUE) # resolve filename if provided as ID
  def_file <- sprintf("%s/%s", work_dir, basename(def_file))
}

if (!file.exists(def_file)) {
  stop(sprintf("--def-file '%s': file does not exist", def_file))
}

vmessage("Parsing endpoint definition file...")

opts <- readLines(def_file)
opts <- gsub("#.*$", "", opts) # strip out comments
opts <- gsub("\\s*$", "", opts) # remove trailing whitespace
opts <- opts[opts != ""] # remove empty lines

# Split into fields and codes
fields <- tolower(gsub(":.*", "", opts))
opts[!grepl(":", opts)] <- TRUE
opts <- gsub(".*:\\s*", "", opts)

# Check for valid fields
valid_options <- c(
  "filter by sex",
  "max follow years", 
  "max follow date", 
  "max follow age",
  "min follow years", 
  "min follow date", 
  "min follow age",
  "use conservative max censor date",
  "no icd-10 lookup in gp records",
  "no opcs-4 lookup in gp records",
  "exclude gp records",
  "filter to primary care sub-cohort",
  "exclude death register",
  "exclude hospital records",
  "exclude cancer register",
  "follow-up from birth",
  "impute follow-up as midpoint",
  "time since prevalent event",
  "most recent prevalent event",
  "do not impute missing event dates",
  "remove missing event dates",
  "allow missing event status",
  "use conservative min censor date"
)

valid_icd10 <- c(
  "icd-10", 
  "excluding icd-10", 
  "icd-10 primary cause only", 
  "excluding icd-10 primary cause only", 
  "prevalent icd-10", 
  "excluding prevalent icd-10", 
  "prevalent icd-10 primary cause only", 
  "excluding prevalent icd-10 primary cause only", 
  "incident icd-10", 
  "excluding incident icd-10", 
  "incident icd-10 primary cause only", 
  "excluding incident icd-10 primary cause only", 
  "non-fatal incident icd-10", 
  "excluding non-fatal incident icd-10", 
  "non-fatal incident icd-10 primary cause only", 
  "excluding non-fatal incident icd-10 primary cause only", 
  "fatal incident icd-10", 
  "excluding fatal incident icd-10", 
  "fatal incident icd-10 primary cause only", 
  "excluding fatal incident icd-10 primary cause only"
)

valid_icd9 <- c(
  "prevalent icd-9", 
  "excluding prevalent icd-9", 
  "prevalent icd-9 primary cause only", 
  "excluding prevalent icd-9 primary cause only"
)

valid_opcs4 <- c(
  "opcs-4", 
  "excluding opcs-4", 
  "opcs-4 primary cause only", 
  "excluding opcs-4 primary cause only", 
  "prevalent opcs-4", 
  "excluding prevalent opcs-4", 
  "prevalent opcs-4 primary cause only", 
  "excluding prevalent opcs-4 primary cause only", 
  "incident opcs-4",
  "excluding incident opcs-4",
  "incident opcs-4 primary cause only", 
  "excluding incident opcs-4 primary cause only"
)

valid_opcs3 <- c(
  "prevalent opcs-3", 
  "excluding prevalent opcs-3", 
  "prevalent opcs-3 primary cause only", 
  "excluding prevalent opcs-3 primary cause only"
)

valid_read3 <- c(
  "read-3", 
  "incident read-3",
  "prevalent read-3"
)

valid_read2 <- c(
  "read-2", 
  "incident read-2",
  "prevalent read-2"
)

valid_self_report <- c(
  "self-report",
  "incident self-report",
  "prevalent self-report"
)

valid <- c(valid_self_report, valid_opcs3, valid_opcs4, valid_icd9, valid_icd10,
           valid_read2, valid_read3, valid_options)

invalid <- setdiff(fields, valid)
if (length(invalid) > 0) {
  warning(sprintf("discarding invalid fields from endpoint file: %s", paste(invalid, collapse=", ")))
}

opts <- opts[which(fields %in% valid)]
fields <- fields[which(fields %in% valid)]
if (length(opts) == 0) {
  stop(sprintf("No valid fields found in --def-file '%s'", def_file))
}

# Get vector of codes for each field. Note fields may be repeated in the def file
opts <- lapply(unique(fields), function(ff) {
  this_opts <- opts[which(fields == ff)]
  if (ff %in% c(valid_icd9, valid_icd10, valid_opcs3, valid_opcs4)) {
    this_opts <- lapply(this_opts, parse_codelist)
    this_opts <- unique(unlist(this_opts))
  } else if (ff %in% c(valid_read2, valid_read3, valid_self_report)) {
    this_opts <- lapply(this_opts, parse_code_value_pairs)
    this_opts <- data.table(
      "codes"=unlist(lapply(this_opts, `[[`, 1)),
      "value_comparison_operators"=unlist(lapply(this_opts, `[[`, 2)),
      "values"=unlist(lapply(this_opts, `[[`, 3))
    )
    this_opts=unique(this_opts)
    if (all(is.na(this_opts$value_comparison_operators))) {
      this_opts[, value_comparison_operators := NULL]
      this_opts[, values := NULL]
    } 
    if (ff %in% valid_self_report) {
      setnames(this_opts, c("fields", "field_code_comparison_operators", "codes"))
    }
    this_opts <- as.list(this_opts)
  } else {
    this_opts <- unique(unlist(this_opts))
  }
  return(this_opts)
})
names(opts) <- unique(fields)

# Check fields with a single allowable entry:
n <- sapply(opts[valid_options], length)
names(n) <- valid_options
bad <- n[n > 1]
if (length(bad) > 0) {
  stop(sprintf("Multiple entries/options found for fields: %s", paste(paste0("'", names(bad), ":'"), collapse=", ")))
}

# Process options with values
if ("filter by sex" %in% names(opts)) {
  optval <- opts[["filter by sex"]]
  if (is.na(optval) || is.na(pmatch(tolower(optval), c("males", "females")))) {
    stop("'filter by sex' must be either \"Males\" or \"Females\"")
  }
  opts[["filter by sex"]] <- c("males", "females")[pmatch(tolower(optval),  c("males", "females"))]
}

if ("max follow years" %in% names(opts)) {
  processed <- as.integer(opts[["max follow years"]])
  if (is.na(processed) || processed != opts[["max follow years"]]) {
    stop("'max follow years:' must be given as an integer")
  } else {
    opts[["max follow years"]] <- processed
  }
}

if ("min follow years" %in% names(opts)) {
  processed <- as.integer(opts[["min follow years"]])
  if (is.na(processed) || processed != opts[["min follow years"]]) {
    stop("'min follow years:' must be given as an integer")
  } else {
    opts[["min follow years"]] <- processed
  }
}

if ("max follow age" %in% names(opts)) {
  processed <- as.integer(opts[["max follow age"]])
  if (is.na(processed) || processed != opts[["max follow age"]]) {
    stop("'max follow age:' must be given as an integer")
  } else {
    opts[["max follow age"]] <- processed
  }
}

if ("min follow age" %in% names(opts)) {
  processed <- as.integer(opts[["min follow age"]])
  if (is.na(processed) || processed != opts[["min follow age"]]) {
    stop("'min follow age:' must be given as an integer")
  } else {
    opts[["min follow age"]] <- processed
  }
}

if ("max follow date" %in% names(opts)) {
  tryCatch({
    opts[["max follow date"]] <- as.IDate(opts[["max follow date"]])
  }, error = function(e) {
    stop("Unable to unambiguously parse value given to 'max follow date:' into a valid date, please use format YYYY-MM-DD")
  })
}

if ("min follow date" %in% names(opts)) {
  tryCatch({
    opts[["min follow date"]] <- as.IDate(opts[["min follow date"]])
  }, error = function(e) {
    stop("Unable to unambiguously parse value given to 'min follow date:' into a valid date, please use format YYYY-MM-DD")
  })
}

# Check conflicts
if ("follow-up from birth" %in% names(opts) && any(names(opts) %like% "prevalent")) {
  stop("Prevalent events cannot be used when 'follow-up from birth' is set as an option")
}

if ("remove missing event dates" %in% names(opts) && !("do not impute missing event dates" %in% names(opts))) {
  warning("'remove missing event dates' implies 'do no impute missing event dates'; adding 'do not impute missing event dates' to options list")
  opts[["do not impute missing event dates"]] <- TRUE
}

if ("use self-report from current visit only" %in% names(opts) && ("incident self-report" %in% names(opts) || "self-report" %in% names(opts))) {
  stop("'use self-report from current visit only' is not compatible with 'self-report:' or 'incident self-report:'")
}

vmessage("Collated and simplified endpoint definition:")
vmessage(show(opts))

########################################################################################
# Determine data sets we actually need
########################################################################################

# Load ICD codes associated with cancer register
if (!("exclude cancer register" %in% names(opts))) {
  cancer_icd10 <- load_from_rap("common/Cancer Register/icd10_codes.csv")
  cancer_icd9 <- load_from_rap("common/Cancer Register/icd9_codes.csv")
}

# Split into prevalent vs. incident for later use as well
need_cancer_register_prevalent <- (
  !("exclude cancer register" %in% names(opts)) && any(names(opts) %in% c(
    "icd-10", "excluding icd-10", "prevalent icd-10", 
    "excluding prevalent icd-10", "prevalent icd-9", 
    "excluding prevalent icd-9"
  )) && (
    !is.null(opts[["icd-10"]]) && length(intersect(opts[["icd-10"]], cancer_icd10) > 0) ||
      !is.null(opts[["excluding icd-10"]]) && length(intersect(opts[["excluding icd-10"]], cancer_icd10) > 0) ||
      !is.null(opts[["prevalent icd-10"]]) && length(intersect(opts[["prevalent icd-10"]], cancer_icd10) > 0) ||
      !is.null(opts[["excluding prevalent icd-10"]]) && length(intersect(opts[["excluding prevalent icd-10"]], cancer_icd10) > 0) ||
      !is.null(opts[["prevalent icd-9"]]) && length(intersect(opts[["prevalent icd-9"]], cancer_icd9) > 0) ||
      !is.null(opts[["excluding prevalent icd-9"]]) && length(intersect(opts[["excluding prevalent icd-9"]], cancer_icd9) > 0) 
  )
)

need_cancer_register_incident <- (
  !("exclude cancer register" %in% names(opts)) && any(names(opts) %in% c(
    "icd-10", "excluding icd-10", "incident icd-10", 
    "excluding incident icd-10", "non-fatal incident icd-10", 
    "excluding non-fatal incident icd-10"
  )) && (
    !is.null(opts[["icd-10"]]) && length(intersect(opts[["icd-10"]], cancer_icd10) > 0) ||
      !is.null(opts[["excluding icd-10"]]) && length(intersect(opts[["excluding icd-10"]], cancer_icd10) > 0) ||
      !is.null(opts[["incident icd-10"]]) && length(intersect(opts[["incident icd-10"]], cancer_icd10) > 0) ||
      !is.null(opts[["excluding incident icd-10"]]) && length(intersect(opts[["excluding incident icd-10"]], cancer_icd10) > 0) ||
      !is.null(opts[["non-fatal incident icd-10"]]) && length(intersect(opts[["non-fatal incident icd-10"]], cancer_icd10) > 0) ||
      !is.null(opts[["excluding non-fatal incident icd-10"]]) && length(intersect(opts[["excluding non-fatal incident icd-10"]], cancer_icd10) > 0) 
  )
)

need_cancer_register <- need_cancer_register_prevalent || need_cancer_register_incident

need_hospital_records_prevalent <- (
  !("exclude hospital records" %in% names(opts)) && any(names(opts) %in% c(
    "icd-10", "excluding icd-10", "icd-10 primary cause only", 
    "excluding icd-10 primary cause only", "prevalent icd-10", 
    "excluding prevalent icd-10", "prevalent icd-10 primary cause only", 
    "excluding prevalent icd-10 primary cause only", "prevalent icd-9", 
    "excluding prevalent icd-9", "prevalent icd-9 primary cause only", 
    "excluding prevalent icd-9 primary cause only"
  ))
)

need_hospital_records_incident <- (
  !("exclude hospital records" %in% names(opts)) && any(names(opts) %in% c(
    "icd-10", "excluding icd-10", "icd-10 primary cause only", 
    "excluding icd-10 primary cause only",  "incident icd-10", 
    "excluding incident icd-10", "incident icd-10 primary cause only", 
    "excluding incident icd-10 primary cause only", "non-fatal incident icd-10", 
    "excluding non-fatal incident icd-10", 
    "non-fatal incident icd-10 primary cause only", 
    "excluding non-fatal incident icd-10 primary cause only"
  ))
)

need_hospital_records <- need_hospital_records_prevalent || need_hospital_records_incident

need_hospital_operations_prevalent <- (
  !("exclude hospital records" %in% names(opts)) && any(names(opts) %in% c(
    "opcs-4", "excluding opcs-4", "opcs-4 primary cause only", 
    "excluding opcs-4 primary cause only", "prevalent opcs-4", 
    "excluding prevalent opcs-4", "prevalent opcs-4 primary cause only", 
    "excluding prevalent opcs-4 primary cause only", "prevalent opcs-3", 
    "excluding prevalent opcs-3", "prevalent opcs-3 primary cause only", 
    "excluding prevalent opcs-3 primary cause only"
  ))
)

need_hospital_operations_incident <- (
  !("exclude hospital records" %in% names(opts)) && any(names(opts) %in% c(
    "opcs-4", "excluding opcs-4", "opcs-4 primary cause only", 
    "excluding opcs-4 primary cause only", "incident opcs-4", 
    "excluding incident opcs-4", "incident opcs-4 primary cause only", 
    "excluding incident opcs-4 primary cause only"
  ))
)

need_hospital_operations <- need_hospital_operations_prevalent || need_hospital_operations_incident

need_death_records <- (
  !("exclude death register" %in% names(opts)) && any(names(opts) %in% c(
    "icd-10", "excluding icd-10", "icd-10 primary cause only", 
    "excluding icd-10 primary cause only", "incident icd-10", 
    "excluding incident icd-10", "incident icd-10 primary cause only", 
    "excluding incident icd-10 primary cause only", "fatal incident icd-10", 
    "excluding fatal incident icd-10", 
    "fatal incident icd-10 primary cause only", 
    "excluding fatal incident icd-10 primary cause only"
  ))
)

need_gp_records_prevalent <- (
  !("exclude gp records" %in% names(opts)) && (
    any(names(opts) %in% c("read-3", "prevalent read-3", "read-2", "prevalent read-2")) ||
    !("no icd-10 lookup in gp records" %in% names(opts)) && any(names(opts) %in% c(
      "icd-10", "prevalent icd-10", "excluding icd-10", "excluding prevalent icd-10"
    )) || !("no opcs-4 lookup in gp records" %in% names(opts)) && any(names(opts) %in% c(
      "opcs-4", "prevalent opcs-4", "excluding opcs-4", "excluding prevalent opcs-4"
    ))
  )
)

need_gp_records_incident <- (
  !("exclude gp records" %in% names(opts)) && (
    any(names(opts) %in% c("read-3", "incident read-3", "read-2", "incident read-2")) ||
    !("no icd-10 lookup in gp records" %in% names(opts)) && any(names(opts) %in% c(
      "icd-10", "incident icd-10", "excluding icd-10", "excluding incident icd-10"
    )) || !("no opcs-4 lookup in gp records" %in% names(opts)) && any(names(opts) %in% c(
      "opcs-4", "incident opcs-4", "excluding opcs-4", "excluding incident opcs-4"
    ))
  )
)

need_gp_records <- need_gp_records_prevalent || need_gp_records_incident

need_self_report_prevalent <- (
  any(names(opts) %in% c("self-report", "prevalent self-report"))
)

need_self_report_incident <- (
  any(names(opts) %in% c("self-report", "incident self-report"))
)

need_self_report <- need_self_report_prevalent || need_self_report_incident 

need_prevalent <- need_cancer_register_prevalent || need_gp_records_prevalent || 
  need_hospital_operations_prevalent || need_hospital_records_prevalent || 
  need_self_report_prevalent

need_incident <- need_cancer_register_incident || need_gp_records_incident ||
  need_hospital_operations_incident || need_hospital_records_incident ||
  need_self_report_incident || need_death_records

vmessage("Using the following data sources:")
if (need_death_records) vmessage("  Death registry data ('common/Deaths/death_causes.csv')")
if (need_hospital_records) vmessage("  Hospital records - Diagnoses ('common/Hospital Records/diagnoses.csv')")
if (need_hospital_operations) vmessage("  Hospital records - Operations & Procedures ('common/Hospital Records/operations.csv')")
if (need_cancer_register) vmessage("  Cancer registry data ('common/Cancer Register/cancer_register.csv')")
if (need_gp_records) vmessage("  Primary care data ('common/Primary Care/gp_clinical_records.csv')")
if (need_self_report) vmessage("  Self-reported medical history ('common/Medical History/medical_history.csv')")

########################################################################################
# Curate information about follow-up time for each participant and assessment
########################################################################################

vmessage("Initializing events and follow-up table...")

vmessage("  Loading demographic data at each UK Biobank assessment...")
demographics <- load_from_rap("common/Demographics/demographics.csv")

if ('filter by sex' %in% names(opts)) {
  if (opts[["filter by sex"]] == "males") {
    vmessage("  Filtering to male participants...")
    demographics <- demographics[sex == "Male"]
  } else {
    vmessage("  Filtering to female participants...")
    demographics <- demographics[sex == "Female"]
  }
}

vmessage("  Loading follow-up information on each participant...")
followup <- load_from_rap("common/Follow-up/follow_up.csv")

if ("follow-up from birth" %in% names(opts)) {
  events <- demographics[visit_index == 0, .(eid, visit_index=-1, assessment_date=approx_birth_date, approx_birth_date)]
} else {
  events <- demographics[,.(eid, visit_index, assessment_date, approx_birth_date)]
}

if ("filter to primary care sub-cohort" %in% names(opts)) {
  vmessage("  Filtering to participants with primary care data...")
  followup <- followup[(linked_gp_records)]
  events <- events[eid %in% followup$eid]
}

if (need_prevalent) {
  vmessage("  Determining minimum censor dates...")
  
  min_followup <- melt(followup, id.vars="eid", measure.vars=c(
    "min_hospital_censor", "min_cancer_censor", "min_gp_censor"
  ), variable.name="id_column", value.name="censor_date", na.rm=TRUE)
  min_followup[, data_source := "common/Follow-up/follow_up.csv"]
  min_followup <- rbind(min_followup, demographics[visit_index == 0, .(eid, 
    data_source="common/Demographics/demographics.csv", 
    id_column="approx_birth_date",
    censor_date=approx_birth_date
  )])
  
  if (!(need_hospital_operations_prevalent) && !(need_hospital_records_prevalent)) min_followup <- min_followup[id_column != "min_hospital_censor"]
  if (!(need_cancer_register_prevalent)) min_followup <- min_followup[id_column != "min_cancer_censor"]
  if (!(need_gp_records_prevalent)) min_followup <- min_followup[id_column != "min_gp_censor"]
  if (!(need_self_report_prevalent)) min_followup <- min_followup[id_column != "approx_birth_date"]

  if ("use conservative min censor date" %in% names(opts)) {
    min_followup <- min_followup[order(-censor_date)][order(eid)]
    min_followup <- min_followup[,.SD[1], by=eid] # maximum censor date per person
  } else {
    min_followup <- min_followup[order(censor_date)][order(eid)]
    min_followup <- min_followup[,.SD[1], by=eid] # minimum censor date per person
  }
  
  events[min_followup, on = .(eid), 
    c("min_censor_date", "min_censor_reason", "min_censor_data_source", "min_censor_id_column") := 
    .(i.censor_date, "minimum censor date", i.data_source, i.id_column)]
  
  if ("min follow years" %in% names(opts)) {
    events[, min_censor_truncation := add_years(assessment_date, -opts[["min follow years"]])]
    events[min_censor_truncation > min_censor_date, 
           c("min_censor_date", "min_censor_reason", "min_censor_data_source", "min_censor_id_column") :=
             .(min_censor_truncation, sprintf("minimum follow years: %s", opts[["min follow years"]]), args[["--def-file"]], NA)]
    events[, min_censor_truncation := NULL]
  } 
} else if ("follow-up from birth" %in% names(opts)) {
  events[,
   c("min_censor_date", "min_censor_reason", "min_censor_data_source", "min_censor_id_column") := 
   .(approx_birth_date, "Date of birth (approximate)", "common/Demographics/demographics.csv", "approx_birth_date")]
} else {
  events[,
   c("min_censor_date", "min_censor_reason", "min_censor_data_source", "min_censor_id_column") := 
   .(assessment_date, "Date of UKB assessment", "common/Demographics/demographics.csv", "assessment_date")]
}

if ("min follow date" %in% names(opts)) {
  events[opts[["min follow date"]] > min_censor_date, 
    c("min_censor_date", "min_censor_reason", "min_censor_data_source", "min_censor_id_column") :=
    .(as.IDate(opts[["min follow date"]]), sprintf("minimum follow date: %s", opts[["min follow date"]]), args[["--def-file"]], NA)]
  if(events[min_censor_date > assessment_date, .N] > 0) {
    vmessage(sprintf("    Removing %s rows where the 'min follow date' occurs after UK Biobank assessment date", events[min_censor_date > assessment_date, .N]))
    events <- events[assessment_date >= min_censor_date]
  }
} 
if ("min follow age" %in% names(opts)) {
  events[, min_censor_truncation := add_years(approx_birth_date, opts[["min follow age"]])]
  events[min_censor_truncation > min_censor_date, 
    c("min_censor_date", "min_censor_reason", "min_censor_data_source", "min_censor_id_column") :=
    .(min_censor_truncation, sprintf("minimum follow age: %s", opts[["min follow age"]]), args[["--def-file"]], NA)]
  events[, min_censor_truncation := NULL]
  if(events[min_censor_date > assessment_date, .N] > 0) {
    vmessage(sprintf("    Removing %s rows where the participant is younger than the 'min follow age' at UK Biobank assessment", events[min_censor_date > assessment_date, .N]))
    events <- events[assessment_date >= min_censor_date]
  }
}

if (need_incident) {
  vmessage("  Determining maximum censor dates...")
  
  max_followup <- melt(followup, id.vars="eid", measure.vars=c(
    "max_death_censor", "max_hospital_censor", "max_cancer_censor", "max_gp_censor"
  ), variable.name="id_column", value.name="censor_date", na.rm=TRUE)
  max_followup[, data_source := "common/Follow-up/follow_up.csv"]
  max_followup <- rbind(max_followup, demographics[, .(
    data_source="common/Demographics/demographics.csv", 
    id_column="assessment_date",
    censor_date=max(assessment_date)
  ), by=eid])
  
  if (!(need_death_records)) max_followup <- max_followup[id_column != "max_death_censor"]
  if (!(need_hospital_operations_incident) && !(need_hospital_records_incident)) max_followup <- max_followup[id_column != "max_hospital_censor"]
  if (!(need_cancer_register_incident)) max_followup <- max_followup[id_column != "max_cancer_censor"]
  if (!(need_gp_records_incident)) max_followup <- max_followup[id_column != "max_gp_censor"]
  if (!(need_self_report_incident)) max_followup <- max_followup[id_column != "assessment_date"]
  
  if ("use conservative max censor date" %in% names(opts)) {
    max_followup <- max_followup[order(censor_date)][order(eid)]
    max_followup <- max_followup[,.SD[1], by=eid] # minimum censor date per person
  } else if (need_hospital_operations || need_hospital_records) {
    max_followup <- max_followup[id_column == "max_hospital_censor"] # Default behavior
  } else {
    # if no hospital record lookup, default to max follow-up date available
    max_followup <- max_followup[order(-censor_date)][order(eid)]
    max_followup <- max_followup[,.SD[1], by=eid] # maximum censor date per person
  }
  
  events[max_followup, on = .(eid), 
    c("max_censor_date", "max_censor_reason", "max_censor_data_source", "max_censor_id_column") := 
    .(i.censor_date, "maximum censor date", i.data_source, i.id_column)]
  
  if ("max follow years" %in% names(opts)) {
    events[, max_censor_truncation := add_years(assessment_date, opts[["max follow years"]])]
    events[max_censor_truncation < max_censor_date, 
      c("max_censor_date", "max_censor_reason", "max_censor_data_source", "max_censor_id_column") :=
      .(max_censor_truncation, sprintf("maximum follow years: %s", opts[["max follow years"]]), args[["--def-file"]], NA)]
    events[, max_censor_truncation := NULL]
  } 
}

if ("max follow date" %in% names(opts)) {
  events[opts[["max follow date"]] < max_censor_date, 
    c("max_censor_date", "max_censor_reason", "max_censor_data_source", "max_censor_id_column") :=
    .(as.IDate(opts[["max follow date"]]), sprintf("maximum follow date: %s", opts[["max follow date"]]), args[["--def-file"]], NA)]
  if(events[max_censor_date < assessment_date, .N] > 0) {
    vmessage(sprintf("    Removing %s rows where the 'max follow date' occurs before UK Biobank assessment date", events[max_censor_date < assessment_date, .N]))
    events <- events[assessment_date <= max_censor_date]
  }
} 
if ("max follow age" %in% names(opts)) {
  events[, max_censor_truncation := add_years(approx_birth_date, opts[["max follow age"]])]
  events[max_censor_truncation < max_censor_date, 
    c("max_censor_date", "max_censor_reason", "max_censor_data_source", "max_censor_id_column") :=
    .(max_censor_truncation, sprintf("maximum follow age: %s", opts[["max follow age"]]), args[["--def-file"]], NA)]
  events[, max_censor_truncation := NULL]
  if(events[max_censor_date < assessment_date, .N] > 0) {
    vmessage(sprintf("    Removing %s rows where the participant is older than the 'max follow age' at UK Biobank assessment", events[max_censor_date < assessment_date, .N]))
    events <- events[assessment_date <= max_censor_date]
  }
}

########################################################################################
# Identify events
########################################################################################
vmessage("Identifying events...")

# Set up incident and prevalent event columns
if (need_incident) {
  events[, inci_censor_date := add_day(max_censor_date)] # Offset max_censor_date by 1 so we can always use < operator.
  events[, inci_censor_data_source := NA_character_]
  events[, inci_censor_id_column := NA_character_]
  events[, inci_censor_index := NA_character_]
}

if (need_prevalent) {
  events[, prev_censor_date := add_day(assessment_date)]
  events[, prev_censor_data_source := NA_character_]
  events[, prev_censor_id_column := NA_character_]
  events[, prev_censor_index := NA_character_]
}

if (need_death_records) {
  vmessage("  Loading death register data...")
  deaths <- load_from_rap("common/Deaths/deaths.csv")
  death_causes <- load_from_rap("common/Deaths/death_causes.csv")
  
  # Filter to people who have consented for health record linkage (and to 
  # the sub-cohort with primary care data if that filter has been applied)
  vmessage("    Filtering to participants consented for electronic health record linkage...")
  deaths <- deaths[eid %in% followup$eid]
  death_causes <- death_causes[eid %in% followup$eid]
  
  # Make sure we take the most recent death certificate where there are multiple
  # (e.g. where cause of death has been refined)
  vmessage("    Extracting most up to date death certificate where multiple have been issued...")
  deaths <- deaths[,.SD[which.max(certificate_number)], by=eid]
  death_causes <- death_causes[dnx_death_id %in% deaths$dnx_death_id]
  death_causes[deaths, on = .(dnx_death_id), date_of_death := i.date_of_death]

  vmessage("    Finding deaths matching endpoint definition...")

  excl_opts <- c(
    "excluding icd-10", 
    "excluding incident icd-10", 
    "excluding fatal incident icd-10"
  )
  for (opt_type in excl_opts) {
    if (opt_type %in% names(opts)) {
      vmessage(sprintf("      Removing codes in '%s'...", opt_type))
      codes <- opts[[opt_type]]
      for (code in codes) {
        vmessage(sprintf("        Removing %s and sub-codes...", code))
        code <- format_code(code)
        death_causes <- death_causes[!(cause_icd10 %like% code)]
      }
    }
  }
  
  excl_opts <- c(
    "excluding icd-10 primary cause only", 
    "excluding incident icd-10 primary cause only",
    "excluding fatal incident icd-10 primary cause only" 
  )
  for (opt_type in excl_opts) {
    if (opt_type %in% names(opts)) {
      vmessage(sprintf("      Removing codes in '%s'...", opt_type))
      codes <- opts[[opt_type]]
      for (code in codes) {
        vmessage(sprintf("        Removing %s and sub-codes...", code))
        code <- format_code(code)
        death_causes <- death_causes[!(cause_icd10 %like% code & cause_number == 1)]
      }
    }
  }

  match_opts <- c("icd-10", "incident icd-10", "fatal incident icd-10")
  for (opt_type in match_opts) {
    if (opt_type %in% names(opts)) {
      vmessage(sprintf("      Matching codes in '%s'...", opt_type))
      codes <- opts[[opt_type]]
      for (code in codes) {
        vmessage(sprintf("        Matching %s and sub-codes...", code))
        code <- format_code(code)
        matches <- death_causes[cause_icd10 %like% code]
        events[matches, on = .(eid, inci_censor_date > date_of_death),
          c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
          .(date_of_death, "common/Deaths/death_causes.csv", "dnx_death_cause_id", dnx_death_cause_id)]
      }
    }
  }
  
  match_opts <- c(
    "icd-10 primary cause only", 
    "incident icd-10 primary cause only", 
    "fatal incident icd-10 primary cause only"
  )
  for (opt_type in match_opts) {
    if (opt_type %in% names(opts)) {
      vmessage(sprintf("      Matching codes in '%s'...", opt_type))
      codes <- opts[[opt_type]]
      for (code in codes) {
        vmessage(sprintf("        Matching %s and sub-codes...", code))
        code <- format_code(code)
        matches <- death_causes[cause_icd10 %like% code & cause_number == 1]
        events[matches, on = .(eid, inci_censor_date > date_of_death),
          c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
          .(date_of_death, "common/Deaths/death_causes.csv", "dnx_death_cause_id", dnx_death_cause_id)]
      }
    }
  }
  
  # Remove tables we no longer need to free up memory
  vmessage("  Unloading death register data...")
  rm(deaths, death_causes)
  invisible(gc())
}

if (need_hospital_records) {
  vmessage("  Loading hospital records diagnosis data...")
  hes <- load_from_rap("common/Hospital Records/diagnoses.csv")
  
  # Filter to people who have consented for health record linkage (and to 
  # the sub-cohort with primary care data if that filter has been applied)
  vmessage("    Filtering to participants consented for electronic health record linkage...")
  hes <- hes[eid %in% followup$eid]

  vmessage("    Sorting events...")
  hes <- hes[order(cause_type)][order(event_start)][order(eid)]
  
  if (need_hospital_records_incident) {
    vmessage("    Finding incident diagnoses in hospital records matching endpoint definition...")
    
    # Create a copy before applying any filters. If no filters are applied, only
    # a pointer is stored, not a full additional copy of 'hes'
    hes_inci <- hes
    
    excl_opts <- c(
      "excluding icd-10", 
      "excluding incident icd-10", 
      "excluding non-fatal incident icd-10"
    )
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          hes_inci <- hes_inci[!(icd_version == 10 & diag_icd %like% code)]
        }
      }
    }
    
    excl_opts <- c(
      "excluding icd-10 primary cause only", 
      "excluding incident icd-10 primary cause only",
      "excluding non-fatal incident icd-10 primary cause only" 
    )
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          hes_inci <- hes_inci[!(icd_version == 10 & diag_icd %like% code & cause_type == "primary")]
        }
      }
    }
    
    match_opts <- c("icd-10", "incident icd-10", "non-fatal incident icd-10")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          first_matches <- hes_inci[icd_version == 10 & diag_icd %like% code, .SD[1], by=eid] 
          events[first_matches, on = .(eid, inci_censor_date > event_start, assessment_date < event_start),
            c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
            .(event_start, "common/Hospital Records/diagnoses.csv", "dnx_hesin_diag_id", dnx_hesin_diag_id)]
        }
      }
    }
    
    match_opts <- c(
      "icd-10 primary cause only", 
      "incident icd-10 primary cause only", 
      "non-fatal incident icd-10 primary cause only"
    )
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          first_matches <- hes_inci[icd_version == 10 & cause_type == "primary" & diag_icd %like% code, .SD[1], by=eid] 
          events[first_matches, on = .(eid, inci_censor_date > event_start, assessment_date < event_start),
            c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
            .(event_start, "common/Hospital Records/diagnoses.csv", "dnx_hesin_diag_id", dnx_hesin_diag_id)]
        }
      }
    }
    
    rm(hes_inci) # table no longer needed
  }
  
  if (need_hospital_records_prevalent) {
    vmessage("    Finding prevalent diagnoses in hospital records matching endpoint definition...")
    
    # Create a copy before applying any filters. If no filters are applied, only
    # a pointer is stored, not a full additional copy of 'hes'
    hes_prev <- hes
    
    excl_opts <- c(
      "excluding icd-10", 
      "excluding prevalent icd-10" 
    )
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          hes_prev <- hes_prev[!(icd_version == 10 & diag_icd %like% code)]
        }
      }
    }
    
    excl_opts <- "excluding prevalent icd-9"
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          hes_prev <- hes_prev[!(icd_version == 9 & diag_icd %like% code)] 
        }
      }
    }
    
    excl_opts <- c(
      "excluding icd-10 primary cause only", 
      "excluding prevalent icd-10 primary cause only"
    )
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          hes_prev <- hes_prev[!(icd_version == 10 & diag_icd %like% code & cause_type == "primary")]
        }
      }
    }
    
    excl_opts <- "excluding prevalent icd-9 primary cause only"
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          hes_prev <- hes_prev[!(icd_version == 9 & diag_icd %like% code & cause_type == "primary")]
        }
      }
    }

    match_opts <- c("icd-10", "prevalent icd-10")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          if ("most recent prevalent event" %in% names(opts)) {
            most_recent <- hes_prev[icd_version == 10 & diag_icd %like% code][
              events, on = .(eid, event_end > min_censor_date, event_end < prev_censor_date),
              mult="last", nomatch=0, .(eid, visit_index=i.visit_index, dnx_hesin_diag_id, event_end=x.event_end)
            ]
            events[most_recent, on = .(eid, visit_index), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(event_end, "common/Hospital Records/diagnoses.csv", "dnx_hesin_diag_id", dnx_hesin_diag_id)]
          } else {
            first_matches <- hes_prev[icd_version == 10 & diag_icd %like% code, .SD[1], by=eid]
            events[first_matches, on = .(eid, min_censor_date < event_start, prev_censor_date > event_start),
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(event_start, "common/Hospital Records/diagnoses.csv", "dnx_hesin_diag_id", dnx_hesin_diag_id)]
          }
        }
      }
    }
    
    match_opts <- c("prevalent icd-9")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          if ("most recent prevalent event" %in% names(opts)) {
            most_recent <- hes_prev[icd_version == 9 & diag_icd %like% code][
              events, on = .(eid, event_end > min_censor_date, event_end < prev_censor_date),
              mult="last", nomatch=0, .(eid, visit_index=i.visit_index, dnx_hesin_diag_id, event_end=x.event_end)
            ]
            events[most_recent, on = .(eid, visit_index), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(event_end, "common/Hospital Records/diagnoses.csv", "dnx_hesin_diag_id", dnx_hesin_diag_id)]
          } else {
            first_matches <- hes_prev[icd_version == 9 & diag_icd %like% code, .SD[1], by=eid]
            events[first_matches, on = .(eid, min_censor_date < event_start, prev_censor_date > event_start),
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(event_start, "common/Hospital Records/diagnoses.csv", "dnx_hesin_diag_id", dnx_hesin_diag_id)]
          }
        }
      }
    }
    
    match_opts <- c("icd-10 primary cause only", "prevalent icd-10 primary cause only")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          if ("most recent prevalent event" %in% names(opts)) {
            most_recent <- hes_prev[icd_version == 10 & diag_icd %like% code & cause_type == "primary"][
              events, on = .(eid, event_end > min_censor_date, event_end < prev_censor_date),
              mult="last", nomatch=0, .(eid, visit_index=i.visit_index, dnx_hesin_diag_id, event_end=x.event_end)
            ]
            events[most_recent, on = .(eid, visit_index), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(event_end, "common/Hospital Records/diagnoses.csv", "dnx_hesin_diag_id", dnx_hesin_diag_id)]
          } else {
            first_matches <- hes_prev[icd_version == 10 & diag_icd %like% code & cause_type == "primary", .SD[1], by=eid]
            events[first_matches, on = .(eid, min_censor_date < event_start, prev_censor_date > event_start),
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(event_start, "common/Hospital Records/diagnoses.csv", "dnx_hesin_diag_id", dnx_hesin_diag_id)]
          }
        }
      }
    }
    
    match_opts <- c("prevalent icd-9 primary cause only")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          if ("most recent prevalent event" %in% names(opts)) {
            most_recent <- hes_prev[icd_version == 9 & diag_icd %like% code & cause_type == "primary"][
              events, on = .(eid, event_end > min_censor_date, event_end < prev_censor_date),
              mult="last", nomatch=0, .(eid, visit_index=i.visit_index, dnx_hesin_diag_id, event_end=x.event_end)
            ]
            events[most_recent, on = .(eid, visit_index), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(event_end, "common/Hospital Records/diagnoses.csv", "dnx_hesin_diag_id", dnx_hesin_diag_id)]
          } else {
            first_matches <- hes_prev[icd_version == 9 & diag_icd %like% code & cause_type == "primary", .SD[1], by=eid]
            events[first_matches, on = .(eid, min_censor_date < event_start, prev_censor_date > event_start),
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(event_start, "common/Hospital Records/diagnoses.csv", "dnx_hesin_diag_id", dnx_hesin_diag_id)]
          }
        }
      }
    }
    
    rm(hes_prev); # table no longer needed
  }
  
  # Remove tables we no longer need to free up memory
  vmessage("  Unloading hospital records diagnosis data......")
  rm(hes)
  invisible(gc())
}

if (need_hospital_operations) {
  vmessage("  Loading hospital records operations and procedures data...")
  hop <- load_from_rap("common/Hospital Records/operations.csv")
  
  # Filter to people who have consented for health record linkage (and to 
  # the sub-cohort with primary care data if that filter has been applied)
  vmessage("    Filtering to participants consented for electronic health record linkage...")
  hop <- hop[eid %in% followup$eid]
  
  vmessage("    Sorting events...")
  hop <- hop[order(procedure_type)][order(imputed_operation_date)][order(eid)] # imputed operation date imputes missing dates as the start of the corresponding hospital admission
  
  if (need_hospital_operations_incident) {
    vmessage("    Finding incident operations and procedures in hospital records matching endpoint definition...")
    
    # Create a copy before applying any filters. If no filters are applied, only
    # a pointer is stored, not a full additional copy of 'hes'
    hop_inci <- hop
    
    excl_opts <- c(
      "excluding opcs-4", 
      "excluding incident opcs-4" 
    )
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          hop_inci <- hop_inci[!(opcs_version == 4 & opcs_code %like% code)]
        }
      }
    }
    
    excl_opts <- c(
      "excluding opcs-4 primary cause only", 
      "excluding incident opcs-4 primary cause only"
    )
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          hop_inci <- hop_inci[!(opcs_version == 4 & opcs_code %like% code & procedure_type == "main")]
        }
      }
    }
    
    match_opts <- c("opcs-4", "incident opcs-4")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          first_matches <- hop_inci[opcs_version == 4 & opcs_code %like% code, .SD[1], by=eid] 
          events[first_matches, on = .(eid, inci_censor_date > imputed_operation_date, assessment_date < imputed_operation_date),
            c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
            .(imputed_operation_date, "common/Hospital Records/operations.csv", "dnx_hesin_oper_id", dnx_hesin_oper_id)]
        }
      }
    }
    
    match_opts <- c(
      "opcs-4 primary cause only", 
      "incident opcs-4 primary cause only"
    )
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          first_matches <- hop_inci[opcs_version == 4 & opcs_code %like% code & procedure_type == "main", .SD[1], by=eid] 
          events[first_matches, on = .(eid, inci_censor_date > imputed_operation_date, assessment_date < imputed_operation_date),
            c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
            .(imputed_operation_date, "common/Hospital Records/operations.csv", "dnx_hesin_oper_id", dnx_hesin_oper_id)]
        }
      }
    }
    
    rm(hop_inci) # table no longer needed
  }
  
  if (need_hospital_operations_prevalent) {
    vmessage("    Finding prevalent operations and procedures in hospital records matching endpoint definition...")
    
    # Create a copy before applying any filters. If no filters are applied, only
    # a pointer is stored, not a full additional copy of 'hes'
    hop_prev <- hop
    
    excl_opts <- c(
      "excluding opcs-4", 
      "excluding prevalent opcs-4" 
    )
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          hop_prev <- hop_prev[!(opcs_version == 4 & opcs_code %like% code)]
        }
      }
    }
    
    excl_opts <- "excluding prevalent opcs-3"
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          hop_prev <- hop_prev[!(opcs_version == 3 & opcs_code %like% code)]
        }
      }
    }
    
    excl_opts <- c(
      "excluding opcs-4 primary cause only", 
      "excluding prevalent opcs-4 primary cause only"
    )
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          hop_prev <- hop_prev[!(opcs_version == 4 & opcs_code %like% code & procedure_type == "main")]
        }
      }
    }
    
    excl_opts <- "excluding prevalent opcs-3 primary cause only"
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          hop_prev <- hop_prev[!(opcs_version == 3 & opcs_code %like% code & procedure_type == "main")]
        }
      }
    }
    
    match_opts <- c("opcs-4", "prevalent opcs-4")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          if ("most recent prevalent event" %in% names(opts)) {
            most_recent <- hop_prev[opcs_version == 4 & opcs_code %like% code][
              events, on = .(eid, imputed_operation_date > min_censor_date, imputed_operation_date < prev_censor_date),
              mult="last", nomatch=0, .(eid, visit_index=i.visit_index, dnx_hesin_oper_id, imputed_operation_date=x.imputed_operation_date)
            ]
            events[most_recent, on = .(eid, visit_index), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(imputed_operation_date, "common/Hospital Records/operations.csv", "dnx_hesin_oper_id", dnx_hesin_oper_id)]
          } else {
            first_matches <- hop_prev[opcs_version == 4 & opcs_code %like% code, .SD[1], by=eid]
            events[first_matches, on = .(eid, min_censor_date < imputed_operation_date, prev_censor_date > imputed_operation_date),
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(imputed_operation_date, "common/Hospital Records/operations.csv", "dnx_hesin_oper_id", dnx_hesin_oper_id)]
          }
        }
      }
    }
    
    match_opts <- c("prevalent opcs-3")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          if ("most recent prevalent event" %in% names(opts)) {
            most_recent <- hop_prev[opcs_version == 3 & opcs_code %like% code][
              events, on = .(eid, imputed_operation_date > min_censor_date, imputed_operation_date < prev_censor_date),
              mult="last", nomatch=0, .(eid, visit_index=i.visit_index, dnx_hesin_oper_id, imputed_operation_date=x.imputed_operation_date)
            ]
            events[most_recent, on = .(eid, visit_index), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(imputed_operation_date, "common/Hospital Records/operations.csv", "dnx_hesin_oper_id", dnx_hesin_oper_id)]
          } else {
            first_matches <- hop_prev[opcs_version == 3 & opcs_code %like% code, .SD[1], by=eid]
            events[first_matches, on = .(eid, min_censor_date < imputed_operation_date, prev_censor_date > imputed_operation_date),
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(imputed_operation_date, "common/Hospital Records/operations.csv", "dnx_hesin_oper_id", dnx_hesin_oper_id)]
          }
        }
      }
    }
    
    match_opts <- c("opcs-4 primary cause only", "prevalent opcs-4 primary cause only")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          if ("most recent prevalent event" %in% names(opts)) {
            most_recent <- hop_prev[opcs_version == 4 & opcs_code %like% code & procedure_type == "main"][
              events, on = .(eid, imputed_operation_date > min_censor_date, imputed_operation_date < prev_censor_date),
              mult="last", nomatch=0, .(eid, visit_index=i.visit_index, dnx_hesin_oper_id, imputed_operation_date=x.imputed_operation_date)
            ]
            events[most_recent, on = .(eid, visit_index), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(imputed_operation_date, "common/Hospital Records/operations.csv", "dnx_hesin_oper_id", dnx_hesin_oper_id)]
          } else {
            first_matches <- hop_prev[opcs_version == 4 & opcs_code %like% code & procedure_type == "main", .SD[1], by=eid]
            events[first_matches, on = .(eid, min_censor_date < imputed_operation_date, prev_censor_date > imputed_operation_date),
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(imputed_operation_date, "common/Hospital Records/operations.csv", "dnx_hesin_oper_id", dnx_hesin_oper_id)]
          }
        }
      }
    }
    
    match_opts <- c("prevalent opcs-3 primary cause only")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          if ("most recent prevalent event" %in% names(opts)) {
            most_recent <- hop_prev[opcs_version == 3 & opcs_code %like% code & procedure_type == "main"][
              events, on = .(eid, imputed_operation_date > min_censor_date, imputed_operation_date < prev_censor_date),
              mult="last", nomatch=0, .(eid, visit_index=i.visit_index, dnx_hesin_oper_id, imputed_operation_date=x.imputed_operation_date)
            ]
            events[most_recent, on = .(eid, visit_index), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(imputed_operation_date, "common/Hospital Records/operations.csv", "dnx_hesin_oper_id", dnx_hesin_oper_id)]
          } else {
            first_matches <- hop_prev[opcs_version == 3 & opcs_code %like% code & procedure_type == "main", .SD[1], by=eid]
            events[first_matches, on = .(eid, min_censor_date < imputed_operation_date, prev_censor_date > imputed_operation_date),
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(imputed_operation_date, "common/Hospital Records/operations.csv", "dnx_hesin_oper_id", dnx_hesin_oper_id)]
          }
        }
      }
    }
    
    rm(hop_prev); # table no longer needed
  }
  
  # Remove tables we no longer need to free up memory
  vmessage("  Unloading hospital records operations and procedures data...")
  rm(hop)
  invisible(gc())
}

if (need_cancer_register) {
  vmessage("  Loading cancer register data...")
  cr <- load_from_rap("common/Cancer Register/cancer_register.csv")
  cr <- cr[!is.na(diagnosis_date)] # All are ICD-9 diagnoses with empty code strings 
  
  # Filter to people who have consented for health record linkage (and to 
  # the sub-cohort with primary care data if that filter has been applied)
  vmessage("    Filtering to participants consented for electronic health record linkage...")
  cr <- cr[eid %in% followup$eid]
  
  vmessage("    Sorting events...")
  cr <- cr[order(diagnosis_date)][order(eid)]
  
  if (need_cancer_records_incident) {
    vmessage("    Finding incident cancer register records matching endpoint definition...")
    
    # Create a copy before applying any filters. If no filters are applied, only
    # a pointer is stored, not a full additional copy of 'hes'
    cr_inci <- cr
    
    excl_opts <- c(
      "excluding icd-10", 
      "excluding incident icd-10", 
      "excluding non-fatal incident icd-10"
    )
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          cr_inci <- cr_inci[!(icd_version == 10 & diag_icd %like% code)]
        }
      }
    }
    
    match_opts <- c("icd-10", "incident icd-10", "non-fatal incident icd-10")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          first_matches <- cr_inci[icd_version == 10 & diag_icd %like% code, .SD[1], by=eid] 
          events[first_matches, on = .(eid, inci_censor_date > diagnosis_date, assessment_date < diagnosis_date),
            c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
            .(diagnosis_date, "common/Cancer Register/cancer_register.csv", "cr_id", cr_id)]
        }
      }
    }
    
    rm(cr_inci) # table no longer needed
  }
  
  if (need_cancer_register_prevalent) {
    vmessage("    Finding prevalent cancer register records matching endpoint definition...")
    
    # Create a copy before applying any filters. If no filters are applied, only
    # a pointer is stored, not a full additional copy of 'hes'
    cr_prev <- cr
    
    excl_opts <- c(
      "excluding icd-10", 
      "excluding prevalent icd-10" 
    )
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          cr_prev <- cr_prev[!(icd_version == 10 & diag_icd %like% code)]
        }
      }
    }
    
    excl_opts <- "excluding prevalent icd-9"
    for (opt_type in excl_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Removing codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Removing %s and sub-codes...", code))
          code <- format_code(code)
          cr_prev <- cr_prev[!(icd_version == 9 & diag_icd %like% code)] 
        }
      }
    }

    match_opts <- c("icd-10", "prevalent icd-10")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          if ("most recent prevalent event" %in% names(opts)) {
            most_recent <- cr_prev[icd_version == 10 & diag_icd %like% code][
              events, on = .(eid, diagnosis_date > min_censor_date, diagnosis_date < prev_censor_date),
              mult="last", nomatch=0, .(eid, visit_index=i.visit_index, cr_id, diagnosis_date=x.diagnosis_date)
            ]
            events[most_recent, on = .(eid, visit_index), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(diagnosis_date, "common/Cancer Register/cancer_register.csv", "cr_id", cr_id)]
          } else {
            first_matches <- cr_prev[icd_version == 10 & diag_icd %like% code, .SD[1], by=eid]
            events[first_matches, on = .(eid, min_censor_date < diagnosis_date, prev_censor_date > diagnosis_date),
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(diagnosis_date, "common/Cancer Register/cancer_register.csv", "cr_id", cr_id)]
          }
        }
      }
    }
    
    match_opts <- c("prevalent icd-9")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        codes <- opts[[opt_type]]
        for (code in codes) {
          vmessage(sprintf("        Matching %s and sub-codes...", code))
          code <- format_code(code)
          if ("most recent prevalent event" %in% names(opts)) {
            most_recent <- cr_prev[icd_version == 9 & diag_icd %like% code][
              events, on = .(eid, diagnosis_date > min_censor_date, diagnosis_date < prev_censor_date),
              mult="last", nomatch=0, .(eid, visit_index=i.visit_index, cr_id, diagnosis_date=x.diagnosis_date)
            ]
            events[most_recent, on = .(eid, visit_index), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(diagnosis_date, "common/Cancer Register/cancer_register.csv", "cr_id", cr_id)]
          } else {
            first_matches <- cr_prev[icd_version == 9 & diag_icd %like% code, .SD[1], by=eid]
            events[first_matches, on = .(eid, min_censor_date < diagnosis_date, prev_censor_date > diagnosis_date),
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(diagnosis_date, "common/Cancer Register/cancer_register.csv", "cr_id", cr_id)]
          }
        }
      }
    }
    
    rm(cr_prev); # table no longer needed
  }
  
  # Remove tables we no longer need to free up memory
  vmessage("  Unlooading cancer register data...")
  rm(cr)
  invisible(gc()) 
}

if (need_gp_records) {
  vmessage("  Loading primary care records...")
  gp <- load_from_rap("/common/Primary Care/gp_clinical_records.csv")
  
  # Filter to people who have consented for health record linkage (and to 
  # the sub-cohort with primary care data if that filter has been applied)
  vmessage("    Filtering to participants consented for electronic health record linkage...")
  gp <- gp[eid %in% followup$eid]
  
  if ("do not impute missing event dates" %in% names(opts)) {
    vmessage("    Setting imputed event dates to NA...")
    gp[(imputed_date), event_date := NA]
    if ("remove missing event dates" %in% names(opts)) {
      vmessage("    Removing events missing event dates...")
      gp <- gp[!is.na(event_date)]
    }
  }
  
  vmessage("    Sorting events...")
  gp <- gp[order(event_date)][order(eid)] # NA dates if present come last within eid group
  
  if (need_gp_records_incident) {
    vmessage("    Finding incident primary care records matching endpoint definition...")
    
    # First, handle Read code matching
    match_opts <- c("read-3", "incident read-3", "read-2", "incident read-2")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        covs <- opts[[opt_type]]
        codes <- covs$codes
        operators <- covs$value_comparison_operators
        values <- covs$values
        this_read_version <- as.integer(gsub("[a-z]| |-", "", opt_type))
        
        # Match each code or code-value combination
        for (ii in seq_along(codes)) {
          this_code <- codes[ii]
          if (is.null(operators)) {
            this_operator <- NA
            this_value <- NA
          } else {
            this_operator <- operators[ii]
            this_value <- values[ii]
          }
          
          if (is.na(this_operator)) {
            vmessage(sprintf("        Matching %s...", this_code))
          } else {
            vmessage(sprintf("        Matching %s %s %s...", this_code, this_operator, this_value))
          }
          
          if (is.na(this_operator)) {
            first_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code, .SD[1], by=eid]
          } else if (this_operator == "=") {
            first_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code & value == this_value, .SD[1], by=eid]
          } else if (this_operator == "!=") {
            first_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code & value != this_value, .SD[1], by=eid]
          } else if (this_operator == ">") {
            first_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code & value > this_value, .SD[1], by=eid]
          } else if (this_operator == ">=") {
            first_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code & value >= this_value, .SD[1], by=eid]
          } else if (this_operator == "<") {
            first_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code & value < this_value, .SD[1], by=eid]
          } else if (this_operator == "<=") {
            first_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code & value <= this_value, .SD[1], by=eid]
          } else {
            stop(sprintf("Unrecognised operator '%s' when matching primary care READ code values", this_operator))
          }
 
          events[first_matches, on = .(eid, inci_censor_date > event_date, assessment_date < event_date),
            c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
            .(event_date, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]            
          
          # Handle NA dates, if applicable
          if (is.na(this_operator)) {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code]
          } else if (this_operator == "=") {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code & value == this_value]
          } else if (this_operator == "!=") {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code & value != this_value]
          } else if (this_operator == ">") {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code & value > this_value]
          } else if (this_operator == ">=") {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code & value >= this_value]
          } else if (this_operator == "<") {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code & value < this_value]
          } else if (this_operator == "<=") {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code & value <= this_value]
          }
          events[na_matches, on = .(eid),
            c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
            .(NA, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
        }
      }
    }
    
    # Next, handle ICD-10 code matching (if relevant) 
    if (!("no icd-10 lookup in gp records" %in% names(options))) {
      # Create a copy before applying any filters. If no filters are applied, only
      # a pointer is stored, not a full additional copy of 'gp'
      gp_inci <- gp
      
      excl_opts <- c("excluding icd-10", "excluding incident icd-10")
      for (opt_type in excl_opts) {
        if (opt_type %in% names(opts)) {
          vmessage(sprintf("      Removing codes in '%s'...", opt_type))
          codes <- opts[[opt_type]]
          for (code in codes) {
            vmessage(sprintf("        Removing %s and sub-codes...", code))
            code <- format_code(code)
            gp_inci <- gp_inci[!(icd10_code %like% code)]
          }
        }
      }
      
      match_opts <- c("icd-10", "incident icd-10")
      for (opt_type in match_opts) {
        if (opt_type %in% names(opts)) {
          vmessage(sprintf("      Matching codes in '%s'...", opt_type))
          codes <- opts[[opt_type]]
          for (code in codes) {
            vmessage(sprintf("        Matching %s and sub-codes...", code))
            code <- format_code(code)
            
            first_matches <- gp_inci[!is.na(event_date) & icd10_code %like% code, .SD[1], by=eid] 
            events[first_matches, on = .(eid, inci_censor_date > event_date, assessment_date < event_date),
              c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
              .(event_date, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
            
            na_matches <- gp_inci[is.na(event_date) & icd10_code %like% code]
            events[na_matches, on = .(eid),
              c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
              .(NA, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
           
          }
        }
      }
    }
    
    # Next, handle OPCS-4 code matching (if relevant) 
    if (!("no opcs-4 lookup in gp records" %in% names(options))) {
      # Create a copy before applying any filters. If no filters are applied, only
      # a pointer is stored, not a full additional copy of 'gp'
      gp_inci <- gp
      
      excl_opts <- c("excluding opcs-4", "excluding incident opcs-4")
      for (opt_type in excl_opts) {
        if (opt_type %in% names(opts)) {
          vmessage(sprintf("      Removing codes in '%s'...", opt_type))
          codes <- opts[[opt_type]]
          for (code in codes) {
            vmessage(sprintf("        Removing %s and sub-codes...", code))
            code <- format_code(code)
            gp_inci <- gp_inci[!(opcs4_code %like% code)]
          }
        }
      }
      
      match_opts <- c("opcs-4", "incident opcs-4")
      for (opt_type in match_opts) {
        if (opt_type %in% names(opts)) {
          vmessage(sprintf("      Matching codes in '%s'...", opt_type))
          codes <- opts[[opt_type]]
          for (code in codes) {
            vmessage(sprintf("        Matching %s and sub-codes...", code))
            code <- format_code(code)
            
            first_matches <- gp_inci[!is.na(event_date) & opcs4_code %like% code, .SD[1], by=eid] 
            events[first_matches, on = .(eid, inci_censor_date > event_date, assessment_date < event_date),
              c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
              .(event_date, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
            
            na_matches <- gp_inci[is.na(event_date) & opcs4_code %like% code]
            events[na_matches, on = .(eid),
              c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
              .(NA, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
           
          }
        }
      }
    }
  }
    
  if (need_gp_records_prevalent) {
    vmessage("    Finding prevalent primary care records matching endpoint definition...")
    
    # First, handle Read code matching
    match_opts <- c("read-3", "prevalent read-3", "read-2", "prevalent read-2")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        covs <- opts[[opt_type]]
        codes <- covs$codes
        operators <- covs$value_comparison_operators
        values <- covs$values
        this_read_version <- as.integer(gsub("[a-z]| |-", "", opt_type))
        
        # Handle specific field=code combinations
        for (ii in seq_along(fields)) {
          this_code <- codes[ii]
          if (is.null(operators)) {
            this_operator <- NA
            this_value <- NA
          } else {
            this_operator <- operators[ii]
            this_value <- values[ii]
          }
          
          if (is.na(this_operator)) {
            vmessage(sprintf("        Matching %s...", this_code))
          } else {
            vmessage(sprintf("        Matching %s %s %s...", this_code, this_operator, this_value))
          }
          
          if (is.na(this_operator)) {
            all_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code]
          } else if (this_operator == "=") {
            all_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code & value == this_value]
          } else if (this_operator == "!=") {
            all_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code & value != this_value]
          } else if (this_operator == ">") {
            all_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code & value > this_value]
          } else if (this_operator == ">=") {
            all_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code & value >= this_value]
          } else if (this_operator == "<") {
            all_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code & value < this_value]
          } else if (this_operator == "<=") {
            all_matches <- gp[!is.na(event_date) & read_version == this_read_version & read_code == this_code & value <= this_value]
          } else {
            stop(sprintf("Unrecognised operator '%s' when matching primary care READ code values", this_operator))
          }
          
          if ("most recent prevalent event" %in% names(opts)) {
            most_recent <- all_matches[events, 
              on = .(eid, event_date > min_censor_date, event_date < prev_censor_date),
              mult="last", nomatch=0, .(eid, gp_id, event_date=x.event_date)]
            events[most_recent, on = .(eid), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(event_date, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
          } else {
            first_matches <- all_matches[events, 
              on = .(eid, event_date > min_censor_date, event_date < prev_censor_date),
              mult="first", nomatch=0, .(eid, sr_id, event_date=x.event_date)]
            events[first_matches, on = .(eid), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(event_date, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
          }
          
          # Handle NA dates, if applicable
          if (is.na(this_operator)) {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code]
          } else if (this_operator == "=") {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code & value == this_value]
          } else if (this_operator == "!=") {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code & value != this_value]
          } else if (this_operator == ">") {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code & value > this_value]
          } else if (this_operator == ">=") {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code & value >= this_value]
          } else if (this_operator == "<") {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code & value < this_value]
          } else if (this_operator == "<=") {
            na_matches <- gp[is.na(event_date) & read_version == this_read_version & read_code == this_code & value <= this_value]
          }
          
          events[na_matches, on = .(eid),
            c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
            .(NA, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
        }
      }
    }
    
    # Next, handle ICD-10 code matching (if relevant) 
    if (!("no icd-10 lookup in gp records" %in% names(options))) {
      # Create a copy before applying any filters. If no filters are applied, only
      # a pointer is stored, not a full additional copy of 'gp'
      gp_prev <- gp
      
      excl_opts <- c("excluding icd-10", "excluding prevalent icd-10")
      for (opt_type in excl_opts) {
        if (opt_type %in% names(opts)) {
          vmessage(sprintf("      Removing codes in '%s'...", opt_type))
          codes <- opts[[opt_type]]
          for (code in codes) {
            vmessage(sprintf("        Removing %s and sub-codes...", code))
            code <- format_code(code)
            gp_prev <- gp_prev[!(icd10_code %like% code)]
          }
        }
      }
      
      match_opts <- c("icd-10", "prevalent icd-10")
      for (opt_type in match_opts) {
        if (opt_type %in% names(opts)) {
          vmessage(sprintf("      Matching codes in '%s'...", opt_type))
          codes <- opts[[opt_type]]
          for (code in codes) {
            vmessage(sprintf("        Matching %s and sub-codes...", code))
            code <- format_code(code)
            
            if ("most recent prevalent event" %in% names(opts)) {
              most_recent <- gp_prev[!is.na(event_date) & icd10_code %like% code][
                events, on = .(eid, event_date > min_censor_date, event_date < prev_censor_date),
                mult="last", nomatch=0, .(eid, visit_index=i.visit_index, gp_id, event_date=x.event_date)
              ]
              events[most_recent, on = .(eid, visit_index), 
                c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
                .(event_date, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
            } else {
              first_matches <- gp_prev[!is.na(event_date) & icd10_code %like% code, .SD[1], by=eid]
              events[first_matches, on = .(eid, min_censor_date < event_date, prev_censor_date > event_date),
                c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
                .(event_date, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
            }
            
            na_matches <- gp_prev[is.na(event_date) & icd10_code %like% code]
            events[na_matches, on = .(eid),
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(NA, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
          }
        }
      }
    }
    
    # Next, handle OPCS-4 code matching (if relevant) 
    if (!("no opcs-4 lookup in gp records" %in% names(options))) {
      # Create a copy before applying any filters. If no filters are applied, only
      # a pointer is stored, not a full additional copy of 'gp'
      gp_prev <- gp
      
      excl_opts <- c("excluding opcs-4", "excluding prevalent opcs-4")
      for (opt_type in excl_opts) {
        if (opt_type %in% names(opts)) {
          vmessage(sprintf("      Removing codes in '%s'...", opt_type))
          codes <- opts[[opt_type]]
          for (code in codes) {
            vmessage(sprintf("        Removing %s and sub-codes...", code))
            code <- format_code(code)
            gp_prev <- gp_prev[!(opcs4_code %like% code)]
          }
        }
      }
      
      match_opts <- c("opcs-4", "prevalent opcs-4")
      for (opt_type in match_opts) {
        if (opt_type %in% names(opts)) {
          vmessage(sprintf("      Matching codes in '%s'...", opt_type))
          codes <- opts[[opt_type]]
          for (code in codes) {
            vmessage(sprintf("        Matching %s and sub-codes...", code))
            code <- format_code(code)
            
            if ("most recent prevalent event" %in% names(opts)) {
              most_recent <- gp_prev[!is.na(event_date) & opcs4_code %like% code][
                events, on = .(eid, event_date > min_censor_date, event_date < prev_censor_date),
                mult="last", nomatch=0, .(eid, visit_index=i.visit_index, gp_id, event_date=x.event_date)
              ]
              events[most_recent, on = .(eid, visit_index), 
                c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
                .(event_date, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
            } else {
              first_matches <- gp_prev[!is.na(event_date) & opcs4_code %like% code, .SD[1], by=eid]
              events[first_matches, on = .(eid, min_censor_date < event_date, prev_censor_date > event_date),
                c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
                .(event_date, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
            }
            
            na_matches <- gp_prev[is.na(event_date) & opcs4_code %like% code]
            events[na_matches, on = .(eid),
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(NA, "common/Primary Care/gp_clinical_records.csv", "gp_id", gp_id)]
          }
        }
      }
    }
  }
  
  # Remove tables we no longer need to free up memory
  vmessage("  Unlooading primary care records...")
  rm(gp)
  invisible(gc())
}

if (need_self_report) {
  vmessage("  Loading self-reported medical history...")
  sr <- load_from_rap("/common/Medical History/medical_history.csv")
  
  vmessage("    Dropping records that don't contain events (e.g. \"No\" answers to multiple-choice questions)...")
  sr <- sr[!(code == 0 & label != "0")] # '0' label indicates a continuous score answer
  sr <- sr[code != -2] # Not applicable - i.e. fields that rely on "Yes" answers to other fields.
  sr <- sr[code != -7] # None of the above - i.e. when asked about multiple disease options
  sr <- sr[!(code == 1 & label == "Not at all")] # Frequency of mental health episodes
  
  if (!("allow missing event status" %in% names(opts))) {
    vmessage("    Removing \"Do not know\" and \"Prefer not to answer\" records to prevent NAs in event column...")
    sr <- sr[!(code == -3 & label == "Prefer not to answer")]
    sr <- sr[!(code == -1 & label == "Do not know")]
  }
  
  if ("do not impute missing event dates" %in% names(opts)) {
    vmessage("    Setting imputed event dates to NA...")
    sr[(imputed), date := NA]
    if ("remove missing event dates" %in% names(opts)) {
      vmessage("    Removing events missing event dates...")
      sr <- sr[!is.na(date)]
    }
  }
  
  vmessage("    Sorting events...")
  sr <- sr[order(date)][order(eid)][order(field_id)] # NA dates if present come last within eid/field_id group

  if (need_self_report_incident) {
    vmessage("    Finding incident self-reported medical history records matching endpoint definition...")
    
    match_opts <- c("self-report", "incident self-report")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        focs <- opts[[opt_type]]
        fields <- focs$fields
        operators <- focs$field_code_comparison_operators
        codes <- focs$codes
        
        # Handle specific field=code combinations
        for (ii in seq_along(fields)) {
          this_field <- fields[ii]
          this_operator <- operators[ii]
          this_code <- codes[ii]
          vmessage(sprintf("        Matching field-id: %s %s %s...", this_field, this_operator, this_code))
          
          if (this_operator == "=") {
            first_matches <- sr[!is.na(date) & field_id == this_field & code == this_code, .SD[1], by=eid]
          } else if (this_operator == "!=") {
            first_matches <- sr[!is.na(date) & field_id == this_field & code != this_code, .SD[1], by=eid]
          } else if (this_operator == ">") {
            first_matches <- sr[!is.na(date) & field_id == this_field & code > this_code, .SD[1], by=eid]
          } else if (this_operator == ">=") {
            first_matches <- sr[!is.na(date) & field_id == this_field & code >= this_code, .SD[1], by=eid]
          } else if (this_operator == "<") {
            first_matches <- sr[!is.na(date) & field_id == this_field & code < this_code, .SD[1], by=eid]
          } else if (this_operator == "<=") {
            first_matches <- sr[!is.na(date) & field_id == this_field & code <= this_code, .SD[1], by=eid]
          } else {
            stop(sprintf("Unrecognised operator '%s' when matching self-report codes", this_operator))
          }
 
          events[first_matches, on = .(eid, inci_censor_date > date, assessment_date < date),
            c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
            .(date, "common/Medical History/medical_history.csv", "sr_id", sr_id)]            
          
          # Handle NA dates, if applicable
          if (this_operator == "=") {
            na_matches <- sr[is.na(date) & field_id == this_field & code == this_code]
          } else if (this_operator == "!=") {
            na_matches <- sr[is.na(date) & field_id == this_field & code != this_code]
          } else if (this_operator == ">") {
            na_matches <- sr[is.na(date) & field_id == this_field & code > this_code]
          } else if (this_operator == ">=") {
            na_matches <- sr[is.na(date) & field_id == this_field & code >= this_code]
          } else if (this_operator == "<") {
            na_matches <- sr[is.na(date) & field_id == this_field & code < this_code]
          } else if (this_operator == "<=") {
            na_matches <- sr[is.na(date) & field_id == this_field & code <= this_code]
          }
          na_matches <- na_matches[order(visit_index)]
          events[na_matches, on = .(eid, visit_index < visit_index), # match where the self-report event is reported at a later visit_index
            c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
            .(NA, "common/Medical History/medical_history.csv", "sr_id", sr_id)]
        }
        
        # Handle NA answers to relevant multiple choice questions that make 
        # absence of case status impossible to determine
        if ("allow missing event status" %in% names(opts)) {
          for (this_field in unique(fields)) {
            na_matches <- sr[field_id == this_field & (
              (code == -1 & label == "Do not know") | 
              (code == -3 & label == "Prefer not to answer")
            )]
            na_matches <- na_matches[order(visit_index)]
            events[na_matches, on = .(eid, visit_index < visit_index), # match where the self-report event is reported at a later visit_index
              c("inci_censor_date", "inci_censor_data_source", "inci_censor_id_column", "inci_censor_index") :=
              .(NA, "common/Medical History/medical_history.csv", "sr_id", sr_id)]
          }
        }
      }
    }
  }
  
  if (need_self_report_prevalent) {
    vmessage("    Finding prevalent self-reported medical history records matching endpoint definition...")
    
    match_opts <- c("self-report", "prevalent self-report")
    for (opt_type in match_opts) {
      if (opt_type %in% names(opts)) {
        vmessage(sprintf("      Matching codes in '%s'...", opt_type))
        focs <- opts[[opt_type]]
        fields <- focs$fields
        operators <- focs$field_code_comparison_operators
        codes <- focs$codes
        
        # Handle specific field=code combinations
        for (ii in seq_along(fields)) {
          this_field <- fields[ii]
          this_operator <- operators[ii]
          this_code <- codes[ii]
          vmessage(sprintf("        Matching field-id: %s %s %s...", this_field, this_operator, this_code))
          
          if (this_operator == "=") {
            all_matches <- sr[!is.na(date) & field_id == this_field & code == this_code]
          } else if (this_operator == "!=") {
            all_matches <- sr[!is.na(date) & field_id == this_field & code != this_code]
          } else if (this_operator == ">") {
            all_matches <- sr[!is.na(date) & field_id == this_field & code > this_code]
          } else if (this_operator == ">=") {
            all_matches <- sr[!is.na(date) & field_id == this_field & code >= this_code]
          } else if (this_operator == "<") {
            all_matches <- sr[!is.na(date) & field_id == this_field & code < this_code]
          } else if (this_operator == "<=") {
            all_matches <- sr[!is.na(date) & field_id == this_field & code <= this_code]
          } else {
            stop(sprintf("Unrecognised operator '%s' when matching self-report codes", this_operator))
          }
          
          if ("most recent prevalent event" %in% names(opts)) {
            if ("use self-report from current visit only" %in% names(opts)) {
              most_recent <- all_matches[events, 
                on = .(eid, visit_index, date > min_censor_date, date < prev_censor_date),
                mult="last", nomatch=0, .(eid, visit_index=i.visit_index, sr_id, date=x.date)]
              events[most_recent, on = .(eid, visit_index), 
                c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
                .(date, "common/Medical History/medical_history.csv", "sr_id", sr_id)]
            } else {
              most_recent <- all_matches[events, 
                on = .(eid, date > min_censor_date, date < prev_censor_date),
                mult="last", nomatch=0, .(eid, visit_index=i.visit_index, sr_id, date=x.date)]
              events[most_recent, on = .(eid, visit_index), 
                c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
                .(date, "common/Medical History/medical_history.csv", "sr_id", sr_id)]
            }
          } else {
            if ("use self-report from current visit only" %in% names(opts)) {
              first_matches <- all_matches[events, 
                on = .(eid, visit_index, date > min_censor_date, date < prev_censor_date),
                mult="first", nomatch=0, .(eid, visit_index=i.visit_index, sr_id, date=x.date)]
              events[first_matches, on = .(eid, visit_index), 
                c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
                .(date, "common/Medical History/medical_history.csv", "sr_id", sr_id)]
            } else {
              first_matches <- all_matches[events, 
                on = .(eid, date > min_censor_date, date < prev_censor_date),
                mult="first", nomatch=0, .(eid, visit_index=i.visit_index, sr_id, date=x.date)]
              events[first_matches, on = .(eid, visit_index), 
                c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
                .(date, "common/Medical History/medical_history.csv", "sr_id", sr_id)]      
            }
          }
          
          # Handle NA dates, if applicable
          if (this_operator == "=") {
            na_matches <- sr[is.na(date) & field_id == this_field & code == this_code]
          } else if (this_operator == "!=") {
            na_matches <- sr[is.na(date) & field_id == this_field & code != this_code]
          } else if (this_operator == ">") {
            na_matches <- sr[is.na(date) & field_id == this_field & code > this_code]
          } else if (this_operator == ">=") {
            na_matches <- sr[is.na(date) & field_id == this_field & code >= this_code]
          } else if (this_operator == "<") {
            na_matches <- sr[is.na(date) & field_id == this_field & code < this_code]
          } else if (this_operator == "<=") {
            na_matches <- sr[is.na(date) & field_id == this_field & code <= this_code]
          }
          na_matches <- na_matches[order(visit_index)]
          
          if ("use self-report from current visit only" %in% names(opts)) {
            events[na_matches, on = .(eid, visit_index), 
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(NA, "common/Medical History/medical_history.csv", "sr_id", sr_id)]
          } else {
            events[na_matches, on = .(eid, visit_index >= visit_index), # match where the self-report event is reported at the current or earlier visit_index
              c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
              .(NA, "common/Medical History/medical_history.csv", "sr_id", sr_id)]
          }
        }
        
        # Handle NA answers to relevant multiple choice questions that make 
        # absence of case status impossible to determine
        if ("allow missing event status" %in% names(opts)) {
          for (this_field in unique(fields)) {
            na_matches <- sr[field_id == this_field & (
              (code == -1 & label == "Do not know") | 
                (code == -3 & label == "Prefer not to answer")
            )]
            na_matches <- na_matches[order(visit_index)]
            
            if ("use self-report from current visit only" %in% names(opts)) {
              events[na_matches, on = .(eid, visit_index), 
                c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
                .(NA, "common/Medical History/medical_history.csv", "sr_id", sr_id)]
            } else {
              events[na_matches, on = .(eid, visit_index >= visit_index), # match where the self-report event is reported at the current or earlier visit_index
                c("prev_censor_date", "prev_censor_data_source", "prev_censor_id_column", "prev_censor_index") :=
                .(NA, "common/Medical History/medical_history.csv", "sr_id", sr_id)]
            }
          }
        }
      }
    }
  }
  
  # Remove tables we no longer need to free up memory
  vmessage("  Unloading self-reported medical history...")
  rm(sr)
  invisible(gc())
}

# Sanity checks - prevalent events should not occur after assessment and vice-versa for incident
# if these errors are thrown, debug the relevant joins above
if (need_prevalent) stopifnot(events[!is.na(prev_censor_index) & prev_censor_date > assessment_date, .N] == 0)
if (need_incident) stopifnot(events[!is.na(inci_censor_index) & inci_censor_date < assessment_date, .N] == 0)

########################################################################################
# Collapse into a single traceable event and follow-up column
########################################################################################

vmessage("Collapsing into a single traceable event and follow-up column...")

ro <- events[,.(eid, visit_index)]

if (need_incident && need_prevalent) {
  events <- rbind(
    # incident events
    events[!is.na(inci_censor_index) & is.na(prev_censor_index), .(
      eid, visit_index, assessment_date, event=TRUE, event_date=inci_censor_date, 
      follow_up=NA_real_, data_source=inci_censor_data_source, 
      data_source_lookup_column=inci_censor_id_column,
      data_source_lookup_id=inci_censor_index,
      min_censor_date)],
    # Prevalent events
    events[!is.na(prev_censor_index), .(
      eid, visit_index, assessment_date, event=TRUE, event_date=prev_censor_date, 
      follow_up=NA_real_, data_source=prev_censor_data_source, 
      data_source_lookup_column=prev_censor_id_column,
      data_source_lookup_id=prev_censor_index,
      min_censor_date)],
    # No events
    events[is.na(inci_censor_index) & is.na(prev_censor_index), .(
      eid, visit_index, assessment_date, event=FALSE, event_date=max_censor_date, 
      follow_up=NA_real_, data_source=max_censor_data_source, 
      data_source_lookup_column=max_censor_id_column,
      data_source_lookup_id=max_censor_reason,
      min_censor_date)]
  ) 
} else if (need_incident && !need_prevalent) {
  events <- rbind(
    # incident events
    events[!is.na(inci_censor_index), .(
      eid, visit_index, assessment_date, event=TRUE, event_date=inci_censor_date, 
      follow_up=NA_real_, data_source=inci_censor_data_source, 
      data_source_lookup_column=inci_censor_id_column,
      data_source_lookup_id=inci_censor_index, min_censor_date)],
    # No events
    events[is.na(inci_censor_index), .(
      eid, visit_index, assessment_date, event=FALSE, event_date=max_censor_date, 
      follow_up=NA_real_, data_source=max_censor_data_source, 
      data_source_lookup_column=max_censor_id_column,
      data_source_lookup_id=max_censor_reason, min_censor_date)]
  ) 
} else if (need_prevalent && !need_incident) {
  events <- rbind(
    # Prevalent events
    events[!is.na(prev_censor_index), .(
      eid, visit_index, assessment_date, event=TRUE, event_date=prev_censor_date, 
      follow_up=NA_real_, data_source=prev_censor_data_source, 
      data_source_lookup_column=prev_censor_id_column,
      data_source_lookup_id=prev_censor_index, min_censor_date)],
    # No events
    events[is.na(prev_censor_index), .(
      eid, visit_index, assessment_date, event=FALSE, event_date=min_censor_date, 
      follow_up=NA_real_, data_source=min_censor_data_source, 
      data_source_lookup_column=min_censor_id_column,
      data_source_lookup_id=min_censor_reason, min_censor_date)]
  ) 
}

# Sanity check that we haven't lost any events
stopifnot(events[,.N] == ro[,.N])
stopifnot(events[,.N,by=.(eid, visit_index)][N > 1, .N] == 0)

# Reorder rows
events <- events[ro, on = .(eid, visit_index)]

# Drop people who haven't consented for health record linkage and have no
# events from the self-reported medical history
events <- events[!is.na(data_source)]

########################################################################################
# Get mid-point dates for computing follow-up as midpoint
########################################################################################

if ("impute follow-up as midpoint" %in% names(opts)) {
  vmessage("Finding closest healthy events required to impute follow-up as midpoint...")
  need_latest_healthy=FALSE
  need_earliest_healthy=FALSE
  if (!("most recent prevalent event" %in% names(opts))) {
    first_inci <- events[(event) & event_date > assessment_date, .(eid, visit_index, latest_healthy=assessment_date, event_date)]
    first_prev <- events[(event) & event_date < assessment_date, .(eid, visit_index, latest_healthy=min_censor_date, event_date)]
    first <- rbind(first_inci, first_prev)
    need_latest_healthy=TRUE
    case_eids <- first$eid
  } else if (!need_incident) {
    last <- events[(event) & event_date < assessment_date, .(eid, visit_index, event_date, earliest_healthy=assessment_date)]
    need_earliest_healthy=TRUE
    case_eids <- last$eid
  } else {
    first <- events[(event) & event_date > assessment_date, .(eid, visit_index, latest_healthy=assessment_date, event_date)]
    last <- events[(event) & event_date < assessment_date, .(eid, visit_index, event_date, earliest_healthy=assessment_date)]
    need_latest_healthy=TRUE
    need_earliest_healthy=TRUE
    case_eids <- c(first$eid, last$eid)
  }
  
  # N.b. death records are not considered here for obvious reasons
  
  if (need_hospital_records) {
    vmessage("  Loading hospital records diagnosis data...")
    hes <- load_from_rap("common/Hospital Records/diagnoses.csv")
    
    # Filter to people who have consented for health record linkage (and to 
    # the sub-cohort with primary care data if that filter has been applied)
    vmessage("    Filtering to cases...")
    hes <- hes[eid %in% case_eids]
    
    vmessage("    Sorting events...")
    hes <- hes[order(event_start)][order(eid)]
    
    # Find records that occur after the most recent known healthy date and before
    # the event date
    if (need_latest_healthy) {
      vmessage("    Finding latest records prior to first events...")
      latest <- hes[first, on = .(eid, event_start > latest_healthy, event_start < event_date),
        mult="last", nomatch=0, .(eid, visit_index=i.visit_index, healthy_date=x.event_start)]
      first[latest, on = .(eid, visit_index), latest_healthy := i.healthy_date]
    }
    # Find records that occur before the most recent known healthy date and after
    # the event date (i.e. where 'most recent prevalent event' is set)
    if (need_earliest_healthy) {
      vmessage("    Finding earliest records after most recent prevalent event...")
      earliest <- hes[last, on = .(eid, event_start > event_date, event_start < earliest_healthy),
        mult="first", nomatch=0, .(eid, visit_index=i.visit_index, healthy_date=x.event_start)]
      last[earliest, on = .(eid, visit_index), earliest_healthy := i.healthy_date]
    }
    
    # Remove tables we no longer need to free up memory
    vmessage("  Unloading hospital records diagnosis data...")
    rm(hes)
    invisible(gc())
  }
  
  if (need_hospital_operations) {
    vmessage("  Loading hospital records operations and procedures data...")
    hop <- load_from_rap("common/Hospital Records/operations.csv")
    
    # Filter to people who have consented for health record linkage (and to 
    # the sub-cohort with primary care data if that filter has been applied)
    vmessage("    Filtering to cases...")
    hop <- hop[eid %in% case_eids]
    
    vmessage("    Sorting events...")
    hop <- hop[order(imputed_operation_date)][order(eid)]
    
    # Find records that occur after the most recent known healthy date and before
    # the event date
    if (need_latest_healthy) {
      vmessage("    Finding latest records prior to first events...")
      latest <- hop[first, on = .(eid, imputed_operation_date > latest_healthy, imputed_operation_date < event_date),
        mult="last", nomatch=0, .(eid, visit_index=i.visit_index, healthy_date=x.imputed_operation_date)]
      first[latest, on = .(eid, visit_index), latest_healthy := i.healthy_date]
    }
    # Find records that occur before the most recent known healthy date and after
    # the event date (i.e. where 'most recent prevalent event' is set)
    if (need_earliest_healthy) {
      vmessage("    Finding earliest records after most recent prevalent event...")
      earliest <- hop[last, on = .(eid, imputed_operation_date > event_date, imputed_operation_date < earliest_healthy),
        mult="first", nomatch=0, .(eid, visit_index=i.visit_index, healthy_date=x.imputed_operation_date)]
      last[earliest, on = .(eid, visit_index), earliest_healthy := i.healthy_date]
    }
    
    # Remove tables we no longer need to free up memory
    vmessage("  Unloading hospital records operations and procedures data...")
    rm(hop)
    invisible(gc())
  }
  
  if (need_cancer_register) {
    vmessage("  Loading cancer register data...")
    cr <- load_from_rap("common/Cancer Register/cancer_register.csv")
    cr <- cr[!is.na(diagnosis_date)] # All are ICD-9 diagnoses with empty code strings 
    
    # Filter to people who have consented for health record linkage (and to 
    # the sub-cohort with primary care data if that filter has been applied)
    vmessage("    Filtering to cases...")
    cr <- cr[eid %in% case_eids]
    
    vmessage("    Sorting events...")
    cr <- cr[order(diagnosis_date)][order(eid)]
    
    # Find records that occur after the most recent known healthy date and before
    # the event date
    if (need_latest_healthy) {
      vmessage("    Finding latest records prior to first events...")
      latest <- cr[first, on = .(eid, diagnosis_date > latest_healthy, diagnosis_date < event_date),
        mult="last", nomatch=0, .(eid, visit_index=i.visit_index, healthy_date=x.diagnosis_date)]
      first[latest, on = .(eid, visit_index), latest_healthy := i.healthy_date]
    }
    # Find records that occur before the most recent known healthy date and after
    # the event date (i.e. where 'most recent prevalent event' is set)
    if (need_earliest_healthy) {
      vmessage("    Finding earliest records after most recent prevalent event...")
      earliest <- cr[last, on = .(eid, diagnosis_date > event_date, diagnosis_date < earliest_healthy),
        mult="first", nomatch=0, .(eid, visit_index=i.visit_index, healthy_date=x.diagnosis_date)]
      last[earliest, on = .(eid, visit_index), earliest_healthy := i.healthy_date]
    }
    
    # Remove tables we no longer need to free up memory
    vmessage("  Unloading cancer register data...")
    rm(cr)
    invisible(gc())
  }
  
  if (need_gp_records) {
    vmessage("  Loading primary care records...")
    gp <- load_from_rap("/common/Primary Care/gp_clinical_records.csv")
    
    # Filter to people who have consented for health record linkage (and to 
    # the sub-cohort with primary care data if that filter has been applied)
    vmessage("    Filtering to cases...")
    gp <- gp[eid %in% case_eids]
    
    if ("do not impute missing event dates" %in% names(opts)) {
      vmessage("    Setting imputed event dates to NA...")
      gp[(imputed_date), event_date := NA]
      if ("remove missing event dates" %in% names(opts)) {
        vmessage("    Removing events missing event dates...")
        gp <- gp[!is.na(event_date)]
      }
    }
    
    vmessage("    Sorting events...")
    gp <- gp[order(event_date)][order(eid)] # NA dates if present come last within eid group
    
    # Find records that occur after the most recent known healthy date and before
    # the event date
    if (need_latest_healthy) {
      vmessage("    Finding latest records prior to first events...")
      latest <- gp[first, on = .(eid, event_date > latest_healthy, event_date < event_date),
        mult="last", nomatch=0, .(eid, visit_index=i.visit_index, healthy_date=x.event_date)]
      first[latest, on = .(eid, visit_index), latest_healthy := i.healthy_date]
      
      na_matches <- gp[is.na(event_date)]
      first[na_matches, on = .(eid), latest_healthy := NA]
    }
    # Find records that occur before the most recent known healthy date and after
    # the event date (i.e. where 'most recent prevalent event' is set)
    if (need_earliest_healthy) {
      vmessage("    Finding earliest records after most recent prevalent event...")
      earliest <- gp[last, on = .(eid, event_date > event_date, event_date < earliest_healthy),
        mult="first", nomatch=0, .(eid, visit_index=i.visit_index, healthy_date=x.event_date)]
      last[earliest, on = .(eid, visit_index), earliest_healthy := i.healthy_date]
      
      na_matches <- gp[is.na(event_date)]
      last[na_matches, on = .(eid), earliest_healthy := NA]
    }
    
    # Remove tables we no longer need to free up memory
    vmessage("  Unloading primary care records...")
    rm(gp)
    invisible(gc())
  }
  
  # N.b. self-reported interpolated ages of events are not considered here
  
  vmessage("  Calculating midpoint dates...")
  if (need_latest_healthy) {
    first[, midpoint := midpoint(latest_healthy, event_date)]
    events[first, on = .(eid, visit_index), event_date := i.midpoint]
  }
  if (need_earliest_healthy) {
    last[, midpoint := midpoint(event_date, earliest_healthy)]
    events[last, on = .(eid, visit_index), event_date := i.midpoint]
  }
}

########################################################################################
# Compute follow-up time and clean up output
########################################################################################

vmessage("Computing follow-up time...")
events[, follow_up := years_between(assessment_date, event_date)]

vmessage("Curating output columns...")
setnames(events, "data_source_lookup_id", "reason_or_id")
setnames(events, "data_source_lookup_column", "lookup_column")

events <- events[,.(eid, visit_index, assessment_date, event, follow_up, censor_date=event_date, reason_or_id, data_source, lookup_column)]

if (!("time since prevalent event" %in% names(opts))) {
  events[censor_date <= assessment_date, follow_up := NA]
  if (!need_incident) {
    events[, follow_up := NULL]
  }
}
if (!("allow missing event status" %in% names(opts))) {
  events[!(event) & assessment_date > censor_date, 
    c("reason_or_id", "follow_up") := .("Assessment after maximum censor date", 0)]
} else {
  events[!(event) & assessment_date > censor_date, 
    c("reason_or_id", "event") := .("Assessment after maximum censor date", NA)]
}
if ("follow-up from birth" %in% names(opts)) {
  events[, visit_index := NULL]
}

########################################################################################
# Write out file to target destination
########################################################################################

vmessage("Writing output...")

if (args[["--local-output"]] || args[["--container-mode"]]) {
  system(sprintf("mkdir -p '%s'", dirname(args[["--output"]])))
  fwrite(events, file=args[["--output"]])
} else {
  fwrite(events, sprintf("%s/%s", work_dir, basename(args[["--output"]])))
  system(sprintf("dx mkdir -p '%s'", dirname(args[["--output"]])), ignore.stdout=!args[["--verbose"]])
  system(sprintf("dx upload '%s/%s' --destination '%s'/", work_dir, basename(args[["--output"]]), dirname(args[["--output"]])), ignore.stdout=!args[["--verbose"]])
  if (args[["--local-input"]] & !(basename(args[["--def-file"]]) %in% system(sprintf("dx ls '%s'", dirname(args[["--output"]])), intern=TRUE))) {
    system(sprintf("dx upload '%s' --destination '%s'/", def_file, dirname(args[["--output"]])), ignore.stdout=!args[["--verbose"]])
  }
}

vmessage("Done.")

