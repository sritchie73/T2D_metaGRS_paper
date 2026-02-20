Endpoints
================================================================================

This folder contains case/control status for various disease endpoints that have
been extracted from linked electronic health record data.

Some key folders:

 - 'CEU curated/': case/control status and time to first event (or disease free
   maximum follow-up for controls) for various cardiovascular endpoints that 
   have been defined based on ICD code lists defined by the CEU Epidemiology 
   working group, i.e. matching the harmonized endpoints that have been defined
   in other cohorts. Follow-up time has been computed relative to each date of 
   UK Biobank assessment. 
   
 - 'Phecodes/': case/control status and time to first event for phecodes version
   X (https://academic.oup.com/bioinformatics/article/39/11/btad655/7335839). 
   For each phecode, four files have been extracted: two based on events in the
   linked hospital episode statistics, and two with case numbers supplemented by
   events in the linked primary care records, which is available for only ~45% 
   of participants. For each of these, there are two files: one with the age at
   first event, and one with follow-up time computed relative to the date of 
   each UK Biobank assessment.
   
 - 'Diabetes (Eastwood Algorithm)'/: diabetes status and follow-up time 
   adjudicated from self-reported medical history, medication status, and linked
   hospital episode statistics following the algorithm described by Eastwood et 
   al. 2016 (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0162388)
 
 - 'UKB curated/': date of first events for UK Biobank's algorithmically defined
   ouctomes (https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=42)
   
As with other files under common/, rows in each dataset are uniquely 
identifiable based on the combination of 'eid' and 'visit_index' columns, which
give the application-specific participant identifier and the UK Biobank 
assessment visit index (0 for baseline assessment, 1 for first repeat 
assessment, 2 for first imaging assessment, 3 for second repeat assessment).

For underlying data, see other folders in common/:

 - 'common/Follow-up/'
 - 'common/Hospital Records/'
 - 'commom/Deaths/'
 - 'common/Cancer Register/'
 - 'common/Primary Care/'
 - 'common/Medical History/'
 - 'common/Medications/'
 
================================================================================
Curation of new endpoints (curate_endpoint.R)
================================================================================

This folder also contains a program, curate_endpoints.R, which can be used to 
determine case/control status and compute follow-up time for new endpoints based
on user-specific combinations of ICD codes, OPCS codes, Read codes, and 
self-report medical history codes, along with various options to control which
data sources are included and how follow-up time is computed.

The overall guiding principle that curate_endpoints.R works on is to try to 
maximise the case numbers across all possible data sources and to maximise the
follow-up time. In practice, this means that data sources that cover only a 
subset of the cohort (e.g. the primary care data, available only in ~45% of 
participants) will supplement the case numbers, rather than restricting the 
cohort size, and the maximum censor date will typically be that of the hospital
records, rather than the data source with the earliest maximum follow-up time.
By default, missing event dates, which applies to many self-reported disease 
events and a handful of primary care records, are imputed (see the respective
data source README files) rather than leading to NA in the follow-up time 
computation. All of these behaviours can be changed with specific options in the
endpoint definition file.

--------------------------------------------------------------------------------
Running curate_endpoint.R
--------------------------------------------------------------------------------

There are several ways you can run the program:

 (1) By clicking on the "curate_endpoint" applet to launch a job via a point 
     and click interface.
     
 (2) As a stand-alone program run on an interactive compute instance (e.g. 
     Posit/Rstudio)
     
 (3) As a part of a job submitted via 'dx run swiss-army-knife'.

In cases 2 and 3 the the program expects the compute instance to be attached to 
the CEU_overarching project.

When run as an applet, you will be asked to supply an endpoint definition file,
a name for the CSV file to save the output to, and the location to save the 
output CSV file.

To run curate_endpoints.R in an interative session, download the file to your 
local compute instance using 'dx download', then run the program from the
terminal window by running:

Rscript curate_endpoints.R --def-file path/to/def_file.R --output path/to/output.csv --verbose

Or from within R by running:

system("Rscript curate_endpoints.R --def-file path/to/def_file.R --output path/to/output.csv --verbose") 

The program requires two arguments: an endpoint definition file (--def-file) 
which contains the list of ICD codes and endpoint definition options (detailed
in a later section), and a file path to save the resulting output. 

The --verbose flag is optional, and tells the program to print out status 
updates as it goes.

By default, the program will look for the endpoint definition file at the given 
path on the CEU_overarching project storage, and save the output file at the 
given path on the CEU_overarching project storage.

To run the program instead on an endpoint defintion file located on your local
compute instance, add the '--local-input' flag, and the save the output file
to local compute without uploading to project storage, add the '--local-output'
flag.

When running on a local compute instance, the program will create a folder
named 'curate_endpoints_workdir/' on the local compute instance which it uses 
to download copies of the underlying data sources and endpoint definition file
needed to curate the endpoint, as well as save the output file for upload to
project storage. To change this folder you can use the --work-dir argument to
give an alternate file path.

To run the program as part of a job using swiss-army-knife, you need to
add the --container-mode flag. Note in this case, the --local-input and 
--local-output are not available.

A template job submission script you can use:

#!/bin/bash

cmd_to_run="Rscript /mnt/project/common/Endpoints/curate_endpoint.R \
  --def-file path/to/endpoint_defintion.txt \
  --container-mode --verbose \
  --output events_and_followup.csv"

dx run swiss-army-knife \
  --icmd $cmd_to_run \
  --destination 'path/to/output/folder/' \
  --instance-type 'mem3_ssd1_v2_x8' \
  --brief --yes 

--------------------------------------------------------------------------------
The endpoint definition file - part 1:
Defining the disease based on ICD, OPCS, Read, and self-report field codes
--------------------------------------------------------------------------------

The curate_endpoint.R takes as input an "endpoint definition file". This is a 
plain text file that contains one or more list of "codes" defining the disease
endpoint of interest, alongside additional "options" that control how follow-up
time is calculated, what data sources are included, and how missing data is
handled (see part 2 below).

A very basic input file looks like this:

-----
# Acute myocardial infarction
ICD-10: I21
-----

This is a file with two lines: one, starting with a #, which is treated as like
a comment in an R script, and ignored by curate_endpoint.R, and a second line
telling curate_endpoint.R to find the first instance of ICD-10 code I21 in the
hospital records, death records, and primary care data for each UK Biobank 
participant.

Endpoints can be defined based on multiple codes, which can be given as a comma
separate list, as a range:

-----
# Coronary heart disease
ICD-10: I20, I21, I22, I23, I24, I25
-----

Or 

-----
# Coronary heart disease
ICD-10: I20-I25
-----

Codes and code lists can also be given on multiple lines to enhance readability
ad documentation:

-----
# Coronary heart disease
ICD-10: I20 # Unstable angina
ICD-10: I21 # Acute myocardial infarction
ICD-10: I22, I23 # You can also mix and match codes, code lists, and code ranges
ICD-10: I24-I25
-----

When specifying codes, curate_endpoint.R will also match all subcodes as well

---
# Acute MI
ICD-10: I21 # Matches I21, I21.0, I21.1, ..., I21.9
----

You can also exclude specific ICD codes when specifying ranges:

----
# Equivalent to I20, I22-I25
ICD-10: I20-I25
Excluding ICD-10: I21
----

There are several other keywords that modify how a code or set of codes is 
looked up in the data: 

 - "Prevalent":            Match the code only where it occurs prior to the 
                           respective assessment time point.
 - "Incident":             Match the code only where it occurs after the 
                           respective assessment time point.
 - "Primary cause only":   Match the code only where it is the primary cause for
                           diagnosis (applies to hospital and death records).
 - "Fatal":                Match the code only where it occurs in the death 
                           records.
 - "Non-fatal:"            Do not look for the code in the death records.
     
For ICD-10 codes, all of these options can be combined in the following ways in
an endpoint definition file:

ICD-10:
Excluding ICD-10:
ICD-10 primary cause only:
Excluding ICD-10 primary cause only:
Prevalent ICD-10:
Excluding prevalent ICD-10:
Prevalent ICD-10 primary cause only:
Excluding prevalent ICD-10 primary cause only:
Incident ICD-10:
Excluding incident ICD-10:
Incident ICD-10 primary cause only:
Excluding incident ICD-10 primary cause only:
Non-fatal incident ICD-10:
Excluding non-fatal incident ICD-10:
Non-fatal incident ICD-10 primary cause only:
Excluding non-fatal incident ICD-10 primary cause only:
Fatal incident ICD-10:
Excluding fatal incident ICD-10:
Fatal incident ICD-10 primary cause only:
Excluding fatal incident ICD-10 primary cause only:

When matching ICD-10 codes, curate_endpoint.R will look in the hospital records,
death records, cancer register, and primary care data. Note only the hospital 
and death records distinguish between primary and secondary causes, so adding 
the "primary cause only" suffix will confine search for those codes to the 
hospital and death records.

For prevalent events, hospital records and cancer records may also be coded with
ICD-9 codes instead of ICD-10 codes. Hospitals in England switched from ICD-9 to 
ICD-10 codes in ~1997, hospitals in Wales switched in ~1999, and hospitals in 
Scotland switched in ~1996. The switch from ICD-9 to ICD-10 codes in the cancer
registers took place in ~1995 (England and Wales) and ~1997 (Scotland).
https://biobank.ndph.ox.ac.uk/showcase/exinfo.cgi?src=Data_providers_and_dates

ICD-9 codes can be supplied in the endpoint definition file using:

Prevalent ICD-9:
Excluding prevalent ICD-9:
Prevalent ICD-9 primary cause only:
Excluding prevalent ICD-9 primary cause only:

Operations and surgical procedures are coded with OPCS-4 codes, or OPCS-3 codes 
for operations taking place prior to ~1997 (hospitals in England), ~1999 
(hospitals in Wales), or ~1996 (hospitals in Scotland). When looking up OPCS-4
codes curate_endpoint.R will look in the hospital records and primary care data.

OPCS-4 codes can be supplied in the endpoint definition file using:

OPCS-4:
Excluding OPCS-4:
OPCS-4 primary cause only:
Excluding OPCS-4 primary cause only:
Prevalent OPCS-4:
Excluding prevalent OPCS-4:
Prevalent OPCS-4 primary cause only:
Excluding prevalent OPCS-4 primary cause only:
Incident OPCS-4:
Excluding incident OPCS-4:
Incident OPCS-4 primary cause only:
Excluding incident OPCS-4 primary cause only:

OPCS-3 codes can be supplied in the endpoint definition file using:

Prevalent OPCS-3:
Excluding prevalent OPCS-3:
Prevalent OPCS-3 primary cause only:
Excluding prevalent OPCS-3 primary cause only:

Primary care records are coded using a combination of version 2 or version 3
Read codes, depending on the general practice recording the event. Read codes
can be supplied in the endpoint definition file using:

Read-2:
Incident Read-2:
Prevalent Read-2:

Read-3:
Incident Read-3:
Prevalent Read-3:

Unlike ICD and OPCS codes, Read codes are matched exactly rather than using 
partial matching, as it is not clear that the read codes inherently have a 
hierarchical structure, so it is not possible to give a range of Read codes.

Some Read codes have values attached to them, in which case it is possible to
filter events to those matching a specific value, or exceeding a specific 
threshold:

----
Read-3: X772q > 6.5 # HbA1c level > 6.5% (warning: field may be a mix of % and mmol/L)
----

When filtering Read codes by value, the following operators can be used: =, !=,
>=, >, <=, <. 

Mappings between Read codes and ICD-10 codes for disease diagnoses and OPCS-4 
codes for surgical procedures have also been curated by UK Biobank. By default,
curate_endpoints.R will look up ICD-10 codes and OPCS-4 codes in the primary 
care records (available for ~45% of participants) to supplement the case
numbers, although this behaviour can be disabled using specific options detailed
below.

Self-reported medical history can be integrated into an endpoint by providing
a combination of UK Biobank field IDs and field-specific disease codes:

----
Prevalent self-report: 20002=1075 # Self-reported heart attack in verbal interview with nurse
----

A full list of self-report fields and associated codes and their meanings that
can be provided can be found in '/common/Medical History/code_labels.csv'. In
general all verbal interview medical questionnaires and all medical and mental
health questions in the touchscreen survey are covered. 

Self-reported medical history fields and codes can be specified using:

Self-report:
Incident self-report:
Prevalent self-report:

Touchscreen questionnaire fields that include scores, e.g. the neuroticism score
in field 20127, can also be filtered to values passing a specific threshold 
using the operators =, !=, >=, >, <=, or <.

--------------------------------------------------------------------------------
The endpoint definition file - part 2:
Options that change the behaviour of curate_endpoint.R file
--------------------------------------------------------------------------------

Options are given on separate lines of the endpoint definition file to change
the default behaviour of curate_endpoint.R

filter by sex: [males/females]
-------------------------------
Subset the cohort to either males or females before curating the endpoint by
adding as a line to the endpoint definition file either "filter by sex: males"
or "filter by sex: females".

max follow years: <N>
----------------------
Restrict the maximum follow-up time for each participant to <N> years rather 
than the maximum available follow-up. E.g. adding "max follow years: 10" to
the endpoint definition file will restrict the maximum follow-up time for both
case identification and disease-free survival for controls to 10 years.

max follow date: <YYYY-MM-DD>
------------------------------
Restrict the maximum follow-up to to a specific date,. E.g. adding 
"max follow date: 2020-01-15" to the endpoint definition file will restrict the
maximum follow-up date to the 15th of January 2020.

max follow age: <N>
--------------------
Restrict the maximum follow-up for each participant to <N> years of age. 

min follow years: <N>
----------------------
Restrict the retrospecitve follow-up time for prevalent case identification to
N years.

min follow date: <YYYY-MM-DD>
------------------------------
Restrict the retrospecitve follow-up time for prevalent case identification to
a specific date.

min follow age: <N>
--------------------
Restrict the retrospecitve follow-up time for each participant to a minimum of
N years of age.

use conservative max censor date
---------------------------------
Adding "use conservative max censor date" as a line to the endpoint definition
file will set the maximum follow-up time to the earliest maximum follow-up time
across all contributing data sources instead of treating these data sources as
supplemental sources for case identification. E.g. if the endpoint contains an
ICD-10 code for a cancer, the maximum follow-up will match the maximum follow
up in the cancer register, which is earlier than the maximum follow-up in the
hospital records.

no icd-10 lookup in gp records 
----------------------------
Adding "no icd-10 lookup in gp records" as a line to the endpoint definition 
file tells curate_endpoints.R not to look in the primary care data when 
identifying events based on ICD-10 codes. Useful if your endpoint has a 
combination of ICD-10 codes and Read-2 or Read-3 codes, and you specifically 
want to use Read codes only when identifying cases in the primary care data.

no opcs-4 lookup in gp records
-----------------------------
Adding "no opcs-4 lookup in gp records" as a line to the endpoint definition 
file tells curate_endpoints.R not to look in the primary care data when 
identifying events based on OPCS-4 codes. Useful if your endpoint has a 
combination of OPCS-4 codes and Read-2 or Read-3 codes, and you specifically 
want to use Read codes only when identifying cases in the primary care data.

exclude gp records 
-------------------
Adding "exclude gp records" as a line to the endpoint definition file tells 
curate_endpoints.R to ignore the primary care data entirely when curating the 
endpoint.

filter to primary care sub-cohort
----------------------------------
Adding "filter to primary care sub-cohort" tells curate_endpoint.R to subset the
cohort to only participants with linked primary care records before curating the
endpoint.

exclude death register 
-----------------------
Adding "exclude death register" as a line to the endpoint definition file tells 
curate_endpoints.R to ignore the death register data entirely when curating the
endpont.

exclude hospital records 
-------------------------
Adding "exclude hospital records" as a line to the endpoint definition file 
tells curate_endpoints.R to ignore the hospital records entirely when curating
the endpont.

exclude cancer register
------------------------
Adding "exclude cancer register" as a line to the endpoint definition file 
tells curate_endpoints.R to ignore the cancer register entirely when curating
the endpont.

follow-up from birth
---------------------
Adding "follow-up from birth" as a line to the endpoint definition file tells
curate_endpoints.R to compute the follow-up time with respect to participant 
birth rather than with respect to each UK Biobank assessment visit date. I.e.
the output will be age of first event, rather than time to event from UK Biobank
assessment.

impute follow-up as midpoint
-----------------------------
Adding "impute follow-up as midpoint" as a line to the endpoint definition file 
tells curate_endpoints.R will find the midpoint between the first disease event 
and the most recent disease-free event prior and use that as the date to compute
the follow-up time. This is useful when defining endpoints for diseases whose 
onset is typically diagnosed earlier than we would expect to find in the 
available data sources. E.g., incident type 2 diabetes from hospital records, 
where type 2 diabetes is diagnosed in primary care and thus the first occurrence
in the hospital record must have occurred at a date later than the actual 
diagnosis.

time since prevalent event
---------------------------
By default, prevalent events have 'NA' for the follow-up time as they occur
prior to the UK Biobank assessment visit. Adding "time since prevalent event" as 
a line to the endpoint definition file tells curate_endpoints.R to compute the
time since the prevalent event in years instead of reporting it as 'NA'. In this
case the time is reported as a negative number of years, to help distinguish the
prevalent events from incident events.

most recent prevalent event 
----------------------------
Adding "most recent prevalent event" as a line to the endpoint definition file 
tells curate_endpoints.R to find the most recent prevalent event when computing
the time since prevalent event, instead of reporting the time since the first
occurring prevalent event.

do not impute missing event dates 
----------------------------------
Adding "do not impute missing event dates" as a line to the endpoint definition 
file tells curate_endpoints.R not to use the dates that have been imputed for
events with no date in the underlying data. Missing event dates are particularly
prevalent in the self-reported medical history, where it often wasn't possible
for the interviewing nurse to determine an age of onset, but there are also a
handful of primary care records missing event dates. When kept as NAs, these
missing data can cascade into the follow-up time in the output, e.g. where it
becomes impossible to determine whether an event is prevalent or incident. For 
details on how missing event dates were imputed, see the 'README.txt' files in 
'common/Medical History/' and 'common/Primary Care/'.

remove missing event dates
---------------------------
Adding "remove missing event dates" as a line to the endpoint definition file 
tells curate_endpoints.R to simply exclude any records with missing event dates
when curating the endpoint.

use conservative min censor date 
---------------------------------
Adding "use conservative min censor date" as a line to the endpoint definition 
file tells curate_endpoints.R to set the minimum censor date based on the 
earliest records available in the linked data sources, rather than treating the
retrospective follow-up as going back to each participant's date of birth.

allow missing event status 
---------------------------
Adding "allow missing event status" as a line to the endpoint definition file
tells curate_endpoints.R to allow the event status (and follow-up) to be NA 
where it is not possible to determine case status. These can arise where the
date of assessment is after the maximum censor date in the linked health records
the endpoint is curated from, or where the endpoint includes self-reported 
medical history and the participant has answered "do not know" or "prefer not to 
answer" and the event status would otherwise be FALSE.

--------------------------------------------------------------------------------
The output file
--------------------------------------------------------------------------------

The output file saved by curate_endpoint.R contains the following columns:

 - 'eid':              Applicant-specific UK Biobank participant identifier.
 - 'visit_index':      Index of the UK Biobank assessment the follow-up time has
                       been calculated with respect to:
                         '0': Baseline assessment (2006-2010).
                         '1': First repeat assessment (2012-2013).
                         '2': First imaging assessment (2014+).
                         '3': Repeat imaging assessment (2019+).
                       This column is omitted if the "follow-up from birth" 
                       option has been set in the endpoint definition file (see
                       option documentation above). 
 - 'assessment_date':  Date of the assessment visit, or date of birth if the 
                       "follow-up from birth" option has been set in the 
                       endpoint definition file (see option documentation 
                       above). 
 - 'event':            TRUE where the participant is a case, FALSE otherwise. NA 
                       where "allow missing non-case status" has been set in the 
                       endpoint definition file (see option documentation above).
 - 'follow_up':        Follow-up time in years (decimal) from the assessment 
                       date of the respective 'visit_index' to the first event 
                       (cases) or maximum follow-up (controls). NA for prevalent 
                       events, i.e. where the first event happened prior to the 
                       assessment date, unless the "time since prevalent event" 
                       option has been set in the endpoint definition file (see 
                       option documentation above), in which case the follow-up 
                       time will be a negative number (years since the prevalent 
                       event). Also NA where it was not possible to calculate 
                       follow-up time, e.g. when the options "do not impute 
                       missing event dates" or "allow missing event status" have 
                       been set in the endpoint definition file (see option 
                       documentation above).
 - 'censor_date':      Date of the event, or maximum censor date, i.e. the date
                       used when computing the 'follow_up' column with respect 
                       to the 'assessment_date' column.
 - 'reason_or_id':     Unique identifier for the event identified in the 
                       underlying health record data that makes the participant 
                       a case, or for controls, the reason for the specific 
                       follow-up time (e.g. "maximum censor date").
 - 'data_source':      Path to the underlying data file that the event was 
                       identified in, or the path to the curated follow-up 
                       information used to determine maximum/minmum censor 
                       dates.
 - 'lookup_column':    Name of the column in the 'data_source' file to use to 
                       lookup the event id ('reason_or_id'), or the name of the
                       column in the follow-up data used to determine the 
                       maximum censor date for that participant.
