Primary Care Records
=====================

This folder contains data recorded by health professionals working at general 
practices. Note that only ~45% of the UK Biobank participants currently have 
data linkage to primary care records. For more information see: 
https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=3000
https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=591
https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=592

Briefly, this folder contains the following files:

---------------------
gp_registrations.csv
---------------------
A table containing dates each UK Biobank participant with primary care linkage
was registered at a specific general practice, and if applicable, de-registered
from that same general practice. Note that information about the specific  
general practices is not provided; participants may have multiple entries that
are distinguishable only by combinations of registration and de-registration
dates; this table is mostly useful to check periods each participant has record
coverage in the other tables.

Special dates used by UK Biobank to code missing or erroneous information (Data
Coding 819) have been replaced with NAs:

Coding	    Meaning
1900-01-01  Code has no event date
1901-01-01  Code has event date before participant's date of birth
1902-02-02  Code has event date matching participant's date of birth
1903-03-03  Code has event date after participant's date of birth and falls in 
            the same calendar year as date of birth
1909-09-09  Code has event date in the future and is presumed to be a place-
            holder or other system default
2037-07-07  Code has event date in the future and is presumed to be a place-
            holder or other system default

----------------------------------------
gp_registrations_column_information.csv
----------------------------------------
Descriptions of the columns found in gp_registrations.csv

------------------------
gp_clinical_records.csv
------------------------
Table of clinical records for each participant. Clinical records are coded using
either Read code version 2 or Read code version 3 depending on the data source.
See https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=591 for further details.

The meaning of each Read code can be obtained from the 'read_v2_lkp' and 
'read_ctv3_lkp' tables for Read codes version 2 and 3 respectively in 
read_code_lookup_maps.xlsx

As with gp_records.csv, special dates used by UK Biobank to code missing or 
erroneous information (Data Coding 819) have been replaced with NAs.

For convenience and to reduce compute in downstream pipelines (e.g. 
curate_endpoint.R) records missing event dates have had those event dates 
imputed as the midpoint in the records available for that person. These can
be filtered out or returned to NAs by filtering on the 'imputed_date' column.

-------------------------------------------
gp_clinical_records_column_information.csv
-------------------------------------------
Descriptions of the columns found in gp_clinical_records.csv

---------------------
gp_prescriptions.csv
---------------------
Table of issued prescriptions for each participant. Drugs may be coded with 
Read codes (version 2), British National Formulary codes, and/or Dictionary of 
Medicines and Devices codes. Drug names and quantities are also provided.
For further details see https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=591

The meaning of each code can be obtained from the 'read_v2_drugs_lkp', 
'bnf_lkp', and 'dmd_lkp' sheets respectively in read_code_lookup_maps.xlsx

As with gp_records.csv, special dates used by UK Biobank to code missing or 
erroneous information (Data Coding 819) have been replaced with NAs.

----------------------------------------
gp_prescriptions_column_information.csv
----------------------------------------
Descriptions of the columns found in gp_clinical_records.csv

----------------
read2_codes.csv
----------------
Table containing information on Read-2 codes from read_code_lookup_maps.xlsx, 
including their description, ICD-10 and OPCS-4 mappings, and a TRUE/FALSE column
indicating whether that Read code occurs in gp_clinical_records.csv

-----------------------------------
read2_codes_column_information.csv
-----------------------------------
Descriptions of the columns found in read2_codes.csv

----------------
read3_codes.csv
----------------
Table containing information on Read-3 codes from read_code_lookup_maps.xlsx, 
including their description, ICD-10 and OPCS-4 mappings, and a TRUE/FALSE column
indicating whether that Read code occurs in gp_clinical_records.csv

-----------------------------------
read3_codes_column_information.csv
-----------------------------------
Descriptions of the columns found in read3_codes.csv

---------------------------
read_code_lookup_maps.xlsx
---------------------------
Excel file provided by UK Biobank containing various lookup maps for read codes 
and how they map to other code types like ICD-10. 

----------------------------------------
read_code_lookup_maps_documentation.pdf
----------------------------------------
Documentation file provided by UK Biobank for read_code_lookup_maps.xlsx
