Death records
==============

This folder contains death registry records for UK Biobank participants. For 
more information see: https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100093

Briefly, this folder contains the following files:

-----------
deaths.csv
-----------
A table containing the date of death and issuing nation/record source for any
UK Biobank participants that have died in England/Wales/Scotland. 

Note that a small number of participants have more than one death record due to
multiple differing death certificates being issued by different data providers,
e.g. following a postmortem. This is indicated by the 'certificate_number'
column.

------------------------------
deaths_column_information.csv
------------------------------
Descriptions of the columns found in deaths.csv

-----------------
death_causes.csv
-----------------
Table containing information about the cause(s) of death for each death record.
The death_causes.csv table can be linked to the deaths.csv table by the 
'dnx_death_id' and 'eid' columns. The cause of death is given as an ICD-10 code,
and information on whether the cause was listed as the primary cause on the 
death certificate (always only one ICD-10 code) or a secondary cause of death
(may be multiple ICD-10 coes).

------------------------------------
death_causes_column_information.csv
------------------------------------
Descriptions of the columns found in death_causes.csv

----------------
icd10_codes.csv
----------------
Descriptions for each ICD-10 code, along with a listing of the full set of 
possible ICD-10 codes used by UK Biobank. Note the majority of these do not
actually appear in any death records; an indicator column is provided to
enable to filtering to only causes of death that are actually present in 
death_causes.csv

-----------------------------------
icd10_codes_column_information.csv
-----------------------------------
Descriptions of the columns found in icd10_codes.csv
