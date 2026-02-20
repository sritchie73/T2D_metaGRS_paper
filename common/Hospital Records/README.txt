Hospital Records
=================

This folder contains hospital inpatient records and diagnosis or operation 
codes for each event. For further details, see: 
https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=2000

Briefly, this folder contains the following files:

----------------------
hospital_episodes.csv
----------------------
Master table of all hospital episodes from different record sources. Each row
corresponds to a single hospital episode for a UK Biobank participant in any
hospitals in England, Wales, or Scotland. A hospital episode is defined as a 
continuous period of admitted patient care administered under one consultant 
within a hospital. Records from hospitals in England and Wales are also grouped
into spells covering a a total continuous stay of a patient in a single hospital 
from admission to discharge. For further details see:
https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/HospitalEpisodeStatistics.pdf

-----------------------------------------
hospital_episodes_column_information.csv
-----------------------------------------
Descriptions of the columns found in hospital_episodes.csv

--------------
diagnoses.csv
--------------
A table of diagnosis codes associated with each hospital record in 
hospital_episodes.csv

---------------------------------
diagnoses_column_information.csv
---------------------------------
Descriptions of the columns found in diagnoses.csv

---------------
operations.csv
---------------
For any hospital stays in hospital_episodes.csv which included an operative 
procedure, provides a table of procedure codes and pre/post operation 
information for each procedure.

----------------------------------
operations_column_information.csv
----------------------------------
Descriptions of the columns found in operations.csv

----------------
icd10_codes.csv
----------------
Descriptions for each ICD-10 code, along with a listing of the full set of 
possible ICD-10 codes used by UK Biobank. Note a large fraction of these do not
actually appear in any hospital records; an indicator column is provided to
enable to filtering to only diagnoses that are actually present in diagnoses.csv

-----------------------------------
icd10_codes_column_information.csv
-----------------------------------
Descriptions of the columns found in icd10_codes.csv

----------------
icd9_codes.csv
----------------
Descriptions for each ICD-9 code, along with a listing of the full set of 
possible ICD-9 codes used by UK Biobank. Note the majority of these do not
actually appear in any hospital records; an indicator column is provided to
enable to filtering to only diagnoses that are actually present in diagnoses.csv

-----------------------------------
icd9_codes_column_information.csv
-----------------------------------
Descriptions of the columns found in icd9_codes.csv

----------------
opcs4_codes.csv
----------------
Descriptions for each OPCS-4 code, along with a listing of the full set of 
possible OPCS-4 codes used by UK Biobank. Note a large fraction of these do not
actually appear in any hospital records; an indicator column is provided to
enable to filtering to only procedures that are actually present in 
operations.csv

-----------------------------------
opcs4_codes_column_information.csv
-----------------------------------
Descriptions of the columns found in opcs4_codes.csv

----------------
opcs3_codes.csv
----------------
Descriptions for each OPCS-3 code, along with a listing of the full set of 
possible OPCS-3 codes used by UK Biobank. Note a large fraction of these do not
actually appear in any hospital records; an indicator column is provided to
enable to filtering to only procedures that are actually present in 
operations.csv

-----------------------------------
opcs3_codes_column_information.csv
-----------------------------------
Descriptions of the columns found in opcs3_codes.csv
