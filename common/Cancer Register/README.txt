Hospital Records
=================

This folder contains data from cancer registers for UK Biobank participants. For
more information see: https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100092

Briefly, this folder contains the following files:

--------------------
cancer_register.csv
--------------------
Master table of all cancer registry records and diagnosis codes from different
record sources. Each row corresponds to a single cancer record - UK Biobank has
made every effort to deduplicate multiple instances of the same cancer recorded
through different sources but cannot exclude the possibility that some of these
may be (pseudo-) duplicates from different providers. For more information see:
https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/CancerLinkage.pdf
 
ICD-O-3 represent five digit codes from the the International Classification of 
Diseases for Oncology, 3rd Revision (ICD-O-3), ranging from M-8000/0 to M-9989/3. 
The first four digits (after the M) code the histology and the fifth digit codes 
the behaviour. Note that ICD-O morphology codes have undergone multiple 
revisions over time, and it is often impossible to determine which version of 
the system has been used. For more details see:
https://iris.who.int/bitstream/handle/10665/96612/9789241548496_eng.pdf
 
-----------------------------------------
cancer_register_column_information.csv
-----------------------------------------
Descriptions of the columns found in cancer_register.csv

----------------
icd10_codes.csv
----------------
Descriptions for each ICD-10 code, along with a listing of the full set of 
possible ICD-10 codes used by UK Biobank. Note a large fraction of these do not
actually appear in any cancer records; an indicator column is provided to
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
actually appear in any cancer records; an indicator column is provided to
enable to filtering to only diagnoses that are actually present in diagnoses.csv

-----------------------------------
icd9_codes_column_information.csv
-----------------------------------
Descriptions of the columns found in icd9_codes.csv
