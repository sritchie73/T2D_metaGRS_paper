Medication data
================

This folder contains information about medications being taken by each UK Biobank
participant collected at each UK Biobank assessment visit. Information on 
medication usage was collected through a combination of simple multiple choice
questions as part of the touchscreen survey, as well as a follow-up verbal 
interview with a trained nurse who collected information on current prescription
medications being taken by the participant. Notably, these data *do not* include
medication history coded in the primary care records.

At the two imaging assessments, information was also collected during the verbal
interview about whether the participant had been taking any antibiotics in the
three months prior, and if so, what type of antibiotics.


A brief overview of the files in this folder:

-------------------------------
curated_medication_classes.csv
-------------------------------
Contains curated information on whether a participant is taking certain types of
medication at each assessment visit. Each column summarises the presence/absence
of various lists of medications, and includes information obtained from both the
verbal interview and touchscreen survey. E.g. the 'lipid_lowering_medication' 
combines information from 18 different prescription medications as well as the
touchscreen question asking each participant whether they have been prescribed
medication for their cholesterol levels. 

Full lists of medications contributing to each column can be found in 
scripts/curate_medication_classes.R

Missing values (NAs) are used where it was not possible to rule out the use of
the medication in that participant based on the combination of touchscreen and
verbal answers given; e.g. for example where a participant said "Prefer not to
answer" to the touchscreen survey and had none of the relevant medications coded
as a result of the verbal interview.

As with other datasets, rows are uniquely identifiable by the 'eid' and 
'visit_index' columns, which contain the application-specific unique identifier 
for each participant and a number indicating the UK Biobank assessment visit 
respectively (0=baseline assessment, 1=first repeat assessment, 2=first imaging 
assessment, 3=second imaging assessment).

-------------------------------------------
curated_medication_classes_column_info.csv
-------------------------------------------
Descriptions of the columns found in curated_medication_classes.csv

----------------------
medication_status.csv
----------------------
Contains TRUE/FALSE or NA values indicating whether, at each timepoint, the 
participant was recorded as taking any prescription medications either during 
the touchscreen survey or verbal interview, and at the imaging assessments, 
whether the participant reported taking any antibiotics in the 3 months prior.

----------------------------------
medication_status_column_info.csv
----------------------------------
Descriptions of the columns found in medication_status.csv

-----------------------------------
touchscreen_survey_medications.csv
-----------------------------------
Contains TRUE/FALSE or NA values for each of the five types of medications 
available as answers on the multiple-choice touchscreen survey question, as well
as TRUE/FALSE or NA values for additional touchscreen survey question asking 
whether they were taking any additional prescription medications not covered by
that first survey question. Note the "hormone_replacement_therapy" and 
"oral_contraceptives" are always set to FALSE in males as they options were not
presented on the touchscreen survey to males. NAs indicate "Do not know", or 
"Prefer not to answer" responses, or an absence of any data for that participant
at the respective assessment.

-----------------------------------------------
touchscreen_survey_medications_column_info.csv
-----------------------------------------------
Descriptions of the columns found in touchscreen_survey_medications.csv

---------------------------------
verbal_interview_medications.csv
---------------------------------
A long format table containing the full list of codes and associated labels for
the set of prescription medications coded by the trained nurse during the verbal
interview for each participant at each UK Biobank assessment.

---------------------------------------------
verbal_interview_medications_column_info.csv
---------------------------------------------
Descriptions of the columns found in verbal_interview_medications.csv

---------------------------------
verbal_interview_antibiotics.csv
---------------------------------
A long format table containing the full list of codes and associated labels for
the set of antibiotics coded by the trained nurse during the verbal interview 
for each participant at the two imaging UK Biobank assessments. Note that these
are antibiotics taken within the three months prior to the assessment visit;
antibiotics still being currently taken at the time of assessment will also 
appear in the verbal_interview_medications.csv table.

---------------------------------------------
verbal_interview_antibiotics_column_info.csv
---------------------------------------------
Descriptions of the columns found in verbal_interview_antibiotics.csv

-----------------------------------------
verbal_interview_medications_summary.csv
-----------------------------------------
A table giving the total number of prescription medications each participant was
recorded as having during the verbal interview at each assessment, and the total 
number of antibiotics each participant was recorded as taking in the three 
months prior to the imaging assessments. 

-----------------------------------------------------
verbal_interview_medications_summary_column_info.csv
-----------------------------------------------------
Descriptions of the columns found in verbal_interview_medications_summary.csv
