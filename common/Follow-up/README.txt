Follow-up
==========

Files in this folder curate the maximum and minimum follow-up dates for each
participant consenting to electronic health record linkage. 

Minimum and maximum follow-up dates are determined based on the record source
(Death registry, hospital records, cancer registry, or primary care records)
based on the earliest and latest dates UK Biobank have recieved complete data
for the given data provider for the given record source. For more information
on the different data providers for each record source and dates complete 
linkage are available see:
https://biobank.ndph.ox.ac.uk/showcase/exinfo.cgi?src=Data_providers_and_dates

Dates of coding system switch overs in the hospital records (ICD-9 to ICD-10 and 
OPCS-3 to OPCS-4) and cancer register (ICD-9 to ICD-10) are also based on 
information provided by UK Biobank at:
https://biobank.ndph.ox.ac.uk/showcase/exinfo.cgi?src=Data_providers_and_dates

Note that records may exist earlier than the minimum follow-up dates or later
than the maximum follow-up dates for some particpants and record sources, as
the records themselves are not truncated to stop or start at the dates where UK
has adjudicated the record sources to be complete across all participants. Other
cases this may occur are records appearing later than UK Biobank reports the 
participant has been lost to follow-up (UKB fields 190 and 191) or earlier than 
the year the participant reports having immigrated to the UK (UKB field 3659).

The different data providers for each record source are typically the different
nations within the UK, which manage their own distinct electronic health care 
record systems. To determine the appropriate minimum and maximum follow-up for
each record source, the nation of residence for each participant was inferred
for each participant based on the closest available record (in time) at each 
minimum or maximum censor date to determine which nation to use when determing
the respective minimum and maximum follow-up dates. While most participants 
reside in the same nation as their baseline assessment location throughout the
linkage period(s), this enables more accurate censoring of controls for any 
given endpoint if the participant has moved between nations and those nations 
have dramatically different linkage periods available.

Minimum and maximum censor dates by record source and data provider are 
detailed in 'available_followup.txt', along with minimum and maximum follow-up
dates reporting the maximum range of follow-up available when also considering
data linkage that has not yet been completely reported for all UK Biobank 
participants.

Minimum follow-up dates for each participant are also restricted to be no 
earlier than the participant's inferred date of birth, or the inferred date of 
immigration to the UK (if relevant). Each participant's date of birth was 
inferred to be the 15th of the month of birth, as exact dates of birth are not
dispensed by UK Biobank to preserve participant privacy. Year of birth and month
of birth were retrieved from UKB fields 34 and 52 respectively. Date of
immigration to the UK was inferred based on self-reported year of immigration to
the UK (UKB field 3659), with the date set at the mid point between the start of 
the year and first record available in that same year, if any. If the 
participant reported immigrating in the same year as their year of birth, the 
mid-point calculation used their inferred date of birth rather than the start of 
the year. In cases where a participant attended multiple UK Biobank assessments, 
and reported different years of immigration at each assessment, the existence of 
electronic health records in intervening or earlier years was used to adjudicate
the most likely year of immigration to the UK. Note that despite this, there are 
still some participants where linked electronic health care record exist in 
years prior to the earliest self-reported year of immigration.

Maximum follow-up dates for each participant were restricted to be no later than
their date of death if a record exists for them in the death registry (or their
death has been reported to UK Biobank from a relative; see UKB Field 191). 
Likewise, the maximum follow-up dates are also restricted to be no later than 
the date the participant was recorded as being lost to follow-up (e.g. due to 
emmigration from the UK; see UKB Fields 190 and 191).

For participants with primary care linkage (~45% of participants), the minimum
censor date was set to be the earliest primary care record existing for that 
person, as the primary care data begins in the late 1940s but is by no means 
retrospectively complete to that date. If the earliest record for the 
participant was a deregistration from a general practice, this record was 
ignored in favour of the next general practice registration or primary care 
record. The maximum censor date was determined based on the inferred nation of
residence at the end of each data provider, with the relevant data provider 
determined based on the most recent primary care record for that participant,
following the same rules as above when handling deaths and loss of follow-up.
The exception to this was if the latest available record for a participant was 
a de-registration from a clinical practice, in which case this was set as the 
maximum follow-up under the assumption that the participant had changed general
practices to one that has not provided data linkage to UK Biobank. Note there 
are a handful of participants with primary care records but without minimum or
maximum censor dates: these were all cases where the participant had only 
records with no date or data provider attached, in which case we deemed it not
possible to calculate follow-up, only presence or absence of those specific read
codes.

In addition to minimum and maximum censor dates for each record source, the file
'follow_up.csv' also reports for each record source the nation the participant 
was inferred to be living in at the given censor date, as well as the type of 
record source this evidence was determined from. The follow-up table also 
reports whether each participant has died, has any hospitalisations, has any
records in the cancer registry, or has primary care linkage (~45% of 
participants). 
