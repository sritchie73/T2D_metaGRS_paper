library(data.table)

# Download curated medication data if not already present on local instance
system("mkdir -p Medication")
system("dx download -f common/Medications/verbal_interview_medications.csv -o Medication/")
system("dx download -f common/Medications/touchscreen_survey_medications.csv -o Medication/")
system("dx download -f common/Medications/medication_status.csv -o Medication/")

verbal <- fread("Medication/verbal_interview_medications.csv")
tchscrn <- fread("Medication/touchscreen_survey_medications.csv")
med_status <- fread("Medication/medication_status.csv")

# Combine touchscreen medications with verbal interview medications into one long format table
tchscrn[, other_prescription_medication := NULL]
tchscrn <- melt(tchscrn, id.vars=c("eid", "visit_index"), variable.name="medication_name")
tchscrn <- tchscrn[(value)]
tchscrn[, medication_code := fcase(
  medication_name == "cholesterol_medication", 1L,
  medication_name == "blood_pressure_medication", 2L, 
  medication_name == "insulin_medication", 3L,
  medication_name == "hormone_replacement_therapy", 4L,
  medication_name == "oral_contraceptive", 5L
)]
tchscrn[, value := NULL]
meds <- rbind(tchscrn, verbal)

# Add in additional entries based on medication status so that all participants
# end up in the resulting table derive from medication combinations
med_status <- med_status[!meds, on = .(eid, visit_index)]
med_status[!(any_current_medications), c("medication_code", "medication_name") := .(0, "No prescription medications")]
med_status[(any_current_medications), c("medication_code", "medication_name") := .(NA, "Other medication not recorded by nurse")] # Becomes NA, rather than FALSE downstream 
med_status[is.na(any_current_medications), c("medication_code", "medication_name") := .(NA, "Uncertain medication usage or status")]
med_status <- med_status[,.(eid, visit_index, medication_code, medication_name)]
meds <- rbind(meds, med_status)

# Curate custom classes of medications
med_classes <- meds[, by=.(eid, visit_index), .(
  # Curate data on lipid lowering medications
  #
  # List of drugs from BNF chapter 2.12: Lipid-regulating drugs
  # https://openprescribing.net/bnf/0212/
  #
  # Note not all drugs listed in BNF are found in the list of UKB
  # medication codes. All drugs have also been cross-checked against
  # drugbank to ensure they are capturing the right condition (i.e.
  # lipid lowering to treat hypercholesterolemia and prevent CVD)
  #
  # E.g., while Ispaghula husk is listed, is primary indication is for constipation,
  # not lipid lowering, so its not included.
  # https://bnf.nice.org.uk/drug/ispaghula-husk.html
  #
  # Conversely, colestyramine is listed, even though it has multiple indications,
  # because one of its main indications is for primary prevention of coronary heart disease
  # https://bnf.nice.org.uk/drug/colestyramine.html
  #
  lipid_lowering_medication = ifelse(
    any(medication_code == 1140861892) |  # acipimox
    any(medication_code == 1141146234) |  # atorvastatin
    any(medication_code == 1140861924) |  # bezafibrate
    any(medication_code == 1141157260) |  # bezafibrate product
    any(medication_code == 1140862026) |  # ciprofibrate
    any(medication_code == 1140888590) |  # colestipol
    any(medication_code == 1140909780) |  # colestyramine
    any(medication_code == 1141180734) |  # colestyramine product
    any(medication_code == 1141180722) |  # colestyramine+aspartame 4g/sachet powder
    any(medication_code == 1141192736) |  # ezetimibe
    any(medication_code == 1140861954) |  # fenofibrate
    any(medication_code == 1140888594) |  # fluvastatin
    any(medication_code == 1140861856) |  # gemfibrozil
    any(medication_code == 1141157262) |  # gemfibrozil product
    any(medication_code == 1140861868) |  # nicotinic acid product
    any(medication_code == 1140888648) |  # pravastatin
    any(medication_code == 1141192410) |  # rosuvastatin
    any(medication_code == 1140861958) |  # simvastatin
    any(medication_code == 1),            # Cholesterol lowering medication reported on touchscreen survey
    TRUE, FALSE),
  
  # Curate data on treated hypertension
  #
  # List of drugs from BNF chapter 2.5: Hypertension and heart failure
  # https://openprescribing.net/bnf/0205/
  #
  # Cross-checked with drug indication at https://bnf.nice.org.uk/drug to ensure
  # those specifically for hypertension (e.g. instead of heart failure) are curated
  #
  # Note not all drugs listed in BNF are found in the list of UKB medication codes.
  #
  hypertension_medication = ifelse(
    any(medication_code == 1141186674) |  # bosentan
    any(medication_code == 1140888686) |  # hydralazine
    any(medication_code == 1140860532) |  # minoxidil
    any(medication_code == 1141168936) |  # sildenafil
    any(medication_code == 1141187810) |  # tadalafil
    any(medication_code == 1140883468) |  # clonidine
    any(medication_code == 1140871986) |  # clonidine hydrochloride 25micrograms tablet
    any(medication_code == 1140860470) |  # methyldopa
    any(medication_code == 1140860562) |  # methyldopa+hydrochlorothiazide 250mg/15mg tablet
    any(medication_code == 1140910606) |  # alpha methyldopa
    any(medication_code == 1140928284) |  # moxonidine
    any(medication_code == 1140888536) |  # guanethidine
    any(medication_code == 1140879778) |  # doxazosin
    any(medication_code == 1140879782) |  # indoramin
    any(medication_code == 1141157490) |  # indoramin product
    any(medication_code == 1140879794) |  # prazosin
    any(medication_code == 1140879798) |  # terazosin
    any(medication_code == 1141156836) |  # candesartan cilexetil
    any(medication_code == 1140860750) |  # captopril
    any(medication_code == 1140860764) |  # captopril+hydrochlorothiazide 25mg/12.5mg tablet
    any(medication_code == 1140860882) |  # cilazapril
    any(medication_code == 1141181186) |  # co-zidocapt 25mg/12.5mg tablet
    any(medication_code == 1140888552) |  # enalapril
    any(medication_code == 1140860790) |  # enalapril maleate+hydrochlorothiazide 20mg/12.5mg tablet
    any(medication_code == 1141171336) |  # eprosartan
    any(medication_code == 1140888556) |  # fosinopril
    any(medication_code == 1141164148) |  # imidapril hydrochloride
    any(medication_code == 1141152998) |  # irbesartan
    any(medication_code == 1141172682) |  # irbesartan+hydrochlorothiazide 150mg/12.5mg tablet
    any(medication_code == 1140860696) |  # lisinopril
    any(medication_code == 1140864952) |  # lisinopril+hydrochlorothiazide 10mg/12.5mg tablet
    any(medication_code == 1140916356) |  # losartan
    any(medication_code == 1141151016) |  # losartan potassium+hydrochlorothiazide 50mg/12.5mg tablet
    any(medication_code == 1140923712) |  # moexipril
    any(medication_code == 1141193282) |  # olmesartan
    any(medication_code == 1140879802) |  # amlodipine
    any(medication_code == 1140888560) |  # perindopril
    any(medication_code == 1141180592) |  # perindopril+indapamide
    any(medication_code == 1140860728) |  # quinapril
    any(medication_code == 1140860806) |  # ramipril
    any(medication_code == 1141165470) |  # felodipine+ramipril
    any(medication_code == 1141166006) |  # telmisartan
    any(medication_code == 1141187788) |  # telmisartan+hydrochlorothiazide 40mg/12.5mg tablet
    any(medication_code == 1140860904) |  # trandolapril
    any(medication_code == 1141153328) |  # trandolapril+verapamil hydrochloride
    any(medication_code == 1140888510) |  # verapamil
    any(medication_code == 1141145660) |  # valsartan
    any(medication_code == 1141201038) |  # valsartan+hydrochlorothiazide 80mg/12.5mg tablet
    any(medication_code == 2),            # Blood pressure medication reported on touchscreen survey
    TRUE, FALSE),

  # Curate data on insulin treatment
  # 
  # From BNF chapter 6.1.1
  # https://openprescribing.net/bnf/060101/
  insulin_medication = ifelse(
    any(medication_code == 1140883066) | # insulin product
    any(medication_code == 3),           # Insulin usage reported on touchscreen survey
    TRUE, FALSE),
  
  # Second generation atypical antipsychotics (used by QDiabetes and QRISK algorithms)
  atypical_antipsychotics = ifelse(
    any(medication_code == 1141153490) |  # amisulpride
    any(medication_code == 1141195974) |  # aripiprazole
    any(medication_code == 1140867420) |  # clozapine
    any(medication_code == 1140928916) |  # olanzapine
    any(medication_code == 1141152848) |  # quetiapine
    any(medication_code == 1140867444) |  # risperidone
    any(medication_code == 1140927956) |  # sertindole
    any(medication_code == 1141169714),   # zotepine
    TRUE, FALSE),
  
  # Routine use of corticosteroids (used by QDiabetes and QRISK algorithms)
  systematic_corticosteroids = ifelse(
    any(medication_code == 1140874790) |  # betamethasone
    any(medication_code == 1141145782) |  # deflazacort
    any(medication_code == 1140874816) |  # dexamethasone
    any(medication_code == 1140874896) |  # hydrocortisone
    any(medication_code == 1140874976) |  # methylprednisolone
    any(medication_code == 1140874930) |  # prednisolone
    any(medication_code == 1141157402) |  # prednisolone product
    any(medication_code == 1140868364) |  # prednisone
    any(medication_code == 1140868426),   # triamcinolone
    TRUE, FALSE),
  
  # Erectile dysfunction treatments, used by QRISK
  erectile_dysfunction = ifelse(
    any(medication_code == 1140869100) | # alprostadil
    any(medication_code == 1140883010) | # papaverine
    any(medication_code == 1141168936) | # sildenafil
    any(medication_code == 1141187810) | # tadalafil
    any(medication_code == 1141192248) | # vardenafil 
    any(medication_code == 1140865136),  # yohimbine/pemoline/methyltestosterone
    TRUE, FALSE)
)]

fwrite(med_classes, "Medication/curated_medication_classes.csv")

info <- rbind(use.names=TRUE, fill=TRUE, 
  data.table(var="eid", name="Application-specific Participant ID"),
  data.table(var="visit_index", name="UK Biobank assessment visit: 0 = Baseline Assessment, 1 = First repeat assessment, 2 = First imaging assessment, 3 = Second imaging assessment"),
  data.table(var="lipid_lowering_medication", name="Self-reporting cholesterol medication on touchscreen survey, or prescribed any lipid regulating drugs listed on the British National Formulary (see scripts/curated_medication_classes.R for full list)"),
  data.table(var="hypertension_medication", name="Self-reporting blood pressure medication on touchscreen survey, or prescribed any hypertension drugs listed on the British National Formulary (see scripts/curated_medication_classes.R for full list)"),
  data.table(var="insulin_medication", name="Self-reporting insulin usage on touchscreen survey, or had an insulin prescription (verbal interview with nurse; UKB Field ID: 20003)"),
  data.table(var="atypical_antipsychotics", name="Prescribed any drugs on the list of atypical antipsychotics used by the QRISK or QDIABETES risk models (see scripts/curated_medication_classes.R for full list)"),
  data.table(var="systematic_corticosteroids", name="Prescribed any drugs on the list of systematic corticosteroids used by the QRISK or QDIABETES risk models (see scripts/curated_medication_classes.R for full list)"),
  data.table(var="erectile_dysfunction", name="Prescribed any drugs on the list of erectile dysfunction treatments used by the QRISK risk models (see scripts/curated_medication_classes.R for full list)")
)
fwrite(info, "Medication/curated_medication_classes_column_info.csv")

# Upload to persistent storage
system("dx upload Medication/curated_medication_classes.csv Medication/curated_medication_classes_column_info.csv --destination common/Medications/")
