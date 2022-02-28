library(data.table)

# Controls TRUE/FALSE/NA values for prevalent/incident T2D for training and testing
case_control_definitions <- rbind(idcol="short_name",
  "prevalent_basic_adjudicated"=data.table("type"="prevalent", "cases"="adjudicated", 
     "missing"="uncertain diabetes", "controls"="others"),
  "prevalent_no_t1d"=data.table("type"="prevalent", "cases"="adjudicated", 
     "missing"="uncertain diabetes or T1D", "controls"="others"),
  "prevalent_no_undiagnosed"=data.table("type"="prevalent", "cases"="adjudicated", 
     "missing"="uncertain or possible undiagnosed diabetes", "controls"="others"),
  "prevalent_no_prediabetes"=data.table("type"="prevalent", "cases"="adjudicated", 
     "missing"="uncertain or possible prediabetes", "controls"="others"),
  "prevalent_no_undiagnosed_or_t1d"=data.table("type"="prevalent", "cases"="adjudicated", 
     "missing"="uncertain or possible undiagnosed diabetes or T1D", "controls"="others"),
  "prevalent_no_prediabetes_or_t1d"=data.table("type"="prevalent", "cases"="adjudicated", 
     "missing"="uncertain or possible prediabetes or T1D", "controls"="others"),
  "prevalent_no_nf_undiagnosed"=data.table("type"="prevalent", "cases"="adjudicated", 
     "missing"="uncertain or possible undiagnosed diabetes (non-fasting glucose)", "controls"="others"),
  "prevalent_no_nf_prediabetes"=data.table("type"="prevalent", "cases"="adjudicated", 
     "missing"="uncertain or possible prediabetes (non-fasting glucose)", "controls"="others"),
  "prevalent_no_nf_undiagnosed_or_t1d"=data.table("type"="prevalent", "cases"="adjudicated", 
     "missing"="uncertain or possible undiagnosed diabetes (non-fasting glucose) or T1D", "controls"="others"),
  "prevalent_no_nf_prediabetes_or_t1d"=data.table("type"="prevalent", "cases"="adjudicated", 
     "missing"="uncertain or possible prediabetes (non-fasting glucose) or T1D", "controls"="others"),

  "incident_basic_adjudicated"=data.table("type"="incident", "cases"="adjudicated", 
     "missing"="prevalent diabetes or uncertain incident diabetes", "controls"="others"),
  "incident_no_t1d"=data.table("type"="incident", "cases"="adjudicated", 
     "missing"="prevalent diabetes or incident uncertain diabetes or T1D", "controls"="others"),
  "incident_no_undiagnosed"=data.table("type"="incident", "cases"="adjudicated, no possible undiagnosed diabetes", 
     "missing"="prevalent diabetes or incident uncertain or possible undiagnosed diabetes", "controls"="others"),
  "incident_no_prediabetes"=data.table("type"="incident", "cases"="adjudicated, no possible prediabetes", 
     "missing"="prevalent diabetes or incident uncertain or possible prediabetes", "controls"="others"),
  "incident_no_undiagnosed_or_t1d"=data.table("type"="incident", "cases"="adjudicated, no possible undiagnosed diabetes", 
     "missing"="prevalent diabetes or incident uncertain or possible undiagnosed diabetes or T1D", "controls"="others"),
  "incident_no_prediabetes_or_t1d"=data.table("type"="incident", "cases"="adjudicated, no possible prediabetes", 
     "missing"="prevalent diabetes or incident uncertain or possible prediabetes or T1D", "controls"="others"),
  "incident_no_nf_undiagnosed"=data.table("type"="incident", "cases"="adjudicated, no possible undiagnosed diabetes", 
     "missing"="prevalent diabetes or incident uncertain or possible undiagnosed diabetes (non-fasting glucose)", "controls"="others"),
  "incident_no_nf_prediabetes"=data.table("type"="incident", "cases"="adjudicated, no possible prediabetes", 
     "missing"="prevalent diabetes or incident uncertain or possible prediabetes (non-fasting glucose)", "controls"="others"),
  "incident_no_nf_undiagnosed_or_t1d"=data.table("type"="incident", "cases"="adjudicated, no possible undiagnosed diabetes", 
     "missing"="prevalent diabetes or incident uncertain or possible undiagnosed diabetes (non-fasting glucose) or T1D", "controls"="others"),
  "incident_no_nf_prediabetes_or_t1d"=data.table("type"="incident", "cases"="adjudicated, no possible prediabetes",
     "missing"="prevalent diabetes or incident uncertain or possible prediabetes (non-fasting glucose) or T1D", "controls"="others"),

  "combined_basic_adjudicated"=data.table("type"="prevalent + incident", "cases"="adjudicated",
    "missing"="uncertain diabetes (prevalent or incident)", "controls"="others"),
  "combined_no_t1d"=data.table("type"="prevalent + incident", "cases"="adjudicated",
    "missing"="uncertain diabetes or T1D (prevalent or incident)", "controls"="others"),
  "combined_no_undiagnosed"=data.table("type"="prevalent + incident", "cases"="adjudicated",
    "missing"="uncertain or possible undiagnosed diabetes (prevalent or incident)", "controls"="others"),
  "combined_no_prediabetes"=data.table("type"="prevalent + incident", "cases"="adjudicated",
    "missing"="uncertain or possible prediabetes (prevalent or incident)", "controls"="others"),
  "combined_no_undiagnosed_or_t1d"=data.table("type"="prevalent + incident", "cases"="adjudicated",
    "missing"="uncertain or possible undiagnosed diabetes or T1D (prevalent or incident)", "controls"="others"),
  "combined_no_prediabetes_or_t1d"=data.table("type"="prevalent + incident", "cases"="adjudicated",
    "missing"="uncertain or possible prediabetes or T1D (prevalent or incident)", "controls"="others"),
  "combined_no_nf_undiagnosed"=data.table("type"="prevalent + incident", "cases"="adjudicated",
    "missing"="uncertain or possible undiagnosed diabetes (non-fasting glucose) (prevalent or incident)", "controls"="others"),
  "combined_no_nf_prediabetes"=data.table("type"="prevalent + incident", "cases"="adjudicated",
    "missing"="uncertain or possible prediabetes (non-fasting glucose) (prevalent or incident)", "controls"="others"),
  "combined_no_nf_undiagnosed_or_t1d"=data.table("type"="prevalent + incident", "cases"="adjudicated",
    "missing"="uncertain or possible undiagnosed diabetes (non-fasting glucose) or T1D (prevalent or incident)", "controls"="others"),
  "combined_no_nf_prediabetes_or_t1d"=data.table("type"="prevalent + incident", "cases"="adjudicated",
    "missing"="uncertain or possible prediabetes (non-fasting glucose) or T1D (prevalent or incident)", "controls"="others"),

  "QDiabetes2018A"=data.table("type"="incident", "cases"="adjudicated", "missing"="missing QDiabetes2018A", controls="others"),
  "QDiabetes2018B_fasting"=data.table("type"="incident", "cases"="adjudicated", "missing"="missing QDiabetes2018B_fasting", controls="others"),
  "QDiabetes2018B_non_fasting"=data.table("type"="incident", "cases"="adjudicated", "missing"="missing QDiabetes2018B_non_fasting", controls="others"),
  "QDiabetes2018C"=data.table("type"="incident", "cases"="adjudicated", "missing"="missing QDiabetes2018C", controls="others"),
  "QDiabetes2013"=data.table("type"="incident", "cases"="adjudicated", "missing"="missing QDiabetes2013", controls="others")
)



# Extract cohorts as described above from provided phenotype data
set_case_control_definition <- function(pheno, short_name) {
  dat <- copy(pheno)
  if (short_name == "prevalent_basic_adjudicated") {
    dat[, T2D := type_2_diabetes]

  } else if (short_name == "prevalent_no_t1d") {
    dat[, T2D := type_2_diabetes]
    dat[(type_1_diabetes), T2D := NA]

  } else if (short_name == "prevalent_no_undiagnosed") {
    dat[, T2D := type_2_diabetes]
    dat[(fasting_undiagnosed_diabetes), T2D := NA] 

  } else if (short_name == "prevalent_no_prediabetes") {
    dat[, T2D := type_2_diabetes]
    dat[(fasting_prediabetes), T2D := NA] 

  } else if (short_name == "prevalent_no_undiagnosed_or_t1d") {
    dat[, T2D := type_2_diabetes]
    dat[(fasting_undiagnosed_diabetes), T2D := NA] 
    dat[(type_1_diabetes), T2D := NA]

  } else if (short_name == "prevalent_no_prediabetes_or_t1d") {
    dat[, T2D := type_2_diabetes]
    dat[(fasting_prediabetes), T2D := NA] 
    dat[(type_1_diabetes), T2D := NA]

  } else if (short_name == "prevalent_no_nf_undiagnosed") {
    dat[, T2D := type_2_diabetes]
    dat[(non_fasting_undiagnosed_diabetes), T2D := NA] 

  } else if (short_name == "prevalent_no_nf_prediabetes") {
    dat[, T2D := type_2_diabetes]
    dat[(non_fasting_prediabetes), T2D := NA] 

  } else if (short_name == "prevalent_no_nf_undiagnosed_or_t1d") {
    dat[, T2D := type_2_diabetes]
    dat[(non_fasting_undiagnosed_diabetes), T2D := NA] 
    dat[(type_1_diabetes), T2D := NA]

  } else if (short_name == "prevalent_no_nf_prediabetes_or_t1d") {
    dat[, T2D := type_2_diabetes]
    dat[(non_fasting_prediabetes), T2D := NA] 
    dat[(type_1_diabetes), T2D := NA]

  } else if (short_name == "incident_basic_adjudicated") {
    dat[, T2D := incident_type_2_diabetes] # prevalent diabetes (of any kind) already NA here

  } else if (short_name == "incident_no_t1d") {
    dat[, T2D := incident_type_2_diabetes] # prevalent diabetes (of any kind) already NA here  
    dat[(incident_type_1_diabetes), T2D := NA]

  } else if (short_name == "incident_no_undiagnosed") {
    dat[, T2D := incident_type_2_diabetes] # prevalent diabetes (of any kind) already NA here  
    dat[(fasting_undiagnosed_diabetes), T2D := NA] # Note may lose some incident cases here, but consistent with limitations of QDiabetes risk scores

  } else if (short_name == "incident_no_prediabetes") {
    dat[, T2D := incident_type_2_diabetes] # prevalent diabetes (of any kind) already NA here  
    dat[(fasting_prediabetes), T2D := NA] # Note may lose some incident cases here, but consistent with limitations of QDiabetes risk scores
    
  } else if (short_name == "incident_no_undiagnosed_or_t1d") {
    dat[, T2D := incident_type_2_diabetes] # prevalent diabetes (of any kind) already NA here  
    dat[(fasting_undiagnosed_diabetes), T2D := NA] # Note may lose some incident cases here, but consistent with limitations of QDiabetes risk scores
    dat[(incident_type_1_diabetes), T2D := NA] 

  } else if (short_name == "incident_no_prediabetes_or_t1d") {
    dat[, T2D := incident_type_2_diabetes] # prevalent diabetes (of any kind) already NA here  
    dat[(fasting_prediabetes), T2D := NA] # Note may lose some incident cases here, but consistent with limitations of QDiabetes risk scores
    dat[(incident_type_1_diabetes), T2D := NA] 

  } else if (short_name == "incident_no_nf_undiagnosed") {
    dat[, T2D := incident_type_2_diabetes] # prevalent diabetes (of any kind) already NA here  
    dat[(non_fasting_undiagnosed_diabetes), T2D := NA] # Note may lose some incident cases here, but consistent with limitations of QDiabetes risk scores

  } else if (short_name == "incident_no_nf_prediabetes") {
    dat[, T2D := incident_type_2_diabetes] # prevalent diabetes (of any kind) already NA here  
    dat[(non_fasting_prediabetes), T2D := NA] # Note may lose some incident cases here, but consistent with limitations of QDiabetes risk scores
    
  } else if (short_name == "incident_no_nf_undiagnosed_or_t1d") {
    dat[, T2D := incident_type_2_diabetes] # prevalent diabetes (of any kind) already NA here  
    dat[(non_fasting_undiagnosed_diabetes), T2D := NA] # Note may lose some incident cases here, but consistent with limitations of QDiabetes risk scores
    dat[(incident_type_1_diabetes), T2D := NA] 

  } else if (short_name == "incident_no_nf_prediabetes_or_t1d") {
    dat[, T2D := incident_type_2_diabetes] # prevalent diabetes (of any kind) already NA here  
    dat[(non_fasting_prediabetes), T2D := NA] # Note may lose some incident cases here, but consistent with limitations of QDiabetes risk scores
    dat[(incident_type_1_diabetes), T2D := NA] 

  } else if (short_name == "combined_basic_adjudicated") {
    dat[, T2D := (type_2_diabetes) | (incident_type_2_diabetes)] 
    dat[(type_2_diabetes), censor_hospital_nation := assessment_nation]

  } else if (short_name == "combined_no_t1d") {
    dat[, T2D := (type_2_diabetes) | (incident_type_2_diabetes)] 
    dat[(incident_type_1_diabetes), T2D := NA] # prevalent T1D cases already set to NA by OR statement above, as it is classed as prevalent diabetes
    dat[(type_2_diabetes), censor_hospital_nation := assessment_nation]

  } else if (short_name == "combined_no_undiagnosed") {
    dat[, T2D := (type_2_diabetes) | (incident_type_2_diabetes)] 
    dat[(fasting_undiagnosed_diabetes) & !(T2D), T2D := NA] # Just want to exclude from controls, not drop incident cases, here
    dat[(type_2_diabetes), censor_hospital_nation := assessment_nation]

  } else if (short_name == "combined_no_prediabetes") {
    dat[, T2D := (type_2_diabetes) | (incident_type_2_diabetes)] 
    dat[(fasting_prediabetes) & !(T2D), T2D := NA] # Just want to exclude from controls, not drop incident cases, here
    dat[(type_2_diabetes), censor_hospital_nation := assessment_nation]

  } else if (short_name == "combined_no_undiagnosed_or_t1d") {
    dat[, T2D := (type_2_diabetes) | (incident_type_2_diabetes)]
    dat[(incident_type_1_diabetes), T2D := NA] # prevalent T1D cases already set to NA by OR statement above, as it is classed as prevalent diabetes
    dat[(fasting_undiagnosed_diabetes) & !(T2D), T2D := NA] # Just want to exclude from controls, not drop incident cases, here
    dat[(type_2_diabetes), censor_hospital_nation := assessment_nation]

  } else if (short_name == "combined_no_prediabetes_or_t1d") {
    dat[, T2D := (type_2_diabetes) | (incident_type_2_diabetes)] 
    dat[(incident_type_1_diabetes), T2D := NA] # prevalent T1D cases already set to NA by OR statement above, as it is classed as prevalent diabetes
    dat[(fasting_prediabetes) & !(T2D), T2D := NA] # Just want to exclude from controls, not drop incident cases, here
    dat[(type_2_diabetes), censor_hospital_nation := assessment_nation]

  } else if (short_name == "combined_no_nf_undiagnosed") {
    dat[, T2D := (type_2_diabetes) | (incident_type_2_diabetes)] 
    dat[(non_fasting_undiagnosed_diabetes) & !(T2D), T2D := NA] # Just want to exclude from controls, not drop incident cases, here
    dat[(type_2_diabetes), censor_hospital_nation := assessment_nation]

  } else if (short_name == "combined_no_nf_prediabetes") {
    dat[, T2D := (type_2_diabetes) | (incident_type_2_diabetes)] 
    dat[(non_fasting_prediabetes) & !(T2D), T2D := NA] # Just want to exclude from controls, not drop incident cases, here
    dat[(type_2_diabetes), censor_hospital_nation := assessment_nation]

  } else if (short_name == "combined_no_nf_undiagnosed_or_t1d") {
    dat[, T2D := (type_2_diabetes) | (incident_type_2_diabetes)]
    dat[(incident_type_1_diabetes), T2D := NA] # prevalent T1D cases already set to NA by OR statement above, as it is classed as prevalent diabetes
    dat[(non_fasting_undiagnosed_diabetes) & !(T2D), T2D := NA] # Just want to exclude from controls, not drop incident cases, here
    dat[(type_2_diabetes), censor_hospital_nation := assessment_nation]

  } else if (short_name == "combined_no_nf_prediabetes_or_t1d") {
    dat[, T2D := (type_2_diabetes) | (incident_type_2_diabetes)] 
    dat[(incident_type_1_diabetes), T2D := NA] # prevalent T1D cases already set to NA by OR statement above, as it is classed as prevalent diabetes
    dat[(non_fasting_prediabetes) & !(T2D), T2D := NA] # Just want to exclude from controls, not drop incident cases, here
    dat[(type_2_diabetes), censor_hospital_nation := assessment_nation]

  } else if (short_name == "QDiabetes2018A") {
    dat[, T2D := incident_type_2_diabetes]
    dat[is.na(QDiabetes2018A), T2D := NA]

  } else if (short_name == "QDiabetes2018B_fasting") {
    dat[, T2D := incident_type_2_diabetes]
    dat[is.na(QDiabetes2018B_fasting), T2D := NA]

  } else if (short_name == "QDiabetes2018B_non_fasting") {
    dat[, T2D := incident_type_2_diabetes]
    dat[is.na(QDiabetes2018B_non_fasting), T2D := NA]

  } else if (short_name == "QDiabetes2018C") {
    dat[, T2D := incident_type_2_diabetes]
    dat[is.na(QDiabetes2018C), T2D := NA]

  } else if (short_name == "QDiabetes2013") {
    dat[, T2D := incident_type_2_diabetes]
    dat[is.na(QDiabetes2013), T2D := NA]

  }
  return(dat[!is.na(T2D)])
}

