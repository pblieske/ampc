#--------------------------------------------------
# Collect and summaries potential confounding variables
#--------------------------------------------------
# We preprocess and aggregate all potential confounders that are observed, e.g. the patient demographics or ICU status.

rm(list=ls())   # Ensure that memory is empty


#--------------------------------------------------
# Paths
#--------------------------------------------------

path="\\\\fs-home/blipi$/Documents/"                  # Root path
path_data=paste0(path, "surveillance_data/")          # path where "raw" surveillance files can be found
path_prep=paste0(path, "surveillance_preprocessed/")  # path where "preprocessed" surveillance files can be found
path_dict=paste0(path, "repo/Antibiotic_resistance_surveillance_eth/pio/dictionaries/")   # path where dictionaries for encoding antibiotics, sample material and pathogens are found


#--------------------------------------------------
# Set-up and load data
#--------------------------------------------------

# Installing and loading packages
# lapply(c("tidyverse", "ggplot2", "data.table", "comorbidity"), install.packages, character.only=T) # Uncomment to install packages
lapply(c("ggplot2", "tidyverse", "data.table", "comorbidity"), library, character.only=T)
source(paste(path, "repo/Antibiotic_resistance_surveillance_eth/pio/chromosomal_AmpC_expression", "/utils.R", sep=""))

# Read in the RDS files
toread=c("df_encounter"="encounter_data_preprocessed", "movement_inpatient"="movemement_inpatient_data_preprocessed")
for (df in names(toread)){
  assign(df, readRDS(file=paste0(path_prep, toread[df], ".RDS")))
}

# Read in the prepossessed test file
df_test=readRDS(paste0(path_prep, "ampc/test_ampc.RDS"))
df_encounter=df_encounter %>% unique()

# Read in the csv files
toread=c("df_info"="info", "df_patient"="patient_data", "df_sofa"="encounters_sofa_score")
for (df in names(toread)){
  assign(df, read.csv(file=paste0(path_data, toread[df], ".csv"), sep=if(df=="df_sofa"){";"}else{","}))
}

# Read the diagnosis data frame and add columns name
df_diagnosis=read.csv(file=paste0(path_data, "diagnosis_non_aggregated.csv"), sep=";", header = F)

# Add the column names
colnames=c("patient_pseudo_id", "encounter_pseudo_id", "Austrittsart", "death_dttm", "weight_value", "weight_unit", "height_value", "height_unit", "encounter_type", "encounter_begin_date", "encounter_discharge_date", "code_system", "code", "code_name", "is_primary")
colnames(df_diagnosis)=colnames

# Names for the aggregated diagnosis information
# df_diagnosis_aggregatd=read.csv(file=paste0(path_data, "diagnosis_aggregated.csv"), sep=";", header = F)
# colnames=c("patient_pseudo_id", "encounter_pseudo_id", "Austrittsart", "death_dttm", "weight", "unit_weight", "height", "unit_height",  "encounter_type", "encounter_begin_date", "encounter_discharge_date", "icd_codes", "icd_names", "primary_icd_code", "primary_icd_name", "chop_codes", "chop_names", "primary_chop_code", "primary_chop_name", "first_chop_date")
# colnames(df_diagnosis_aggregatd)=colnames

# Prepare data frames containing the test dates of only the stationary patients
df_test_s = df_test %>% left_join(df_encounter, by="encounter_pseudo_id") %>% 
  filter(encounter_type=="S") %>% select(-encounter_type)


#--------------------------------------------------
# Patient information that are constant, including ages, sex, weight, and height. 
#--------------------------------------------------

# Get age and sex from the patient info data frame
get_const_confounders=function(df){
  # Help function for extracting the constant confounders from the patient data frame.
  df %>% select(encounter_pseudo_id, patient_pseudo_id, encounter_begin_date, encounter_discharge_date) %>% unique() %>% left_join(df_patient %>% select(patient_pseudo_id, sex, birth_date) %>% unique(), by="patient_pseudo_id") %>%
    mutate(birth_date=as.Date(birth_date)) %>% mutate(age = time_length(difftime(as.Date(encounter_begin_date), birth_date), "years")) %>% select(birth_date, encounter_pseudo_id, patient_pseudo_id, sex, age, encounter_begin_date, encounter_discharge_date) %>%
    mutate(sex=as.factor(sex), age=if_else(age>=0, age, 0))
}

# Obtain the constant confounder for the stationary encounters
df_con=get_const_confounders(df_test_s)

# Get height and weight from the diagnosis data frame
df_height_weight=df_diagnosis %>% select(encounter_pseudo_id, patient_pseudo_id, weight_value, height_value, weight_unit, height_unit) %>% filter(!(height_value=="NULL" & weight_value=="NULL")) %>% unique() 

# Add the height to the constant confounder data frame
df_height=df_height_weight %>% select(encounter_pseudo_id, patient_pseudo_id, height_value, height_unit) %>% filter(height_value!="NULL") %>% mutate(height_value=as.numeric(height_value))  %>% unique() %>%
  mutate(height_cm=if_else(height_unit %in% c("cm", "NULL"), height_value, ifelse(height_unit=="m", height_value/100, NA)))

# Add the heights to the confounder data frame
df_con = df_con %>% left_join(df_height %>% select(encounter_pseudo_id, height_cm), by="encounter_pseudo_id")

# Add still missing values by going over patient_id
df_con= df_con  %>% left_join(df_height %>% select(patient_pseudo_id, height_cm) %>% rename(height_update=height_cm) %>% summarise(height_update=max(height_update), .by="patient_pseudo_id"), by="patient_pseudo_id") %>%
  mutate(height_cm=if_else(is.na(height_cm), height_update, height_cm)) %>% select(-height_update)  %>%
  mutate(height_cm=if_else(height_cm>250, NA, height_cm))

# Add the weight to the constant confounder data frame
df_weight=df_height_weight %>% select(encounter_pseudo_id, patient_pseudo_id, weight_value, weight_unit) %>% filter(weight_value!="NULL") %>% mutate(weight_value=as.numeric(weight_value))  %>% unique() %>%
  mutate(weight_kg=if_else(weight_value<=500, weight_value, weight_value/1000)) %>% summarise(weight_kg=max(weight_kg), .by=c("patient_pseudo_id", "encounter_pseudo_id"))

# Add the to the df_con
df_con = df_con %>% left_join(df_weight %>% select(encounter_pseudo_id, weight_kg), by="encounter_pseudo_id")

# Add still missing values by going over patient_id
df_con = df_con %>% left_join(df_weight %>% select(patient_pseudo_id, weight_kg) %>% rename(weight_update=weight_kg) %>% summarise(weight_update=max(weight_update), .by="patient_pseudo_id"), by="patient_pseudo_id") %>%
  mutate(weight_kg=if_else(is.na(weight_kg), weight_update, weight_kg)) %>% select(-weight_update) %>% mutate(weight_kg=if_else(weight_kg>0, weight_kg, NA))

# Add the bmi
df_con=df_con %>% mutate(bmi=weight_kg/(height_cm/100)^2) 


#--------------------------------------------------
# Time varying patient information: ICU (yes/no), and movement
#--------------------------------------------------

# help function to deal with missing values
conv_date=function(date){ 
  date_conv=tryCatch(as.POSIXlt(date, format = c("%Y-%m-%d %T")), error=function(e) NA)
  indx_na=is.na(date_conv)
  date_conv[indx_na]=tryCatch(as.POSIXlt(date[indx_na], format = c("%Y-%m-%d")), error=function(e) NA)
  return(date_conv)
}

dttm_fact=c("movement_start", "movement_end")
movement_inpatient[ , dttm_fact]=apply(movement_inpatient[ , dttm_fact], 2, FUN=conv_date)

# Make a data frame containing only the movement between wards (not room and bed)
df_departement=movement_inpatient %>% mutate(remove=difftime(movement_end, movement_start)<=0) %>% filter(!remove) %>% 
  select(encounter_pseudo_id, movement_start, departement_abbreviation) %>%
  group_by(encounter_pseudo_id) %>% arrange(movement_start) %>% mutate(changed=(departement_abbreviation!=lag(departement_abbreviation))) %>%
  filter(is.na(changed) | changed) %>% select(-changed) %>% ungroup() %>% rename(entry_dttm=movement_start)

# Construct dictionary with abbreviations for the ICU's
dict_departement=movement_inpatient %>% select(departement_name, departement_abbreviation) %>% unique()
idx_ips=str_detect(tolower(dict_departement$departement_name), "ips")
dict_ips=dict_departement[idx_ips, ]

# Add the ICU information to the dataframe
df_departement=df_departement %>% mutate(ips=departement_abbreviation %in% dict_ips$departement_abbreviation)

# Add the discharge date, ward entered and left
df_departement= movement_inpatient %>% summarise(entry_dttm=max(movement_end), .by="encounter_pseudo_id") %>% mutate(departement_abbreviation="DISCHARGED", ips=F) %>%
  add_row(df_departement) %>% group_by(encounter_pseudo_id) %>% arrange(entry_dttm) %>% mutate(dep_entry=departement_abbreviation, dep_left=lag(departement_abbreviation)) %>%
  ungroup()


#--------------------------------------------------
# SOFA score
#--------------------------------------------------

# Clean SOFA score data, remove values out of range (0,24)
df_sofa=df_sofa %>% select(patient_pseudo_id, Value, Time) %>% rename(sofa=Value, dttm=Time) %>%
  filter(0<=sofa & sofa<=24)

df_sofa = df_encounter %>% select(encounter_pseudo_id, patient_pseudo_id) %>% unique() %>%
  full_join(df_sofa, by="patient_pseudo_id", multiple = "all", relationship = "many-to-many") %>%
  filter(!is.na(sofa)) %>% mutate(dttm=as.POSIXct(dttm))


#--------------------------------------------------
# AmpC hyperexpression on each ward
#--------------------------------------------------
# number of pseudomonas aerugionsa with resistance or MDR on ward within the last (two??) week
# or even simpler, count the number of patient with MDR on ward that were present at the ward
# Here we only work with the date and not the exact time as the discovery of a positive case can be 
# left or right censored and testing can take up to 3 days.

# Read the sequence data frame with the ampc information
df_seq=readRDS(paste0(path_prep, "ampc/test_ampc_seq.RDS")) 
df_seq[ , "ampc"]=df_seq$ampc$ampc

# Convert to wide fromat
df_ampc=df_seq %>% filter(ampc) %>% mutate(date=as.Date(dttm_min)) %>% 
  group_by(encounter_pseudo_id, bacteria, date) %>% summarise(ampc=sum(ampc)>0) %>% ungroup() %>%
  pivot_wider(names_from=bacteria, values_from=ampc)

# Initialize data frame to track AmpC-like on each ward, containing all combinations of dates and wards
departements=levels(as.factor(movement_inpatient$departement_abbreviation))
n=nlevels(departements)
dates=seq(as.Date("2017-07-01"), as.Date("2024-06-01"), by="days")
species=df_seq$bacteria %>% unique()
df_ampc_ward=data.frame(departement_abbreviation=rep(departements, each=length(dates)), date=dates)
df_ampc_ward[ ,species]=0

# Map the positive results to the wards that the patients have visited
df_test_temp = df_ampc %>% group_by(encounter_pseudo_id) %>% arrange(date) %>% mutate(date=as.POSIXct(date)) %>% ungroup()
df_dep_temp = df_departement %>% group_by(encounter_pseudo_id) %>% arrange(entry_dttm) %>% 
  rename(date=entry_dttm) %>% ungroup() %>% filter(!is.na(date))

# Encounters that have been tested positive at least once
encounters_test= df_test_temp %>%  filter(encounter_pseudo_id %in% df_dep_temp$encounter_pseudo_id) %>%  
  select(encounter_pseudo_id) %>% unique()

df_test_temp=df_test_temp %>% filter(encounter_pseudo_id %in% encounters_test$encounter_pseudo_id)
df_dep_temp=df_dep_temp %>% filter(encounter_pseudo_id %in% encounters_test$encounter_pseudo_id)

# Combine the test and movement data
df_ampc_ward_encounter= df_test_temp  %>% mutate(ind="test", departement_abbreviation=NA) %>%
  add_row(df_dep_temp %>% select(date, encounter_pseudo_id, departement_abbreviation) %>% mutate(ind="move"))%>% arrange(date) 

# Combine everything to one data frame containing current department and AmpC status for every date
df_ampc_ward_encounter=df_ampc_ward_encounter %>% 
  mutate(across(all_of(species), ~ifelse(!is.na(departement_abbreviation), ifelse(departement_abbreviation=="DISCHARGED", F, .x),  .x))) %>%  mutate(Date=as.Date(date)) %>%
  group_by(encounter_pseudo_id) %>% arrange(date) %>% complete(Date=seq(first(Date), last(Date), by="day")) %>%
  mutate(date=if_else(!is.na(date), date, as.POSIXct(Date))) %>% select(-Date) %>% ungroup() %>%
  add_row(df_dep_temp %>% select(date, encounter_pseudo_id, dep_left) %>% filter(!is.na(dep_left)) %>% rename(departement_abbreviation=dep_left) %>% mutate(date=date-seconds(1))) %>%
  group_by(encounter_pseudo_id) %>% arrange(date) %>% fill(c(species, "departement_abbreviation"), .direction="down") %>% ungroup() %>%
  mutate(date=as.Date(date)) %>% unique() %>% filter(!is.na(departement_abbreviation)) %>% mutate(across(all_of(species), ~if_else(is.na(.x), F, .x))) %>%
  filter(departement_abbreviation!="DISCHARGED")

rm(list=c("df_test_temp", "df_dep_temp"))

# Summarize the number of encounters with mdr for each ward and date
df_ampc_ward = df_ampc_ward %>% rows_update(df_ampc_ward_encounter %>% summarise(across(all_of(species), ~sum(.x)), .by=c("departement_abbreviation", "date")), by=c("departement_abbreviation", "date"))


#--------------------------------------------------
# Exposure to AmpC expressing species
#--------------------------------------------------

# Extract all the stationary encounters
encounters_s=df_test_s$encounter_pseudo_id %>% unique()

# Subsetting df for stationary encounters
df_departement_s=df_departement %>% filter(encounter_pseudo_id %in% encounters_s) %>%
  group_by(encounter_pseudo_id) %>% arrange(entry_dttm) %>% mutate(discharge_dttm=lead(entry_dttm)) %>%
  ungroup() %>% mutate(discharge_date=as.Date(discharge_dttm), entry_date=as.Date(entry_dttm)) %>%
  select(departement_abbreviation, encounter_pseudo_id, entry_date, discharge_date)

# Allocate the memory for the exposure data
df_exposure=setDT(df_test_s %>% select(encounter_pseudo_id, encounter_begin_date, encounter_discharge_date) %>% unique())[ , .(date=as.Date(seq(encounter_begin_date, encounter_discharge_date+days(1), by="day"))), by=c("encounter_pseudo_id")] %>%
  mutate(date=as.Date(date))
df_exposure[, species]=0

# Convert the logical to numerical to be able to sum them up later
df_ampc_ward_encounter=df_ampc_ward_encounter %>% mutate(across(all_of(species), ~as.numeric(.x)))

# Go through every encounter and determine the number of people on the same ward with AmpC hyperexpression
# Unfortunately, this can take up to 2 hours.

# Prepare the joining argument and progress bar
by=join_by(departement_abbreviation==departement_abbreviation, date>=entry_date, date<=discharge_date)
pb = txtProgressBar(min = 0, max = length(encounters_s), initial = 0, style=3) 

# Track all potential contacts back
for(i in 1:length(encounters_s)){
  encounter=encounters_s[i]
  df_dep_temp=df_departement_s %>% filter(encounter_pseudo_id==encounter) %>% mutate(matched=T)
  df_add=df_ampc_ward_encounter %>% filter(encounter_pseudo_id!=encounter) %>% left_join(df_dep_temp %>% mutate(matched=T), by=by) %>% filter(matched) %>% summarise(across(all_of(species), ~sum(na.omit(.x))), .by="date")
  df_exposure=df_exposure %>% rows_update(df_add %>% mutate(encounter_pseudo_id=encounter), by=c("encounter_pseudo_id", "date"))
  setTxtProgressBar(pb,i)   # Update progress bar
}

close(pb)


#--------------------------------------------------
# Summaries the diagnosis
#--------------------------------------------------

# Add the in-hospital mortality and date of death to df_con
df_con=df_diagnosis %>% filter(encounter_pseudo_id %in% encounters_s) %>% select(encounter_pseudo_id, death_dttm, Austrittsart)  %>% unique() %>%
  filter(Austrittsart=="05") %>% rowwise() %>% mutate(in_hospital_mortality=T, death_dttm=tryCatch(as.POSIXct(death_dttm),  error=function(e) NA)) %>% ungroup() %>% select(-Austrittsart) %>% 
  right_join(df_con, by="encounter_pseudo_id") %>% mutate(in_hospital_mortality=if_else(is.na(in_hospital_mortality), F, in_hospital_mortality))

# Make a summary data frame with the diagnosis
df_diagnoses=df_diagnosis %>% filter(encounter_pseudo_id %in% encounters_s) %>% select(encounter_pseudo_id, patient_pseudo_id, code_system, code, code_name, is_primary) %>%
  filter(code_system=="ICD") %>% select(-code_system) %>% unique() %>% mutate(is_primary=is_primary==1, code=as.factor(code))

# Compute the Elixhauser and Charlson Comorbidity Index
df_elixhauser=comorbidity(x = df_diagnoses, id = "encounter_pseudo_id", code = "code", map = "elixhauser_icd10_quan", assign0 = T)
df_charlson=comorbidity(x = df_diagnoses, id = "encounter_pseudo_id", code = "code", map = "charlson_icd10_quan", assign0 = T)

df_charlson$score=score(df_charlson, weights="charlson", assign0 = T)
df_elixhauser$score=score(df_elixhauser, assign0 = T)

# Add the summary score to the confounder data frame
df_con=df_con %>% left_join(df_elixhauser %>% select(encounter_pseudo_id, score) %>% rename(elixhauser=score), by="encounter_pseudo_id") %>%
  left_join(df_charlson %>% select(encounter_pseudo_id, score) %>% rename(charlson=score), by="encounter_pseudo_id")

# Add the discharge codes
df_con=df_diagnosis %>% filter(encounter_pseudo_id %in% df_con$encounter_pseudo_id) %>% select(encounter_pseudo_id, Austrittsart) %>% unique() %>%
  rename(discharge_code=Austrittsart) %>% right_join(df_con, by="encounter_pseudo_id")

# After a sanity check, discharge encoded by "8" could be both, death or discharge.
# We corrected them by checking them by hand in KISIM.
df_correct=read.csv(paste0(path_prep, "ampc/df_check_discharge_dw.csv"), sep=";", header = T)
df_con=df_con %>% left_join(df_correct %>% select(encounter_pseudo_id, in_hospital_mortality, death_dttm) %>% rename(m_c=in_hospital_mortality, dttm_c=death_dttm), by="encounter_pseudo_id") %>%
  mutate(discharge_code=if_else(!is.na(m_c), if_else(m_c, "05", "01"), discharge_code), death_dttm=if_else(!is.na(dttm_c), as.Date(dttm_c), death_dttm)) %>%
  select(-c("m_c", "dttm_c")) %>% 
  mutate(in_hospital_mortality=if_else(!is.na(discharge_code), discharge_code=="05", F)) %>% 
  mutate(death_dttm=if_else(in_hospital_mortality & is.na(death_dttm), encounter_discharge_date, death_dttm))


#--------------------------------------------------
# Save the data frames
#--------------------------------------------------

# Save everything to the path_prep/ampc
tosave=c("df_con", "df_departement", "df_exposure", "df_diagnoses", "df_charlson", "df_elixhauser", "df_sofa")
invisible(lapply(tosave, FUN = function(x){saveRDS(get(x), file=paste0(path_prep,"ampc/", x, ".RDS"))}))
