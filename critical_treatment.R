#--------------------------------------------------
# Preprocessing treatment data
#--------------------------------------------------
# We subset and summarize the inappropriate treatments for AmpC

rm(list=ls())   # Ensure that memory is empty


#--------------------------------------------------
# Paths
#--------------------------------------------------

path="\\\\fs-home/blipi$/Documents/"  
path_prep=paste0(path, "surveillance_preprocessed/")  # path where "preprocessed" surveillance files can be found
path_dict=paste0(path, "repo/Antibiotic_resistance_surveillance_eth/pio/dictionaries/")   # path where dictionaries for encoding antibiotics, sample material and pathogens are found


#--------------------------------------------------
# Set-up and load data
#--------------------------------------------------

# loading packages
lapply(c("tidyverse", "data.table"), library, character.only=T)

# Read files with resistance data and dictionaries
antiinfektiva=readRDS(paste0(path_prep, "antiinfektiva_data_preprocessed_2018.RDS"))
df_test=readRDS(paste0(path_prep, "ampc/test_ampc.RDS"))
dict_anti_ampc=read.csv(paste0(path_dict, "anti_ampc.csv"), sep=";")


#--------------------------------------------------
# Compute a treatment data frame, containing the encounters ID, date, ATCKEY and indicator if they were treated
#--------------------------------------------------
# We summarize the treament per day, wrong entries of treatment that was not given is corrected by a subscription of a negative amount.

# Make a data frame containing all antibiotic of interest
df_treat=antiinfektiva %>% filter(ATCKEY %in% dict_anti_ampc$ATCKEY & encounter_pseudo_id %in% df_test$encounter_pseudo_id) %>% 
  mutate(calculated_amount=as.numeric(calculated_amount)) %>% 
  summarise(calculated_amount=sum(calculated_amount, na.rm = T), .by=c("calculated_unit", "encounter_pseudo_id", "date", "ATCKEY")) %>%
  mutate(calculated_amount=pmax(calculated_amount, 0)) %>% summarise(treated=any(calculated_amount>0), .by=c("encounter_pseudo_id", "date", "ATCKEY")) %>%
  left_join(dict_anti_ampc %>% select(ATCKEY, ATCKEY_Description) %>% rename(anti_name=ATCKEY_Description), by="ATCKEY")

# Summarize the rest by others
df_treat=antiinfektiva %>% filter(!(ATCKEY %in% dict_anti_ampc$ATCKEY) & encounter_pseudo_id %in% df_test$encounter_pseudo_id) %>% 
  mutate(calculated_amount=as.numeric(calculated_amount)) %>% 
  summarise(calculated_amount=sum(calculated_amount, na.rm = T), .by=c("calculated_unit", "encounter_pseudo_id", "date", "ATCKEY")) %>%
  mutate(calculated_amount=pmax(calculated_amount, 0)) %>% summarise(treated=any(calculated_amount>0), .by=c("encounter_pseudo_id", "date", "ATCKEY")) %>%
  mutate(anti_name="others") %>%
  add_row(df_treat)


#--------------------------------------------------
# Construct a treatment dataframae with higher temporal resolution
#--------------------------------------------------
# Resolution of 1 hour.

# Make a data frame containing all antibiotic of interest
df_treat_temp=antiinfektiva %>% filter(ATCKEY %in% dict_anti_ampc$ATCKEY & encounter_pseudo_id %in% df_test$encounter_pseudo_id) %>% 
  mutate(calculated_amount=as.numeric(calculated_amount)) %>% mutate(dttm=round(start_time, unit="hours")) %>%
  summarise(calculated_amount=sum(calculated_amount, na.rm = T), .by=c("calculated_unit", "encounter_pseudo_id", "dttm", "ATCKEY")) %>%
  mutate(calculated_amount=pmax(calculated_amount, 0)) %>% summarise(treated=any(calculated_amount>0), .by=c("encounter_pseudo_id", "dttm", "ATCKEY")) %>%
  left_join(dict_anti_ampc %>% select(ATCKEY, ATCKEY_Description) %>% rename(anti_name=ATCKEY_Description), by="ATCKEY")

# Summarize the rest by others
df_treat_temp=antiinfektiva %>% filter(!(ATCKEY %in% dict_anti_ampc$ATCKEY) & encounter_pseudo_id %in% df_test$encounter_pseudo_id) %>% 
  mutate(calculated_amount=as.numeric(calculated_amount)) %>% mutate(dttm=round(start_time, unit="hours")) %>%
  summarise(calculated_amount=sum(calculated_amount, na.rm = T), .by=c("calculated_unit", "encounter_pseudo_id", "dttm", "ATCKEY")) %>%
  mutate(calculated_amount=pmax(calculated_amount, 0)) %>% summarise(treated=any(calculated_amount>0), .by=c("encounter_pseudo_id", "dttm", "ATCKEY")) %>%
  mutate(anti_name="others") %>%
  add_row(df_treat_temp)


#--------------------------------------------------
# Save the data
#--------------------------------------------------

saveRDS(df_treat_temp, file=paste0(path_prep, "ampc/treatment_ampc_temp.RDS"))
saveRDS(df_treat, file=paste0(path_prep, "ampc/treatment_ampc.RDS"))
