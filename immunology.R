#--------------------------------------------------
# Preprocessing immunology data
#--------------------------------------------------
# We preprocess 

rm(list=ls())   # Ensure that memory is empty


#--------------------------------------------------
# Load the data
#--------------------------------------------------

# Paths
path="\\\\fs-home/blipi$/Documents/"    
path_data=paste0(path, "surveillance_data/")          # path where "raw" surveillance files can be found
path_prep=paste0(path, "surveillance_preprocessed/")  # path where "preprocessed" surveillance files can be found
path_dict=paste0(path, "repo/Antibiotic_resistance_surveillance_eth/pio/dictionaries/")   # path where dictionaries for encoding antibiotics, sample material and pathogens are found

# Load the packages
invisible(lapply(c("tidyverse", "lubridate", "dplyr", "gtsummary"), library, character.only=T))
source(paste(path, "repo/Antibiotic_resistance_surveillance_eth/pio/chromosomal_AmpC_expression/utils.R", sep=""))
select=dplyr::select   # Make sure that select from the dplyr package is not overwritten by the MASS package

# Read the data files and dictionaries
dict_imm=read.csv(paste0(path_dict, "dict_imm_final.csv"), sep=";")
df_imm=read.csv(paste0(path_data, "immunology_20260311/20260311_lab_values_chosen_ids.csv"), sep=";", header=T) %>% 
  rename(obs_name=ObservationIdentifierId, obs_val=ObservationValue, obs_units=Units)


#--------------------------------------------------
# Format conversion and remove missing
#--------------------------------------------------

# Mapping of the name to the codes
df_imm=df_imm %>% left_join(dict_imm, by="obs_name")

# Low PCT values are not measurable, replace them with 0
df_imm = df_imm %>% mutate(obs_val=if_else(cat_code=="pct" & obs_val=="Nicht messbar", "0", obs_val))

# Remove results that have been canceled
results_rm=c("Massenhaft Erythrozyten, Urinstatus nicht durchführbar", "Präparat archiviert aber nicht differenziert.", "#nrm", "NULL", "siehe Sediment, mikroskopisch-visuell", "s.Text", "Durch Regelwerk nachverordnet. Das Präparat wird archiviert.", "Für die Beurteilung des Sediments zu wenig Material.", "keine Angaben", "Material geronnen. Neue EDTA-Blutentnahme angefordert.", "Material geronnen. Neue Probe angefordert.", "falsches Material siehe UZL \"Analysen-Auskunftsystem\", www.uzl.usz.ch", "leer", "misslungen", "nicht beurteilbar", "Fehlerhafte Blutentnahme.", "Probe verdünnt.", "Zu wenig Blut im Röhrchen (ergibt Verdünnungsfehler).", "Material geronnen. Neue EDTA-Bluentnahme angefordert.", "falsches Mat.", "Kein Mat.", "keine Berechnung möglich", "Analyse annulliert.", "Analyse annulliert", "gelöscht", "Labortechnische Panne. Auftrag nicht ausgeführt.", "Kein EDTA-Röhrchen erhalten.", "Nicht messbar", "alsches Material siehe UZL \"Analysen-Auskunftsystem\", www.uzl.usz.ch", "k.Mat.", "Verwechslung (falscher Patient gestochen)", "geronnen", "Follows", "STARKE Störung d. Hämolyse, neue Probe abnehmen!", "Analyse nicht durchführbar, zu wenig Material erhalten.", "annulliert", "Nicht durchgeführt - Verordnung wurde durch Auftraggeber gelöscht")
df_imm=df_imm %>% filter(!(obs_val %in% results_rm))

# When an interval is given, compute the mean
df_imm=df_imm %>% mutate(obs_val=str_replace_all(obs_val, pattern=c(">"="", "<"=""))) %>%
  mutate(loc=regexpr(pattern="-", obs_val)) %>%
  mutate(obs_val=ifelse(loc!=-1, as.character(0.5*(as.numeric(substr(obs_val, 0, loc-1))+as.numeric(substr(obs_val, loc+1, nchar(obs_val))))), obs_val)) %>%
  select(-loc)

# Correct the categorization of the leukocytes
df_imm = df_imm %>% mutate(cat_code=if_else(cat_code=="leukozyten" & is.na(as.numeric(obs_val)), "leukozyten_q", cat_code))


#--------------------------------------------------
# Format conversion and remove missing
#--------------------------------------------------

convert_units=function(obs_val, obs_units, cat_code){
  if(cat_code %in% c("lymphozyten", "leukozyten") & obs_units %in% c("Zellen/µl", "/µl")){
      obs_val=obs_val/10^3
      obs_units="G/l"
  }else if(cat_code %in% c("il-6") & obs_units=="pg/ml"){
        obs_val=obs_val
        obs_units="ng/l"
  }
  return(tibble(obs_val=obs_val, obs_units=obs_units))
}

# Convert the units to the target unit
df_sum=df_imm %>% filter(!is.na(cat_code)) %>%
  mutate(obs_val=as.numeric(obs_val)) %>%
  rowwise() %>% mutate(convert_units(obs_val = obs_val, obs_units=obs_units, cat_code=cat_code)) %>%
  ungroup()

# Tidy and clean the dataframe
df_sum=df_sum %>% filter(target_unit==obs_units) %>% 
  select(c(encounter_pseudo_id, SampleDrawDtTm, ObservationEndDtTm, ObservationDtTm, obs_val, obs_units, cat_code, obs_name)) %>%
  rename(dttm=SampleDrawDtTm, imm_code=cat_code)

df_sum %>%
  group_by(imm_code) %>%
  summarise(
    Min    = min(obs_val, na.rm = TRUE),
    Q1     = quantile(obs_val, 0.25, na.rm = TRUE),
    Median = median(obs_val, na.rm = TRUE),
    Q3     = quantile(obs_val, 0.75, na.rm = TRUE),
    Max    = max(obs_val, na.rm = TRUE),
    .groups = "drop") %>%
  gt() %>%
  fmt_number(columns = -imm_code, decimals = 2) %>%
  cols_label(
    imm_code = "IMM Code",
    Min = "Min",
    Q1 = "Q1",
    Median = "Median",
    Q3 = "Q3",
    Max = "Max")

ggplot(df_sum %>% filter(imm_code=="crp"), aes(x=obs_val)) +
  geom_histogram()



#--------------------------------------------------
# Save the results
#--------------------------------------------------

saveRDS(df_sum, file=paste0(path_prep, "ampc/df_imm.RDS"))
