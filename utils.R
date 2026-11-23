#--------------------------------------------------
# Utils
#--------------------------------------------------
# This files contains all help functions for the AmpC project.
# Run this file to get the functions.

# Load the required packages
packages=c("survival", "ggplot2", "broom", "dplyr", "extrafont", "cobalt", "latex2exp")
lapply(packages, library, character.only=T)


#--------------------------------------------------
# Graphical parameters
#--------------------------------------------------

# Colorscheme used for the plots
col = c( "#dc267f", "#785ef0",  "#ffb000")
col_long=c("#648fff", "#785ef0", "#dc267f", "#fe6100", "#ffb000", "#000000", "#808080")
col_lab=c("AmpC-vulnerable"= "#dc267f",  "AmpC-stable"="#785ef0", "adequate"="#785ef0"  , "Non-beta lactam"="#ffb000", "meropenem"="#785ef0", "cefepime"="#ffb000")
line_lab=c("AmpC-vulnerable"= "dashed",  "AmpC-stable"="solid", "adequate"="solid"  , "Non-beta lactam"="solid", "meropenem"="solid", "cefepime"="solid")

# Choose the font
loadfonts(device="all", quiet=T)
par(family="LM Roman 10")
font="Times"
font_sz=16

# Plot size in cm
width=25
height=15

# Vector for adjusting the plot labels
translate=c("TRUE"="", "arm"="", "_"=" ", "sofa"="SOFA", "ips"="ICU", "ind"="inducers", "anti"="", "inad before"="AmpC-vulnerable empirical", "other treat"="other treatment", "leukozyten"="leukocytes", "t risk"="time at risk", "neutrophile"="neutrophils", "lymphozyten"="lymphocytes", "polym"="polymicrobial", "poly bsi"="polymicrobial BSI", "crp"="CRP", "pct"="PCT", "bmi"="BMI")


#--------------------------------------------------
# Help function for determining the treatment arm
#--------------------------------------------------
# This function is similar to get_arm_outcome, with the difference that we look at cabapenems vs. cefepime as arms.

# Vector with all occurring carbapenems antibiotics
carbapenems=c("meropenem", "etrapenem", "imipenem and enzyme inhibitor")

get_anti_outcome=function(t_start=hours(48), t_duration=hours(24), t_max=hours(72), t_zero="start"){
  # This function returns the unique antibiotic treatment under which each encounter was between t_start and t_max after the first positive test.
  # Arguments:
  #     t_start: duration after which we apply the treatment criteria
  #     t_max: maximal duration after the test at which the treatment can start
  #     t_duration: every 24h during t_duration, the treatment has to be administered and no other of the three arms
  # Return: dataframe containing the arm and the ICU status at the start of the treatment


  #--------------------------------------------------
  # Determine the antibiotic treatment
  #--------------------------------------------------
  carbapenems=c("meropenem", "etrapenem", "imipenem and enzyme inhibitor")

  # Collect the treatment given after dttm_test+t_start
  df_treat_after=df_treat %>% left_join(df_empirical %>% select(encounter_pseudo_id, dttm_test), by=c("encounter_pseudo_id")) %>%
    filter(dttm>=dttm_test+t_start & dttm<=dttm_test+t_max+t_duration & anti_name %in% c("cefepime", carbapenems))

  # Get the first treatment for each encounter
  df_treat_first=df_treat_after %>% filter(arm!="4" & dttm<=dttm_test+t_max) %>% summarise(dttm_first=min(dttm), .by="encounter_pseudo_id") %>%
    left_join(df_treat_after, by=join_by(encounter_pseudo_id==encounter_pseudo_id, dttm_first==dttm)) %>% filter(arm!="4") %>%
    select(dttm_first, encounter_pseudo_id, anti_name, arm) %>% unique() %>% mutate(dttm_check=dttm_first+t_duration) %>% rename(anti_def=anti_name)

  # Check if other treatments were given
  enc_mult= df_treat_after %>% filter(arm!="4") %>%
    left_join(df_treat_first, by=join_by(encounter_pseudo_id==encounter_pseudo_id, dttm>=dttm_first, dttm<=dttm_check)) %>%
    filter(anti_name!=anti_def) %>% pull(encounter_pseudo_id) %>% unique()

  # Exclude encounter that have received any other treatment arm
  df_treat_first=df_treat_first %>% filter(!(encounter_pseudo_id %in% enc_mult))
  print(paste0("Encounters excluded due to non unique treatment arm: ", as.character(length(enc_mult))))


  # Check if the treatment was given at least every 24h during t_duration
  n_int=round(as.numeric(t_duration, unit="hours") / 24, digits=0)
  enc_cont=df_treat_first %>% pull(encounter_pseudo_id) %>% unique()

  for(i in 1:max(n_int, 1)){
    # Get all encounters that received the treatment again in 24h
    enc_cont_up=df_treat_after %>% left_join(df_treat_first, by="encounter_pseudo_id") %>%
      filter(dttm>=dttm_first+(1-i)*hours(24) & dttm<=dttm_first+i*hours(24) & anti_name==anti_def) %>%
      pull(encounter_pseudo_id)
    # Keep only the encounters that have received the treatment again
    enc_cont=enc_cont[enc_cont %in% enc_cont_up]
  }

  # Keep only the encounters that have continued the treatment
  df_arm=df_treat_first %>% filter(encounter_pseudo_id %in% enc_cont) %>%
    rename(dttm_start=dttm_first) %>% rename(anti=anti_def)

  # Check if AmpC-vulnerable treatment has been given before
  enc_inad=df_treat %>% filter(arm=="1") %>% left_join(df_empirical %>% select(encounter_pseudo_id, dttm_test), by=c("encounter_pseudo_id")) %>%
    left_join(df_arm %>% select(encounter_pseudo_id, dttm_start), by="encounter_pseudo_id") %>%
    filter(dttm>=dttm_test & dttm<dttm_start) %>% pull(encounter_pseudo_id) %>% unique()

  # Check if other antibiotics have been administered
  enc_4= df_treat %>% filter(arm=="4") %>%
    left_join(df_treat_first, by=join_by(encounter_pseudo_id==encounter_pseudo_id, dttm>dttm_first, dttm<dttm_check)) %>%
    pull(encounter_pseudo_id) %>% unique()

  # Add the information
  df_arm=df_arm %>% mutate(inad_before=(encounter_pseudo_id %in% enc_inad), other_treat=(encounter_pseudo_id %in% enc_4))

  #--------------------------------------------------
  # Compute the outcome, incl. ICU admissions and current ICU status
  #--------------------------------------------------

  df_icu=df_departement %>% filter(encounter_pseudo_id %in% df_arm$encounter_pseudo_id) %>%
    group_by(encounter_pseudo_id) %>% arrange(entry_dttm) %>% mutate(discharge_dttm=lead(entry_dttm)) %>% ungroup() %>%
    filter(departement_abbreviation!="DISCHARGE")

  # Movements between ICU and normal stations
  by=join_by(encounter_pseudo_id==encounter_pseudo_id, dttm_start>=entry_dttm, dttm_start<discharge_dttm)
  df_arm_cov=df_arm %>% left_join(df_icu, by) %>% left_join(df_empirical, by=c("encounter_pseudo_id")) %>% replace_na(list(ips=F))

  # New ICU admission afterwards
  df_icu=df_departement %>% left_join(df_arm_cov %>% select(encounter_pseudo_id, dttm_start), by="encounter_pseudo_id") %>%
    filter(entry_dttm>=dttm_start) %>% filter(ips) %>% summarise(entry_ips_dttm=min(entry_dttm), .by="encounter_pseudo_id") %>%
    mutate(entry_ips_dttm=if_else(is.finite(entry_ips_dttm), entry_ips_dttm, NA))

  df_arm_cov=df_arm_cov %>% left_join(df_icu, by="encounter_pseudo_id")

  # Filter for 48h before the test
  df_imm_48=df_imm %>% left_join(df_arm_cov %>% select(encounter_pseudo_id, dttm_start), by="encounter_pseudo_id") %>%
    filter(dttm<dttm_start, dttm>dttm_start - hours(48))

  # Keep the last test, convert to wide format
  df_imm_add=df_imm_48 %>% summarise(dttm_last=max(dttm), .by=c("encounter_pseudo_id", "imm_code")) %>%
    rename(dttm=dttm_last) %>%
    right_join(df_imm_48, by=c("encounter_pseudo_id", "imm_code", "dttm")) %>%
    select(encounter_pseudo_id, imm_code, obs_val) %>%
    pivot_wider(values_from = obs_val, names_from = imm_code, values_fn=mean)

  # Add the immunology data
  df_arm_cov=df_arm_cov %>% left_join(df_imm_add, by="encounter_pseudo_id")


  # Return the summary data frame
  return(df_arm_cov)
}



#--------------------------------------------------
# Help function for plotting a weighted non-parametric estimator
#--------------------------------------------------
# Help function for plotting weighted Aalen-Johansen curves

plot_surv=function(df, title, xmax=25, legend_sum=F, caption=F, carbfep=F, km=F){
  # Rename the treatment arms
  if(carbfep){
    df=df  %>% mutate(strata=factor(ifelse(strata=="arm=cefepime", "cefepime", "meropenem"), levels = c("meropenem", "cefepime")))
  }else{
    if(legend_sum){
      df=df %>% mutate(strata=factor(ifelse(strata=="arm=1", "inadaquate", "adqaute"), levels = c("inadaquate", "adqaute")))
    }
    else{
      df=df %>% mutate(strata=factor(ifelse(strata=="arm=1", "AmpC-vulnerable", ifelse(strata=="arm=2", "AmpC-stable", "Non-beta lactam")), levels=c("AmpC-vulnerable", "AmpC-stable", "Non-beta lactam")))
    }
  }  
  
  # Add the last observation such that the we have the wanted observation interval (0, x_max).
  df_add=df %>% left_join(df %>% filter(time<=xmax) %>% summarise(min_diff=min(xmax-time), .by="strata"), by="strata") %>%
    filter(xmax-time==min_diff) %>% mutate(time=xmax) %>% select(-min_diff)
  
  df=df %>% add_row(df_add)
  
  # Make the ggplot
  if(km){
    p=ggplot(df, aes(time, estimate))+
      geom_line(aes(color=strata), linewidth=1)+
      geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=strata), alpha=0.2)+
      scale_fill_manual(values=col_lab) +
      scale_color_manual(values = col_lab)+
      xlim(0,xmax)+
      labs(y="probability", x="time (days)", title=title, fill="treatment arm", color="treatment arm")+
      theme_light()+
      theme(text=element_text(family=font, size=font_sz), axis.title=element_text(size=font_sz*3/4), axis.text.x = element_text(size=font_sz*3/4), axis.text.y=element_text(size=font_sz*3/4), plot.title = element_text(hjust=0.5, family=font, size=font_sz), plot.caption=element_text(hjust=0, size=font_sz, family=font), strip.background = element_rect(color="gray", fill="white"), strip.text = element_text(size=font_sz, color="black", family=font), legend.text = element_text(size=font_sz, family=font))
  }else{
    p=ggplot(df %>% filter(state!="(s0)"), aes(time, estimate))+
      geom_line(aes(color=as.factor(strata), linetype = as.factor(strata)), linewidth=1, show.legend = T)+
      geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=strata), alpha=0.2, show.legend = T)+
      scale_fill_manual(values=col_lab, breaks=c("AmpC-vulnerable", "AmpC-stable", "Non-beta lactam", "meropenem", "cefepime"), labels=c("AmpC-vulnerable", "AmpC-stable", "Non-beta lactam", "meropenem", "cefepime"), drop=F, name=NULL) +
      scale_color_manual(values = col_lab, breaks=c("AmpC-vulnerable", "AmpC-stable", "Non-beta lactam", "meropenem", "cefepime"), labels=c("AmpC-vulnerable", "AmpC-stable", "Non-beta lactam", "meropenem", "cefepime"), drop=F, name=NULL)+
      scale_linetype_manual(values=line_lab, breaks=c("AmpC-vulnerable", "AmpC-stable", "Non-beta lactam", "meropenem", "cefepime"), labels=c("AmpC-vulnerable", "AmpC-stable", "Non-beta lactam", "meropenem", "cefepime"), drop=F, name=NULL)+
      xlim(0,xmax)+
      labs(y="probability", x="time (days)", title=title, fill="treatment arm", color="treatment arm")+
      facet_grid(rows=vars(state), scales="free_y")+
      theme_light()+
      theme(text=element_text(family=font, size=font_sz),  axis.title=element_text(size=font_sz*3/4), axis.text.x = element_text(size=font_sz*3/4), axis.text.y=element_text(size=font_sz*3/4), plot.title = element_text(hjust=0.5, family=font, size=font_sz), plot.caption=element_text(hjust=0, size=font_sz, family=font), strip.background = element_rect(color="gray", fill="white"), strip.text = element_text(size=font_sz, color="black", family=font), legend.text = element_text(size=font_sz, family=font))
    
    }
    
  if(caption){
    p=p+  labs(caption="Aalen-Johansen estimators for the competing risks discharge, ICU admission, and death adjusted for \nconfounding by patient demographics, comorbidities, and infection sites using inverse probability \nweighting (IPWT). IPWT leads to a pairwise balanced sample with absolute standardized mean \ndifferences <0.1. There is a clear difference for the discharge and ICU admission probabilities for \nDNA inhibitor antibiotics compared to carbapenems and cefepime.")
  }
  
  return(p)
}


#--------------------------------------------------
# Help function for plotting glms results
#--------------------------------------------------

plot_glm=function(df_plot, title=""){
  # This function produces a ggplot with the coefficient and confidence intervals of a glm.
  # Arguments:
  #     df_plot: data frame containing the columns var with variablenames, coef the coefficients and confidence interval coef_l, coef_u.
  # Returns:
  #     p: ggplot
  
  # Rename the variables
  df_plot$var=df_plot$var %>% str_replace_all(pattern=translate)
  
  p=ggplot(df_plot, aes(y=coef, x=factor(var, levels=var %>% rev()))) + 
    # theme(text=element_text(family=font), axis.title.x =element_text(margin=margin(0,20,0,0), family = font))+
    # theme(axis.title.x = element_blank())+
    geom_hline(yintercept=1, lty=2, lwd=1, colour=col[1]) +
    geom_errorbar(aes(ymin=coef_l, ymax=coef_u), 
                  lwd=1, width=0.5,  position = position_dodge(width = 0.9), colour=col[2]) + 
    geom_point(size=3, pch=16, position = position_dodge(width = 0.9), colour=col[2]) + 
    geom_text(data = subset(df_plot, pvalue < 0.05), aes(label = "*"), position = position_dodge(width = 0.9), vjust = 0.1, hjust=-1, color = "black", size = 4)+
    coord_flip(ylim=c(min(df_plot$coef_l), max(df_plot$coef_u))) +
    theme_light()+
    theme(text=element_text(family=font, size=font_sz),  axis.title=element_text(size=font_sz*3/4), axis.text.x = element_text(size=font_sz*3/4), axis.text.y=element_text(size=font_sz*3/4), plot.title = element_text(hjust=0.5, family=font))+
    labs(y="odds ratio", x="", title = title) 
  
  return(p)
  
}


#--------------------------------------------------
# Help function for plotting glms: adjusted and unadjusted
#--------------------------------------------------

plot_glm_grp=function(df_plot, title="Treatment Effects", xlab="odds ratio"){
  # This function produces a ggplot with the coefficient and confidence intervals of a glm.
  # Arguments:
  #     df_plot: data frame containing the columns var with variable names, coef the coefficients and confidence interval coef_l, coef_u, adj: T/F.
  # Returns:
  #     p: ggplot
  
  
  p=ggplot(df_plot, aes(y=coef, x=factor(group, levels=group %>% rev()), colour = factor(adj, levels = c(FALSE, TRUE)))) + 
    # theme(text=element_text(family=font), axis.title.x =element_text(margin=margin(0,20,0,0), family = font))+
    # theme(axis.title.x = element_blank())+
    geom_hline(yintercept=1, lwd=1) + 
    geom_linerange( aes(ymin = coef_l,ymax = coef_u, linetype = factor(adj, levels = c(FALSE, TRUE)) ), linewidth = 1, position = position_dodge(width = 0.9), key_glyph = "path") +
    geom_crossbar( aes(ymin = coef_l, ymax = coef_l,y = coef_l),width = 0.5,linewidth = 1,linetype = "solid",fatten = 0,position = position_dodge(width = 0.9), show.legend = FALSE) +
    geom_crossbar(aes(ymin = coef_u, ymax = coef_u,y = coef_u ), width = 0.5,linewidth = 1,linetype = "solid",fatten = 0,position = position_dodge(width = 0.9), show.legend = FALSE)+
    geom_point(aes(shape=factor(adj, levels = c(FALSE, TRUE))), size=3, position = position_dodge(width = 0.9)) + 
    geom_text(data = subset(df_plot, pvalue < 0.05), aes(label = "*"), position = position_dodge(width = 0.9), vjust = 0.1, hjust=-1, color = "black", size = 4)+
    scale_x_discrete(labels=c("discharge ( unadjusted )"="discharge", "discharge ( adjusted )"="discharge", "death ( unadjusted )"="death", "death ( adjusted )"="death", "icu ( unadjusted )"="icu", "icu ( adjusted )"="icu"))+
    coord_flip(ylim=c(min(df_plot$coef_l), max(df_plot$coef_u))) +
    theme_light()+
    theme(,text=element_text(family=font, size=font_sz),  axis.title=element_text(size=font_sz*3/4), axis.text.x = element_text(size=font_sz*3/4), axis.text.y=element_text(size=font_sz*3/4), plot.title = element_text(hjust=0.5, family=font, size=font_sz))+
    scale_color_manual(name = "",values = c("FALSE" = col[1], "TRUE" = col[2]), labels = c("FALSE" = "unadjusted", "TRUE" = "adjusted"))+
    scale_shape_manual(name="", values=c("FALSE"=18, "TRUE"=19), , labels = c("FALSE" = "unadjusted", "TRUE" = "adjusted"))+
    scale_linetype_manual(name="", values=c("FALSE"="dashed", "TRUE"="solid") , labels = c("FALSE" = "unadjusted", "TRUE" = "adjusted"))+
    labs(y=xlab, x="", title = title) 
  
  return(p)
  
}


#--------------------------------------------------
# Help function for plotting SMD
#--------------------------------------------------

plot_smd=function(df_ipw, title="Sample balance", match=F, abs=F, cont=F){
  # Function for plotting the standardized mean difference for the weighted and unweighted sample
  # Arguments:
  #   df_ipw: output the matching or weighting
  #   title: string for the title of the plot
  #   match: indicator for preprocessed dataframe obtained by matching
  #   abs: indicator if absolute value should be plotted
  # Return:
  #   p: ggplot
                    
  facet=c("CRP"="severity", "leukocytes"="severity", "charlson"="severity",   "ICU"="severity",  "BMI"="demographics", "age"="demographics", "sex W"="demographics", 
          "klebsiella aerogenes"="species", "providencia spp"="species", "enterobacter cloacae"="species", "citrobacter freundii"="species","morganella morganii"="species",  "serratia marcescens"="species",
          "installation"="infection site", "pus"="infection site", "invasive not specified"="infection site",  "bones joints"="infection site", "upper respiratory tract"="infection site", "abdomen"="infection site", "lower respiratory tract"="infection site", "blood"="infection site",  "urine"="infection site", "wound"="infection site", "other treatment"="infection site",
          "time at risk"="", "PCT"="severity", "neutrophils"="severity", "SOFA"="severity", "lymphocytes"="severity", "il 6"="severity", "polymicrobial"="species", "polymicrobial BSI"="infection site", "other"="infection site", "AmpC-vulnerable empirical"="demographics")
  
  # Prepare the dataframe for plotting
  if(!match & !cont){
    df_plot=bal.tab(df_ipw, un=T)
    if(length(df_plot)==3){
      df_plot=bal.tab(df_ipw, un=T)$Balance %>% rename(unadjusted='Diff.Un', adjusted='Diff.Adj')
    }else{
      df_plot=bal.tab(df_ipw, un=T)$Balance %>% rename(unadjusted='Max.Diff.Un', adjusted='Max.Diff.Adj')
    }
  }else if(cont){
    df_plot=bal.tab(df_ipw, un=T)$Balance %>% rename(unadjusted='Corr.Un', adjusted='Corr.Adj')
  }
  else{
    df_plot=df_ipw
  }
  
  # Renaming and filtering of rows
  df_plot$var_name=rownames(df_plot)
  df_plot=df_plot %>% select(-Type) %>% pivot_longer(cols=c(unadjusted, adjusted), names_to="group", values_to="SMD") %>%
    filter(!(var_name %in% c("prop.score", "distance")))
  
  # Adjust the variable names and add the groups
  df_plot$var_name=df_plot$var_name %>% str_replace_all(pattern=translate) 
  df_plot=df_plot %>% mutate(var_group=facet[var_name])
  
  # Absolute values
  x_tresh=c(-0.1, 0.1)
  if(abs){
    df_plot=df_plot %>% mutate(SMD=abs(SMD))
    x_tresh=c(0.1)
  }
  
  # Plotting
  p=ggplot(df_plot, aes(SMD, var_name))+
    geom_point(aes(colour=factor(group), shape=factor(group)), size=3)+
    geom_vline(xintercept = c(0))+
    geom_vline(xintercept = x_tresh, linetype="dashed")+
    facet_grid(
      rows = vars(var_group),
      scales = "free_y",
      space  = "free_y",
      switch = "y"
    ) +
    theme_light()+
    theme(text=element_text(family=font, colour = "black", size=font_sz*3/4),axis.text.x = element_text(size=font_sz*3/4), axis.text.y=element_text(size=font_sz*3/4), plot.title = element_text(hjust=0.5, family=font, size=font_sz),
          strip.background = element_blank(),  axis.title=element_text(size=font_sz*3/4), strip.text = element_text(colour="black", size=font_sz*3/4, face="bold"), strip.placement = "outside", 
          panel.spacing=unit(0, "lines"), strip.clip="off")+
    scale_color_manual(values=c(col[2], col[1]))+
    scale_shape_manual(values=c(19, 18))+
    labs(x=if(!cont){TeX("SMD = \\tiny $(m_{adequate} - m_{AmpC-vulnerable})/sd$")}else{"Correlation"} , y="", title=title, colour="", shape="")
  
  return(list(plot=p, df=df_plot))
}



#--------------------------------------------------
# Help function for plotting SMD
#--------------------------------------------------

get_finegray=function(data, weights=NULL, arm=c("2"), etypes=c("discharge", "death"), imm_var_rest=F){
  # Function for fitting the Fine-Gray model, adjusted and unadjusted
  # Arguments:
  #   data: dataframe with data on which the model should be fitted
  #   weights: inverse probability weights, e.g. from the WeightIt function
  #   etypes: vector with event types
  # Return:
  #   df: dataframe with the estimation and confidence intervals
  
  # Formula for the Fine-Gray model
  var_exp=c(materials[materials %in% colnames(data)], species, patient_encounter[patient_encounter!="ips"], ifelse(imm_var_rest, imm_var_gen, imm_var), "other_treat", "polym", "poly_bsi", "inad_before")
  formula_fg=as.formula(paste0("Surv(fgstart, fgstop, fgstatus)~arm+", paste(var_exp, collapse="+")))
  
  # Intialize dataframe
  df=data.frame('coef'=double(), 'lower.95'=double(), 'upper.95'=double(), 'adj'=logical(), 'event'=character(), 'pvalue'=double(), 'arm'=character())
  
  for(evn in etypes){
    for(adj in c(F, T)){
      # Fit an adjusted Fine-Gray model 
      fdata <- finegray(as.formula(paste0("Surv(time, event) ~ arm+encounter_pseudo_id+", paste(var_exp, collapse = "+"))), 
                        data = if(adj & any(weights==0)){data[as.logical(weights), ]}else{data}, etype=evn, weight=if(adj & !any(weights==0)){weights}else{NULL})
      fgfit <- coxph(formula_fg, data = fdata, weight = fgwt, id=encounter_pseudo_id, control=coxph.control(iter.max=200), robust=T)
      
      # Add the results to the dataframe df
      for(a in arm){
        df_add=summary(fgfit)$conf.int[paste0("arm", a), c("exp(coef)", "lower .95", "upper .95")] %>%
          t() %>% as.data.frame() %>% mutate(event=evn, adj=adj)  %>% rename(coef='exp(coef)') 
        df_add=df_add %>% mutate(pvalue=summary(fgfit)$coefficients[paste0("arm", a), "Pr(>|z|)"], arm=paste0("arm", a))
        colnames(df_add)=colnames(df_add) %>% gsub(" ", "", .)
        df=df %>% add_row(df_add)
      }
    }
  }
  
  return(df)
}



#--------------------------------------------------
# Help function for determine the AmpC like resistances
#--------------------------------------------------
 
# AmpC=function(cro, caz){
#   # Function for summarizing the AmpC expression
#   
#   anti_codes=c("cro", "caz")
#   for(code in anti_codes){
#     assign(code, case_when(is.na(get(code)) ~  "S", T ~  get(code)))
#   }
#   
#   ampc_n=(cro=="R")+(caz=="R")
#   ampc=ampc_n>=1
#   
#   return(tibble(ampc_n=ampc_n, ampc=ampc))
# }



get_arm_outcome_2=function(t_start=hours(48), t_duration=hours(24), t_max=hours(72), t_zero="start", route="all"){
  # This function returns the unique treatment arm under which each encounter was between t_start and t_max after the first positive test.
  # Arguments:
  #     t_start: duration after which we apply the treatment criteria
  #     t_max: maximal duration after the test at which the treatment can start
  #     t_duration: every 24h during t_duration, the treatment has to be administered and no other of the three arms
  # Return: dataframe containing the arm, the ICU status and SOFA score at the start of the treatment
  
  
  #--------------------------------------------------
  # Determine the treatment arm
  #--------------------------------------------------
  
  # Collect the treatment given after dttm_test+t_start
  df_treat_after=df_treat %>% left_join(df_empirical %>% select(encounter_pseudo_id, dttm_test), by=c("encounter_pseudo_id")) %>%
    filter(dttm>=dttm_test+t_start & dttm<=dttm_test+t_max+t_duration)
  
  # Get the first treatment for each encounter
  df_treat_first=df_treat_after %>% filter(arm!="4" & dttm<=dttm_test+t_max) %>% summarise(dttm_first=min(dttm), .by="encounter_pseudo_id") %>% 
    left_join(df_treat_after, by=join_by(encounter_pseudo_id==encounter_pseudo_id, dttm_first==dttm)) %>% filter(arm!="4") %>%
    select(dttm_first, encounter_pseudo_id, arm) %>% unique() %>% mutate(dttm_check=dttm_first+t_duration) %>% rename(arm_def=arm)
  
  # Check if other treatments were given
  df_check= df_treat_after %>% filter(arm!="4") %>%
    left_join(df_treat_first, by=join_by(encounter_pseudo_id==encounter_pseudo_id, dttm>=dttm_first, dttm<=dttm_check)) 
  
  enc_mult= df_check %>%
    filter(arm!=arm_def) %>% pull(encounter_pseudo_id) %>% unique()
  
  # Check if the route is correct
  if(route=="oral"){
    enc_route=df_check %>% filter(arm==arm_def) %>% filter(!oral) %>%
      pull(encounter_pseudo_id) %>% unique()
  }else if(route=="iv"){
    enc_route=df_check %>% filter(arm==arm_def) %>% filter(oral)  %>%
      pull(encounter_pseudo_id) %>% unique()
  }else{
    enc_route=c()
  }
  
  # Exclude encounter that have received any other treatment arm
  df_treat_first=df_treat_first %>% filter(!(encounter_pseudo_id %in% enc_mult) & !(encounter_pseudo_id %in% enc_route))
  print(paste0("Encounters excluded due to non unique treatment arm: ", as.character(length(enc_mult)), " \n Encounters exclude due to administration route: ", as.character(length(enc_route))))
  
  # Check if the treatment was given at least every 24h during t_duration
  n_int=round(as.numeric(t_duration, unit="hours") / 24, digits=0)
  enc_cont=df_treat_first %>% pull(encounter_pseudo_id) %>% unique()
  
  for(i in 1:max(n_int, 1)){
    # Get all encounters that received the treatment again in 24h
    enc_cont_up=df_treat_after %>% left_join(df_treat_first, by="encounter_pseudo_id") %>% 
      filter(dttm>=dttm_first+(1-i)*hours(24) & dttm<=dttm_first+i*hours(24) & arm==arm_def) %>% 
      pull(encounter_pseudo_id)
    # Keep only the encounters that have received the treatment again
    enc_cont=enc_cont[enc_cont %in% enc_cont_up]
  }
  
  # Keep only the encounters that have continued the treatment
  df_arm=df_treat_first %>% filter(encounter_pseudo_id %in% enc_cont) %>%
    rename(dttm_start=dttm_first) %>% rename(arm=arm_def)
  
  
  
  # Check if AmpC-vulnerable treatment has been given before
  enc_inad=df_treat %>% filter(arm=="1") %>% left_join(df_empirical %>% select(encounter_pseudo_id, dttm_test), by=c("encounter_pseudo_id")) %>%
    left_join(df_arm %>% select(encounter_pseudo_id, dttm_start), by="encounter_pseudo_id") %>%
    filter(dttm>=dttm_test & dttm<dttm_start) %>% pull(encounter_pseudo_id) %>% unique()
  
  # Check if other antibiotics have been administered
  enc_4= df_treat_after %>% filter(arm=="4") %>%
    left_join(df_treat_first, by=join_by(encounter_pseudo_id==encounter_pseudo_id, dttm>dttm_first, dttm<dttm_check)) %>%
    pull(encounter_pseudo_id) %>% unique()
  
  # Add the information
  df_arm=df_arm %>% mutate(inad_before=(encounter_pseudo_id %in% enc_inad), other_treat=(encounter_pseudo_id %in% enc_4))
  
  
  #--------------------------------------------------
  # Compute the outcome, incl. ICU admissions and current ICU status
  #--------------------------------------------------
  
  df_icu=df_departement %>% filter(encounter_pseudo_id %in% df_arm$encounter_pseudo_id) %>%
    group_by(encounter_pseudo_id) %>% arrange(entry_dttm) %>% mutate(discharge_dttm=lead(entry_dttm)) %>% ungroup() %>%
    filter(departement_abbreviation!="DISCHARGE")
  
  # Movements between ICU and normal stations
  by=join_by(encounter_pseudo_id==encounter_pseudo_id, dttm_start>=entry_dttm, dttm_start<discharge_dttm)
  df_arm_cov=df_arm %>% left_join(df_icu, by) %>% left_join(df_empirical, by=c("encounter_pseudo_id")) %>% replace_na(list(ips=F))
  
  # New ICU admission afterwards
  df_icu=df_departement %>% left_join(df_arm_cov %>% select(encounter_pseudo_id, dttm_start), by="encounter_pseudo_id") %>%
    filter(entry_dttm>=dttm_start) %>% filter(ips) %>% summarise(entry_ips_dttm=min(entry_dttm), .by="encounter_pseudo_id") %>%
    mutate(entry_ips_dttm=if_else(is.finite(entry_ips_dttm), entry_ips_dttm, NA))
  
  # Combine the dataframes
  df_arm_cov=df_arm_cov %>% left_join(df_icu, by="encounter_pseudo_id")
  
  
  #--------------------------------------------------
  # Add the immunology data
  #--------------------------------------------------
  
  # Filter for 72h before the test
  df_imm_48=df_imm %>% left_join(df_arm_cov %>% select(encounter_pseudo_id, dttm_start), by="encounter_pseudo_id") %>%
    filter(dttm<dttm_start, dttm>dttm_start - hours(72))
  
  # Keep the last test, convert to wide format
  df_imm_add=df_imm_48 %>% summarise(dttm_last=max(dttm), .by=c("encounter_pseudo_id", "imm_code")) %>%
    rename(dttm=dttm_last) %>%
    right_join(df_imm_48, by=c("encounter_pseudo_id", "imm_code", "dttm")) %>% 
    select(encounter_pseudo_id, imm_code, obs_val) %>% 
    pivot_wider(values_from = obs_val, names_from = imm_code, values_fn=mean) 
  
  # Add the immunology data
  df_arm_cov=df_arm_cov %>% left_join(df_imm_add, by="encounter_pseudo_id")
  
  
  #--------------------------------------------------
  # Add the SOFA score
  #--------------------------------------------------
  
  # Filter for 72h before the test
  df_sofa_48=df_sofa %>% left_join(df_arm_cov %>% select(encounter_pseudo_id, dttm_start), by="encounter_pseudo_id") %>%
    filter(dttm<=dttm_start+hours(2), dttm>dttm_start - days(14))
  
  # Keep the last value
  df_sofa_add=df_sofa_48 %>% summarise(dttm_last=max(dttm), .by=c("encounter_pseudo_id")) %>%
    rename(dttm=dttm_last) %>%
    right_join(df_sofa_48, by=c("encounter_pseudo_id", "dttm")) %>% 
    select(encounter_pseudo_id, sofa) %>% summarise(sofa=mean(sofa), .by="encounter_pseudo_id")

  # Add the sofa score to the dataframe
  df_arm_cov=df_arm_cov %>% left_join(df_sofa_add, by="encounter_pseudo_id")
  
  
  # Return the summary data frame
  return(df_arm_cov)
}


get_arm_outcome_def=function(t_grace=hours(24), t_duration=hours(48), t_max=days(30), route="all"){
  # This function returns the unique treatment arm under which each encounter was between t_start and t_max after the first positive test.
  # Arguments:
  #     t_grace: how much time after the susceptibility test is allowed to switch treatment
  #     t_max: maximal duration after the test at which the treatment can start
  #     t_duration: every 24h during t_duration, the treatment has to be administered and no other of the three arms
  # Return: dataframe containing the arm, the ICU status and SOFA score at the start of the treatment
  
  
  #--------------------------------------------------
  # Determine the treatment arm
  #--------------------------------------------------
  
  # Collect the treatment given after dttm_result
  df_treat_after=df_treat %>% left_join(df_empirical %>% select(all_of(c("encounter_pseudo_id", "dttm_test", "dttm_result"))), by=c("encounter_pseudo_id")) %>%
    filter(dttm>=dttm_result & dttm<=dttm_result+t_grace+t_duration & dttm<dttm_test+t_max)
  
  # Get the first treatment for each encounter
  df_treat_first=df_treat_after %>% filter(arm!="4" & dttm<=dttm_test+t_max) %>% summarise(dttm_first=min(dttm), .by="encounter_pseudo_id") %>% 
    left_join(df_treat_after, by=join_by(encounter_pseudo_id==encounter_pseudo_id, dttm_first==dttm)) %>% filter(arm!="4") %>%
    select(dttm_first, encounter_pseudo_id, arm) %>% unique() %>% mutate(dttm_check=dttm_first+t_duration) %>% rename(arm_def=arm)
  
  # Check if other treatments were given
  df_check= df_treat_after %>% filter(arm!="4") %>%
    left_join(df_treat_first, by=join_by(encounter_pseudo_id==encounter_pseudo_id, dttm>=dttm_first, dttm<=dttm_check)) 
  
  enc_mult= df_check %>%
    filter(arm!=arm_def) %>% pull(encounter_pseudo_id) %>% unique()
  
  # Check if the route is correct
  if(route=="oral"){
    enc_route=df_check %>% filter(arm==arm_def) %>% filter(!oral) %>%
      pull(encounter_pseudo_id) %>% unique()
  }else if(route=="iv"){
    enc_route=df_check %>% filter(arm==arm_def) %>% filter(oral)  %>%
      pull(encounter_pseudo_id) %>% unique()
  }else{
    enc_route=c()
  }
  
  # Exclude encounter that have received any other treatment arm
  df_treat_first=df_treat_first %>% filter(!(encounter_pseudo_id %in% enc_mult) & !(encounter_pseudo_id %in% enc_route))
  print(paste0("Encounters excluded due to non unique treatment arm: ", as.character(length(enc_mult)), " \n Encounters exclude due to administration route: ", as.character(length(enc_route))))
  
  # Check if the treatment was given at least every 24h during t_duration
  n_int=round(as.numeric(t_duration, unit="hours") / 24, digits=0)
  enc_cont=df_treat_first %>% pull(encounter_pseudo_id) %>% unique()
  
  for(i in 1:max(n_int, 1)){
    # Get all encounters that received the treatment again in 24h
    enc_cont_up=df_treat_after %>% left_join(df_treat_first, by="encounter_pseudo_id") %>% 
      filter(dttm>=dttm_first+(1-i)*hours(24) & dttm<=dttm_first+i*hours(24) & arm==arm_def) %>% 
      pull(encounter_pseudo_id)
    # Keep only the encounters that have received the treatment again
    enc_cont=enc_cont[enc_cont %in% enc_cont_up]
  }
  
  # Keep only the encounters that have continued the treatment
  df_arm=df_treat_first %>% filter(encounter_pseudo_id %in% enc_cont) %>%
    rename(dttm_start=dttm_first) %>% rename(arm=arm_def)
  
  
  
  # Check if AmpC-vulnerable treatment has been given before
  enc_inad=df_treat %>% filter(arm=="1") %>% left_join(df_empirical %>% select(encounter_pseudo_id, dttm_test), by=c("encounter_pseudo_id")) %>%
    left_join(df_arm %>% select(encounter_pseudo_id, dttm_start), by="encounter_pseudo_id") %>%
    filter(dttm>=dttm_test & dttm<dttm_start) %>% pull(encounter_pseudo_id) %>% unique()
  
  # Check if other antibiotics have been administered
  enc_4= df_treat_after %>% filter(arm=="4") %>%
    left_join(df_treat_first, by=join_by(encounter_pseudo_id==encounter_pseudo_id, dttm>dttm_first, dttm<dttm_check)) %>%
    pull(encounter_pseudo_id) %>% unique()
  
  # Add the information
  df_arm=df_arm %>% mutate(inad_before=(encounter_pseudo_id %in% enc_inad), other_treat=(encounter_pseudo_id %in% enc_4))
  
  
  #--------------------------------------------------
  # Compute the outcome, incl. ICU admissions and current ICU status
  #--------------------------------------------------
  
  df_icu=df_departement %>% filter(encounter_pseudo_id %in% df_arm$encounter_pseudo_id) %>%
    group_by(encounter_pseudo_id) %>% arrange(entry_dttm) %>% mutate(discharge_dttm=lead(entry_dttm)) %>% ungroup() %>%
    filter(departement_abbreviation!="DISCHARGE")
  
  # Movements between ICU and normal stations
  by=join_by(encounter_pseudo_id==encounter_pseudo_id, dttm_start>=entry_dttm, dttm_start<discharge_dttm)
  df_arm_cov=df_arm %>% left_join(df_icu, by) %>% left_join(df_empirical, by=c("encounter_pseudo_id")) %>% replace_na(list(ips=F))
  
  # New ICU admission afterwards
  df_icu=df_departement %>% left_join(df_arm_cov %>% select(encounter_pseudo_id, dttm_start), by="encounter_pseudo_id") %>%
    filter(entry_dttm>=dttm_start) %>% filter(ips) %>% summarise(entry_ips_dttm=min(entry_dttm), .by="encounter_pseudo_id") %>%
    mutate(entry_ips_dttm=if_else(is.finite(entry_ips_dttm), entry_ips_dttm, NA))
  
  # Combine the dataframes
  df_arm_cov=df_arm_cov %>% left_join(df_icu, by="encounter_pseudo_id")
  
  
  #--------------------------------------------------
  # Add the immunology data
  #--------------------------------------------------
  
  # Filter for 72h before the test
  df_imm_48=df_imm %>% left_join(df_arm_cov %>% select(encounter_pseudo_id, dttm_start), by="encounter_pseudo_id") %>%
    filter(dttm<dttm_start, dttm>dttm_start - hours(72))
  
  # Keep the last test, convert to wide format
  df_imm_add=df_imm_48 %>% summarise(dttm_last=max(dttm), .by=c("encounter_pseudo_id", "imm_code")) %>%
    rename(dttm=dttm_last) %>%
    right_join(df_imm_48, by=c("encounter_pseudo_id", "imm_code", "dttm")) %>% 
    select(encounter_pseudo_id, imm_code, obs_val) %>% 
    pivot_wider(values_from = obs_val, names_from = imm_code, values_fn=mean) 
  
  # Add the immunology data
  df_arm_cov=df_arm_cov %>% left_join(df_imm_add, by="encounter_pseudo_id")
  
  
  #--------------------------------------------------
  # Add the SOFA score
  #--------------------------------------------------
  
  # Filter for 72h before the test
  df_sofa_48=df_sofa %>% left_join(df_arm_cov %>% select(encounter_pseudo_id, dttm_start), by="encounter_pseudo_id") %>%
    filter(dttm<=dttm_start+hours(2), dttm>dttm_start - days(14))
  
  # Keep the last value
  df_sofa_add=df_sofa_48 %>% summarise(dttm_last=max(dttm), .by=c("encounter_pseudo_id")) %>%
    rename(dttm=dttm_last) %>%
    right_join(df_sofa_48, by=c("encounter_pseudo_id", "dttm")) %>% 
    select(encounter_pseudo_id, sofa) %>% summarise(sofa=mean(sofa), .by="encounter_pseudo_id")
  
  # Add the sofa score to the dataframe
  df_arm_cov=df_arm_cov %>% left_join(df_sofa_add, by="encounter_pseudo_id")
  
  
  # Return the summary data frame
  return(df_arm_cov)
}


