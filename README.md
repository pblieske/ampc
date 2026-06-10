# Inadequate treatment of chromosomal AmpC Enterobacterales infections increases the risk of phenotypic resistance despite similar clinical outcomes

This is the GitHub repository of the paper [title](url), containing all the R scripts used to obtain the presented results. The data can not be provided due to ethical reasons.

## Abstract

**Background:** Exposure to certain β-lactams can induce overexpression of chromosomally-encoded AmpC β-lactamase in ESCPM bacteria (Enterobacter cloacae, Klebsiella aerogenes, Serratia marcescens, Citrobacter freundii, Providencia spp., and Morganella morganii), decreasing susceptibility. While this mechanism is well understood, clinical impact and extent to which treatment drives resistance remains unclear due to lack of randomized controlled trials and confounding in observational studies.
**Methods:** We emulate target trials using electronic health records of adult inpatients with ESCPM infections from a tertiary-care center, comparing three strategies: (i) inadequate antibiotics (by AmpC hydrolysable β-lactams), (ii) carbapenems/cefepime, and (iii) DNA inhibitors. Outcomes include length of stay, in-hospital mortality, microbiological relapse, and phenotypic resistance. We adjusted for confounding by patient demographics, comorbidities, infection site, and severity using inverse probability weighting (IPW).
**Results:** Among 4’228 eligible encounters, 961 with repeated cultures were included in the resistance analysis. After IPW, covariate balance was achieved (standardized mean differences < 0.1). Inadequate therapy was associated with increased odds of new AmpC-like resistance (adjusted odds ratio (aOR) 2.03; 95% CI, 1.12 – 3.73) and persistence of resistance (aOR 1.85, 95% CI 1.02 – 3.36). In a stratified analysis, meropenem was associated with lower odds of reduced susceptibility compared with cefepime (aOR 0.20, 95% CI 0.06 – 0.63). Clinical outcomes – mortality, ICU admission, and length of stay – did not differ substantially, and microbiological relapse was not associated with treatment strategy.
**Conclusion:** The findings suggest that while short-term clinical outcomes may be similar, inappropriate use of β-lactams is associated with decreased susceptibility, emphasizing the need of antimicrobial stewardship.

## Usage

The main statistical analysis is performed in the Quarto file target_trial_emulation.qmd. The utils.R contains all needed help functions. The remaining files were used for preprocessing the data after collecting it from the hospital servers.