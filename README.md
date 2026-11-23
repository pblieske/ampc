# AmpC-vulnerable β-lactam therapy of chromosomal AmpC Enterobacterales infections increases the risk of phenotypic resistance despite similar clinical outcomes: Emulation of target trials

This is the GitHub repository contains all the R scripts used to obtain the presented results. The data can not be provided due to ethical reasons.

## Abstract

**Background:** Exposure to certain β-lactams can induce overexpression of chromosomally-encoded AmpC β-lactamases in ESCPM bacteria (*Enterobacter cloacae*, *Klebsiella aerogenes*, *Serratia marcescens*, *Citrobacter freundii*, *Providencia spp.*, and *Morganella morganii*), decreasing susceptibility. While this mechanism is well understood, clinical impact and extent to which treatment drives resistance remains unclear due to lack of randomized controlled trials and confounding in observational studies.

**Methods:** We emulate target trials using electronic health records of adult inpatients with ESCPM infections from a tertiary-care centre, comparing three strategies: (i) AmpC-vulnerable β-lactams ( hydrolysable by AmpC), (ii) AmpC-stable β-lactams, and (iii) non-β-lactam antibiotics. Outcomes include length of stay, in-hospital mortality, microbiological relapse, and phenotypic resistance. We adjusted for confounding by demographics, comorbidities, infection site, severity, and testing using inverse probability weighting (IPW).

**Results:**  Among 3’962 eligible encounters, 1’170 had follow-up cultures within 30 days. IPW achieved covariate balance (maximum standardized mean differences 0.14). AmpC-vulnerable therapy was associated with increased odds of new AmpC-compatible resistance (adjusted odds ratio (aOR) 4.74, 95% CI: 2.32 – 9.66), but not with persistence of resistance (aOR 0.78, 95% CI 0.36 – 1.72). In a stratified analysis, meropenem was associated with lower odds of reduced susceptibility compared with cefepime (aOR 0.01, 95% CI: 0.004 − 0.043). Clinical outcomes – mortality, ICU admission, and length of stay – did not differ substantially, and microbiological relapse was not associated with treatment strategy.

**Conclusion:** The findings suggest that while short-term clinical outcomes may be similar, inappropriate use of β-lactams is associated with decreased susceptibility, emphasising the need of antimicrobial stewardship.


## Usage

The main statistical analysis is performed in the Quarto file target_trial_emulation.qmd. The utils.R contains all needed help functions. The remaining files were used for preprocessing the data after collecting it from the hospital servers.
