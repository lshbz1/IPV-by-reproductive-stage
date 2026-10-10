This repository contains the R code used for the analyses and figures examining intimate partner violence (IPV) across female reproductive stages.

The analytical workflow includes data preparation, pooled individual-level analyses, parity-specific analyses, sensitivity analyses using more finely defined reproductive stages, analyses of factors associated with postpartum IPV burden, and generation of the main figures.

## Repository structure
The repository contains seven R scripts:

1_Main analysis.R: Performs data preparation and the primary individual-level pooled analysis of IPV across reproductive stages. This script constructs the main analytical variables, defines reproductive-stage groups and IPV outcomes, and fits the primary regression models. 
2_Analysis by parity.R: Conducts parity-specific analyses to examine whether associations between reproductive stage and IPV differ according to the number of previous children. 
3_Sensitivity analysis.R: Conducts sensitivity analyses using more finely classified reproductive-stage and time-since-childbirth categories to assess the robustness and temporal patterning of the main findings. 
4_Risk factor analysis.R: Examines individual- and country-level factors associated with postpartum IPV burden using multivariable regression models. 
Fig 2.R: Generates Figure 2 from the corresponding analysis outputs. 
Fig 3.R: Generates Figure 3 from the corresponding analysis outputs. 
Fig 4.R: Generates Figure 4 from the corresponding analysis outputs. 

## Data source
The analyses use data from the Demographic and Health Surveys (DHS).
DHS microdata are not included in this repository because access is subject to the DHS Program's data-use requirements. Researchers wishing to reproduce the analyses should request access to the relevant surveys directly from the DHS Program.
