---
title: ""
output: html_document
date: "2024-01-31"
author: "Hélène Cecilia"
---

====================

This repository contains code used in the following paper.

Cecilia H, Althouse BM, Azar SR, Moehn B, Yun R, Rossi SL, Vasilakis N, Hanley KA (2024) ***Aedes albopictus* is not an arbovirus aficionado - Impacts of sylvatic flavivirus infection in vectors and hosts on mosquito engorgement on non-human primates** *bioRxiv*

All code contained within this repository is released under the [CRAPL v0.1 License](http://matt.might.net/articles/crapl/). Data provided in this repository are sufficient to rerun all analyses.

The analyses were all performed on a desktop computer (Ubuntu 22.04)

====================

### Set up of code base:

-   `script` contains code to analyze the data, perform calculations, and generate output and figures
-   `data` contains experimental data described in the manuscript and used by scripts
-   `output` contains results from analyses and figures

### Script folder

*Note : some scripts are in the quarto format (.qmd), which is similar to R Markdown (see [quarto.org](https://quarto.org/)). It is for long scripts that benefit from being executed in chunks, but the programming language used is still R.*

-   `Data_Viz.R` produces **Figures 4 and S1**, related to host body temperature, engorgement, and time of day. 

-   `Day0_Engorgement_Vector_Infection_Status.R` tests the effect of vector infection status on engorgement rates. This script produces **Figure 2**. The selected model is saved in `output/engorgement_models/vector_infection.rds`. 

-   `Engorgement_Host_Infection_Status.qmd` tests the effect of host infection status (what we call the simple model) and other biological and experimental variables (what we call the complete model) on engorgement rates. It produces **Figures 3, 5, and S2A**. Model objects are saved in `output/engorgement_models/host_infection_simple.rds` and `output/engorgement_models/host_infection_complete.rds`. 

-   `Preliminary_Tests.R` performs tests that justify the exclusion of data from day 28 of the experiments, as well as interaction terms in the complete model. 

-   `Approach_Data_Manaus.R` uses data from a field study in Manaus, Brazil, measuring approach rates of mosquitoes on different occasions, with different collectors, and different levels of urbanization and times of day (Hernandez Acosta, 2023). It computes the coefficient of variation in what we define as optimum conditions. It produces **Figure S2B**. 

-   `Cytokines_Analysis.R` runs generalized linear mixed effect models, as well as generalized additive models, and performs a selection, to describe the possible associations between cytokine concentrations in control monkeys, and the short-term and long-term exposure to uninfected mosquito bites. It produces **Figures 6 and Figures S3-S8**.

### Data folder

Files ending with `data_engorgement` contain data to derive engorgement rates, for all arms of the experiments, either at day 0 (`Day0`) for infected mosquitoes, or other days for uninfected mosquitoes.

Files starting with `Hanley2024` correspond to data described in the paper by Hanley *et al.* 2024 (see reference below).

Files starting with `Temperatures`provide host body temperature from a transponder every 15 minutes.

`Cytokines_and_bites_data.csv` contains cytokine concentrations over time along with short-term and long-term exposure to mosquito bites.

`Control_monkeys_cumulative_bites.csv` contains various variables related to bite exposure over time (those used in the analysis as well as other cumulative metrics). It includes days 11, 14, and 21, as mosquito exposure took place on those days, even though cytokine concentrations were not measured.

`Approach_rate_Hernandez_Acosta.csv` contains approach data from Hernandez Acosta, 2023.

### Output folder

#### Figure subfolder

Contains figures from the manuscript and the supplementary information, numbered accordingly.

The `outlines` subfolder contains monkey and human images used in figures. The human image was created by Soremba and downloaded from [flaticon.com](https://www.flaticon.com/). The monkey images are licensed from Shutterstock.

#### Engorgement_models subfolder

Contains selected models to describe observed engorgement rates.

#### Cytokines_analysis

Contains reports from models applied to cytokines. These reports were used to produce **Table S12** (provided as Supplementary Information, not in the present repository).

Any questions about this code base can be directed at [helene.cecilia3\@gmail.com](mailto:helene.cecilia3@gmail.com)

### References 
Hernandez Acosta, E. *Some like it hot : How urban microclimate across a tropical city impacts the capacity of Aedes mosquitoes to transmit flaviviruses*. PhD thesis (New Mexico State University, **2023**).

Hanley KA, Cecilia H, Azar SR, Moehn B, Gass J, Oliveira da Silva NI, Yu W, Yun R, Althouse BM, Vasilakis N, Rossi SL. **2024** Trade-offs shaping transmission of sylvatic dengue and Zika viruses in native and novel monkey hosts. *bioRxiv* 10.1101/2023.06.30.547187 now accepted in *Nature Communications*
