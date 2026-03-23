# Cost-effectiveness of Frontline Blinatumomab for Pediatric B-cell Acute Lymphoblastic Leukemia

## Overview

Repository for the manuscript "**Cost-effectiveness of frontline blinatumomab in children with standard- and high-risk B-cell acute lymphoblastic leukemia in Ontario, Canada**". The repository describes a cost-effectiveness analysis (CEA) of using frontline blinatumomab in children diagnosed with NCI standard- and high-risk B-cell acute lymphoblastic leukemia (ALL). This analysis adopts a previously published "ALL Policy" microsimulation model ([Development of a policy model for pediatric acute lymphoblastic leukemia to facilitate economic evaluation](https://doi.org/10.1093/jnci/djaf024)) by Pechlivanoglou _et al_ (2025) that was constructed using real-world evidence from children diagnosed and treated from 2002 to 2017 in the universal healthcare system of Ontario, Canada. With the ALL-Policy model as a basis, we have augmented it to incorporate an economic evaluation of blinatumomab as a frontline therapy using treatment efficacy informed by the [AALL1731 trial](https://clinicaltrials.gov/study/NCT03914625) results ([Gupta _et al_ NEJM 2025](https://doi.org/10.1056/NEJMoa2411680)).

The objectives of this manuscript are to assess the cost-effectiveness of frontline blinatumomab for:

- i) standard-risk B-ALL (defined as age 1-9 at diagnosis with white blood cell count < 50,000/&micro;L)
- ii) high-risk B-ALL (white blood cell count &ge; 50,000 /&micro;L or age 10-17 at diagnosis)

## Usage

This project contains several wrapper scripts to conduct the CEA analyses. The main scripts are `analysis/SR_analysis.R` or `analysis/HR_analysis.R` as entry points, depending on the target risk group. These scripts execute the microsimulation pipeline and save the results in `model-out/`, in the respective target folders based on the population being studied. The `economic_evaluation.R` script is then run once the microsimulation has completed and saves aggregated output into the `report/` folder which is used for manuscript production.

## Project Structure

```plaintext
.
├── README.md
├── analysis/
│   ├── scenario_analyses/          # Scenario analysis scripts
│   │   ├── HR_bsc-in-relapse.R     #   Chemotherapy only for all relapses in high-risk B-ALL
│   │   └── SR_bsc-in-relapse.R     #   Chemotherapy only for all relapses in standard risk B-ALL
│   ├── 00_setup.R                  # Initializes libraries, globals, and helper setup
│   ├── 01_input_params.R           # Sets up model parameters and sources functions
│   ├── 02_descriptive.R            # Initiates synthetic data generation
│   ├── 03_model_run.R              # Main microsimulation run script
│   ├── 04_save_results.R           # Saves model output and calculates costs/effects
│   ├── SR_analysis.R               # Standard-risk analysis wrapper
│   └── HR_analysis.R               # High-risk analysis wrapper
├── blina-trial-calibration/        # Trial calibration code and data for blinatumomab first-line use
│   ├── analysis/
│   ├── data/
│   ├── data-raw/
│   ├── figs/
│   └── R/
├── manuscript/                     # Manuscript-related outputs
├── model-out/                      # Microsimulation model outputs (costs, effects, state histories, etc.)
│   ├── blin_1stline_relapse_blin/  #   Main results SR B-ALL: output from `SR_analysis.R`
│   ├── blin_1stline_nci-high-blin/ #   Main results HR B-ALL: output from `HR_analysis.R`
│   ├── blin_1stline_nci-high-bsc/  #   Scenario analysis (output from `HR_bsc-in-relapse.R`)
│   └── blin_1stline_relapse_bsc/   #   Scenario analysis (output from `HR_bsc-in-relapse.R`)
├── model-parameters/               # Final model inputs estimated from data
│   ├── cohort/                     # Cohort generation parameters and distributions
│   ├── costs/                      # Cost regression estimates (betas)
│   ├── events/                     # Event/transition model regression estimates (betas)
│   └── raw-data/                   # Underlying data used to estimate parameters (e.g., life tables, growth charts)
├── R/
│   ├── functions/                  # Functions needed for the project
│   ├── manuscript/                 # Scripts used for manuscript table generation
│   ├── economic_evaluation.R       # Cost-effectiveness evaluation and summary outputs
│   ├── generate-synthetic-cohort.R # Creates synthetic cohort for microsimulation
│   ├── microsim_functions_auto.R   # Sets up microsimulation environment (loads models, costs, etc.)
│   └── validate_blina_1stline.R    # Validation plots and combined survival outputs
└── report/                         # Processed outputs for reporting
```

## Requirements

```r
install.packages(c(
  "admiral", "broom", "cmprsk", "corpcor", "data.table", "DAAG", "dedupewider", "demography",
  "dplyr", "doParallel", "doRNG", "etm", "flexsurv", "flexsurvcure", "foreach", "fst", "furrr", "gems",
  "ggplot2", "gmodels", "gridExtra", "MASS", "matrixStats", "mgcv", "mice", "mstate", "msm",
  "muhaz", "mvtnorm", "nnet", "openxlsx", "pacman", "plyr", "purrr", "rlang", "reshape",
  "splines", "survminer", "survival", "tableone", "tidyr", "tidyverse", "zoo", "zscorer"
))

# darthtools
pacman::p_load_gh("DARTH-git/dampack")
pacman::p_load_gh("DARTH-git/darthtools")
```
