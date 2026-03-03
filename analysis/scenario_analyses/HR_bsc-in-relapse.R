#' -----------------------------------------------------------------------------
#' Title:
#'   Microsimulation Analysis for Blinatumomab First-Line in NCI High-Risk
#'   B-cell ALL, BSC for Relapsed Disease
#' 
#' Population:
#'   High-risk (NCI) pediatric B-cell ALL
#' 
#' Description:
#'   Two-arm microsimulation comparing chemo vs chemo+blina in first line NCI
#'   high-risk B-cell ALL. For any relapsed disease, all patients receive
#'   best supportive care.
#'
#' Sourced Scripts:
#'   analysis/00_setup.R
#'   analysis/01_input_params.R
#'   analysis/02_descriptive.R
#'   analysis/03_model_run.R
#'   analysis/04_save_results.R
#'   R/validate_blina_1stline.R
#'   R/economic_evaluation.R
#' 
#' Outputs:
#'   - Model results saved under model-out/blin_1stline_nci-high-bsc/
#'   - Summary statistics and cost-effectiveness results
#'
#' -----------------------------------------------------------------------------

# Load the setup script to initialize the environment
source("analysis/00_setup.R")

# Load specific inputs for HR's
bHr.norm <- exp(rnorm(n_sim, -0.356, 0.2)) # from Brown et al
b1_HR_PRE_REL1 <- read.csv(paste0(getwd(), "/blina-trial-calibration/data/Gupta_bootHR_PFS_PSA.csv"))

# Best supportive care arm: chemo alone ----------------------------------------
set.seed(SEED)
n_i <- 2000000
n_i_i <- n_i_s_all
n_i_s <- n_i_i

n_cat <- "NCI_high"
v_cat <- NA
subgroup <- ""

path.out <- "/model-out/blin_1stline_nci-high-bsc/bsc/"
path.out2 <- "/model-out/blin_1stline_nci-high-bsc/"

filepath <- getwd()

v_M_Init_nii <- rep("PRE", n_i_i)

source(paste0(filepath, "/analysis/01_input_params.R"))

# First line treatment flags
bsc_1stline <- TRUE
blin_1stline <- blin_1stline_CEres <- FALSE

# Relapse specific flags
blin_1stline_relapse_bsc <- TRUE
blin_1stline_relapse_blin <- FALSE

source(paste0(filepath, "/analysis/02_descriptive.R"))
source(paste0(filepath, "/analysis/03_model_run.R"), local = T)
source(paste0(filepath, "/analysis/04_save_results.R"))

# Intervention arm: standard of care chemo + blina -----------------------------

set.seed(SEED)
n_i <- 2000000
n_i_i <- n_i_s_all
n_i_s <- n_i_i

n_cat <- "NCI_high"
v_cat <- NA
subgroup <- ""

path.out <- "/model-out/blin_1stline_nci-high-bsc/blin/"
path.out2 <- "/model-out/blin_1stline_nci-high-bsc/"

filepath <- getwd()

v_M_Init_nii <- rep("PRE", n_i_i)

source(paste0(filepath, "/analysis/01_input_params.R"))

# First line treatment flags
bsc_1stline <- FALSE
blin_1stline <- blin_1stline_CEres <- TRUE

# Relapse specific flags
blin_1stline_relapse_bsc <- TRUE
blin_1stline_relapse_blin <- FALSE

source(paste0(filepath, "/analysis/02_descriptive.R"))
source(paste0(filepath, "/analysis/03_model_run.R"), local = T)
source(paste0(filepath, "/analysis/04_save_results.R"))

# Health outcomes and economic evaluation --------------------------------------

source(paste0(filepath, "/R/validate_blina_1stline.R"))
source(paste0(filepath, "/R/economic_evaluation.R"))
