#' -----------------------------------------------------------------------------
#' Title:
#'   Setup for Microsimulation Analysis
#'
#' Description:
#'   This script initializes the environment for conducting a microsimulation
#'   analysis. It handles package loading, workspace cleanup, and
#'   cohort‐specification settings for downstream probabilistic sensitivity
#'   analyses.
#'
#' Sections:
#'   01 Load Libraries        - Installs (if needed) and loads CRAN packages
#'   02 Clean Environment     - Clears all objects
#'   03 Cohort Specifications - Sets random seed, PSA iterations, and
#'                              sample sizes
#'
#' Key Parameters:
#'   SEED      : Integer seed for reproducibility of random sampling (2022)
#'   n_sim     : Number of probabilistic sensitivity analysis iterations
#'   n_i       : Total number of individuals generated initially
#'   n_i_s_all : Number of individuals sampled per subpopulation
#'
#' Usage:
#'   Source this file at the start of analysis workflows to ensure a consistent
#'   computing environment and cohort setup.
#'
#' -----------------------------------------------------------------------------

## 01 Load libraries ------------------------------------------------

# Install packages from CRAN (if required) and load packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(
  "survival", "flexsurv", "DAAG", "msm", "demography",
  "gems", "cmprsk", "MASS", "openxlsx", "reshape", "etm",
  "rlang", "mgcv", "mvtnorm", "corpcor", "gridExtra", "data.table", "plyr",
  "survminer", "foreach", "doParallel", "splines", "gmodels", "flexsurvcure",
  "matrixStats", "muhaz", "fst", "doRNG",
  "tableone", "zoo", "dedupewider", "zscorer", "admiral", "furrr", "tidyverse"
)

# Load custom packages
p_load_gh("DARTH-git/darthtools")
p_load_gh("DARTH-git/dampack")

# Load custom functions
for (file in list.files("R/functions", pattern = "\\.R$", full.names = TRUE)) {
  source(file)
}

## 02 Cohort specifications ------------------------------------------------

SEED <- 2022 # set the seed for the probabilistic analysis
set.seed(SEED)

# Number of Probabilistic Sensitivity Analyses to capture parameter uncertainty
n_sim <- 500
# Number of individuals to generate from which a subset will be used in analysis
n_i <- 200000
# Number of individuals to be randomly sampled from n_i by (sub)population.
n_i_s_all <- 30000
