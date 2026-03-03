################################################################################
# Microsimulation Input Setup
#
# This script defines and samples the core input parameters used by the
# microsimulation model, including state definitions, time settings, treatment
# assumptions, adverse event inputs, utilities, costs, and output paths.
################################################################################

set.seed(SEED) # ensure seed is set to ensure reproducibility each time this script is run

# Specify paths to model inputs ------------------------------------------------
# - paths where the model-parameters derived from real-world data are stored

event.path <- paste0(getwd(), "/model-parameters/events/") # where event inputs are stored
costs.path <- paste0(getwd(), "/model-parameters/costs/") # where cost inputs are stored

# Model parameters -------------------------------------------------------------

v_n      <- c("PRE","REL1","REL2","BMT1","BMT1_REL","DEAD")     # model states (EVENTS)
n_s      <- length(v_n)                                         # number of states
n_t      <- 100                                                  # number of years to run
c_l      <- 30/365                                              # cycle length (30 days) 
step     <- 30                                                  # Time step relative to days 
times    <- seq(from = 0, to = n_t * 365, by = step)            # sequence of times to be considered in the model
cycles   <- length(times)                                       # number of cycles

################################################################################
# Blina first-line specific inputs
################################################################################

# Time (days) to start blina to align with the trial randomization time
START_BLIN <- 30
PRICE_BLIN <- 2978.26 # price per vial of blina (CADTH report for Blincyto; https://www.ncbi.nlm.nih.gov/books/NBK616209/)

# Create adverse events (grade >=3) and uncertainty parameters for inputs -----
# Notes:
#    - combining the SR-average and SR_high groups
#    - Table 2 of the NEJM paper (https://doi.org/10.1056/NEJMoa2411680) for grade 3+ AEs
#    - Table S5 in appendix for grade 4 AEs

# Define the NEJM trial demoninators in the blina + chem and chemo groups
denom_AALL1731 <- c(
    (351 + 273), # blina + chemo arm
    (376 + 277) # chemo arm
)

# Cytokine Release Syndrome
events_crs <- c(
    (1 + 1), # blina + chemo arm
    (0 + 0) # chemo arm
)

p_CRS <- sim_risk_diff(events_crs, denom_AALL1731, n_sim)

# Febrile Neutropenia (fever) ---------------------------------
events_fever <- c(
    (165 + 156), # blina + chemo arm
    (149 + 140) # chemo arm
)
events_fever_gr4 <- c(
    (4 + 1), # blina + chemo arm
    (1 + 0) # chemo arm
)

p_FN <- sim_risk_diff(events_fever, denom_AALL1731, n_sim)

dur_FN <- sim_weighted_ae_day(
    n_g3p = events_fever,
    dur_g3 = 2, # duration of hospitalization if FN (Source: Teuffel 2011)
    n_g4 = events_fever_gr4,
    dur_g4 = 6, # grade 4 fevers requiring hospitalization (Source: Teuffel 2011)
    n_sim = n_sim
) # create duration of weighted averages for febrile neutropenia days

dis_FN <- 0.15 # disutility for febrile neutropenia (Source: LLoyd 2008)

# Sepsis ---------------------------------------------
events_sepsis <- c(
    (52 + 57), # blina + chemo arm
    (19 + 47) # chemo arm
)
events_sepsis_gr4 <- c(
    (3 + 8), # blina + chemo arm
    (0 + 5) # chemo arm
)

p_SEP <- sim_risk_diff(events_sepsis, denom_AALL1731, n_sim)

dur_SEP <- sim_weighted_ae_day(
    n_g3p = events_sepsis,
    dur_g3 = 2, # duration of hosp days with grade <4 sepsis (Assumption)
    n_g4 = events_sepsis_gr4,
    dur_g4 = 18, # duration of hosp days with grade 4 sepsis (Basu et al 2006 JCO)
    n_sim = n_sim
)

dis_SEP <- 0.218 # disutility for sepsis (Source: Stein 2018)

# IVIG Scenario analysis -------------------------------------------------------
# - assume IVIG is dosed at 0.7g/kg per infusion
# - assume each patient receives 2 cycles of IVIG to align with blina cycles
# - assume IVIG reduces the risk of FN by 28% as seen in RCT (Thus et al, 2025 Haemtologica)

# IVIG Price per gram
PRICE_IVIG <- 23.88 # CAD per gram of IVIG (in 2018 dollars)
IVIG_REDUC <- (1 - 0.28) # reduction in risk of FN with IVIG

# ------------------------------------------------------------------------------

# Flags for the analysis scenarios ---------------------------------------------

blin_2ndline <- bsc_2ndline <- FALSE # blina as second line treatment
bsc_1stline <- blin_1stline <- FALSE # blina as first line treatment
blin_1stline_CEres <- blin_2ndline_CEres <- FALSE # blina cost-effectiveness

# Discounting  -----------------------------------------------------------------

d_r <- 0.015 # discount rate
v_dwe <- v_dwc <- 1 / ((1 + d_r)^(times / step / 12)) # discount weight

# Health utilities  ------------------------------------------------------------

# from Furlong et al. Table 3, HUI3 2012 by treatment phase (Mean, SD and N):
# - add uncertainty for each PSA iteration set by `n_sim`
sim_u__PRE <-((71 * rnorm(n_sim, 0.87, 0.21 / sqrt(121)) +     # Continuation
               30 * rnorm(n_sim, 0.79, 0.25 / sqrt(124)) +     # Intensification
                3 * rnorm(n_sim, 0.75, 0.23 / sqrt(137)) +     # CNS
                4 * rnorm(n_sim, 0.67, 0.32 / sqrt(196)))/104) # Induction

# Health utility values for each health state
u_PRE    <- sim_u__PRE                         # Post-diagnosis (no relapse hx) short-term (if <2 years in state) [Furlong et al., Table 3 HUI3]
u_PRE_LT <- rnorm(n_sim, 0.90, 0.18/sqrt(171)) # Post-diagnosis (no relapse hx) long-term  (if >2 years in state) [Furlong et al., Table 3 HUI3]
u_REL    <- rnorm(n_sim, 0.75, 0.01)           # Post-relapse short-term (if <5 years in state) [Kelly et al. Table 1, EQ-5D]
u_REL_LT <- rnorm(n_sim, 0.87, 0.02)           # Post-relapse long-term (if >5 years in state) [Lawitschka et al, Table 3, PedsQL]
u_BMT    <- rnorm(n_sim, 0.73, 0.01)           # Post-BMT short-term (if <5 years in state) [Kurosawa et al., Table 5, EuroQOL5D]
u_BMT_LT <- rnorm(n_sim, 0.87, 0.02)           # Post-BMT long-term (if >5 years in state)  [Lawitschka et al., Table 3, PedsQL]
u_DEAD   <- 0

# Cost model parameters --------------------------------------------------------

costs.trans  <- "log"    # transformation used on dependent variable in cost model
costs.spline <- "3"      # df of splines


# Specify output paths ---------------------------------------------------------
# - calls the `path.out` and `path.out2` variables to create a directory to
#   store model outputs

filepath <- getwd()
output.path <- paste0(filepath, path.out)
dir.create(output.path, recursive = TRUE)

if (is_TRUE("blin_2ndline") || is_TRUE("blin_1stline")) {
    output.path2 <- paste0(filepath, path.out2)
}

# Source functions and MSM outputs ---------------------------------------------
# - source the functions that are used in the model to setup the microsim
source(file = paste0(filepath, "/R/microsim_functions_auto.R"))
