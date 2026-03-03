#' Calculate Chemotherapy and Blinatumomab Treatment Costs
#'
#' @description
#'   This function estimates the total costs associated with chemotherapy (BSC)
#'   or Blinatumomab (Blina) treatment for pediatric patients, based on body
#'   surface area (BSA) and probability of adverse events 
#'   (eg. febrile neutropenia, sepsis, cytokine release syndrome). It
#'   includes drug costs, nurse time, equipment, and hospitalization.
#'
#' @note
#'   Costs are to be provided in 2018 CADs since that is what the model was
#'   built using. Any inflation adjustments are made in post-processing.
#'
#' @param bsa Numeric or vector. Patient body surface area in m^2.
#'   Default is 1.09 (average height and weight for a 4 year old male/female).
#' @param blina_price Numeric. Price per vial of Blinatumomab.
#'   Default is 2978.26 (Source: CADTH report for Blincyto; https://www.ncbi.nlm.nih.gov/books/NBK616209/)
#' @param p_CRS Numeric. Probability of experiencing Cytokine Release
#'   Syndrome. Default is 0.42 (Source: Brown et al. trial for relapsed disease)
#' @param p_FN Numeric. Probability of experiencing febrile neutropenia.
#'   Default is 0 (since this is not a cost included in the base case, but can be added in sensitivity analyses).
#' @param dur_FN Numeric. Duration of febrile neutropenia in days.
#'   Default is 0 (since this is not a cost included in the base case, but can be added in sensitivity analyses).
#' @param p_SEP Numeric. Probability of experiencing sepsis.
#'   Default is 0 (since this is not a cost included in the base case, but can be added in sensitivity analyses).
#' @param dur_SEP Numeric. Duration of sepsis in days.
#'   Default is 0 (since this is not a cost included in the base case, but can be added in sensitivity analyses).
#' @param trt_output Character. Specify which treatment cost to return:
#'   - "bsc" for best supportive care chemotherapy
#'   - "blin" for Blinatumomab.
#'   Default is "bsc"
#'
#' @return Numeric.
#'   Total cost for the specified treatment over the modeled period. Will return
#'   a vector of costs per individual is a vector of `bsa` is provided, or a
#'   single numeric value if `bsa` is a single number.
#'
#' @details
#' - Costs are calculated per cycle and summed over the treatment period.
#' - Includes microcosting for chemotherapy drugs and administration.
#' - Blinatumomab costs include drug, nurse time, IV pump, and hospitalization
#'   (including CRS).
#' - Hospitalization costs are based on average pediatric stays.
#'
#' @examples
#' chemo_costs(bsa = 1.2, p_CRS = 0.42, trt_output = "blin")
#' chemo_costs(bsa = c(0.89, 1.2), p_CRS = 0.000035, trt_output = "bsc")
#'
chemo_costs <- function(bsa = 1.09, blina_price = 2978.26, p_CRS = 0.42, p_FN = 0, dur_FN = 0, p_SEP = 0, dur_SEP = 0, trt_output = "bsc") {
  
  # Wrap this in an `mapply` to vectorized operations
  mapply(function(bsa_i, p_CRS_i, p_FN_i, dur_FN_i, p_SEP_i, dur_SEP_i, blina_price_i) {

    # General --------------------------------------------------------
    model_cycle  <- 7 # model's cycle length (days)
    c_nurse_hour <- 39.16 * 1.25 # salary of pediatric nurse * 25% (Source: Employment and Social Development Canada Job Bank)

    # Blina specific  --------------------------------------------------------

    c_Blina           <- blina_price # cost of blina per vial

    c_pump            <- 5000    # cost of IV pump (Source: Assumption)
    le_pump           <- 10      # life expectancy of pump (Source: Assumption)
    Blina_cycle       <- 28      # length of one Blina cycle (days)
    Blina_unit        <- 38.5    # mcg per vial
    Blina_day_dose    <- 15      # 15 mcg/m2/day (max 28 mcg/day)
    hosp1_days        <- 2       # average days of in-patient admission with each Blina cycle (Source: Assumption)
    nurse_blina_mins  <- 30      # 30 mins of nurse time assumed for each blina administration (Source: Assumption)
    visits_blina_days <- 4       # assuming a 96-hour bag: bag refill required every 4 days. (28 day cycle / 4 = 7 outpatient visits)
    visit_week        <- model_cycle / visits_blina_days # 7 visits in a 28 day cycle is 1.75 on average, weekly
    c_pump_week       <- c_pump / (le_pump * 52)
    c_nurse_week      <- c_nurse_hour * nurse_blina_mins / 60 * visit_week

    # Calculate blina treatment costs
    c_Blina_week <- model_cycle * Blina_day_dose / Blina_unit * c_Blina * bsa_i
    c_trt <- sum(c_pump_week + c_nurse_week + c_Blina_week)

    # Hospitalization associated costs (additional) ----------------------------
    
    c_hosp_total <- 24931 # average cost of hospital stay for child aged 1-7 (Source: CIHI Patient Cost Estimator. CMG: 625 - Acute Leukemia except Myeloid, Jurisdiction: Canada)
    hosp_days    <- 11.7 # average length of stay for child aged 1-7 (Source: CIHI Patient Cost Estimator. CMG: 625 - Acute Leukemia except Myeloid, Jurisdiction: Canada)
    c_hosp_day   <- c_hosp_total / hosp_days # cost per day of hospitalization

    # Cytokine Release Syndrome (CRS)
    hosp_CRS_days <- 9 # duration of hospitalization if CRS (average of 2 + additional 7 = 9) (Source: Assumption)

    # Febrile Neutropenia (FN)
    hosp_FN_days <- 2 # duration of hospitalization if FN (Source: Teuffel 2011)
    c_FN_perEpisode <- 5579 * 1.17 # avg cost per person per FN episode inflated from 2009 in the paper to 2018 when all other model costs are in (Source: Teuffel 2011)

    # Sepsis
    hosp_SEP_days <- 18 # duration of hospitalization if Sepsis for children with cancer and grade 4 sepsis (Basu et al 2006 JCO)

    # Determine additional costs associated with comorbidities -----------------
    c_hosp <- c_hosp_day * hosp1_days # baseline hospitalization costs for Blina

    c_hosp <- c_hosp + (c_hosp_day * p_CRS_i * hosp_CRS_days) # CRS costs
    #c_hosp <- c_hosp + (c_hosp_day * p_FN_i * dur_FN_i) # FN costs (using cihi estimate per day)
    c_hosp <- c_hosp + (p_FN_i * c_FN_perEpisode) # FN costs via the Teuffel 2011 estimate
    c_hosp <- c_hosp + (c_hosp_day * p_SEP_i * dur_SEP_i) # Sepsis costs

    # Add additional hospitalization costs to blina costs ----------------------
    c_Blina_total <- c(rep(c_trt, Blina_cycle / model_cycle), 0, rep(c_trt, Blina_cycle / model_cycle))
    c_Blina_total[1] <- c_Blina_total[1] + c_hosp

    # Total costs for Blina treatment
    trtcost_blin <- sum(c_Blina_total)

    # Chemotherapy Microcosting  -----------------------------------------------
    # - These costs concern chemotherapy administration for relapsed disease
    # - Assuming a treatment protocol similar to the SOC arm in Brown et al.
    #
    # Sources for drug costs:
    # - ODBF = Ontario Drug Benefit Formulary/Comparative Drug Index (Source: https://www.formulary.health.gov.on.ca/formulary/)
    # - HQO = Health Quality Ontario, report on ALL MRD (Source: PMID27099644)
    #
    
    Dex_Unit   <- 5
    c_Vin      <- 30.60    # per mg      (Source: ODBF search Vincristine)
    c_Peg      <- 5106.98  # 2500 units  (Source: HQO, Table A5.8, Pegasparaginase lower range)
    c_Aspa     <- 6500     # 25000 units (Source: ?, Asparaginase Erwinia)
    c_Cyt      <- 9.48     # 100 mg      (Source: HQO, Table A5.8, Cytarabine)
    Cyt_unit   <- 100
    c_Eto      <- 0.75     # per mg      (Etoposide, Source: CADTH report for Atezolizumab)
    c_Leuc     <- 7        # per 10mg/ml (Source: ODBF search Leucovorin)
    c_MTX      <- 8        # 25mg        (Source: ODBF search Methotrexate)
    c_Cycl     <- 0.4740   #             (Source: ODBF search Cyclophosphamide)
    Cycl_unit  <- 50
    MTX_unit   <- 25

    bsc_Dex_dose  <- 3      # dose in micrograms/m2  dose twice daily PO or IV  days 1 - 5
    MTX_dose      <- 1000   # mg
    Vin_dose      <- 1.5    # dose in mg/m2 (max 2mg) IV  day 1
    bsc_meth_dose <- 1000   # mg/m2 dose IV over 36 hours day 8
    Leuc_dose     <- 15     # mg/m2 dose every 6 hours IV or PO  days 10, 11
    Peg_dose      <- 2500   # IU/m2 dose IV  days 9 or 10
    Cycl_dose     <- 440    # mg /m2 dose IV days 15 - 19
    Eto_dose      <- 100    # mg/m2 dose IV days  15 - 19

    Cyt_dose      <- 1000   # mg/m2 dose every 12 hours IV over 3 hours days 1, 2, 8, 9
    Aspa_dose     <- 25000  # IU/m2 dose IM or IV days 2, 4, 9, 11, 23

    ## Working on the chemotherapy treatment

    # block 2
    c_Vin_week2   <- Vin_dose * bsa_i
    c_MTX_week2   <- c_MTX * MTX_dose / MTX_unit * bsa_i
    c_Peg_Week2   <- c_Peg * bsa_i
    c_Cycl_week2  <- c_Cycl * Cycl_dose / Cycl_unit * bsa_i * 4 # 4 days of therapy
    c_Eto_week2   <- c_Eto * Eto_dose * 4 * bsa_i               # 4 days of therapy

    c_block2      <- sum(c_Vin_week2, c_MTX_week2, c_Peg_Week2, c_Cycl_week2, c_Eto_week2)
    c_block2      <- c_block2 * model_cycle / Blina_cycle

    # block 3
    c_Vin_week3 <- Vin_dose * bsa_i
    c_MTX_week3 <- c_MTX * MTX_dose / MTX_unit * bsa_i
    c_Asp_Week3 <- c_Aspa * bsa_i * 5
    c_Cyt_week3 <- c_Cyt * Cyt_dose / Cyt_unit * bsa_i * 4 * 2 # 4 days  of therapy 2x daily
    c_Eto_week3 <- c_Eto * Eto_dose * 4 * bsa_i # 4 days of therapy

    c_block3 <- sum(c_Vin_week3, c_MTX_week3, c_Asp_Week3, c_Cyt_week3)
    c_block3 <- c_block3 * model_cycle / Blina_cycle

    # Put these costs into a vector representing weekly costs
    c_BSC_total <- c(
      rep(c_block2, Blina_cycle / model_cycle),
      rep(c_block3, Blina_cycle / model_cycle), 0
    )

    # Total costs for BSC treatment
    trtcost_bsc <- sum(c_BSC_total)

    if (trt_output == "bsc") {
      return(trtcost_bsc)
    } else if (trt_output == "blin") {
      return(trtcost_blin)
    } else {
      stop("Invalid treatment output specified. Use 'bsc' or 'blin'.")
    }
  }, bsa, p_CRS, p_FN, dur_FN, p_SEP, dur_SEP, blina_price)
}
