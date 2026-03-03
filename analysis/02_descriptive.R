
###
### Structure of the `bl.dat` dataset
#
#' Variables:
#' - age_Dx: age at ALL diagnosis (numeric, years)
#' - sexM: sex
#'    - female = 0
#'    - male = 1
#' - SES: income quintile
#'    - If all below vars are 0, then the individual is part of SES1
#'    - SES2
#'    - SES3
#'    - SES4
#'    - SES5
#' - ruralY: rurality
#'    - 0 = urban residence
#'    - 1 = rural residence
#' - distanceshort: distance to treating pediatric centre
#'    - 0 = long, 1 = short
#' - dxyear: year of ALL diagnosis
#'    - If all below vars are 0, then dxyear is 2002
#'    - 2002
#'    - 2003
#'    - 2004
#'    - 2005
#'    - 2006
#'    - 2007
#'    - 2008
#'    - 2009
#'    - 2010
#'    - 2011
#'    - 2012
#' - lineageT: type of ALL
#'    - 0 = B-cell
#'    - 1 = T-cell
#' - protocolDFCI: protocol type
#'    - 0 = DFCI (Dana-Farber Cancer Institute)
#'    - 1 = COG (Children's Oncology Group)
#' - cytogenetics: genetic features
#'    - if all below vars are 0, then cytogenetics category is 9;22 (Ph+)
#'        catFavorable
#'        catUnfavorable
#'        catNeutral
#' - wbc_cat>=50: white blood cell (wbc) count at diagnosis
#'    - 0 = <50 
#'    - 1 = >=50
#' - leukemia_risk_groupstandard: NCI risk group
#'    - 0 = High (all others)
#'    - 1 = Standard (aged 1.0-9.99 at diagnosis, WBC of < 50)
#' - cns_status: CNS disease at diagnosis
#'    - If all below vars are 0, then CNS status is 1
#'    - CNS2
#'    - CNS3
#' - mrd: Minimal residual disease status. If all are 0, then mrd is 0.01%-0.1%
#'        mrdNegative
#'        mrd>0.1%-1%
#'        mrd>1%
#'        mrdNot measured
#' - dxyear: year of ALL diagnosis (numeric, years)
#'    - ref is 2002. So 0 = 2002, 1 = 2003, etc.
#' - age>=15: age at diagnosis category
#'    - 0 = aged 14 or lower 
#'    - 1 = 15 or higher
#' - AALL1731_elig: AALL1731 eligibility (Fig S3; DOI: 10.1056/NEJMoa2411680)
#'    - 0 = not eligible
#'    - 1 = eligible
#' - AALL1731_risk: trial risk stratification
#'    - "Fav" = favourable
#'    - "Avg" = average
#'    - "High" = high

# Set random seed and generate synthetic cohort --------------------------------
set.seed(SEED) # uses the SEED defined in the calling script

# Baseline covariate generation ------------------------------------------------
source(file = paste0(filepath, "/R/generate-synthetic-cohort.R"))

bl.dat$Age <- bl.dat$age_Dx

# Subset the cohort based on the analysis being run ----------------------------
# - `n_cat` is the scenario identifier used to subset cohorts for analyses

# if MRD, remove mrd not measured - the synthetic cohort has complete info
if (n_cat == "MRD" & sum(bl.dat$`mrdNot measured` == 1) != 0) {
  bl.dat <- bl.dat[-which(bl.dat$`mrdNot measured` == 1), ]
} 

## Blina (first line therapy)
# - subset the cohort to only those who would be eligible AALL1731 trial
if (n_cat == "blin_1stline") {
  bl.dat <- bl.dat[bl.dat$AALL1731_elig == 1 & bl.dat$AALL1731_risk != "Fav", ]
}

## NCI high risk cohort subset
# - subset to age >= 1 and < 18 years
if (n_cat == "NCI_high") {

  # filter to NCI high, and 1-17
  bl.dat <- bl.dat[bl.dat$NCI_high == 1 & bl.dat$age_Dx >= 1 & bl.dat$age_Dx < 18, ]

  # force all NCI high risk to not be standard risk (for policy model)
  bl.dat$leukemia_risk_groupstandard <- 0
}

# Overall Population
if (n_cat %in% c("overall", "blin_1stline", "NCI_high")) {
  b.dat.01 <- bl.dat[sample(nrow(bl.dat), n_i_i), ]
} else if (subgroup == TRUE || subgroup == "Bcell") { # Subgroup Analysis

  for (i in 1:length(v_cat)) {
    if (i == 1) {
      temp.dat <- data.frame("cat" = rep(0, nrow(bl.dat)))   #create temp dataframe equal to size of # of subgroups (Cat) as many rows as bl.dat, sets everything to 0 in that cat
    } else {
      temp.dat$cat[bl.dat[, v_cat[i]] == TRUE] <- i - 1      #set to 1,2, etc for the categories you want (replace the 0s)
    }
  }

  bl.dat$cat <- temp.dat$cat

  b.dat.01 <- c()
  for (i in 1:length(v_cat)) {
    bl.help <- bl.dat %>%
      mutate(id = rownames(.)) %>%
      filter(temp.dat$cat == i - 1) %>%  #check on no use of brackets here: implement bracket for [i-1]?
      sample_n(n_i_s, replace = FALSE) %>%  # sample n_i_s individuals with v_cat[i]
      ungroup()

    bl.help <- bl.help[, names(bl.help) != "id"]

    assign(paste0("bl.dat.", i), bl.help)
    b.dat.01 <- rbind(b.dat.01, get(paste0("bl.dat.", i)))   #stacking happens here (find the individuals who have the cat(s) you want)
  }

} else { # Counter-factual
  for (i in 1:length(v_cat)) {
    if (i == 1) {
      temp.dat <- data.frame("cat" = rep(0, nrow(bl.dat)))
    } else {
      temp.dat$cat[bl.dat[, v_cat[i]] == TRUE] <- i - 1
    }
  }

  b.dat.01 <- c()

  for (i in 1:length(v_cat)) {
    if (i == 1) {
      assign(
        paste0("bl.dat.", i), bl.dat %>%
          mutate(id = rownames(.)) %>%
          filter(temp.dat$cat == i - 1) %>%
          sample_n(n_i_s, replace = FALSE) %>%
          ungroup() %>%
          dplyr::select(-c(id)) %>%
          mutate(cat = i)
      )
    } else {
      temp.bl.dat <- bl.dat.1
      temp.bl.dat[, v_cat[i]] <- 1
      temp.bl.dat$cat <- i
      assign(paste0("bl.dat.", i), temp.bl.dat)
    }

    b.dat.01 <- rbind(b.dat.01, get(paste0("bl.dat.", i)))
  }
}

# Calculate other variables of interest for descriptives and cost model --------

# Variables used in model not in cohort gen
b.dat.01[, 'I(age_Dx^2)'] <- b.dat.01$age_Dx^2
b.dat.01[, "REL1.t"] <- NA_real_

b.dat.01 <- b.dat.01 %>%
  mutate(
    `age>=15` = case_when(
      age_Dx >= 15 ~ 1,
      TRUE ~ 0
    ),
    mrdPositive = case_when(
      `mrd0.1%-1%` == 1 ~ 1,
      `mrd>1%` == 1 ~ 1,
      mrdNegative == 0 & `mrdNot measured` == 0 & `mrd0.1%-1%` == 0 & `mrd>1%` == 0 ~ 1,
      TRUE ~ 0
    ),
    prev_REL1 = case_when(
      v_M_Init_nii %chin% c("REL1", "REL2", "BMT1_REL") ~ 1,
      TRUE ~ 0
    ),
    prev_BMT = case_when(
      v_M_Init_nii == "BMT1_REL" ~ 1,
      TRUE ~ 0
    )
  )

# Body surface area (BSA) calculation
# - use WHO data to get height and weight for BSA calculation
WHO_age_height <- read.csv("model-parameters/raw-data/WHOCan_hgt_wgt.csv") %>%
  mutate(
    sexM = ifelse(sex == "male", 1, 0)
  )

# Get height and weight for BSA calculation
temp_data <- b.dat.01 %>%
  transmute(age = round(age_Dx), sexM) %>%
  left_join(WHO_age_height, by = c("age", "sexM")) %>%
  dplyr::select(age, sex = sexM, height, weight)

# Calculate BSA using Mosteller formula
b.dat.01$bsa <- admiral::compute_bsa(
  height = temp_data$height,
  weight = temp_data$weight,
  method = "Mosteller"
)

rm(temp_data, WHO_age_height) # remove temporary data objects

# Baseline for costs
bl.dat.c <- b.dat.01
bl.dat.c$interval <- rep(step/30, n_i_i)
bl.dat.c$`I(interval^2)` <- rep((step/30)^2, n_i_i)
bl.dat.c$`(Intercept)` <- 1
bl.dat.c$dxyear <- bl.dat.c$dxyear + 2002

# Descriptive table of the cohort ----------------------------------------------

# Specify variables that are categorical for descriptive table
catvars <- c("sexM", "ruralY", "dxyear","wbc_cat>=50", "lineageT","leukemia_risk_groupstandard")

# Create descriptive table of the cohort
myvars  <- c("age_Dx", "cat", catvars)

if (all(is.na(v_cat)) || subgroup == "Bcell") {
  # Table for blina first line
  if (n_cat == "blina_first_risk") {
    myvars <- c(myvars, "AALL1731_risk")
    catvars <- c(catvars, "AALL1731_risk")
    tab2 <- CreateTableOne(vars = myvars, data = b.dat.01, factorVars = catvars, strata = "cat")
  } else {
    tab2 <- CreateTableOne(vars = myvars, data = b.dat.01, factorVars = catvars)
  }
} else {
  tab2 <- CreateTableOne(vars = myvars, data = b.dat.01, factorVars = catvars, strata = "cat")
}

# Save descriptive table and cohort info to the output folder
write.csv(x = print(tab2), file = paste0(filepath, path.out, "tableone.csv"))

# Save BSA data if blina or bsc analysis
# - save the BSA information since it is used to calculate chemo costs later
if (is_TRUE(blin_2ndline) || is_TRUE(bsc_2ndline) || is_TRUE(bsc_1stline) || is_TRUE(blin_1stline)) {
  save_df <- b.dat.01 |> dplyr::select(age = age_Dx, sexM, bsa)
  write.csv(x = save_df, file = paste0(filepath, path.out, "cohort_info.csv"), row.names = FALSE)
  rm(save_df)
}
