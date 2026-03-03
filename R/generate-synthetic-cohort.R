################################################################################
# Generate Synthetic Patient Cohort
#
# Purpose:
#   Generates a synthetic baseline cohort of `n_i` patients for use in the
#   microsimulation model. Patient characteristics are sampled sequentially
#   using regression models estimated from real-world ALL registry data from
#   Ontario, Canada (see Pechlivanoglou et al. 2025 JNCI).
#
# What it does:
#   1. Reads a dummy dataset to extract the variable names and factor levels
#      for all baseline covariates.
#   2. Samples age at diagnosis from a Gamma distribution fit to the cohort.
#   3. Sequentially samples each covariate (sex, SES, rurality, distance,
#      diagnosis year, lineage, protocol, cytogenetics, WBC category, risk
#      group, CNS status, and MRD status) conditioning on previously sampled 
#      covariates.
#   4. Create deterministic variables including:
#        - Numeric diagnosis year
#        - Age >= 15 flag
#        - MRD-positive flag
#        - NCI high-risk classification (B-cell, age, or WBC criteria)
#        - AALL1731 trial eligibility
#        - AALL1731 risk stratification (SR-Average, SR-High, SR-Favourable)
#
# Inputs:
#   - n_i                  : number of patients to simulate (set upstream)
#   - filepath             : root path to model-parameters/ directory
#   - model-parameters/cohort/dummy.dat.1000.csv    : reference dataset for levels
#   - model-parameters/cohort/model_<var>.csv       : regression coefficients
#   - model-parameters/cohort/vcov_<var>.csv        : variance-covariance matrices
#
# Output:
#   - bl.dat : a data frame of simulated patients (`n_i` rows) with all sampled 
#              and derived baseline covariates
################################################################################

# Specify all variables
vars <- c(
  "age_Dx", "sex", "SES",
  "rural", "distance", "dxyear",
  "lineage", "protocol", "cytogenetics_cat", "wbc_cat",
  "leukemia_risk_group", "cns_status", "mrd"
)

# Read in the dummy dataset to get variable names and levels for sampling
dummy <- read.csv(paste0(filepath, "/model-parameters/cohort/dummy.dat.1000.csv"))

# Remove the first column (age_Dx) and grab the variable names and levels
dummy <- dummy[, vars[-1]]
varnames <- list()
for (i in 1:ncol(dummy)) {
  varnames[[i]] <- levels(as.factor(dummy[, i]))
  varnames[[i]] <- varnames[[i]][varnames[[i]] != "Unknown"]
  names(varnames)[i] <- colnames(dummy)[i]
}

# Overwrite the variable names for `mrd`
varnames[["mrd"]] <- c("<0.1%", "0.1%-1%", ">1%", "Negative", "Not measured")

# Average age and SD of age in the ALL cohort
age.g <- dampack::gamma_params(5.88, 4.52)

# Initialize the data frame with the intercept and age, so that synthetic data
# can be generated
bl.dat <- data.frame(
  int = rep(1, n_i), # intercept
  age_Dx = rgamma(n_i, shape = age.g$shape, scale = age.g$scale) # age
)

# Simulate patients sequentially -----------------------------------------------

# Loop through the columns (besides age which was previously defined)
# - j indexes the variable to sample, and starts at position 3
# - for each variable to simulate, read in its corresponding regression coeff
#   and variance-covariance matrix, then sample it

jj <- 3 # index for bl.dat column (starting at 3 since age is already defined)

for (i in 1:length(vars[-1])) {
  # BINARY VARAIBLES
  vcov.m <- read.csv(paste0(filepath, "/model-parameters/cohort/vcov_", vars[i + 1], ".csv"), row.names = 1)
  est.m <- read.csv(paste0(filepath, "/model-parameters/cohort/model_", vars[i + 1], ".csv"), row.names = 1)

  # if one column then GLM with logit link
  if (ncol(est.m) == 1) {
    odds1 <- exp(as.matrix(bl.dat[, 1:nrow(est.m)]) %*% as.matrix(est.m))
    prob1 <- odds1 / (1 + odds1)
    dat.bin <- rbinom(n_i, 1, prob = prob1)
    bl.dat[[jj]] <- dat.bin
    colnames(bl.dat)[jj] <- paste0(vars[i + 1], varnames[[i]][2])
    jj <- jj + 1
    print(jj)
  } else {
    odds2 <- matrix(NA, nrow = n_i, ncol = nrow(est.m) + 1)
    odds2[, 1] <- 1

    for (j in 1:nrow(est.m)) {
      odds2[, j + 1] <- exp(as.matrix(bl.dat[, 1:length(est.m[j, ])]) %*% as.matrix(t(est.m[j, ])))
    }

    probs.m2 <- odds2 / rowSums(odds2)
    colnames(probs.m2) <- varnames[[i]]
    data.bin <- darthtools::samplev(probs.m2)

    for (l in 1:(length(varnames[[i]]) - 1)) {
      bl.dat[[jj]] <- (data.bin == varnames[[i]][l + 1]) * 1
      colnames(bl.dat)[jj] <- paste0(vars[i + 1], varnames[[i]][l + 1])
      jj <- jj + 1
      print(jj)
    }
  }
}

# Determininstic variable creation ---------------------------------------------
bl.dat <- bl.dat %>%
  mutate(
    dxyear = case_when(
      dxyear2003 == 1 ~ 1,
      dxyear2004 == 1 ~ 2,
      dxyear2005 == 1 ~ 3,
      dxyear2006 == 1 ~ 4,
      dxyear2007 == 1 ~ 5,
      dxyear2008 == 1 ~ 6,
      dxyear2009 == 1 ~ 7,
      dxyear2010 == 1 ~ 8,
      dxyear2011 == 1 ~ 9,
      dxyear2012 == 1 ~ 10,
      TRUE ~ 0
    ),
    `age>=15` = case_when(
      age_Dx >= 15 ~ 1,
      TRUE ~ 0
    ),
    mrdPositive = case_when(
      `mrd0.1%-1%` == 1 ~ 1,
      `mrd>1%` == 1 ~ 1,
      `mrd0.1%-1%` == 0 & `mrd>1%` == 0 & `mrdNegative` == 0 & `mrdNot measured` == 0 ~ 1,
      TRUE ~ 0
    ),

    # NCI high risk creation (for B-ALL)
    # define NCI high risk using the following criteria
    # - WBC count >= 50k
    #        or
    # - dx age < 1 or dx age >= 10
    NCI_high = case_when(
      lineageT == 1             ~ 0, # only B-cell eligible
      `wbc_cat>=50` == 1        ~ 1,
      age_Dx < 1 | age_Dx >= 10 ~ 1,
      TRUE ~ 0
    ),

    # Determine AALL1731 eligibility from NEJM Fig S3 contained in supplement
    AALL1731_elig = case_when(
      # exclusions
      lineageT == 1 ~ 0, # only B-cell eligible
      `cns_statusCNS3` == 1 ~ 0, # CNS3 not eligible
      `mrdNot measured` == 1 ~ 0, # MRD not measured
      rowSums(dplyr::across(dplyr::starts_with("cytogenetics"))) == 0 ~ 0, # excludes 9;22 translocation since that is the reference group and dummy columns for the variable will sum to 0

      # eligibility
      (`age_Dx` >= 1 & `age_Dx` < 10) & (`wbc_cat>=50` == 0) ~ 1, # age 1-9.9 with wbc<50
      TRUE ~ 0 # all other cases not eligible
    ),

    # Determine AALL1731 risk stratifications
    AALL1731_risk = case_when(
      # empty for not applicables
      AALL1731_elig == 0 ~ NA_character_,

      # SR-high
      (`cytogenetics_catUnfavorable` == 1) |
        (`mrd>1%` == 1 | `mrd0.1%-1%` == 1) |
        (`cytogenetics_catNeutral` == 1 & `cns_statusCNS2` == 1) ~ "High",

      # SR-favourable
      (`cytogenetics_catFavorable` == 1) &
        (`cns_statusCNS2` == 1 | (`cns_statusCNS2` == 0 & `cns_statusCNS3` == 0)) &
        (`mrdNegative` == 1) ~ "Fav",

      # SR-avg
      TRUE ~ "Avg"
    ),

    # Make dummy variables for AALL1731_risk
    # AALL1731_riskFav = ifelse(AALL1731_risk == "Fav", 1, 0), # exclude reference group
    AALL1731_riskHigh = ifelse(AALL1731_risk == "High", 1, 0),
    AALL1731_riskAvg = ifelse(AALL1731_risk == "Avg", 1, 0)
  )
