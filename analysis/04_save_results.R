#' -----------------------------------------------------------------------------
#' Calculate total costs and QALYs for the ALL model
#'
#' This script calculates the total costs and quality-adjusted life years
#' (QALYs) for the Acute Lymphoblastic Leukemia (ALL) model.
#' It processes simulation data, computes costs based on treatment type,
#' and summarizes results for different treatment scenarios.
#' It also handles specific cases for Blinatumomab as a second- or first-line
#' treatment, and calculates life expectancy metrics.
#'
#'
#' Dependencies:
#' - `chemo_costs()`: Calculates BSC chemotherapy and Blinatumomab costs
#' - `admiral::compute_bsa()`: Computes body surface area (BSA)
#'
#' -----------------------------------------------------------------------------

# Get body surface area information if evaluating Blinatumomab or BSC since it
# is needed for calculating the costs of the treatment
if (is_TRUE(blin_2ndline) || is_TRUE(bsc_2ndline) || is_TRUE(bsc_1stline) || is_TRUE(blin_1stline) || is_TRUE(blin_CEres)) {

  # Load age at diagnosis, sex, and BSA info at baseline
  # - save as vectors
  cohort_info <- read.csv(paste0(filepath, path.out, "cohort_info.csv"))
  age_base <- cohort_info$age
  male_base <- cohort_info$sexM
  bsa <- cohort_info$bsa

  # Load WHO height and weight data
  age_height <- read.csv("model-parameters/raw-data/WHOCan_hgt_wgt.csv")
  age_height$sex_num <- 0
  age_height$sex_num[age_height$sex == "male"] <- 1

  # Add height and weight to cohort_info so IVIG costs can be calculated
  cohort_info <- left_join(
    cohort_info |> mutate(age = floor(age)),
    age_height,
    by = c("age" = "age", "sexM" = "sex_num")
  )
  weight_base <- cohort_info$weight
}

# Initialize output objects ----------------------------------------------------

tc.i <- c(); tc.cat.i <- list()
le_i <- le_5i <-le_10i <- le_20i <- le_30i <- le_40i <- c()
le_iREL <- le_5iREL <- le_10iREL <- le_20iREL <- le_30iREL <- le_40iREL <- c()
le_iBMT <- le_5iBMT <- le_10iBMT <- le_20iBMT <- le_30iBMT <- le_40iBMT <- c()
m_bsa <- c()

tc.i_p10 <- tc.i_m10 <- tc.i_ivig <- c() # discounted costs +/- 10% blina price
tc.i.ud <- c(); tc.cat.i.ud <- list() # undiscounted costs

# Store event probabilities -------

time_points <- 1:10  # time period to calculate event probabilities
col_names <- c(
  "sim",
  "relapse_any", "bmt_any",
  paste0("relapse_", time_points, "yr"),
  paste0("bmt_", time_points, "yr"),
  paste0("death_", time_points, "yr")
)

event_df <- matrix(
  nrow = n_sim,
  ncol = length(col_names),
  dimnames = list(NULL, col_names)
)

# 1. Total costs ---------------------------------------------------------------

cat("\n")
for (i in 1:n_sim) {

  cat("\rCalculating total costs on simulation:  ", i, "\r")
  gc()

  data.i  <- readRDS(paste0(output.path, "data.", i, ".rds"))
  costs.i <- readRDS(paste0(output.path, "cost.", i, ".rds"))

  ## Calculate life expectancy -------------------------------------------------

  se        <- rowSums(data.i != "DEAD") * c_l
  le_i[i]   <- mean(se )
  le_5i[i]  <- mean(se > 5 )
  le_10i[i] <- mean(se > 10 )
  le_20i[i] <- mean(se > 20 )
  le_30i[i] <- mean(se > 30 )
  le_40i[i] <- mean(se > 40 )

  # Handle life expectancy for relapse
  id.REL1 <- rowSums(data.i == "REL1") > 0
  id.REL2 <- rowSums(data.i == "BMT1_REL") > 0
  id.REL22 <- (id.REL1 == F & id.REL2 == T)

  data.i2  <- data.i[id.REL1, ]
  data.i22 <- data.i[id.REL22, ]
  data.i2[data.i2   == "PRE"]  = NA
  data.i2[data.i2   == "PRE"]  = NA
  data.i22[data.i22 == "PRE"]  = NA
  data.i22[data.i22 == "BMT1"] = NA

  data.i2 <- rbind(data.i2,data.i22)


  seREL        <- rowSums(data.i2 != "DEAD", na.rm = T) * c_l
  le_iREL[i]   <- mean(seREL )
  le_5iREL[i]  <- mean(seREL > 5 )
  le_10iREL[i] <- mean(seREL > 10 )
  le_20iREL[i] <- mean(seREL > 20 )
  le_30iREL[i] <- mean(seREL > 30 )
  le_40iREL[i] <- mean(seREL > 40 )

  # Handle life expectancy for BMT
  id.BMT  <- rowSums(data.i == "BMT1") > 0
  data.i3 <- data.i[id.BMT, ]

  data.i3[data.i3 == "PRE"] <- NA

  seBMT         <- rowSums (data.i3 != "DEAD", na.rm = T) * c_l
  le_iBMT[i]   <- mean(seBMT )
  le_5iBMT[i]  <- mean(seBMT > 5 )
  le_10iBMT[i] <- mean(seBMT > 10 )
  le_20iBMT[i] <- mean(seBMT > 20 )
  le_30iBMT[i] <- mean(seBMT > 30 )
  le_40iBMT[i] <- mean(seBMT > 40 )

  # calculate costs depending on the treatment
  # - blina 2nd line after relapse vs. bsc
  # - blina 1st line vs. bsc
  # - all other treatments
  if (is_TRUE(blin_2ndline) || is_TRUE(bsc_2ndline)) { # Blinatumomab as 2nd line (JNCI paper) ------------------------

    when.REL11 <- rowSums(data.i == "REL1") > 0
    when.REL12 <- rowSums(data.i == "PRE") < 36
    # sum(when.REL11 * when.REL12)

    id.REL1 <- when.REL11 #when.REL12    # identify those that progressed in the first 3 years
    #print( sum(id.REL1))
    #print( sum(when.REL11))
    #print( sum(when.REL12))

    age.REL1 <- age_base[id.REL1 ]
    male.REL1 <- male_base[id.REL1 ]
    costs.i.REL1 <- costs.i[id.REL1,]
    data.i.REL1  <- data.i[id.REL1,]
    costs.i.REL1[data.i.REL1 == "PRE"] <- NA  ###Because Brown trial 2nd line is related to post-relapse events
    costs.i.REL1   <- t(na_move(data.frame(t(costs.i.REL1)), direction = "bottom")) # add the NA moving here
  
    # Calculate age at relapse
    age.tREL1 <- age.REL1 + (rowSums(is.na(costs.i.REL1))*30)/365.25
    age.tREL1 <- round(age.tREL1 )
    agesex.REL <- data.frame(age = age.tREL1, sex = male.REL1)

    # Calculate BSA at the time of relapse
    aswh <- left_join(agesex.REL, age_height, join_by(age == age, sex == sex_num))
    aswh <- aswh[, c("age", "sex", "height", "weight")]
    bsa.REL1 <- admiral::compute_bsa(aswh$height, aswh$weight, method = "Mosteller")
    m_bsa[i] <- mean(bsa.REL1, na.rm = T)

    costs.i.REL1[is.na(costs.i.REL1)] <- 0      # set all NAs to 0 cost
    costs.i.REL1 <- (as.matrix(costs.i.REL1))

    # Determine discounted and undiscounted costs ------------------------------

    data.tc.i <- data.frame(
      cost_d  = as.matrix(costs.i.REL1) %*% (v_dwc), # discount costs
      cost_ud = rowSums(costs.i.REL1)                # undiscounted costs
    )

    # Calculate total costs for blina or bsc using bsa at relapse `bsa.REL1`
    if (is_TRUE(blin_2ndline)) {
      data.tc.i$cost_d <- data.tc.i$cost_d + chemo_costs(bsa.REL1, trt_output = "blin")
      data.tc.i$cost_ud <- data.tc.i$cost_ud + chemo_costs(bsa.REL1, trt_output = "blin")
    } else if (is_TRUE(bsc_2ndline)) {
      data.tc.i$cost_d <- data.tc.i$cost_d + chemo_costs(bsa.REL1, trt_output = "bsc")
      data.tc.i$cost_ud <- data.tc.i$cost_ud + chemo_costs(bsa.REL1, trt_output = "bsc")
    }

    # Remove extreme outliers -------------
    # - we remove individuals with extreme costs outside the 99% CI's

    data.tc.i <- trim_outliers(data.tc.i, prob = c(0.01, 0.99))

  } else if (is_TRUE(blin_1stline) || is_TRUE(bsc_1stline)) { # First-line blina -------------------------

    # Determine discounted and undiscounted costs ------------------------------

    data.tc.i <- data.frame(
      cost_d  = as.matrix(costs.i) %*% (v_dwc), # discount costs
      cost_ud = rowSums(costs.i)                # undiscounted costs
    )

    #* Add additional Blina costs for 1st line treatment to the total costs
    #* - we dont need to do this for SoC because it is already captured in the
    #*   costing estimates from ICES
    #*
    #* Note:
    #*  - we incorporate uncertainty in the p_* parameters by sampling from a
    #*    uniform distribution for each simulation using trial estimates
    if (is_TRUE(blin_1stline)) {

      # Determine the blina treatment costs for each individual ----------------

      blina_trt_cost <- chemo_costs(
        bsa = bsa,
        blina_price = PRICE_BLIN,
        p_CRS = p_CRS[i],
        p_FN = p_FN[i],
        dur_FN = dur_FN[i],
        p_SEP = p_SEP[i],
        dur_SEP = dur_SEP[i],
        trt_output = "blin"
      )

      # Add blina-associated costs to total discounted and undiscounted costs
      # - blina is applied in the first cycle, so no discounting needed
      data.tc.i$cost_d  <- data.tc.i$cost_d  + blina_trt_cost
      data.tc.i$cost_ud <- data.tc.i$cost_ud + blina_trt_cost

      # Scenarios --------------------------------------------------------------
      # 1) +/- 10% blina price
      # 2) IVIG implementation

      # 1) +/- 10% of blina price
      blina_trt_cost_p10 <- chemo_costs(
        bsa = bsa,
        blina_price = PRICE_BLIN * 1.1, # +10%
        p_CRS = p_CRS[i],
        p_FN = p_FN[i],
        dur_FN = dur_FN[i],
        p_SEP = p_SEP[i],
        dur_SEP = dur_SEP[i],
        trt_output = "blin"
      )

      blina_trt_cost_m10 <- chemo_costs(
        bsa = bsa,
        blina_price = PRICE_BLIN * 0.9, # -10%
        p_CRS = p_CRS[i],
        p_FN = p_FN[i],
        dur_FN = dur_FN[i],
        p_SEP = p_SEP[i],
        dur_SEP = dur_SEP[i],
        trt_output = "blin"
      )

      # Save costs for +/- 10% price scenarios to the discounted costs
      # - subtract away `blina_trt_cost` to avoid double counting
      data.tc.i$cost_d_p10 <- data.tc.i$cost_d  + (blina_trt_cost_p10 - blina_trt_cost)
      data.tc.i$cost_d_m10 <- data.tc.i$cost_d  + (blina_trt_cost_m10 - blina_trt_cost)

      # 2) IVIG implementation
      # - 2 cycles of IVIG at 0.7g/kg per infusion
      ivig_cost <- (weight_base * 0.7 * PRICE_IVIG) * 2
      blina_trt_cost_ivig <- chemo_costs(
        bsa = bsa,
        blina_price = PRICE_BLIN,
        p_CRS = p_CRS[i],
        p_FN = p_FN[i] * IVIG_REDUC, # reduced risk of FN with IVIG
        dur_FN = dur_FN[i],
        p_SEP = p_SEP[i],
        dur_SEP = dur_SEP[i],
        trt_output = "blin"
      ) + ivig_cost # tack on IVIG costs

      # Save costs for IVIG scenario
      data.tc.i$cost_d_ivig <- data.tc.i$cost_d  + (blina_trt_cost_ivig - blina_trt_cost)

    }

    if (is_TRUE(bsc_1stline)) {
      # For Scenario analysis need to copy over the base costs
      data.tc.i$cost_d_p10  <- data.tc.i$cost_d
      data.tc.i$cost_d_m10  <- data.tc.i$cost_d
      data.tc.i$cost_d_ivig <- data.tc.i$cost_d
    }

    # Relapse scenarios when blina first line is given -------------------------

    if (is_TRUE(blin_1stline_relapse_blin) || is_TRUE(blin_1stline_relapse_bsc)) {

      # Calculate information at the time of relapse
      id.REL1 <- rowSums(data.i == "REL1") > 0

      # Determine when relapse occurs
      relapse_time_idx <- apply(data.i[id.REL1, ], 1, function(x) {
        idx <- which(x %in% "REL1")
        idx[1] # First occurrence of relapse
        #times[idx[1]] # First occurrence of relapse time in days
      })

      # Determine the occurrence of the first relapse time in years
      relapse_time <- times[relapse_time_idx] / 30 / 12  # Convert to years
      
      # Get the discount factor weights to be applied at the time of relapse
      v_dwc_rel <- v_dwc[relapse_time_idx]

      # Calculate age at relapse
      age.tREL1  <- age_base[id.REL1] + relapse_time

      # Set maximum age to match available in height/weight data
      age.tREL1[age.tREL1 > max(age_height$age)] <- max(age_height$age)

      # Age/sex for relapsed individuals
      agesex.REL <- data.frame(
        age = round(age.tREL1),
        sex = male_base[id.REL1]
      )

      # Calculate BSA at the time of relapse
      aswh <- left_join(agesex.REL, age_height, join_by(age == age, sex == sex_num))
      aswh <- aswh[, c("age", "sex", "height", "weight")]
      bsa.REL1 <- admiral::compute_bsa(aswh$height, aswh$weight, method = "Mosteller")

      # Calculate additional treatment costs at the time of relapse
      # - blina at relapse
      # - bsc at relapse
      if (is_TRUE(blin_1stline_relapse_blin)) {

        rel_cost_v_ud <- chemo_costs(bsa.REL1, blina_price = PRICE_BLIN, trt_output = "blin")  # undiscounted
        rel_cost_v_d  <- rel_cost_v_ud * v_dwc_rel       # discounted

        # Scenarios:
        # 1) +/- 10% of blina price in relapse --------------------------

        # Determine the relapse costs for +/- 10% blina price scenarios
        rel_cost_v_d_p10 <- chemo_costs(bsa.REL1, trt_output = "blin", blina_price = PRICE_BLIN * 1.1) * v_dwc_rel # +10%
        rel_cost_v_d_m10 <- chemo_costs(bsa.REL1, trt_output = "blin", blina_price = PRICE_BLIN * 0.9) * v_dwc_rel # -10%

        # Add to the discounted costs for +/- 10% blina price scenarios
        data.tc.i$cost_d_p10[id.REL1] <- data.tc.i$cost_d_p10[id.REL1] + rel_cost_v_d_p10
        data.tc.i$cost_d_m10[id.REL1] <- data.tc.i$cost_d_m10[id.REL1] + rel_cost_v_d_m10

        # 2) IVIG scenario -----------------------------------------------------

        # Calculate IVIG costs at relapse (discounted)
        ivig_cost <- ((weight_base[id.REL1] * 0.7 * PRICE_IVIG) * 2) * v_dwc_rel
        rel_cost_v_d_ivig <- rel_cost_v_d + ivig_cost

        # Add to the discounted costs for IVIG scenario
        data.tc.i$cost_d_ivig[id.REL1] <- data.tc.i$cost_d_ivig[id.REL1] + rel_cost_v_d_ivig

      } else if (is_TRUE(blin_1stline_relapse_bsc)) {

        rel_cost_v_ud <- chemo_costs(bsa.REL1, trt_output = "bsc") # undiscounted
        rel_cost_v_d <- rel_cost_v_ud * v_dwc_rel                  # discounted

      }

      # Update total costs for individuals who experienced relapse
      data.tc.i$cost_d[id.REL1]  <- data.tc.i$cost_d[id.REL1]  + rel_cost_v_d
      data.tc.i$cost_ud[id.REL1] <- data.tc.i$cost_ud[id.REL1] + rel_cost_v_ud

    }

    # Remove extreme outliers -------------
    # - we remove individuals with extreme costs outside the 99% CI's
    data.tc.i <- trim_outliers(data.tc.i, prob = c(0.01, 0.99))

  } else { # all other treatment scenarios -------------------------------------

    data.tc.i <- data.frame(
      cost_d  = as.matrix(costs.i) %*% (v_dwc), # discount costs
      cost_ud = rowSums(costs.i)                # undiscounted costs
    )

    # Remove extreme outliers -------------
    # - we remove individuals with extreme costs outside the 99% CI's
    data.tc.i <- trim_outliers(data.tc.i, prob = c(0.01, 0.99))

  }

  # Get total costs overall ----------------------------------------------------

  tc.i[i] <- (data.tc.i %>% summarise(mean = mean(cost_d, na.rm=T))) # total costs discounted
  tc.i.ud[i] <- (data.tc.i %>% summarise(mean = mean(cost_ud, na.rm=T))) # total costs undiscounted

  # Scenario analysis total costs
  if (is_TRUE(blin_1stline) || is_TRUE(bsc_1stline)) {
    tc.i_p10[i] <- (data.tc.i %>% summarise(mean = mean(cost_d_p10, na.rm=T))) # total costs discounted +10% blina price
    tc.i_m10[i] <- (data.tc.i %>% summarise(mean = mean(cost_d_m10, na.rm=T))) # total costs discounted -10% blina price
    tc.i_ivig[i] <- (data.tc.i %>% summarise(mean = mean(cost_d_ivig, na.rm=T))) # total costs discounted IVIG scenario
  }

  # Calculate event counts for reporting ---------------------------------------
  # - loop through the times of interest to calculate event proportions
  for (t in c(0, time_points)) {
    event_df[i, "sim"] <- i

    # time = 0 captures any event during the entire simulation
    if (t == 0) {
      event_df[i, "relapse_any"] <- event_risk(data = data.i, event = c("REL1"), event_time = NULL)
      event_df[i, "bmt_any"]     <- event_risk(data = data.i, event = c("BMT1", "BMT1_REL"), event_time = NULL)
    } else {
      event_df[i, paste0("relapse_", t, "yr")] <- event_risk(data = data.i, event = c("REL1"), event_time = t)
      event_df[i, paste0("bmt_", t, "yr")]     <- event_risk(data = data.i, event = c("BMT1", "BMT1_REL"), event_time = t)
      event_df[i, paste0("death_", t, "yr")]   <- event_risk(data = data.i, event = c("DEAD"), event_time = t)
    }
  }

}

# Life expectancy overall ------------------------------------------------------

les <- list(le_i = le_i,le_5i = le_5i,le_10i =le_10i ,le_20i = le_20i,le_30i = le_30i,le_40i= le_40i)
lesREL <- list(le_iREL = le_iREL,le_5iREL = le_5iREL,le_10iREL =le_10iREL ,le_20iREL = le_20iREL,le_30iREL = le_30iREL,le_40iREL= le_40iREL)
lesBMT <- list(le_iBMT = le_iBMT,le_5iBMT = le_5iBMT,le_10iBMT =le_10iBMT ,le_20iBMT = le_20iBMT,le_30iBMT = le_30iBMT,le_40iBMT= le_40iBMT)


# Total costs overall ----------------------------------------------------------

tc.i <- data.frame(mean = unlist(tc.i)) # per iteration
tc <- summarize_sim_output(tc.i) # summarized

tc.i.ud <- data.frame(mean = unlist(tc.i.ud))
tc.ud <- summarize_sim_output(tc.i.ud)

# Scenario analysis:
if (is_TRUE(blin_1stline) || is_TRUE(bsc_1stline)) {
  tc.i_p10 <- data.frame(mean = unlist(tc.i_p10))
  tc.p10 <- summarize_sim_output(tc.i_p10)

  tc.i_m10 <- data.frame(mean = unlist(tc.i_m10))
  tc.m10 <- summarize_sim_output(tc.i_m10)

  tc.i_ivig <- data.frame(mean = unlist(tc.i_ivig))
  tc.ivig <- summarize_sim_output(tc.i_ivig)
}

if (is_TRUE(blin_2ndline)) {

  saveRDS(tc.i, paste0(output.path, "tc.i_new_blin",  ".rds"), compress =T)
  saveRDS(tc.i.ud, paste0(output.path, "tc.i.ud_new_blin",  ".rds"), compress =T)
  saveRDS(tc, paste0(output.path, "tc_hat_blin", n_t, "year.rds"), compress =T)
  saveRDS(tc.ud, paste0(output.path, "tc_hat_ud_blin", n_t, "year.rds"), compress =T)

} else if (is_TRUE(bsc_2ndline)) {

  saveRDS(tc.i, paste0(output.path, "tc.i_new_bsc",  ".rds"), compress =T)
  saveRDS(tc.i.ud, paste0(output.path, "tc.i.ud_new_bsc",  ".rds"), compress =T)
  saveRDS(tc, paste0(output.path, "tc_hat_bsc", n_t, "year.rds"), compress =T)
  saveRDS(tc.ud, paste0(output.path, "tc_hat_ud_bsc", n_t, "year.rds"), compress =T)

} else if (is_TRUE(blin_1stline)) {

  saveRDS(tc.i,  paste0(output.path, "tc.i_new_blin",  ".rds"), compress =T)
  saveRDS(tc.i.ud, paste0(output.path, "tc.i.ud_new_blin",  ".rds"), compress =T)
  saveRDS(tc,    paste0(output.path, "tc_hat_blin", n_t, "year.rds"), compress =T)
  saveRDS(tc.ud,    paste0(output.path, "tc_hat_ud_blin", n_t, "year.rds"), compress =T)
  
  ### Scenario analysis:
  # +/- 10% blina price
  saveRDS(tc.i_p10, paste0(output.path, "tc.i_new_blin_p10",  ".rds"), compress =T)
  saveRDS(tc.p10, paste0(output.path, "tc_hat_blin_p10_", n_t, "year.rds"), compress =T)
  saveRDS(tc.i_m10, paste0(output.path, "tc.i_new_blin_m10",  ".rds"), compress =T)
  saveRDS(tc.m10, paste0(output.path, "tc_hat_blin_m10_", n_t, "year.rds"), compress =T)
    # IVIG implementation
  saveRDS(tc.i_ivig, paste0(output.path, "tc.i_new_blin_ivig",  ".rds"), compress =T)
  saveRDS(tc.ivig, paste0(output.path, "tc_hat_blin_ivig_", n_t, "year.rds"), compress =T)
  
  saveRDS(les,   paste0(output.path, "les",  ".rds"), compress =T) 
  saveRDS(lesREL,paste0(output.path, "lesREL",  ".rds"), compress =T) 
  saveRDS(lesBMT,paste0(output.path, "lesBMT",  ".rds"), compress =T) 
  saveRDS(event_df, paste0(output.path, "event_probs.rds"), compress =T)

} else if (is_TRUE(bsc_1stline)) {

  saveRDS(tc.i,  paste0(output.path, "tc.i_new_bsc",  ".rds"), compress =T)
  saveRDS(tc.i.ud, paste0(output.path, "tc.i.ud_new_bsc",  ".rds"), compress =T)
  saveRDS(tc,    paste0(output.path, "tc_hat_bsc", n_t, "year.rds"), compress =T)
  saveRDS(tc.ud, paste0(output.path, "tc_hat_ud_bsc", n_t, "year.rds"), compress =T)

  ### Scenario analysis:
  # +/- 10% blina price
  saveRDS(tc.i_p10, paste0(output.path, "tc.i_new_bsc_p10",  ".rds"), compress =T)
  saveRDS(tc.p10, paste0(output.path, "tc_hat_bsc_p10_", n_t, "year.rds"), compress =T)
  saveRDS(tc.i_m10, paste0(output.path, "tc.i_new_bsc_m10",  ".rds"), compress =T) 
  saveRDS(tc.m10, paste0(output.path, "tc_hat_bsc_m10_", n_t, "year.rds"), compress =T)
  # IVIG implementation
  saveRDS(tc.i_ivig, paste0(output.path, "tc.i_new_bsc_ivig",  ".rds"), compress =T)
  saveRDS(tc.ivig, paste0(output.path, "tc_hat_bsc_ivig_", n_t, "year.rds"), compress =T)

  saveRDS(les,   paste0(output.path, "les",  ".rds"), compress =T) 
  saveRDS(lesREL,paste0(output.path, "lesREL",  ".rds"), compress =T) 
  saveRDS(lesBMT,paste0(output.path, "lesBMT",  ".rds"), compress =T)
  saveRDS(event_df, paste0(output.path, "event_probs.rds"), compress =T)

} else {

  saveRDS(tc.i, paste0(output.path, "tc.i_new",  ".rds"), compress =T)
  saveRDS(tc.i.ud, paste0(output.path, "tc.i.ud_new",  ".rds"), compress =T)
  saveRDS(tc, paste0(output.path, "tc_hat_", n_t, "year.rds"), compress =T)
  saveRDS(tc.ud, paste0(output.path, "tc_hat_ud_", n_t, "year.rds"), compress =T)

  saveRDS(les,   paste0(output.path, "les",  ".rds"), compress =T) 
  saveRDS(lesREL,paste0(output.path, "lesREL",  ".rds"), compress =T) 
  saveRDS(lesBMT,paste0(output.path, "lesBMT",  ".rds"), compress =T) 
}


# 2. Total effects -------------------------------------------------------------

te.i <- le.i <- c(); te.cat.i <- le.cat.i <- list()

te.i.ud <- le.i.ud <- c(); te.cat.i.ud <- le.cat.i.ud <- list() # undiscounted effects

te.i_ivig <- le.i_ivig <- c() # IVIG scenario

cat("\n")
for (i in 1:n_sim) {

  cat("\rCalculating effects on simulation:  ", i, "\r")

  effect.i <- readRDS(paste0(output.path, "effect.", i, ".rds"))

  if (is_TRUE(blin_2ndline) || is_TRUE(bsc_2ndline)) {
    # subset to first relapse people, and only calculate the effect of those in the first relapse, 
    data.i  <- readRDS(paste0(output.path, "data.", i, ".rds"))
    when.REL11 <- rowSums(data.i == "REL1") > 0

    id.REL1 <- when.REL11
    effect.i.REL1 <- effect.i[id.REL1, ]
    data.i.REL1 <- data.i[id.REL1, ]
    effect.i.REL1[data.i.REL1=="PRE"] <- NA                                 # set to NA time spent in PRE 
    effect.i.REL1<- t(na_move(data.frame(t(effect.i.REL1)), direction = "bottom")) # move all NA's to bottom
    effect.i.REL1[is.na(effect.i.REL1)] = 0                                 # all the time in PRE set it to NA

    data.te.i <- effect.i.REL1 %*% (v_dwe) # discounted QALYs
    data.te.i.ud <- rowSums(effect.i.REL1) # undiscounted QALYs

    data.le.i <- rowSums(as.matrix(effect.i.REL1)>0) * 30/365.25
    data.le.i.ud <- data.le.i

  } else {

    data.te.i <- as.matrix(effect.i) %*% (v_dwe) # discounted QALYs
    data.te.i.ud <- rowSums((effect.i)) # undiscounted QALYs

    data.le.i <- rowSums(as.matrix(effect.i)>0) * 30/365.25
    data.le.i.ud <- data.le.i

    # Scenario analysis for IVIG implementation
    # - create copies of the base case effects so i can save them later
    data.te.i_ivig <- data.te.i
    data.le.i_ivig <- data.le.i

    # Incorporate disutilities for added adverse events in the blina arm
    # - only if blina is being evaluated as first line treatment
    if (is_TRUE(blin_1stline)) {
      data.te.i <- data.te.i - (dis_FN * p_FN[i] * (dur_FN[i]/365.25)) # disutility from FN
      data.te.i <- data.te.i - (dis_SEP * p_SEP[i] * (dur_SEP[i]/365.25)) # disutility from sepsis

      data.te.i.ud <- data.te.i.ud - (dis_FN * p_FN[i] * (dur_FN[i]/365.25)) # disutility from FN
      data.te.i.ud <- data.te.i.ud - (dis_SEP * p_SEP[i] * (dur_SEP[i]/365.25)) # disutility from sepsis

      # Scenario: IVIG implementation
      data.te.i_ivig <- data.te.i_ivig - (dis_FN * p_FN[i] * IVIG_REDUC * (dur_FN[i]/365.25)) # disutility from FN
      data.te.i_ivig <- data.te.i_ivig - (dis_SEP * p_SEP[i] * (dur_SEP[i]/365.25)) # disutility from sepsis
    }
  }

  data.te.i <- data.frame(V1 = data.te.i) # discounted effects
  data.te.i.ud <- data.frame(V1 = data.te.i.ud) # undiscounted effects

  data.le.i <- data.frame(V1 = data.le.i)
  data.le.i.ud <- data.frame(V1 = data.le.i.ud)

  # Remove extreme outliers -------------
  
  data.te.i <- trim_outliers(data.te.i, prob = c(0.01, 0.99))
  data.te.i.ud <- trim_outliers(data.te.i.ud, prob = c(0.01, 0.99))
  data.le.i <- trim_outliers(data.le.i, prob = c(0.01, 0.99))
  data.le.i.ud <- trim_outliers(data.le.i.ud, prob = c(0.01, 0.99))

  if (is_TRUE(blin_1stline) || is_TRUE(bsc_1stline)) {

    data.te.i_ivig <- data.frame(V1 = data.te.i_ivig) # discounted effects IVIG scenario
    data.le.i_ivig <- data.frame(V1 = data.le.i_ivig) # discounted le IVIG scenario
    data.te.i_ivig <- trim_outliers(data.te.i_ivig, prob = c(0.01, 0.99))
    data.le.i_ivig <- trim_outliers(data.le.i_ivig, prob = c(0.01, 0.99))

    te.i_ivig[i] <- data.te.i_ivig %>% summarise(mean = mean(V1, na.rm = T)) # discounted effects IVIG scenario
    
  }
  
  # Get total effects overall
  te.i[i] <- data.te.i %>% summarise(mean = mean(V1, na.rm = T)) # discounted effects
  te.i.ud[i] <- data.te.i.ud %>% summarise(mean = mean(V1, na.rm = T)) # undiscounted effects

}

# Total effects overall
te.i <- data.frame(mean = unlist(te.i)) # discounted
te.i.ud <- data.frame(mean = unlist(te.i.ud)) # undiscounted

te <- summarize_sim_output(te.i)
te.ud <- summarize_sim_output(te.i.ud)

# Scenario analysis IVIG implementation
if (is_TRUE(blin_1stline) || is_TRUE(bsc_1stline)) {
  te.i_ivig <- data.frame(mean = unlist(te.i_ivig)) # discounted effects
  te.ivig <- summarize_sim_output(te.i_ivig)
}

if (is_TRUE(blin_2ndline)) {

  saveRDS(te.i, paste0(output.path, "te.i_new_blin",  ".rds"), compress =T)
  saveRDS(te.i.ud, paste0(output.path, "te.i.ud_new_blin",  ".rds"), compress =T)
  saveRDS(te,  paste0(output.path, "te_hat_blin", n_t, "year.rds"), compress =T)
  saveRDS(te.ud,  paste0(output.path, "te_hat_ud_blin", n_t, "year.rds"), compress =T)

} else if (is_TRUE(bsc_2ndline)) {

  saveRDS(te.i, paste0(output.path, "te.i_new_bsc",  ".rds"), compress =T)
  saveRDS(te.i.ud, paste0(output.path, "te.i.ud_new_bsc",  ".rds"), compress =T)
  saveRDS(te, paste0(output.path, "te_hat_bsc", n_t, "year.rds"), compress =T)
  saveRDS(te.ud, paste0(output.path, "te_hat_ud_bsc", n_t, "year.rds"), compress =T)

} else if (is_TRUE(blin_1stline)) {

  saveRDS(te.i, paste0(output.path, "te.i_new_blin",  ".rds"), compress =T)
  saveRDS(te.i.ud, paste0(output.path, "te.i.ud_new_blin",  ".rds"), compress =T)
  saveRDS(te,  paste0(output.path, "te_hat_blin", n_t, "year.rds"), compress =T)
  saveRDS(te.ud,  paste0(output.path, "te_hat_ud_blin", n_t, "year.rds"), compress =T)

  saveRDS(te.i_ivig, paste0(output.path, "te.i_ivig_new_blin",  ".rds"), compress =T)
  saveRDS(te.ivig, paste0(output.path, "te_hat_ivig_blin_", n_t, "year.rds"), compress =T)

} else if (is_TRUE(bsc_1stline)) {

  saveRDS(te.i, paste0(output.path, "te.i_new_bsc",  ".rds"), compress =T)
  saveRDS(te.i.ud, paste0(output.path, "te.i.ud_new_bsc",  ".rds"), compress =T)
  saveRDS(te, paste0(output.path, "te_hat_bsc", n_t, "year.rds"), compress =T)
  saveRDS(te.ud, paste0(output.path, "te_hat_ud_bsc", n_t, "year.rds"), compress =T)

  saveRDS(te.i_ivig, paste0(output.path, "te.i_ivig_new_bsc",  ".rds"), compress =T)
  saveRDS(te.ivig, paste0(output.path, "te_hat_ivig_bsc_", n_t, "year.rds"), compress =T)

} else {

  saveRDS(te.i, paste0(output.path, "te.i_new",  ".rds"), compress =T)
  saveRDS(te.i.ud, paste0(output.path, "te.i.ud_new",  ".rds"), compress =T)
  saveRDS(te, paste0(output.path, "te_hat_", n_t, "year.rds"), compress =T)
  saveRDS(te.ud, paste0(output.path, "te_hat_ud_", n_t, "year.rds"), compress =T)

}

cat("\n")
