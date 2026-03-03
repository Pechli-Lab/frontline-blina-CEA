#' Update transition probabilities based on health state occupied at j
#'
#' @param M_it matrix with individuals health state at j
#' @param X    dataframe of covariate values from individuals
#'
#' @return transiton probabilities
#' 
#' @note Global variables used (not passed in as arguments):
#'   - n_s: number of states
#'   - n_ind: number of individuals
#'   - v_n: vector of state names
#'   - t: current time
#'   - k: current iteration of the simulation
#'   - p_death_f: probability of death for females (list)
#'   - p_death_m: probability of death for males (list)
#'   - bHr.norm: hazard ratio values for blinatumomab 2nd line analysis
#'   - b1_HR vars: hazard ratio values for blinatumomab 1st line analysis
#'   - p.PRE_DEAD:  trans prob from PRE --> DEAD
#'   - p.PRE_REL1:  trans prob from PRE --> REL1
#'   - p.PRE_BMT1:  trans prob from PRE --> BMT1
#'   - p.REL1_DEAD: trans prob from REL1 --> DEAD
#'   - p.REL1_REL2: trans prob from REL1 --> REL2
#'   - p.REL1_BMT1: trans prob from REL1 --> BMT1
#'   - p.REL2_DEAD: trans prob from REL2 --> DEAD
#'   - p.REL2_BMT1: trans prob from REL2 --> BMT1
#'   - p.BMT1_DEAD: trans prob from BMT1 --> DEAD
#'   - p.BMT1_BMT1_REL: trans prob from BMT1 --> BMT1_REL
#'   - p.BMT1_REL_DEAD: trans prob from BMT1_REL --> DEAD
#'
Probs <- function(M_it, X = X, blin_2ndline = FALSE, blin_1stline = FALSE){

  p.it <- matrix(0, nrow = n_s, ncol = n_ind)   # create vector of state transition probabilities
  rownames(p.it) <- v_n                         # name the rows


  for (i in 1:length(v_n[v_n != "DEAD"])){
    v_n.temp <- v_n[v_n != "DEAD"]
    assign(paste0("pos.", v_n.temp[i]), M_it == v_n[i])
  }

  # * Update age
  temp.age <- floor(X$age_Dx + t/30/12)

  # probability of dying from other causes
  p_DEATH_add <- p_death_f[[k]][temp.age + 1]
  p_DEATH_add[X$sexM == 1] <- p_death_m[[k]][temp.age + 1][X$sexM == 1]


  # HR of dying for blinatumomab (2nd line in relapse/refractory disease)
  # - sample from log hazard dist (normal) of Brown et al paper
  # - These are based on a HR in the Brown paper that related to PFS
  # - adjusts the probability of dying or relapsing (again)
  if (is_TRUE(blin_2ndline) || is_TRUE(blin_1stline_relapse_blin)) {
   # p.PRE_DEAD      <- 1 - exp(log(1 - p.PRE_DEAD)      * bHr.norm[k])
    p.REL1_DEAD      <- 1 - exp(log(1 - p.REL1_DEAD)     * bHr.norm[k])
    p.REL1_REL2      <- 1 - exp(log(1 - p.REL1_REL2)     * bHr.norm[k])
    p.REL1_BMT1      <- 1 - exp(log(1 - p.REL1_BMT1)     )
   # p.BMT1_DEAD     <- 1 - exp(log(1 - p.BMT1_DEAD)     * bHr.norm[k])
   # p.BMT1_REL_DEAD <- 1 - exp(log(1 - p.BMT1_REL_DEAD) * bHr.norm[k])
  }
  

  # Blinatumomab (1st line therapy)
  # - time-varying hazard ratios based on AALL1731 trial (S1 Fig, Gupta et al)
  #   (extracted via SurvdigitizeR)
  # - need to apply time-dependent HR at the current time 't'
  # - time of randomization (`START_BLIN`) to mimic the original trial design,
  #  which was after induction chemo
  if (is_TRUE(blin_1stline)) {
    if ( (t / START_BLIN < ncol(b1_HR_PRE_REL1) ) && (t >= START_BLIN) ) {

      # Convert `t` (days) to index position for hazard ratio dataframe
      # - columns in `b1_HR_PRE_REL1` are in 30d intervals after randomization
      idx <- (t / START_BLIN)

      # Update probabilities based on hazard ratios from trial
      # p.PRE_DEAD <- 1 - exp(log(1 - p.PRE_DEAD) * b1_HR_PRE_DEAD[k, idx]) # trial did not show any major difference in OS, so not applying HR to death probability
      p.PRE_REL1 <- 1 - exp(log(1 - p.PRE_REL1) * b1_HR_PRE_REL1[k, idx]) # transition probability to relapse affected by blina in first line, so adjusting the HR for it
    }
  }
    

  # Add addition probability of death
  # to do - sum rates instead of probabilities below
  
  p.PRE_DEAD[pos.PRE]           <- p.PRE_DEAD[pos.PRE]   + p_DEATH_add[pos.PRE]
  p.REL2_DEAD[pos.REL2]         <- p.REL2_DEAD[pos.REL2] + p_DEATH_add[pos.REL2]
  p.REL1_DEAD[pos.REL1]         <- p.REL1_DEAD[pos.REL1] + p_DEATH_add[pos.REL1]
  p.BMT1_DEAD[pos.BMT1]         <- p.BMT1_DEAD[pos.BMT1] + p_DEATH_add[pos.BMT1]
  p.BMT1_REL_DEAD[pos.BMT1_REL] <- p.BMT1_REL_DEAD[pos.BMT1_REL] + p_DEATH_add[pos.BMT1_REL]




  # Update p.it with appropriate transition probabilities
  # p.it has the transition probabilities for each individual
  p.it[, M_it == "PRE"]    <- rbind((1 - (p.PRE_REL1[pos.PRE] + p.PRE_BMT1[pos.PRE] + p.PRE_DEAD[pos.PRE])), #PRE to PRE
                                    (p.PRE_REL1[pos.PRE]),                                                   #PRE to REL1
                                    0,                                                                       #PRE to REL2
                                    p.PRE_BMT1[pos.PRE],                                                     #PRE to BMT1
                                    0,                                                                       #PRE to BMT1_REL
                                    p.PRE_DEAD[pos.PRE])                                                     #PRE to DEAD

  p.it[, M_it == "REL1"]   <- rbind(0,	                                                                             #REL1 to PRE,
                                    (1 - (p.REL1_REL2[pos.REL1] + p.REL1_BMT1[pos.REL1] + p.REL1_DEAD[pos.REL1])),	 #REL1 to REL1
                                    p.REL1_REL2[pos.REL1],                                                           #REL1 to REL2
                                    p.REL1_BMT1[pos.REL1],                                                           #REL1 to BMT1
                                    0,                                                                               #REL1 to BMT1_REL
                                    p.REL1_DEAD[pos.REL1])                                                           #REL1 to DEAD

  p.it[, M_it =="REL2"]   <- rbind(0,  0,	           	                                         #REL2 to PRE, REL1
                                    (1 - (p.REL2_BMT1[pos.REL2] + p.REL2_DEAD[pos.REL2])),       #REL2 to REL2
                                    p.REL2_BMT1[pos.REL2],                                       #REL2 to BMT1
                                    0,                                                           #REL2 to BMT1_REL
                                    p.REL2_DEAD[pos.REL2])                                       #REL2 to DEAD

  p.it[, M_it == "BMT1"]   <- rbind(0,  0,	0,                                                   #BMT1 to PRE, REL1, REL2
                                    (1 - (p.BMT1_BMT1_REL[pos.BMT1] + p.BMT1_DEAD[pos.BMT1])),   #BMT1 to BMT1
                                    p.BMT1_BMT1_REL[pos.BMT1],                                   #BMT1 to BM1REL
                                    p.BMT1_DEAD[pos.BMT1])                                       #BMT1 to DEAD

  p.it[, M_it == "BMT1_REL"] <- rbind(0,  0,	 0,	 0,                                            #BMT1_REL to PRE, REL1, REL2, BMT2
                                     (1 - (p.BMT1_REL_DEAD[pos.BMT1_REL])),                      #BMT1_REL to BMT1_REL
                                     p.BMT1_REL_DEAD[pos.BMT1_REL])                              #BMT1_REL to DEAD

  p.it[, M_it == "DEAD"] <- c(0,	0, 0, 0, 0, 1)


  temp.age <- X$age_Dx + t/30/12  # updated age - not rounded
  if(sum(temp.age > 110) > 1) p.it[, which(temp.age > 110)] <- c(0,0,0,0,0,1) # if over 100, die

  return(t(p.it))# return the transition probabilities
}
