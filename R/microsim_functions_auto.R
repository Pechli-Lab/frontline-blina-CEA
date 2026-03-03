###   Functions for ALL Simulation

options(scipen=999) # prevents scientific notation


# IMPORT MULTISTATE MODEL OUTPUTS -------------------------------------------------
load(paste0(event.path, "m_t.Rdata")) # Model matrix

# * MSM Estimates ----
MSM_est <- read.csv(paste0(event.path, "l_selected_models_by_transition_est.csv"))
# MSM_est <- MSM_est[MSM_est$transition != "model_TRT",]  # remove from MSM_est

MSM_est$model <- as.character(MSM_est$model)
MSM_est$covariate <- as.character(MSM_est$covariate)
MSM_est$transition <- factor(MSM_est$transition, levels=unique(MSM_est$transition)) # factor transition to preserve order
l_MSM_est <- split(MSM_est, f = MSM_est$transition)                                 # convert to list of dataframes based on trans


# * MSM Covariance Matrices ----
MSM_cov <- list()

for (i in 1:sum(!is.na(m_t))){    # for each transition number
  transition <- which(m_t == i)   # get transition index in m_t

  s.from <- transition %% length(v_n)     # from state = row of the transition matrix, or remainder when divided by # states
  s.to <- ceiling(transition/length(v_n)) # to state = column of the transition matrix, or rounded up division by # states

  # assign tryCatch value to x: reads in CSV if it's not empty and takes cov matrix, otherwise takes value of NA
  x <- tryCatch({
    read.csv(paste0(event.path, "cov_selected_model_", v_n[s.from], "_", v_n[s.to], ".csv"), skip = 0)
    },
    error = function(e){
      NA}
  )

  MSM_cov[[i]] <- x
  names(MSM_cov)[i] <- paste0("model_", v_n[s.from], "_", v_n[s.to])
}


# * Spline Knot Locations ----
MSM_knots <- read.csv(paste0(event.path, "df_model_knots.csv"))
MSM_knots$transition <- as.character(MSM_knots$transition)
MSM_knots$transition <- factor(MSM_knots$transition, levels = unique(MSM_knots$transition))
l_MSM_knots <- lapply(levels(MSM_knots$transition), function(x) MSM_knots[which(MSM_knots$transition == x),])


# * Add Uncertainty to event coefficients ----
# Requires transition matrix (m_t) to run but automatically runs all the transitions

norm.mat <- list()
for (i in 1:sum(!is.na(m_t))){    # for each transition number
  transition <- which(m_t == i)   # get transition index in m_t
  s.from <- transition %% length(v_n)     # from state = row of the transition matrix, or remainder when divided by # states
  s.to <- ceiling(transition/length(v_n)) # to state = column of the transition matrix, or rounded up division by # states

  assign(paste0(v_n[s.from], "_", v_n[s.to], ".norm.mat"), model.rmvnorm.f(MSM_est, MSM_cov, index = i, n_sim = n_sim))
  norm.mat[[i]] <- get(paste0(v_n[s.from], "_", v_n[s.to], ".norm.mat"))
}


# * Mortality rates from lifetables ----
mx <- read.table("model-parameters/raw-data/Mx_1x1.txt", header = T)
mx_year <- max(mx$Year)
mx <- mx[mx$Year == mx_year,]

mx$Year <- as.numeric(mx$Year)
mx$Female <- as.numeric(mx$Female)
mx$Male   <- as.numeric(mx$Male)
m_r_f  <- (1 + mx$Female)^(1/12) - 1  # calculate monthly mortality rates for latest year
m_r_m  <- (1 + mx$Male)  ^(1/12) - 1


# * Relative Risks from Yeh JM et al.----
Yeh.RR <- read.csv("model-parameters/raw-data/Yeh.RR.csv") # 10 year RR for 30, 40, 50, and 60 y.o., extracted using WebPlotDigitizer

# add CI (hazard with log normal distribution)
RR.norm <- data.frame(rep(0, n_sim))
for (i in 1:nrow(Yeh.RR)){
  RR.norm[,i] <- exp(rnorm(n_sim, log(Yeh.RR[i,2]), Yeh.RR[i,5]))
}
names(RR.norm) <- Yeh.RR$X

# * Calculate probability of death ----
# * * Treating RR as HRs - Get close results if probability of dying is small but difference increases with probability of dying
RR  <- p_death_f <- p_death_m <- list()

for(i in 1:n_sim){
  RR[[i]]  <- c(rep(1, 30),                                            # 0 to 29
                rep(RR.norm[i,1],10),                                  #30 to 39
                rep(RR.norm[i,2],10),                                  #40 to 49
                rep(RR.norm[i,3],10),                                  #50 to 59
                rep(RR.norm[i,4],length(mx$Female) - (30+10+10+10))) #60+

  p_death_f[[i]] <- (1 - exp(-m_r_f * RR[[i]])) # add CI assume normal on log scale
  p_death_m[[i]] <- (1 - exp(-m_r_m * RR[[i]]))
}


# IMPORT COST MODEL OUTPUTS ------------------------------------------------

# * Cost Estimates ----
cost_est <- read.csv(paste0(costs.path, "/", costs.trans, "_df_est_ns", costs.spline, ".csv"))


# * Initialized vectors ----
cost_reff <- cost_cov <- vector(mode = "list", length = n_s)
names(cost_cov) <- names(cost_reff) <-v_n


# * Covariance Matrix - Fixed Effects ----
for (i in 1:n_s){
  # assign tryCatch value to x: reads in CSV if it's not empty and takes cov matrix, otherwise takes value of NA
  x <- tryCatch({
    read.csv(paste0(costs.path, "/", costs.trans, "_cov_ns", costs.spline, "_", v_n[i], ".csv"), skip = 0)
  },
  error = function(e){
    NA}
  )

  cost_cov[[i]] <- x
}


# * Random effects - Covariance Matrix ----
for(i in 2:(n_s-3)){
  cost_reff[[i]] <- read.csv(paste0(costs.path, "/", costs.trans, "_r.eff_ns", costs.spline, "_", v_n[i], ".csv"))
}

# mixed model >0
l_cost_reff <- lapply(cost_reff, function(x) length(x))


# * Smearing Estimates  ----
cost_smear <- read.csv(paste0(costs.path, "/", costs.trans, "_smear_hetero_ns", costs.spline, ".csv"))

# * Spline Knot Locations ----
cost_knots <- read.csv(paste0(costs.path, "/", costs.trans, "_df_knots", costs.spline, ".csv"))

# * List of covariate names for each state  ----
l_costcovs <- list()

for (i in 1:n_s){
  l_costcovs[[i]] <- as.vector(cost_est$covariate[which(cost_est$state == v_n[i])])
  names(l_costcovs)[i] <- v_n[i]
}
names(l_costcovs) <- v_n

l_costcovs_ns <- lapply(l_costcovs, function(x) x[grep("ns", x)])


# * Add Uncertainty and Random Effects to Costs -----------------------------

# * * Uncertainty ----
norm.mat.cost <- c()

for (i in 1:n_s){
  est <- as.matrix(cost_est$estimate[cost_est$state == v_n[i]]) # Select coefficient estimates for each state
  norm.mat.cost[[i]] <- rmvnorm(n_sim, est, as.matrix(cost_cov[[i]][-1])) # sample from multivariate normal with mean = estimates, sigma = covairance matrix
  colnames(norm.mat.cost[[i]]) <- cost_est$covariate[cost_est$state == v_n[i]] # set column names to covariate names
}


# * * Random Effects ----
RE <- c()

for (i in 2:(n_s-3)){
  RE[[i]] <- rmvnorm(n_i_i*n_sim, rep(0, nrow(cost_reff[[i]])), data.matrix(cost_reff[[i]])) # random effect sample from multivariate normal
  colnames(RE[[i]]) <- rownames(cost_reff[[i]])
}


# * * Natural Cubic Splines - Get Basis Matrices
nspline <- list()
for (i in 1:(n_s-1)){
  nspline[[i]] <- ns(1:cycles, knots = c(cost_knots[i,]$I.K1, cost_knots[i,]$I.K2),  # ns df=3, get basis matrix
                     Boundary.knots = c(cost_knots[i,]$B.K1, cost_knots[i,]$B.K2))
}


# * * Calculate costs for DEAD (interval prior to death) for each sim - DEAD state costs are an empty lm model
c_dead <- exp(norm.mat.cost[[6]] + 1/2*cost_smear[[6]])


# # SMEARING ESTIMATE(homoskedastic - log)
# # exp(rowSums(as.matrix(temp.mat) %*% as.matrix(coef(model)[-1])) * smear) # Calculate adjusted cost predictions
# tempcost <- exp(tempcost * cost_smear[i])
#
# # SMEARING ESTIMATE(homoskedastic - sqrt)
# # exp(rowSums(as.matrix(temp.mat) %*% as.matrix(coef(model)[-1])) * smear) # Calculate adjusted cost predictions
# tempcost <- rowSums(tempcost)^2 + cost_smear[i]
