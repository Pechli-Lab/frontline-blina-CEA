################################################################################
#                                   MICROSIMULATION START                      #
################################################################################
p <- Sys.time() # start time
set.seed(2)

# Set up parallel
Cores <- detectCores() -3
cl <- makeCluster(Cores)
registerDoParallel(cl)
registerDoRNG(seed = 123)

################################################################
## OPEN LOOP (1) - number of simulations
################################################################
for(k in 1:n_sim){

  # Change seed for each simulation
  set.seed(SEED + k)

  gc() # clears the memory

  print(paste(format(Sys.time()), "-- k =", k, sep = " ")) # print current time and simulation number

  results.sim.par <- 0                                                        # Initialize storage for parallel loop
  inds <- split(seq_len(n_i_i), sort(rep_len(seq_len(Cores), n_i_i)))         # Split number of individuals into groups for parallel loop

  ################################################################
  ## OPEN LOOP (2) - number of individuals (parallel)
  ################################################################
  results.sim.par <- foreach(
    l = seq_along(inds),
    .combine = rbind,
    .packages = c(
      "mvtnorm", "corpcor", "rlang", "survival", "flexsurv", "MASS", "msm", "gems",
      "dplyr", "DAAG", "cmprsk", "reshape", "splines", "survival", "doParallel", "data.table",
      "plyr", "ggplot2", "purrr", "survminer", "gridExtra", "foreach", "flexsurvcure", "etm"
    )
  ) %dopar% {

    print(l)
    n_ind <- length(inds[[l]])        # Number of individuals in this group
    X     <- b.dat.01[inds[[l]],]     # Assign baseline characteristics
    X_c   <- bl.dat.c[inds[[l]],]     # Assign baseline costs characteristics
    v_M_Init <- v_M_Init_nii[inds[[l]]] # Initial health states for individuals in this group


    # INITIALIZE  
    v_init <- c("PRE", "REL1", "REL2", "BMT1", "BMT1_REL") # vector of non-absorbing states

    for (i in 1:length(v_init)) {
      assign(paste0("v_", v_init[i], "_Init"), as.numeric(v_M_Init == v_init[i]) * step)
    }

    # m_M: health state for each patient at each cycle
    m_M <- m_C <- m_E <- matrix(
      nrow = n_ind, 
      ncol = length(times),
      dimnames = list(
        paste("ind", 1:n_ind, sep = " "),
        times
      )
    )


    # CREATE ATTRIBUTE MATRICES 'm_(state)' to track time in states that aren't absorbing
    m_PRE <- m_REL1 <- m_REL2 <- m_BMT1 <- m_BMT1_REL <- matrix(0,
      nrow = n_ind,
      ncol = length(times),
      dimnames = list(
        paste("ind", 1:n_ind, sep = " "),
        paste("cycle", times, sep = " ")
      )
    )

    # Initialize health-state matrices: m_
    m_M[,1] <- v_M_Init          # Health states of individuals at cycle 0
    for (i in 1:length(v_init)){
      temp <- '[<-'(eval(as.name(paste0("m_", v_init[i]))), i = , j = 1, get(paste0("v_", v_init[i], "_Init")))
      assign(paste0("m_", v_init[i]), temp)
    }

    # m_int:    tracks interval count for each cycle
    # m_age:    updates individuals age for cost models when individual changes state
    # m_dxyear: updates individuals dxyear for cost models when individual changes state 
    m_int <- m_age <- m_dxyear <- matrix(0,
      nrow = n_ind, 
      ncol = length(times),
      dimnames = list(
        paste("ind", 1:n_ind, sep = " "), # rows
        paste("cycle", times, sep = " ")  # cols
      )
    )

    # Initialize m_int, m_age and m_dxyear
    m_int[, 1] <- 1
    m_age[, 1] <- X_c$age_Dx
    m_dxyear[, 1] <- X_c$dxyear

    # Index variable for cycles
    j <- 0

    # Calculate initial costs and effects at cycle 0
    m_C[, j + 1] <- pred.cost(M_t = m_M, df = X_c, interval = m_int, age = m_age, year = m_dxyear)
    m_E[, j + 1] <- Effs(M_t = m_M[, j + 1], k = k)


    ################################################################
    ## OPEN LOOP (3) - time (number of cycles)
    ################################################################
    for (t in times[-length(times)]) {
      j = j + 1 # update cycle index (ex. first cycle (cycle 0) is at index 1)

      # TRANSITION PROBABILITIES for n_sim
      for (trans in 1:sum(!is.na(m_t))) {
        transition <- which(m_t == trans, arr.ind = T)   # transition is the index in the transition matrix for transition = i

        s.from <- transition[1] # from state = row of the transition matrix, or remainder when divided by # states
        s.to <- transition[2] # to state = column of the transition matrix, or rounded up division by # states

        # List of covariate names in model for transition 'trans'
        cov.trans <- as.character(l_MSM_est[[trans]]$covariate)

        # Calculate transition probabilities
        p.trans <- t(model.dist.f(
          dist.v = l_MSM_est[[trans]]$model[[1]], # model distribution (exp, weibull...)
          d.data = norm.mat[[trans]][k, ], # multivariate normal estimates
          dat.x = as.matrix(X[, cov.trans[cov.trans %chin% names(X)]]), # matrix of baseline covariates for covariates used in model, ordered by how they were used in the model
          t = get(paste0("m_", v_n[s.from]))[, j], # time spent in states
          trans = trans # transition number
        ))
        assign(paste0("p.", v_n[s.from], "_", v_n[s.to]), p.trans)
      }

      # Get transition probabilities based on current health state and time since illness
      m_Probs <- Probs(m_M[, j], X = X, blin_2ndline = is_TRUE(blin_2ndline), blin_1stline = is_TRUE(blin_1stline))

      # Sample the  health state at next time interval based on transition probabilities p
      m_M[, j + 1] <- darthtools::samplev(m_Probs, m = 1)

      # Update time in states, time to relapse, previous relapses, interval, age and dxyear
      m_PRE      [,j + 1] <- ifelse(m_M[,j + 1] == "PRE",      m_PRE      [,j] + step, 0)
      m_REL1     [,j + 1] <- ifelse(m_M[,j + 1] == "REL1",     m_REL1     [,j] + step, 0)
      m_REL2     [,j + 1] <- ifelse(m_M[,j + 1] == "REL2",     m_REL2     [,j] + step, 0)
      m_BMT1     [,j + 1] <- ifelse(m_M[,j + 1] == "BMT1",     m_BMT1     [,j] + step, 0)
      m_BMT1_REL [,j + 1] <- ifelse(m_M[,j + 1] == "BMT1_REL", m_BMT1_REL [,j] + step, 0)

      X$REL1.t    <- ifelse(m_M[,j+1] == "REL1" & m_M[,j] != "REL1", t + step, X$REL1.t)

      # Whether you relapsed or not as captured in prev_REL1 is agnostic to previous BMT
      X$prev_REL1 <- ifelse(m_M[, j + 1] == "REL1" | m_M[, j + 1] == "BMT1_REL", 1, X$prev_REL1)
      X$prev_BMT  <- ifelse(m_M[, j + 1] == "BMT1" | m_M[, j + 1] == "BMT1_REL", 1, X$prev_BMT)

      # Update interval, age and dxyear
      m_int[, j+1]  <- ifelse(m_M[,j] == m_M[,j+1], m_int[,j] + 1, 1)
      m_age[, j+1]  <- ifelse(m_M[,j] == m_M[,j+1], m_age[,j], (m_age[,j] + m_int[,j] * step/365.25))
      m_dxyear[,j+1] <- ifelse(m_M[,j] == m_M[,j+1], m_dxyear[,j], m_dxyear[,j] + m_int[,j] * step/365.25)

      # Calculate costs and effects
      m_C[, j + 1] <- pred.cost(M_t = m_M, df = X_c, interval = m_int, age = m_age, year = m_dxyear)
      m_E[, j + 1] <- Effs(M_t = m_M[, j + 1], k = k)
    } # CLOSE LOOP (cycles)


    # COSTS
    # Capping costs - look at diff to see when it's positive (costs increase)

    for (i in c("PRE", "REL1", "REL2", "BMT1")) { #BMT1_REL didn't explode
      l_time <- apply(m_M, 1, function(x) which(x == i))

      if (length(l_time) != 0) {

        indx <- which(sapply(l_time, function(x) length(x) != 0))   # get individuals who visit state x
        l_time <- l_time[indx]

        for (j in 1:length(indx)) {
          dif <- diff(m_C[indx[j], l_time[[j]]], differences = 1)
          dif_sign <- dif/abs(dif)

          coln <- names(which(dif_sign == 1)[1]) 

          y_min <- which(colnames(m_C)==coln)
          y_max <- which(colnames(m_C)==tail(names(m_C[indx[j], l_time[[j]]]),1))

          if (length(coln != 0)) { # not NULL
            if (!is.na(coln)) {
              m_C[indx[j], y_min:y_max] <- m_C[indx[j], which(colnames(m_C)==coln)-1]
            }
          }
        }
      }
    }

    # Add costs of death for interval prior to death
    # Get first time for "DEAD" in each row, if none then NA, change costs for cycle prior to death to DEAD costs 
    j_dead <- apply(m_M, 1, function(x) {
      ifelse(
        length(as.numeric(names(head(x[x == "DEAD"], 1)))) == 0,
        NA,
        as.numeric(names(head(x[x == "DEAD"], 1)))
      )
    })
    x <- which(!is.na(j_dead))
    y <- j_dead[which(!is.na(j_dead))] / step # column prior to death

    m_C[cbind(x, y)] <- c_dead[k] # For each individual that dies, change costs of interval prior to death to DEAD costs 

    m_C <- round(m_C, 2)      # round to 2 decimal places


    tc <- m_C %*% v_dwc  # total (discounted) cost per individual
    te <- m_E %*% v_dwe  # total (discounted) QALYs per individual
    tc_hat <- mean(tc, na.rm = TRUE)   # average (discounted) cost
    te_hat <- mean(te, na.rm = TRUE)   # average (discounted) QALY

    return(list(events = m_M, costs = m_C, effects = m_E, tc = tc, te = te, tc_hat = tc_hat, te_hat = te_hat, tREL = X$prev_REL1, tBMT= X$prev_BMT))
  } # CLOSE LOOP (individuals)
  
  
  # Save resutls from all cores back together
  events <- costs <- effects <- tc <- te <- tc_hat <- te_hat <-t_REL <- t_BMT <- list()

  # Extract the results
  for (i in 1:dim(results.sim.par)[1]){
    events  <- `if`(i == 1, results.sim.par[[i,1]], rbind(events,  results.sim.par[[i,1]]))
    costs   <- `if`(i == 1, results.sim.par[[i,2]], rbind(costs,   results.sim.par[[i,2]]))
    effects <- `if`(i == 1, results.sim.par[[i,3]], rbind(effects, results.sim.par[[i,3]]))
    tc      <- `if`(i == 1, results.sim.par[[i,4]], rbind(tc,      results.sim.par[[i,4]]))
    te      <- `if`(i == 1, results.sim.par[[i,5]], rbind(te,      results.sim.par[[i,5]]))
    tc_hat  <- `if`(i == 1, results.sim.par[[i,6]], rbind(tc_hat,  results.sim.par[[i,6]]))
    te_hat  <- `if`(i == 1, results.sim.par[[i,7]], rbind(te_hat,  results.sim.par[[i,7]]))
    t_REL  <- `if`(i == 1, results.sim.par[[i,8]], rbind(t_REL,  results.sim.par[[i,8]]))
    t_BMT  <- `if`(i == 1, results.sim.par[[i,9]], rbind(t_BMT,  results.sim.par[[i,9]]))
  }

  # Save the various results to associated .rds files
  saveRDS(as.data.frame(events),  paste0(output.path, "data.",   k, ".rds"))
  saveRDS(as.data.frame(costs),   paste0(output.path, "cost.",   k, ".rds"))
  saveRDS(as.data.frame(effects), paste0(output.path, "effect.", k, ".rds"))
  saveRDS(as.data.frame(tc)     , paste0(output.path, "tc.",     k, ".rds"))
  saveRDS(as.data.frame(te)     , paste0(output.path, "te.",     k, ".rds"))
  saveRDS(as.data.frame(tc_hat) , paste0(output.path, "tc_hat.", k, ".rds"))
  saveRDS(as.data.frame(te_hat) , paste0(output.path, "te_hat.", k, ".rds"))
  saveRDS(as.data.frame(t_REL) , paste0(output.path, "tREL.", k, ".rds"))
  saveRDS(as.data.frame(t_BMT) , paste0(output.path, "tBMT.", k, ".rds"))
  rm(list = c("events", "costs", "effects", "tc", "te", "tc_hat", "te_hat", "t_REL", "t_BMT"))
} # CLOSE LOOP (n_sim)

comp.time = Sys.time() - p
print(n_i_i)
print(n_sim)
print(comp.time)
# Print the time taken (in mins) to a txt file
write(paste0("Time taken for ", n_sim, " simulations with ", n_i_i, " individuals: ", round(comp.time, 2), " hours"), file = paste0(output.path, "simulation_time.txt"))
print(Cores)
stopCluster(cl)
