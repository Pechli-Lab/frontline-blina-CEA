#' Predict Costs
#'
#' @param M_t       matrix of individuals states for that cycle
#' @param df        dataframe of individuals covariate values
#' @param interval  matrix of interval counts for each individual
#' @param year      matrix of dxyear for that state
#' @param age       matrix of ages when individuals entered their current state
#'
#' @return vector of costs for individuals during cycle `j`
#'
pred.cost <- function(M_t, df, interval, year, age){


  newdata <- df
  newdata$interval <- interval[ , j+1]
  newdata$`I(interval^2)` <- newdata$interval^2

  newdata$age_Dx <- age[, j+1]
  newdata$dxyear <- year[, j+1]


  cost <- if_else(M_t[,j+1] == "PRE",  Costcycle("PRE",  newdata = newdata),
          if_else(M_t[,j+1] == "REL1", Costcycle("REL1", newdata = newdata),
          if_else(M_t[,j+1] == "REL2", Costcycle("REL2", newdata = newdata),
          if_else(M_t[,j+1] == "BMT1", Costcycle("BMT1", newdata = newdata),
          if_else(M_t[,j+1] == "BMT1_REL", Costcycle("BMT1_REL", newdata = newdata),
          if_else(M_t[,j+1] == "DEAD", 0, NA_real_))))))

  return(cost)
}


#' Costcycle
#'
#' @param state which state to calculate costs for
#' @param newdata individual covariates
#' @param inds.   list of individuals index's that were split for parallelization (see `model_run.R` where created)
#' @param k.      current iteration of the psa simulation
#' @param l.      position in the list of individuals index's (l = 1, 2, 3, ...)
#'
#' @return vector of costs for individuals during cycle `j`
#'
#' @note Global variables used (not passed as arguments):
#'   - j: current time/cycle index
#'   - v_n: vector of state names
#'   - l_costcovs_ns: 
#'   - nspline: 
#'   - l_costcovs: 
#'   - norm.mat.cost: 
#'   - l_cost_reff: 
#'   - RE: random effects 
#'   - cost_smear: 
#'   - costs.trans: 
#'
Costcycle <- function(state, newdata, inds. = inds, k. = k, l. = l){

  i <- which(v_n == state)

  if (length(l_costcovs_ns[[i]]) >= 1) {
    newdata[, l_costcovs_ns[[i]]] <- nspline[[i]][newdata$interval, ]
  }

  m.newdata <- as.matrix(newdata[, l_costcovs[[i]]])

  temp.coeff <- t(as.matrix(norm.mat.cost[[i]][k., ]))

  # Add random effects
  if (l_cost_reff[i] != 0) {
    tempnames <- names(norm.mat.cost[[i]][k., ])[names(norm.mat.cost[[i]][k., ]) %chin% colnames(RE[[i]][inds[[l.]] * k., ])]
    temp.coeff <- temp.coeff[rep(1, length(inds.[[l.]])), ]
    temp.coeff[, tempnames] <- temp.coeff[, tempnames] + RE[[i]][inds.[[l.]] * k., ]

    tempcost <- m.newdata * temp.coeff
  } else {
    tempcost <- sweep(m.newdata, 2, temp.coeff, "*")
  }

  # SMEARING (heteroskedastic)
  if (costs.trans == "log") {
    tempcost <- exp(rowSums(tempcost) + 1/2 * cost_smear[[i]])
  } else if (costs.trans == "sqrt") {
    tempcost <- (rowSums(tempcost) + 1/2 * cost_smear[[i]])^2
  }

  return(tempcost)
}
