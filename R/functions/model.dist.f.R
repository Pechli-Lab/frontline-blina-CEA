#' Get transition probabilities
#'
#' @param dist.v character: model distribution
#' @param d.data vector: model estimates (flexsurvreg(Surv(t, event) ~ x1 + x2 + x3, dist = "...")$res[,1])
#' @param dat.x  matrix: individuals covariates
#' @param t      vector: time spent in state
#' @param model  character: name of model
#'
#' @return transition probability for this cycle, for each individual
#'
#' @note Global variables used:
#'   - step: Time step size used in calculating transition probabilities
#'   - l_MSM_knots: List of knots for spline models, accessed by transition number in spline models
#'   - n_ind and n_i: Number of individuals, used in spline models
#' 
model.dist.f <- function(dist.v = NA, d.data = NA, dat.x = NA, t = NA, trans = NA){

  # Get cumulative survival probability for each individual
  ############################ Gamma ############################
  if(dist.v == "Gamma"){

    pars <- d.data[1:2]
    beta <- if(length(d.data) > 2){d.data[3:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]

    for (i in 1:length(pars)) {
      if(i == 2){

        # beta affects the rate
        if(length(beta) == 0){beta.raw = 0 # if no covariates
        }else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms
      }
      else{print("gamma")}
    }
    p.res   <- pgamma(t,        shape = pars[1], rate = pred, lower.tail = F)
    p.res.1 <- pgamma(t - step, shape = pars[1], rate = pred, lower.tail = F)


    ## PP : Perhaps consider changing to
    #H.res   <- Hgamma(t,        shape = pars[1], rate = pred)
    #H.res.1 <- Hgamma(t - step, shape = pars[1], rate = pred)

    # And then at the end of the function do:
    #H_diff <- Hres - Hres1
    #tp <- 1- exp(-H_diff)

  ############################ Gamma Cure ############################
  } else if (dist.v == "Gamma.c"){

    pars <- d.data[1:3]
    beta <- if(length(d.data) > 3){d.data[4:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]

    for (i in 1:length(pars)){
      if (i == 3){

        # beta affects the rate
        if(length(beta) == 0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)           # fit$dlist$inv.transforms
      }
      else{print("gamma")}
    }
    p.res   <- pmixsurv(pgamma, t,        theta = pars[1], shape = pars[2], rate = pred, lower.tail = F)
    p.res.1 <- pmixsurv(pgamma, t - step, theta = pars[1], shape = pars[2], rate = pred, lower.tail = F)

  ############################ Exponential ############################
  } else if (dist.v == "Exponential") {

    pars <- d.data[1]
    beta <- if(length(d.data) > 1){d.data[2:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]

    for (i in 1:length(pars)) {

      # beta affects the rate
      pars <- d.data[1]
      beta <- if(length(d.data) > 1){d.data[2:length(d.data)]} else{NA}
      beta <- beta[!is.na(beta)]

      if(length(beta) ==0){beta.raw = 0}
      else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

      pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
      pred     <- exp(pred.raw)              # fit$dlist$inv.transforms
    }
    p.res   <- pexp(t,        rate = pred, lower.tail = F)
    p.res.1 <- pexp(t - step, rate = pred, lower.tail = F)

  ############################ Exponential Cure ############################
  } else if (dist.v == "Exponential.c") {

    pars <- d.data[1:2]
    beta <- if(length(d.data) > 2){d.data[3:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]

    for (i in 1:length(pars)) {
      if (i == 2){
        print("rate") # Beta affects the rate

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms
      }
      else{print("exponential cure")}
    }
    p.res   <- pmixsurv(pexp, t,        theta = pars[1], rate = pred, lower.tail = F)
    p.res.1 <- pmixsurv(pexp, t - step, theta = pars[1], rate = pred, lower.tail = F)

  ############################ Weibull ############################
  } else if (dist.v == "Weibull") {

    pars <- d.data[1:2]
    beta <- if(length(d.data) > 2){d.data[3:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]

    for (i in 1:length(pars)) {
      if(i == 2) {
        # Beta affects the scale

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms

      }
      else{print("weibull")}
    }
    p.res   <- pweibull(t,        shape = pars[1], scale = pred, lower.tail = F)
    p.res.1 <- pweibull(t - step, shape = pars[1], scale = pred, lower.tail = F)

  ############################ Weibull Cure ############################
  } else if (dist.v == "Weibull.c") {

    pars <- d.data[1:3]
    beta <- if(length(d.data) > 3){d.data[4:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]

    for (i in 1:length(pars)) {
      if(i == 3) {
        # Beta affects the scale

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms

      }
      else{print("weibull")}
    }
    p.res   <- pmixsurv(pweibull, t,        theta = pars[1], shape = pars[2], scale = pred, lower.tail = F)
    p.res.1 <- pmixsurv(pweibull, t - step, theta = pars[1], shape = pars[2], scale = pred, lower.tail = F)

  ############################ LogNormal ############################
  } else if (dist.v == "LogNormal") {

    pars <- d.data[1:2]
    beta <- if(length(d.data) > 2){d.data[3:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]

    for (i in 1:length(pars)) {
      if(i == 1) {
        #meanlog

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- pars[i] + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- pred.raw              # fit$dlist$inv.transforms
      }
      else{print("lognormal")}
    }
    p.res   <- plnorm(t,        meanlog = pred, sdlog = pars[2], lower.tail = F)
    p.res.1 <- plnorm(t - step, meanlog = pred, sdlog = pars[2], lower.tail = F)

  ############################ LogNormal Cure ############################
  } else if (dist.v == "LogNormal.c") {

    pars <- d.data[1:3]
    beta <- if(length(d.data) > 3){d.data[4:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]

    for (i in 1:length(pars)) {
      if(i == 2) {
        #meanlog

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- pars[i] + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- pred.raw              # fit$dlist$inv.transforms
        }
      else{print("lognormal cure")}
    }

    p.res   <- pmixsurv(plnorm, t       , theta = pars[1], meanlog = pred, sdlog = pars[3], lower.tail = F)
    p.res.1 <- pmixsurv(plnorm, t - step, theta = pars[1], meanlog = pred, sdlog = pars[3], lower.tail = F)
  ############################ Gompertz ############################
  } else if (dist.v == "Gompertz") {

    pars <- d.data[1:2]
    beta <- if(length(d.data) > 2){d.data[3:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]

    for (i in 1:length(pars)) {
      if(i == 2) {

        #rate
        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms

      }
      else{print("gompertz")}
    }
    p.res   <- pgompertz(t,        shape = pars[1], rate = pred, lower.tail = F)
    p.res.1 <- pgompertz(t - step, shape = pars[1], rate = pred, lower.tail = F)

  ############################ Gompertz Cure ############################
  } else if (dist.v == "Gompertz.c") {

    pars <- d.data[1:3]
    beta <- if(length(d.data) > 3){d.data[4:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]

    for (i in 1:length(pars)) {
      if(i == 3) {
        #rate

        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms
      }
      else{print("gompertz cure")}
    }
    p.res   <- pmixsurv(pgompertz, t,        theta = pars[1], shape = pars[2], rate = pred, lower.tail = F)
    p.res.1 <- pmixsurv(pgompertz, t - step, theta = pars[1], shape = pars[2], rate = pred, lower.tail = F)

  ############################ LogLogistic ############################
  } else if (dist.v == "LogLogistic") {

    pars <- d.data[1:2]
    beta <- if(length(d.data) > 2){d.data[3:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]


    for (i in 1:length(pars)) {
      if(i == 2){

        #Scale
        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms
      }
      else{print("loglogistic")}
    }
    p.res   <- pllogis(t,        shape = pars[1], scale = pred, lower.tail = F)
    p.res.1 <- pllogis(t - step, shape = pars[1], scale = pred, lower.tail = F)

  ############################ LogLogistic Cure ############################
  } else if (dist.v == "LogLogistic.c") {

    pars <- d.data[1:3]
    beta <- if(length(d.data) > 3){d.data[4:length(d.data)]} else{NA}
    beta <- beta[!is.na(beta)]

    for (i in 1:length(pars)){
      if(i == 3){

        #Scale
        if(length(beta) ==0){beta.raw = 0}
        else{beta.raw <- t(as.matrix(beta)) %*% t(dat.x)}

        pred.raw <- log(pars[i]) + beta.raw    # fit$dlist$transforms + beta.raw
        pred     <- exp(pred.raw)              # fit$dlist$inv.transforms
      }
      else{print("loglogistic cure")}
    }

    p.res   <- pmixsurv(pllogis, t,        theta = pars[1], shape = pars[2], scale = pred, lower.tail = F)
    p.res.1 <- pmixsurv(pllogis, t - step, theta = pars[1], shape = pars[2], scale = pred, lower.tail = F)

  ############################ Spline ############################
  } else if (dist.v == "Spline") {
    beta  <- d.data[!is.na(d.data)]  # remove extra NA value
    p.tot <- length(beta)            # total no. of pars and betas

    pars <- d.data[1:(p.tot - (dim(dat.x)[2]))]                   # parameter estimates (gamma0, gamma1, ...)

    if (p.tot > length(pars)) {
      beta <- beta[(length(pars) + 1) : p.tot]                    # if covariates are included, extract coefficient estimates
    }else {beta <- vector()}                                      # otherwise beta is an empty vector

    pred <- matrix(NA, nrow = nrow(dat.x), ncol = length(pars))   # matrix to store pred values

    for (i in 1:length(pars)) {
      if(i == 1) { # Beta affects the first param
        beta.raw <- t(as.matrix(beta)) %*% t(dat.x)

        pred.raw <- pars[i] + beta.raw
        pred[, i] <- pred.raw
      }
      if(i > 1) {
        #pred[, i] <- rep(pars[i], n_i)   # use if not parallel
        pred[, i] <- rep(pars[i], n_ind) # use if in parallel
      }
    }

    # return results for spline
      p.res   <- t(psurvspline(t,        gamma = pred,  knots = l_MSM_knots[[trans]]$knots, lower.tail = F))
      p.res.1 <- t(psurvspline(t - step, gamma = pred,  knots = l_MSM_knots[[trans]]$knots, lower.tail = F))

  ############################ No Distribution ############################
  } else {
    print("no distribution found")
    p.res <- NA
    p.res.1 <- NA
  }

  return(1 - p.res/p.res.1)   # return transition probability at timepoint
}
