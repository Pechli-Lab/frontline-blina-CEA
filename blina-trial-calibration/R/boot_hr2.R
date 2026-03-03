boot_hr2 <- function (surv_model1 = NULL, surv_model2 = NULL, rx1 = NULL, 
          rx2 = NULL, rx = F, surv_model_rx = NULL, times, B = 100) 
{
  if (rx == F) {
    boot.haz1 = boot.haz(x = surv_model1, B = B, t = times)
    boot.haz1[boot.haz1 == 0] = 0.01
    boot.haz1 <- boot.haz1[1, , ]
    boot.log.haz1 <- log(boot.haz1)
    boot.haz2 = boot.haz(x = surv_model2, B = B, t = times)
    boot.haz2[boot.haz2 == 0] = 0.01
    boot.haz2 <- boot.haz2[1, , ]
    boot.log.haz2 <- log(boot.haz2)
    boot.log.hr <- boot.log.haz1 - boot.log.haz2
    HR <- apply(boot.log.hr, 2, function(x) exp(quantile(x, 
                                                         probs = c(0.025, 0.5, 0.975))))
    HR <- as.data.frame(t(HR))
    HR$times <- times
    colnames(HR) <- c("lcl", "med", "ucl", "time")
    return(list(HR = HR,
                boot.log.hr = exp(boot.log.hr)))
  }
  else {
    mod.rx <- surv_model_rx
    distt <- mod.rx$call$dist
    t <- times
    hr.est <- summary(mod.rx, t = t, type = "hazard", newdata = data.frame(rx = rx1), 
                      ci = FALSE)[[1]][, "est"]/summary(mod.rx, t = t, 
                                                        type = "hazard", newdata = data.frame(rx = rx2), 
                                                        ci = FALSE)[[1]][, "est"]
    pars <- normboot.flexsurvreg(mod.rx, B = B, newdata = data.frame(rx = c(rx2, 
                                                                            rx1)))
    hr <- matrix(nrow = B, ncol = length(t))
    for (i in seq_along(t)) {
      haz.rx1.rep <- do.call(mod.rx$dfns$h, c(list(t[i]), 
                                              as.data.frame(pars[[2]])))
      haz.rx2.rep <- do.call(mod.rx$dfns$h, c(list(t[i]), 
                                              as.data.frame(pars[[1]])))
      hr[, i] <- haz.rx1.rep/haz.rx2.rep
    }
    hr <- apply(hr, 2, quantile, c(0.025, 0.975))
    HR <- data.frame(lcl = hr[1, ], med = hr.est, ucl = hr[2, 
    ], time = times)
    plot(t, hr.est, type = "l", ylim = c(0, 2), col = "red", 
         xlab = "time", ylab = paste0("Hazard ratio", rx1, 
                                      "/ ", rx2), lwd = 2)
    lines(t, hr[1, ], col = "red", lwd = 1, lty = 2)
    lines(t, hr[2, ], col = "red", lwd = 1, lty = 2)
    return(HR)
  }
}
