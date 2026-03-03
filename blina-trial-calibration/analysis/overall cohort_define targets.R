##### OVERALL COHORT ####

##### BLINA ARM #####

# 01 Initial Setup --------------------------------------------------------

## 01.01 Clean environment ------------------------------------------------
rm(list = ls())
gc()

## 01.02 Load libraries ----------------------------------------------------

# Install packages from CRAN (if required) and load packages
if (!require('pacman')) install.packages('pacman'); library(pacman)
p_load("ggplot2", "dplyr", "scales", "data.table", "flexsurv",
       "patchwork", "matrixStats", "tidyverse", "IPDfromKM")

p_load_gh("Pechli-Lab/SurvdigitizeR")
# library(devtools)
# devtools::install_github("Pechli-Lab/SurvdigitizeR") #unhide if installing SurvdigitizeR package for the first time

## 01.03 Load functions ----------------------------------------------------

source(file = "R/00_digitize_KM_curves_functions.R")
source(file = "R/boot_hr2.R")

# 02 Digitize printed survival curves -------------------------------------

## 02.01 Apply SurvdigitizeR ----------------------------------------------

# Specify image paths and intervention

OS_img_path  <- "figs/AALL1731_overall cohort_OS.png"
PFS_img_path <- "figs/AALL1731_overall cohort_PFS.png"

intervention <- "blina_overall"

# Digitize the plotted OS KM curve
os_curve_dig <- SurvdigitizeR::survival_digitize(
  img_path = OS_img_path,
  num_curves = 2,
  x_start = 0,
  x_end = 4,
  x_increment = 1,
  y_start = 0,
  y_end = 1,
  y_increment = 0.20,
  censoring = F,
  y_text_vertical = T,
  line_censoring = F)

#Visualize
plot(y = os_curve_dig$St[os_curve_dig$curve == 1],
     x = os_curve_dig$time[os_curve_dig$curve == 1], type = "l", ylim = c(0,1))
lines(y = os_curve_dig$St[os_curve_dig$curve == 2],
      x = os_curve_dig$time[os_curve_dig$curve == 2], col = 2)

# Digitize the plotted PFS KM curve
pfs_curve_dig <- SurvdigitizeR::survival_digitize(
  img_path = PFS_img_path,
  num_curves = 2,
  x_start = 0,
  x_end = 4,
  x_increment = 1,
  y_start = 0,
  y_end = 1,
  y_increment = 0.20,
  censoring = F,
  y_text_vertical = T,
  line_censoring = F)

#Visualize
plot(y = pfs_curve_dig$St[pfs_curve_dig$curve == 1],
     x = pfs_curve_dig$time[pfs_curve_dig$curve == 1], type = "l", ylim = c(0,1))
lines(y = pfs_curve_dig$St[pfs_curve_dig$curve == 2],
      x = pfs_curve_dig$time[pfs_curve_dig$curve == 2], col = 2)


## 02.02 Filter digitized data sets ----------------------------------------

os_curve_dig <- os_curve_dig %>% filter(curve == 2) #blina arm

pfs_curve_dig <- pfs_curve_dig %>% filter(curve == 2) #blina arm

# Keep only the columns related to time and survival probabilities

# OS
os_curve_dig <- os_curve_dig %>% 
  select(time, St)

# PFS
pfs_curve_dig <- pfs_curve_dig %>% 
  select(time, St)


## 02.03 Save digitized KM curves ------------------------------------------

#* Save digitized data coming from plots of the digitized KM curves

# write.csv(os_curve_dig,  file = paste0("data-raw/", "OS_digitized","_",intervention,".csv"),  row.names = FALSE)
# write.csv(pfs_curve_dig,  file = paste0("data-raw/", "PFS_digitized", "_", intervention, ".csv"),  row.names = FALSE)


# 03 Reconstruct IPD from digitized KM curves -----------------------------------------

## 03.01 Specify general parameters -----------------------------------------

# Number of people at risk in the OS arm (from the trial publication)

OS_nrisk <- c(718, 628, 465, 249, 50)
OS_trisk <-  c(0,1,2,3,4)             

# Number of people at risk in the PFS arm (from the trial publication)
PFS_nrisk <- c(718, 623, 457, 241, 50)
PFS_trisk <- c(0,1,2,3,4)              

#* Specify if max value for y axis in digitized data is 100 or 1.0
maxy <- round(max(os_curve_dig$St),1)

## 03.02 Generate IPD --------------------------

# Reconstruct IPD for OS
OS_target  <- generate_IPDfromKM(path = paste0("data-raw/", "OS_digitized","_",intervention,".csv"),
                                 nrisk = OS_nrisk,
                                 trisk = OS_trisk,
                                 maxy = maxy)

# Reconstruct IPD for PFS
PFS_target <- generate_IPDfromKM(path = paste0("data-raw/", "PFS_digitized","_",intervention,".csv"),
                                 nrisk = PFS_nrisk,
                                 trisk = PFS_trisk,
                                 maxy = maxy)


#* Join IPD from both OS/PFS into one data frame
df_IPD_all <- dplyr::bind_rows(OS  = OS_target$IPD, 
                               PFS = PFS_target$IPD, 
                               .id = "type")


df_IPD_blina <- df_IPD_all


# write.csv(x = df_IPD_blina,
#           file = paste0("data-raw/df_IPD_blina.csv"),
#           row.names = FALSE)

##### REGENERATE OS/DFS CURVES FOR VALIDATION #####

## 03.03 Fit survival model -----------------------------------

#* Fit a KM model using the reconstructed IPDs from the digitized data. 

df_IPD_blina$time <- df_IPD_blina$time * 365 #convert timescale to days

sfit_reconstructed_KM <- survival::survfit(formula = Surv(time, status) ~ type, 
                                           data    = df_IPD_blina)

times <- seq(0, 1650, 30) #from 0 to 1650 days (4.5 yearS), in months

#* Obtain survival estimates 
rec_KM_surv_curves <- summary(
  object = sfit_reconstructed_KM, 
  scale = 1,
  times = times, 
  extend = TRUE)

# Construct target dataset with reconstructed KM
df_rec_KM_surv_curves <- tibble(
  target      = "Survival",
  type        = stringr::str_sub(rec_KM_surv_curves$strata, 6, -1),
  time        = rec_KM_surv_curves$time,
  time_scale  = "days",
  n.risk      = rec_KM_surv_curves$n.risk,
  n.event     = rec_KM_surv_curves$n.event,
  survival    = rec_KM_surv_curves$surv,
  se          = rec_KM_surv_curves$std.err,
  lower_95_CI = rec_KM_surv_curves$lower,
  upper_95_CI = rec_KM_surv_curves$upper)

# Compute the variance from the survival estimates
df_rec_KM_surv_curves <- df_rec_KM_surv_curves %>% 
  mutate(var = se^2, .after = "se")


write.csv(x = df_rec_KM_surv_curves,
          file = paste0("data-raw/", "KM_reconstructed_", intervention, ".csv"),
          row.names = FALSE)


# ############################################################################################ #

##### SOC ARM #####

# 01 Initial Setup --------------------------------------------------------

## 01.01 Clean environment ------------------------------------------------
rm(list = ls())
gc()

## 01.02 Load libraries ----------------------------------------------------

# # Install packages from CRAN (if required) and load packages
# if (!require('pacman')) install.packages('pacman'); library(pacman)
# p_load("ggplot2", "dplyr", "scales", "IMIS", "data.table", "flexsurv",
#        "patchwork", "matrixStats", "tidyverse", "SurvdigitizeR", "IPDfromKM")
# 
# p_load_gh("Pechli-Lab/SurvdigitizeR")
# library(devtools)
# devtools::install_github("Pechli-Lab/SurvdigitizeR") #unhide if installing SurvdigitizeR package for the first time

## 01.03 Load functions ----------------------------------------------------

source(file = "R/00_digitize_KM_curves_functions.R")

# 02 Digitize printed survival curves -------------------------------------

## 02.01 Apply SurvdigitizeR ----------------------------------------------

# Specify image paths and intervention

OS_img_path  <- "figs/AALL1731_overall cohort_OS.png"
PFS_img_path <- "figs/AALL1731_overall cohort_PFS.png"

intervention <- "soc_overall"

# Digitize the plotted OS KM curve
os_curve_dig <- SurvdigitizeR::survival_digitize(
  img_path = OS_img_path,
  num_curves = 2,
  x_start = 0,
  x_end = 4,
  x_increment = 1,
  y_start = 0,
  y_end = 1,
  y_increment = 0.20,
  censoring = F,
  y_text_vertical = T,
  line_censoring = F)

#Visualize
plot(y = os_curve_dig$St[os_curve_dig$curve == 1],
     x = os_curve_dig$time[os_curve_dig$curve == 1], type = "l", ylim = c(0,1))
lines(y = os_curve_dig$St[os_curve_dig$curve == 2],
      x = os_curve_dig$time[os_curve_dig$curve == 2], col = 2)

# Digitize the plotted PFS KM curve
pfs_curve_dig <- SurvdigitizeR::survival_digitize(
  img_path = PFS_img_path,
  num_curves = 2,
  x_start = 0,
  x_end = 4,
  x_increment = 1,
  y_start = 0,
  y_end = 1,
  y_increment = 0.20,
  censoring = F,
  y_text_vertical = T,
  line_censoring = F)

#Visualize
plot(y = pfs_curve_dig$St[pfs_curve_dig$curve == 1],
     x = pfs_curve_dig$time[pfs_curve_dig$curve == 1], type = "l", ylim = c(0,1))
lines(y = pfs_curve_dig$St[pfs_curve_dig$curve == 2],
      x = pfs_curve_dig$time[pfs_curve_dig$curve == 2], col = 2)

## 02.02 Filter digitized data sets ----------------------------------------

os_curve_dig <- os_curve_dig %>% filter(curve == 1) #soc arm

pfs_curve_dig <- pfs_curve_dig %>% filter(curve == 1) #soc arm

# Keep only the columns related to time and survival probabilities

# OS
os_curve_dig <- os_curve_dig %>% 
  select(time, St)

# PFS
pfs_curve_dig <- pfs_curve_dig %>% 
  select(time, St)


## 02.03 Save digitized KM curves ------------------------------------------

#* Save digitized data coming from plots of the digitized KM curves

# write.csv(os_curve_dig,  file = paste0("data-raw/", "OS_digitized","_",intervention,".csv"),  row.names = FALSE)
# write.csv(pfs_curve_dig,  file = paste0("data-raw/", "PFS_digitized", "_", intervention, ".csv"),  row.names = FALSE)


# 03 Reconstruct IPD from digitized KM curves -----------------------------------------

## 03.01 Specify general parameters -----------------------------------------

# Number of people at risk in the OS arm (from the trial publication)

OS_nrisk <- c(722, 636, 466, 245, 50)
OS_trisk <-  c(0,1,2,3,4)

# Number of people at risk in the PFS arm (from the trial publication)
PFS_nrisk <- c(722, 627, 449, 221, 44)
PFS_trisk <- c(0,1,2,3,4)

#* Specify if max value for y axis in digitized data is 100 or 1.0
maxy <- round(max(os_curve_dig$St),1)

## 03.02 Generate IPD --------------------------

# Reconstruct IPD for OS
OS_target  <- generate_IPDfromKM(path = paste0("data-raw/", "OS_digitized","_",intervention,".csv"),
                                 nrisk = OS_nrisk,
                                 trisk = OS_trisk,
                                 maxy = maxy)

# Reconstruct IPD for PFS
PFS_target <- generate_IPDfromKM(path = paste0("data-raw/", "PFS_digitized","_",intervention,".csv"),
                                 nrisk = PFS_nrisk,
                                 trisk = PFS_trisk,
                                 maxy = maxy)


#* Join IPD from both OS/PFS into one data frame
df_IPD_all <- dplyr::bind_rows(OS  = OS_target$IPD, 
                               PFS = PFS_target$IPD, 
                               .id = "type")

df_IPD_soc <- df_IPD_all

# write.csv(x = df_IPD_soc,
#           file = paste0("data-raw/df_IPD_soc.csv"),
#            row.names = FALSE)


##### REGENERATE OS/DFS CURVES FOR VALIDATION #####

## 03.03 Fit survival model -----------------------------------

#* Fit a KM model using the reconstructed IPDs from the digitized data. 

df_IPD_soc$time <- df_IPD_soc$time * 365 #convert timescale to days

sfit_reconstructed_KM <- survival::survfit(formula = Surv(time, status) ~ type, 
                                           data    = df_IPD_soc)

times <- seq(0, 1650, 30) #from 0 to 1650 days (4.5 yearS), in months

#* Obtain survival estimates 
rec_KM_surv_curves <- summary(
  object = sfit_reconstructed_KM, 
  scale = 1,
  times = times, 
  extend = TRUE)

# Construct target dataset with reconstructed KM
df_rec_KM_surv_curves <- tibble(
  target      = "Survival",
  type        = stringr::str_sub(rec_KM_surv_curves$strata, 6, -1),
  time        = rec_KM_surv_curves$time,
  time_scale  = "years",
  n.risk      = rec_KM_surv_curves$n.risk,
  n.event     = rec_KM_surv_curves$n.event,
  survival    = rec_KM_surv_curves$surv,
  se          = rec_KM_surv_curves$std.err,
  lower_95_CI = rec_KM_surv_curves$lower,
  upper_95_CI = rec_KM_surv_curves$upper)

# Compute the variance from the survival estimates
df_rec_KM_surv_curves <- df_rec_KM_surv_curves %>% 
  mutate(var = se^2, .after = "se")


write.csv(x = df_rec_KM_surv_curves,
          file = paste0("data-raw/", "KM_reconstructed_", intervention, ".csv"),
          row.names = FALSE)




##### ESTIMATION OF A RELATIVE HAZARD FOR DFS #####

## 01 Fit parametric survival models ------------------------------------------

library(darthtools)
library(survHE)

df_IPD_blina <- read.csv(file = "data-raw/df_IPD_blina.csv")
df_IPD_soc <- read.csv(file = "data-raw/df_IPD_soc.csv")

df_IPD_blina2 <- df_IPD_blina %>% filter(type == "PFS") %>% mutate(time = time*365.25) #convert time from years to days
df_IPD_soc2 <- df_IPD_soc %>% filter(type == "PFS") %>% mutate(time = time*365.25)     #convert time from years to days

mods <- c("exponential", "weibull", "weibullPH", "loglogistic","lognormal", "gamma", "gompertz", "rps")

fit.survHE.DFS.blina <- survHE::fit.models( formula = Surv(time, status) ~ 1,
                                  data = df_IPD_blina2, distr = mods)
plot(fit.survHE.DFS.blina)

fit.survHE.DFS.soc <- survHE::fit.models( formula = Surv(time, status) ~ 1,
                                             data = df_IPD_soc2, distr = mods)
plot(fit.survHE.DFS.soc)

aicres.blina <- data.frame(
  model = names(fit.survHE.DFS.blina$models),
  AIC = fit.survHE.DFS.blina$model.fitting$aic,
  BIC = fit.survHE.DFS.blina$model.fitting$bic)

aicres.blina <- aicres.blina %>% arrange(AIC)

aicres.soc <- data.frame(
  model = names(fit.survHE.DFS.soc$models),
  AIC = fit.survHE.DFS.soc$model.fitting$aic,
  BIC = fit.survHE.DFS.soc$model.fitting$bic) 

aicres.soc <- aicres.soc %>% arrange(AIC)

best.DFS.blina <- fit.survHE.DFS.blina$models[[aicres.blina$mod[1]]]
best.DFS.soc <- fit.survHE.DFS.soc$models[[aicres.soc$mod[1]]]

## 02 Bootstrap and extrapolate the hazard ------------------------------------------

times <- seq(0,7300, by = 30)  #20 year extrapolation

# df_HR_PFS <- darthtools::boot_hr(surv_model1 = best.DFS.blina,
#                                  surv_model2 = best.DFS.soc,
#                                  times = times, B = 500)

df_HR_PFS2 <- boot_hr2(surv_model1 = best.DFS.blina,
                                 surv_model2 = best.DFS.soc,
                                 times = times, B = 500)

#View(df_HR_PFS2[["boot.log.hr"]]) 
#* rows: 1 of 500 iterations
#* cols: 244 monthly cycles (i.e., 0-7300 days)

plt <- df_HR_PFS2[["HR"]] %>% ggplot(
  aes(x = time, y = med)) + geom_line(linewidth = 0.8, color = "pink"
  ) + geom_ribbon(
    aes(ymin=lcl, ymax=ucl), colour = NA, fill = "orange2", alpha=0.1) +
  #scale_x_continuous(breaks = c(1:max(times[-1]))) +
  labs(#title = "Relative hazard ratio for DFS",
       #subtitle = paste0("blina = ", aicres.blina$mod[1], " / ", "soc = ", aicres.soc$mod[1]),
       x = "Time in days", y = "Relative hazard ratio (95% CI)") +
  theme_bw() + scale_y_continuous(limits = c(0,1.4)) + geom_hline(yintercept = 1)

plt

# ggsave(plot = plt,
#        filename = paste0("figs/Gupta_HR_PFS_Jan6.png"),
#        height = 4,
#        width = 5)


trial_end <- floor(1650/30)
max_time  <- ceiling(3650/30)

bootHR <- df_HR_PFS2[["boot.log.hr"]]

v_time <- c(trial_end:max_time)

for(i in 1:length(v_time)){
  
bootHR[,v_time[i]] <- bootHR[,v_time[i]] + seq(0,1,length.out = length(v_time))[i]

bootHR[,v_time[i]] <- if_else(bootHR[,v_time[i]] > 1, 1, (bootHR[,v_time[i]]))
}

bootHR[, max_time:ncol(bootHR)] <- 1

bootHR2 = bootHR[,-1]

bootHR2 = cbind(bootHR2, bootHR2[,ncol(bootHR2)])

# write.csv(bootHR2,  file = paste0("data/Gupta_bootHR_PFS_PSA.csv"),  row.names = FALSE)

#Plot again to illustrate the relative hazard's U shape

bootHR_U <- t(bootHR2)

#colnames(bootHR_U) <- times

bootHR_U <- apply(bootHR_U, 1, quantile, c(0.025,0.5,0.975), na.rm = TRUE)
bootHR_U <- as.data.frame(bootHR_U)
bootHR_U <- t(bootHR_U)
bootHR_U <- as.data.frame(bootHR_U)
bootHR_U <- bootHR_U %>% 
  mutate(time = times) %>% 
  rename(med = "50%",
         lcl = "2.5%",
         ucl = "97.5%")


plt_U <- bootHR_U %>% ggplot(
  aes(x = time, y = med)) + geom_line(linewidth = 0.8, color = "pink"
  ) + geom_ribbon(
    aes(ymin=lcl, ymax=ucl), colour = NA, fill = "orange2", alpha=0.1) +
  labs(#title = "Relative hazard ratio for DFS",
    x = "Time in days", y = "Relative hazard ratio (95% CI)") +
  scale_x_continuous(breaks = seq(0, max(bootHR_U$time), by = 1000)) +
  theme_bw()

plt_U

ggsave(plot = plt_U,
       filename = paste0("figs/Gupta_HR_PFS_Ushape.png"),
       height = 4,
       width = 5)

