##IPDfromKM_gen-----
###Generates IPD data using the IPDfromKM function
generate_IPDfromKM <- function(path, 
                               tot.e = NULL, 
                               nrisk,
                               trisk,
                               maxy){
  
  # Arguments:
  # path: the file path to .csv with times, overall survival probabilities
  # tot.e: the total events that occurred
  # nrisk: the no. of patients at risk at given times points
  # trisk: the time intervals at which the nrisk are given
  
  # Returns:
  # IPDfromKM object
  
  #Read in data
  KM   <- read.csv(path)
  
  
  #Preprocess coordinates
  #Intervention group
  preprocess_KM <- IPDfromKM::preprocess(dat   = KM, 
                                         trisk = trisk, 
                                         nrisk = nrisk, 
                                         maxy  = maxy)
  
  #Get IPDs and plot them
  #Overall survival, intervention
  est_surv          <- IPDfromKM::getIPD(prep       = preprocess_KM, 
                                         tot.events = tot.e)
  
  list <- list(IPD = est_surv$IPD,
               unit = est_surv)
  return(list)
  
}