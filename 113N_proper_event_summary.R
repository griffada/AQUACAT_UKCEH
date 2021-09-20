######
# Adam Griffin, 2020-04-27
#
# Summarising size and time, and extracting daily PoE of extreme events.
# 
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-04-27
# Pipeline version ABG 2020-09-07
#
# Now includes clusters and peakyness.
#
#
# OUTPUTS: eventdf_***.csv: dataframe of each event summarised. One event per row.
#
#####

if(interactive()){commandArgs <- function(...){c("15","future")}}
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

#if(settings$OBSsumm){
#  stop("OBS summary exists. stopping 113N.")
#}

print(paste("Running for threshold", thresh1, "at", ws1, "minimum spread."))

print("Loading in data")
threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                            RCM, suffix, ".rds"))

load(paste0(data_wd, subfold, "eventListsMar_RCM", RCM, suffix, ".RDa"))

param_table <- readdf(paste0(data_wd, subfold, "paramtableG_POT2_RCM",
                            RCM, suffix, ".csv"))

obs_events <- nc_open(paste0(data_wd,subfold, "eventOBSMar_",
                            thresh1,"_", ws1,"_RCM", RCM, suffix, ".nc"))

eventNo <- ncvar_get(obs_events, "eventNo")
NN <- sum(eventNo > 0)

eventSummaryFrame <- data.frame(eventNumber=numeric(NN),
                                eventDay=numeric(NN),
                                eventLength = numeric(NN),
                                area = numeric(NN),
                                peakA = numeric(NN),
                                peakD = numeric(NN),
                                season = character(NN),
                                nclusters = numeric(NN),
                                peakyness = numeric(NN),
                                stringsAsFactors = FALSE)

# (D %/% 90) %% 4 for season 0,1,2,3: DJF, MAM, JJA, SON
season <- function(D) switch(((D %/% 90) %% 4) + 1, "DJF", "MAM", "JJA", "SON")

nclustevent <- function(v, radii=150) {
  # This function estimates the number of "seperate clusters" of event going
  # on during a widespread event. 
  ne0 <- sum(v)
  event1 <- rbind(as.matrix(rn[v,]), matrix(c(
    700, 0, 700000, 1000000
  ), nrow = 1))
  d <- dist(event1[, 1:2])
  d2 <- dist(event1[1:ne0, 1:2])
  if ((max(d2) > radii) & (ne0 > 1)) {
    hh <- hclust(d, method = "single")
    ch <- max(cutree(hh, h = 0.8 * sort(hh$height, decreasing = T)[2])) - 1
    if(sort(hh$height,decreasing=T)[2] < 0.2*max(hh$height)){
      ch <- 1
    }
    ch <- min(ch, 4)
  } else{
    ch <- 1
  }
  ch
}

peakyfun <- function(v, thr) {
  mn <- min(v)
  if (all(is.na(v))| mn > thr) {
    return(0)
  }else{
    mn50 <- mean(c(mn, thr))
    return(sum(v < mn50) / sum(v < thr))
  }
}

rn <- as.data.frame(rn)
ST <- Sys.time()
ST0 <- Sys.time()
print(NN)
for (i in 1:NN) {
  if ((i < 10) | (i %% 100 == 0)) { # time recording
    print(i)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((NN-i)/NN ,2)))
    print(paste("Time remaining", round((NN-i)/i * I0,2)))
    print(paste("Since last readout:", round(I,2)))
    ST <- Sys.time()
    print(ST) 
  }
  try({
  eventSummaryFrame[i,] <- eventSummaryLine(obs_events, i, threshMat[,jT])
  })
}

print("saving OBS Summary")
readr::write_csv(eventSummaryFrame, paste0(data_wd,subfold, "eventSumm_OBS_",
                                  thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))

settings$OBSsumm <- TRUE
write_yaml(settings, settingspath)

print("113N done.")