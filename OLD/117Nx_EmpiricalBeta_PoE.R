#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2021-05-04
#
# The Empirical Copula functions for simulating new events from given data.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-15
# Pipeline version ABG 2020-09-07
# Version 2 using empirical beta copulas
# netcdf version ABG 2021-07-07
#
# OUTPUTS: NewEventPresentEC_***.nc: data table of events, one event per row.
#
#~~~~~~~~~~~~~~~~~~~~~~~

if(interactive()){commandArgs <- function(...){c("13","present")}}
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

# if(settings$EC2ape & settings$EC2dpe){
#  print("EC probabilities exist. Stopping 117N.") 
# }#else{

##### SETUP #####---------------------------------------------------
suppressPackageStartupMessages({
library(extRemes)
library(dplyr)
library(fitdistrplus)
library(parallel)
library(foreach)
library(ilaprosUtils)
library(lmomco)
})

readdf <- function(...){as.data.frame(data.table::fread(...))}
### DATA --------------------------

print(paste("Running for threshold", thresh1, "at", ws1, "minimum spread."))
print(paste("RCM", RCM, "period", period))

suffix_pres <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix_pres, "_out.nc")
ncin_pres <- nc_open(ncpres)

threshMat <- readRDS(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"))

obs_events  <- nc_open(paste0(data_wd,subfold,"eventOBS_",thresh1, "_", ws1,
                             "_RCM", RCM, suffix, ".nc"))
NE <- obs_events$dim$event$size
NH <- nrow(rn)

partable    <- readdf(paste0(data_wd,subfold,"paramtableH_",thresh1,
                              "_RCM", RCM, suffix, ".csv"))
partable_pres <- readdf(paste0(data_wd,subfold_pres, 
                  "paramtableH_", thresh1, "_RCM", RCM, suffix_pres, ".csv"))

ec_events <- nc_open(paste0(data_wd, subfold, "eventECD_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".nc"), write=T)

# ec_events_old <- nc_open(paste0(data_wd, subfold, "eventECB_",
#                         thresh1,"_", ws1, "_RCM", RCM, suffix, ".nc"))

##### DPE/APE CALCULATION #####---------------------------------------------

gpa_tracker <- c()
gpa_worst <- c()

EN <- sum(ncvar_get(ec_events, "eventNo")>0)
ST0 <- Sys.time()
ST <- Sys.time()
print("loop start")
for(h in 1:20){
  #print(h)
  if((h < 10) | (h %% 1000 == 0)){ # time recording
    print(h)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((NH-h)/NH ,2)))
    print(paste("Time remaining", round((NH-h)/h * I0,2)))
    print(paste("Since last readout:", round(I,2)))
    ST <- Sys.time()
    print(ST) 
  }
  
  thr <- threshMat[h,jT]
  meanInt <- partable$meanint[h]
  thresholdH <- partable$threshold[h]
  locH <- partable$loc[h]
  scaleH <- partable$sca[h]
  shapeH <- partable$shape[h]
  
  vals <- ncvar_get(ncin_pres, "dmflow",
                        start=c(rn$row[h], rn$col[h], 1),
                        count=c(1, 1, -1))
  
  obs_flow <- ncvar_get(ec_events, "flow",
                          start=c(h,1), count=c(1,EN))

  da <- dpeApeComputer(h,
                 vals,
                 obs_flow,
                 ncin=ec_events,
                 pars=vec2par(c(locH,scaleH,shapeH), type='gpa'),
                 thresh_val=thresholdH)

  gpa_tracker[h] <- da[1]
  gpa_worst[h] <- da[2]
  
}

print(quantile(gpa_tracker, probs=seq(0,1,by=0.1)))
print(quantile(gpa_worst, probs=seq(0,1,by=0.1)))

#### SAVE OUTPUTS ####--------------------------------------------------------
print("Saving outputs")
nc_close(obs_events)
nc_close(ec_events)
#nc_close(ncin_pres)

settings$EC2dpe <- TRUE
settings$EC2ape <- TRUE
write_yaml(settings, settingspath)
#}
print(Sys.time())
print("117N done.")