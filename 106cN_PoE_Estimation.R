#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-15
#
# Estimating probability of exceedence along time series using ecdf and
# plotting positions.
# Uses inputs from 102_Threshold_Extract, 104_Event_Extract and 
# 105_Event_Summary
#
# Note that for FUTURE estimation of return periods, make use of 
# the present day data (tSlice) to compute the ecdf.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-15
# Pipeline version ABG 2020-09-07
# NetCDF version ABG 2021-07-06
#
# OUTPUTS: present_returnlevels***.nc
#
#~~~~~~~~~~~~~~~~~~~~~~~
print("running 106N")
if(interactive()){commandArgs <- function(...){c("01","present")}}
##### SETUP #####------------------------------------------------------------
suppressMessages({
  library(extRemes)
  library(lmomco)
})

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

subfold <- paste0("RCM", RCM, suffix, "/")

if(settings$OBSdpe & settings$OBSape){
  stop("OBS probs already exist. stopping 106N.")
}

##### DATA #####--------------------------------------------------------------

suffix_pres  <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- ncoriginal <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix_pres, "_out.nc") 
ncin_pres <- nc_open(ncpres)

# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                            RCM, suffix, ".rds"))

# eventLList length of event L, 
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))

# timewise maxima at each cell for each event (NE x NH)
nc_events <- nc_open(paste0(data_wd,subfold, "eventOBS_",
                          thresh1,"_", ws1,"_RCM", RCM, suffix, ".nc"),
                     write=TRUE)

NE <- sum(ncvar_get(nc_events, "eventNo") > 0)

partable <- readdf(paste0(data_wd,subfold_pres, 
                  "paramtableG_",thresh1, "_RCM", RCM, suffix_pres, ".csv"))

gpa_tracker <- c()
gpa_worst   <- c()
ST0 <- Sys.time()
ST  <- Sys.time()

print("loop start")
for (h in 1:NH) {
  #print(h)
  if ((h < 10) | (h %% 200 == 0)) { # time recording
    print(h)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100 * round((NH - h)/NH ,2)))
    ST <- Sys.time()
    print(ST) 
  }
  
  thr        <- threshMat[h,jT]
  meanInt    <- partable$meanint[h]
  thresholdH <- partable$threshold[h]
  locH       <- partable$loc[h]
  scaleH     <- partable$sca[h]
  shapeH     <- partable$shape[h]

    # Full timeseries for location
  vals <- ncvar_get(ncin_pres, "dmflow",
                         start=c(rn$row[h], rn$col[h], 1),
                         count=c(1, 1, -1))
  
    # Only event summaries at location
  obs_events <- ncvar_get(nc_events, "flow",
                          start=c(h, 1), count=c(1, -1))

    # Compute and store Daily and Annual PoE
  da <- dpeApeComputer(h=h,
                       vals=vals,
                       obs_events=obs_events,
                       ncin=nc_events,
                       pars=partable[h,],
                       thresh_val=thresholdH)

  gpa_tracker[h] <- da[1]
  gpa_worst[h] <- da[2]
}

nc_close(ncin_pres)
nc_close(nc_events)
print(paste("Number of extreme events:", sum(gpa_tracker > 0)))
print("Quantiles of highest return period AT LOCATION:")
print(quantile(gpa_worst, probs=seq(0, 1, by=0.1)))

settings$OBSdpe <- TRUE
settings$OBSape <- TRUE
settings$PoE_Estimation <- "106cN"
write_yaml(settings, settingspath)

print(Sys.time())
print("106N complete.")