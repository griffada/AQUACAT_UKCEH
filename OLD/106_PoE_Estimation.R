#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-15
#
# Estimating probability of exceedence along time series using ecdf and
# plotting positions.
# Uses inputs from 102_Threshold_Extract, 104_Event_Extract and 105_Event_Summary
#
# Note that for FUTURE estimation of return periods, make use of the present day
# data (tSlice) to compute the ecdf.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-15
# Pipeline version ABG 2020-09-07
#
# OUTPUTS: present_returnlevels***.csv
#
#~~~~~~~~~~~~~~~~~~~~~~~

##### SETUP #####------------------------------------------------------------
suppressMessages({
library(extRemes)
})

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

# thresh1 <- "POT2"
# ws1 <- "pc05"
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)
print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

# if(file.exists(paste0(data_wd, subfold, "eventdpe_OBS_",
#                       thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))){
#   stop("eventdpe_OBS_ already exists. ending 106.")
# }

subfold <- paste0("RCM", RCM, suffix, "/")


##### FUNCTIONS #####---------------------------------------

logit <- function(x){log(x/(1-x))}
invlogit <- function(y){1/(1 + exp(-1*y))}

gringorten <- function(v){
  ((length(v) + 1 - rank(v)) - 0.44)/(length(v) + 0.12)
}

weibull <- function(v){
  (length(v) + 1 - rank(v))/(length(v) + 1)
}

##### DATA #####--------------------------------------------------------------

print(ST <- Sys.time())

suffix_pres <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- ncoriginal <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix_pres, "_out.nc") 

# ncin <- nc_open(ncoriginal) # This file is ~2.5GB on the linux server.
# print(ncin)
# print(floor(Sys.time() - ST))

ncin_pres <- nc_open(ncpres)

# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                            RCM, suffix, ".rds"))
#dim(threshMat) #19914 x 5

#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))

# timewise maxima at each cell for each event (NE x NH)

obs_events  <- as.data.frame(readr::read_csv(paste0(data_wd,subfold, "eventflow_OBS_",
                          thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"),
                          col_types=cols(.default = col_double())))[,-(1:4)]

NE <- ncol(obs_events)

partable <- as.data.frame(readr::read_csv(paste0(data_wd,subfold_pres, 
                            "paramtable_",thresh1, "_RCM", RCM, suffix_pres, ".csv"),
                            col_types=cols(.default= col_double())))

colnames(partable)[1] <- "meanint"

##### PROB CALCULATION #####-------------------------------------------------

### prealloc -----------------------
eventDpeFrame <- matrix(NA, ncol=NE, nrow=NH)
eventApeFrame <- matrix(NA, ncol=NE, nrow=NH)

ST0 <- Sys.time()
ST <- Sys.time()
print("loop start")
for(h in 1:NH){
  
  if((h < 10) | (h %% 200 == 0)){ # time recording
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
  scaleH <- partable$scale[h]
  shapeH <- partable$shape[h]
  
  # vals <- ncvar_get(ncin, "dmflow",
  #                   start=c(rn$row[h], rn$col[h], 1),
  #                   count=c(1, 1, -1))
  # if(period == "future"){
  vals <- ncvar_get(ncin_pres, "dmflow",
                         start=c(rn$row[h], rn$col[h], 1),
                         count=c(1, 1, -1))
  #}
  #fv <- fevd(x=vals, threshold = threshMat[h,jT], type="GP")$results$par
  # DPoE per cell for event maxima based on whole present time series.
  #if(period == "present"){
  ecd <- ecdf(c(-1,vals,1e8))
  #}else{
  #  ecd <- ecdf(c(-1,vals_pres,1e8))
  #}
  valsdpe <- 1 - ecd(unlist(obs_events[h,]))
  
  # APoE per cell for event maxima based on either ecdf or gpa
  wh_ext <- (valsdpe < threshVal[jT])
  gpa_poe <- (1 - pevd(as.numeric(unlist(obs_events[h,])),
                       threshold=thr, scale=scaleH, shape=shapeH, type='GP'))
  if(any(gpa_poe < 1e-3)){
    print(paste("****", h))
    print(paste("peaksonly:", scaleH, shapeH, 1/min(gpa_poe)))
    print(paste("obs beyond 1000 yr", obs_events[h, gpa_poe < 1e-3]))
    gpa_poe[gpa_poe < 1e-3] <- 1e-3
    print(paste("maxGPA = ", thr - scaleH/shapeH))
  }
  #gpa_poe2 <- (1 - pevd(as.numeric(unlist(obs_events[h,])),
  #                      threshold=thr, scale=fv[1], shape=fv[2], type='GP'))
  
  #print(paste("num exceed: should be 60", sum(vals > thr)))
  #print(paste("allexceed:", fv[1], fv[2], 1/min(gpa_poe2)))
  #print(paste("peaksonly:", scaleH, shapeH, 1/min(gpa_poe)))
  
  #at_site_dpe[wh_ext] <- (2/360)*gpa_poe
  
  valsape <- 1 - exp(-valsdpe*360)
  valsape[wh_ext] <- 1 - exp(-gpa_poe[wh_ext]/meanInt)
                    # 1 - exp(-gpa_poe/meanInt), #gpa scaled to year
                    # 1 - exp(-valsdpe*360)) #dpoe scaled to year
  
  eventDpeFrame[h,]  <- valsdpe
  eventApeFrame[h,]  <- valsape
}

print("writing new files")
eventDpeFrame <- cbind(rn, round(as.data.frame(eventDpeFrame),8))
colnames(eventDpeFrame)[-(1:4)] <- paste0("E",1:(ncol(eventDpeFrame)-4))

eventApeFrame <- cbind(rn, round(as.data.frame(eventApeFrame),8))
colnames(eventApeFrame)[-(1:4)] <- paste0("E",1:(ncol(eventApeFrame)-4))

readr::write_csv(eventDpeFrame, path=paste0(data_wd,subfold,
                "eventdpe_OBS_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
readr::write_csv(eventApeFrame, path=paste0(data_wd,subfold,
                "eventape_OBS_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
#nc_close(ncin)
print(Sys.time())
######
#
# CONVERTION FROM PoE IN DAYS (p) TO PoE IN YEARS (b): b = 1 - (1-p)^(360)
#
# p = 1- (1-b)^(1/360)
#
#####