#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-15
#
# Estimating probability of exceedence along time series using ecdf and
# plotting positions.
# Uses inputs from 102_Threshold_Extract, 104_Event_Extract and 105_Event_Summary

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
library(extRemes)

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

thresh1 <- "POT2"
ws1 <- "pc05"
print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

if(file.exists(paste0(data_wd, subfold, "returnlevels_",
                      thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))){
  stop("returnlevels already exists. ending 106.")
}

subfold <- paste0("RCM", RCM, suffix, "/")

jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)
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

ncin <- nc_open(paste0(data_wd, subfold_pres, "dmflow_copy_RCM",
                       RCM, suffix_pres, ".nc")) # This file is ~2.5GB on the linux server.
print(ncin)
print(floor(Sys.time() - ST))

# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                            RCM, suffix, ".rds"))
#dim(threshMat) #19914 x 5

#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))

# timewise maxima at each cell for each event (NE x NH)

obs_events  <- readr::read_csv(paste0(data_wd,subfold, "eventflow_OBS_",thresh1,"_", ws1,
                                  "_RCM", RCM, suffix, ".csv"),
                           col_types=cols(.default = col_double()))

NE <- ncol(obs_events) - 4

partable <- readr::read_csv(paste0(data_wd,subfold, 
                            "paramtable_",thresh1, "_RCM", RCM, suffix, ".csv"))

colnames(partable)[1] <- "meanint"

##### PROB CALCULATION #####-------------------------------------------------

### prealloc -----------------------
# rarityDF <- data.frame(eventNo = numeric(),
#                        loc = numeric(),
#                        Easting = numeric(),
#                        Northing = numeric(), 
#                        thresh = numeric(),
#                        DayS = numeric(),
#                        val = numeric(),
#                        gpa_apoe = numeric(),
#                        rp_years = numeric())
eventDpeFrame <- matrix(NA, ncol=ncol(obs_events), nrow=NH)
eventApeFrame <- matrix(NA, ncol=ncol(obs_events), nrow=NH)

ST0 <- proc.time()
ST <- proc.time()
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
  
  thr <- thresMat[h,jT]
  meanInt <- partable$meanint[h]
  scaleH <- partable$scale[h]
  shapeH <- partable$shape[h]
  
  vals <- ncvar_get(ncin, "dmflow",
                    start=c(rn$row[h], rn$col[h], 1),
                    count=c(1, 1, -1))
  
  # DPoE per cell for event maxima based on whole time series.
  ecd <- ecdf(vals)
  valsdpe <- 1 - ecd(obs_events[h,])
  
  # APoE per cell for event maxima based on either ecdf or gpa
  wh_ext <- (valsdpe < (2/360))
  gpa_poe <- (1 - pevd(as.numeric(obs_events[h,]),
                       threshold=thr, scale=scaleH, shape=shapeH, type='GP'))
  
  at_site_dpe[wh_ext] <- (2/360)*gpa_poe
  
  valsape <- ifelse(wh_ext,
                    1 - exp(-gpa_poe/meanInt), #gpa scaled to year
                    1 - exp(-valsdpe/360)) #dpoe scaled to year
  
  eventDpeFrame[h,]  <- valsdpe
  eventApeFrame[h,]  <- valsape
}

# eventDpeFrame  <- cbind(rn, eventDpeFrame)
# eventApeFrame  <- cbind(rn, eventApeFrame)
# 
# readr::write_csv(x=rarityDF,
#                  path=paste0(data_wd, subfold, "returnlevels_",
#                              thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

readr::write_csv(as.data.frame(eventDpeFrame), path=paste0(data_wd,subfold,
                "eventdpe_OBS_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
readr::write_csv(as.data.frame(eventApeFrame), path=paste0(data_wd,subfold,
                "eventape_OBS_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
nc_close(ncin)
print(Sys.time())
######
#
# CONVERTION FROM PoE IN DAYS (p) TO PoE IN YEARS (b): p = 1 - (1-b)^(360)
#
# b = 1- (1-p)^(1/360)
#
#####