#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-15
#
# Estimating probability of exceedence along time series using ecdf and
# plotting positions.
# Uses inputs from 104_Event_Extract and 105_timeDF_compile.

# Note that for FUTURE estimation of return periods, make use of the present day
# data (tSlice) to compute the ecdf.
#
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-15
# Pipeline version ABG 2020-09-07
#
# Outputs: present_returnlevels***.csv
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
### Functions ###---------------------------------------

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

#library(data.table)
#edf <- fread(paste0(wd,"/Data/eventdf_POT2_pc2.csv"), colClasses=rep("numeric",287))

obs_events  <- readr::read_csv(paste0(data_wd,subfold, "eventdf_",thresh1,"_", ws1,
                                  "_RCM", RCM, suffix, ".csv"),
                           col_types=cols(.default = col_double()))

NE <- ncol(obs_events) - 4

thresMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM", RCM, suffix, ".rds"))

partable <- readr::read_csv(paste0(data_wd,subfold, 
                                   "paramtable_",thresh1, "_RCM", RCM, suffix, ".csv"))

colnames(partable)[1] <- "meanint"

##### PROB CALCULATION #####-------------------------------------------------

### prealloc -----------------------
rarityDF <- data.frame(eventNo = numeric(),
                       loc = numeric(),
                       Easting = numeric(),
                       Northing = numeric(), 
                       thresh = numeric(),
                       DayS = numeric(),
                       val = numeric(),
                       gpa_apoe = numeric(),
                       rp_years = numeric())

ST0 <- proc.time()
ST <- proc.time()
print("loop start")
for(h in 1:NH){
  thr <- thresMat[h,jT]
  meanInt <- partable$meanint[h]
  scaleH <- partable$scale[h]
  shapeH <- partable$shape[h]
  
  obs_events_h <- obs_events[h,-(1:4)]
  
  at_site_poeE <- 1 - pevd(as.numeric(obs_events_h), threshold=thr,
                           scale=scaleH, shape=shapeH, type='GP')
  
  rp_ver <- ifelse(at_site_poeE > 1-(1e-8), NA, meanInt/at_site_poeE)  
  # expected rate per year given POT and GPA PoE.
  
  at_site_apoe <- ifelse(is.na(rp_ver1E), NA, 1 - exp(-at_site_poeE/meanInt))  
  # Poisson assumption
  
  rarityTemp <- data.frame(eventNo = 1:NE,
                           loc = h,
                           Easting = rn[h, 1],
                           Northing = rn[h, 2], 
                           thresh = thr,
                           DayS = eventDayList[[jT]][[jW]],
                           val = obs_events_h,
                           gpa_apoe = at_site_apoe,
                           rp_years = rp_ver)
  
  rarityDF <- rbind(rarityDF, rarityTemp)
}

readr::write_csv(x=rarityDF,
                 path=paste0(data_wd, subfold, "returnlevels_",
                             thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
nc_close(ncin)
print(Sys.time())
######
#
# CONVERTION FROM PoE IN DAYS (p) TO PoE IN YEARS (b): p = 1 - (1-b)^(1/365.25)
#
# b = 1- (1-p)^(365.25)
#
#####

