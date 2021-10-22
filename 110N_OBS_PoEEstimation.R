#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-17
#
# Estimating probability of exceedence for Heffernan-Tawn events.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-17
# Pipeline version 2020-09-07
# netcdf version 2021-07-07
#
# OUTPUTS: HTraritydf_***.csv
#
#~~~~~~~~~~~~~~~~~~~~~~~
if(interactive()){
  commandArgs <- function(...){c("01","present","NW")}
}
#### SETUP ####----------------------------------------------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

# if(settings$HTdpe & settings$HTape){
  # print("HT probs exist. Stopping 110N.")
# }else
{

library(ilaprosUtils)
library(lmomco)
library(extRemes)
library(pastecs)
library(dplyr)

# ws1 <- "pc05"
# thresh1 <- "POT2"
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)
# print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE", "SEV", "SSC",
             "SW", "THA", "TRE", "WAL")

if(length(args)==3){
  RCM <- sprintf("%02d", as.numeric(args[1]))
  period <- args[2]
  if(period=="present"){
    suffix <- "_198012_201011"
  }else if (period=="future"){
    suffix <- "_205012_208011"
  }
  REG <- args[3]
  if(!any(REG == regions)){
    stop(paste("incorrect call: Rscript 110_HT_PoEEstimation.R gcm period region \n",
               "- Region must be one of: ANG, ESC, NE, NSC, NW, SE, SEV, SSC, SW, THA, TRE, WAL."))
  }
}

# suffix_pres <- "_198012_201011"
# subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix, "_out.nc") 

##### DATA #####-------------------------------------------------------

#5 lists of NH lists, one for each threshold

threshMat <- readdf(paste0(data_wd, subfold, "threshMat_RCM",
                           RCM, suffix, ".csv"))

ST <- Sys.time()
ncin <- nc_open(ncpres)
print(Sys.time() - ST)

partable <- readdf(paste0(data_wd,subfold, 
                        "paramtableG_",thresh1, "_RCM", RCM, suffix, ".csv"))

ht_events <- nc_open(paste0(data_wd,subfold, REG, "/eventOBS_region_",
                  REG,"_RCM", RCM, suffix, ".nc"), write=T)

NE <- sum(ncvar_get(ht_events, "eventNo") > 0)
rn_regions$locnum <- 1:nrow(rn_regions)
rn_reg <- data.frame(rn_regions) %>% dplyr::filter(REGION == REG)
r1 <- which(rn_regions$REGION == REG)
NH1 <- length(r1)
partable <- partable[r1,]

### Prealloc ###--------------------------------------------------------------

##### DPoE Calculation using GPA #####-----------------------------------------
ST <- Sys.time()
ST0 <- Sys.time()
for(h in 1:NE){
  
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
  
  
  H <- rn_reg$locnum[h]
  thr  <- threshMat[H,jT]
  locH <- partable$loc[h]
  scaleH <- partable$sca[h]
  shapeH <- partable$shape[h]
  thresholdH <- partable$threshold[h]
  
  vals <- ncvar_get(ncin, "dmflow",
                    start=c(rn$row[H], rn$col[H], 1),
                    count=c(1, 1, 10800))
  
  eventflow_HT <- ncvar_get(ht_events, "flow",
                    start=c(1,h), count=c(NE,1))

  da <- dpeApeComputer(h,
                 vals,
                 obs=eventflow_HT,
                 ncin=ht_events,
                 pars=partable[h,],
                 thresh_val=thresholdH)
}

##### OUTPUTS #####-----------------------------------------------------------

nc_close(ht_events)
nc_close(ncin)

settings$OBSdpe <- TRUE
settings$OBSape <- TRUE
settings$HTPoE_Estimation <- "110cN"
write_yaml(settings, settingspath)

}
print("110N complete.")
print(Sys.time())
