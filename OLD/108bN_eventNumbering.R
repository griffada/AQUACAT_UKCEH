#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-07-31
#
# Splitting events by region defined by Hydrological Area. 
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-07-31
# Pipeline version ABG 2020-09-07
# netcdf version ABG 2021-07-07
#
# OUTPUTS: regionalEvents_***.csv,
#          eventdf_region_***.csv
#
#~~~~~~~~~~~~~~~~~~~~~~~

if(interactive()){commandArgs <- function(...){c("04","future")}}

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}
regions <- c("NW")
print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

suppressPackageStartupMessages({
library(readr)
library(dplyr)
library(ncdf4)
library(extRemes)
library(reshape2)
library(tidyverse)
})

### DATA ###--------------------------------------------------------------

# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"))
thresh0 <- unlist(threshMat[,which(threshName==thresh1)], use.names=FALSE)
#dim(threshMat) #19914 x 5


obs_national <- nc_open(paste0(data_wd, subfold, "eventOBS_", thresh1,
                                     "_", ws1, "_RCM", RCM, suffix, ".nc"))
### OBS EVENT SPLITTING ###---------------------------------------------------

for(REG in c("NW")){
    
  if(!dir.exists(paste0(data_wd, subfold, "/", REG))){
    dir.create(paste0(data_wd, subfold, "/", REG))
  }
  
  print(paste("splitting", REG))
  
  r1 <- which(rn_regions$REGION == REG) # length = 1437
  thresh_region <- threshMat[r1, jT]
  
  savepath <- paste0(data_wd,subfold, "/", REG, "/",
                  "eventOBS_region_", REG, "_RCM", RCM, suffix, ".nc")
  
  rn0 <- rn
  rn <- rn[r1,]
  
  cdfPrimer(RCM=RCM, period=period, method="OBS", NE=250, NH=length(r1),
            thresh1=thresh1, ws1=ws1, rn=rn, savepath=savepath)
  obs_region <- nc_open(savepath, write=T)
  
  rn <- rn0
  
  NEreg <- 0
  thresh_region <- threshMat[r1, jT]
  NE <- obs_national$dim$event$len
  WE <- c()
  obs_event_number <- ncvar_get(obs_national, "eventNo")
  ST <- Sys.time()
  ST0 <- Sys.time()
  for(i in 1:NE){
    if((i < 10) | (i %% 200 == 0)){ # time recording
      print(i)
      I <- difftime(Sys.time(), ST, units="secs")
      I0 <- difftime(Sys.time(), ST0, units="secs")
      print(paste("Percent remaining", 100*round((NE-i)/NE ,2)))
      print(paste("Time remaining", round((NE-i)/i * I0,2)))
      print(paste("Since last readout:", round(I,2)))
      ST <- Sys.time()
      print(ST) 
    }
    
    obs_events_region <- ncvar_get(obs_national, "flow", start=c(1,i),
                                   count=c(-1,1))[r1]
    if(any(obs_events_region > thresh_region)){
      WE <- c(WE, obs_event_number[i])
    }
  }
  ncvar_put(obs_region, "eventNo", WE, start=1, count=length(WE))
  print(paste(NE, "OBS events saved for", REG, "region."))
}