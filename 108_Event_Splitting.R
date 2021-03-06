#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-07-31
#
# Splitting events by region defined by Hydrological Area. 
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-07-31
# Pipeline version ABG 2020-09-07
#
# OUTPUTS: regionalEvents_***.csv,
#          eventdf_region_***.csv
#
#~~~~~~~~~~~~~~~~~~~~~~~

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE", "SEV", "SSC", "SW", "THA", "TRE", "WAL")

# if(length(args)==3){
#   RCM <- sprintf("%02d", as.numeric(args[1]))
#   period <- args[2]
#   if(period=="present"){
#     suffix <- "_198012_201011"
#   }else if (period=="future"){
#     suffix <- "_205012_201011"
#   }
#   REG <- args[3]
#   if(!any(REG == regions)){
#     stop("incorrect call: Rscript 108_eventSplitting.R gcm period region - Region must be one of: ANG, ESC, NE, NW, NSC, SE, SEV, SSC, SW, THA, TRE, WAL.")
#   }
# }

thresh1 <- "POT2" # Important constants to select.
ws1 <- "pc05"
print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

jV <- which(threshName==thresh1)
jI <- which(wsName == ws1)

library(readr)
library(dplyr)
library(ncdf4)
library(extRemes)
library(reshape2)
library(tidyverse)


### DATA ###--------------------------------------------------------------

#Hydrometric areas map
HA <- readOGR(dsn=paste0(data_wd,"hydrometricAreas"), layer="hyd_areas",
              stringsAsFactors=FALSE)
HA@data$HA_NUM <- as.numeric(HA@data$HA_NUM)



# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS( paste0(data_wd, subfold, "threshDayExcList_RCM",
                                    RCM, suffix,".rds"))


# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"))
thresh0 <- unlist(threshMat[,which(threshName==thresh1)], use.names=FALSE)
#dim(threshMat) #19914 x 5


#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))
NE <- length(eventDayList[[jV]][[jI]]) # POT2, 2% inun.

# timewise maxima at each cell for each event ((NE + 2) x NH)
obs_events <- readr::read_csv(paste0(data_wd, subfold, "eventflow_OBS_", thresh1,
                                     "_", ws1, "_RCM", RCM, suffix, ".csv"))[,-(1:4)]

obs_dpe <- readr::read_csv(paste0(data_wd, subfold, "eventdpe_OBS_", thresh1,
                                       "_", ws1, "_RCM", RCM, suffix, ".csv"))[,-(1:4)]
  
obs_ape <- readr::read_csv(paste0(data_wd, subfold, "eventape_OBS_", thresh1,
                                       "_", ws1, "_RCM", RCM, suffix, ".csv"))[,-(1:4)]

# # PoE under different computations with extra data. Tidy format.
# present <- readr::read_csv(paste0(data_wd, subfold, "returnlevels_",
#                               thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))



### EVENT SPLITTING ###---------------------------------------------------

for(REG in regions){
  
  print(paste("splitting", REG))
  
  r1 <- which(rn_regions$REGION == REG) # length = 1437
  obs_events_region <- obs_events[r1, ]
  obs_dpe_region <- obs_dpe[r1, ]
  obs_ape_region <- obs_ape[r1, ]
  #present_region <- present %>% subset(loc %in% r1) 
  thresh_region <- threshMat[r1, ]
  
  KEEP <- apply(obs_events_region, 2, function(x){any(x > thresh_region)})
  
  obs_dpe_region <- obs_dpe_region[,KEEP]
  obs_dpe_region <- cbind(rn, obs_dpe_region)
  
  obs_events_region <- obs_events_region[,KEEP]
  obs_events_region <- cbind(rn, obs_events_region)
  
  
  obs_ape_region <- obs_ape_region[,KEEP]
  obs_ape_region <- cbind(rn, obs_ape_region)
  
  # p_a_events <- unique(present_above$eventNo) #79 events
  #   print(length(p_a_events))
  # present_above2 <- present_region %>% subset(eventNo %in% p_a_events)
  # 
  # event_region <- event_region[, (p_a_events+2)]
  
  ### SAVE OUTPUTS ###------------------------------------------------
  
  print("saving outputs")
  
  if(!dir.exists(paste0(data_wd, subfold, "/", REG))){
  dir.create(paste0(data_wd, subfold, "/", REG))
  }
  
  write_csv(obs_events_region, path=paste0(data_wd, subfold, "/", REG, "/eventflow_OBS_region_",
                                        REG,"_RCM", 
                                        RCM, suffix,".csv"))
  
  write_csv(obs_dpe_region, path=paste0(data_wd, subfold,  "/", REG, "/eventdpe_OBS_region_",
                                        REG,"_RCM", 
                                        RCM, suffix,".csv"))
  
  write_csv(obs_ape_region, path=paste0(data_wd, subfold,  "/", REG, "/eventape_OBS_region_",
                                      REG,"_RCM", 
                                      RCM, suffix,".csv"))
}
print("108 done.")