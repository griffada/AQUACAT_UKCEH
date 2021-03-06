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
#~~~~~~~~~~~~~~~~~~~~~~~

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

regions <- c("ANG", "NE", "NW", "SCO", "SE", "SEV", "SW", "THA", "TRE", "WAL")
REG <- NA
if(length(args)==3){
  RCM <- sprintf("%02d", args[1])
  period <- args[2]
  if(period=="present"){
    suffix <- "_198012_201011"
  }else if (period=="future"){
    suffix <- "_205012_201011"
  }
  REG <- args[3]
  if(!any(REG == regions)){
    stop(paste("incorrect call: Rscript 109_HeffTawn_Modelling.R gcm period region \n",
               "- Region must be one of: ANG, NE, NW, SCO, SE, SEV, SW, THA, TRE, WAL."))
  }
}

thresh1 <- "POT2" #!#!#!#!# Important constants to select.
ws1 <- "pc05"
print(paste0("Running for threshold ", thresh1, " at ", ws1, " minimum spread."))
jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)
library(readr)
library(dplyr)
library(ncdf4)
library(extRemes)
library(reshape2)
library(tidyverse)


### DATA ###--------------------------------------------------------------

HA <- readOGR(dsn=paste0(data_wd,"hydrometricAreas"), layer="hyd_areas",
              stringsAsFactors=FALSE)
HA@data$HA_NUM <- as.numeric(HA@data$HA_NUM)

# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshMat <- read_csv(paste0(wd_id, "threshMat2.csv"))
thresh0 <- unlist(threshMat['X2'], use.names=F)
#dim(threshMat)  =  19914 x 5

# eventLList length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(wd_id, "eventLists03.RDa")) 
NE <- length(eventDayList[[jT]][[jW]]) # POT2, 2% inun.

# timewise maxima at each cell for each event ((NE + 2) x NH)
eventDF <- readr::read_csv(paste0(data_wd,"TestData/eventdf_POT2_pc05.csv"))

# PoE under different computations with extra data. Tidy format.
present <- readr::read_csv(paste0(data_wd,"TestData/present_returnlevels_",
                                  thresh1,"_",ws1,".csv"))



### EVENT SPLITTING ###---------------------------------------------------

REG <- "ANG"

for(REG in regions){
  r1 <- which(rn_regions$REGION == REG) # length = 1437
  event_region <- eventDF[r1,]
  present_region <- present %>% subset(loc %in% r1) 
  thresh_region <- threshMat[r1,]

  # Method 3: Only keep events where inundation is >1% (reduced inundation for 
  # smaller area.)
  
  p_a_event_summ <- present_region %>%
    mutate(aboveThresh = (val > thresh)) %>% 
    dplyr::select(eventNo, aboveThresh) %>%
    group_by(eventNo) %>%
    summarise(inun_event = mean(aboveThresh)) %>%
    mutate(inund = inun_event >= 0.01)
  
  w1 <- p_a_event_summ$eventNo[which(p_a_event_summ$inund)] #44 events
  
  present_above3 <- present_region %>% subset(eventNo %in% w1)
  
  write_csv(present_above3, path=paste0(wd_id, "regionalEvents_example",REG,"_0.csv"))

}
event_region <- event_region[, (p_a_events+2)]