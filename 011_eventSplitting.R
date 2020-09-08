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

regions <- c("ANG", "NE", "NW")

if(length(args)==3){
  RCM <- sprintf("%02d", args[1])
  period <- args[2]
  if(period=="present"){
    suffix <- "_198012_201011"
  }else if (period=="future"){
    suffix <- "_205012_201011"
  }
  REG <- args[3]
  if(!any(REG == c()))
}

thresh1 <- "POT2" #!#!#!#!# Important constants to select.
ws1 <- "pc05"
print("Running for threshold", POT2, "at ", ws1, "minimum spread.")

library(readr)
library(dplyr)
library(ncdf4)
library(extRemes)
library(reshape2)
library(tidyverse)


### DATA ###--------------------------------------------------------------

HA <- readOGR(dsn="~/FEH_C/HydrometricAreas/temp", layer="hyd_areas",
              stringsAsFactors=FALSE)
HA@data$HA_NUM <- as.numeric(HA@data$HA_NUM)

NE <- length(eventDayList[[jV]][[jI]]) # POT2, 2% inun.

# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS( paste0(data_wd, subfold, "threshDayExcList_RCM",
                                    RCM, suffix,".rds"))


# matrix of threshold value (col) at a given cell (row)
threshMat <- read.csv(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"),
                      stringsAsFactors=FALSE)
#dim(threshMat) #19914 x 5


#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))
NE <- length(eventDayList[[jV]][[jI]]) # POT2, 2% inun.


# timewise maxima at each cell for each event ((NE + 2) x NH)
eventDF <- readr::read_csv(paste0(wd,"Data/eventdf_POT2_pc05.csv"))

# PoE under different computations with extra data. Tidy format.
present <- readr::read_csv(paste0(wd,"Data/present2_returnlevels.csv"))



### EVENT SPLITTING ###---------------------------------------------------

REG <- "NW"

r1 <- which(rn_regions$REGION == REG) # length = 1437
event_region <- eventDF[r1,]
present_region <- present %>% subset(loc %in% r1) 
thresh_region <- threshMat[r1,]

# Method 1: keep all events.

present_above <- present_region
write_csv(present_above, path=paste0(wd_id, "regionalEvents_exampleNW_1.csv"))


# Method 2: Only keep events where the threshold is passed somewhere.

present_above <- present_region %>% subset(val > thresh)
p_a_events <- unique(present_above$eventNo) #79 events
  print(length(p_a_events))
present_above2 <- present_region %>% subset(eventNo %in% p_a_events)

write_csv(present_above2, path=paste0(wd_id, "regionalEvents_exampleNW_2.csv"))


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

write_csv(present_above3, path=paste0(wd_id, "regionalEvents_exampleNW_3.csv"))





