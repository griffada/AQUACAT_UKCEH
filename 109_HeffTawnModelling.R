#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-17
#
# Using Heffernan and Tawn model for spatial coherence to generate new events.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-17
# Pipeline version 2020-09-07
#
#~~~~~~~~~~~~~~~~~~~~~~~

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

##### SETUP #####-----------------------------------------------------------
library(texmex)
library(dplyr)
library(extRemes)
library(reshape2)
library(tidyverse)
  
ws1 <- "pc05"
thresh1 <- "POT2"

print("Running for threshold", POT2, "at ", ws1, "minimum spread.")

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



# timewise maxima at each cell for each event (NE x NH)

#library(data.table)
#edf <- fread(paste0(wd,"/Data/eventdf_POT2_pc2.csv"), colClasses=rep("numeric",287))

eventDF <- readr::read_csv(paste0(data_wd,subfold, "eventdf_",thresh1,"_", ws1,
                                  "_RCM", RCM, suffix, ".csv"))
  
  
# PoE under different computations with extra data. Tidy format.
present <- readr::read_csv(paste0(data_wd, subfold, "present_returnlevels_",
                                  thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
 
print("SETUP DONE.")
##### HEFFTAWN CODING #####------------------------------------------------




