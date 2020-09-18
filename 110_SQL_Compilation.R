#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-09-10
#
# SQL database construction for HT/EC modelled events.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-09-10
#
# Outputs:
#    ***.sqlite
#
#~~~~~~~~~~~~~~~~~~~~~~~



##### SETUP -------------------------------------------------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

regions <- c("ANG", "NE", "NW", "SCO", "SE", "SEV", "SW", "THA", "TRE", "WAL")

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


library(dbplyr)
library(RSQLite)
library(reshape2)

##### DATA-------------------------------------------------------------

PRES <- readr::read_csv(paste0(data_wd,"present_returnlevels_POT2_pc05.csv"))

NEWEC <- readr::read_csv(paste0(data_wd,"NewEventPresentEC_POT2_pc05.csv"))

load(paste0(wd_id, "eventLists03.RDa"))


