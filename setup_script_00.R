#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-02-21
#
# A proper setup file which can be run whenever needed.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-02-21, edited 2020-05-04
#
#~~~~~~~~~~~~~~~~~~~~~~~~

##### PARSE COMMAND LINE ARGUMENTS FIRST ##### -------------------------------
{
args <- commandArgs(trailingOnly=TRUE)
if(length(args) > 3){
  stop("incorrect call. Rscript 10X_*****.R [RCM] [period]. ")
}
if(length(args)==0){
  RCM <- "01"
  period <- "present"
  suffix <- "_198012_201011"
}else if(length(args)==1){
  RCM <- sprintf("%02d",as.numeric(args[1]))
  period <- "present"
  suffix <- "_198012_201011"
}else if(length(args)>=2){
  RCM <- sprintf("%02d", as.numeric(args[1]))
  period <- args[2]
  if(period=="present"){
    suffix <- "_198012_201011"
  }else if (period=="future"){
    suffix <- "_205012_208011"
  }else{  
    stop("incorrect call: Rscript 10X_*****.R RCM [period]. Period should be 'present' or 'future'.")
  }
}
if(as.numeric(RCM) < 0 | as.numeric(RCM) > 16){
  stop("correct call: Rscript 10X_*****.R RCM [period]. RCM should be between 1 and 16.")
}
}
#### CORE PACKAGES ---------------------------------------------------
print(Sys.time())
library(ncdf4)
library(raster)
library(fields)
library(rgdal)
library(readr)

### SET WORKING DIRECTORY ###
if (substr(osVersion,1,3) == "Win") {
  
  wd <- "S:/CodeABG/"
  wd_id <- "S:/CodeABG/InterimData/"
  data_wd <- "S:/Data/"
  g2g_wd <- "S:/run_hmfg2g/outputs/"
  
}else{
  
  wd <- "/prj/aquacat/CodeABG/"
  wd_id <- "/prj/aquacat/CodeABG/InterimData/"
  data_wd <- "/prj/aquacat/Data/"
  g2g_wd <- "/prj/aquacat/run_hmfg2g/outputs/"
  
}

if(!dir.exists(paste0(data_wd, "RCM", RCM, suffix))){
  dir.create(paste0(data_wd, "RCM", RCM, suffix))
}
subfold <- paste0("RCM", RCM, suffix, "/")


ncoriginal <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix, "_out.nc") 

ncname <- paste0(data_wd, subfold, "dmflow_copy_RCM", RCM, suffix, ".nc")

# threshold for inundation
threshVal <- c(5/365, 2/365, 1/365, 0.2/365, 0.1/365)
threshName <- c("POT5", "POT2", "POT1", "Q5", "Q10")
NT <- length(threshVal)

# cut-off for widespread event
wsCutoff <- c(0.05, 0.02, 0.01, 0.005, 0.001)
wsName <- c("pc5", "pc2", "pc1", "pc05", "pc01")
NW <- length(wsCutoff)

# location of river network
rn <- read_csv(paste0(wd_id,"hasData3.csv"))
#rn <- read_csv(paste0(data_wd, "hasData_primary.csv"))
NH <- nrow(rn)

rn_regions <- try({read_csv(paste0(wd_id, "hasData_Regions.csv"))})

ND <- 10800 # Number of days
