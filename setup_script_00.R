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
# if(interactive()){
#   commandArgs <- function(...){c("10","present","NW")}
# }
##### PARSE COMMAND LINE ARGUMENTS FIRST ##### -------------------------------

args <- commandArgs(trailingOnly=TRUE)
FF_flag <- FALSE
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

regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE", "SEV",
             "SSC", "SW", "THA", "TRE", "WAL")


if(length(args)==3){
  
  if(args[3] == "FF" & period == "future"){
    FF_flag <- TRUE
  }else{
    REG <- args[3]
    if(!(REG %in% regions)){
      stop(paste("incorrect call: Rscript 109_HeffTawn_Modelling.R gcm period [region] \n",
      "- Region must be one of: ANG, ESC, NE, NSC, NW, SE, SEV, SSC, SW, THA, TRE, WAL."))
    }
  }
}
if(as.numeric(RCM) < 0 | as.numeric(RCM) > 16){
  stop("correct call: Rscript 10X_*****.R RCM [period]. RCM should be between 1 and 16.")
}
}
#### CORE PACKAGES ---------------------------------------------------
print(Sys.time())
options("rgdal_show_exportToProj4_warnings"="none")
options(warn=1)
suppressMessages({
library(ncdf4)
library(raster)
library(fields)
library(rgdal)
library(readr)
library(yaml)
})

readdf <- function(...){as.data.frame(data.table::fread(...))}

### SET WORKING DIRECTORY ###
if (substr(osVersion,1,3) == "Win") {
  
  wd <- "S:/CodeABG/"
  wd_id <- "S:/CodeABG/InterimData/"
  data_wd <- "S:/Data/"
  g2g_wd <- "S:/run_hmfg2g/outputs/"
  
}else if (substr(osVersion,1,3) == "Fed"){
  
  wd <- "/prj/aquacat/CodeABG/"
  wd_id <- "/prj/aquacat/CodeABG/InterimData/"
  data_wd <- "/prj/aquacat/Data/"
  g2g_wd <- "/prj/aquacat/run_hmfg2g/outputs/"
  
}else{
  wd <- "~/AQUACAT/CodeABG/"
  wd_id <- "~/AQUACAT/CodeABG/InterimData/"
  data_wd <- "~/AQUACAT/Data/"
  g2g_wd <- "~/AQUACAT/run_hmfg2g/outputs/"
}


if(FF_flag){
  subfold <- paste0("RCM", RCM, suffix, "_FF/")
}else{
  subfold <- paste0("RCM", RCM, suffix, "/")
}

if(!dir.exists(paste0(data_wd, subfold))){
  dir.create(paste0(data_wd, subfold))
}

ncoriginal <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix, "_out.nc") 

ncname <- paste0(data_wd, subfold, "dmflow_copy_RCM", RCM, suffix, ".nc")

# threshold for inundation
threshVal <- c(5/360, 2/360, 1/360, 0.2/360, 0.1/360)
threshName <- c("POT5", "POT2", "POT1", "Q5", "Q10")
NT <- length(threshVal)

# cut-off for widespread event
wsCutoff <- c(0.05, 0.02, 0.01, 0.005, 0.001)
wsName <- c("pc5", "pc2", "pc1", "pc05", "pc01")
NW <- length(wsCutoff)

# location of river network
rn <- read_csv(paste0(data_wd,"hasData_primary.csv"), 
                         col_types=cols(
                           row = col_double(),
                           col = col_double(),
                           east = col_double(),
                           nor = col_double()
                         ))
# print("*********** THIS IS USING A SHORT VERSION OF RN, CHANGE BACK **********************")
# print("*********** THIS IS USING A SHORT VERSION OF RN, CHANGE BACK **********************")
# rn <- rn[17000:18999,]
#rn <- read_csv(paste0(data_wd, "hasData_primary.csv"))
NH <- nrow(rn)

rn_regions <- try({
  read_csv(paste0(data_wd,"hasData_Regions.csv"),
    col_types = cols(
      row = col_double(),
      col = col_double(),
      HA_NUM = col_double(),
      HA_NAME = col_character(),
      REGION = col_character()
    )
  )
})

ND <- 10800 # Number of days



# Msims <- 25000
# nSampleHT <- 400
# RECOMPUTE_FLAG <- FALSE

source(paste0(wd,"extra_functions.R"))

settingspath <- paste0(data_wd,subfold,"settings.yaml")
if(file.exists(settingspath)){
  settings <- read_yaml(settingspath)

thresh1 <- settings$thresh
ws1 <- settings$wsname
Msims <- settings$Msims
nSampleHT <- settings$nSampleHT
# Only using POT2 and 2% inundation minimums.
jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)
print(paste("Setup for threshold", thresh1, "at", ws1, "minimum spread."))
}

