#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2021-01-13
#
# Script for testing outputs from scripts 102-110 for plausibility.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2021-01-13
#
# Outputs:
#
#~~~~~~~~~~~~~~~~~~~~~~~

#### Packages -----------------------------------------------------------------

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE",
             "SEV", "SSC", "SW", "THA", "TRE", "WAL")

thresh1 <- "POT2" # Important constants to select.
ws1 <- "pc05"
print("Running 112_output_testing.R")
print(paste("Testing for threshold", thresh1, "at ", ws1, "minimum spread."))

jV <- which(threshName==thresh1)
jI <- which(wsName == ws1)

suppressPackageStartupMessages({
  library(readr)
  library(ncdf4)
  library(extRemes)
  library(reshape2)
  library(tidyverse)
})

#### THRESHMAT ----------------------------------------------------------------

threshMat <- readRDS(
                paste0(data_wd, subfold,"threshMat_RCM", RCM, suffix, ".rds"))

if(any(dim(threshMat) != c(NH, NT))){
  print("threshMat does not match expected dimensions: NH locations, NT thresholds")
  print(str(threshMat))
}
if(any(is.na(as.vector(threshMat)))){
  print(paste0("There are ", sum(is.na(as.vector(threshMat))),
               "missing values in threshMat."))
  print(which(is.na(threshMat), arr.ind=TRUE))
}
if(any((threshMat[,2:5] - threshMat[,1:4]) < 0)){
  print("thresholds not in ascending order. Check for problems in extraction.")
}


#### REGIONS ------------------------------------------------------------------

rn_regions <- readr::read_csv(paste0(data_wd, "hasData_Regions.csv"),
                              col_types=
                                cols(
                                  row = col_double(),
                                  col = col_double(),
                                  HA_NUM = col_double(),
                                  HA_NAME = col_character(),
                                  REGION = col_character()
                                ))

rn_regions$REGNUM <- as.numeric(factor(rn_regions$REGION))

genericMapPlotter(rn_regions[,1:4], rn_regions$REGNUM, gbgrid=gbgrid)


#### EVENTS -------------------------------------------------------------------

eventflow <- readr::read_csv(
                  paste0(data_wd,subfold, "eventflow_OBS_",thresh1,"_", ws1,
                         "_RCM", RCM, suffix, ".csv"),
                  col_types=cols(.default=col_double()))
