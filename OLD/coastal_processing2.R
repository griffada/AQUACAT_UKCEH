#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2021-04-29
#
# Event threshold extraction from daily flow for various thresholds.
# A GPA distribution is fitted to the peaks over these thresholds.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-02-21, edited 2020-05-04
# Pipeline version 2020-09-07
# Version b uses GLO 2021-03-23
#
# OUTPUTS: threshDayExcList.rda,
#          thresMat.rda,
#          paramtable.csv
#
#~~~~~~~~~~~~~~~~~~~~~~~
print("running 102_sea")
if(interactive()){commandArgs <- function(...){c("01","present")}}
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

library(ilaprosUtils)
library(extRemes)
library(lmomco)
library(data.table)
library(tidyverse)
L <- list.files(paste0(data_wd,"SeaLevels"), full.names=T)

LN <- length(L)

ND_true <- 10957

obs_table <- matrix(NA, nrow=ND, ncol=LN)
days_removed <- floor(seq(0,10958, length.out=159)[2:158])
days_kept <- (1:ND_true)[-days_removed]
daynames_removed <- as.Date("1980-11-30") + days_removed
daynames_kept <- as.Date("1980-11-30") + days_kept

for(l in 1:length(L)){
  print(l)
  df <- as.data.frame(fread(L[l]))
  colnames(df) <- c("No","day","flow")
  df$day <- as.Date(df$day)
  df <- df %>% dplyr::filter(day > as.Date("1980-11-30") & 
                             day < as.Date("2010-12-01") &
                           (day %in% daynames_kept))
  obs_table[which((daynames_kept %in% df$day)),l] <- df$flow
}

obs_table[obs_table < -98] <- NA

threshs <- apply(obs_table, 2, function(x){quantile(x, probs=c(1 - (2/360)), na.rm=T)})








