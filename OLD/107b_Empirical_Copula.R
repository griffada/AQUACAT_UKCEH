#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-15
#
# The Empirical Copula functions for simulating new events from given data.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-15
# Pipeline version ABG 2020-09-07
#
# OUTPUTS: NewEventPresentEC_***.csv: data table of events, one event per row.
#
#~~~~~~~~~~~~~~~~~~~~~~~

if(interactive()){commandArgs <- function(...){c("10","present")}}
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  #source("C:/Users/adagri/Documents/AQUACAT_C/CodeABG/setup_script_00.R")
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

##### SETUP #####---------------------------------------------------

suppressPackageStartupMessages({
library(extRemes)
library(dplyr)
library(fitdistrplus)
library(lmomco)
})

### DATA --------------------------

# thresh1 <- "POT2" #!#!#!#!# Important constants to select.
# ws1 <- "pc05"
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)

#print(paste("Running for threshold", thresh1, "at", ws1, "minimum spread."))
M <- Msims
print(paste("> > > > > >M = ", M))
suffix_pres <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix_pres, "_out.nc") 



# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS(paste0(data_wd, subfold, "threshDayExcList_RCM",
                                    RCM, suffix,".rds"))

# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"))
#dim(threshMat) #19914 x 5

#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))

NE <- length(eventDayList[[jT]][[jW]]) # POT2, 0.5% inun.

NH <- nrow(rn)

partable <- as.data.frame(readr::read_csv(paste0(data_wd,subfold, 
                           "paramtableG_",thresh1, "_RCM", RCM, suffix, ".csv"),
                            col_types=cols(.default= col_double())))

partable_pres <- as.data.frame(readr::read_csv(paste0(data_wd,subfold_pres, 
                        "paramtableG_",thresh1, "_RCM", RCM, suffix_pres, ".csv"),
                            col_types=cols(.default= col_double())))

# colnames(partable)[1] <- "meanInt"
# colnames(partable_pres)[1] <- "meanInt"

newEventDpe <- as.data.frame(readr::read_csv(paste0(data_wd, subfold, "eventdpe_EC_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"),
                        col_types=cols(.default= col_double()))[,-(1:4)])

returnLevelsEC <- function(eventSimTable, ncin, paramtable,
                           NH=19914, pthr=(2/360)){
  # Gets empirical quantiles for simulated PoE, and GPA fitted values of flow
  # for simulated PoE below threshold (e.g 2/360 for POT2).
  #
  # eventSimTable     matrix of DPoE, one event per column, one location per row.
  # ncin              netCDF object showing daily flow
  # paramtable        GPA parameters table for POTs
  # NH                number of locations
  # pthr              probability for threshold exceedence, default = POT2
  
  ST0 <- Sys.time()
  ST <- Sys.time()
  eventFlow <- eventSimTable
  for(i in 1:NH){
    
    if((i < 10) | (i %% 200 == 0)){ # time recording
      print(i)
      I <- difftime(Sys.time(), ST, units="secs")
      I0 <- difftime(Sys.time(), ST0, units="secs")
      print(paste("Percent remaining", 100*round((NH-i)/NH ,2)))
      print(paste("Time remaining", round((NH-i)/i * I0,2)))
      print(paste("Since last readout:", round(I,2)))
      ST <- Sys.time()
      print(ST) 
    }
    
    nor <- rn$row[i]
    eas <- rn$col[i]
    tSlice <- ncvar_get(ncin, varid="dmflow",
                        start=c(nor, eas,  1),
                        count=c(1, 1, -1))
    rare <- (eventSimTable[i,] < pthr)
    #print(sum(rare))
    quant <- quantile(tSlice, probs=1-eventSimTable[i,])
    # fit exceedences using obs if under threshold
    #print(i)
    # use GPA if over threshold (note divide by pthr to get P[X | X>thr])
    if(sum(rare)>0){

    o <- try(quaglo(1 - pmax(1e-7,eventSimTable[i,rare]),
                        vec2par(c(paramtable$loc[i],
                                  paramtable$sca[i],
                                  paramtable$shape[i]), type="glo")))
    
    if(inherits(o, "try-error")){
      print(i)
      print(range(eventSimTable[i,rare]))
      if(interactive()) browser()
      warning("Error in returnLevelsEC")
    }else{
      quant[rare] <- o
    }
    }else{
      print("length rare = 0")
    }
    eventFlow[i,] <- quant
  }
  eventFlow
}

### CONVERSION TO FLOW AND APoE ###------------------------------------------
print("Conversion to flow")
# Conversion to flow
ncin <- nc_open(ncpres)
newEventFlow <- returnLevelsEC(as.matrix(newEventDpe),
                               ncin=ncin,
                               NH=NH,
                               paramtable=partable_pres,
                               pthr=threshVal[jT])
summary(warnings())
colnames(newEventFlow) <- colnames(newEventDpe)

readr::write_csv( cbind(rn, round(as.data.frame(newEventFlow),8)),
                 paste0(data_wd, subfold, "eventflow_EC_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))


print(Sys.time())

#### SAVE OUTPUTS ####--------------------------------------------------------
print("Saving outputs")

# 
# W <- paste0("E",1:(ncol(newEventFlow)-4))
# W1 <- sapply(1:length(W), function(i){paste0(W[i],".",sum(W[1:i]==W[i]))})
# colnames(newEventApe) <- W1
# newEventApe <- cbind(rn, round(as.data.frame(newEventApe),8))
# readr::write_csv(newEventApe,
#                  paste0(data_wd, subfold, "eventape_EC_",
#                         thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

print("107 done.")