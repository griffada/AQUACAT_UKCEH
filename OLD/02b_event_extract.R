######
# Adam Griffin, 2020-04-27
#
# Event extraction from daily flow. Practiced on one 30-year period of data.
# Rewritten for linux to reduce read-write times for netCDF.
# 
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-04-27
#
#####

##### SETUP #####------------------------------------------------------------

library(ncdf4)
library(raster)

threshVal <- c(5/365, 2/365, 1/365, 0.2/365, 0.1/365)
threshName <- c("POT5", "POT2", "POT1", "Q5", "Q10")
NT <- length(threshVal)


##### DATA #####------------------------------------------------------------
ncnameLin <-"/prj/ResilRiskInds/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
ncnameWin <- "V:/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"  
# This file is ~36GB on the linux server.

if(substr(osVersion,1,3)=="Win"){
  ncname <- ncnameWin
}else{
  ncname <- ncnameLin
}
# This file is ~36GB on the linux server.
print(ST <- Sys.time())
ncin <- nc_open(ncname)
print(ncin)
print(Sys.time() - ST)

ND <- 10957 # Number of days


hasData <- read.csv("./InterimData/hasData.csv", stringsAsFactors=FALSE)
NH <- nrow(hasData)


threshDayExcList <- loadRDS(file="./InterimData/threshDayExcList.rds") #5 lists of NH lists
thresMat <- loadRDS(file="./InterimData/threshMat.rds")


### Find days with most exceedences ###-----------------------------------

# Number of sites exceeded for each day at each threshold
inunMat <- matrix(0, nrow=length(NT), ncol=ND)

for(j in 1:NT){  # for each threshold value
  for(k in 1:NH){ # at each river network gridcell
    tde_jk <- threshDayExcList[[j]][[k]]
    
    for(n in 1:length(tde_jk)){
      tde <- tde_jk[n] # which days was it exceeded
      inunMat[j,tde] <- exceed[j,tde] + 1
    }
    
  }
}



# Percentage of inundated cells
wsBound <- c(0.05, 0.02, 0.01, 0.005, 0.001)
NW <- length(wsBound)

widespreadArray <- array(FALSE, dim=c(NT, ND, NW))
for(j in 1:NW){
  # Is there a widespread event on this day (at bound j)?
  widespreadArray[j,,] <- inunMat >= ND*wsBound[j] 
}


### Extract data at timepoints for HT/EC ###-----------------------------



