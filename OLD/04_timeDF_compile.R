######
# Adam Griffin, 2020-04-27
#
# Summarising size and time of extreme events. Practiced on one 30-year period of data.
# 
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-04-27
#
#####

##### SETUP #####------------------------------------------------------------

library(ncdf4)
library(raster)

if(substr(osVersion,1,3)=="Win"){
  ncname <- "V:/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
  wd <- "S:"
}else{
  ncname <- "/prj/ResilRiskInds/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
  wd <- "/prj/aquacat"
}

# threshold for inundation
threshVal <- c(5/365, 2/365, 1/365, 0.2/365, 0.1/365)
threshName <- c("POT5", "POT2", "POT1", "Q5", "Q10")
NT <- length(threshVal)

# cut-off for widespread event
wsCutoff <- c(0.05, 0.02, 0.01, 0.005, 0.001)
NW <- length(wsCutoff)


##### DATA #####------------------------------------------------------------

print(ST <- Sys.time())
	ncin <- nc_open(ncname) # This file is ~36GB on the linux server.
	print(ncin)
print(Sys.time() - ST)

ND <- 10957 # Number of days

#cells on the ruver network
rn <- read.csv(paste0(wd,"/CodeABG/InterimData/hasData.csv",
				stringsAsFactors=FALSE)
NH <- nrow(rn)


threshDayExcList <- readRDS("/CodeABG/InterimData/threshDayExcList.rds")
  # 5 lists of NH lists
thresMat <- read.csv("/CodeABG/InterimData/threshMat.csv", stringsAsFactors=FALSE)



### Find days with most exceedences ###-----------------------------------

# Number of sites exceeded for each day at each threshold
inunMat <- array(0, dim=c(NT,ND)) # thres x days
dimnames(inunMat) <- list(threshold=threshName,
                          day=1:ND)

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
wsCutoff <- c(0.05, 0.02, 0.01, 0.005, 0.001)
wsCutoffNames <- c("pc5", "pc2", "pc1", "pc0.5", "pc0.1")
NW <- length(wsCutoff)

widespreadArray <- array(FALSE, dim=c(NT, NW, ND)) # thres x cutoff x days
dimnames(widespreadArray) <- list(threshold=threshName,
                                  cutoff=wsCutoffNames,
                                  days=1:ND)
for(j in 1:NW){
  # Is there a widespread event on this day (at bound j)?
  widespreadArray[,j,] <- inunMat >= ND*wsBound[j] 
}

# How many events for each thres/cutoff combo
eventCount <- apply(widespreadArray, c(1,2), function(x){sum(1*x)})
saveRDS(eventCount, paste0(wd,"/CodeABG/InterimData/eventCount.rds"))

### Extract data at timepoints for HT/EC ###-----------------------------

eventArrayList <- vector("list", 5)
names(eventArrayList) <- wsCutoffNames

for(i in 1:NW){  # for each ws cutoff
  
  eventArrayList[[i]] <- vector("list",5)
  names(eventArrayList[[i]]) <- threshName
  
  for(j in 1:NT){  #for each inun threshold
    
    eventArray <- matrix(NA, nrow=eventCount[i,j], ncol=NH)
    
    eDays <- which(widespreadArray[i,j,] == TRUE)
    rownames(eventArray) <- eDays
    
    for(k in 1:eventCount[i,j]){ # for each event
      ncwide <- ncvar_get(ncin,
                          "dmflow",
                          start=c(1,1,eDays[k]),
                          count=c(-1, -1, 1))
      for(l in 1:NH){
        # add the event to the event array
        eventArray[k, l] <- rn[hd[l, 1], rn[l, 2]]
      }
    }
    
    eventArrayList[[i]][[j]] <- eventArray
    
  }
}
saveRDS(eventArrayList, paste0(wd,"/Code/InterimData/eventArrayList.rds"))
nc.close(ncin)