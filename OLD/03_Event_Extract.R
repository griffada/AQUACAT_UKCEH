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
library(fields)
library(rgdal)

threshVal <- c(5/365, 2/365, 1/365, 0.2/365, 0.1/365)
threshName <- c("POT5", "POT2", "POT1", "Q5", "Q10")
NT <- length(threshVal)


##### DATA #####------------------------------------------------------------
if(substr(osVersion,1,3) == "Win"){
  wd <- "S:/"
}else{
  wd <- "/prj/aquacat/"
}

ncname <- "run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc" #2GB

ND <- 10800 # Number of days


rn <- read.csv(paste0(wd,"CodeABG/InterimData/hasData.csv"),
               stringsAsFactors=FALSE)
NH <- nrow(rn)


threshDayExcList <- readRDS(paste0(wd,
                            "CodeABG/InterimData/threshDayExcList2.rds"))
#5 lists of NH lists

# Number of events
S <- sapply(threshDayExcList, function(l){length(l[[1]])})
for(i in 1:5){
  print(paste(threshName[i], "events per grid-cell:", S[i]))
}

thresMat <- readRDS(paste0(wd,"CodeABG/InterimData/threshMat2.rds"))


### Find days with most exceedences ###-----------------------------------

# Number of sites exceeded for each day at each threshold
inunMat <- matrix(0, nrow=NT, ncol=ND)
inunDays <- matrix(0, nrow=NT, ncol=NH)

for(j in 1:NT){  # for each threshold value
  #print(threshName[j])
  for(k in 1:NH){ # at each river network gridcell
    tde_jk <- threshDayExcList[[j]][[k]]
    
    for(n in 1:length(tde_jk)){
      tde <- tde_jk[n] # which days was it exceeded
      inunMat[j,tde] <- inunMat[j,tde] + 1
      
    }
    
  }
}

rownames(inunMat) <- threshName

EventSizeSumm <- apply(inunMat, 1, function(v){
  quantile(v[v>0], probs=
             c(0.99,0.9,0.8,0.5,0.4, 0.3,0.2,0.1,0.01))
})

### Extract data at timepoints for HT/EC ###-----------------------------

# Percentage of inundated cells
wsBound <- c(0.05, 0.02, 0.01, 0.005, 0.001)
wsName <- c("pc5", "pc2", "pc1", "pc05", "pc01")
NW <- length(wsBound)

widespreadArray <- array(FALSE, dim=c(NT, ND, NW))
for(j in 1:NW){
  # Is there a widespread event on this day (at bound j)?
  widespreadArray[,,j] <- inunMat >= NH*wsBound[j] 
}
dimnames(widespreadArray) <- list(threshold = threshName,
                                  days=1:ND,
                                  inundation = c("pc5", "pc2", "pc1",
                                                 "pc05", "pc01"))

widespreadCount <- apply(widespreadArray, c(1,3), sum)

##### event lengths and event starts #####-------------------------------

eventLList <- vector("list",5)
names(eventLList) <- threshName
eventDayList <- vector("list", 5)
names(eventDayList) <- threshName

for(j in 1:5){
  eventLList[[j]] <- vector("list",5)
  names(eventLList[[j]]) <- wsName
  eventDayList[[j]] <- vector("list", 5)
  names(eventDayList[[j]]) <- wsName
  for(k in 1:5){
    eventGo <- 0
    eventL <- c()
    eventD <- c()
    for(i in 1:ND){
      if (widespreadArray[j, i, k]) {
        eventGo <- eventGo + 1
      }else{
        if (eventGo > 0){
          eventD <- c(eventD, i)
          eventL <- c(eventL, eventGo)
          #print(eventGo)
        }  
        eventGo <- 0
      }
    }
    eventLList[[j]][[k]] <- eventL
    eventDayList[[j]][[k]] <- eventD
  }
}
save(eventLList, eventDayList,
     file=paste0(wd,"CodeABG/InterimData/eventLists03.RDa"))
##### PLOTTING EXAMPLES #####--------------------------------------------------
load(paste0(wd, "CodeABG/InterimData/eventLists03.RDa"))

qv <- c()
quan <-  quantile(inunMat[4,inunMat[4,]>0],
                  probs=c(0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95))
for(qi in quan){
  qv <- c(qv,which.min(abs(inunMat[4,] - qi)))
}