######
# Adam Griffin, 2020-04-27
#
# Event extraction from daily flow, no test for independence. Using outputs
# from 102_Threshold_Extract.
# 
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-04-27
# Pipeline version ABG 2020-09-07
# 
# Outputs:
#   eventLists***.Rda = eventLList, eventDayList
#
#####

##### SETUP #####------------------------------------------------------------


if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

##### DATA #####------------------------------------------------------------


threshDayExcList <- readRDS(paste0(data_wd, subfold,
                            "threshDayExcList_RCM", RCM, suffix, ".rds"))
#5 lists of NH lists

thresMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                           RCM, suffix, ".rds"))


# Number of events
S <- sapply(threshDayExcList, function(l){length(l[[1]])})
for(i in 1:5){
  print(paste(threshName[i], "events per grid-cell:", S[i]))
}



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

# Percentage of inundated cell

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

#### OUTPUTS ----------------------------------------------------
save(eventLList, eventDayList,
     file=paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))
