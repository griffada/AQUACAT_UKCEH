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
# OUTPUTS: eventLists***.Rda = eventLList, eventDayList
#
#####

##### SETUP #####------------------------------------------------------------


if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

if(file.exists(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))){
  stop("eventLists alread exist for 104. Proceeding to next job.")
}else{
  print("Proceeding to job.")
}

##### DATA #####------------------------------------------------------------


threshDayExcList <- readRDS(paste0(data_wd, subfold,
                            "threshDayExcList_RCM", RCM, suffix, ".rds"))
#5 lists of NH lists, one for each threshold

threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                           RCM, suffix, ".rds"))


# Number of events
S <- sapply(threshDayExcList, function(l){length(l[[1]])})
for(i in 1:5){
  print(paste(threshName[i], "events per grid-cell:", S[i]))
}

### Find days with most exceedences ###-----------------------------------

# Number of sites exceeded for each day at each threshold
inunMat <- matrix(0, nrow=NT, ncol=ND)
#inunDays <- matrix(0, nrow=NT, ncol=NH)

for(j in 1:NT){  # for each threshold value j
  #print(threshName[j])
  for(k in 1:NH){ # at each river network gridcell k
    tde_jk <- threshDayExcList[[j]][[k]]
    
    for(n in 1:length(tde_jk)){
      tde <- tde_jk[n] # which days was k exceeded past j (tde \subset 1:NH)
      inunMat[j,tde] <- inunMat[j,tde] + 1 # add 1 to all those locations.
    }
  }
}

rownames(inunMat) <- threshName

### Extract data at timepoints for HT/EC ###-----------------------------

widespreadArray <- array(FALSE, dim=c(NT, ND, NW))
for(j in 1:NW){
  for(k in 1:NT){
  # Is there a widespread event on this day (at bound j)?
  widespreadArray[k, , j] <- inunMat[k,] >= (NH*wsCutoff[j]) 
  }
}
dimnames(widespreadArray) <- list(threshold = threshName,
                                  days = 1:ND,
                                  inundation = c("pc5", "pc2", "pc1",
                                                 "pc05", "pc01"))

widespreadCount <- apply(widespreadArray, c(1, 3), sum)

##### event lengths and event starts #####-------------------------------

eventLList <- vector("list", 5)
names(eventLList) <- threshName
eventDayList <- vector("list", 5)
names(eventDayList) <- threshName

for (j in 1:5) {
  print(j)
  eventLList[[j]]          <- vector("list", 5)
  names(eventLList[[j]])   <- wsName
  eventDayList[[j]]        <- vector("list", 5)
  names(eventDayList[[j]]) <- wsName
  for (k in 1:5) {
    print(k)
    eventGo <- 0
    eventSt <- 1
    eventL  <- c()
    eventD  <- c()
    for (i in 1:ND) {
      if (widespreadArray[j, i, k] & (eventSt < 15)) { #if the event is still going
        eventGo <- eventGo + 1 
      }else{ #count up the length of the event
        if (eventGo > 0) {
          eventD <- c(eventD, eventSt)
          eventL <- c(eventL, eventGo)
          #print(eventGo)
        }  
        # then reset and continue
        eventGo <- 0
        eventSt <- i+1
      }
    }
    eventLList[[j]][[k]] <- eventL
    eventDayList[[j]][[k]] <- eventD
  }
}

#### OUTPUTS ----------------------------------------------------
save(eventLList, eventDayList,
     file=paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))
