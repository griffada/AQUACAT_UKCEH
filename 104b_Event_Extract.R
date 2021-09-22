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

if(interactive()){commandArgs <- function(...){c("01","present")}}

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

if(settings$eventlist){
   stop("eventLists alread exist for 104. Proceeding to next job.")
}

##### DATA #####------------------------------------------------------------


threshDayExcList <- readRDS(paste0(data_wd, subfold,
                            "threshDayExcList_RCM", RCM, suffix, ".rds"))
#5 lists of NH lists, one for each threshold

threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                           RCM, suffix, ".rds"))
# thresh1 <- "POT2"
# ws1 <- "pc05"
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)


ST   <-  Sys.time()
ncin <- nc_open(ncoriginal, readunlim=FALSE)
  # this is a huge file, do not open without reason.
print(Sys.time() - ST)
#print(ncin)

##### Find days with most exceedences #####--------------------------------
  ST0 <- Sys.time()
  ST <- Sys.time()
# Number of sites exceeded for each day at each threshold
dvec <- matrix(NA, ncol=5, nrow=ND)
for (i in 1:ND) {
  if ((i < 10) | (i %% 200 == 0)) { # time recording
    print(paste(i, "out of", ND))
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((ND - i)/ND, 2)))
    print(paste("Time remaining", round((ND - i)/i * I0, 2)))
    print(paste("Since last readout:", round(I, 2)))
    ST <- Sys.time()
    print(ST) 
  }
  tSlice <- ncvar_get(ncin, varid="dmflow",
                       start=c(1, 1,  i),
                       count=c(-1, -1, 1))
  
  # find quantile for threshold
  vvec <- c()

  for (k in 1:NH) vvec[k] <- tSlice[rn$row[k], rn$col[k]]

  for (j in 1:5) dvec[i,j] <- sum(vvec > threshMat[,j])
}

saveRDS(dvec, paste0(data_wd, subfold, "dvec_RCM", RCM, suffix, ".RDs"))
dvec <- readRDS(paste0(data_wd, subfold, "dvec_RCM", RCM, suffix, ".RDs"))

##### Event lengths and event starts #####-------------------------------

eventLList   <- vector("list", 5)
names(eventLList)   <- threshName
eventDayList <- vector("list", 5)
names(eventDayList) <- threshName
maxInunList  <- vector("list", 5)
names(maxInunList)  <- threshName

for (j in jT) {
  eventLList[[j]]          <- vector("list", 5)
  names(eventLList[[j]])   <- wsName
  eventDayList[[j]]        <- vector("list", 5)
  names(eventDayList[[j]]) <- wsName
  maxInunList[[j]] <- vector("list", 5)
  names(maxInunList[[j]]) <- wsName
  for (k in jW) {
    eventGo <- 0
    eventSt <- 1
    eventL  <- c()
    eventD  <- c()
    for (i in 1:ND) {
      if (dvec[i,j] > NH*wsCutoff[k]) { #if the event is still going
        eventGo <- eventGo + 1 
      }else{ #count up the length of the event
        if (eventGo > 0) {
          eventD <- c(eventD, eventSt)
          eventL <- c(eventL, eventGo)
        }  
        # then reset and continue
        eventGo <- 0
        eventSt <- i+1
      }
    }
    maxInun <- apply(cbind(eventL, eventD),
                     1,
                     function(x){max(dvec[x[2] + (0:(x[1] - 1)), j])})
    xstart <- apply(do.call(cbind, list(eventL, eventD, maxInun)), 1,
       function(x){max(x[2], x[2] + which(dvec[x[2] + (0:(x[1] - 1))]==x[3]) - 7)})
    eventD[eventL > 14] <- xstart[eventL > 14]
    eventL[eventL > 14] <- 14
    eventLList[[j]][[k]]   <- eventL
    eventDayList[[j]][[k]] <- eventD
    maxInunList[[j]][[k]]  <- maxInun
  }
}

#season <- function(D) switch(((D %/% 90) %% 4) + 1, "DJF", "MAM", "JJA", "SON")
season2 <- function(D) c("DJF","MAM","JJA","SON")[((D %/% 90) %% 4) + 1]

# Initial summary of events, possibly to be used in summaries.
summInitial <- data.frame(Start=eventDayList[[jT]][[jW]],
                          Length=eventLList[[jT]][[jW]],
                          MaxInun=maxInunList[[jT]][[jW]],
                          Season=season2(eventDayList[[jT]][[jW]]))
                          

#### OUTPUTS ##### ----------------------------------------------------
print("saving outputs")
save(eventLList,
     eventDayList,
     dvec,
     maxInunList,
     file=paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))
readr::write_csv(summInitial,
          paste0(data_wd, subfold, "initialSummary_RCM", RCM, suffix, ".csv"))
settings$paramtable <- TRUE
settings$eventlist <- TRUE
settings$event_extract <- "104b"
write_yaml(settings, settingspath)

print(Sys.time())
print("104b Done")
