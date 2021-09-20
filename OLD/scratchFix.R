#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-02-21
#
# Event threshold extraction from daily flow for various thresholds.
# A GPA distribution is fitted to the peaks over these thresholds.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-02-21, edited 2020-05-04
# Pipeline version 2020-09-07
#
# OUTPUTS: threshDayExcList.rda,
#          thresMat.rda,
#          paramtable.csv
#
#~~~~~~~~~~~~~~~~~~~~~~~

#### SETUP ####----------------------
commandArgs <- function(...){"01"}
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~ /AQUACAT/CodeABG/setup_script_00.R")
}

library(ilaprosUtils)
library(extRemes)
library(lmomco)

thresh1 <- "POT2"
ws1 <- "pc05"
# Only using POT2 and 2% inundation minimums.
jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)

print(paste0(data_wd, subfold, "threshMat_RCM",
                           RCM, suffix, ".rds"))

threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                           RCM, suffix, ".rds"))

threshDayExcList <- readRDS(paste0(data_wd, subfold,
                            "threshDayExcList_RCM", RCM, suffix, ".rds"))

ST <-  Sys.time()
ncin <- nc_open(ncoriginal, readunlim=FALSE)
  # this is a huge file, do not open without reason.
print(Sys.time() - ST)
print(ncin)
NH <- nrow(rn)
ND <- 10800

dvec <- c()
for(n in 1:ND){
  if(n < 11 | n%%50==0){
    print(n)
  }
    tSlice <- ncvar_get(ncin, varid="dmflow",
                         start=c(1, 1,  n),
                         count=c(-1, -1, 1))
    
    # find quantile for threshold
  vvec <- c()
  for(k in 1:NH){
    vvec[k] <- tSlice[rn$row[k], rn$col[k]]
  }
  dvec[n] <- sum(vvec > threshMat[,jT])
}

eventLList <- vector("list", 5)
names(eventLList) <- threshName
eventDayList <- vector("list", 5)
names(eventDayList) <- threshName

for(j in jT){
  print(j)
  eventLList[[j]] <- vector("list", 5)
  names(eventLList[[j]]) <- wsName
  eventDayList[[j]] <- vector("list", 5)
  names(eventDayList[[j]]) <- wsName
  for(k in jW){
    print(k)
    eventGo <- 0
    eventL <- c()
    eventD <- c()
    eventSt <- 1
    for(i in 1:ND){
      if (dvec[i] > (NH*0.005)) {
        eventGo <- eventGo + 1
      }else{
        if (eventGo > 0){
          eventD <- c(eventD, eventSt)
          eventL <- c(eventL, eventGo)
          #print(eventGo)
        }  
        eventGo <- 0
        eventSt <- i+1
      }
    }
    md <- apply(cbind(eventL, eventD), 1, function(x){max(dvec[x[2]+(0:(x[1]-1))])})
    xstart <- apply(df, 1, function(x){max(x[2], x[2]+which(dvec[x[2] + (0:(x[1]-1))]==x[3])-7)})
    eventD0 <- eventD
    eventL0 <- eventL
    eventD[eventL > 14] <- xstart[eventL > 14]
    eventL[eventL > 14] <- 14
    eventLList[[j]][[k]] <- eventL
    eventDayList[[j]][[k]] <- eventD
  }
}
VVV <- c()
save(VVV, eventLList, eventDayList, dvec, file=paste0(data_wd,"trials.RDa"))

inunMat <- matrix(0, nrow=NT, ncol=ND)
#inunDays <- matrix(0, nrow=NT, ncol=NH)

for(j in jT){  # for each threshold value j
  #print(threshName[j])
  for(k in 1:NH){ # at each river network gridcell k
    tde_jk <- threshDayExcList[[j]][[k]]
    
    for(n in 1:length(tde_jk)){
      tde <- tde_jk[n] # which days was k exceeded past j (tde \subset 1:NH)
      inunMat[j,tde] <- inunMat[j,tde] + 1 # add 1 to all those locations.
    }
  }
}

VVV <- c()
for(i in 1:length(eventLList[[jT]][[jW]])){
  L <- eventLList[[jT]][[jW]][i]
  D <- eventDayList[[jT]][[jW]][i]
 VVV[i] <- (max(inunMat[jT,D:(D+L-1)]))
}

save(VVV, eventLList, eventDayList, dvec, file=paste0(data_wd,"trials.RDa"))



####################################################
commandArgs <- function(...){"01"}
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~ /AQUACAT/CodeABG/setup_script_00.R")
}

thresh1 <- "POT2"
ws1 <- "pc05"
# Only using POT2 and 2% inundation minimums.
jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)

load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))

eventDayList0 <- eventDayList
eventLList0 <- eventLList

load("S:/Data/trials.RDa")

dvt <- dvec > NH*0.005
RL <- rle(dvt)
RL1 <- RL$lengths[RL$values]
wdvt <- which(dvt[-length(dvt)] < dvt[-1]) + 1
df <- data.frame(RL1, wdvt)
dvec[6156+(-1:52)]
df$md <- apply(df, 1, function(x){max(dvec[x[2]+(0:(x[1]-1))])})
plot(md)
df$xstart <- apply(df, 1, function(x){max(x[2], x[2]+which(dvec[x[2] + (0:(x[1]-1))]==x[3])-7)})