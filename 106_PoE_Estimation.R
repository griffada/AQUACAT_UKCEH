#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-15
#
# Estimating probability of exceedence along time series using ecdf and
# plotting positions.
# Uses inputs from 104_Event_Extract and 105_timeDF_compile.

# Note that for FUTURE estimation of return periods, make use of the present day
# data (tSlice) to compute the ecdf.
#
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-15
# Pipeline version ABG 2020-09-07
#
# Outputs: present_returnlevels***.csv
#
#~~~~~~~~~~~~~~~~~~~~~~~

##### SETUP #####------------------------------------------------------------
library(extRemes)

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

thresh1 <- "POT2"
ws1 <- "pc05"
print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

if(file.exists(paste0(data_wd, subfold, "returnlevels_",
                      thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))){
  stop("returnlevels already exists. ending 106.")
}




subfold <- paste0("RCM", RCM, suffix, "/")

jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)
### Functions ###---------------------------------------

logit <- function(x){log(x/(1-x))}
invlogit <- function(y){1/(1 + exp(-1*y))}

gringorten <- function(v){
  ((length(v) + 1 - rank(v)) - 0.44)/(length(v) + 0.12)
}

weibull <- function(v){
  (length(v) + 1 - rank(v))/(length(v) + 1)
}

##### DATA #####--------------------------------------------------------------

print(ST <- Sys.time())

suffix_pres <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")


ncin <- nc_open(paste0(data_wd, subfold_pres, "dmflow_copy_RCM",
                       RCM, suffix_pres, ".nc")) # This file is ~2.5GB on the linux server.
print(ncin)
print(floor(Sys.time() - ST))


# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS( paste0(data_wd, subfold, "threshDayExcList_RCM",
                                    RCM, suffix,".rds"))


# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                            RCM, suffix, ".rds"))
#dim(threshMat) #19914 x 5


#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))



# timewise maxima at each cell for each event (NE x NH)

#library(data.table)
#edf <- fread(paste0(wd,"/Data/eventdf_POT2_pc2.csv"), colClasses=rep("numeric",287))

eventDF <- readr::read_csv(paste0(data_wd,subfold, "eventdf_",thresh1,"_", ws1,
                                  "_RCM", RCM, suffix, ".csv"),
                           col_types=cols(
                             .default = col_double()
                           ))

NE <- ncol(eventDF) - 4


##### PROB CALCULATION #####-------------------------------------------------


### prealloc -----------------------
rarityDF <- expand.grid("eventNo" = 1:NE,
                        "loc" = 1:NH)

rarityDF$Easting <- rn[rarityDF[, 2], 1]
rarityDF$Northing <- rn[rarityDF[, 2], 2]
rarityDF$thresh <- NA
rarityDF$DayS <- eventDayList[[jT]][[jW]][rarityDF[, 1]]
rarityDF$val <- NA
rarityDF$gpp <- NA
rarityDF$ecdf <- NA
rarityDF$gev <- NA


n <- 1

edl <- eventDayList[[jT]][[jW]]

ST0 <- proc.time()
ST <- proc.time()
print("loop start")
for(n in 1:NH){
  if(n %% 200 == 0){
    print(paste(n, "out of", NH))
  }
  
  
  i <- rn[n,1]
  j <- rn[n,2]
  #!#!#!#!# Pull out spaceslice from present day to get present day PoE
  tSlice <- ncvar_get(ncin, varid="dmflow",
                      start=c(i, j,  1),
                      count=c(1, 1, -1))
  
  #event magnitudes from correct period for a specific location(present or future)
  tSliceEvent <- unlist(eventDF[n, -(1:4)], use.names=FALSE)
  
  threshval <- threshMat[n, jT] # POT2 column
  
  
  
  rarityDF$thresh[which(rarityDF$loc == n)] <- threshval
  rarityDF$val[which(rarityDF$loc == n)] <- tSliceEvent

  # get ecdf and estimate PoE
  
  ecdfSlice <- ecdf(tSlice)
  poeEvent <- 1 - ecdfSlice(tSliceEvent)
  
  
  # Weibull or Gringorten plotting position
  
  
  if(period=="future"){
    grEvent <- sapply(tSliceEvent, function(x){ gringorten(c(x, tSlice))[1] })
  }else{
    grSlice <- gringorten(tSlice)
    grEvent <- sapply(tSliceEvent, function(x){ grSlice[which(tSlice == x)[1]] })
  }
  
  rarityDF$gpp[which(rarityDF$loc == n)] <- grEvent
  rarityDF$ecdf[which(rarityDF$loc == n)] <- poeEvent

  # GEV fitted to whole spaceslice
  FFGEV <- fevd(x=tSlice,
                type='GEV')$results$par
  QFGEV <- 1- pevd(q=tSliceEvent,
                   loc=FFGEV[1],
                   scale=FFGEV[2],
                   shape=FFGEV[3],
                   type='GEV')
  
  rarityDF$gev[which(rarityDF$loc==n)] <- QFGEV
  
  if(n %% 200 == 0){
    print(floor(proc.time() - ST0)[1:3])
    print(floor(proc.time() -  ST)[1:3])
  }
  ST <- proc.time()
}



readr::write_csv(x=rarityDF,
                 path=paste0(data_wd, subfold, "returnlevels_",
                             thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
nc_close(ncin)
print(Sys.time())
######
#
# CONVERTION FROM PoE IN DAYS (p) TO PoE IN YEARS (b): p = 1 - (1-b)^(1/365.25)
#
# b = 1- (1-p)^(365.25)
#
#####

