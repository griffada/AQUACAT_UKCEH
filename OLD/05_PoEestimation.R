#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-15
#
# Estimating PoE along time series using ecdf and plotting positions.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-15
#
#~~~~~~~~~~~~~~~~~~~~~~~

##### SETUP #####------------------------------------------------------------

print(Sys.time())
library(ncdf4)
library(raster)
library(fields)
library(extRemes)

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

if(substr(osVersion,1,3) == "Win"){
  ncname <- "S:/CodeABG/InterimData/dmflow_copy.nc"
  ncname2 <- "S:/CodeABG/InterimData/dmflow_trial.nc"
  wd <- "S:/"
}else{
  ncname <- "/prj/aquacat/CodeABG/InterimData/dmflow_copy.nc"
  wd <- "/prj/aquacat/"
}

ND <- 10800 # Number of days

# river network
rn <- read.csv(paste0(wd,"CodeABG/InterimData/hasData2.csv"),
               stringsAsFactors=FALSE)
NH <- nrow(rn)

# threshold for inundation
threshVal <- c(5/365, 2/365, 1/365, 0.2/365, 0.1/365)
threshName <- c("POT5", "POT2", "POT1", "Q5", "Q10")
NT <- length(threshVal)

# cut-off bound for widespread event
wsBound <- c(0.05, 0.02, 0.01, 0.005, 0.001)
wsName <- c("pc5", "pc2", "pc1", "pc05", "pc01")
NW <- length(wsBound)




print(ST <- Sys.time())
ncin <- nc_open(ncname) # This file is ~2.5GB on the linux server.
print(ncin)
print(floor(Sys.time() - ST))




# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS(
                      paste0(wd, "CodeABG/InterimData/threshDayExcList2.rds"))



# matrix of threshold value (col) at a given cell (row)
threshMat <- read.csv(paste0(wd, "CodeABG/InterimData/threshMat2.csv"),
                      stringsAsFactors=FALSE)
dim(threshMat) #19914 x 5



# eventLList length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(wd, "CodeABG/InterimData/eventLists03.RDa")) 
NE <- length(eventDayList[[2]][[4]])



# timewise maxima at each cell for each event (NE x NH)

#library(data.table)
#edf <- fread(paste0(wd,"/Data/eventdf_POT2_pc2.csv"), colClasses=rep("numeric",287))

eventDF <- readr::read_csv(paste0(wd,"Data/eventdf_POT2_pc05.csv"))




##### PROB CALCULATION #####-------------------------------------------------


### prealloc -----------------------
rarityDF <- expand.grid("eventNo" = 1:NE,
                        "loc" = 1:NH)

rarityDF$Easting <- rn[rarityDF[, 2], 1]
rarityDF$Northing <- rn[rarityDF[, 2], 2]
rarityDF$thresh <- NA
rarityDF$DayS <- eventDayList[[2]][[2]][rarityDF[, 1]]
rarityDF$val <- NA
rarityDF$gpp <- NA
rarityDF$ecdf <- NA
rarityDF$gev <- NA


n <- 1

edl <- eventDayList[[2]][[2]]

ST0 <- proc.time()
ST <- proc.time()
print("loop start")
for(n in 1:NH){
  if(n %% 200 == 0){
    print(paste(n, "out of", NH))
  }
  
  
  i <- rn[n,1]
  j <- rn[n,2]
  # Pull out spaceslice
  tSlice <- ncvar_get(ncin, varid="dmflow",
                      start=c(i, j,  1),
                      count=c(1, 1, -1))
  
  tSliceEvent <- unname(unlist(eventDF[n, -(1:2)]))
  
  threshval <- threshMat[n, 2] # POT2 column
  
  
  
  rarityDF$thresh[which(rarityDF$loc == n)] <- threshval
  rarityDF$val[which(rarityDF$loc == n)] <- tSliceEvent
  
  # get ecdf and estimate PoE
  
  ecdfSlice <- ecdf(tSlice)
  poeEvent <- 1 - ecdfSlice(tSliceEvent)
  
  
  # Weibull or Gringorten plotting position
  
  grSlice <- gringorten(tSlice)
  grEvent <- sapply(tSliceEvent,
                            function(x){grSlice[which(tSlice == x)[1]]})
  
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
                 path=paste0(wd,"Data/present_returnlevels_POT2_pc05.csv"))
nc_close(ncin)
print(Sys.time())
######
#
# CONVERTION FROM PoE IN DAYS (p) TO PoE IN YEARS (b): p = 1 - (1-b)^(1/365.25)
#
# b = 1 - (1-p)^(365.25)
#
#####

