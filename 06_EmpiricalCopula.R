#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-15
#
# The Empirical Copula functions
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-15
#
#~~~~~~~~~~~~~~~~~~~~~~~

##### SETUP #####---------------------------------------------------

library(ncdf4)
library(extRemes)
library(dplyr)
library(readr)
library(fitdistrplus)

if (substr(osVersion,1,3) == "Win") {
  ncname <- "S:/CodeABG/InterimData/dmflow_copy.nc"  # rechunked for spaceslices
  ncname2 <- "S:/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc" 
    # pre-rechunk
  wd <- "S:/"
  wd_id <- "S:/CodeABG/InterimData/"
  wd_cd <- "S:/CodeABG/"
  
} else {
  ncname <- "/prj/aquacat/CodeABG/InterimData/dmflow_copy.nc"
  ncname2 <- "S:/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc"
  wd <- "/prj/aquacat/"
  wd_id <- "/prj/aquacat/CodeABG/InterimData/"
  wd_cd <- "/prj/aquacat/CodeABG/"
}

### DATA --------------------------

# print(ST <- Sys.time())
# ncin <- nc_open(ncname) # This file is ~2.5GB on the linux server.
# print(ncin)
# print(floor(Sys.time() - ST))

ND <- 10800 # Number of days


# river network
rn <- read_csv(paste0(wd_id, "hasData2.csv"))
NH <- nrow(rn)


# threshold for inundation
threshVal <- c(5/365, 2/365, 1/365, 0.2/365, 0.1/365)
threshName <- c("POT5", "POT2", "POT1", "Q5", "Q10")
NT <- length(threshVal)


# cut-off bound for widespread event
wsBound <- c(0.05, 0.02, 0.01, 0.005, 0.001)
wsName <- c("pc5", "pc2", "pc1", "pc05", "pc01")
NW <- length(wsBound)


# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS(paste0(wd_id, "threshDayExcList2.rds"))



# matrix of threshold value (col) at a given cell (row)
threshMat <- read_csv(paste0(wd_id, "threshMat2.csv"))
#dim(threshMat)  =  19914 x 5



# eventLList length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(wd_id, "eventLists03.RDa")) 
NE <- length(eventDayList[[2]][[4]]) # POT2, 0.5% inun.



# timewise maxima at each cell for each event ((NE + 2) x NH)
eventDF <- readr::read_csv(paste0(wd,"Data/eventdf_POT2_pc05.csv"))

# PoE under different computations with extra data. Tidy format.
present <- readr::read_csv(paste0(wd,"Data/present2_returnlevels.csv"))

dummy <- readr::read_csv(paste0(wd, "Data/dummy1_returnlevels.csv"))
dummy <- present %>% filter(loc < 4)


##### COPULA FUNCTIONS #####---------------------------------------------------

# Get event set

# select event

tailsGen <- function(n, low=1, pool=present$gpp, maxit=1000, betapar=c(1,1)){
  # draw new probabilities of exceedence from pool, allowing minumum lower bound
  # to be enforced. Technically drawing from the conditional distribution 
  # P[X=x | X_i < low for some i].
  #
  # n         number of PoEs to draw
  # low       lower target on PoE
  # pool      existing PoE to take quantiles from.
  # maxit     number of trials to get an event of the required probability
  # betapar   vector of two parameters for the Beta distribution
  
  QV <- rep(2, n)
  
  #TODO needs the right distribution to extrapolate from pool to reach tails.
  
  z <- 0
  while(min(QV) > low){
    V <- runif(n)
    z <- z+1
    
    QV <- qbeta(V, betapar[1], betapar[2])
    if(z >= maxit){
      print("maxit exceeded")
      break
    }
  }
  if(z>1){print(paste(z, "iterations."))}
  unname(QV)
  
 # print(QV)
  
  return(unname(QV))
  
  
}




generateNewEvent <- function(eventSet=present, NE=285, NH=19914,
                             maxit=10000, low=1, pool=present$gpp, betapar=c(1,1)){
  # generates a new widespread event based on a given event library
  #
  # eventSet  object like present1_return_levels.csv, needs the right col names
  # NE        number of events in the set
  # NH        number of gridpoints involved (~20000 in standard set)
  
  U <- sample.int(NE, 1)
  
  #print(paste(U,NE))
  
  eventSubset <-  eventSet %>% filter(eventNo == U)
  
  #print(paste("gpp",eventSubset$gpp))
  rankNew <- rank(eventSubset$gpp, ties.method="random")
  
  eventMags <- tailsGen(NH, low=low, pool=pool, maxit=maxit, betapar=betapar)
  
 #print(paste("eventMags",length(eventMags)))
  
  # add the magnitudes according to rank
  magsNew <- sort(eventMags)[rankNew]
  
  
  #print(paste("magsNew",length(magsNew)))
  
    #print(paste("rankNew",rankNew))
  
  df <- eventSubset[,c("loc", "Northing", "Easting")]
  #print(paste("eventSubset",dim(eventSubset), collapse=" "))
  #print(paste("df", dim(df)))
  df$rank <- rankNew
  df$gpp <- magsNew
  df$low <- rep(low, nrow(df))
  
  df
  
  
}

print("get started")

FF <- fitdist(dummy$gpp, "beta", method="mle")

T3 <- tailsGen(1, pool=dummy$gpp, maxit=1000, betapar=FF$estimate)

M3 <- generateNewEvent(eventSet=dummy, NE=NE, NH=3,
                       maxit=100, pool=dummy$gpp,
                       low=0.01, betapar=FF$estimate)

print(M3)


# Get ranking according to PoE per event
NH <- 3
M <- 1000 # number of new events
newEventMat <- matrix(NA, nrow=NH, ncol=M)

print("dummyset")

for(m in 1:M){
  
  newEventMat[,m] <- generateNewEvent(eventSet=dummy, NE=NE, NH=3,
                                      maxit=100, pool=dummy$gpp,
                                      low=1, betapar=FF$estimate)$gpp
   
}

readr::write_csv(data.frame(newEventMat),
                 paste0(wd,"Data/NewEventTrialEC_POT2_pc05.csv"))

FF <- fitdist(present$gpp, "beta", method="mle")

print("fullset")
NH <- nrow(rn)
newEventMat <- matrix(NA, nrow=NH, ncol=M)
for(m in 1:M){
  #if(m %% 50 == 0){print(m)}
  newEventMat[,m] <- generateNewEvent(eventSet=present, NE=NE, NH=nrow(rn),
                                      maxit=100, pool=present$gpp,
                                      low=1, betapar=FF$estimate)$gpp
}

readr::write_csv(data.frame(newEventMat),
                 paste0(wd,"Data/NewEventPresentEC_POT2_pc05.csv"))
print(Sys.time())
# Roll new tails

# Apply to ordered event

# save as new event