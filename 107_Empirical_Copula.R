#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-15
#
# The Empirical Copula functions for simulating new events from given data.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-15
# Pipeline version ABG 2020-09-07
#
# Outputs: NewEventPresentEC_***.csv: data table of events, one event per row.
#   
#
#~~~~~~~~~~~~~~~~~~~~~~~

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

##### SETUP #####---------------------------------------------------

library(extRemes)
library(dplyr)
library(fitdistrplus)

### DATA --------------------------

thresh1 <- "POT2" #!#!#!#!# Important constants to select.
ws1 <- "pc05"
print("Running for threshold", POT2, "at ", ws1, "minimum spread.")

jV <- which(threshName==thresh1)
jI <- which(wsName == ws1)

# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS( paste0(data_wd, subfold, "threshDayExcList_RCM",
                                    RCM, suffix,".rds"))


# matrix of threshold value (col) at a given cell (row)
threshMat <- read.csv(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"),
                      stringsAsFactors=FALSE)
#dim(threshMat) #19914 x 5


#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))

# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))

NE <- length(eventDayList[[jV]][[jI]]) # POT2, 0.5% inun.




# timewise maxima at each cell for each event ((NE + 2) x NH)
eventDF <- readr::read_csv(paste0(data_wd,subfold, "eventdf_",
                                  thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

# PoE under different computations with extra data. Tidy format.
present <- readr::read_csv(paste0(data_wd, subfold, "returnlevels_",
                                  thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

#dummy <- present %>% filter(loc < 4)


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


print("Determining tail distribution")

FF <- fitdist(present$gpp, "beta", method="mle")

T3 <- tailsGen(1, pool=present$gpp, maxit=1000, betapar=FF$estimate)

# Number of new events to simulate
M <- 500

print("Simulating new events")

newEventMat <- matrix(NA, nrow=NH, ncol=M)
for(m in 1:M){
  #if(m %% 50 == 0){print(m)}
  newEventMat[,m] <- generateNewEvent(eventSet=present, NE=NE, NH=nrow(rn),
                                      maxit=100, pool=present$gpp,
                                      low=1, betapar=FF$estimate)$gpp
}

readr::write_csv(data.frame(newEventMat),
                 paste0(data_wd, subfold, "NewEventEC_",
                            thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
print(Sys.time())
# Roll new tails

# Apply to ordered event

# save as new event