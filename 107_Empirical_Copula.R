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
# OUTPUTS: NewEventPresentEC_***.csv: data table of events, one event per row.
#
#~~~~~~~~~~~~~~~~~~~~~~~

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

##### SETUP #####---------------------------------------------------

library(extRemes)
library(dplyr)
library(fitdistrplus)

### DATA --------------------------

thresh1 <- "POT2" #!#!#!#!# Important constants to select.
ws1 <- "pc05"

print(paste("Running for threshold", thresh1, "at", ws1, "minimum spread."))

suffix_pres <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix_pres, "_out.nc") 

jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)

# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS(paste0(data_wd, subfold, "threshDayExcList_RCM",
                                    RCM, suffix,".rds"))

# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"))
#dim(threshMat) #19914 x 5

#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))

NE <- length(eventDayList[[jT]][[jW]]) # POT2, 0.5% inun.

NH <- nrow(rn)

# timewise maxima at each cell for each event ((NE + 2) x NH)
obs_events  <- readr::read_csv(paste0(data_wd,subfold, "eventflow_OBS_",thresh1,
                                      "_", ws1, "_RCM", RCM, suffix, ".csv"),
                               col_types=cols(.default = col_double()))[,-(1:4)]

obs_dpoe  <- readr::read_csv(paste0(data_wd,subfold, "eventdpe_OBS_",thresh1,
                                    "_", ws1, "_RCM", RCM, suffix, ".csv"),
                             col_types=cols(.default = col_double()))[,-(1:4)]

obs_apoe  <- readr::read_csv(paste0(data_wd,subfold, "eventape_OBS_",thresh1,
                                    "_", ws1, "_RCM", RCM, suffix, ".csv"),
                             col_types=cols(.default = col_double()))[,-(1:4)]

partable <- readr::read_csv(paste0(data_wd,subfold, 
                           "paramtable_",thresh1, "_RCM", RCM, suffix, ".csv"),
                            col_types=cols(.default= col_double()))

partable_pres <- readr::read_csv(paste0(data_wd,subfold_pres, 
                        "paramtable_",thresh1, "_RCM", RCM, suffix_pres, ".csv"),
                            col_types=cols(.default= col_double()))

colnames(partable)[1] <- "meanint"





### COPULA FUNCTIONS ###---------------------------------------------------

tailsGen <- function(n, low=1, maxit=1000, betapar=c(1,1)){
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
  return(unname(QV))
}

generateNewEvent <- function(eventSet=obs_dpoe, NE=285, NH=19914,
                             maxit=10000, low=1, betapar=c(1,1)){
  # generates a new widespread event based on a given event library
  #
  # eventSet  grid of daily PoE, one location per row, one event per col.
  # NE        number of events in the set
  # NH        number of gridpoints involved (~20000 in standard set)
  
  U <- sample.int(NE, 1)
  eventSubset <-  eventSet[,U]
  rankNew <- rank(eventSubset, ties.method="random")
  eventMags <- tailsGen(NH, low=low, maxit=maxit, betapar=betapar)
  # add the magnitudes according to rank
  magsNew <- sort(eventMags)[rankNew]
  df <- eventSubset#[,c("loc", "Northing", "Easting")]
  df$rank <- rankNew
  df$dpoe <- magsNew
  df$low <- rep(low, nrow(df))
  df$U <- U
  df
}

returnLevelsEC <- function(eventSimTable, ncin, paramtable,
                           NH=19914, pthr=(2/360)){
  # Gets empirical quantiles for simulated PoE, and GPA fitted values of flow
  # for simulated PoE below threshold (e.g 2/360 for POT2).
  #
  # eventSimTable     matrix of DPoE, one event per column, one location per row.
  # ncin              netCDF object showing daily flow
  # paramtable        GPA parameters table for POTs
  # NH                number of locations
  # pthr              probability for threshold exceedence, default = POT2
  
  ST0 <- Sys.time()
  ST <- Sys.time()
  eventFlow <- eventSimTable
  for(i in 1:NH){
    
    if((i < 10) | (i %% 200 == 0)){ # time recording
      print(i)
      I <- difftime(Sys.time(), ST, units="secs")
      I0 <- difftime(Sys.time(), ST0, units="secs")
      print(paste("Percent remaining", 100*round((NH-i)/NH ,2)))
      print(paste("Time remaining", round((NH-i)/i * I0,2)))
      print(paste("Since last readout:", round(I,2)))
      ST <- Sys.time()
      print(ST) 
    }
    
    nor <- rn$row[i]
    eas <- rn$col[i]
    tSlice <- ncvar_get(ncin, varid="dmflow",
                        start=c(nor, eas,  1),
                        count=c(1, 1, -1))
    rare <- (eventSimTable[i,] < pthr)
    quant <- quantile(tSlice, probs=1-eventSimTable[i,])
    # fit exceedences using obs if under threshold
    
    # use GPA if over threshold (note divide by pthr to get P[X | X>thr])
    quant[rare] <- qevd(eventSimTable[i,rare]/pthr,
                        thr=paramtable$threshold[i],
                        scale=paramtable$scale[i],
                        shape=paramtable$shape[i], type="GP")
    eventFlow[i,] <- quant
  }
  eventFlow
}






### FITTING TAILS ###------------------------------------------------------

print("Determining tail distribution")

u <- unlist(obs_dpoe, use.names=FALSE)
w <- which(u < 1e-15)
u[w] <- 1e-15
ww <- which((1-u) < 1e-15)
u[ww] <- 1 - 1e-15
m_x <- mean(u, na.rm = TRUE)
s_x <- sd(u, na.rm = TRUE)
alpha <- m_x*((m_x*(1 - m_x)/s_x^2) - 1)
beta <- (1 - m_x)*((m_x*(1 - m_x)/s_x^2) - 1)



FF <- c(alpha,beta)
try({
FF <- fitdist(u, "beta", method="mle",
              start=list(shape1=alpha, shape2=beta))$estimate
#get tails from events
})
rm(u)



### SIMULATION OF NEW EVENTS ###--------------------------------------------
# Number of new events to simulate
print("Simulating new events")

M <- 5000
newEventDraw <- rep(NA, M)
newEventDpe <- matrix(NA, nrow=NH, ncol=M)
for(m in 1:M){
  if(m %% 50 == 0){print(m)}
  gne <-  generateNewEvent(eventSet=obs_dpoe, NE=NE, NH=NH,
                                    maxit=100, low=1, betapar=FF)
  newEventDpe[,m] <-gne$dpoe
  newEventDraw[m] <- gne$U
}


### CONVERSION TO FLOW AND APoE ###------------------------------------------
print("Conversion to flow")
# Conversion to flow
ncin <- nc_open(ncpres)
newEventFlow <- returnLevelsEC(newEventDpe, ncin=ncin, paramtable=partable_pres)

# Conversion to APoE
newEventApe <- 1 - exp(-newEventDpe*360)

fn1 <- function(h){
  wh_ext <- (newEventDpe < (2/360))
  gpa_poe <- (1 - pevd(newEventFlow[h,],
                       threshold=threshMat[h,jT], scale=partable$scale[h],
                       shape=partable$shape[h], type='GP'))

  valsape <- ifelse(wh_ext,
                    1 - exp(-gpa_poe/partable$meanInt[h]), #gpa scaled to year
                    1 - exp(-newEventDpe[h,]*360)) #dpoe scaled to year
  valsape
}

newEventApe <- sapply(1:NH, fn1)

print(Sys.time())

#### SAVE OUTPUTS ####--------------------------------------------------------
print("Saving outputs")

newEventFlow <- cbind(rn, round(as.data.frame(NewEventFlow),8))
W <- paste0("E",1:(ncol(newEventFlow)-4))
W1 <- sapply(1:length(W), function(i){paste0(W[i],".",sum(W[1:i]==W[i]))})
colnames(newEventFlow)[-(1:4)] <- W1

newEventApe <- cbind(rn, round(as.data.frame(NewEventApe),8))
W <- paste0("E",1:(ncol(newEventFlow)-4))
W1 <- sapply(1:length(W), function(i){paste0(W[i],".",sum(W[1:i]==W[i]))})
colnames(newEventFlow)[-(1:4)] <- W1

newEventDpe <- cbind(rn, round(as.data.frame(NewEventDpe),8))
W <- paste0("E",1:(ncol(newEventFlow)-4))
W1 <- sapply(1:length(W), function(i){paste0(W[i],".",sum(W[1:i]==W[i]))})
colnames(newEventDpe)[-(1:4)] <- W1

readr::write_csv(newEventFlow,
                 paste0(data_wd, subfold, "eventflow_EC_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

readr::write_csv(newEventDpe,
                 paste0(data_wd, subfold, "eventdpe_EC_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

readr::write_csv(newEventApe,
                 paste0(data_wd, subfold, "eventape_EC_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

print("107 done.")