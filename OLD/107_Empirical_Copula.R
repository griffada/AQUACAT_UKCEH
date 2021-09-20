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

if(interactive()){commandArgs <- function(...){c("06","future")}}
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

##### SETUP #####---------------------------------------------------
suppressPackageStartupMessages({
library(extRemes)
library(dplyr)
library(fitdistrplus)
})

# if(file.exists(paste0(data_wd, subfold, "eventdpe_EC_",
#                         thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))){
#   stop("EC events already exist. Finishing 107.")}


### DATA --------------------------

# thresh1 <- "POT2" #!#!#!#!# Important constants to select.
# ws1 <- "pc05"
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)

print(paste("Running for threshold", thresh1, "at", ws1, "minimum spread."))

suffix_pres <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix_pres, "_out.nc") 



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
obs_events  <- as.data.frame(readr::read_csv(paste0(data_wd,subfold, "eventflow_OBS_",thresh1,
                                      "_", ws1, "_RCM", RCM, suffix, ".csv"),
                               col_types=cols(.default = col_double()))[,-(1:4)])

obs_dpoe  <-  as.data.frame(readr::read_csv(paste0(data_wd,subfold, "eventdpe_OBS_",thresh1,
                                    "_", ws1, "_RCM", RCM, suffix, ".csv"),
                             col_types=cols(.default = col_double()))[,-(1:4)])

obs_apoe  <-  as.data.frame(readr::read_csv(paste0(data_wd,subfold, "eventape_OBS_",thresh1,
                                    "_", ws1, "_RCM", RCM, suffix, ".csv"),
                             col_types=cols(.default = col_double()))[,-(1:4)])

partable <-  as.data.frame(readr::read_csv(paste0(data_wd,subfold, 
                           "paramtableG_",thresh1, "_RCM", RCM, suffix, ".csv"),
                            col_types=cols(.default= col_double())))

partable_pres <-  as.data.frame(readr::read_csv(paste0(data_wd,subfold_pres, 
                        "paramtableG_",thresh1, "_RCM", RCM, suffix_pres, ".csv"),
                            col_types=cols(.default= col_double())))

# colnames(partable)[1] <- "meanInt"
# colnames(partable_pres)[1] <- "meanInt"





### COPULA FUNCTIONS ###---------------------------------------------------

tailsGen <- function(n, eventSet, threshVal, nabove=0, low=1, maxit=1000){
  # draw new probabilities of exceedence from pool, allowing minumum lower bound
  # to be enforced. Technically drawing from the conditional distribution 
  # P[X=x | X_i < low for some i].
  #
  # n         number of PoEs to draw
  # low       lower target on PoE
  # maxit     number of trials to get an event of the required probability
  # betapar   vector of two parameters for the Beta distribution
  
  # QV <- rep(2, n)
  # #TODO needs the right distribution to extrapolate from pool to reach tails.
  # z <- 0
  # while(min(QV) > low){
  #   V <- runif(n)
  #   z <- z+1
  #   QV <- qbeta(V, betapar[1], betapar[2])
  #   QV <- QV*(1-2e-5) + 1e-5
  #   if(z >= maxit){
  #     print("maxit exceeded")
  #     break
  #   }
  # }
  # if(z>1){print(paste(z, "iterations."))}
  # return(unname(QV))
  N <- 1
  Y0 <- rep(2, nabove)
  Y1 <- rep(2, n-nabove)
  UL <- log(unlist(eventSet, use.names=F))
  logt <- log(threshVal)
  
  while(any(Y0 > threshVal | Y0 < -14) & N < maxit){
    UL0 <- UL[UL < (logt+0)]
    NY <- sum(Y0 > threshVal | Y0 < -14)
    #print(NY)
    Z <- sample(UL0, NY, replace=T) + rnorm(NY,0,0.06)
    Y0[Y0 > threshVal | Y0 < -14] <- exp(Z)
    N <- N+1
  }
  while(any((Y1 < threshVal) | (Y1 > 1)) & N < maxit){
    UL0 <- UL[UL > (logt-0)]
    NY <- sum(Y1 < threshVal | Y1 > 1)
    #print(NY)
    Z <- sample(UL0, NY, replace=T) + rnorm(NY,0,0.0008)
    Y1[Y1 < threshVal | Y1 > 1] <- exp(Z)
    N <- N+1
  }
  if(N == maxit){
    print("Hit maxit.")
  }
  
  return(pmax(0,pmin(1,c(Y0,Y1))))
  
  
}

generateNewEvent <- function(eventSet=obs_dpoe, NE=285, NH=19914,
                             maxit=10000, low=1, threshVal=(2/360)){
  # generates a new widespread event based on a given event library
  #
  # eventSet  grid of daily PoE, one location per row, one event per col.
  # NE        number of events in the set
  # NH        number of gridpoints involved (~20000 in standard set)
  
  U <- sample.int(NE, 1)
  V <- sample.int(NE, 5, replace=T)
  #print(U)
  #print(V)
  eventSubset <-  eventSet[,U]
  #print(range(eventSet[,V]))
  su <- sum(eventSubset < threshVal)
  su <- max(0, min(NH, floor(su*rnorm(1, 1, 0.05)))) # number of rare points.
  rankNew <- rank(eventSubset, ties.method="random")
  #eventMags <- exp(pmin(1e-4,sample(log(unlist(eventSet)), NH, replace=T)+rnorm(NH,0,0.6)))
  eventMags <- tailsGen(n=NH, nabove=su, threshVal=threshVal, eventSet=eventSet[,V])
  # add the magnitudes according to rank
  magsNew <- sort(eventMags)[rankNew]
  df <- as.data.frame(eventSubset)#[,c("loc", "Northing", "Easting")]
  
  df$rank <- unname(rankNew)
  df$dpoe <- magsNew
  df$low <- rep(low, nrow(df))
  df$U <- U
  df$su <- su
  df
}


### FITTING TAILS ###------------------------------------------------------

print("Determining tail distribution")

# FFX <- matrix(NA, nrow=ncol(obs_dpoe), ncol=2)
# for(i in 1:ncol(obs_dpoe)){
# #u <- unlist(obs_dpoe, use.names=FALSE)
# u <- unlist(obs_dpoe[,i])
# #print(range(u))
# w <- which(u < 1e-15)
# u[w] <- 1e-15
# ww <- which((1-u) < 1e-15)
# u[ww] <- 1 - 1e-15
# m_x <- mean(u, na.rm = TRUE)
# s_x <- sd(u, na.rm = TRUE)
# alpha <- m_x*((m_x*(1 - m_x)/s_x^2) - 1)
# beta <- (1 - m_x)*((m_x*(1 - m_x)/s_x^2) - 1)
# 
# DD <- density(log(obs_dpoe))
# SC <- smooth.spline(DD$x,DD$y, nknots=8)
# 
# FF <- c(alpha,beta)
# try({
# FF <- fitdist(u, "beta", method="mle",
#               start=list(shape1=alpha, shape2=beta))$estimate
# #get tails from events
# })
# rm(u)
# FFX[i,] <- FF
# 
# }
# 
# CL <- cov(log(FFX))
# mu <- c(mean(log(FFX[,1])), mean(log(FFX[,2])))

### SIMULATION OF NEW EVENTS ###--------------------------------------------
# Number of new events to simulate
print("Simulating new events")

M <- Msims
print(paste("> > > > > > M = ", M))
newEventDraw <- rep(NA, M)
newEventDpe <- matrix(NA, nrow=NH, ncol=M)
newEventArea <- c()

areas <- apply(obs_dpoe, 2, function(x){sum(x<(2/360))})

for(m in 1:M){
  if ( m < 11 |(m %% 50 == 0)) {print(m)}
  # FFa <- c(-1,-1)
  # while(any(FFa < 0)){
  # FFa <- exp(mvrnorm(1, mu=mu, Sigma=CL))
  # }
  gne <-  generateNewEvent(eventSet=obs_apoe,  NE=NE, NH=NH,
                                    maxit=100, low=1, threshVal=(1-exp(-2))
                           #, betapar=FFa
                           )
  if(m < 3)str(gne)
  newEventDpe[,m] <-unlist(gne$dpoe)
  #print(range(gne$dpoe))
  newEventArea[m] <- gne$su[1]
  newEventDraw[m] <- gne$U[1]
}
# 1/range(apply(obs_apoe,2,min))
# 
# 1/range(apply(newEventDpe, 2, min))
# 
# apply(newEventDpe, 2, function(x){sum(x<(1-exp(-2)))})
# summary(warnings())

newEventDpe0 <- -log(pmax(1-newEventDpe))/360 #The new dpe grid Annual > Daily
newEventDpe0[is.na(newEventDpe0)] <- 1

newEventDpe <- cbind(rn, round(as.data.frame(newEventDpe),8))
W <- paste0("E",newEventDraw)
W1 <- sapply(1:length(W), function(i){paste0(W[i],".",sum(W[1:i]==W[i]))})
colnames(newEventDpe)[-(1:4)] <- W1

newEventDpe0 <- cbind(rn, round(as.data.frame(newEventDpe0),8))
W <- paste0("E",newEventDraw)
W1 <- sapply(1:length(W), function(i){paste0(W[i],".",sum(W[1:i]==W[i]))})
colnames(newEventDpe0)[-(1:4)] <- W1


#### SAVE OUTPUTS ####--------------------------------------------------------
print("Saving outputs")
readr::write_csv(newEventDpe,
                 paste0(data_wd, subfold, "eventape_EC_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
readr::write_csv(newEventDpe0,
                 paste0(data_wd, subfold, "eventdpe_EC_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))


print(Sys.time())
print("107 done.")