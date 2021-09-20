#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2021-05-04
#
# The Empirical Copula functions for simulating new events from given data.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-15
# Pipeline version ABG 2020-09-07
# Version 2 using empirical beta copulas
#
# OUTPUTS: NewEventPresentEC_***.csv: data table of events, one event per row.
#
#~~~~~~~~~~~~~~~~~~~~~~~

if(interactive()){commandArgs <- function(...){c("05","present")}}
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
library(parallel)
library(foreach)
library(ilaprosUtils)
library(lmomco)
})

### DATA --------------------------

# thresh1 <- "POT2" #!#!#!#!# Important constants to select.
# ws1 <- "pc05"
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)

print(paste("Running for threshold", thresh1, "at", ws1, "minimum spread."))
print(paste("RCM", RCM, "period", period))

suffix_pres <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix_pres, "_out.nc") 

readdf <- function(...){as.data.frame(data.table::fread(...))}

# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"))
#dim(threshMat) #19914 x 5

NH <- nrow(rn)

# timewise maxima at each cell for each event ((NE + 4) x NH)
obs_events  <- readdf(paste0(data_wd,subfold,"eventflow_OBS_",thresh1, "_", ws1,
                             "_RCM", RCM, suffix, ".csv"))[,-(1:4)]

obs_dpoe    <- readdf(paste0(data_wd,subfold,"eventdpe_OBS_",thresh1, "_", ws1,
                            "_RCM", RCM, suffix, ".csv"))[,-(1:4)]

obs_apoe    <- readdf(paste0(data_wd,subfold,"eventape_OBS_",thresh1, "_", ws1,
                              "_RCM", RCM, suffix, ".csv"))[,-(1:4)]

partable    <- readdf(paste0(data_wd,subfold,"paramtableG_",thresh1,
                              "_RCM", RCM, suffix, ".csv"))

partable_pres <- readdf(paste0(data_wd,subfold_pres, 
                  "paramtableG_", thresh1, "_RCM", RCM, suffix_pres, ".csv"))

NE <- ncol(obs_events)

suffix_pres <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- ncoriginal <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix_pres, "_out.nc") 

# ncin <- nc_open(ncoriginal) # This file is ~2.5GB on the linux server.
# print(ncin)
# print(floor(Sys.time() - ST))

ncin_pres <- nc_open(ncpres)

### COPULA FUNCTIONS ###---------------------------------------------------

rank_events <- apply(obs_events, 1, rank, ties.method = "random")

flow_glo <- function(v, params, OBS){
 thresh0 <- params$threshold
 threshquan <- params$threshquan
 qout <- rep(0, length(v))
 up <- v > threshquan
 if(any(up)){
 vsub <- v[up]
 qout[up] <- qglo(1 - (1-vsub)/(1-threshquan),
                  loc=params$loc, scale=params$sca, sh=params$shape)
 }
 if(any(!up)){
 vsub <- v[!up]
 qout[!up] <- quantile(OBS, probs=vsub)
 }
 qout[qout < 0] <- 1e-3
 unlist(qout)
}


### SIMULATION OF NEW EVENTS ###--------------------------------------------
# Number of new events to simulate
print("Simulating new events")
#Msims <- 20
dM <- 5 #number of draws per tester; keep small.
M <- Msims - (Msims%%dM)
print(paste("> > > > > > M = ", M))
#newEventDraw <- rep(NA, M)
newEventDpe <- matrix(NA, nrow=NH, ncol=M)

newEventDraw <- rep(NA, M)

areas <- apply(obs_dpoe, 2, function(x){sum(x<(2/360))})

cl <- parallel::makeCluster(3, outfile = "")
doParallel::registerDoParallel(cl)
Vnew <- foreach(m = 1:(M/dM) , .combine = cbind) %dopar% {
    #print(m*dM)
    U <- matrix(runif(NH * NE), nrow = NE, ncol = NH)
    V <- apply(U, 2, sort)  # This is the slow bit.
    w <- sample.int(NE, dM)
    newEventDraw[1:5 + (m-1)*dM] <- w
    Vnew0 <- matrix(NA, nrow=dM, ncol=NH)
    for (i in 1:NH) {
      Vnew0[,i] <- V[rank_events[w, i], i]
    }
    t(Vnew0)
}
stopCluster(cl)

saveRDS(Vnew, paste0(data_wd, subfold, "Vnew_temp.RDS"))
#Vnew <- readRDS(paste0(data_wd, subfold, "Vnew_temp.RDS"))
#print(M)
#print(dim(Vnew))
newEventFlow <- matrix(NA, nrow=NH, ncol=ncol(Vnew))

for(i in 1:NH){
  if(i < 10 | i %% 1000 == 0){print(i)}
  o <- try(
  flow_glo(Vnew[i, ], params = partable_pres[i, ], OBS = unlist(obs_events[i, ]))
  )
  #print(dim(o))
  if(inherits(o, c("Error", "try-error"))){
   print(i)
   print(o)
  }else{
   newEventFlow[i,] <- unlist(o)
  }
}

readr::write_csv(cbind(rn, signif(as.data.frame(newEventFlow),6)),
                 paste0(data_wd, subfold, "eventflow_EC2_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"), append=T)

# newEventFlow <- data.matrix(newEventFlow)

### GETTING PoE FOR NEW EVENTS FROM FLOW ###----------------------------------
# newEventFlow <- readdf(paste0(data_wd, subfold, "eventflow_EC2_",
#                         thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

newEventDpe <- matrix(NA, ncol=M, nrow=NH)
eventDpeFrame <- matrix(NA, ncol=M, nrow=NH)
newEventApe <- matrix(NA, ncol=M, nrow=NH)
eventApeFrame <- matrix(NA, ncol=M, nrow=NH)

glo_tracker <- c()
glo_worst <- c()

ST0 <- Sys.time()
ST <- Sys.time()
print("loop start")
for(h in 1:NH){
  #print(h)
  if((h < 10) | (h %% 1000 == 0)){ # time recording
    print(h)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((NH-h)/NH ,2)))
    print(paste("Time remaining", round((NH-h)/h * I0,2)))
    print(paste("Since last readout:", round(I,2)))
    ST <- Sys.time()
    print(ST) 
  }
  
  thr <- threshMat[h,jT]
  meanInt <- partable_pres$meanint[h]
  thresholdH <- partable_pres$threshold[h]
  locH <- partable_pres$loc[h]
  scaleH <- partable_pres$sca[h]
  shapeH <- partable_pres$shape[h]

  vals <- ncvar_get(ncin_pres, "dmflow",
                         start=c(rn$row[h], rn$col[h], 1),
                         count=c(1, 1, -1))

  ecd <- ecdf(c(-1,vals,1e8))

  valsdpe <- 1 - ecd(newEventFlow[h,])
  
  # APoE per cell for event maxima based on either ecdf or glo
  wh_ext <- (valsdpe < threshVal[jT])

  ub <- ifelse(shapeH < 0, thresholdH - (scaleH/shapeH), -9999)

  glo_poe <- (1 - lmomco::cdfglo(unlist(newEventFlow[h,]),
                       vec2par(c(locH,scaleH,shapeH), type='glo')))
  glo_poe[is.na(glo_poe)] <- 0
  
  glo_tracker[h] <- sum(glo_poe < 1e-3)
  glo_worst[h] <- 1/min(glo_poe)
  # if(any(glo_poe < 1e-3)){
  #   print(paste("****", h))
  #   #print(paste("peaksonly:", locH, scaleH, 1/min(glo_poe)))
  #   print(paste("obs beyond 1000 yr", paste(unlist(newEventFlow[h,])[glo_poe < 1e-3],
  #                                           collapse=" ")))
  # }

  valsape <- 1 - exp(-valsdpe*360)
  valsape[wh_ext] <- 1 - exp(-glo_poe[wh_ext]/meanInt)
  
  eventDpeFrame[h,]  <- valsdpe
  eventApeFrame[h,]  <- valsape
}

print(quantile(glo_tracker, probs=seq(0,1,by=0.1)))
print(quantile(glo_worst, probs=seq(0,1,by=0.1)))

newEventDpe <- cbind(rn, signif(as.data.frame(eventDpeFrame),6))
W <- paste0("E",newEventDraw)
W1 <- sapply(1:length(W), function(i){paste0(W[i],".",sum(W[1:i]==W[i]))})
colnames(newEventDpe)[-(1:4)] <- W1

newEventDpe0 <- cbind(rn, signif(as.data.frame(eventApeFrame),6))
colnames(newEventDpe0)[-(1:4)] <- W1


#### SAVE OUTPUTS ####--------------------------------------------------------
print("Saving outputs")
readr::write_csv(as.data.frame(eventDpeFrame),
                 paste0(data_wd, subfold, "eventdpe_EC2_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"), append=T)
readr::write_csv(as.data.frame(eventApeFrame),
                 paste0(data_wd, subfold, "eventape_EC2_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"), append=T)


file.remove(paste0(data_wd, subfold, "Vnew_temp.RDS"))

print(Sys.time())
print("EC2 done.")