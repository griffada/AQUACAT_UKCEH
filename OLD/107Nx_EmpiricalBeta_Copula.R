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

if(interactive()){commandArgs <- function(...){c("13","present")}}
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
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

# if(settings$EC2flow){
#   print("EC flow exists already. Stopping 107N.")
# }else{

print(paste("Running for threshold", thresh1, "at", ws1, "minimum spread."))
print(paste("RCM", RCM, "period", period))

suffix_pres <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
#ncpres <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix_pres, "_out.nc") 

readdf <- function(...){as.data.frame(data.table::fread(...))}

# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"))
#dim(threshMat) #19914 x 5

NH <- nrow(rn)

# timewise maxima at each cell for each event ((NE + 4) x NH)
obs_events  <- nc_open(paste0(data_wd,subfold,"eventOBS_",thresh1, "_", ws1,
                             "_RCM", RCM, suffix, ".nc"))

partable    <- readdf(paste0(data_wd,subfold,
                  "paramtableH_",thresh1, "_RCM", RCM, suffix, ".csv"))

partable_pres <- readdf(paste0(data_wd,subfold_pres, 
                  "paramtableH_", thresh1, "_RCM", RCM, suffix_pres, ".csv"))

NE <- obs_events$dim$event$len

savepath <- paste0(data_wd, subfold, "eventECD_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".nc")
cdfPrimer(RCM, period, "EC2", NE=Msims, NH, thresh1, ws1, rn, savepath, chunks=F)
ec_events <- nc_open(savepath, write=T)

### COPULA FUNCTIONS ###---------------------------------------------------

flow_glo <- function(v, params, eventthresh, OBS, pot2){
    flow <- rep(0, length(v))
    up <- v < 1-eventthresh
    setpars <- pargpa(lmoms(OBS[OBS > params$threshold]),
                      xi=params$threshold)

    if(any(up)) flow[up] <- quagpa(1 - (v[up])/(1-eventthresh), setpars)
    if(any(!up)) flow[!up] <- quantile(OBS, 1-v[!up])

    flow[flow < 0] <- 1e-5
    #print(sum(flow > pot2))
    flow
}

betacop <- function(NE, NH, dM, rankevents, thresholds=(2/360), mincov=0.002){
    U <- matrix(runif(NH * NE), nrow = NH, ncol = NE) #NHxNE
    V <- apply(U, 1, sort) # This is the slow bit.
    V <- sapply(1:NH, function(i){V[rankevents[i, ], i]})
    ap <- apply(V, 1, function(v){sum(v < thresholds)})
    w0 <- which(ap >= (mincov*NH))
    if(length(w0) < dM){w0 <- c(rep(w0, ceiling(dM/length(w0))))}
    w <- sample(w0, dM, replace=F)
    Vnew0 <- t(V[w, ])
    return(list(newEventDraw=w, Vnew0=Vnew0))
}

### SIMULATION OF NEW EVENTS ###--------------------------------------------
# Number of new events to simulate
print("Simulating new events")
Msims <- 20
dM <- 5 #number of draws per tester; keep small.
M <- Msims - (Msims%%dM)
print(paste("> > > > > > M = ", M))

newEventDraw <- c()

obs_slice <- ncvar_get(obs_events, "flow") #NH x NE

obs_thresh <- 1 - sapply(1:NH,
                         function(i){ecdf(obs_slice[i,])(threshMat[i,jT])})

rank_events <- t(apply(-obs_slice, 1, rank, ties.method = "random")) 

# cl <- parallel::makeCluster(3, outfile = "")
# doParallel::registerDoParallel(cl)
Vnew1 <- foreach(m = 1:(M/dM)) %do% {
  if(m %% 200 == 0){print(paste0(m,"/",M))}
  B <- betacop(NE, NH, dM, rank_events, thresholds=obs_thresh)
  newEventDraw <<- c(newEventDraw, B$newEventDraw)
  B
}
#stopCluster(cl)
Vnew <- do.call(cbind, lapply(Vnew1, function(v){v$Vnew0}))
newEventDraw <- do.call(c, lapply(Vnew1, function(v){v$newEventDraw}))
# V to flow
for(h in 1:NH){
  if(h < 10 | h %% 1000 == 0){print(h)}
  obs_location <- obs_slice[h,]
  o <- try(
    flow_glo(Vnew[h, ],
             params = partable_pres[h, ],
             eventthresh = ecdf(obs_location)(partable_pres$threshold[h]),
             OBS = obs_location,
             pot2 = threshMat[h,jT])
  )
  if(inherits(o, c("Error", "try-error"))){
   print(h)
   print(o)
  }else{
    o <- unlist(o)
    ncvar_put(ec_events, "flow", unlist(o), start=c(h,1), count=c(1,Msims))
  }
}

W <- newEventDraw
W1 <- sapply(1:length(W), function(i){paste0(W[i],".",sum(W[1:i]==W[i]))})

ncvar_put(ec_events, "eventNo", W, start=1, count=Msims)

nc_close(obs_events)
nc_close(ec_events)

settings$EC2flow <- TRUE
write_yaml(settings, settingspath)
}
print(Sys.time())
print("EC2 done")
