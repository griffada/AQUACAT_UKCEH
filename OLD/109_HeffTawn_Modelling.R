#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-17
#
# Using Heffernan and Tawn model for spatial coherence to generate new events.
#
# For aquaCAT, Project 07441.
#
# OUTPUTS: NewEventHT_***.csv,
#          coefficients.rds, 
#          zscores.rds, 
#          depStruct.rds
# 
# Created ABG 2020-06-17
# Pipeline version 2020-09-07
#
#~~~~~~~~~~~~~~~~~~~~~~~
if(interactive()){
  commandArgs <- function(...){c("09","present","NW")}
}
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE", "SEV", "SSC", "SW", "THA", "TRE", "WAL")
if(length(args)==3){
  RCM <- sprintf("%02d",as.numeric(args[1]))
  period <- args[2]
  if(period=="present"){
    suffix <- "_198012_201011"
  }else if (period=="future"){
    suffix <- "_205012_208011"
  }
  REG <- args[3]
  if(!any(REG == regions)){
    stop(paste("incorrect call: Rscript 109_HeffTawn_Modelling.R gcm period region \n",
    "- Region must be one of: ANG, ESC, NE, NSC, NW, SE, SEV, SSC, SW, THA, TRE, WAL."))
  }
}

##### SETUP #####-----------------------------------------------------------
library(texmex)
library(dplyr)
library(extRemes)
library(reshape2)
library(tidyverse)
source(paste0(wd,"07c_texmex_slimline.R"))

nSample <- nSampleHT #<<<<<<<<< see setup_script_00
if(RECOMPUTE_FLAG){print("RECOMPUTING TEXMEX")}
ST <- Sys.time()  
# ws1 <- "pc05"
# thresh1 <- "POT2"
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)
print(paste0("Running for threshold ", thresh1, " at ", ws1, " minimum spread."))

# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS( paste0(data_wd, subfold, "threshDayExcList_RCM",
                                    RCM, suffix,".rds"))


# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"))
#dim(threshMat) #19914 x 5


#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))



# timewise maxima at each cell for each event (NE x NH)

#library(data.table)
#edf <- fread(paste0(wd,"/Data/eventdf_POT2_pc2.csv"), colClasses=rep("numeric",287))

eventFlow <- as.data.frame(readr::read_csv(paste0(data_wd,subfold,  "/", REG,"/eventflow_OBS_region_",
                                REG, "_RCM", RCM, suffix, ".csv"),
                             col_types=cols(.default=col_double()))[,-(1:4)])

eventDpe <- as.data.frame(readr::read_csv(paste0(data_wd,subfold,  "/", REG,"/eventdpe_OBS_region_",
                                   REG, "_RCM", RCM, suffix, ".csv"),
                            col_types=cols(.default=col_double()))[,-(1:4)])

eventApe <- as.data.frame(readr::read_csv(paste0(data_wd,subfold,  "/", REG,"/eventape_OBS_region_",
                                   REG, "_RCM", RCM, suffix, ".csv"),
                            col_types=cols(.default=col_double()))[,-(1:4)])

REG <- REG
r1 <- which(rn_regions$REGION == REG)
NREG <- length(r1)
thresh0 <- unlist(threshMat[r1,jT], use.names=FALSE)

print("SETUP DONE.")
print(Sys.time() - ST)

##### HEFFTAWN CODING #####------------------------------------------------

### KEY ARGUMENTS ---------------------------------
ST <- Sys.time()
DATA <- t(as.data.frame(eventFlow))
mqu <- 0.75
dqu <- 0.8


### STEP 1 ### MARGINALS ###--------------------------------------------------
    marginfns_temp <- list("laplace",
                      p2q = function(p) ifelse(p <  0.5, log(2 * p), -log(2 * (1 - p))),
                      q2p = function(q) ifelse(q <  0, exp(q)/2, 1 - 0.5 * exp(-q)))
### STEP 2 ### TRANSFORM TO LAPLACE ###----------------------------------------
if(!RECOMPUTE_FLAG & file.exists(paste0(data_wd, subfold, "/", REG, "/transformedData.rds"))){
  MODELS <- readRDS(paste0(data_wd, subfold, REG, "/marginal_models.rds"))
  mth_in <- MODELS[[2]]
  mqu_in <- MODELS[[3]]
  MODELS <- MODELS[[1]]
  print("Step 1 completed: Marginal models.")
  TRANSFORMED <- readRDS(paste0(data_wd, subfold, "/", REG, "/transformedData.rds"))
}else{
  marginals <- migpd_slim(mqu=mqu, penalty="none", verbose=F)
    # str(step1_test1, max.level=1)
  mth_in <- marginals$mth
  mqu_in <- marginals$mqu

    # Marginal models using Generalised Pareto distribution.
  MODELS <- marginals$models

  
#}
  print("Step 1 completed: Marginal models.")

    #Transform DATA into Laplace distribution
  TRANSFORMED <- mexTransform_slim(marginfns=marginfns_temp, mth=mth_in,
                                      method="mixture", r=NULL)$TRANSFORMED
  
 
  saveRDS(list(marginals$models, mth_in, mqu_in),
          file=paste0(data_wd, subfold, REG, "/marginal_models.rds"))
  rm(marginals)
  saveRDS(TRANSFORMED, file=paste0(data_wd, subfold, REG, "/transformedData.rds"))
}


print("Step 2 completed: Data transform.")
print(Sys.time() - ST)
ST <- Sys.time()
ST0 <- ST
### STEP 3 ### DEPT STRUCT ### ---------------------
if(!RECOMPUTE_FLAG & file.exists(paste0(data_wd, subfold, REG,  "/depStruct", REG, "_",
                                thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))){
    COEFFS <- readRDS(paste0(data_wd, subfold, REG,  "/coefficients_", REG, "_",
                                thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))
    DEPENDENCE <- readRDS(paste0(data_wd, subfold, REG,  "/depStruct", REG, "_",
                                thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))
    Z <- readRDS(paste0(data_wd, subfold, REG,  "/zScores", REG, "_",
                                thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))
}else{
    COEFFS <- array(NA, dim=c(6,NREG-1,NREG))  # parameters for dep struct
    Z <- vector("list", NREG)  # Z-scores for dependence structure.
    DEPENDENCE <- vector("list", NREG)
    # Compute parametric dependence structure
    for(k in 1:NREG){
      
        if((k < 10) | (k %% 100 == 0)){ # time recording
          print(k)
          I <- difftime(Sys.time(), ST, units="secs")
          I0 <- difftime(Sys.time(), ST0, units="secs")
          print(paste("Percent remaining", 100*round((NREG-k)/NREG ,2)))
          print(paste("Time remaining", round((NREG-k)/k * I0,2)))
          print(paste("Since last readout:", round(I,2)))
          ST <- Sys.time()
          print(ST) 
        }
      
      o <- try(mexDependence_slim(dqu=dqu, mth=mth_in, whch=k,
                                  marginsTransformed=TRANSFORMED))
      
      if(inherits(o, "try-error")){
        if(interactive()) browser()
        warning("Error in mexDependence")
      }
      #str(o)
      if((length(o) > 1) & o$errcode > 0){
        print(paste0("Error code ", o$errcode, " on iter ", k, "."))
      }else{
        DEPENDENCE[[k]] <- o
      }
    }

saveRDS(DEPENDENCE, file=paste0(data_wd, subfold, REG,  "/depStruct", REG, "_",
                                thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))

saveRDS(COEFFS, file=paste0(data_wd, subfold, REG, "/coefficients_", REG, "_",
                            thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))

saveRDS(Z, file=paste0(data_wd, subfold, REG, "/zScores", REG, "_",
                       thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))
}
print("Step 3 Complete: Dependence Structure computation.")
print(Sys.time() - ST)
ST <- Sys.time()

### STEP 4 ### SIMULATION ###-----------------------------------------

nSample <- nSample
mult <- 2
dS <- 50


MCEVENTS <- mexMonteCarlo_slim(mexList=DEPENDENCE,
                                  marginfns=marginfns_temp,
                                  mth=mth_in,
                                  mqu=mqu_in,
                                  nSample=min(nSample,dS), mult=mult)
str(MCEVENTS, max.level=1)

MCS <- t(MCEVENTS$MCsample)

W <- colnames(MCS)
W1 <- sapply(1:length(W), function(i){paste0(W[i],".",sum(W[1:i]==W[i]))})
MCS <- cbind(rn[r1,], MCS)
colnames(MCS)[-(1:4)] <- W1

readr::write_csv(MCS,
                 paste0(data_wd, subfold, REG, "/eventflow_HT_",REG,"_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
saveRDS(MCS,paste0(data_wd, subfold, REG, "/eventflow_HT_",REG,"_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, "_0.RDS"))

if(nSample > dS){
    sampleStep <- seq(0,nSample,by=dS)
    if(max(sampleStep) != nSample){ sampleStep <- c(sampleStep, nSample) }
    for(i in seq_len(length(sampleStep)-1)){
      MCEVENTS <- mexMonteCarlo_slim(mexList=DEPENDENCE,
                                  marginfns=marginfns_temp,
                                  mth=mth_in,
                                  mqu=mqu_in,
                                  nSample=(sampleStep[i+1] - sampleStep[i]),
                                  mult=mult)
      
      MCS <- t(MCEVENTS$MCsample)
      
      W <- colnames(MCS)
      W1 <- sapply(1:length(W), function(i){paste0(W[i], ".", sum(W[1:i]==W[i]))})
      MCS <- cbind(rn[r1,], MCS)
      colnames(MCS)[-(1:4)] <- W1
      
      saveRDS(MCS,paste0(data_wd, subfold, REG, "/eventflow_HT_", REG, "_",
                          thresh1, "_", ws1, "_RCM", RCM, suffix, "_", i, ".RDS"))
  }
}

# MCS <- cbind(rn[r1, ], t(MCS))

### SAVE OUTPUTS ###---------------------------------------
print("Step 4 Complete: New events simulated.")
print(Sys.time() - ST)
ST <- Sys.time()
