#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-17
#
# Using Heffernan and Tawn model for spatial coherence to generate new events.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-17
# Pipeline version 2020-09-07
#
#~~~~~~~~~~~~~~~~~~~~~~~

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE", "SEV", "SSC", "SW", "THA", "TRE", "WAL")

if(length(args)==3){
  RCM <- sprintf("%02d", args[1])
  period <- args[2]
  if(period=="present"){
    suffix <- "_198012_201011"
  }else if (period=="future"){
    suffix <- "_205012_201011"
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
  
ws1 <- "pc05"
thresh1 <- "POT2"
jV <- which(threshName==thresh1)
jI <- which(wsName == ws1)
print("Running for threshold", POT2, "at ", ws1, "minimum spread.")

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

eventDF <- readr::read_csv(paste0(data_wd,subfold, "eventdf_region_",REG,"_",thresh1,"_", ws1,
                                  "_RCM", RCM, suffix, ".csv"))
  
  
# PoE under different computations with extra data. Tidy format.
poeTable <- readr::read_csv(paste0(data_wd, subfold, "regionalEvents_",REG,"_RCM", 
                                  RCM, suffix,".csv"))
 
REG <- REG
r1 <- which(rn_regions$REGION == REG)
NREG <- length(r1)
thresh0 <- unlist(threshMat[r1,jV], use.names=FALSE)

print("SETUP DONE.")
print(Sys.time() - ST)
ST <- Sys.time()
##### HEFFTAWN CODING #####------------------------------------------------


poeTable_melt <- dcast(poeTable, eventNo~loc, value.var="val")
rownames(poeTable_melt) <- paste0("E",poeTable$eventNo)
colnames(poeTable_melt) <- paste0("L", 0:NREG)

### KEY ARGUMENTS ---------------------------------
DATA <- poeTable_melt[,-1]
mqu <- 0.7
dqu <- 0.7
nSample <- 101
mult <- 10



marginals <- migpd_slim(mqu=mqu, penalty="none")
# str(step1_test1, max.level=1)
mth_in <- marginals$mth
mqu_in <- marginals$mqu

# Marginal models using Generalised Pareto distribution.
MODELS <- marginals$models
saveRDS(MODELS, file=paste0(data_wd, subfold, "marginal_models.rds"))

rm(marginals)
print("Step 1 completed: Marginal models.")






marginfns_temp <- list("laplace",
                     p2q = function(p) ifelse(p <  0.5, log(2 * p), -log(2 * (1 - p))),
                     q2p = function(q) ifelse(q <  0, exp(q)/2, 1 - 0.5 * exp(-q)))

#Transform DATA into Laplace distribution
TRANSFORMED <- mexTransform_slim(marginfns=marginfns_temp, mth=mth_in,
                                  method="mixture", r=NULL)$TRANSFORMED
saveRDS(TRANSFORMED, file=paste0(data_wd, subfold, "transformedData.rds"))
print("Step 2 completed: Data transform.")
print(Sys.time() - ST)
ST <- Sys.time()






COEFFS <- array(NA, dim=c(6,NREG-1,NREG))  # parameters for dep struct
Z <- array(NA, dim=c(24,NREG-1,NREG))  # Z-scores for dependence structure.
DEPENDENCE <- vector("list", NREG)
# Compute parametric dependence structure
for(k in 1:NREG){
  if (k%%50==0){
    print(paste0("Dependence at", round(100*k/NREG,2), "%"))
  } 
  o <- try(mexDependence_slim(dqu=dqu, mth=mth_in, which=k,
                              marginsTransformed=TRANSFORMED))
  
  if(inherits(o, "try-error")){
    if(interactive()) browser()
    warning("Error in mexDependence")
  }
  #str(o)
  if(o$errcode > 0){
    print("Error code ", o$errcode, " on iter ", k, ".")
  }else{
    DEPENDENCE[[k]] <- o
  }
}

saveRDS(DEPENDENCE, file=paste0(data_wd, subfold, "depStruct.rds"))
saveRDS(COEFFS, file=paste0(data_wd, subfold, "coefficients.rds"))
saveRDS(Z, file=paste0(data_wd, subfold, "zScores.rds"))

print("Step 3 Complete: Dependence Structure computation.")
print(Sys.time() - ST)
ST <- Sys.time()


MCEVENTS <- mexMonteCarlo_slim(mexList=DEPENDENCE,
                                  marginfns=marginfns_temp,
                                  mth=mth_in,
                                  mqu=mqu_in,
                                  nSample=nSample, mult=mult)

str(MCEVENTS, max.level=1)

MCS <- MCEVENTS$MCsample
colnames(MCS) <- paste0("L", 1:NREG)

readr::write_csv(MCEVENTS$MCsample,
                 paste0(data_wd, subfold, "NewEventHT_",REG,"_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"), 
                 append=T)
print("Step 4 Complete: New events simulated.")
print(Sys.time() - ST)
ST <- Sys.time()


