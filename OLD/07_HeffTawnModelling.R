#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-17
#
# Using HT model for spatial coherence.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-17
#
#~~~~~~~~~~~~~~~~~~~~~~~

##### SETUP #####-----------------------------------------------------------
if(T){
  #setwd("S:/")
  library(texmex)
  library(readr)
  library(dplyr)
  library(ncdf4)
  library(extRemes)
  library(reshape2)
  library(tidyverse)
  
  if (substr(osVersion,1,3) == "Win") {
    ncname <- "S:/CodeABG/InterimData/dmflow_timechunks.nc"  # rechunked for spaceslices
    ncname2 <- "S:/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc" 
    # pre-rechunk
    wd <- "S:/"
    wd_id <- "S:/CodeABG/InterimData/"
    wd_cd <- "S:/CodeABG/"
    
  } else {
    ncname <- "/prj/aquacat/CodeABG/InterimData/dmflow_timechunks.nc"
    ncname2 <- "S:/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc"
    wd <- "/prj/aquacat/"
    wd_id <- "/prj/aquacat/CodeABG/InterimData/"
    wd_cd <- "/prj/aquacat/CodeABG/"
  }
  
  source(paste0(wd_cd,"07b_HTfunctions.R"))
  
##### DATA #####-----------------------------------------------------------
  
  ND <- 10800 # Number of days
  
  # river network
  rn <- read_csv(paste0(wd_id, "hasData2.csv"))
  
  rn_regions <- read_csv(paste0(wd_id, "hasData_Regions.csv"))
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
  thresh0 <- unname(unlist(threshMat['X2']))
  #dim(threshMat)  =  19914 x 5
  
  
  # eventLList length of event L, NT lists (by threshold) of NW lists 
  # (by inun cutoff))
  # eventDayList start of event L, NT lists of NW lists
  load(paste0(wd_id, "eventLists03.RDa")) 
  NE <- length(eventDayList[[2]][[4]]) # POT2, 2% inun.
  
  
  # timewise maxima at each cell for each event ((NE + 2) x NH)
  eventDF <- readr::read_csv(paste0(wd,"Data/eventdf_POT2_pc05.csv"))
  
  
  # PoE under different computations with extra data. Tidy format.
}  
print("SETUP DONE.")
##### HEFFTAWN CODING #####------------------------------------------------

#### TRIAL FIT ####-------------------------------------------
if(FALSE){
  dummy_melt <- dcast(dummy, eventNo~loc, value.var="val")
  
  rownames(dummy_melt) <- paste0("E",dummy_melt$eventNo)
  dummy_melt <- dummy_melt[,-1]
  colnames(dummy_melt) <- paste0("L", 1:ncol(dummy_melt))
  
  ST <- Sys.time()
  # fits everything else conditional on location 1.
  mmex <- mex(dummy_melt, mqu=.9, penalty="none", dqu=.7, which=1)
  print(Sys.time() - ST)
  
  ST <- Sys.time()
  myboot <- bootmex(mmex, R=100)  # R=10 is quite low.
  print(Sys.time() - ST)
  plottexmex(myboot, plots="dependence")
  
  ST <- Sys.time()
  mypred <- predict.mex(myboot,  pqu=.99)
  print(Sys.time() - ST)
  plotpredmex(mypred)
  
  
  
  ST <- Sys.time()
  rp <-  rep(0.7,9)
  myAll <- vector("list",9)
  for(i in 1:9){
    myAll[[i]] <- mex(dummy_melt, mth=thresh0[1:9], penalty="none",
                      dqu=0.7, which=i)
  }
  names(myAll) <- names(dummy_melt)
  oldClass(myAll) <- "mexList"
  myMCMC <- mexMonteCarlo(500, myAll)
  print(Sys.time() - ST)
  
  
  ab_mmex <- mmex$dependence$coefficients
  marg_mmex <- sapply(mmex$margins$models, function(L){L$coefficients})
  
  
  myMCMC_dummy <- mexMonteCarlo2(500, myAll)
  print(myMCMC_dummy$nR)
  
  POE_dummy <- 1 - myMCMC_dummy$MCsample_P[myMCMC_dummy$whichMaxAboveThresh,] 
  #PoExceedence
}
### PRESENT FIT ###-----------------------------------------

# Get region.
REG <- "NW"
r1 <- which(rn_regions$REGION == REG) # length = 1437
NREG <- length(r1)
event_region <- eventDF[r1,]
thresh_region <- threshMat[r1,]

# Get events from region.
present_region <- readr::read_csv(paste0(wd_id,"regionalEvents_exampleNW_2.csv"))

thresh0 <- unname(unlist(thresh_region['X2']))

#melt to reshape dataframe
present_melt <- dcast(present_region, eventNo~loc, value.var="val")
rownames(present_melt) <- paste0("E",present_melt$eventNo)
colnames(present_melt) <- paste0("L", 0:NREG)
present_melt <- present_melt[,-1]

D <- 60

present_migpd_list <- vector("list", ceiling(NREG/D))

for(j in 1:ceiling(NREG/D)){
  
  suppressWarnings({
    present_migpd_list[[j]] <- migpd(present_melt[,((j-1)*D+1):(min(NREG, j*D))],
                               mqu=0.7, penalty="none")
  })
                               
}
print("DONE MIGPD STEP")

present_migpd_comb <- present_migpd_list[[1]]
present_migpd_comb$models <- do.call(c, lapply(present_migpd_list, function(x)x$models))
present_migpd_comb$mth <- do.call(c, lapply(present_migpd_list, function(x)x$mth))
present_migpd_comb$data <- present_melt

save(present_migpd_comb, file=paste0(wd_id, "present_migpd_",REG,"_",threshName[2],"_",wsName[4],".RDa"))

print("SAVED INTERIM STEP")
present_mexdep <- vector("list", NREG)

#load(paste0(wd_id, "present_migpd_",REG,"_",threshName[2],"_",wsName[4],".RDa"))

m <- 1
for(k in 1:NREG){
  if(k %% 50 == 0){print(k)}
  suppressWarnings({
    present_mexdep_temp <- mexDependence(present_migpd_comb, dqu=0.7, which=k)
  })
  present_mexdep[[m]] <- present_mexdep_temp
  m <- m+1
  if(m > 99){
    names(present_mexdep) <- names(present_melt[(k-99):k])
    save(present_mexdep, file=paste0(wd_id, "/present_mexdep/present_mexdep_list",REG,"_",threshName[2],"_",wsName[4],"_K",k,".RDa"))
    m <- 1
  }
}
#names(present_mexdep) <- names(present_melt)
oldClass(present_mexdep) <- "mexList"
print("SAVED_PARTIAL_MEXDEPS")
#save(present_mexdep, file=paste0(wd_id, "present_mexdep"))

#present_MCMC <- mexMonteCarloBIG(500, present_mexdep)

#save(present_MCMC, file=paste0(wd,"./Data/HT_MCMC_",REG,"_",threshName[2],"_",wsName[4],".RDa"))
