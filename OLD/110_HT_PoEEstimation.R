#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-06-17
#
# Estimating probability of exceedence for Heffernan-Tawn events.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-06-17
# Pipeline version 2020-09-07
#
# OUTPUTS: HTraritydf_***.csv
#
#~~~~~~~~~~~~~~~~~~~~~~~
if(interactive()){
  commandArgs <- function(...){c("04","present","NW")}
}
#### SETUP ####----------------------------------------------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

library(ilaprosUtils)
library(lmomco)
library(extRemes)
library(pastecs)
library(dplyr)

# ws1 <- "pc05"
# thresh1 <- "POT2"
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)
# print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE", "SEV", "SSC", "SW", "THA", "TRE", "WAL")

if(length(args)==3){
  RCM <- sprintf("%02d", as.numeric(args[1]))
  period <- args[2]
  if(period=="present"){
    suffix <- "_198012_201011"
  }else if (period=="future"){
    suffix <- "_205012_208011"
  }
  REG <- args[3]
  if(!any(REG == regions)){
    stop(paste("incorrect call: Rscript 110_HT_PoEEstimation.R gcm period region \n",
               "- Region must be one of: ANG, ESC, NE, NSC, NW, SE, SEV, SSC, SW, THA, TRE, WAL."))
  }
}

suffix_pres <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix_pres, "_out.nc") 

##### DATA #####-------------------------------------------------------

#5 lists of NH lists, one for each threshold

threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                           RCM, suffix, ".rds"))

ST <- Sys.time()
ncin <- nc_open(ncpres)
print(Sys.time() - ST)

partable <- readr::read_csv(paste0(data_wd,subfold, 
                            "paramtableG_",thresh1, "_RCM", RCM, suffix, ".csv"),
                     col_types=cols(
                       .default = col_double()
                     ))

eventflow_HT <- readr::read_csv(paste0(data_wd,subfold, REG, "/eventflow_HT_",
                                REG,"_",thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"),
                           col_types=cols(.default = col_double()))[,-(1:4)]

NE <- ncol(eventflow_HT)

rn_regions$locnum <- 1:nrow(rn_regions)

rn_reg <- data.frame(rn_regions) %>% dplyr::filter(REGION == REG)

r1 <- which(rn_regions$REGION == REG)

NH1 <- length(r1)

partable <- partable[r1,]

### Prealloc ###--------------------------------------------------------------
eventflow_HT <- as.matrix(eventflow_HT)

#rarityDF <- eventflow_HT

eventDpeFrame <- matrix(NA, ncol=ncol(eventflow_HT), nrow=NH1)
eventApeFrame <- matrix(NA, ncol=ncol(eventflow_HT), nrow=NH1)

##### DPoE Calculation using GPA #####-----------------------------------------
ST <- Sys.time()
ST0 <- Sys.time()
for(h in 1:nrow(eventflow_HT)){
  
  if((h < 10) | (h %% 200 == 0)){ # time recording
    print(h)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((NH-h)/NH ,2)))
    print(paste("Time remaining", round((NH-h)/h * I0,2)))
    print(paste("Since last readout:", round(I,2)))
    ST <- Sys.time()
    print(ST) 
  }
  
  H <- rn_reg$locnum[h]
  
  # thr <- thresMat[H,jT]
  # 
  # meanInt <- partable$meanint[h]
  # scaleH <- partable$scale[h]
  # shapeH <- partable$shape[h]
  
  thr <- threshMat[H,jT]
  meanInt <- partable$meanint[h]
  thresholdH <- partable$threshold[h]
  scaleH <- partable$sca[h]
  shapeH <- partable$shape[h]

  
  vals <- ncvar_get(ncin, "dmflow",
                    start=c(rn$row[H], rn$col[H], 1),
                    count=c(1, 1, -1))
  
  # DPoE per cell for event maxima based on whole time series.
  ecd <- ecdf(c(-1,vals,1e8))
  valsdpe <- 1 - ecd(eventflow_HT[h,])
  
  thr <- threshMat[H,jT]
  meanInt <- partable$meanint[h]
  thresholdH <- partable$threshold[h]
  locH <- partable$loc[h]
  scaleH <- partable$sca[h]
  shapeH <- partable$shape[h]

  # APoE per cell for event maxima based on either ecdf or gpa
  wh_ext <- (valsdpe < threshVal[jT])
  wh_ext[is.na(wh_ext)] <- FALSE
  
  ub <- ifelse(shapeH < 0, thresholdH - (scaleH/shapeH), -9999)

  gpa_poe <- (1 - cdfglo(as.numeric(unlist(eventflow_HT[h,])),
                       vec2par(c(locH,scaleH,shapeH), type='glo')))
  gpa_poe[is.na(gpa_poe)] <- 1
  
  #gpa_tracker[h] <- sum(gpa_poe < 1e-3)
  #gpa_worst[h] <- 1/min(gpa_poe)
  if(any(gpa_poe < 1e-3)){
    print(paste("****", h))
    #print(paste("peaksonly:", locH, scaleH, 1/min(gpa_poe)))
    print(paste("obs beyond 1000 yr", paste(eventflow_HT[h, gpa_poe < 1e-3],
                                            collapse=" ")))
  #gpa_poe[gpa_poe < 1e-5] <- 1e-5
  #print(paste("maxGPA = ", ub))
  }
  #valsdpe <- gpa_poe
  valsape <- 1 - exp(-valsdpe*360)
  valsape[wh_ext] <- 1 - exp(-gpa_poe[wh_ext]/meanInt)
  
  # valsape <- ifelse(wh_ext,
  #                   1 - exp(-gpa_poe/meanInt), #gpa scaled to year
  #                   1 - exp(-valsdpe/360)) #dpoe scaled to year
  
  eventDpeFrame[h,]  <- valsdpe
  eventApeFrame[h,]  <- valsape
  
}

eventDpeFrame  <- cbind(rn[r1,], eventDpeFrame)
eventApeFrame  <- cbind(rn[r1,], eventApeFrame)

##### OUTPUTS #####-----------------------------------------------------------

readr::write_csv(as.data.frame(eventDpeFrame), path=paste0(data_wd,subfold, REG,
                  "/eventdpe_HT_",REG,"_", thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
readr::write_csv(as.data.frame(eventApeFrame), path=paste0(data_wd,subfold, REG,
                  "/eventape_HT_",REG, "_", thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

