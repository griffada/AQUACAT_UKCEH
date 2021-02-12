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

ws1 <- "pc05"
thresh1 <- "POT2"
jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)
print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

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
    stop(paste("incorrect call: Rscript 110_HT_PoEEstimation.R gcm period region \n",
               "- Region must be one of: ANG, ESC, NE, NSC, NW, SE, SEV, SSC, SW, THA, TRE, WAL."))
  }
}


##### DATA #####-------------------------------------------------------

#5 lists of NH lists, one for each threshold

thresMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                           RCM, suffix, ".rds"))

print(ST <- Sys.time())
ncin <- nc_open(ncoriginal)
print(Sys.time() - ST)

partable <- readr::read_csv(paste0(data_wd,subfold, 
                            "paramtable_",thresh1, "_RCM", RCM, suffix, ".csv"),
                     col_types=cols(
                       .default = col_double()
                     ))

eventflow_HT <- readr::read_csv(paste0(data_wd,subfold, "/", REG, "/eventflow_HT_",
                                REG,"_",thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"),
                           col_types=cols(.default = col_double()))[,-(1:4)]

NE <- ncol(eventflow_HT) - 4

rn_regions$locnum <- 1:nrow(rn_regions)

rn_reg <- data.frame(rn_regions) %>% dplyr::filter(REGION == REG)

r1 <- which(rn_regions$REGION == REG)

partable <- partable[r1,]

### Prealloc ###--------------------------------------------------------------
eventflow_HT <- as.matrix(eventflow_HT)

#rarityDF <- eventflow_HT

eventDpeFrame <- matrix(NA, ncol=ncol(obs_events), nrow=NH)
eventApeFrame <- matrix(NA, ncol=ncol(obs_events), nrow=NH)

##### DPoE Calculation using GPA #####-----------------------------------------
for(h in 1:nrow(eventflow_HT)){
  
  if((h < 10) | (h %% 200 == 0)){ # time recording
    print(n)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((NH-h)/NH ,2)))
    print(paste("Time remaining", round((NH-h)/h * I0,2)))
    print(paste("Since last readout:", round(I,2)))
    ST <- Sys.time()
    print(ST) 
  }
  
  H <- rn_reg$locnum[h]
  
  thr <- thresMat[H,jT]
  
  meanInt <- partable$meanint[h]
  scaleH <- partable$scale[h]
  shapeH <- partable$shape[h]
  
  vals <- ncvar_get(ncin, "dmflow",
                    start=c(rn$row[H], rn$col[H], 1),
                    count=c(1, 1, -1))
  
  # DPoE per cell for event maxima based on whole time series.
  ecd <- ecdf(vals)
  valsdpe <- 1 - ecd(eventflow_HT[h,])
  
  # APoE per cell for event maxima based on either ecdf or gpa
  wh_ext <- (eventflow_HT[h,] > thr)
  gpa_poe <- (1 - pevd(as.numeric(eventflow_HT[h,]),
                       threshold=thr, scale=scaleH, shape=shapeH, type='GP'))
  
  valsdpe[wh_ext] <- (2/360)*gpa_poe
  
  valsape <- ifelse(wh_ext,
                    1 - exp(-gpa_poe/meanInt), #gpa scaled to year
                    1 - exp(-valsdpe/360)) #dpoe scaled to year
  
  eventDpeFrame[h,]  <- valsdpe
  eventApeFrame[h,]  <- valsape
  
}

eventDpeFrame  <- cbind(rn, eventDpeFrame)
eventApeFrame  <- cbind(rn, eventApeFrame)

##### OUTPUTS #####-----------------------------------------------------------

readr::write_csv(as.data.frame(eventDpeFrame), path=paste0(data_wd,subfold, "/", REG,
                  "/eventdpe_HT_",REG,"_", thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
readr::write_csv(as.data.frame(eventApeFrame), path=paste0(data_wd,subfold, "/", REG,
                  "/eventape_HT_",REG, "_", thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

