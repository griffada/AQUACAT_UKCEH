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
    stop(paste("incorrect call: Rscript 109b_HT_PoEEstimation.R gcm period region \n",
               "- Region must be one of: ANG, ESC, NE, NSC, NW, SE, SEV, SSC, SW, THA, TRE, WAL."))
  }
}


##### DATA #####-------------------------------------------------------

#5 lists of NH lists, one for each threshold

thresMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                           RCM, suffix, ".rds"))

print(ST <- Sys.time())
ncin <- nc_open(ncname)
print(ncin)
print(Sys.time() - ST)

partable <- read_csv(paste0(data_wd,subfold, 
                            "paramtable_",thresh1, "_RCM", RCM, suffix, ".csv"))


eventDF <- readr::read_csv(paste0(data_wd,subfold, "NewEventHT_",thresh1,"_", ws1,
                                  "_RCM", RCM, suffix, ".csv"),
                           col_types=cols(
                             .default = col_double()
                           ))



### Prealloc ###--------------------------------------------------------------
eventDF <- matrix(unlist(eventDF), nrow=nrow(eventDF), ncol=ncol(eventDF))

rarityDF <- eventDF

NE <- ncol(eventDF)-4

rn_regions$locnum <- 1:nrow(rn_regions)

rn_reg <- data.frame(rn_regions) %>% dplyr::filter(REGION == REG)


##### PoE Calculation using GPA #####-----------------------------------------
for(h in 1:nrow(eventDF)){
  if((h %% 200) == 0){
    print(paste(h, "of", nrow(eventDF)))
  }
  
  thr <- thresMat[rn_reg$locnum[h],jT]
  
  rps <- as.numeric(pevd(eventDF[h,-1], scale=partable$scale[h],
                         shape=partable$shape[h],
                         threshold = thr, type='GP'))
  
  rps[rps < 1e-10] <- NA
  
  rarityDF[h,-1] <- 1 - rps
  
}


##### OUTPUTS #####-----------------------------------------------------------

rarityDF <- data.frame(rarityDF)
colnames(rarityDF) <- 
  c("row", "col","east","nor", paste0("E", 1:(ncol(eventDF)-4)))

write_csv(data.frame(rarityDF), 
          paste0(data_wd,subfold, "HTraritydf_", thresh1, "_", ws1, "_RCM", RCM, suffix, ".csv"))
