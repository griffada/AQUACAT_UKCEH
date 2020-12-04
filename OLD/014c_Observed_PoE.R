# Get annual return periods, and 


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

ws1 <- "pc05"
thresh1 <- "POT2"
jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)
print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

obs_events <- readr::read_csv(paste0(data_wd,subfold,
                      "eventdf_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

NE <- ncol(obs_events) - 4

thresMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM", RCM, suffix, ".rds"))

load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))

partable <- readr::read_csv(paste0(data_wd,subfold, 
                           "paramtable_",thresh1, "_RCM", RCM, suffix, ".csv"))

colnames(partable)[1] <- "meanint"
# one event per column, from col 5.

rarityDF <- data.frame(eventNo = numeric(),
                         loc = numeric(),
                         Easting = numeric(),
                         Northing = numeric(), 
                         thresh = numeric(),
                         DayS = numeric(),
                         val = numeric(),
                         gpa_apoe = numeric(),
                         rp_years = numeric())

for(h in 1:NH){
  thr <- thresMat[h,jT]
  meanInt <- partable$meanint[h]
  scaleH <- partable$scale[h]
  shapeH <- partable$shape[h]

  obs_events_h <- obs_events[h,-(1:4)]
  
  at_site_poeE <- 1 - pevd(as.numeric(obs_events_h), threshold=thr,
                           scale=scaleH, shape=shapeH, type='GP')
  
  rp_ver <- ifelse(at_site_poeE > 1-(1e-8), NA, meanInt/at_site_poeE)  
            # expected rate per year given POT and GPA PoE.
  
  at_site_apoe <- ifelse(is.na(rp_ver1E), NA, 1 - exp(-at_site_poeE/meanInt))  
                # Poisson assumption
  
  rarityTemp <- data.frame(eventNo = 1:NE,
                           loc = h,
                           Easting = rn[h, 1],
                           Northing = rn[h, 2], 
                           thresh = thr,
                           DayS = eventDayList[[jT]][[jW]],
                           val = obs_events_h,
                           gpa_apoe = at_site_apoe,
                           rp_years = rp_ver)
    
  rarityDF <- rbind(rarityDF, rarityTemp)
}


saveRDS(rarityDF, paste0(wd_id, "observed_rps.rds"))