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

ec_events <- readr::read_csv(paste0(data_wd,subfold,
                      "NewEventEC_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

NE <- ncol(ec_events) - 4

thresMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                           RCM, suffix, ".rds"))

partable <- readRDS(paste0(wd_id, "parTable1.RDa"))
colnames(partable)[1] <- "meanint"
# one event per column, from col 5.

rarityDF <- data.frame(eventNo = numeric(),
                         loc = numeric(),
                         Easting = numeric(),
                         Northing = numeric(), 
                         thresh = numeric(),
                         val = numeric(),
                         gpa_apoe = numeric(),
                         rp_years = numeric())

for(h in 1:NH){
  
  thr <- thresMat[h,jT]
  meanInt <- partable$meanint[h]
  scaleH <- partable$scale[h]
  shapeH <- partable$shape[h]

  ec_events_h <- ifelse(ec_events[h, ] <= 2/360, ec_events[h,]*360/2, NA)

  at_site_levelE <- qevd(as.numeric(ec_events_h), threshold=thr,
                           scale=scaleH, shape=shapeH, type='GP', lower.tail=F)
  
  rp_ver <- meanInt/ec_events_h
            # expected rate per year given POT and GPA PoE.
  
  at_site_apoe <- ifelse(is.na(rp_ver), NA, 1 - exp(-ec_events_h/meanInt))  
                # Poisson assumption
  
  rarityTemp <- data.frame(eventNo = 1:NE,
                           loc = h,
                           Easting = rn[h, 1],
                           Northing = rn[h, 2], 
                           thresh = thr,
                           val = at_site_levelE,
                           gpa_apoe = at_site_apoe,
                           rp_years = rp_ver,
                           ec_poe = ec_events_h)
    
  rarityDF <- rbind(rarityDF, rarityTemp)

}

readr::write_csv(x=rarityDF,
                 path=paste0(data_wd, subfold, "returnlevelsEC_",
                             thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
nc_close(ncin)
print(Sys.time())