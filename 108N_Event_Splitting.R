#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-07-31
#
# Splitting events by region defined by Hydrological Area. 
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-07-31
# Pipeline version ABG 2020-09-07
# netcdf version ABG 2021-07-07
#
# OUTPUTS: regionalEvents_***.csv,
#          eventdf_region_***.csv
#
#~~~~~~~~~~~~~~~~~~~~~~~


if(interactive()){commandArgs <- function(...){c("04","present")}}

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

#regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE", "SEV", "SSC", "SW", "THA", "TRE", "WAL")
regions <- c("NW")

if(settings$HTsplit){
 stop("Split events already exist. Stopping 108.") 
}

print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

suppressPackageStartupMessages({
library(readr)
library(dplyr)
library(ncdf4)
library(extRemes)
library(reshape2)
library(tidyverse)
})

### DATA ###--------------------------------------------------------------

#Hydrometric areas map
# HA <- readOGR(dsn=paste0(data_wd,"hydrometricAreas"), layer="hyd_areas",
#               stringsAsFactors=FALSE)
# HA@data$HA_NUM <- as.numeric(HA@data$HA_NUM)

param_table <- readdf(paste0(data_wd, subfold, "paramtableG_POT2_RCM",
                             RCM, suffix, ".csv"))

# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS( paste0(data_wd, subfold, "threshDayExcList_RCM",
                                    RCM, suffix,".rds"))


# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"))
thresh0 <- unlist(threshMat[,which(threshName==thresh1)], use.names=FALSE)
#dim(threshMat) #19914 x 5


#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
# load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))
# NE <- length(eventDayList[[jT]][[jW]]) # POT2, 2% inun.

# timewise maxima at each cell for each event ((NE + 2) x NH)

obs_national <- nc_open(paste0(data_wd, subfold, "eventOBS_", thresh1,
                                     "_", ws1, "_RCM", RCM, suffix, ".nc"))



### OBS EVENT SPLITTING ###---------------------------------------------------

#for(REG in c("NW")){
  REG="NW"  
  if(!dir.exists(paste0(data_wd, subfold, "/", REG))){
    dir.create(paste0(data_wd, subfold, "/", REG))
  }
  
  print(paste("splitting", REG))
  
  r1 <- which(rn_regions$REGION == REG) # length = 1437
  thresh_region <- threshMat[r1, jT]
  
  savepath <- paste0(data_wd,subfold, REG, "/",
                  "eventOBS_region_", REG, "_RCM", RCM, suffix, ".nc")
  
  rn0 <- rn
  rn <- rn[r1,]
  param_table <- param_table[r1,]
  cdfPrimer(RCM=RCM, period=period, method="OBS", NE=2, NH=length(r1),
            thresh1=thresh1, ws1=ws1, rn=rn, savepath=savepath, chunks=F)
  obs_region <- nc_open(savepath, write=T)
  
  rn <- rn0
  
  NEreg <- 0
  thresh_region <- threshMat[r1, jT]
  NE <- obs_national$dim$event$len
  WE <- c()
  obs_event_number <- ncvar_get(obs_national, "eventNo")
  ST <- Sys.time()
  ST0 <- Sys.time()
  for(i in 1:NE){
    #print(i)
    
    # if((i < 10) | (i %% 200 == 0)){ # time recording
    #   print(i)
    #   I <- difftime(Sys.time(), ST, units="secs")
    #   I0 <- difftime(Sys.time(), ST0, units="secs")
    #   print(paste("Percent remaining", 100*round((NE-i)/NE ,2)))
    #   print(paste("Time remaining", round((NE-i)/i * I0,2)))
    #   print(paste("Since last readout:", round(I,2)))
    #   ST <- Sys.time()
    #   print(ST) 
    # }
    # obs_dpe_region <- ncvar_get(obs_national, "dpe",
    #                         start=c(1,i), count=c(-1,1))[r1]
    # if(sum(is.na(obs_dpe_region) > 0)){
    #   print(paste0("*", sum(is.na(obs_dpe_region))))
    # }
    
    obs_flow_region <- ncvar_get(obs_national, "flow",
                            start=c(1,i), count=c(-1,1))[r1]
    if(any(obs_flow_region > param_table$threshold)){
      #print(i)
      #next
      NEreg <- NEreg + 1
      obs_events_region <- ncvar_get(obs_national, "flow",
                               start=c(1,i), count=c(-1,1))[r1]
      ncvar_put(obs_region, "flow", obs_events_region,
                start=c(1, NEreg), count=c(-1,1))
      
      obs_dpe_region <- ncvar_get(obs_national, "dpe",
                        start=c(1,i), count=c(-1,1))[r1]
      ncvar_put(obs_region, "dpe", obs_dpe_region,
                start=c(1, NEreg), count=c(-1,1))
      
      obs_ape_region <- ncvar_get(obs_national, "ape",
                                  start=c(1,i), count=c(-1,1))[r1]
      ncvar_put(obs_region, "ape", obs_ape_region,
                start=c(1, NEreg), count=c(-1,1))
      
      WE <- c(WE, obs_event_number[i])
    }
    #if(NEreg > 3){break}
  }
  ncvar_put(obs_region, "event", vals=(1:NEreg), start=1, count=NEreg)
  ncvar_put(obs_region, "eventNo", vals=WE, start=1, count=length(WE))
  print(paste(NEreg, "OBS events saved for", REG, "region."))
  
  ### SAVE OUTPUTS ###------------------------------------------------
  
  print("saving OBS outputs")

  nc_close(obs_region)
  nc_close(obs_national)
  rn <- rn0
#}

#stop("Only doing OBS")

 
#### EC SPLITTING ####---------------------------------------------------------
ec_national <- nc_open(paste0(data_wd, subfold, "eventEC_", thresh1,
                              "_", ws1, "_RCM", RCM, suffix, ".nc"))

# # PoE under different computations with extra data. Tidy format.
# present <- readr::read_csv(paste0(data_wd, subfold, "returnlevels_",
#                               thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))



### EVENT SPLITTING ###---------------------------------------------------
ST <- Sys.time()
ST0 <- Sys.time()
#for(REG in c("NW")){
  REG = "NW"
  print(paste("splitting", REG))
  
  if(!dir.exists(paste0(data_wd, subfold, "/", REG))){
    dir.create(paste0(data_wd, subfold, "/", REG))
  }
  
  #print(paste("splitting", REG))
  
  r1 <- which(rn_regions$REGION == REG) # length = 1437
  thresh_region <- threshMat[r1, jT]
  rn0 <- rn
  rn <- rn[r1,]  
  
  savepath <- paste0(data_wd,subfold, REG, "/",
                  "eventEC_region_", REG, "_RCM", RCM, suffix, ".nc")
  
  cdfPrimer(RCM, period, "EC2", NE=250, NH=length(r1), thresh1, ws1, rn, savepath)
  
  ec_region <- nc_open(savepath, write=T)
  
  rn <- rn0
  
  NE <- ec_national$dim$event$len
  NEreg <- 0
  thresh_region <- threshMat[r1, jT]
  ec_event_numbers <- ncvar_get(ec_national, "eventNo")
  WE <- c()
  for(i in 1:NE){
    
    if((i < 10) | (i %% 200 == 0)){ # time recording
    print(i)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((NE-i)/NE ,2)))
    print(paste("Time remaining", round((NE-i)/i * I0,2)))
    print(paste("Since last readout:", round(I,2)))
    ST <- Sys.time()
    print(ST) 
    }
    
    ec_events_region <- ncvar_get(ec_national, "flow",
                              start=c(1,i), count=c(-1,1))[r1]
    if(any(ec_events_region < param_table$threshold)){
      NEreg <- NEreg + 1
    
      ncvar_put(ec_region, "flow", ec_events_region,
          start=c(1, NEreg), count=c(-1,1))
      
      ec_dpe_region <- ncvar_get(ec_national, "dpe",
                       start=c(1,i), count=c(-1,1))[r1]
      ncvar_put(ec_region, "dpe", ec_dpe_region,
                start=c(1, NEreg), count=c(-1,1))
      
      ec_ape_region <- ncvar_get(ec_national, "ape",
                                 start=c(1,i), count=c(-1,1))[r1]
      ncvar_put(ec_region, "ape", ec_ape_region,
                start=c(1, NEreg), count=c(-1,1))
      WE <- c(WE, ec_event_numbers[i])
    }
    
  }
  ncvar_put(obs_region, "eventNo", vals=WE, start=1, count=length(WE))
#}


nc_close(ec_region)
nc_close(ec_national)

settings$HTsplit <- TRUE
write_yaml(settings, settingspath)

print(Sys.time())
print("108N done.")