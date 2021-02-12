######
# Adam Griffin, 2020-04-27
#
# Summarising size and time, and extracting daily PoE of extreme events.
# 
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-04-27
# Pipeline version ABG 2020-09-07
#
#
# OUTPUTS: eventdf_***.csv: dataframe of each event summarised. One event per row.
#
#####

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}
regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE", "SEV", "SSC", "SW", "THA", "TRE", "WAL")


thresh1 <- "POT2"
ws1 <- "pc05"
# Only using POT2 and 2% inundation minimums.
jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)
print(paste("Running for threshold", thresh1, "at", ws1, "minimum spread."))

if(file.exists(paste0(data_wd,subfold,
                      "eventdf_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))){
  stop("eventdf_ exists. finishing 105.")
}

# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                            RCM, suffix, ".rds"))
#dim(threshMat) #19914 x 5


#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))



for(REG in REGIONS){

obs_event_flow  <- readr::read_csv(paste0(data_wd, subfold, "/", REG, "/eventflow_OBS_region_",
                                        REG,"_RCM", RCM, suffix,".csv"),
                               col_types=cols(.default = col_double()))[,-(1:4)]

obs_event_ape  <- readr::read_csv(paste0(data_wd, subfold, "/", REG, "/eventape_OBS_region_",
                                        REG,"_RCM", RCM, suffix,".csv"),
                               col_types=cols(.default = col_double()))[,-(1:4)]

obs_event_dpe <- readr::read_csv(paste0(data_wd, subfold, "/", REG, "/eventdpe_OBS_region_",
                                        REG,"_RCM", RCM, suffix,".csv"),
                                 col_types=cols(.default = col_double()))[,-(1:4)]
NN <- ncol(obs_event_flow)-4
eventSummaryFrame <- data.frame(eventNumber=numeric(NN),
                                eventDay=numeric(NN),
                                eventLength = numeric(NN),
                                area = numeric(NN),
                                peak = numeric(NN),
                                season = character(NN),
                                stringsAsFactors = FALSE)

# (D %/% 90) %% 4 for season 0,1,2,3: DJF, MAM, JJA, SON
season <- function(D) switch(((D %/% 90) %% 4) + 1, "DJF", "MAM", "JJA", "SON")
r1 <- which(rn_regions$REGION == REG)

for(i in 1:(ncol(obs_event_flow))){
  ni <- substr(colnames(obs_event_flow)[i],2,6)
  vvec <- obs_event_flow[,i]
  avec <- obs_event_ape[,i]
  D <- eventDayList[[jT]][[jW]][ni]
  L <- eventLList[[jT]][[jW]][ni]
  eventSummaryFrame[i,] <- list(i, D, L, sum(vvec > threshMat[r1, jT]), min(avec), season(D))
}

readr::write_csv(eventSummaryFrame, paste0(data_wd,subfold, "/", REG, "/eventSumm_OBS_region_",
                                        thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))


#### EC SUMMARY ####

obs_event_flow  <- readr::read_csv(paste0(data_wd,subfold, "/", REG, "/eventflow_HT_",
                                      thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"),
                               col_types=cols(.default = col_double()))[,-(1:4)]

obs_event_ape  <- readr::read_csv(paste0(data_wd,subfold, "/", REG, "/eventape_HT_",
                                      thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"),
                               col_types=cols(.default = col_double()))[,-(1:4)]

obs_event_dpe <- readr::read_csv(paste0(data_wd,subfold, "/", REG, "/eventdpe_HT_",
                                        thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"),
                                 col_types=cols(.default = col_double()))[,-(1:4)]
NN <- ncol(obs_event_flow)-4
eventSummaryFrame <- data.frame(eventNumber=numeric(NN),
                                eventDay=numeric(NN),
                                eventLength = numeric(NN),
                                area = numeric(NN),
                                peak = numeric(NN),
                                season = character(NN),
                                stringsAsFactors = FALSE)

# (D %/% 90) %% 4 for season 0,1,2,3: DJF, MAM, JJA, SON
season <- function(D) switch(((D %/% 90) %% 4) + 1, "DJF", "MAM", "JJA", "SON")

for(i in 1:(ncol(obs_event_flow)-4)){
  vvec <- obs_event_flow[,i]
  avec <- obs_event_ape[,i]
  D <- NA
  L <- NA
  eventSummaryFrame[i,] <- list(ni, D, L, sum(vvec > threshMat[r1,jT]), min(avec), season(D))
}

readr::write_csv(eventSummaryFrame, paste0(data_wd,subfold, "/", REG, "/eventSumm_EC_",
                                        thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))
}