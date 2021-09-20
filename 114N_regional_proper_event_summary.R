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

if(interactive()){commandArgs <- function(...){c("01","present", "NW")}}

#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

if(settings$HTsumm){
  stop("Regional summaries exist. Stopping 114N.")
}

regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE", "SEV", "SSC", "SW",
             "THA", "TRE", "WAL")
REGIONS <- c("NW")

rn <- as.data.frame(rn)
r1 <- which(rn_regions$REGION == REGIONS[1])

nclustevent <- function(v, radii=50) {
    # This function estimates the number of "seperate clusters" of event going on during
    # a widespread event. 
    ne0 <- sum(v, na.rm=T)
    event1 <- rbind(as.matrix(rn[r1,][v,]), matrix(c(
      700, 0, 700000, 1000000
    ), nrow = 1))
    d <- dist(event1[, 1:2])
    d2 <- dist(event1[1:ne0, 1:2])
    if ((max(d2) > radii) & (ne0 >= 10)) {
      hh <- hclust(d, method = "single")
      ch <- max(cutree(hh, h = 0.8 * sort(hh$height, decreasing = T)[2])) - 1
      if(sort(hh$height,decreasing=T)[2] < 0.2*max(hh$height)){
        ch <- 1
      }
      ch <- min(ch, 4)
    } else{
      ch <- 1
    }
    ch
  }
  
peakyfun <- function(v, thr) {
  mn <- min(v, na.rm=T)
  if (mn > thr | mn > 1e8) {
    return(0)
  }else{
    mn50 <- mean(c(mn, thr))
    return(sum(v < mn50) / sum(v < thr))
  }
}

# thresh1 <- "POT2"
# ws1 <- "pc01"
# # Only using POT2 and 2% inundation minimums.
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)
print(paste("Running for threshold", thresh1, "at", ws1, "minimum spread."))

eventSummaryLine <- function(obs_events, i, threshold, r1=1:19914){
  ni <- eventNo[i]
  vvec <- sum(
    (ncvar_get(obs_events, "flow", start=c(1,i), count=c(-1,1)) > threshold),
    na.rm=T)
  avec <- ncvar_get(obs_events, "ape", start=c(1,i), count=c(-1,1))
  min_dvec <- min(
    ncvar_get(obs_events, "dpe", start=c(1,i), count=c(-1,1)),
    na.rm=T)
  D <- eventDayList[[jT]][[jW]][ni]
  L <- eventLList[[jT]][[jW]][ni]
  ncl <- nclustevent(avec < (1-exp(-2)))
  pk <- peakyfun(avec, (1-exp(-2)))

  return(list(i, D, L, vvec, min(avec), min_dvec, season(D), ncl, pk))
}

threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                            RCM, suffix, ".rds"))
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))

param_table <- readdf(paste0(data_wd, subfold, "paramtableG_POT2_RCM",
                            RCM, suffix, ".csv"))


# for(REG in c("NW")){
REG <- "NW"
  ##### REGIONAL OBS SUMMARY #####-----------------------------------------------
  
r1 <- which(rn_regions$REGION == REG)

param_table <- param_table[r1,]
  
obs_events <- nc_open(paste0(data_wd,subfold, REG, "/",
                    "eventOBS_region_", REG, "_RCM", RCM, suffix, ".nc"))
  
NN <- obs_events$dim$event$len
  
eventSummaryFrame <- data.frame(eventNumber=numeric(NN),
                                  eventDay=numeric(NN),
                                  eventLength = numeric(NN),
                                  area = numeric(NN),
                                  peakA = numeric(NN),
                                  peakD = numeric(NN),
                                  season = character(NN),
                                  nclusters = numeric(NN),
                                  peakyness = numeric(NN),
                                  stringsAsFactors = FALSE)
  
  # (D %/% 90) %% 4 for season 0,1,2,3: DJF, MAM, JJA, SON
season <- function(D) switch(((D %/% 90) %% 4) + 1,
                             "DJF", "MAM", "JJA", "SON")
ST <- Sys.time()
ST0 <- Sys.time()

eventNo <- ncvar_get(obs_events, "eventNo")

for(i in 1:NN){
    if((i < 10) | (i %% 40 == 0)){ # time recording
    print(i)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((NN-i)/NN ,2)))
    print(paste("Time remaining", round((NN-i)/i * I0,2)))
    print(paste("Since last readout:", round(I,2)))
    ST <- Sys.time()
    print(ST) 
    }
  try({
  eventSummaryFrame[i,] <- eventSummaryLine(obs_events, i, threshMat[r1,jT])
  })
}

readr::write_csv(eventSummaryFrame, 
                 paste0(data_wd,subfold, REG, "/eventSumm_OBS_",REG,"_",
                        thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))

print(paste("OBS regional summary done."))














### EC SUMMARY ####-------------------------------------------------------------
obs_events <- nc_open(paste0(data_wd,subfold, "/", REG, "/",
                  "eventEC_region_", REG, "_RCM", RCM, suffix, ".nc"))

NN <- obs_events$dim$event$len
eventSummaryFrame <- data.frame(eventNumber=numeric(NN),
                                eventDay=numeric(NN),
                                eventLength = numeric(NN),
                                area = numeric(NN),
                                peakA = numeric(NN),
                                peakD = numeric(NN),
                                season = character(NN),
                                nclusters = numeric(NN),
                                peakyness = numeric(NN),
                                stringsAsFactors = FALSE)

# (D %/% 90) %% 4 for season 0,1,2,3: DJF, MAM, JJA, SON
season <- function(D) switch(((D %/% 90) %% 4) + 1, "DJF", "MAM", "JJA", "SON")
r1 <- which(rn_regions$REGION == REG)

ST <- Sys.time()
ST0 <- Sys.time()
for(i in 1:NN){
  if((i < 10) | (i %% 40 == 0)){ # time recording
    print(i)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((NN-i)/NN ,2)))
    print(paste("Time remaining", round((NN-i)/i * I0,2)))
    print(paste("Since last readout:", round(I,2)))
    ST <- Sys.time()
    print(ST) 
  }
  #print(i)
  try({
  eventSummaryFrame[i,] <- eventSummaryLine(obs_events, i, threshMat[r1,jT])
  })
}

readr::write_csv(eventSummaryFrame, paste0(data_wd,subfold, REG,
                "/eventSumm_EC_",REG,"_", thresh1,"_", ws1,
                "_RCM", RCM, suffix, ".csv"))

print(paste("EC regional summary done."))










#### HT SUMMARY ####------------------------------------------------------------
if(settings$HTflow){
  
  obs_events <- nc_open(paste0(data_wd,subfold, REG, "/eventHT_",
                        REG,"_",thresh1,"_", ws1,"_RCM", RCM, suffix, ".nc"))
  
  NN <- obs_events$dim$event$len
  eventSummaryFrame <- data.frame(eventNumber=numeric(NN),
                                  eventDay=numeric(NN),
                                  eventLength = numeric(NN),
                                  area = numeric(NN),
                                  peakA = numeric(NN),
                                  peakD = numeric(NN),
                                  season = character(NN),
                                  nclusters = numeric(NN),
                                  peakyness = numeric(NN),
                                  stringsAsFactors = FALSE)
  
  # (D %/% 90) %% 4 for season 0,1,2,3: DJF, MAM, JJA, SON
  season <- function(D) switch(((D %/% 90) %% 4) + 1, "DJF", "MAM", "JJA", "SON")
  
  ST <- Sys.time()
  ST0 <- Sys.time()
  for(i in 1:NN){
    if((i < 10) | (i %% 40 == 0)){ # time recording
      print(i)
      I <- difftime(Sys.time(), ST, units="secs")
      I0 <- difftime(Sys.time(), ST0, units="secs")
      print(paste("Percent remaining", 100*round((NN-i)/NN ,2)))
      print(paste("Time remaining", round((NN-i)/i * I0,2)))
      print(paste("Since last readout:", round(I,2)))
      ST <- Sys.time()
      print(ST) 
      }
    #print(i)
    try({
    eventSummaryFrame[i,] <- eventSummaryLine(obs_events, i, threshMat[r1,jT])
    })
  }
  readr::write_csv(eventSummaryFrame,
                   paste0(data_wd,subfold, REG, "/eventSumm_HT_",REG,"_",
                          thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))
  print("HT regional summary done.")
}
print(paste(REG, "region done."))


settings$HTsumm <- TRUE
write_yaml(settings, settingspath)

print("114N done.")