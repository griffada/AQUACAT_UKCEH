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

if(interactive()){commandArgs <- function(...){c("04","present", "NW")}}
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}
regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE", "SEV", "SSC", "SW", "THA", "TRE", "WAL")
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

# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                            RCM, suffix, ".rds"))
#dim(threshMat) #19914 x 5


#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))
for(REG in REGIONS){
  
  r1 <- which(rn_regions$REGION == REG)
  



obs_event_flow  <- as.data.frame(readr::read_csv(paste0(data_wd, subfold, REG,
                                                        "/eventflow_OBS_region_",
                                        REG,"_RCM", RCM, suffix,".csv"),
                               col_types=cols(.default = col_double()))[,-(1:4)])

obs_event_ape  <- as.data.frame(readr::read_csv(paste0(data_wd, subfold, REG,
                                                       "/eventape_OBS_region_",
                                        REG,"_RCM", RCM, suffix,".csv"),
                               col_types=cols(.default = col_double()))[,-(1:4)])

obs_event_dpe <- as.data.frame(readr::read_csv(paste0(data_wd, subfold,
                                                      REG, "/eventdpe_OBS_region_",
                                        REG,"_RCM", RCM, suffix,".csv"),
                                 col_types=cols(.default = col_double()))[,-(1:4)])
NN <- ncol(obs_event_flow)
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
NC <- ncol(obs_event_flow)
for(i in 1:(ncol(obs_event_flow))){
    if((i < 10) | (i %% 40 == 0)){ # time recording
    print(i)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((NC-i)/NC ,2)))
    print(paste("Time remaining", round((NC-i)/i * I0,2)))
    print(paste("Since last readout:", round(I,2)))
    ST <- Sys.time()
    print(ST) 
  }
  ni <- as.numeric(strsplit(colnames(obs_event_flow)[i],split="[E.]")[[1]][2])
  vvec <- unlist(obs_event_flow[,i])
  avec <- obs_event_ape[,i]
  dvec <- obs_event_dpe[,i]
  D <- eventDayList[[jT]][[jW]][ni]
  L <- eventLList[[jT]][[jW]][ni]
  ncl <- nclustevent(dvec < threshVal[jT])
  pk <- peakyfun(dvec, threshVal[jT])
  eventSummaryFrame[i,] <- list(i, D, L, sum(vvec > threshMat[r1,jT]),
                                min(avec), min(dvec), season(D), ncl, pk)
}

readr::write_csv(eventSummaryFrame, 
                 paste0(data_wd,subfold, REG, "/eventSumm_OBS_",REG,"_",
                        thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))



print(paste("OBS regional summary done."))














### EC SUMMARY ####
obs_event_flow  <- as.data.frame(readr::read_csv(paste0(data_wd, subfold, REG,
                                                        "/eventflow_EC_region_",
                                        REG,"_RCM", RCM, suffix,".csv"),
                               col_types=cols(.default = col_double()))[,-(1:4)])

obs_event_ape  <- as.data.frame(readr::read_csv(paste0(data_wd, subfold, REG,
                                                       "/eventape_EC_region_",
                                        REG,"_RCM", RCM, suffix,".csv"),
                               col_types=cols(.default = col_double()))[,-(1:4)])

obs_event_dpe <- as.data.frame(readr::read_csv(paste0(data_wd, subfold,
                                                      REG, "/eventdpe_EC_region_",
                                        REG,"_RCM", RCM, suffix,".csv"),
                                 col_types=cols(.default = col_double()))[,-(1:4)])
NN <- ncol(obs_event_flow)
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
NC <- ncol(obs_event_flow)
for(i in 1:(ncol(obs_event_flow))){
    if((i < 10) | (i %% 40 == 0)){ # time recording
    print(i)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((NC-i)/NC ,2)))
    print(paste("Time remaining", round((NC-i)/i * I0,2)))
    print(paste("Since last readout:", round(I,2)))
    ST <- Sys.time()
    print(ST) 
  }
  #print(i)
  ni <- as.numeric(substr(colnames(obs_event_flow)[i],2,6))
  vvec <- unlist(obs_event_flow[,i])
  avec <- unlist(obs_event_ape[,i])
  dvec <- unlist(obs_event_dpe[,i])
  D <- eventDayList[[jT]][[jW]][ni]
  L <- eventLList[[jT]][[jW]][ni]
  ncl <- nclustevent(avec < 1-exp(-360*threshVal[jT]))
  pk <- peakyfun(avec, 1-exp(-360*threshVal[jT]))
  eventSummaryFrame[i,] <- list(i, D, L, sum(vvec > threshMat[r1,jT]),
                                min(avec), min(dvec), season(D), ncl, pk)
}

readr::write_csv(eventSummaryFrame, paste0(data_wd,subfold, REG, "/eventSumm_EC_",REG,"_",
                                        thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))


print(paste("EC regional summary done."))






#### HT SUMMARY ####
if(file.exists(paste0(data_wd,subfold, REG, "/eventflow_HT_",REG,"_",
                                      thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))){

obs_event_flow  <- as.data.frame(readr::read_csv(paste0(data_wd,subfold, REG, "/eventflow_HT_",REG,"_",
                                      thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"),
                               col_types=cols(.default = col_double()))[,-(1:4)])

obs_event_ape  <- as.data.frame(readr::read_csv(paste0(data_wd,subfold, REG, "/eventape_HT_",REG,"_",
                                      thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"),
                               col_types=cols(.default = col_double()))[,-(1:4)])[1:length(r1),]

obs_event_dpe <- as.data.frame(readr::read_csv(paste0(data_wd,subfold, REG, "/eventdpe_HT_",REG,"_",
                                        thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"),
                                 col_types=cols(.default = col_double()))[,-(1:4)])[1:length(r1),]
NN <- ncol(obs_event_flow)
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
NC <- ncol(obs_event_flow)
for(i in 1:(ncol(obs_event_flow))){
    if((i < 10) | (i %% 40 == 0)){ # time recording
    print(i)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((NC-i)/NC ,2)))
    print(paste("Time remaining", round((NC-i)/i * I0,2)))
    print(paste("Since last readout:", round(I,2)))
    ST <- Sys.time()
    print(ST) 
    }
  #print(i)
  ni <- as.numeric(strsplit(colnames(obs_event_flow)[i],split="[E.]")[[1]][2])
  vvec <- unlist(obs_event_flow[,i])
  avec <- unlist(obs_event_ape[,i])
  dvec <- unlist(obs_event_dpe[,i])
  D <- eventDayList[[jT]][[jW]][ni]
  L <- eventLList[[jT]][[jW]][ni]
  ncl <- nclustevent(avec < 1-exp(-360*threshVal[jT]))
  pk <- peakyfun(avec, 1-exp(-360*threshVal[jT]))
  eventSummaryFrame[i,] <- list(i, D, L, sum(vvec > threshMat[r1,jT]),
                                min(avec), min(dvec), season(D), ncl, pk)
  
}
readr::write_csv(eventSummaryFrame, paste0(data_wd,subfold, REG, "/eventSumm_HT_",REG,"_",
                                        thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))
print("HT regional summary done.")
}
print(paste(REG, "region done."))
}