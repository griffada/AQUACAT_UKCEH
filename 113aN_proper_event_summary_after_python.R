if(interactive()){commandArgs <- function(...){c("01","future")}}
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

season <- function(D) switch(((D %/% 90) %% 4) + 1, "DJF", "MAM", "JJA", "SON")

load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))

eventSummaryFrame <- readdf(paste0(data_wd,subfold, "eventSumm_OBS_",
                                  thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))

eventSummaryFrame$eventDay <- eventDayList[[jT]][[jW]]
eventSummaryFrame$eventLength <- eventLList[[jT]][[jW]]
eventSummaryFrame$season <- sapply(eventDayList[[jT]][[jW]], season)

print("added Day/Length/Season to OBS_Summ")

readr::write_csv(eventSummaryFrame, paste0(data_wd,subfold, "eventSumm_OBS_",
                                  thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))

ecSummaryFrame <- readdf(paste0(data_wd,subfold, "eventSumm_EC_",
                                  thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))

ecSummaryFrame$eventDay <- sapply(ecSummaryFrame$eventNumber,
                                  function(i){eventDayList[[jT]][[jW]][i]})
ecSummaryFrame$eventLength <- sapply(ecSummaryFrame$eventNumber,
                                     function(i){eventLList[[jT]][[jW]][i]})
ecSummaryFrame$season <- sapply(ecSummaryFrame$eventDay, season)

readr::write_csv(ecSummaryFrame, paste0(data_wd,subfold, "eventSumm_EC_",
                                  thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))

print("added Day/Length/Season to EC_Summ")
print("113a done.")