######
# Adam Griffin, 2021-04-19
#
# Summarising size and time into single documents for full analysis.
# 
# For aquaCAT, Project 07441.
#
# OUTPUTS: eventdfALL_***.csv: dataframe of each event summarised. One event per row.
#
#####

RCMS <- c("01","04","05","06","07","08","09","10","11","12","13","15")
PERIODS <- c("present", "future")

commandArgs <- function(...){c(RCMS[1], PERIODS[1])}

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

notfirst <- FALSE
for(RCM in RCMS){
  for(period in PERIODS){
    
    suffix <- ifelse(period=="present","_198012_201011","_205012_208011")
    subfold <- paste0("RCM", RCM, suffix, "/")
    
    df_temp <- readr::read_csv(paste0(data_wd,subfold, "eventSumm_OBSB_",
                                        thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"), 
                          col_types=cols(
                            eventNumber = col_double(),
                            eventDay = col_double(),
                            eventLength = col_double(),
                            area = col_double(),
                            peakA = col_double(),
                            peakD = col_double(),
                            season = col_character(),
                            nclusters = col_double(),
                            peakyness = col_double()
                          ))
    
    df_temp$period <- rep(period, nrow(df_temp))
    df_temp$rcm <- rep(RCM, nrow(df_temp))
    if(notfirst){
      df <- rbind(df, df_temp)
    }else{
      df <- df_temp
    }
    notfirst <- TRUE
  }
}
    readr::write_csv(df, paste0(data_wd, "eventSumm_OBSB_", thresh1,
                                "_", ws1, "_ALL.csv"))
    
notfirst <- FALSE
for(RCM in RCMS){
  for(period in PERIODS){
    
    suffix <- ifelse(period=="present","_198012_201011","_205012_208011")
    subfold <- paste0("RCM", RCM, suffix, "/")    
    df_temp <- readr::read_csv(paste0(data_wd,subfold, "eventSumm_ECB_",
                                        thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"),
                          col_types=cols(
                            eventNumber = col_double(),
                            eventDay = col_double(),
                            eventLength = col_double(),
                            area = col_double(),
                            peakA = col_double(),
                            peakD = col_double(),
                            season = col_character(),
                            nclusters = col_double(),
                            peakyness = col_double()
                          ))

    df_temp$period <- rep(period, nrow(df_temp))
    df_temp$rcm <- rep(RCM, nrow(df_temp))
    if(notfirst){
      df <- rbind(df, df_temp)
    }else{
      df <- df_temp
    }
    notfirst <- TRUE
    

  }
}
    readr::write_csv(df,paste0(data_wd, "eventSumm_ECB_", thresh1, "_",
                           ws1, "_ALL.csv"))


print("National eventSumm compiled for OBS and EC.")
# 
# notfirst <- FALSE
# for(RCM in RCMS){
#   for(period in PERIODS){
#     REG <- "NW"
#     suffix <- ifelse(period=="present","_198012_201011","_205012_208011")
#     subfold <- paste0("RCM", RCM, suffix, "/")
#     if(file.exists(paste0(data_wd,subfold, REG, "/eventSumm_OBS_",REG,"_",
#                                           thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))){
#     df <- readr::read_csv(paste0(data_wd,subfold, REG, "/eventSumm_OBS_",REG,"_",
#                                           thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))
# 
#     df$period <- rep(period, nrow(df))
#     df$rcm <- rep(RCM, nrow(df))
#     df$region <- rep(REG , nrow(df))
# 
#     readr::write_csv(df,
#       paste0(data_wd, "eventSumm_OBS_region_", REG, "_", thresh1, "_", ws1, "_ALL.csv"),
#              append=notfirst, col_names=!notfirst)
# 
#       df <- readr::read_csv(paste0(data_wd,subfold, REG, "/eventSumm_EC_",REG,"_",
#                                           thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))
# 
#     df$period <- rep(period, nrow(df))
#     df$rcm <- rep(RCM, nrow(df))
#     df$region <- rep(REG , nrow(df))
# 
#     readr::write_csv(df,
#       paste0(data_wd, "eventSumm_EC_region_", REG, "_", thresh1, "_", ws1, "_ALL.csv"),
#              append=notfirst, col_names=!notfirst)
#     
#     notfirst <- TRUE
#     }
#   }
# }
# 
# print("Regional eventSumm compiled for OBS and EC.")
# 
# notfirst <- FALSE
# for(RCM in RCMS){
#   for(period in PERIODS){
#     REG <- "NW"
#     suffix <- ifelse(period=="present","_198012_201011","_205012_208011")
#     subfold <- paste0("RCM", RCM, suffix, "/")
#     if(file.exists(paste0(data_wd,subfold, REG, "/eventSumm_HT_",REG,"_",
#                                         thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))){
#     df <- readr::read_csv(paste0(data_wd,subfold, REG, "/eventSumm_HT_",REG,"_",
#                                         thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))
# 
#     df$period <- rep(period, nrow(df))
#     df$rcm <- rep(RCM, nrow(df))
#     df$region <- rep(REG , nrow(df))
# 
#     readr::write_csv(df,
#       paste0(data_wd, "eventSumm_HT_", REG, "_", thresh1, "_", ws1, "_ALL.csv"),
#              append = notfirst, col_names = !notfirst)
#     
#     notfirst <- TRUE
#     }
#   }
# }
# 
# print("Regional eventSumm compiled for HT.")