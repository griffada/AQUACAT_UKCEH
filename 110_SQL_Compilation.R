#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-09-10
#
# SQL database construction for HT/EC modelled events.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-09-10
#
# Outputs:
#    ***.sqlite
#
#~~~~~~~~~~~~~~~~~~~~~~~



##### SETUP -------------------------------------------------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

regions <- c("ANG", "NE", "NW", "SCO", "SE", "SEV", "SW", "THA", "TRE", "WAL")
REG <- NA
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
    stop(paste("incorrect call: Rscript 109_HeffTawn_Modelling.R gcm period region \n",
               "- Region must be one of: ANG, NE, NW, SCO, SE, SEV, SW, THA, TRE, WAL."))
  }
}


library(dbplyr)
library(RSQLite)
library(reshape2)

##### DATA-------------------------------------------------------------
dir.create(paste0(data_wd,subfold, "/SQLite/"))

if(!is.na(REG)){
  PRES_REG <- readr::read_csv(paste0(data_wd,subfold, "eventdf_region_",
                          REG,"_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
  dimnames(PRES_GB) <- list(paste0("LR", 1:ncol(PRES_REG)),
                            paste0("E", 1:nrow(PRES_REG)))
  
  NEW_HT <- readr::read_csv(paste0(data_wd, subfold, "NewEventHT_",REG,"_",
                                   thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
  dimnames(NEW_HT) <- list(paste0("LR", 1:ncol(NEW_HT)),
                           paste0("D", 1:nrow(NEW_HT)))
  

  rl_db_file <- paste0(data_wd,subfold, "/SQLite/eventdb_HT_region_", REG, "_RCM",
                       RCM, suffix, ".sqlite")
  
  rl_db <- src_sqlite(rl_db_file, create=TRUE)
  
  copy_to(rl_db,
          PRES_REG,
          temporary=FALSE,
          overwrite=TRUE)
  
  copy_to(rl_db,
          NEW_HT,
          temporary=FALSE,
          overwrite=TRUE)
  
  
}else{
  
  PRES_GB <- readr::read_csv(paste0(data_wd,subfold, "eventdf_",
                                    thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
  dimnames(PRES_GB) <- list(paste0("L", 1:NREG),
                            paste0("E", 1:nrow(PRES_GB)))
  
  NEW_EC <- readr::read_csv(paste0(data_wd, subfold, "NewEventEC_",
                                   thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
  dimnames(NEW_EC) <- list(paste0("L", 1:ncol(NEW_EC)),
                           paste0("D", 1:nrow(NEW_EC)))
  
  rl_db_file <- paste0(data_wd, subfold, "/SQLite/eventdb_EC_RCM", 
                       RCM, suffix, ".sqlite")
  
  rl_db <- src_sqlite(rl_db_file, create=TRUE)
  
  copy_to(rl_db,
          PRES_GB,
          temporary=FALSE,
          overwrite=TRUE)
  
  copy_to(rl_db,
          NEW_EC,
          temporary=FALSE,
          overwrite=TRUE)
}


