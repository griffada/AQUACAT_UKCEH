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
#    eventdb_EC_***.sqlite,
#    eventdb_HT_***.sqlite
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
    stop(paste("incorrect call: Rscript 111_SQL_Compilation.R gcm period region \n",
               "- Region must be one of: ANG, NE, NW, SCO, SE, SEV, SW, THA, TRE, WAL."))
  }
}


library(dbplyr)
library(RSQLite)
library(reshape2)

##### DATA-------------------------------------------------------------
dir.create(paste0(data_wd,subfold, "/SQLite/"))

if(!is.na(REG)){
  
  rl_db_file <- paste0(data_wd,subfold, "/SQLite/eventdb_HT_region_", REG, "_RCM",
                       RCM, suffix, ".sqlite")
  
  rl_db <- src_sqlite(rl_db_file, create=TRUE)
  
  rn_regions$locnum <- 1:nrow(rn_regions)
  
  rn_reg <- data.frame(rn_regions) %>% dplyr::filter(REGION == REG)
  
  LNAMES <- c("eventflow_OBS_region", "eventdpe_OBS_region", "eventape_OBS_region",
              "eventflow_HT", "eventdpe_HT", "eventape_HT")
  
  for(l in LNAMES){
    x <- read_csv(paste0(data_wd, subfold, l, REG,"_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
    dimnames(x) <- list(paste0("L", rn_reg$locnum),
                        paste0("E", 1:ncol(x)))
    copy_to(rl_db, x, name = l, temporary=FALSE, overwrite=TRUE)
  }
  
  MODELS <- readRDS(paste0(data_wd, subfold, "marginal_models", REG, "_",
                           thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))
  COEFFS <- readRDS(paste0(data_wd, subfold, "coefficients_", REG, "_",
                           thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))
  ZSCORES <- readRDS(paste0(data_wd, subfold, "zScores", REG, "_",
                            thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))
  
  COEFF_tall <- data.frame(from=numeric(),
                           to=numeric(),
                           a=numeric(),
                           b=numeric(),
                           c=numeric(),
                           d=numeric(),
                           m=numeric(),
                           s=numeric())
  for(i in 1:dim(COEFFS)[3]){
    if(i %% 25 == 0){print(i)}
    for(j in 1:dim(COEFFS)[2]){
      row <- COEFFS[,j,i]
      row <- c(i, ifelse(i <= j, j+1, j), row)
      COEFF_tall[nrow(COEFF_tall)+1,] <- row
    }
  }
  
  MARGINALS <- t(sapply(MODELS, function(x){c(x$par, x$threshold)}))
  colnames(MARGINALS) <- c("scale", "shape", "threshold")
  
  copy_to(rl_db, MARGINALS, name='HT_marginals', temporary=FALSE, overwrite=TRUE)
  copy_to(rl_db, COEFF_tall, name="COEFF_tall", temporary=FALSE, overwrite=TRUE)
  
}else{
  
  rl_db_file <- paste0(data_wd, subfold, "/SQLite/eventdb_EC_RCM", 
                       RCM, suffix, ".sqlite")
  
  rl_db <- src_sqlite(rl_db_file, create=TRUE)
  
  LNAMES <- c("eventflow_OBS", "eventdpe_OBS", "eventape_OBS",
              "eventflow_EC", "eventdpe_EC", "eventape_EC")
  
  for(l in LNAMES){
    x <- read_csv(paste0(data_wd, subfold, l, "_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
    dimnames(x) <- list(paste0("L", 1:NH),
                        paste0("E", 1:ncol(x)))
    copy_to(rl_db, x, name = l, temporary=FALSE, overwrite=TRUE)
  }
}


