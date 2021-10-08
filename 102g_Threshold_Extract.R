#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-02-21
#
# Event threshold extraction from daily flow for various thresholds.
# A GPA distribution is fitted to the peaks over these thresholds.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-02-21, edited 2020-05-04
# Pipeline version 2020-09-07
# Version b uses GLO 2021-03-23
# Revert to GPA, 2021-09-03
#
# OUTPUTS: threshDayExcList.rda,
#          thresMat.rda,
#          paramtable.csv
#
#~~~~~~~~~~~~~~~~~~~~~~~
print("running 102")
if(interactive()){commandArgs <- function(...){c("01","present")}}
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

library(ilaprosUtils)
library(extRemes)
library(lmomco)

if(FALSE){
  print("thresh* files already exist for 102. Proceeding to next job.")
}else{

suffix_pres  <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")

# thresh1 <- "POT2"
# ws1 <- "pc01"
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)

#prealloc.
threshDayExcList <- lapply(seq_len(length(threshVal)),
                           function(i){vector("list", NH)})


ST <-  Sys.time()
ncin <- nc_open(ncoriginal, readunlim=FALSE)
  # this is a huge file, do not open without reason.
print(Sys.time() - ST)
print(ncin)


##### THRESHOLD EXTRACT #####------------------------------------
## Get grid of river network ##-------------------
tStart <- 1
deltaT <- 1
NH <- nrow(rn)

##### Get quantiles for thresholds #####--------------------------------
print("loop start")
ST <- Sys.time()
print(ST)
ST0 <- ST
print(paste("NH =", NH))
threshMat <- readdf(paste0(data_wd, subfold_pres,
                          "threshMat_RCM", RCM, suffix_pres,".csv"))

if(TRUE){ ### PRESENT ###----------------------------------
  for(h in 1:NH){
    if((h < 10) | (h %% 200 == 0)){ # time recording
      print(h)
      I <- difftime(Sys.time(), ST, units="secs")
      I0 <- difftime(Sys.time(), ST0, units="secs")
      print(paste("Percent remaining", 100 * round((NH - h)/NH , 2)))
      print(paste("Time remaining", round((NH - h)/h * I0, 2)))
      print(paste("Since last readout:", round(I, 2)))
      ST <- Sys.time()
      print(ST) 
    }
    
    i <- rn$row[h]
    j <- rn$col[h]
    tSlice <- as.vector(ncvar_get(ncin, varid="dmflow",
                         start=c(i, j,  1),
                         count=c(1, 1, -1)))
    
    k <- jT
    # # save which days cell n was exceeded.
    threshDayExcList[[k]][[h]] <- which(tSlice > threshMat[h,jT])

  }
}

nc_close(ncin)

partable <- readdf(paste0(data_wd,subfold, 
                  "paramtableG_",thresh1, "_RCM", RCM, suffix, ".csv"))

partable$POT2 <- threshMat[,jT]

readr::write_csv(partable, paste0(data_wd,subfold, 
                          "paramtableG_",thresh1, "_RCM", RCM, suffix, ".csv"))

saveRDS(threshMat, file=paste0(data_wd, subfold,
                               "threshMat_RCM", RCM, suffix,".rds"))
readr::write_csv(data.frame(threshMat), paste0(data_wd, subfold,
                              "threshMat_RCM", RCM, suffix,".csv"))

## Save outputs ##-----------------------------------------------------
# Note this will make copies of some objects in the Future case,
# but this is safer
saveRDS(threshDayExcList, file=paste0(data_wd, subfold,
                          "threshDayExcList_RCM", RCM, suffix,".rds"))

settings$paramtable <- TRUE
settings$threshold_extract <- "102g"
write_yaml(settings, settingspath)
}
warnings()
print(Sys.time())
print("102 Complete")