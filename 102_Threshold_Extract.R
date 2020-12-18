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
#
# OUTPUTS: threshDayExcList.rda,
#          thresMat.rda,
#          paramtable.csv
#
#~~~~~~~~~~~~~~~~~~~~~~~

#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~ /AQUACAT/CodeABG/setup_script_00.R")
}

library(ilaprosUtils)
library(extRemes)
library(lmomco)

thresh1 <- "POT2"
ws1 <- "pc05"
# Only using POT2 and 2% inundation minimums.
jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)


## Does this job need doing? ##
if(file.exists(paste0(data_wd, subfold, "threshDayExcList_RCM", RCM, suffix,".rds")) &&
   file.exists(paste0(data_wd, subfold, "threshMat_RCM", RCM, suffix,".rds")) &&
   file.exists(paste0(data_wd, subfold, "threshMat_RCM", RCM, suffix,".csv")) &&
   file.exists(paste0(data_wd, subfold, "paramtable_", thresh1, "_RCM", RCM, suffix, ".csv"))
   ){
    stop("thresh* files already exist for 102. Proceeding to next job.")
}else{
   print("Proceeding to job.")
}

#prealloc.
threshDayExcList <- vector("list", length(threshVal))
names(threshDayExcList) <- threshName

ST <-  Sys.time()
ncin <- nc_open(ncoriginal, readunlim=FALSE)
  # this is a huge file, do not open without reason.
print(Sys.time() - ST)
print(ncin)


#### THRESHOLD EXTRACT ####------------------------------------
## Get grid of river network ##-------------------
tStart <- 1
deltaT <- 1
###### Only uncomment if needed #####
# ncwide <- ncvar_get(ncin, "dmflow",
# 					start=c(1, 1, tStart),
# 					count=c(-1, -1, deltaT))
# 
# rn <- which(apply(ncwide, c(1,2),
#                  function(v){sum(v[!is.na(v)] > -1) == deltaT}), arr.ind=T)
# 
# rn <- data.frame(row=rn[,1], col=rn[,2])
# east <- ncvar_get(ncin, "Easting")
# north <- ncvar_get(ncin, "Northing")
# rn$east <- east[rn[,1]]
# rn$nor <- north[rn[,2]]
################################### 

NH <- nrow(rn)
## Initialise dfs and lists ##-------------------------

threshGrid <- ncvar_get(ncin, "dmflow",
                        start=c(1, 1, tStart),
                        count=c(-1, -1, 1))

threshMat <- matrix(NA, ncol=NT, nrow=NH)
# Reset grid to easily spot things
threshGrid[!is.na(threshGrid) & threshGrid > -1] <- 0

eastRange <- range(ncvar_get(ncin, "Easting"))
norRange <- range(ncvar_get(ncin, "Northing"))

partable <- data.frame(loc=numeric(), threshold=numeric(),
                       scale=numeric(), shape=numeric())

## Get quantiles for thresholds ##-------------------------------------
print("loop start")
ST <- Sys.time()
print(ST)
ST0 <- ST
print(paste("NH =", NH))

if(period=="present"){ ### PRESENT ###----------------------------------
  for(n in 1:NH){

    if((n < 10) | (n %% 200 == 0)){ # time recording
      print(n)
      I <- difftime(Sys.time(), ST, units="secs")
      I0 <- difftime(Sys.time(), ST0, units="secs")
      print(paste("Percent remaining", 100*round((NH-n)/NH ,2)))
      print(paste("Time remaining", round((NH-n)/n * I0,2)))
      print(paste("Since last readout:", round(I,2)))
      ST <- Sys.time()
      print(ST) 
    }
    
    i <- rn[n,1]
    j <- rn[n,2]
    tSlice <- as.vector(ncvar_get(ncin, varid="dmflow",
                         start=c(i, j,  1),
                         count=c(1, 1, -1)))
    
    # find quantile for threshold
    thresh <- quantile(tSlice, prob=c(1 - threshVal), na.rm=T) #vec
    threshMat[n,] <- thresh
    
    for(k in 1:NT){

      # save which days cell n was exceeded.
      threshDayExcList[[k]][[n]] <- which(tSlice > thresh[k])
      
      # fit GPA to peaks above the threshold
      if(k == jT){
        o <- try({
          
          ep1 <- (extractPeaks(vecObs=tSlice, mintimeDiff=7) == 1) & (tSlice > thresh[k])
          peak_vals <- tSlice[ep1]
          whattime <- which(ep1)
          meanint <- mean(whattime[-1] - whattime[-length(whattime)])/360
          at_site_gpa <- fevd(x=peak_vals,
                              threshold = thresh[k], type="GP")$results$par
          
          partable[nrow(partable)+1,] <- list(meanint, thresh[k], at_site_gpa[1], at_site_gpa[2])
          
          1
        })
        if(inherits(o, "try-error")){
          peak_vals <- blockmaxxer(data.frame(tSlice),blocks=rep(1:30,each=360))
          at_site_gpa <- pargpa(lmoms(peak_vals))
          
          partable[nrow(partable)+1,] <- list(1, at_site_gpa$para[1],
                               at_site_gpa$para[2], at_site_gpa$para[3])
        }
      }
    }
  }
}

if(period=="future"){  #### FUTURE ####-------------------------------------
  suffix_pres <- "_198012_201011"
  subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
  #use the Present day threshold matrix
  threshMat <- readRDS(paste0(data_wd, subfold_pres,
                          "threshMat_RCM", RCM, suffix_pres,".rds"))
  # Just load in threshGridList to save to the Future folder
  # threshGridList <- readRDS(file=paste0(data_wd, subfold_pres,
  #                               "threshGridList_RCM", RCM, suffix_pres,".rds"))
  
  for(n in 1:3){
    
    if((n < 10) | (n %% 200 == 0)){ # time recording
      print(n)
      I <- difftime(Sys.time(), ST, units="secs")
      I0 <- difftime(Sys.time(), ST0, units="secs")
      print(paste("Percent remaining", 100*round((NH-n)/NH ,2)))
      print(paste("Time remaining", round((NH-n)/n * I0,2)))
      print(paste("Since last readout:", round(I,2)))
      ST <- Sys.time()
      print(ST) 
    }
    
    
    i <- rn[n,1]
    j <- rn[n,2]
    tSlice <- as.vector(ncvar_get(ncin, varid="dmflow",
                        start=c(i,j,1), count=c(1,1,-1)))
    
    for(k in 1:NT){
      # Get exceedences
      threshDayExcList[[k]][[n]] <- which(tSlice > threshMat[n,k])
      if(k == jT){
        o <- try({
          # Take peaks from exceedences
          ep1 <- (extractPeaks(vecObs=tSlice, mintimeDiff=7) == 1) & 
            (tSlice > threshMat[n,k])
          
          peak_vals <- tSlice[ep1]
          whattime <- which(ep1)
          meanint <- mean(whattime[-1] - whattime[-length(whattime)])/360
          
          # fit GPA to peaks
          at_site_gpa <- fevd(x=peak_vals,
                              threshold = threshMat[n,k], type="GP")$results$par
          
          partable[nrow(partable)+1,] <- list(meanint, threshMat[n,k],
                               at_site_gpa[1], at_site_gpa[2])
          
          1
        })
        if(inherits(o, "try-error")){
          # If it breaks, fit to the AMAX series.
          peak_vals <- blockmaxxer(data.frame(tSlice),blocks=rep(1:30,each=360))
          
          at_site_gpa <- pargpa(lmoms(peak_vals))
          
          partable[nrow(partable)+1,] <- list(1, at_site_gpa$para[1],
                               at_site_gpa$para[2], at_site_gpa$para[3])
        }
      }
    }
  }
}
colnames(partable) <- c("meanint", "threshold", "scale", "shape")
nc_close(ncin)

## Save outputs ##-----------------------------------------------------

# Note this will make copies of some objects in the Future case, but this is safer
saveRDS(threshDayExcList, file=paste0(data_wd, subfold,
                          "threshDayExcList_RCM", RCM, suffix,".rds"))
saveRDS(threshMat, file=paste0(data_wd, subfold,
                               "threshMat_RCM", RCM, suffix,".rds"))
readr::write_csv(data.frame(threshMat), paste0(data_wd, subfold,
                              "threshMat_RCM", RCM, suffix,".csv"))
readr::write_csv(partable, paste0(data_wd,subfold, 
                              "paramtable_",thresh1, "_RCM", RCM, suffix, ".csv"))

warnings()
print("102 Complete")