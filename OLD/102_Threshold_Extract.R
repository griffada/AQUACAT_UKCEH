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
print("running 102")
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("C:/Users/adagri/Documents/AQUACAT_C/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

library(ilaprosUtils)
library(extRemes)
library(lmomco)

rn <- rn[17000:18999,]

# thresh1 <- "POT2"
# ws1 <- "pc01"
# # Only using POT2 and 2% inundation minimums.
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)


# Does this job need doing? ##
# if(file.exists(paste0(data_wd, subfold, "threshDayExcList_RCM", RCM, suffix,".rds")) &&
#    file.exists(paste0(data_wd, subfold, "threshMat_RCM", RCM, suffix,".rds")) &&
#    file.exists(paste0(data_wd, subfold, "threshMat_RCM", RCM, suffix,".csv")) &&
#    file.exists(paste0(data_wd, subfold, "paramtable_", thresh1, "_RCM", RCM, suffix, ".csv"))
#    ){
#     stop("thresh* files already exist for 102. Proceeding to next job.")
# }else{
#    print("Proceeding to job.")
# }

#prealloc.
#threshDayExcList <- vector("list", length(threshVal))
threshDayExcList <- lapply(seq_len(length(threshVal)), function(i){vector("list", NH)})
#names(threshDayExcList) <- threshName

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

medianThresh <- function(vals, threshVal=(2/360)){
  vq <- c()
  for(i in 1:30){
    v <- vals[1:(1*360) + (i-1)*360]
    vq[i] <- quantile(v, probs=1 - threshVal)
  }
  return(median(vq))
}



NH <- nrow(rn)
## Initialise dfs and lists ##-------------------------

# threshGrid <- ncvar_get(ncin, "dmflow",
#                         start=c(1, 1, 1),
#                         count=c(-1, -1, 1))

threshMat <- matrix(NA, ncol=NT, nrow=NH)
# Reset grid to easily spot things
#threshGrid[!is.na(threshGrid) & threshGrid > -1] <- 0

eastRange <- range(ncvar_get(ncin, "Easting"))
norRange <- range(ncvar_get(ncin, "Northing"))

partable <- data.frame(loc=numeric(), threshold=numeric(),
                       scale=numeric(), shape=numeric(),
                       threshquan=numeric())

## Get quantiles for thresholds ##-------------------------------------
print("loop start")
ST <- Sys.time()
print(ST)
ST0 <- ST
print(paste("NH =", NH))

if(period=="present"){ ### PRESENT ###----------------------------------
  for(h in 1:NH){
    if((h < 10) | (h %% 200 == 0)){ # time recording
      print(h)
      I <- difftime(Sys.time(), ST, units="secs")
      I0 <- difftime(Sys.time(), ST0, units="secs")
      print(paste("Percent remaining", 100*round((NH-h)/NH ,2)))
      print(paste("Time remaining", round((NH-h)/h * I0,2)))
      print(paste("Since last readout:", round(I,2)))
      ST <- Sys.time()
      print(ST) 
    }
    
    i <- rn$row[h]
    j <- rn$col[h]
    tSlice <- as.vector(ncvar_get(ncin, varid="dmflow",
                         start=c(i, j,  1),
                         count=c(1, 1, -1)))
    
    # find quantile for threshold
    thresh <- quantile(tSlice, prob=c(1 - threshVal), na.rm=T) #vec
    threshMat[h,] <- thresh
    
    for(k in 1:NT){
      # save which days cell n was exceeded.
      threshDayExcList[[k]][[h]] <- which(tSlice > thresh[k])
      
      # fit GPA to peaks above the threshold
      if(k == jT){ # No point doing it for all 5 thresholds.
        #o <- try({
          # get Peaks
        # o1 <- try({
        #   ep1 <- (ilaprosUtils::extractPeaks(vecObs=tSlice, mintimeDiff=7) == 1) &
        #     (tSlice > thresh[k])
        #   
        # })
        # if(inherits(o1, "try-error") | sum(ep1) < 30){
        #   warning(paste("using Median POT2 threshold for location", n))
        #   # If there are not enough peaks, use median thresh.
        o1 <- try({
          ep1 <- (ilaprosUtils::extractPeaks(vecObs=tSlice, mintimeDiff=7) == 1)
        })
        if(inherits(o1, "try-error")){
          print(paste(n, " Using UKFE2 and medianThresh"))
          mt <- medianThresh(tSlice, threshVal=0.05)
          thresh1b <- ecdf(tSlice)(mt)
          df <- data.frame(Date=1:10800, Flow=tSlice)
          
          o2 <- try({ep1 <- POT_extract_UKFE2(df,
                              div=quantile(tSlice[tSlice < mt], probs=0.67),
                              thresh=thresh1b, Plot=FALSE)
          })
          if(inherits(o2, "try-error")){
            print(n)
            stop("This still doesn't work.")
          }
          ep1 <- 1:10800 %in% ep1$Date
        }
        thresh0 <- sort(tSlice[which(ep1)], decreasing=T)[min(60, sum(ep1))]
        ep1 <- ep1 & (tSlice > thresh0)
          
        #}
        peak_vals <- tSlice[ep1]
        whattime <- which(ep1)
          # get inter-arrival time
        meanint <- mean(whattime[-1] - whattime[-length(whattime)])/360
          # fit a GPA distribution
        o <- try({
            at_site_gum <- pargum(lmoms(peak_vals))
            at_site_gpa <- pargpa(lmoms(peak_vals), xi=thresh0)
            a <- 1/(1- cdfgum(peak_vals, at_site_gum))
            b <- 1/(1- cdfgpa(peak_vals, at_site_gpa))
            invisible({cdfgum(peak_vals, at_site_gum)})
            at_site_gum <- at_site_gum$para
            #lmomco uses opposite shape to fevd
        })
        if(inherits(o, "try-error")){
          warning(paste("ML used for site",h))
          at_site_gum <- fevd(x=peak_vals,
                              type="Gumbel")$results$par
        }
          
        threshquan <- ecdf(tSlice)(thresh0)
          
        partable[nrow(partable)+1,] <- 
                    list(meanint, thresh0, at_site_gum[1], at_site_gum[2], threshquan)
          
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
  partable <- as.data.frame(readr::read_csv(paste0(data_wd,subfold_pres, 
                              "paramtable_",thresh1, "_RCM", RCM, suffix_pres, ".csv")))
  str(threshMat)
  # Just load in threshGridList to save to the Future folder
  # threshGridList <- readRDS(file=paste0(data_wd, subfold_pres,
  #                               "threshGridList_RCM", RCM, suffix_pres,".rds"))
  
  for (n in 1:NH) {
    if ((n < 10) | (n %% 200 == 0)) { # time recording
      print(n)
      I <- difftime(Sys.time(), ST, units="secs")
      I0 <- difftime(Sys.time(), ST0, units="secs")
      print(paste("Percent remaining", 100*round((NH-n)/NH ,2)))
      print(paste("Time remaining", round((NH-n)/n * I0,2)))
      print(paste("Since last readout:", round(I,2)))
      ST <- Sys.time()
      print(ST) 
    }
    
    i <- rn$row[n]
    j <- rn$col[n]
    tSlice <- as.vector(ncvar_get(ncin, varid="dmflow",
                        start=c(i,j,1), count=c(1,1,-1)))
    thresh <- threshMat[n,]
    for(k in 1:NT){
      # save which days cell n was exceeded.
      threshDayExcList[[k]][[n]] <- which(tSlice > thresh[k])
    }
  }
}
print("saving outputs")
colnames(partable) <- c("meanint", "threshold", "scale", "shape", "threshquan")
str(partable)
str(data.frame(threshMat))
nc_close(ncin)

## Save outputs ##-----------------------------------------------------
# Note this will make copies of some objects in the Future case, but this is safer
# saveRDS(threshDayExcList, file=paste0(data_wd, subfold,
#                           "threshDayExcList_RCM", RCM, suffix,".rds"))
# saveRDS(threshMat, file=paste0(data_wd, subfold,
#                                "threshMat_RCM", RCM, suffix,".rds"))
# readr::write_csv(data.frame(threshMat), paste0(data_wd, subfold,
#                               "threshMat_RCM", RCM, suffix,".csv"))
readr::write_csv(partable, paste0(data_wd,subfold, 
                              "paramtableG_",thresh1, "_RCM", RCM, suffix, ".csv"))

warnings()
print("102 Complete")