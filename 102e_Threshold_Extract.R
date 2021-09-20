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

if(settings$paramtable){
  print("thresh* files already exist for 102. Proceeding to next job.")
}else{

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
threshMat <- matrix(NA, ncol=NT, nrow=NH)

partable <- data.frame(meanint=numeric(),
                        threshold=numeric(),
                        loc=numeric(),
                        scale=numeric(),
                        shape=numeric(),
                        threshquan=numeric(),
                        pot2=numeric())
NV <- 60
##### Get quantiles for thresholds #####--------------------------------
print("loop start")
ST <- Sys.time()
print(ST)
ST0 <- ST
print(paste("NH =", NH))

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
    
    # find quantile for threshold
    thresh <- quantile(tSlice, prob=c(1 - threshVal[jT]), na.rm=T) #vec
    mt <- medianThresh(tSlice, threshVal=(2/360))
    threshMat[h,jT] <- thresh
    
    k <- jT
    #for(k in 1:NT){
    # save which days cell n was exceeded.
    threshDayExcList[[k]][[h]] <- which(tSlice > threshMat[h,jT])
      
    # fit GPA to peaks above the threshold
    #if(k == jT){ # No point doing it for all 5 thresholds.
        
    thresh1b <- ecdf(tSlice)(mt)
        
    o1 <- try({
      ep1 <- (ilaprosUtils::extractPeaks(vecObs=tSlice, mintimeDiff=7) == 1)
    
      min_ep <- ecdf(tSlice)(
        sort(tSlice[which(ep1)], decreasing=T)[min(NV, sum(ep1))])
      #print(min_ep)
    })
    
    if(inherits(o1, "try-error") | min_ep < 0.8){
      print(paste("min_ep =", min_ep))
      print(paste(h, " Using UKFE2 and medianThresh"))

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
    thresh0 <- 0.9999 * 
      sort(tSlice[which(ep1)], decreasing=T)[min(NV, sum(ep1))]
    
    ep1 <- ep1 & (tSlice >= thresh0)
      
    #}
    peak_vals <- tSlice[ep1]
    meanint <- 30/sum(ep1)
      # fit a GPA distribution
    o <- try({
        LM <- lmoms(peak_vals) 
        LM$lambdas[1] <- median(peak_vals)
        LM$ratios[2] <- LM$lambdas[2]/LM$lambdas[1]
        at_site_glo <- parglo(LM)
        invisible({cdfglo(peak_vals, at_site_glo)})
        at_site_glo <- at_site_glo$para
        #lmomco uses opposite shape to fevd
    })
    if(inherits(o, "try-error")){
      warning(paste("ML used for site",h))
      at_site_gum <- fevd(x=peak_vals,
                          type="Gumbel")$results$par
    }
      
    threshquan <- ecdf(tSlice)(thresh0)
      
    partable[nrow(partable)+1,] <- 
                list(meanint,
                     thresh0,
                     at_site_glo[1], at_site_glo[2], at_site_glo[3], 
                     threshquan, threshMat[h,k])
  }
}

# 

print("saving outputs")
colnames(partable) <- c("meanint", "threshold", "loc", "sca", "shape",
                        "threshquan", "pot2")
str(partable)
str(data.frame(threshMat))
nc_close(ncin)

## Save outputs ##-----------------------------------------------------
# Note this will make copies of some objects in the Future case,
# but this is safer
saveRDS(threshDayExcList, file=paste0(data_wd, subfold,
                          "threshDayExcList_RCM", RCM, suffix,".rds"))
saveRDS(threshMat, file=paste0(data_wd, subfold,
                               "threshMat_RCM", RCM, suffix,".rds"))
readr::write_csv(data.frame(threshMat), paste0(data_wd, subfold,
                              "threshMat_RCM", RCM, suffix,".csv"))
readr::write_csv(partable, paste0(data_wd,subfold, 
                          "paramtableG_",thresh1, "_RCM", RCM, suffix, ".csv"))

settings$paramtable <- TRUE
write_yaml(settings, settingspath)
}
warnings()
print(Sys.time())
print("102 Complete")