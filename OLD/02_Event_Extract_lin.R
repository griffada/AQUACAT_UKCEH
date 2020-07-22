#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-02-21
#
# Event threshold extraction from daily flow.
# Practiced on one 30-year period of data.
# Rewritten for linux to reduce read-write times for netCDF.
#
# This takes about 9 hours to run on linux if you're lucky.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-02-21
#
#~~~~~~~~~~~~~~~~~~~~~~~

#### SETUP ####----------------------
library(ncdf4)
library(raster)

# threshold for inundation
threshVal <- c(5/365, 2/365, 1/365, 0.2/365, 0.1/365)
threshName <- c("POT5", "POT2", "POT1", "Q5", "Q10")
NT <- length(threshVal)

# cut-off for widespread event
wsCutoff <- c(0.05, 0.02, 0.01, 0.005, 0.001)
NW <- length(wsCutoff)


threshGridList <- list()

threshDayExcList <- vector("list", length(threshVal))

#### DATA ####-----------------------
ncnameLin <-"/prj/ResilRiskInds/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
ncnameWin <- "V:/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"  
# This file is ~36GB on the linux server.

if(substr(osVersion,1,3)=="Win"){
  ncname <- ncnameWin
}else{
  ncname <- ncnameLin
}
# This file is ~36GB on the linux server.
print(ST <- Sys.time())
ncin <- nc_open(ncname)
print(ncin)
print(Sys.time() - ST)


#### EVENT EXTRACT ####------------------------------------
## Get grid of river network ##-------------------
tStart <- 1
deltaT <- 5
ncwide <- ncvar_get(ncin, "dmflow", start=c(1,1,tStart), count=c(-1, -1, deltaT))

rn <- read.csv("./InterimData/hasData.csv")

#rn <- which(apply(ncwide, c(1,2), 
#                       function(v){sum(v[!is.na(v)] > -1) == deltaT}), arr.ind=T)
#readr::write_csv(data.frame(rn), "./InterimData/hasData.csv")

NH <- nrow(rn)
## Initialise dfs and lists ##-------------------------

threshGrid <- ncvar_get(ncin, "dmflow", start=c(1,1,tStart), count=c(-1, -1, 1))
threshMat <- matrix(NA, ncol=NT, nrow=NH)
# Reset grid to easily spot things
threshGrid[!is.na(threshGrid) & threshGrid > -1] <- 0

eastRange <- range(ncvar_get(ncin, "Easting"))
norRange <- range(ncvar_get(ncin, "Northing"))
threshGrid <- raster(threshGrid,
                     ymn=eastRange[1] - 500,
                     ymx=eastRange[2] + 500,
                     xmn=norRange[1] - 500,
                     xmx=norRange[2] + 500
                     )

threshGridList <- lapply(1:length(threshVal), function(i){threshGrid})

## Get quantiles for thresholds ##-------------------------------------
print("loop start")
ST <- Sys.time()
for(n in 1:NH){
  
  # running time diagnosis
  if(n %% 250 == 0){ print(Sys.time()-ST); ST <- Sys.time() }
  
  i <- rn[n,1]
  j <- rn[n,2]
  tSlice <- ncvar_get(ncin, varid="dmflow",
                       start=c(i, j,  1),
                       count=c(1, 1, -1))
  # find quantile for threshold
  thresh <- quantile(as.vector(tSlice), prob=c(1 - threshVal), na.rm=T)
  threshMat[n,] <- thresh
  for(k in 1:NT){
    threshGridList[[k]][i,j] <- thresh[k]
    # save which days cell n was exceeded.
    threshDayExcList[[k]][[n]] <- which(tSlice > thresh[k])
  }
}


## Save outputs ##-----------------------------------------------------
for(i in 1:NT){
  writeRaster(threshGridList[[i]],
              filename=paste0("./InterimData/",threshName[i],"_threshGrid.asc"),
              format="ascii", overwrite=TRUE)
}

saveRDS(threshDayExcList, file="./InterimData/threshDayExcList.rds")
saveRDS(threshMat, file="./InterimData/threshMat.rds")
readr::write_csv(x=data.frame(threshMat), path="./InterimData/threshMat.csv")
