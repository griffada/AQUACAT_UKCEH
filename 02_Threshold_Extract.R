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
# Created ABG 2020-02-21, edited 2020-05-04
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
if(substr(osVersion,1,3) == "Win"){
  #data_wd0 <- "S:/run_hmfg2g/outputs/"
  data_wd <- "S:/CodeABG/InterimData/"
  wd <- "S:/CodeABG/"
}else{
  #data_wd <- "/prj/aquacat/run_hmfg2g/outputs/"
  data_wd <- "/prj/aquacat/CodeABG/InterimData/"
  wd <- "/prj/aquacat/CodeABG/"
}

#ncname <- "dmflow_RCM01_198012_201011_out.nc" #2GB
ncname <- "dmflow_copy.nc"
ST <-  Sys.time()
ncin <- nc_open(paste0(data_wd,ncname), readunlim=FALSE) 
# this is a huge file, do not open without reason.
print(Sys.time() - ST)
print(ncin)


#### THRESHOLD EXTRACT ####------------------------------------
## Get grid of river network ##-------------------
tStart <- 1
deltaT <- 5
ncwide <- ncvar_get(ncin, "dmflow",
					start=c(1,1,tStart),
					count=c(-1, -1, deltaT))

#rn <- read.csv(paste0(wd,"InterimData/hasData.csv"))

rn <- which(apply(ncwide, c(1,2),
                 function(v){sum(v[!is.na(v)] > -1) == deltaT}), arr.ind=T)
readr::write_csv(data.frame(rn), paste0(wd,"InterimData/hasData2.csv"))

NH <- nrow(rn)
## Initialise dfs and lists ##-------------------------

threshGrid <- ncvar_get(ncin, "dmflow", start=c(1,1,tStart), count=c(-1, -1, 1))
threshMat <- matrix(NA, ncol=NT, nrow=NH)
# Reset grid to easily spot things
threshGrid[!is.na(threshGrid) & threshGrid > -1] <- 0

eastRange <- range(ncvar_get(ncin, "Easting"))
norRange <- range(ncvar_get(ncin, "Northing"))
#!#!#!
threshGrid <- raster(threshGrid,  # threshgrid might need transposing?
                     ymn=eastRange[1] - 500,
                     ymx=eastRange[2] + 500,
                     xmn=norRange[1] - 500,
                     xmx=norRange[2] + 500
                     )

threshGridList <- lapply(1:length(threshVal), function(i){threshGrid})

## Get quantiles for thresholds ##-------------------------------------
print("loop start")
ST <- Sys.time()
ST0 <- ST
for(n in 1:NH){
  print(n)
  # running time diagnosis
  if(n %% 50 == 0){
    I <- Sys.time() - ST
    I0 <- Sys.time() - ST0
    print(paste("Percent remaining", (NH-n)/NH))
    print(paste("Time remaining", (NH-n)/n * I0))
    print(paste("Since last readout:", I)); ST <- Sys.time() 
  }
  
  i <- rn[n,1]
  j <- rn[n,2]
  tSlice <- ncvar_get(ncin, varid="dmflow",
                       start=c(i, j,  1),
                       count=c(1, 1, -1))
  # find quantile for threshold
  thresh <- quantile(as.vector(tSlice), prob=c(1 - threshVal), na.rm=T) #vec
  threshMat[n,] <- thresh
  for(k in 1:NT){
	# add threshold k at position (i,j) to raster k
    threshGridList[[k]][i,j] <- thresh[k]
    # save which days cell n was exceeded.
    threshDayExcList[[k]][[n]] <- which(tSlice > thresh[k])
  }
}


## Save outputs ##-----------------------------------------------------
for(i in 1:NT){
  writeRaster(threshGridList[[i]],
              filename=paste0(wd, "InterimData/", threshName[i], 
                              "_threshGrid2.asc"),
              format="ascii",
              overwrite=TRUE)
}

saveRDS(threshDayExcList,
        file=paste0(wd,"InterimData/threshDayExcList2.rds"))
saveRDS(threshMat,
        file=paste0(wd,"InterimData/threshMat2.rds"))
readr::write_csv(x=data.frame(threshMat),
                 path=paste0(wd,"InterimData/threshMat2.csv"))

nc_close(ncin)