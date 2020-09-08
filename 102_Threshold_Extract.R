#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-02-21
#
# Event threshold extraction from daily flow.
#
# This takes about 9 hours to run on linux if you're lucky.
# This should only need running for one G2G output if everything goes as 
# expected.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-02-21, edited 2020-05-04
# Pipeline version 2020-09-07
#
# OUTPUTS: threshDayExcList.rda,
#          thresMat.rda,
#          threshGrid2.asc
#
#~~~~~~~~~~~~~~~~~~~~~~~

#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

threshGridList <- list()

threshDayExcList <- vector("list", length(threshVal))

#### PARSE COMMAND LINE ARGS --------------------------------------------------
# if(length(args)==0){
#   RCM <- "01"
#   period <- "present"
#   suffix <- "_198012_201011"
# }else if(length(args)==1){
#   RCM <- sprintf("%02d",args[1])
#   period <- "present"
#   suffix <- "_198012_201011"
# }else if(length(args)==2){
#   RCM <- sprintf("%02d", args[1])
#   period <- args[2]
#   if(period=="present"){
#     suffix <- "_198012_201011"
#   }else if (period=="future"){
#     suffix <- "_205012_201011"
#   }else{  
#     stop("correct call: 102_Threshold_Extract.R RCM [period]. Period should be 'present' or 'future'.")
#   }
# }
# if(as.numeric(RCM) < 0 | as.numeric(RCM) > 16){
#   stop("correct call: 102_Threshold_Extract.R RCM [period]. RCM should be between 1 and 16.")
# }
# 
# ncname <- paste0(wd, "run_hmfg2g/outputs/dmflow_RCM", RCM, suffix, "_out.nc") 
# 
# nccopy <- paste0(wd, "Data/dmflow_copy_RCM", RCM, suffix, ".nc")


ST <-  Sys.time()
ncin <- nc_open(nccopy, readunlim=FALSE) 
# this is a huge file, do not open without reason.
print(Sys.time() - ST)
print(ncin)


#### THRESHOLD EXTRACT ####------------------------------------
## Get grid of river network ##-------------------
tStart <- 1
deltaT <- 1
ncwide <- ncvar_get(ncin, "dmflow",
					start=c(1,1,tStart),
					count=c(-1, -1, 1))

rn <- which(apply(ncwide, c(1,2),
                 function(v){sum(v[!is.na(v)] > -1) == deltaT}), arr.ind=T)

rn <- data.frame(row=rn[,1], col=rn[,2])
east <- ncvar_get(ncin, "Easting")
north <- ncvar_get(ncin, "Northing")
rn$east <- east[rn[,1]]
rn$nor <- north[rn[,2]]

readr::write_csv(data.frame(rn), paste0(data_wd, "hasData_primary.csv"))

NH <- nrow(rn)
## Initialise dfs and lists ##-------------------------

threshGrid <- ncvar_get(ncin, "dmflow", start=c(1,1,tStart),
                        count=c(-1, -1, 1))
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
print(ST)
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


saveRDS(threshGridList, file=paste0(data_wd, subfold,
                          "threshGridList_RCM", RCM, suffix,".rds"))

## Save outputs ##-----------------------------------------------------
# for(i in 1:NT){
#   writeRaster(threshGridList[[i]],
#               filename=paste0(wd, "InterimData/", threshName[i], 
#                               "_threshGrid2.asc"),
#               format="ascii",
#               overwrite=TRUE)
# }
# 
saveRDS(threshDayExcList, file=paste0(data_wd, subfold,
                          "threshDayExcList_RCM", RCM, suffix,".rds"))
saveRDS(threshMat, file=paste0(data_wd, subfold,
                          "threshMat_RCM", RCM, suffix,".rds"))
readr::write_csv(x=data.frame(threshMat), path=paste0(data_wd, subfold,
                          "threshMat_RCM", RCM, suffix,".csv"))

nc_close(ncin)