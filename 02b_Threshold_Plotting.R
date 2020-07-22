#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-05-04
#
# Event extraction from daily flow. Practiced on one 30-year period of data.
# Rewritten for linux to reduce read-write times for netCDF.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-05-04
#
#~~~~~~~~~~~~~~~~~~~~~~~

print(Sys.time())
##### SETUP #####---------------------------------------------
library(ncdf4)
library(raster)
library(fields)

# threshold for inundation
threshVal <- c(5/365, 2/365, 1/365, 0.2/365, 0.1/365)
threshName <- c("POT5", "POT2", "POT1", "Q5", "Q10")
NT <- length(threshVal)

# cut-off for widespread event
wsCutoff <- c(0.05, 0.02, 0.01, 0.005, 0.001)
NW <- length(wsCutoff)


##### DATA #####--------------------------------------------

if(substr(osVersion,1,3) == "Win"){
  wd <- "S:/"
}else{
  wd <- "/prj/aquacat/"
}

ncname <- "run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc" #2GB
# ST <-  Sys.time()
# ncin <- nc_open(paste0(data_wd,ncname))  # this is a huge file, do not open without reason.
# print(Sys.time() - ST)
# print(ncin)

# # This file is ~36GB on the linux server.
# if(substr(osVersion,1,3)=="Win"){
#   ncname <- ncnameWin
#   wd <- "S:/"
# }else{
#   ncname <- ncnameLin
#   wd <- "/prj/aquacat/"
# }

ncname <- paste0()

rn <- read.csv(paste0(wd,"CodeABG/InterimData/hasData.csv"))
NH <- nrow(rn) #19914

#matrix of pointwise threshold, matches rn rows.
threshMat <- read.csv(paste0(wd, "CodeABG/InterimData/threshMat2.csv"),
					  stringsAsFactors=FALSE)
					  
# threshDayExcList is a list of NT lists of NH vectors:
#     days when gridpoint (i,j) exceeded threshold T.
threshDayExcList <- readRDS(paste0(wd,"CodeABG/InterimData/threshDayExcList2.rds"))

# print(ST <- Sys.time())
# ncin <- nc_open(ncname)
# print(ncin)
# print(Sys.time() - ST)

# threshold rasters
threshGridList <- vector("list", NT)
brks <- c(-1.01,seq(0, 1800, length.out=21))

for(i in 1:NT){
  threshGridList[[i]] <- raster(paste0(
                    wd,"CodeABG/InterimData/",threshName[i],"_threshGrid2.asc"))
  focal1 <- focal(threshGridList[[i]], w=matrix(1,3,3), fun=max, na.rm=T)
  png(paste0(wd,"/Plots/", threshName[i], "_threshMap2.png"),
      width=100, height=120, units='mm', res=200, pointsize=8)
	  
      plot(#t(threshGridList[[i]]),
           t(focal1),zlim=c(-2,1800),
                 #x=0:700, y=0:1250, ylim=c(1250,0), xlim=c(0,700),
                 asp=1, col=c("grey90", heat.colors(20)), breaks=brks
				 )
				 
  dev.off()
  print(paste("plot", i, "complete."))
}

#nc_close(ncin)