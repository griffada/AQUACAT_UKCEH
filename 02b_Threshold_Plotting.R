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
# OUTPUTS: /Plots/***_threshMap2.png
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

rn <- read.csv(paste0(wd,"CodeABG/InterimData/hasData.csv"))
NH <- nrow(rn) #19914

#matrix of pointwise threshold, matches rn rows.
threshMat <- read.csv(paste0(wd, "CodeABG/InterimData/threshMat2.csv"),
					  stringsAsFactors=FALSE)
					  
# threshDayExcList is a list of NT lists of NH vectors:
#     days when gridpoint (i,j) exceeded threshold T.
threshDayExcList <- readRDS(paste0(wd,"CodeABG/InterimData/threshDayExcList2.rds"))

# threshold rasters
threshGridList <- vector("list", NT)
brks <- c(-1.01,seq(0, 1800, length.out=21))
Lbrks <- c(-1, signif(exp(brks[-1]-1),1))

for(i in 1:NT){
  tG <- raster(paste0(
    wd,"CodeABG/InterimData/",threshName[i],"_threshGrid2.asc"))
  tG[tG>0] <- log(tG[tG>0]+1)
  focal1 <- suppressWarnings(focal(tG, w=matrix(1,3,3), fun=max, na.rm=T))
  png(paste0(wd,"/Plots/", threshName[i], "_threshMap2.png"),
      width=100, height=120, units='mm', res=200, pointsize=8)
  
    image.plot(
      as.matrix(focal1)[,1000:1], zlim=c(-2,10),
      x=0:700, y=0:1000, ylim=c(1000,0), xlim=c(0,700),
      xlab="Easting (thousands)", ylab="Northing (thousands)",
      asp=1, col=c("grey90", rev(heat.colors(10))), breaks=brks,
      lab.breaks=Lbrks
    )
				 
  dev.off()
  print(paste("plot", i, "complete."))
}

#nc_close(ncin)