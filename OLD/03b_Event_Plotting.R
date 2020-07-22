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


# netCDF test file
ncnameLin <-"/prj/ResilRiskInds/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
ncnameWin <- "V:/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"  
# This file is ~36GB on the linux server.
if(substr(osVersion,1,3)=="Win"){
  ncname <- ncnameWin
  wd <- "S:"
}else{
  ncname <- ncnameLin
  wd <- "/prj/aquacat"
}

rn <- read.csv(paste0(wd,"/CodeABG/InterimData/hasData.csv"))
NH <- nrow(rn) #19914

#matrix of pointwise threshold, matches rn.
threshMat <- readRDS(paste0(wd,"/CodeABG/InterimData/threshMat.rds"))
# threshDayExcList is a list of NT lists of NH vectors:
#     days when gridpoint (i,j) exceeded threshold T.
threshDayExcList <- readRDS(paste0(wd,"/CodeABG/InterimData/threshDayExcList.rds"))

# print(ST <- Sys.time())
# ncin <- nc_open(ncname)
# print(ncin)
# print(Sys.time() - ST)

# threshold rasters
threshGridList <- vector("list", NT)
#brks <- c(-1.01,seq(0, 800, length.out=21))

for(i in 1:NT){
  threshGridList[[i]] <- raster(paste0(
                    wd,"/CodeABG/InterimData/",threshName[i],"_threshGrid.asc"))
  png(paste0(wd,"/Plots/", threshName[i], "_threshMap.png"),
      width=70, height=150, units='mm', res=200, pointsize=8)
	  
      plot(threshGridList[[i]], 
                 #x=0:700, y=0:1250, ylim=c(1250,0), xlim=c(0,700),
                 asp=1, col=c("grey90", heat.colors(20))#, breaks=brks
				 )
				 
  dev.off()
  print(paste("plot", i, "complete."))
}

#nc_close(ncin)