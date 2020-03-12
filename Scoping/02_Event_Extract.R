######
# Adam Griffin, 2020-02-13
#
# Event extraction from daily flow. Practiced on one 30-year period of data.
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-02-13
#
#####

##### SETUP ##### --------------------------------------------------------

#### Packages -------------------------------------------------------
library(ncdf4)
library(fields)
library(raster)


#### Data -----------------------------------------------------------


## Read in netCDF ------

ncname_lin <-"/prj/ResilRiskInds/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
ncname_win <- "K:/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"  
# This file is ~36GB on the linux server.

if(substr(osVersion,1,3)=="Win"){
  ncname <- ncname_win
}else{
  ncname <- ncname_lin
}
ncin <- nc_open(ncname, readunlim=FALSE)
print("ncin open")

# Subset for testing
istart <- 1
Di <- 50 
jstart <- 1
Dj <- 50
tstart <- 1

## Restrict to river network ---------

# `NA` is sea pixels, `-1` is not on the river network.
Dt <- 1
ncwide <- ncvar_get(ncin, "dmflow",
                    start=c(1,1,tstart), count=c(-1, -1, Dt))
print("ncwide get")

HasData <- which(apply(ncwide, c(1,2),
                       function(v){sum(v[!is.na(v)] > -1) == Dt}),
                 arr.ind=T)
NH <- nrow(HasData)
# HasData is important for gluing the events back in at the other end 
# of the pipeline.


### EXTRA IMPORTANT SHAPES -------------------------------------------------
HD <- raster(x=ncname_win, layer=1)  # gets extent etc. automatically
HD[HD < 0] <- NA  # sets "land" to "ocean"
HD[HD >= 0] <- 1 # sets "rivers" to 1

writeRaster(HD, "./Scoping/InterimData/RiverNetwork.asc", format='ascii')

HX <- raster(x=ncname_win, layer=1)  # gets extent etc. automatically
HX[HX >= -1] <- -1  # sets rivers to "land"
plot(HX)
writeRaster(HX, "./Scoping/InterimData/GBshape.asc", format="ascii")


## Event extraction ------------
qpot5_grid <- raster(x=ncname_win, layer=1)

qpot5_grid[!is.na(qpot5_grid) & qpot5_grid > -1] <- 0

qpot_extract <- function(n){
  i <- HasData[n,1]
  j <- HasData[n,2]
  tslice <- ncvar_get(ncin, varid="dmflow",
                      start=c(i, j,  1),
                      count=c(1, 1, 10957))
  qpot5 <- quantile(as.vector(tslice), prob=c(1 - (5/365)))
  qpot5_grid[i,j] <- qpot5
}

qpot_extract_full <- function(ncfile, varid, quant, step=50){
  
  qpot_grid <- raster(x=ncfile, layer=1)
  di <- seq(1, dim(qpot_grid)[1], by=step)
  dj <- seq(1, dim(qpot_grid)[2], by=step)
  pb <- txtProgressBar(0, length(di)*length(dj))
  k <- 0
  for(i in di){
    for(j in dj){
      tslice <- ncvar_get(ncin, varid="dmflow",
                          start=c(i, j,  1),
                          count=c(step, step, -1))
      qpot <- apply(X=tslice, c(1,2), FUN=quantile, probs=quant, na.rm=T,
                    names=FALSE)
      qpot_grid[i+(0:(step-1)), j+(0:(step-1))] <- qpot
      k <- k+1
      setTxtProgressBar(pb, k)
    }
  }
  close(pb)
  return(qpot_grid)
  
}

QP99 <- qpot_extract_full(ncname_win, "dmflow", 0.99, step=50)



for(n in 1:NH){
  qpot_extract(n)
  if(n %% 1 == 0){
    print("Expected time to end = ")
    print(((NH-n)/n)*difftime(Sys.time(),ST.old, unit="hours"))
  }
}

save(qpot5_grid, file="./InterimData/qpot5_grid.RDa")