######
# Adam Griffin, 2020-03-10
#
# Event plotting from daily flow. Practiced on one 30-year period of data.
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-03-10
#
#####

library(here)
library(ncdf4)
library(fields)
library(readr)

##### SETUP ----------------------------------------------------
ncname_lin <-"/prj/ResilRiskInds/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
ncname_win <- "K:/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"  
# This file is ~36GB on the linux server.

if(substr(osVersion,1,3)=="Win"){
  ncname <- ncname_win
}else{
  ncname <- ncname_lin
}



##### FUNCTIONS ------------------------------------------------

exceedencePlot <- function(day, netcdf, var="dmflow", threshgrid){
  flow <- raster(ncname_win, layer=day)
  #ncvar_get(netcdf, var, start=c(1,1,day), count=c(-1, -1, 1))
  #flowM <- as.matrix(flow)
  threshgrid <- raster(threshgrid)
  
  s <- sum(flow[!is.na(flow)] > 0)
  
  WW <- which(flow <= threshgrid, arr.ind=T)
  sw <- round(length(which(flowM > threshgrid))/s, 2) * 100 # % inundation
  
  for(i in 1:nrow(WW)){
    if(is.na(flowM[WW[i,1], WW[i,2]])) next
    if(!(flowM[WW[i,1], WW[i,2]] == -1)){
      flowM[WW[i,1], WW[i,2]] <- -0.5
    }
  }
  
  
  
  brks <- c(-1.01, -0.51, seq(0, max(flowM, na.rm=T), length.out=20))
  
  image.plot(flowM, x=0:700, y=0:1250, ylim=c(1250,0),
             col=c("grey90", "grey70", heat.colors(19)), breaks=brks, asp=1,
             xlab="Easting", ylab="Northing")
  
  pu <- par('usr')
  
  text(pu[1] + (pu[2]-pu[1])*0.9, pu[3] + (pu[4] - pu[3])*0.9,
       paste0("POT5 threshold\n", sw, "% inundation"),
       adj=c(1,0.5))
  return(flowM)
}


#### TESTING --------------------------------------------------------
if(FALSE){
  ncin <- nc_open(ncname, readunlim=FALSE)
  
  filepath <- paste0("/Scoping/InterimData/qPOT_", events_per_year,
                     "epy_sans_mask.csv")
  qPOT5 <- unname(as.matrix(read_csv(here::here(filepath),
                                     col_names=TRUE, na="NA")))

  w <- extremes$Day[250]
  TAQ <- t(as.matrix(qPOT5))
  FM <- exceedencePlot(day=w, netcdf=ncin, var="dmflow", threshgrid=TAQ)
  
  nc_close(ncin)
}