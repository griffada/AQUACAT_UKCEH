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
library(raster)

##### SETUP ----------------------------------------------------
ncname_lin <-"/prj/ResilRiskInds/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
ncname_win <- "V:/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"  
# This file is ~36GB on the linux server.

if(substr(osVersion,1,3)=="Win"){
  ncname <- ncname_win
}else{
  ncname <- ncname_lin
}



##### FUNCTIONS ------------------------------------------------

#exceedencePlot <- function(day, netcdf, var="dmflow", threshgrid){
threshgrid <- as.matrix(TAQ)  
for(k in 1:length(w)){
  day <- extremes$Day[w[k]]
  #flow <- raster(ncname_win, layer=day)
  fl <- ncvar_get(netcdf, var, start=c(1,1,day), count=c(-1, -1, 1))
  flowM <- as.matrix(fl)[,1250:1]
  
  #TAQ <- TAQ[,1250:1]
  
  
  s <- sum(flowM[!is.na(flowM)] > 0)
  
  Wdown <- which(flowM > threshgrid, arr.ind=T)
  print(dim(Wdown))
  WW <- which(flowM <= threshgrid, arr.ind=T)
  sw <- round(length(which(flowM > threshgrid))/s, 2) * 100 # % inundation
  
  for(i in 1:nrow(WW)){
    if(is.na(flowM[WW[i,1], WW[i,2]])) next
    if(!(flowM[WW[i,1], WW[i,2]] == -1)){
      flowM[WW[i,1], WW[i,2]] <- -0.5
    }
  }
  
  
  
  brks <- c(-1.01, -0.51, seq(0, max(1500, na.rm=T), length.out=20))
  
  png(paste0("./exampleevent_",extremes$Day[w[k]],"_",extremes$Inun[w[k]], ".png"),
      res=300, width=70, height=70, units='mm', pointsize=9)
  par(mar=c(3,3,1,1), mgp=c(2,1,0))
  image.plot(flowM, x=0:700, y=0:1250, ylim=c(0,1250), zlim=c(-2,1500),
             col=c("grey70", "grey70", heat.colors(19)), breaks=brks, asp=1,
             xlab="Easting", ylab="Northing")
  
  
  pu <- par('usr')
  
  text(pu[1] + (pu[2]-pu[1])*0.9, pu[3] + (pu[4] - pu[3])*0.9,
       paste0("POT5 threshold\n", sw, "% inundation"),
       adj=c(1,0.5))
  dev.off()
}
#  return(flowM)
#}


#### TESTING --------------------------------------------------------
if(FALSE){
  ncin <- nc_open(ncname, readunlim=FALSE)
  
  filepath <- paste0("/Scoping/InterimData/qPOT_", events_per_year,
                     "epy_sans_mask.csv")
  qPOT5 <- unname(as.matrix(read_csv(here::here(filepath),
                                     col_names=TRUE, na="NA")))

  w <- extremes$Day[250]
  TAQ <- t(as.matrix(qPOT5))
  w <- c(250, 903:908, 986, 803, 1975, 779, 1273)
  
  raster <- 
  
  for(i in 1:length(w)){
  png(paste0("./exampleevent_",extremes$Day[w[i]],"_",extremes$Inun[w[i]], ".png"),
      res=300, width=70, height=110, units='mm', pointsize=10)
    FM <- exceedencePlot(day=w[i], netcdf=ncin, var="dmflow", threshgrid=TAQ)
  dev.off()
  }
  nc_close(ncin)
}