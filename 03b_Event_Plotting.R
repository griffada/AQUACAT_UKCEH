######
# Adam Griffin, 2020-07-22
#
# Event plotting from daily flow. Practiced on one 30-year period of data.
# 
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-07-22
#
#####

##### SETUP #####------------------------------------------------------------

library(fields)
library(ncdf4)
library(raster)
library(knitr)
library(kableExtra)
library(dplyr)


if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

#### FUNCTIONS ####--------------------------------------

eventPlotter <- function(flow, threshold, filename=NULL, loc_only=FALSE){
  
  NH <- sum(flow[!is.na(flow)] > 0) # number of river network cells
  
  inunLocs <- which(flow < threshold, arr.ind=T)
  iL2 <- which(flow > threshold, arr.ind=T)
  inunLevel <- round((NH - nrow(inunLocs))/NH, 2) * 100 # % inundation
  
  flow2 <- flow
  flow2[inunLocs] <- -0.5
  if (loc_only) flow2[iL2] <- 1
  flow2 <- as.matrix(focal(raster(flow2), w=matrix(1,3,3), 
                           fun=function(x){ifelse(all(is.na(x)), NA, max(x, na.rm=T))}))
  
  MF <- ceiling(max(flow2[iL2]))
  if (loc_only) {
    colpal <- c("grey90", "grey90", "red", "blue")
    brks <- c(-1.01, -0.51, -0.01, 0.99, 1.99)
  } 
  else{
    colpal <- c("grey90", "grey90", heat.colors(19))
    brks <- c(-1.01, -0.51, seq(0, MF, length.out=20))
  }
  
  if (!is.null(filename)) {
    png(filename, res=300, width=70, height=70, units='mm', pointsize=9)
    par(mar=c(3,3,1,1), mgp=c(2,1,0))
    image.plot(flow2[1000:1,], x=0:700, y=0:1000, ylim=c(0,1000),
               zlim=c(-2,MF),
               col=colpal,
               breaks=brks, asp=1,
               xlab="Easting", ylab="Northing")
    
    pu <- par('usr')
    text(pu[1] + (pu[2]-pu[1])*0.9, pu[3] + (pu[4] - pu[3])*0.9,
         paste0("POT5 threshold\n", inunLevel, "% inundation"),
         adj=c(1,0.5))
    dev.off()
  }
  par(mar=c(3, 3, 1, 1), mgp=c(2, 1, 0))
  image.plot(flow2[1000:1,], x=0:700, y=0:1000, ylim=c(0,1000),
             zlim=c(-2, MF),
             col=colpal,
             breaks=brks, asp=1,
             xlab="Easting", ylab="Northing")
  
  pu <- par('usr')
  text(pu[1] + (pu[2] - pu[1])*0.9, pu[3] + (pu[4] - pu[3])*0.9,
       paste0("POT5 threshold\n", inunLevel, "% inundation"),
       adj=c(1,0.5))
  
  return(flow2)
}


eventPlotterNoMap <- function(flow_tab, threshold_tab, base_map,
                              filename=NULL, loc_only=FALSE){
  
  NH <- nrow(flow_tab) # number of river network cells
  if(nrow(flow_tab) != nrow(threshold_tab)){
    stop("flow cells should match threshold cells")
  }
  
  rn <- flow_tab[,1:2] # Eastings and Northings
  
  base_map <- as.matrix(base_map)
  if (nrow(base_map) != 1000) base_map <- t(base_map)
  
  flow_tab[flow_tab < threshold_tab] <- -0.5
  
  
  
  
  if (loc_only){
    flow_tab[flow_tab > threshold_tab] <- 1
    colpal <- c("grey90", "grey90", "red", "blue")
    brks <- c(-1.01, -0.51, -0.01, 0.99, 1.99)
    MF <- 1
  }else{
    MF <- ceiling(max(flow_tab))
    colpal <- c("grey90", "grey90", heat.colors(19))
    brks <- c(-1.01, -0.51, seq(0, MF, length.out=20))
  }
  
  for(n in 1:NH){
    base_map[rn[n,1], rn[n,2]] <- flow_tab[n,3]
  }  
  
  flow2 <- as.matrix(focal(raster(base_map), w=matrix(1,3,3), 
                           fun=function(x){ifelse(all(is.na(x)), NA, max(x, na.rm=T))}))
  
  if (!is.null(filename)) {
    png(filename, res=300, width=70, height=70, units='mm', pointsize=9)
    par(mar=c(3,3,1,1), mgp=c(2,1,0))
    image.plot(flow2, x=0:700, y=0:1000, ylim=c(0,1000),
               zlim=c(-2,MF),
               col=colpal,
               breaks=brks, asp=1,
               xlab="Easting", ylab="Northing")
    
    pu <- par('usr')
    text(pu[1] + (pu[2]-pu[1])*0.9, pu[3] + (pu[4] - pu[3])*0.9,
         paste0("POT5 threshold\n", inunLevel, "% inundation"),
         adj=c(1,0.5))
    dev.off()
  }
  par(mar=c(3, 3, 1, 1), mgp=c(2, 1, 0))
  image.plot(flow2, x=0:700, y=0:1000, ylim=c(0,1000),
             zlim=c(-2, MF),
             col=colpal,
             breaks=brks, asp=1,
             xlab="Easting", ylab="Northing")
  
  pu <- par('usr')
  text(pu[1] + (pu[2] - pu[1])*0.9, pu[3] + (pu[4] - pu[3])*0.9,
       paste0("POT5 threshold\n", inunLevel, "% inundation"),
       adj=c(1,0.5))
  
  return(flow2)
}


### TESTING ### --------------------------
if(FALSE){
  
  threshGrid <- raster(x=paste0(wd_id,
                                threshName[2], "_threshGrid2.asc"), layer=1)
  netcdf <- nc_open(ncname)
  V <- values(threshGrid)
  tg <- t(matrix(V, nrow=1000, ncol=700))
  var <- "dmflow"
  fl <- ncvar_get(netcdf, var, start=c(1,1,5369), count=c(-1, -1, 1))
  
  #debug(eventPlotter)
  ep <- eventPlotter(flow=fl, threshold=tg, filename=NULL, loc_only=FALSE)
  
  image(ep[,1000:1])
  
}