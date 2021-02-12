#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-10-16
#
# Function for plotting extracted values in terms of return period.
#
# For aquaCAT, Project 07441.
# 
# OUTPUTS: .png plots
#
#~~~~~~~~~~~~~~~~~~~~~~~
library(ncdf4)
library(raster)
library(fields)
library(readr)
library(dplyr)

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

gbgrid <- readRDS("S:/data/Misc/gbshape.RDS")

genericMapPlotter <- function(rn, vals, threshold=NULL, 
                              label="Value", gbgrid,
                              filename=NULL){
  # rn          data frame of river network: row, col, E, N. 
  #               (should be loaded by setup_script_00.R)
  # vals        vector of values to plot at corresponding points of rn.
  # gbgrid      matrix of values for GB shape (-1 for land, NA for sea)
  # threshold   vector of values for thresholds corresponding to points of rn.
  #               if supplied, vals below threshold are greyed out.
  # filename    if supplied, .png is saved to this filepath.
  #
  # label       label for legend describing value plotted.
  
  NH <- nrow(rn)
  if(length(vals) != NH){
    stop("Vals not correct length.")
  }
  if(is.null(threshold)){
    threshold <- rep(-1e6, NH)
  }else{
    if(length(threshold) != NH){
      stop("Threshold not correct length.")
    }
  }
  vals[vals < threshold] <- -0.5
  vals[is.na(vals)] <- -1
  vmax <- max(vals, na.rm=T)
  
  M <- gbgrid
  for(i in 1:NH){
    M[rn$row[i], rn$col[i]]<-vals[i]
  }
  
  Mfocal <- raster::focal(raster(M[,ncol(M):1]),
                          w=matrix(1, 3, 3),
                          fun=function(x){
                            if(!all(is.na(x))){max(x, na.rm=T)}else{return(NA)}
                          })
  brks <- c(-1.01, -.51, seq(0, max(vals, na.rm=T), length.out=20))
  
  f <- function(){
    par(mar=c(1,1,1,1), mgp=c(1,1,0))
    image.plot(as.matrix(Mfocal),
               x=0:700, y=0:1000, ylim=c(0,1000), zlim=c(-2,max(M)),
               col=c("darkseagreen1", "darkseagreen2", topo.colors(19)),
               breaks=brks, asp=1,
               xlab="", ylab="", axes=F)
  }
  
  f()
  if(!is.null(filename)){
    png(filename, res=300, width=100, height=100, units='mm', pointsize=10)
    f()
    dev.off()
  }
}



returnLevelMapPlotter <- function(RL_df, rn_matrix, rn, eventNumber,
                                  poe_days=TRUE,
                                  filename=NULL){
  # plots a map of return periods (in terms of years) and saves it as .png
  #
  # RL_df       event dataframe: should include Easting, Northing,
  #                            eventNo, gpp, val, thresh.
  # rn_matrix   matrix showing location of river network in GB. 
  #               (can be taken straight from ncdf4::ncvar_get)
  # rn          data frame of river network: row, col, E, N. 
  #               (should be loaded by setup_script_00.R)
  # eventNumber number from RL_df of desired event.
  # poe_days    if TRUE, converts from POE in days to POE in years.
  # filename    if provided, a .png is saved to this filepath.
  
  
  RL_day <- RL_df %>% dplyr::filter(eventNo == eventNumber)
  
  day <- RL_day$DayS[1] # day number, you should infer from the event number.
  
  s <- sum(rn_matrix[!is.na(rn_matrix)]>0)
  
  w_above_thresh <- RL_day$val > RL_day$thresh
  
  Wdown <- RL_day[which(w_above_thresh),  ]
  WW    <- RL_day[which(!w_above_thresh), ]
  sw    <- round(sum(w_above_thresh)/s, 2) * 100

  rn_matrix[rn[,1:2]] <- -0.5
  
  if(poe_days){
    rn_matrix[as.matrix(Wdown[,c("Easting","Northing")])] <- 
        1 / (1 - (1-Wdown$gpp)^360)
  }
  flowR <- raster(rn_matrix[,rev(seq_len(ncol(flowR)))])
  
  flowN <- focal(flowR, w=matrix(1, 3, 3),
                 fun=function(x){
                    if(!all(is.na(x))){max(x, na.rm=T)}else{return(NA)}
                 })
  
  brks <- c(-1.01, -0.51, seq(0, max(70, na.rm=T), length.out=20))
  
  f <- function(){
    par(mar=c(1,1,1,1), mgp=c(1,1,0))
    image.plot(as.matrix(flowN),
               x=0:700, y=0:1000, ylim=c(0,1000), zlim=c(-2,200),
               col=c("darkseagreen1", "darkseagreen2", topo.colors(19)),
               breaks=brks, asp=1,
               xlab="", ylab="", axes=F)
    pu <- par('usr')
    text(pu[1] + (pu[2] - pu[1])*0.9, pu[3] + (pu[4] - pu[3])*0.9,
         paste0("POT2 threshold\n", sw, "% inundation"),
         adj=c(1,0.5), cex=0.8)
  }
  
  f()
  
  if(!is.null(filename)){
    png(filename, res=300, width=100, height=100, units='mm', pointsize=10)
      f()
    dev.off()
  }
}

#### DEBUGGING ####
if(FALSE){
  present <- readr::read_csv(paste0(data_wd,
                                    "TestData/present_returnlevels_POT2_pc05.csv"),
                             col_types = cols( .default = col_double()))
  
  netcdf <- nc_open(paste0(wd_id,"dmflow_timechunks.nc"))
  fl <- ncvar_get(netcdf, "dmflow", start=c(1,1,1), count=c(-1, -1, 1))
  eventNo <- 17
  
  returnLevelMapPlotter(RL_df=present, rn_matrix=fl, rn,
                        eventNumber=eventNo, filename="test0.png")
}