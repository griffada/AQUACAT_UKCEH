#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-07-31
#
# Splitting GB rover network points by Hydrological Areas.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-07-31
# Pipeline version 2020-09-07
#
#~~~~~~~~~~~~~~~~~~~~~~~

library(rgeos)

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

#### DATA ------------------------------------------------------------
HA <- readOGR(dsn=paste0(wd_id, "hydrometricAreas"), layer="hyd_areas",
              stringsAsFactors=FALSE)
HA@data$HA_NUM <- as.numeric(HA@data$HA_NUM)
#plot(HA[80+1:4,], axes=T)

coordinates(rn) <- ~ east + nor
# Set the projection of the SpatialPointsDataFrame using the projection of
# the shapefile.
proj4string(rn) <- proj4string(HA)

pts <- over(rn, HA)

w <- which(is.na(pts[,1]))

nearw <- rep(0,NR)

pts$REGION <- rep(NA, NR)


#### PROCESSING -----------------------------------------------
for(i in w){  # find hydrometric area field
  
  tempDist <- gDistance(rn[i,], HA, byid=TRUE)
  nearw[i] <- which.min(tempDist)
  pts$HA_NUM[i] <-  HA$HA_NUM[which.min(tempDist)]
  try({
  pts$HA_NAME[i] <- HA$HA_NAME[which.min(tempDist)]
  })
  
}

for(i in 1:length(pts$HA_NUM)){
  if(is.na(pts$HA_NUM[i])){ 
    next 
  }else if(pts$HA_NUM[i] >= 21 & pts$HA_NUM[i] < 28){
    pts$REGION[i] <- "NE"
  }else if(pts$HA_NUM[i] == 28){
    pts$REGION[i] <- "TRE"
  }else if(pts$HA_NUM[i] > 28 & pts$HA_NUM[i] < 38){
    pts$REGION[i] <- "ANG"
  }else if(pts$HA_NUM[i] > 39 & (pts$HA_NUM[i] < 43 | pts$HA_NUM[i] == 101)){
    pts$REGION[i] <- "SE"
  }else if(pts$HA_NUM[i] > 42 & pts$HA_NUM[i] < 54){
    pts$REGION[i] <- "SW"
  }else if(pts$HA_NUM[i] == 54){
    pts$REGION[i] <- "SEV"
  }else if(pts$HA_NUM[i] > 54 & (pts$HA_NUM[i] < 68 | pts$HA_NUM[i] == 102)){
    pts$REGION[i] <- "WAL"
  }else if(pts$HA_NUM[i] > 67 & (pts$HA_NUM[i] < 78 | pts$HA_NUM[i] == 103)){
    pts$REGION[i] <- "NW"
  }else if(pts$HA_NUM[i] > 37 & pts$HA_NUM[i] < 40){
    pts$REGION[i] <- "THA"
  }else{
    pts$REGION[i] <- "SCO"
  }
}

rn2 <- cbind(rn@data,pts)

#### OUTPUT ----------------------------------------------------------
write_csv(rn2, path=paste0(wd_id, "hasData_Regions.csv"))
write_csv(rn2, path=paste0(data_wd, "hasData_Regions.csv"))
