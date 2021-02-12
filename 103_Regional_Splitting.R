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
# OUTPUTS: hasData_Regions.csv
#
#~~~~~~~~~~~~~~~~~~~~~~~

library(rgeos)

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}


# if(file.exists(paste0(data_wd, "hasData_Regions.csv"))){
#     stop("hasData_Regions already exists for 103. Proceeding to next job.")
# }else{
#    print("Proceeding to job.")
# }

if(period == "future"){
  stop("Regions computed for present, not required for future (103). Proceeding to next job.")
}
#### DATA ------------------------------------------------------------
HA <- readOGR(dsn=paste0(data_wd, "hydrometricAreas"), layer="hyd_areas",
              stringsAsFactors=FALSE)
HA@data$HA_NUM <- as.numeric(HA@data$HA_NUM)
#plot(HA[80+1:4,], axes=T)

coordinates(rn) <- ~ east + nor
# Set the projection of the SpatialPointsDataFrame using the projection of
# the shapefile.
proj4string(rn) <- proj4string(HA)

pts <- over(rn, HA)

w <- which(is.na(pts[,1]))

nearw <- rep(0,NH)

pts$REGION <- rep(NA, NH)


#### PROCESSING -----------------------------------------------
for(i in w){  # find hydrometric area field
  
  tempDist <- gDistance(rn[i,], HA, byid=TRUE)
  nearw[i] <- which.min(tempDist)
  pts$HA_NUM[i] <-  HA$HA_NUM[which.min(tempDist)]
  try({
  pts$HA_NAME[i] <- HA$HA_NAME[which.min(tempDist)]
  })
  
}
HA$REGION <- 0
for(i in 1:length(HA$HA_NUM)){
  if(is.na(HA$HA_NUM[i])){ 
    next 
  }else if(HA$HA_NUM[i] >= 21 & HA$HA_NUM[i] < 28){
    HA$REGION[i] <- 1
  }else if(HA$HA_NUM[i] == 28){
    HA$REGION[i] <- 2
  }else if(HA$HA_NUM[i] > 28 & HA$HA_NUM[i] < 38){
    HA$REGION[i] <- 3
  }else if(HA$HA_NUM[i] > 39 & (HA$HA_NUM[i] < 43 | HA$HA_NUM[i] == 101)){
    HA$REGION[i] <- 4
  }else if(HA$HA_NUM[i] > 42 & HA$HA_NUM[i] < 54){
    HA$REGION[i] <- 5
  }else if(HA$HA_NUM[i] == 54){
    HA$REGION[i] <- 6
  }else if(HA$HA_NUM[i] > 54 & (HA$HA_NUM[i] < 68 | HA$HA_NUM[i] == 102)){
    HA$REGION[i] <- 7
  }else if(HA$HA_NUM[i] > 67 & (HA$HA_NUM[i] < 78 | HA$HA_NUM[i] == 103)){
    HA$REGION[i] <- 8
  }else if(HA$HA_NUM[i] > 37 & HA$HA_NUM[i] < 40){
    HA$REGION[i] <- 9
  }else if(HA$HA_NUM[i] %in% c(1:9, 90:97, 105:108)){
    HA$REGION[i] <- 10
  }else if(HA$HA_NUM[i] %in% 10:20){
    HA$REGION[i] <- 11
  }else if(HA$HA_NUM[i] %in% c(78:89, 104)){
    HA$REGION[i] <- 12
  }else{
    HA$REGION[i] <- NA
  }
}
# 
# vvv<- viridis(12)[sample.int(12,12)]
# 
# png("./regionmap.png", width=80, height=100, units="mm", res=300, pointsize=10)
# par(mar=c(1,0,1,1), mgp=c(2,1,0))
# plot(HA, col=vvv[HA$REGION], border=ifelse(is.na(HA$REGION),NA,1))
# #image.plot(z=HA$REGION, zlim=c(1,12), levels=12, col=vvv, legend.only=T)
# dev.off()

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
  }else if(pts$HA_NUM[i] %in% c(1:9, 90:97, 105:108)){
    pts$REGION[i] <- "NSC"
  }else if(pts$HA_NUM[i] %in% 10:20){
    pts$REGION[i] <- "ESC"
  }else if(pts$HA_NUM[i] %in% c(78:89, 104)){
    pts$REGION[i] <- "SSC"
  }else{
    pts$REGION[i] <- "UNK"
  }
}

rn2 <- cbind(rn@data,pts)

#### OUTPUT ----------------------------------------------------------
write_csv(rn2, path=paste0(wd_id, "hasData_Regions.csv"))
write_csv(rn2, path=paste0(data_wd, "hasData_Regions.csv"))
