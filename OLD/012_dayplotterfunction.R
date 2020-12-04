#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-08-11
#
# Plotting example events in terms of return period.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-08-11
#
#~~~~~~~~~~~~~~~~~~~~~~~
if(substr(osVersion,1,3) == "Win"){
  ncname <- "S:/CodeABG/InterimData/dmflow_timechunks.nc"
  wd <- "S:/"
}else{
  ncname <- "/prj/aquacat/CodeABG/InterimData/dmflow_timechunks.nc"
  wd <- "/prj/aquacat/"
}



present <- readr::read_csv(paste0(wd,"Data/present_returnlevels_POT2_pc05.csv"))

present0 <- present %>%
  mutate(flood = val > thresh) %>%
  group_by(eventNo, DayS) %>%
  summarise(size = sum(flood)/19914) %>%
  ungroup() %>%
  arrange(desc(size))

netcdf <- nc_open(ncname)

rn <- readr::read_csv(paste0(wd,"CodeABG/InterimData/hasData_regions.csv"))

dayplotter <- function(present, netcdf, event=NULL, day=NULL, plotpath){
  
  if(is.null(event) & is.null(day)){
    stop("Supply day or event or both.")
  }
  
  if(!is.null(event) & is.null(day)){
    day <- present$DayS[(present$eventNo==event)[1]]
  }else if(is.null(event) & !is.null(day)){
    event <- present$eventNo[(present$DayS==day)[1]]
  }
  
  present_day <- present %>% dplyr::filter(eventNo == event)
  
  fl <- ncvar_get(netcdf, var, start=c(1,1,day), count=c(-1, -1, 1))
  
  flowM <- as.matrix(fl)
  
  s <- sum(flowM[!is.na(flowM)] > 0)
  
  Wdown <- present_day[which(present_day$val > present_day$thresh), 
                       c("loc","Easting","Northing")]
  print(dim(Wdown))
  WW <- present_day[which(present_day$val <= present_day$thresh), 
                    c("loc","Easting","Northing")]
  sw <- round(nrow(Wdown)/s, 2)*100
  
  for (i in 1:nrow(WW)) {
    if(is.na(flowM[WW$Easting[i], WW$Northing[i]])) next
    if(!(flowM[WW$Easting[i], WW$Northing[i]] == -1)){
      flowM[WW$Easting[i], WW$Northing[i]] <- -0.5
    }
  }
  for (j in 1:nrow(rn)) {
    flowM[rn[j,1],rn[j,2]] <- -0.5
  }
  for (i in 1:nrow(Wdown)) {
    flowM[WW$Easting[i], WW$Northing[i]] <- 
      1 / present_day$gpp[present_day$loc == WW$loc[i]]
  }
  
  flowN <- focal(raster(flowM[, 1000:1]),
                 w=matrix(1, 3, 3),
                 fun=function(x){ifelse(!all(is.na(x)),max(x, na.rm=T),NA)})
  
  brks <- c(-1.01, -0.51, seq(0, max(200, na.rm=T), length.out=20))
  
  png(plotpath,
      res=300, width=100, height=100, units='mm', pointsize=10)
  par(mar=c(1,1,1,1), mgp=c(1,1,0))
  image.plot(as.matrix(flowN), x=0:700, y=0:1000, ylim=c(0,1000), zlim=c(-2,200),
             col=c("darkseagreen1", "darkseagreen1", topo.colors(19)),
             breaks=brks, asp=1,
             xlab="", ylab="", axes=F)
  dev.off()
  
  par(mar=c(1,1,1,1), mgp=c(1,1,0))
  image.plot(as.matrix(flowN), x=0:700, y=0:1000, ylim=c(0,1000), zlim=c(-2,200),
             col=c("darkseagreen1", "darkseagreen1", topo.colors(19)),
             breaks=brks, asp=1,
             xlab="", ylab="", axes=F)
  
  return(flowN)
}

for(n in 1:5){
  event <- present0$eventNo[n]
  day <- present0$DayS[n]
  
  dayplotter(present, netcdf, event, day,
             plotpath = paste0("S:/Plots/example_event",event,"_day",day,".png")
  
}