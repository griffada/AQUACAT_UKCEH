######
# Adam Griffin, 2020-04-27
#
# Summarising size and time of extreme events. Practiced on one 30-year period of data.
# 
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-04-27
#
#####

##### SETUP #####------------------------------------------------------------

library(ncdf4)
library(raster)
library(fields)


if(substr(osVersion,1,3)=="Win"){
  wd <- "S:"
}else{
  wd <- "/prj/aquacat"
}
ncname <- paste0(wd, "/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc")

# threshold for inundation
threshVal <- c(5/365, 2/365, 1/365, 0.2/365, 0.1/365)
threshName <- c("POT5", "POT2", "POT1", "Q5", "Q10")
NT <- length(threshVal)

# cut-off bound for widespread event
wsBound <- c(0.05, 0.02, 0.01, 0.005, 0.001)
wsName <- c("pc5", "pc2", "pc1", "pc05", "pc01")
NW <- length(wsBound)


##### DATA #####------------------------------------------------------------

print(ST <- Sys.time())
	ncin <- nc_open(ncname) # This file is ~36GB on the linux server.
	print(ncin)
print(Sys.time() - ST)

ND <- 10800 # Number of days

#cells on the ruver network
rn <- read.csv(paste0(wd,"/CodeABG/InterimData/hasData2.csv"),
				stringsAsFactors=FALSE)
NH <- nrow(rn)


threshDayExcList <- readRDS(paste0(wd,"/CodeABG/InterimData/threshDayExcList2.rds"))
  # 5 lists of NH lists
thresMat <- read.csv(paste0(wd,"/CodeABG/InterimData/threshMat2.csv"),
                     stringsAsFactors=FALSE)

load(paste0(wd,"/CodeABG/InterimData/eventLists03.RDa"))



##### GET EVENT #####
jV <- which(threshName=="POT2")
jI <- which(wsName == "pc2")

event2_2 <- eventLList$POT2$pc2
days2_2 <- eventDayList$POT2$pc2

eventDataFrame <- matrix(NA, nrow=length(event2_2), ncol=NH)

for(i in 1:length(event2_2)){
  
  if(i %% 25 == 0){print(i)}
  L <- event2_2[i]
  D <- days2_2[i]
  vals <- ncvar_get(ncin, "dmflow", start=c(1,1,D), count=c(-1,-1,L))

  valsmax <- apply(vals, c(1, 2),
                   function(x){ifelse(all(is.na(x)),NA,max(x, na.rm=T))})
  
  vvec <- rep(NA,NH)
  for(k in 1:NH){
    vvec[k] <- valsmax[rn[k,1],rn[k,2]]
  }
  eventDataFrame[i,] <- vvec
  
}

eventDataFrame <- cbind(rn, eventDataFrame)

readr::write_csv(eventDataFrame, file=paste0(wd,"/CodeABG/InterimData/eventdf04.csv"))




