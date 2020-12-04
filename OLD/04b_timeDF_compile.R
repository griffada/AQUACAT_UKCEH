######
# Adam Griffin, 2020-04-27
#
# Summarising size and time of extreme events. 
# Practiced on one 30-year period of data.
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

#cells on the river network
rn <- read.csv(paste0(wd,"/CodeABG/InterimData/hasData2.csv"),
				stringsAsFactors=FALSE)
NH <- nrow(rn)


# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS(
  paste0(wd,"/CodeABG/InterimData/threshDayExcList2.rds"))
  

# matrix of threshold value (col) at a given cell (row)
threshMat <- read.csv(paste0(wd,"/CodeABG/InterimData/threshMat2.csv"),
                     stringsAsFactors=FALSE)
dim(threshMat) #19914 x 5


#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(wd,"/CodeABG/InterimData/eventLists03.RDa"))


##### GET EVENT #####----------------------------------------------------

# Only using POT2 and 2% inundation minimums.
jV <- which(threshName=="POT2")
jI <- which(wsName == "pc05")

event2_2 <- eventLList[[jV]][[jI]]
days2_2 <- eventDayList[[jV]][[jI]]


# prealloc
eventDataFrame <- matrix(NA, ncol=length(event2_2), nrow=NH)

for(i in 1:length(event2_2)){
  
  if(i %% 25 == 0){print(i)}
  
  L <- event2_2[i]
  D <- days2_2[i]
  
  #space-slice
  vals <- ncvar_get(ncin,
                    "dmflow",
                    start=c(1, 1, days2_2[i]),
                    count=c(-1, -1, event2_2[i]))

  #maxima per cell for each event
  valsmax <- apply(vals, c(1, 2),
                   function(x){ifelse(all(is.na(x)),NA,max(x, na.rm=T))})
  
  vvec <- rep(NA,NH)
  for(k in 1:NH){
    # swap from array to vec
    vvec[k] <- valsmax[rn[k,1],rn[k,2]]
  }
  eventDataFrame[,i] <- vvec
  
}


##### OUTPUT #####----------------------------------------------------------

eventDataFrame <- cbind(rn, eventDataFrame)

readr::write_csv(eventDataFrame, path=paste0(wd,"/Data/eventdf_POT2_pc05.csv"))
