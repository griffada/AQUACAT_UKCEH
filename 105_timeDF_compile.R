######
# Adam Griffin, 2020-04-27
#
# Summarising size and time of extreme events.
# 
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-04-27
# Pipeline version ABG 2020-09-07
#
#
# Outputs:
#
# eventdf_***.csv: dataframe of each event summarised. One event per row.
#
#####

##### SETUP #####------------------------------------------------------------

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

thresh1 <- "POT2" #!#!#!#!# Important constants to select.
ws1 <- "pc05"
print("Running for threshold", POT2, "at ", ws1, "minimum spread.")

##### DATA #####------------------------------------------------------------

print(ST <- Sys.time())
	ncin <- nc_open(ncname) # This file is ~36GB on the linux server.
	print(ncin)
print(Sys.time() - ST)

# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS( paste0(data_wd, subfold, "threshDayExcList_RCM",
                                    RCM, suffix,".rds"))
  

# matrix of threshold value (col) at a given cell (row)
threshMat <- read.csv(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"),
                     stringsAsFactors=FALSE)
#dim(threshMat) #19914 x 5


#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))


##### GET EVENT #####----------------------------------------------------

# Only using POT2 and 2% inundation minimums.
jV <- which(threshName==thresh1)
jI <- which(wsName == ws1)

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

readr::write_csv(eventDataFrame, path=paste0(data_wd,subfold,
                  "eventdf_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
