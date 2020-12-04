######
# Adam Griffin, 2020-04-27
#
# Summarising size and time, and extracting daily PoE of extreme events.
# 
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-04-27
# Pipeline version ABG 2020-09-07
#
#
# OUTPUTS: eventdf_***.csv: dataframe of each event summarised. One event per row.
#
#####

##### SETUP #####------------------------------------------------------------

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}



thresh1 <- "POT2"
ws1 <- "pc05"
# Only using POT2 and 2% inundation minimums.
jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)
print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

if(file.exists(paste0(data_wd,subfold,
                      "eventdf_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))){
  stop("eventdf_ exists. finishing 105.")
}


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
threshMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                            RCM, suffix, ".rds"))
#dim(threshMat) #19914 x 5


#eventLList (length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))


paramtable <- readr::read_csv(paste0(data_wd,subfold, 
                                     "paramtable_",thresh1, "_RCM", RCM, suffix, ".csv"))


##### GET EVENT #####----------------------------------------------------

eventz <- eventLList[[jT]][[jW]]
dayz <- eventDayList[[jT]][[jW]]

# prealloc
eventDataFrame <- matrix(NA, ncol=length(eventz), nrow=NH)


# Can be more parallel for timewise maxima
for(i in seq_len(length(eventz))){
  
  if(i %% 25 == 0){print(i)}
  
  L <- eventz[i]
  D <- dayz[i]
  
  #space-slice
  vals <- ncvar_get(ncin,
                    "dmflow",
                    start=c(1, 1, D),
                    count=c(-1, -1, L))
  
  vals_event <- vals[D+seq_len(L)-1]

  #maxima per cell for each event
  valsmax <- apply(vals, c(1, 2),
                   function(x){ifelse(all(is.na(x)),NA,max(x, na.rm=T))})
  
  vvec <- rep(NA, NH)
  for(k in 1:NH){
    # swap from array to vec
    vvec[k] <- valsmax[rn$row[k], rn$col[k]]
  }
  eventDataFrame[,i] <- vvec

}

##### SAVE OUTPUTS #####----------------------------------------------------------

# eventDataFrame <- cbind(rn, eventDataFrame)

readr::write_csv(as.data.frame(eventDataFrame), path=paste0(data_wd,subfold,
                  "eventflow_OBS_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))

nc_close(ncin)