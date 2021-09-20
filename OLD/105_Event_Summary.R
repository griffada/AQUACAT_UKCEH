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

print("running 105")
if(interactive()){commandArgs <- function(...){c("01","present")}}
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}



# thresh1 <- "POT2"
# ws1 <- "pc05"
# # Only using POT2 and 2% inundation minimums.
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)
print(paste("Running for threshold", thresh1, "at", ws1, "minimum spread."))
# 
# if(file.exists(paste0(data_wd,subfold,"eventflow_OBS_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))){
#   stop("eventflow_OBS_ exists. finishing 105.")
# }


##### DATA #####------------------------------------------------------------
suffix_pres <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- ncoriginal <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix, "_out.nc") 

print(ST <- Sys.time())
	ncin <- nc_open(ncoriginal) # This file is ~36GB on the linux server.
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

# paramtable <- as.data.frame(readr::read_csv(paste0(data_wd,subfold_pres, 
#                             "paramtable_",thresh1, "_RCM", RCM, suffix_pres, ".csv"),
#                             col_types=cols(.default= col_double())))
# 
# colnames(partable)[1] <- "meanint"


##### GET EVENT #####----------------------------------------------------

eventz <- eventLList[[jT]][[jW]]
dayz <- eventDayList[[jT]][[jW]]

maxInun <- apply(cbind(eventz, dayz), 1, function(x){max(dvec[x[2]+(0:(x[1]-1))])})

# prealloc
eventDataFrame <- matrix(NA, ncol=length(eventz), nrow=NH)
  ST0 <- Sys.time()
  ST <- Sys.time()
  
fullInun <- c()
# Can be more parallel for timewise maxima
#for(i in seq_len(length(eventz))){
for(i in 1:length(eventz)){ 
  if((i < 10) | (i %% 100 == 0)){ # time recording
    print(i)
    I <- difftime(Sys.time(), ST, units="secs")
    I0 <- difftime(Sys.time(), ST0, units="secs")
    print(paste("Percent remaining", 100*round((length(eventz)-i)/length(eventz) ,2)))
    print(paste("Time remaining", round((length(eventz)-i)/i * I0,2)))
    print(paste("Since last readout:", round(I,2)))
    ST <- Sys.time()
    print(ST) 
  }
  L <- eventz[i]
  D <- dayz[i]
  
  #space-slice
  vals <- ncvar_get(ncin,
                    "dmflow",
                    start=c(1, 1, D),
                    count=c(-1, -1, L), collapse_degen=FALSE)

  #maxima per cell for each event
  # valsmax <- apply(vals, c(1, 2),
  #                  function(x){ifelse(all(is.na(x)),NA,max(x, na.rm=T))})
  vvec <- rep(NA, NH)

  for(k in 1:NH){
    # swap from array to vec
    vvec[k] <- max(vals[rn$row[k], rn$col[k],], na.rm=T)
  }

  eventDataFrame[,i] <- vvec
  fullInun[i] <- sum(vvec > threshMat[,jT])
}

##### SAVE OUTPUTS #####----------------------------------------------------------
eventDataFrame <- cbind(rn, round(as.data.frame(eventDataFrame),8))
colnames(eventDataFrame) <- c(colnames(rn), paste0("E",1:(ncol(eventDataFrame)-4)))
readr::write_csv(eventDataFrame, path=paste0(data_wd,subfold,
                  "eventflow_OBS_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
saveRDS(fullInun, file=paste0(data_wd,subfold,
                  "eventInun_OBS_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".RDS"))

nc_close(ncin)