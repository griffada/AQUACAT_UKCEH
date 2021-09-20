######
# Adam Griffin, 2020-04-27
#
# Summarising size and time, and extracting daily PoE of extreme events.
# 
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-04-27
# Pipeline version ABG 2020-09-07
# NetCDF version 2021-07-06
#
#
# OUTPUTS: eventdf_***.nc: dataframe of each event summarised.
#
#####

##### SETUP #####------------------------------------------------------------

print("running 105")
if (interactive()) {commandArgs <- function(...){c("15","present")}}
#### SETUP ####----------------------
if (substr(osVersion,1,3) == "Win") {
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}
print(paste("Running for threshold", thresh1, "at", ws1, "minimum spread."))


if(settings$OBSflow){
 print("OBS_flow already exists. Finishing 105N.") 
}else{


##### DATA #####--------------------------------------------------------------
suffix_pres  <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- ncoriginal <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix, "_out.nc") 

ncin <- nc_open(ncpres) # This file is ~36GB on the linux server.

# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
threshDayExcList <- readRDS(
        paste0(data_wd, subfold, "threshDayExcList_RCM", RCM, suffix,".rds"))
  
# matrix of threshold value (col) at a given cell (row)
threshMat <- readRDS(
        paste0(data_wd, subfold, "threshMat_RCM", RCM, suffix, ".rds"))

#eventLList length of event L,
#eventDayList start of event L, NT lists of NW lists
load(paste0(data_wd, subfold, "eventLists_RCM", RCM, suffix, ".RDa"))

# initialSumm <- readdf(file=paste0(data_wd, subfold, "initialSummary_RCM",
#                                   RCM, suffix, ".csv"))
##### GET EVENTS #####-------------------------------------------------------
eventz  <- eventLList[[jT]][[jW]]
dayz    <- eventDayList[[jT]][[jW]]
maxInun <- maxInunList[[jT]][[jW]]

##### SET UP NETCEF #####----------------------------------------------------

savepath <- paste0(data_wd,subfold,
                  "eventOBS_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".nc")
cdfPrimer(RCM, period, method="OBS", NE=length(eventz), NH, thresh1,
          ws1, rn, savepath)
obs_events <- nc_open(savepath, write=TRUE)

##### GET EVENTS #####--------------------------------------------------------
# Get pointwise maxima across events.
fullInun <- c()
# Can be more parallel for timewise maxima
NEV <- length(eventz)
for (i in 1:length(eventz)) { 
  if ((i < 10) | (i %% 100 == 0)) { # time recording
    print(i)
    print(paste("Percent remaining", 100 * round((NEV - i)/NEV, 2)))
    ST <- Sys.time()
    print(ST) 
  }
  L <- eventz[i]
  D <- dayz[i]
  
  #space-slice
  vals <- ncvar_get(ncin,
                    "dmflow",
                    start=c(1, 1, D),
                    count=c(-1, -1, L),
                    collapse_degen=FALSE)
  vvec <- rep(NA, NH)
  #vmat <- matrix(NA, NH, L)
  for(k in 1:NH){
    vvec[k] <- max(vals[rn$row[k], rn$col[k],], na.rm=T)
    #vmat[k,] <- vals[rn$row[k], rn$col[k],]
  }

  ncvar_put(obs_events, "flow",
            vals=signif(vvec,6),
            start=c(1,i), count=c(-1,1))

  fullInun[i] <- sum(vvec > threshMat[,jT])
}
ncvar_put(obs_events, "eventNo", 1:length(eventz))

saveRDS(fullInun,
        file=paste0(data_wd,subfold, "eventInun_OBS_",thresh1,"_", ws1, 
                    "_RCM", RCM, suffix, ".RDS"))

nc_close(obs_events)
nc_close(ncin)

settings$OBSflow <- TRUE
write_yaml(settings, settingspath)
}
print(Sys.time())
print("finished 105N.")