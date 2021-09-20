rcms <- c("01","04","05","06","07","08","09","10","11","12","13","15")
REG <- "NW"

commandArgs <- function(...){c("06","present","NW")}

#### SETUP ####----------------------------------------------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

for (RCM in rcms){
  for(period in c("future")){
    suffix <- ifelse(period=="present","_198012_201011","_205012_208011")
    subfold <- paste0("RCM", RCM, suffix, "/")
    
    savepath <- paste0(data_wd,subfold, REG, "/eventHT_",
                        REG,"_",thresh1,"_", ws1,"_RCM", RCM, suffix, ".nc")
    
    eventflow <- readdf(paste0(data_wd,subfold, REG, "/eventflow_HT_",
                        REG,"_",thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"))
    
    r1 <- which(rn_regions$REGION == REG)
    NH1 <- length(r1)
    
    cdfPrimer(RCM=RCM, period=period, method="HT", NE=ncol(eventflow)-4,
              NH=NH1, thresh1, ws1, rn[r1,], savepath, chunks=F)
    ht_events <- nc_open(savepath, write=T)
     
    ncvar_put(ht_events, "flow", as.matrix(eventflow[,-(1:4)]))
    
    cn <- sapply(colnames(eventflow)[-(1:4)], function(x){strsplit(x,"[E\\.]")[[1]][2]})
    ncvar_put(ht_events, "eventNo", as.numeric(cn))
    
    nc_close(ht_events)
  }
}
