
#### SETUP ####----------------------------------------------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

library(ilaprosUtils)
library(lmomco)
library(extRemes)
library(pastecs)

ws1 <- "pc05"
thresh1 <- "POT2"
jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)
print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))


##### DATA #####-------------------------------------------------------

#5 lists of NH lists, one for each threshold

# Step 1, get the GPA estimates of return period
# 1a) Get all values above threshold.
# 1b) Extract peaks (this may have to be reversed)

# 1c) Fit GPA using L-moments.
# 1d) Compute prob of exceedence (1-CDF) in 6 months.
# 1e) Compute annual prob of exceedence (2 periods).

# Here we use the at-site exceedences, not the widespread events.

# Observed Threshold from 102, Threshold Extract.

thresMat <- readRDS(paste0(data_wd, subfold, "threshMat_RCM",
                           RCM, suffix, ".rds"))

print(ST <- Sys.time())
ncin <- nc_open(ncname)
print(ncin)
print(Sys.time() - ST)

partable <- read_csv(paste0(data_wd,subfold, 
                      "paramtable_",thresh1, "_RCM", RCM, suffix, ".csv"))


eventDF <- readr::read_csv(paste0(data_wd,subfold, "eventdf_",thresh1,"_", ws1,
                                  "_RCM", RCM, suffix, ".csv"),
                           col_types=cols(
                             .default = col_double()
                           ))

eventDF <- matrix(unlist(eventDF), nrow=nrow(eventDF), ncol=ncol(eventDF))

rarityDF <- eventDF

NE <- ncol(eventDF)-4

for(h in 1:nrow(eventDF)){
  if((h %% 200) == 0){
    print(paste(h, "of", nrow(eventDF)))
  }
  
  thr <- thresMat[h,jT]
  
  rps <- as.numeric(pevd(eventDF[h,-(1:4)], scale=partable$scale[h],
                         shape=partable$shape[h],
              threshold = thr, type='GP'))
  
  rps[rps < 1e-10] <- NA
  
  rarityDF[h,-(1:4)] <- 1 - rps
}

rarityDF <- data.frame(rarityDF)
colnames(rarityDF) <- c("row", "col","east","nor", paste0("E", 1:(ncol(eventDF)-4)))

write_csv(data.frame(rarityDF), paste0(data_wd,subfold, "OBSraritydf_",thresh1,"_", ws1,
                 "_RCM", RCM, suffix, ".csv"))
