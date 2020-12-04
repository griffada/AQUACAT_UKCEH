
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

partable <- data.frame(loc=numeric(), threshold=numeric(),
                       scale=numeric(), shape=numeric())

for(h in WWW){
  
  if(h %% 20 == 0){print(paste(h, "of", NH))}
  
  ei <- rn[h,1]
  ni <- rn[h,2]
  
  thr <- thresMat[h,jT]
  
  vals <- ncvar_get(ncin,
                    "dmflow",
                    start=c(ei, ni, 1),
                    count=c(1, 1, -1))
  o <- try({
  peak_vals <- vals[(extractPeaks(vecObs=vals, mintimeDiff=7) == 1) &
                      (vals > thr)]

  whattime <- which((extractPeaks(vecObs=vals, mintimeDiff=7) == 1) &
                      (vals > thr))
  meanint <- mean(whattime[-1] - whattime[-length(whattime)])/360

  #at_site_gpa <- pargpa(lmoms(peak_vals))

  at_site_gpa <- fevd(x=peak_vals, threshold = thr, type="GP")$results$par

  partable[h,] <- list(meanint, thr, at_site_gpa[1], at_site_gpa[2])

    1
  })
  if(inherits(o, "try-error")){
    peak_vals <- blockmaxxer(data.frame(vals),blocks=rep(1:30,each=360))
    
    at_site_gpa <- pargpa(lmoms(peak_vals))
    
    partable[h,] <- list(1, at_site_gpa$para[1],
                         at_site_gpa$para[2], at_site_gpa$para[3])
  }
  print(partable[h,])
  
  # We will assume the POTs to occur according to a Poisson process.
  # For events that occur according to a Poisson process, the return period is 
  # NOT equal to 1/(Annual Prob of Exc), because the events could happen multiple
  # times per year.
  
  # Instead we have RP = 2*(RP in half-years).
  # And we have AnnProbOfExc = 1/(1 - exp(-2*ProbOfExc)), where the Prob of Exc 
  # is derived from the GPA distribution.
}


##### EXPORT #####----------------------------- 
readr::write_csv(partable, paste0(data_wd,subfold, 
                                  "paramtable_",thresh1, "_RCM", RCM, suffix, ".csv"))
warnings()

