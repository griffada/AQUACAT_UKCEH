#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin, 2020-02-21
#
# Testing of texmex slimline functions from 07c_texmex_slimline.R
# 
# This takes an unknown amount of time to run.
#
# For aquaCAT, Project 07441.
# 
# Created ABG 2020-02-21, edited 2020-05-04
#
# OUTPUTS: InterimData/slimline/step1_test1.rds,
#          InterimData/slimline/step2_test1.rds,
#          InterimData/slimline/step3_test1.rds,
#          InterimData/slimline/step4_test1.rds,
#          InterimData/slimline/step5_test1.rds,
#          InterimData/slimline/step6_test1.rds,
#          InterimData/slimline/step7_test1.rds
#
#~~~~~~~~~~~~~~~~~~~~~~~
#setwd("S:")
library(texmex)
library(dplyr)
library(extRemes)
library(reshape2)
library(tidyverse)

print(ST <- Sys.time())

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

source(paste0(wd, "07b_HTfunctions.R"))
source(paste0(wd, "07c_texmex_slimline.R"))


##### DATA #####-----------------------------------------------------------

# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
#threshDayExcList <- readRDS(paste0(wd_id, "threshDayExcList2.rds"))

# matrix of threshold value (col) at a given cell (row)
threshMat <- read_csv(paste0(wd_id, "threshMat2.csv"))
thresh0 <- unname(unlist(threshMat['X2']))
#dim(threshMat)  =  19914 x 5

# eventLList length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(wd_id, "eventLists03.RDa")) 
NE <- length(eventDayList[[2]][[4]]) # POT2, 2% inun.

# timewise maxima at each cell for each event ((NE + 2) x NH)
eventDF <- readr::read_csv(paste0(data_wd,"eventdf_POT2_pc05.csv"))

# PoE under different computations with extra data. Tidy format.
print("SETUP DONE.")
print(Sys.time() - ST)
ST <- Sys.time()

REG <- "NW"
r1 <- which(rn_regions$REGION == REG) # length = 1437
NREG <- length(r1)
event_region <- eventDF[r1,]
thresh_region <- threshMat[r1,]

present_region <- readr::read_csv(paste0(wd_id,"regionalEvents_exampleNW_2.csv"))
thresh0 <- unname(unlist(thresh_region['X2']))
#melt to reshape dataframe
present_melt <- dcast(present_region, eventNo~loc, value.var="val")
rownames(present_melt) <- paste0("E",present_melt$eventNo)
colnames(present_melt) <- paste0("L", 0:NREG)
present_melt <- present_melt[,-1]

D <- 200; j <- 1

### Need a "data" object
DATA <- present_melt[,((j-1)*D+1):(min(NREG, j*D))]
DATA0 <- present_melt
mqu <- 0.7
step1_test1 <- migpd_slim(mqu=mqu, penalty="none")

str(step1_test1, max.level=2)

#saveRDS(step1_test1, file=paste0(wd_id, "slimline/step1_D.rds"))
print("STEP 1 COMPLETE")
print(Sys.time() - ST)
ST <- Sys.time()

#step1_test1 <- readRDS(paste0(wd_id, "slimline/step1.rds"))
#MODELS <- step1_test1$models #

#saveRDS(MODELS, file=paste0(wd_id, "slimline/MODELS_D.rds"))

margins_temp <- list("laplace",
                p2q = function(p) ifelse(p <  0.5, log(2 * p), -log(2 * (1 - p))),
                q2p = function(q) ifelse(q <  0, exp(q)/2, 1 - 0.5 * exp(-q)))

# step2_test1 <- mexTransform_slim(marginfns=margins_temp, mth=step1_test1$mth,
#                                  method="mixture", r=NULL)
# 
# TRANSFORMED <- step2_test1$transformed
# 
# saveRDS(TRANSFORMED, file=paste0(wd_id, "slimline/TRANSFORMED_D.rds"))
# 
# str(step2_test1, max.level=2)
# 
# saveRDS(step2_test1, file=paste0(wd_id, "slimline/step2_D.rds"))
print("STEP 2 COMPLETE")
print(Sys.time() - ST)
ST <- Sys.time()
NREG <- D
COEFFS <- array(NA, dim=c(6,NREG-1,NREG))
Z <- array(NA, dim=c(24,NREG-1,NREG))

# k <- 1
# step3_test1 <- mexDependence_slim(dqu=0.7, mth=step1_test1$mth, which=k,
#                                   marginsTransformed=TRANSFORMED[,1:D])
# str(step3_test1, max.level=2)
# 
# k <- 2
# step3_test2 <- mexDependence_slim(dqu=0.7, mth=step1_test1$mth, which=k,
#                                   marginsTransformed=TRANSFORMED[,1:D])
# str(step3_test2, max.level=2)
# 
# saveRDS(step3_test1, file=paste0(wd_id, "slimline/step3_D.rds"))
# print("STEP 3 COMPLETE")
# print(Sys.time() - ST)
# ST <- Sys.time()


dqu <- 0.7

step1_test1 <- readRDS(paste0(wd_id, "slimline/step1_D.rds"))
MODELS <-      readRDS(paste0(wd_id, "slimline/MODELS_D.rds"))
TRANSFORMED <- readRDS(paste0(wd_id, "slimline/TRANSFORMED_D.rds"))
step2_test1 <- readRDS(paste0(wd_id, "slimline/step2_D.rds"))
step3_test1 <- readRDS(paste0(wd_id, "slimline/step3_D.rds"))


#library(doParallel)
# cl <- makeCluster(3, outfile=paste0("cl_out_D.txt"))
# registerDoParallel(cl)
# step4_test1 <- foreach(k=1:D) %dopar% {
#   
#   if((k < 10) | (k %% 50 == 0)){ print(paste0("Dependence at ", 100*k/NREG, "%")) }
#   
#   try(mexDependence_slim(dqu=dqu, mth=step1_test1$mth, which=k))
#   
# }
# 
# stopCluster(cl)
if(file.exists(paste0(wd_id, "slimline/step4_1_D.rds"))){
  step4_test1 <- readRDS(file=paste0(wd_id, "slimline/step4_1_D.rds"))
  COEFFS <- readRDS(file=paste0(wd_id, "slimline/COEFFS_D.rds"))
  Z <- readRDS(file=paste0(wd_id, "slimline/Z_D.rds"))
  print(step4_test1$errcode)
  print("STEP 4 COMPLETE")
}else{
  step4_test1 <- lapply(1:NREG,
                        function(k){print(paste0("Dependence at", 100*k/NREG, "%"))
                          try(mexDependence_slim(dqu=dqu, mth=step1_test1$mth, which=k))})
  saveRDS(step4_test1, file=paste0(wd_id, "slimline/step4_1_D.rds"))
  saveRDS(COEFFS, file=paste0(wd_id, "slimline/COEFFS_D.rds"))
  saveRDS(Z, file=paste0(wd_id, "slimline/Z_D.rds"))
  
  
  print("STEP 4 COMPLETE")
  print(Sys.time() - ST)
  ST <- Sys.time()
}
str(step4_test1, max.level=2)





# 
nSample <- 11  # still a very low number of events.
d <- NREG
mult <- 2
k <- 1

step5_test1 <- predict.mex_slim(which=k, referenceMargin=NULL, marginfns=margins_temp,
                                constrain=step3_test1$dependence$constrain,
                                coeffs_in=COEFFS[,,k], z_in=Z[,,k],
                                mth=step1_test1$mth, mqu=step1_test1$mqu, pqu = dqu,
                                nsim = nSample*d*mult,
                                d=d)
str(step5_test1, max.level=2)
# 
saveRDS(step5_test1, file=paste0(wd_id, "slimline/step5_D.rds"))
print("STEP 5 COMPLETE")
print(Sys.time() - ST)
ST <- Sys.time()

#step4_test1 <- readRDS(paste0(wd_id, "/slimline/step4.rds"))
#library(profvis)
#profvis({
step6_test1 <- mexMonteCarlo_slim(mexList=step4_test1,
                                  marginfns=step3_test1$dependence$marginfns,
                                  mth=step1_test1$mth,
                                  mqu=step1_test1$mqu,
                                  nSample=20, mult=2)
str(step6_test1, max.level=2)
#})
#step6_test1 <- readRDS(file=paste0(wd_id, "slimline/step6_D.rds"))

write_csv(step6_test1$MCsample, paste0(wd_id, "slimline/step7_MCSample.csv"), append=T)
saveRDS(step6_test1, file=paste0(wd_id, "slimline/step6_D.rds"))

print("STEP 6 COMPLETE.")
print(Sys.time() - ST)
ST <- Sys.time()
print("END")
print(Sys.time())