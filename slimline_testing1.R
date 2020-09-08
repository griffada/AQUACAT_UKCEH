
#setwd("S:")
library(texmex)
library(readr)
library(dplyr)
library(ncdf4)
library(extRemes)
library(reshape2)
library(tidyverse)

print(ST <- Sys.time())

if (substr(osVersion,1,3) == "Win") {
  ncname <- "S:/CodeABG/InterimData/dmflow_timechunks.nc"  # rechunked for spaceslices
  ncname2 <- "S:/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc" 
  # pre-rechunk
  wd <- "S:/"
  wd_id <- "S:/CodeABG/InterimData/"
  wd_cd <- "S:/CodeABG/"
  
} else {
  ncname <- "/prj/aquacat/CodeABG/InterimData/dmflow_timechunks.nc"
  ncname2 <- "S:/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc"
  wd <- "/prj/aquacat/"
  wd_id <- "/prj/aquacat/CodeABG/InterimData/"
  wd_cd <- "/prj/aquacat/CodeABG/"
}

source(paste0(wd_cd,"07b_HTfunctions.R"))

##### DATA #####-----------------------------------------------------------

ND <- 10800 # Number of days

# river network
rn <- read_csv(paste0(wd_id, "hasData2.csv"))

rn_regions <- read_csv(paste0(wd_id, "hasData_Regions.csv"))
NH <- nrow(rn)


# threshold for inundation
threshVal <- c(5/365, 2/365, 1/365, 0.2/365, 0.1/365)
threshName <- c("POT5", "POT2", "POT1", "Q5", "Q10")
NT <- length(threshVal)



# cut-off bound for widespread event
wsBound <- c(0.05, 0.02, 0.01, 0.005, 0.001)
wsName <- c("pc5", "pc2", "pc1", "pc05", "pc01")
NW <- length(wsBound)



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
eventDF <- readr::read_csv(paste0(wd,"Data/eventdf_POT2_pc05.csv"))


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
source(paste0(wd_cd,"07c_texmex_slimline.R"))
thresh0 <- unname(unlist(thresh_region['X2']))
#melt to reshape dataframe
present_melt <- dcast(present_region, eventNo~loc, value.var="val")
rownames(present_melt) <- paste0("E",present_melt$eventNo)
colnames(present_melt) <- paste0("L", 0:NREG)
present_melt <- present_melt[,-1]

D <- 60; j <- 1

### Need a "data" object
data0 <- present_melt[,((j-1)*D+1):(min(NREG, j*D))]
data <- present_melt
mqu <- 0.7
# step1_test1 <- migpd_slim(mqu=mqu, penalty="none")
# 
# str(step1_test1, max.level=2)
# 
# saveRDS(step1_test1, file=paste0(wd_id, "/slimline/step1.rds"))
# print("STEP 1 COMPLETE")
# print(Sys.time() - ST)
# ST <- Sys.time()

step1_test1 <- readRDS(paste0(wd_id, "/slimline/step1.rds"))
models <- step1_test1$models #

# margins_temp <- list("laplace",
#                 p2q = function(p) ifelse(p <  0.5, log(2 * p), -log(2 * (1 - p))),
#                 q2p = function(q) ifelse(q <  0, exp(q)/2, 1 - 0.5 * exp(-q)))

# step2_test1 <- mexTransform_slim(marginfns=margins_temp, mth=step1_test1$mth,
#                                  method="mixture", r=NULL)
# 
# str(step2_test1, max.level=2)
# 
# saveRDS(step2_test1, file=paste0(wd_id, "/slimline/step2.rds"))
# print("STEP 2 COMPLETE")
# print(Sys.time() - ST)
# ST <- Sys.time()

k <- 1
step3_test1 <- mexDependence_slim(dqu=0.7, mth=step1_test1$mth, which=k)
str(step3_test1, max.level=2)

#saveRDS(step3_test1, file=paste0(wd_id, "/slimline/step3.rds"))
print("STEP 3 COMPLETE")
print(Sys.time() - ST)
ST <- Sys.time()


dqu <- 0.7
# step4_test1 <- lapply(1:NREG, function(k){mexDependence_slim(dqu=dqu, mth=step1_test1$mth, which=k)})
# 
# saveRDS(step4_test1, file=paste0(wd_id, "/slimline/step4.rds"))
# print("STEP 4 COMPLETE")
# print(Sys.time() - ST)
# ST <- Sys.time()
# 
nSample <- 50
d <- NREG
mult <- 10
# step5_test1 <- predict.mex_slim(object= step3_test1,
#                                 mth=step1_test1$mth, mqu=step1_test1$mqu, pqu = dqu,
#                                 nsim = nSample * d * mult)
# str(step5_test1, max.level=2)
# 
# saveRDS(step5_test1, file=paste0(wd_id, "/slimline/step5.rds"))
# print("STEP 5 COMPLETE")
# print(Sys.time() - ST)
# ST <- Sys.time()

step4_test1 <- readRDS(paste0(wd_id, "/slimline/step4.rds"))

step6_test1 <- mexMonteCarlo_slim(marginfns=step3_test1$dependence$marginfns,
                                  mth=step1_test1$mth,
                                  mqu=step1_test1$mqu,
                                  nSample=nSample, mexList=step4_test1, mult=mult)
str(step6_test1, max.level=2)

saveRDS(step6_test1, file=paste0(wd_id, "/slimline/step6.rds"))
print("STEP 6 COMPLETE.")
print(Sys.time() - ST)
ST <- Sys.time()
print("END")
print(Sys.time())