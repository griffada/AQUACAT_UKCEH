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
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

#source(paste0(wd, "07b_HTfunctions.R"))
source(paste0(wd, "07c_texmex_slimline.R"))


##### DATA #####-----------------------------------------------------------
jI <- which(threshName=="POT2")
jV <- which(wsName=="pc05")
# lists of which days different thresholds were exceeded at different points
# NT lists of NW lists
#threshDayExcList <- readRDS(paste0(wd_id, "threshDayExcList2.rds"))

# matrix of threshold value (col) at a given cell (row)
threshMat <- readr::read_csv(paste0(wd_id, "threshMat2.csv"), col_types=cols(
  X1 = col_double(),
  X2 = col_double(),
  X3 = col_double(),
  X4 = col_double(),
  X5 = col_double()
))
thresh0 <- unname(unlist(threshMat['X2']))
#dim(threshMat)  =  19914 x 5

# eventLList length of event L, NT lists (by threshold) of NW lists 
# (by inun cutoff))
# eventDayList start of event L, NT lists of NW lists
load(paste0(wd_id, "eventLists03.RDa")) 
NE <- length(eventDayList[[jI]][[jV]]) # POT2, 2% inun.

# timewise maxima at each cell for each event ((NE + 2) x NH)
eventDF <- readr::read_csv(paste0(data_wd,"TestData/eventdf_POT2_pc05.csv"),
                           col_types=paste0(rep("c",431),collapse="")
                           )

# PoE under different computations with extra data. Tidy format.


REG <- "NW"
r1 <- which(rn_regions$REGION == REG) # length = 1437
NREG <- length(r1)
event_region <- eventDF[r1,]
thresh_region <- threshMat[r1,]

present_region <- readr::read_csv(paste0(wd_id,"regionalEvents_exampleNW_2.csv"),
                                  col_types=cols(
                                    eventNo = col_double(),
                                    loc = col_double(),
                                    Easting = col_double(),
                                    Northing = col_double(),
                                    thresh = col_double(),
                                    DayS = col_double(),
                                    val = col_double(),
                                    gpp = col_double(),
                                    ecdf = col_double(),
                                    gev = col_double()
                                  ))
thresh0 <- unname(unlist(thresh_region['X2']))
#melt to reshape dataframe
present_melt <- dcast(present_region, eventNo~loc, value.var="val")
rownames(present_melt) <- paste0("E",present_melt$eventNo)
colnames(present_melt) <- paste0("L", 0:NREG)
present_melt <- present_melt[,-1]

D <- 200; j <- 1

### Need a "data" object
DATA0 <- present_melt[,((j-1)*D+1):(min(NREG, j*D))]
DATA <- present_melt

print("SETUP DONE.")
print(Sys.time() - ST)
ST <- Sys.time()

mqu <- 0.7

dqu <- 0.7
if(file.exists(paste0(wd_id, "slimline/MODELS.rda"))){
  load(paste0(wd_id, "slimline/MODELS.rda"))
}else{
  step1_test1 <- migpd_slim(mqu=mqu, penalty="none")
  str(step1_test1, max.level=2)
  MODELS <- step1_test1$models
  mth <- step1_test1$mth
  mqu <- step1_test1$mth
  save(MODELS, mth, file=paste0(wd_id, "slimline/MODELS.rda"))
}
mqu <- rep(mqu, length(mth))[1:length(mth)]
#step1_test1 <- readRDS(paste0(wd_id, "slimline/step1.rds"))


print("STEP 1 COMPLETE")
print(Sys.time() - ST)
ST <- Sys.time()


#MODELS<- readRDS(paste0(wd_id, "slimline/MODELS.rds"))
#TRANSFORMED <- readRDS(file=paste0(wd_id, "slimline/TRANSFORMED.rds"))

margins_temp <- list("laplace",
                     p2q = function(p) ifelse(p <  0.5, log(2 * p), -log(2 * (1 - p))),
                     q2p = function(q) ifelse(q <  0, exp(q)/2, 1 - 0.5 * exp(-q)))

if(file.exists(paste0(wd_id, "slimline/TRANSFORMED.rds"))){
  TRANSFORMED <- readRDS(paste0(wd_id, "slimline/TRANSFORMED.rds"))
}else{
  step2_test1 <- mexTransform_slim(marginfns=margins_temp, mth=mth,
                                   method="mixture", r=NULL)
  
  TRANSFORMED <- step2_test1$transformed
  
  saveRDS(TRANSFORMED, file=paste0(wd_id, "slimline/TRANSFORMED.rds"))
  
  str(step2_test1, max.level=2)
  
  saveRDS(step2_test1, file=paste0(wd_id, "slimline/step2.rds"))
}

print("STEP 2 COMPLETE")
print(Sys.time() - ST)
ST <- Sys.time()
#NREG <- D

if(file.exists(paste0(wd_id, "slimline/step3.rds"))){
  step3_test1 <- readRDS(paste0(wd_id, "slimline/step3.rds"))
  
}else{
  COEFFS <- array(NA, dim=c(6,NREG-1,NREG))
  Z <- array(NA, dim=c(24, NREG-1, NREG))
  
  k <- 1
  step3_test1 <- mexDependence_slim(dqu=0.7, mth=mth, which=k,
                                    marginsTransformed=TRANSFORMED)
  str(step3_test1, max.level=2)
  saveRDS(step3_test1, file=paste0(wd_id, "slimline/step3.rds"))
}


# COEFFS <- array(NA, dim=c(6,NREG-1,NREG))
# Z <- array(NA, dim=c(24,NREG-1,NREG))
# 
# k <- 2
# step3_test2 <- mexDependence_slim(dqu=0.7, mth=step1_test1$mth, which=k, 
#                                   marginsTransformed=TRANSFORMED, zspot=TRUE)
# str(step3_test2, max.level=2)
# 

COEFFS <- array(NA, dim=c(6,NREG-1,NREG))
Z <- array(NA, dim=c(step3_test1$zspot[1], NREG-1, NREG))

#step3_test1 <- readRDS(file=paste0(wd_id, "slimline/step3.rds"))
if(file.exists(paste0(wd_id, "slimline/step4.rds"))){
  step4_test1 <- readRDS(file=paste0(wd_id, "slimline/step4.rds"))
  COEFFS <- readRDS(file=paste0(wd_id, "slimline/COEFFS.rds"))
  Z <- readRDS(file=paste0(wd_id, "slimline/Z.rds"))
  print(step4_test1$errcode)
  print("STEP 4 COMPLETE")
  print(Sys.time() - ST)
  ST <- Sys.time()
}else{
step4_test1 <- lapply(1:NREG,
                function(k){
                  #print(k)
                  if (k%%25==0) print(paste0("Dependence at", round(100*k/NREG,2), "%"))
                  o <- try(mexDependence_slim(dqu=dqu, mth=mth, which=k,
                                                    marginsTransformed=TRANSFORMED))

                  if(inherits(o, "try-error")){
                    print(k)
                    if(interactive()) browser()
                    warning(paste0("Error in mexDependence loop, iteration", k))
                  }
                  else if(o$errcode > 0){
                    print("Error code ", o$errcode, " on iter ", k, ".")
                  }
                  if(k < 5){print(str(o))}
                  o
              })

print("STEP 4 COMPLETE")
print(Sys.time() - ST)
ST <- Sys.time()
#str(step4_test1, max.level=1)

# 
saveRDS(step4_test1, file=paste0(wd_id, "slimline/step4.rds"))

saveRDS(COEFFS, file=paste0(wd_id, "slimline/COEFFS.rds"))

saveRDS(Z, file=paste0(wd_id, "slimline/Z.rds"))
}


#Z <- readRDS(paste0(wd_id, "slimline/Z.rds"))

# print("STEP 4 COMPLETE")
# print(Sys.time() - ST)
# ST <- Sys.time()
# 
# # 
# step4_test1 <- readRDS(paste0(wd_id, "slimline/step4.rds"))
# 
# COEFFS <- readRDS(paste0(wd_id, "slimline/COEFFS.rds"))
# 
# Z <- readRDS(paste0(wd_id, "slimline/Z.rds"))

nSample <- 20  # still a very low number of events.
d <- NREG
mult <- 2
k <- 1

margins_temp <- list("laplace",
                     p2q = function(p) ifelse(p <  0.5, log(2 * p), -log(2 * (1 - p))),
                     q2p = function(q) ifelse(q <  0, exp(q)/2, 1 - 0.5 * exp(-q)))

step5_test1 <- predict.mex_slim(which=k, referenceMargin=NULL,
                                marginfns=margins_temp,
                                constrain=step3_test1$dependence$constrain,
                                coeffs_in=COEFFS[,,k], z_in=Z[,,k],
                                mth=mth, mqu=mqu,
                                pqu = dqu,
                                nsim = nSample * d * mult,
                                printout = TRUE,
                                d = d)
str(step5_test1, max.level=2)
# 
saveRDS(step5_test1, file=paste0(wd_id, "slimline/step5.rds"))
print("STEP 5 COMPLETE")
print(Sys.time() - ST)
ST <- Sys.time()

#step4_test1 <- readRDS(paste0(wd_id, "/slimline/step4.rds"))

step6_test1 <- mexMonteCarlo_slim(mexList=step4_test1,
                                  marginfns=margins_temp,
                                  mth=mth,
                                  mqu=mqu,
                                  nSample=nSample, mult=mult)
str(step6_test1, max.level=1)
#step6_test1 <- readRDS(file=paste0(wd_id, "slimline/step6.rds"))

write.csv(step6_test1$MCsample, paste0(wd_id, "slimline/step7_MCSample2.csv"))
saveRDS(step6_test1, file=paste0(wd_id, "slimline/step6.rds"))

# print(ST <- Sys.time())
# ncin <- nc_open(ncname) # This file is ~2.5GB on the linux server.
# print(ncin)
# print(floor(Sys.time() - ST))
# 
# rn_here <- cbind(rn,rn_regions[,3:5]) %>%
#   dplyr::filter(REGION=="NW") %>%
#   dplyr::select(row,col,east,nor)
# 
# NE <- nrow(step6$MCsample)
# NH <- ncol(step6$MCsample)
# 
# rarityDF_HT <- expand.grid("eventNo" = 1:NE,
#                         "loc" = 1:NH)
# 
# rarityDF_HT$Easting <- rn_here[rarityDF_HT[, 2], 3]
# rarityDF_HT$Northing <- rn_here[rarityDF_HT[, 2], 4]
# rarityDF_HT$thresh <- NA
# rarityDF_HT$val <- NA
# rarityDF_HT$gpp <- NA
# rarityDF_HT$ecdf <- NA
# rarityDF_HT$gev <- NA
# 
# 
# logit <- function(x){log(x/(1-x))}
# invlogit <- function(y){1/(1 + exp(-1*y))}
# 
# gringorten <- function(v){
#   ((length(v) + 1 - rank(v)) - 0.44)/(length(v) + 0.12)
# }
# 
# weibull <- function(v){
#   (length(v) + 1 - rank(v))/(length(v) + 1)
# }
# 
# 
# ST0 <- proc.time()
# ST <- proc.time()
# print("loop start")
# for(n in 1:NH){ # for eacl location N
#   if(n < 10 |(n %% 200 == 0)){
#     print(paste(n, "out of", NH))
#   }
#   
#   i <- rn_here[n,1]
#   j <- rn_here[n,2]
#   # Pull out spaceslice
#   tSlice <- ncvar_get(ncin, varid="dmflow",
#                       start=c(i, j,  1),
#                       count=c(1, 1, -1))
#   
#   tSliceEvent <- unname(unlist(step6$MCsample[n, -(1:4)]))
#   
#   W <- which(rn_regions$REGION=="NW")[1:NH]
#   
#   thresh1 <- thresh0[W]
#   
#   threshval <- thresh1[n] # POT2 column
#   
#   WU <- which(rarityDF_HT$loc == W[n])
#   
#   rarityDF_HT$thresh[WU] <- threshval
#   rarityDF_HT$val[WU] <- tSliceEvent
#   
#   # get ecdf and estimate PoE
#   
#   ecdfSlice <- ecdf(tSlice)
#   poeEvent <- 1 - ecdfSlice(tSliceEvent)
#   
#   
#   # Weibull or Gringorten plotting position
#   
#   grSlice <- gringorten(tSlice)
#   grEvent <- sapply(tSliceEvent,
#                     function(x){grSlice[which.min(abs(tSlice - x))]})
#   
#   rarityDF_HT$gpp[WU] <- grEvent
#   rarityDF_HT$ecdf[WU] <- poeEvent
#   
#   # GEV fitted to whole spaceslice
#   FFGEV <- fevd(x=tSlice,
#                 type='GEV')$results$par
#   QFGEV <- 1- pevd(q=tSliceEvent,
#                    loc=FFGEV[1],
#                    scale=FFGEV[2],
#                    shape=FFGEV[3],
#                    type='GEV')
#   
#   rarityDF_HT$gev[WU] <- QFGEV
#   
#   if(n %% 200 == 0){
#     print(floor(proc.time() - ST0)[1:3])
#     print(floor(proc.time() -  ST)[1:3])
#   }
#   ST <- proc.time()
# }


print("STEP 6 COMPLETE.")
print(Sys.time() - ST)
ST <- Sys.time()
print("END")
print(Sys.time())