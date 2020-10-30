library(tidyverse)
library(dbplyr)
library(RSQLite)
library(reshape2)
library(extRemes)

logit <- function(x){log(x/(1-x))}
invlogit <- function(y){1/(1 + exp(-1*y))}

gringorten <- function(v){
  ((length(v) + 1 - rank(v)) - 0.44)/(length(v) + 0.12)
}

weibull <- function(v){
  (length(v) + 1 - rank(v))/(length(v) + 1)
}


if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

thresh1 <- "POT2" #!#!#!#!# Important constants to select.
ws1 <- "pc05"
print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

threshMat <- readRDS(paste0(wd_id, "threshMat2.rds"))

jV <- which(threshName==thresh1)
jI <- which(wsName == ws1)

inEvents <- read_csv(paste0(data_wd, "RCM01_198012_201011/",
                            "eventdf_POT2_pc05_", "RCM01_198012_201011.csv"),
                     col_types=cols(.default = col_double()))

inRL <- read_csv(paste0(data_wd, "RCM01_198012_201011/",
                        "returnlevels_POT2_pc05_","RCM01_198012_201011.csv"),
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

ECevents <- read_csv(paste0(data_wd, "RCM01_198012_201011/",
                      "NewEventEC_POT2_pc05_","RCM01_198012_201011.csv"),
                     col_types=cols(.default = col_double()))
ECevents <- cbind(rn, ECevents)
colnames(ECevents)[-(1:4)] <- paste0("E",seq_len(ncol(ECevents)-4))
rownames(ECevents) <- NULL
# pivotRL <- tidyr::pivot_wider(inRL, names_from=eventNo, values_from=gpp,
#                               names_prefix="E")


M <- matrix(NA, ncol=429, nrow=NH)
for(i in 1:429){
  M[,i] <- 1 - (1 - inRL$gpp[inRL$eventNo == i])^360
}
M <- cbind(rn, M)
write.csv(M, paste0(data_wd, "RCM01_198012_201011/",
                 "wideEventDF_POT2_pc05_","RCM01_198012_201011.csv"))
colnames(M)[-(1:4)] <- paste0("E",1:429)





HTevents <- read_csv(paste0(wd_id, "slimline/step7_MCsample.csv"))
rn_here <- do.call(cbind,list(rn,rn_regions[, 3:5],LocNum=1:nrow(rn))) %>%
          dplyr::filter(REGION=="NW") %>%
          dplyr::select(row,col,east,nor,LocNum)
HTevents <- cbind(rn_here[1:200, ], t(HTevents))
colnames(HTevents)[-(1:5)] <- paste0("E", 1:(ncol(HTevents) - 5))

NE <- ncol(HTevents) - 5
NH <- 200


NH_reg <- nrow(rn_here)
rarityDF <- expand.grid("eventNo" = seq_len(NE),
                        "loc_reg" = seq_len(length(unique(rn_here$LocNum))))

UU <- unique(rn_here$LocNum)

rarityDF$locNum <- UU[rarityDF$loc_reg]
rarityDF$Easting <- sapply(1:nrow(rarityDF),
                           function(x){rn_here$east[rarityDF$loc_reg[x]]})
rarityDF$Northing <- sapply(1:nrow(rarityDF),
                           function(x){rn_here$nor[rarityDF$loc_reg[x]]})
rarityDF$thresh <- NA
rarityDF$val <- NA
rarityDF$gpp <- NA
rarityDF$ecdf <- NA
rarityDF$gev <- NA

print(ST <- Sys.time())
ncin <- nc_open(ncname) # This file is ~2.5GB on the linux server.
print(ncin)
print(floor(Sys.time() - ST))
ST0 <- proc.time()
ST <- proc.time()
print("loop start")
for(n in 1:NH){
  #if(n %% 200 == 0){
    print(paste(n, "out of", NH))
  #}
  
  i <- rn_here[n,1]
  j <- rn_here[n,2]
  # Pull out spaceslice
  tSlice <- ncvar_get(ncin, varid="dmflow",
                      start=c(i, j,  1),
                      count=c(1, 1, -1))
  
  tSliceEvent <- unname(unlist(HTevents[n, -(1:5)]))
  
  threshval <- threshMat[UU[n], jV] # POT2 column
  
  
  
  rarityDF$thresh[which(rarityDF$loc_reg == n)] <- threshval
  rarityDF$val[which(rarityDF$loc_reg == n)] <- tSliceEvent
  
  # get ecdf and estimate PoE
  
  ecdfSlice <- ecdf(tSlice)
  poeEvent <- 1 - ecdfSlice(tSliceEvent)
  
  
  # Weibull or Gringorten plotting position
  
  grSlice <- gringorten(tSlice)
  grEvent <- sapply(tSliceEvent,
                    function(x){grSlice[which(tSlice == x)[1]]})
  
  rarityDF$gpp[which(rarityDF$loc_reg == n)] <- grEvent
  rarityDF$ecdf[which(rarityDF$loc_reg == n)] <- poeEvent
  
  # GEV fitted to whole spaceslice
  FFGEV <- fevd(x=tSlice,
                type='GEV')$results$par
  QFGEV <- 1- pevd(q=tSliceEvent,
                   loc=FFGEV[1],
                   scale=FFGEV[2],
                   shape=FFGEV[3],
                   type='GEV')
  
  rarityDF$gev[which(rarityDF$loc_reg == n)] <- QFGEV
  
  if(n %% 10 == 0){
    print(floor(proc.time() - ST0)[1:3])
    print(floor(proc.time() -  ST)[1:3])
  }
  ST <- proc.time()
}



rl_db_file <- "./Data/present_return_levels3.sqlite"
rl_db <- src_sqlite(rl_db_file, create=TRUE)

PRES <- readr::read_csv("./Data/present2_returnlevels.csv")

load(paste0(wd_id,"eventLists03.RDa"))

generated_events <- readr::read_csv("./Data/NewEventPresentEC_POT2_pc05.csv")
colnames(generated_events) <- paste0("L", 1:19914)
rownames(generated_events) <- paste0("D", 1:nrow(generated_events))

colnames(PRES)[6] <- "startDay"

locationThresh <- PRES %>%
  dplyr::filter(eventNo == 1) %>%
  dplyr::select(loc, Easting, Northing, thresh)

eventDuration <- data.frame(eventNo     = 1:285,
                            startDay    = eventDayList[[2]][[2]],
                            eventLength = eventLList[[2]][[2]])

presmin <- PRES %>%
  dplyr::select(eventNo, loc, val, gpp, ecdf, gev)

copy_to(rl_db, locationThresh, temporary=FALSE, overwrite=TRUE)
copy_to(rl_db, eventDuration,  temporary=FALSE, overwrite=TRUE)
copy_to(rl_db, presmin,        temporary=FALSE, overwrite=TRUE)

rl_db

eventVloc <- acast(presmin, eventNo ~ loc, value.var="gpp")

colnames(eventVloc) <- paste0("L",1:19914)
rownames(eventVloc) <- paste0("E",1:285)

eventMags <- data.frame(t(eventVloc))

copy_to(rl_db,
        eventMags,
        temporary=FALSE,
        overwrite=TRUE)

copy_to(rl_db,
        generated_events,
        temporary=FALSE,
        overwrite=TRUE)


