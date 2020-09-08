
source("/CodeABG/setup_script_00.R")

folderN <- "JulyV1"
if(!dir.exists(paste0(wd, "CodeABG/InterimData/",folderN))){
  dir.create(paste0(wd, "CodeABG/InterimData/",folderN))
}

tStart <- 1

ncin <- nc_open(nc_rechunk, readunlim=FALSE)

threshGrid <- ncvar_get(ncin, "dmflow", start=c(1,1,tStart),
                        count=c(-1, -1, 1))

threshGridList <- lapply(1:NT, function(i){threshGrid})

print("loop start")
ST <- Sys.time()
print(ST)
ST0 <- ST
for(n in 1:NH){
  print(n)
  # running time diagnosis
  if(n %% 50 == 0){
    I <- Sys.time() - ST
    I0 <- Sys.time() - ST0
    print(paste("Percent remaining", (NH-n)/NH))
    print(paste("Time remaining", (NH-n)/n * I0))
    print(paste("Since last readout:", I)); ST <- Sys.time() 
  }
  
  i <- rn[n,1]
  j <- rn[n,2]
  tSlice <- ncvar_get(ncin, varid="dmflow",
                      start=c(i, j,  1),
                      count=c(1, 1, -1))
  # find quantile for threshold
  thresh <- quantile(as.vector(tSlice), prob=c(1 - threshVal), na.rm=T) #vec
  threshMat[n,] <- thresh
  for(k in 1:NT){
    # add threshold k at position (i,j) to raster k
    threshGridList[[k]][i,j] <- thresh[k]
    # save which days cell n was exceeded.
    threshDayExcList[[k]][[n]] <- which(tSlice > thresh[k])
  }
}


saveRDS(threshGridList, file=paste0(wd, "CodeABG/InterimData/",folderN, 
                                    "/threshGridList.rds"))