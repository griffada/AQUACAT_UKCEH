print("running 102")
if(interactive()){commandArgs <- function(...){c("04","present")}}
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

library(ilaprosUtils)
library(extRemes)
library(lmomco)

ncinp <- nc_open(ncoriginal)
print(ncinp)

if(interactive()){commandArgs <- function(...){c("04","future")}}
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

ncinf <- nc_open(ncoriginal)
print(ncinf)

tStart <- 1
deltaT <- 1

medianThresh <- function(vals, threshVal=(2/360)){
  vq <- c()
  for(i in 1:30){
    v <- vals[1:(1*360) + (i-1)*360]
    vq[i] <- quantile(v, probs=1 - threshVal)
  }
  return(median(vq))
}
NH <- nrow(rn)
NV <- 30
h=900
i <- rn$row[h]
j <- rn$col[h]
tSlicep <- as.vector(ncvar_get(ncinp, varid="dmflow",
                         start=c(i, j,  1),
                         count=c(1, 1, -1)))

tSlicef <- as.vector(ncvar_get(ncinf, varid="dmflow",
                         start=c(i, j,  1),
                         count=c(1, 1, -1)))

threshMat <- matrix(NA, ncol=NT, nrow=NH)

for(tSlice in list(tSlicep, tSlicef)){
  thresh <- quantile(tSlice, prob=c(1 - threshVal[jT]), na.rm=T) #vec
  mt <- medianThresh(tSlice, threshVal=(2/360))
  threshMat[h,jT] <- thresh
  
  thresh1b <- ecdf(tSlice)(mt)
  
  ep1 <- (ilaprosUtils::extractPeaks(vecObs=tSlice, mintimeDiff=7) == 1)
      
  min_ep <- ecdf(tSlice)(
    sort(tSlice[which(ep1)], decreasing=T)[min(NV, sum(ep1))])
  
  
  thresh0 <- 0.9999 * 
        sort(tSlice[which(ep1)], decreasing=T)[min(NV, sum(ep1))]
      
  ep1 <- ep1 & (tSlice >= thresh0)
  
  peak_vals <- tSlice[ep1]
  meanint <- 30/sum(ep1)
  
  PV <- sort(peak_vals, decreasing=T)[-1]
  LM <- lmoms(PV) 
  
  at_site_gpa <- pargpa(LM, xi=thresh0)
  print(1/range(1-cdfgpa(peak_vals, at_site_gpa)))
  #at_site_gpa <- at_site_gpa$para
  
  d2a <- function(x){
    y=(1-x)*(1-ecdf(tSlice)(thresh0))
    p=1 - exp(-360*y)
    p
  }
  
  anngpa <- d2a(cdfgpa(peak_vals, at_site_gpa))
  
  
  LM$lambdas[1] <- median(PV)
  LM$ratios[2] <- LM$lambdas[2]/LM$lambdas[1]
  at_site_glo <- parglo(LM)
  print(1/range(1-cdfglo(peak_vals, at_site_glo)))
  
  annglo <- d2a(cdfglo(peak_vals, at_site_glo))
  #at_site_glo <- at_site_glo$para
  threshquan <- ecdf(tSlice)(thresh0)
  print(c(meanint,thresh0, at_site_glo[1], at_site_glo[2], at_site_glo[3],
          at_site_gpa[1], at_site_gpa[2], at_site_gpa[3],
                     threshquan, threshMat[h,jT]))
}
plot(peak_vals, 1/annglo)
points(peak_vals, 1/anngpa, col=2)
