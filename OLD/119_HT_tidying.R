if(interactive()){
  commandArgs <- function(...){c("04","present","NW")}
}
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

library(stringr)

regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE", "SEV", "SSC", "SW", "THA", "TRE", "WAL")
if(length(args)==3){
  RCM <- sprintf("%02d",as.numeric(args[1]))
  period <- args[2]
  if(period=="present"){
    suffix <- "_198012_201011"
  }else if (period=="future"){
    suffix <- "_205012_208011"
  }
  REG <- args[3]
  if(!any(REG == regions)){
    stop(paste("incorrect call: Rscript 109_HeffTawn_Modelling.R gcm period region \n",
    "- Region must be one of: ANG, ESC, NE, NSC, NW, SE, SEV, SSC, SW, THA, TRE, WAL."))
  }
}

r1 <- which(rn_regions$REGION == REG)
NREG <- length(r1)
threshMat <- readRDS(paste0(data_wd, subfold,"threshMat_RCM", 
                             RCM, suffix,".rds"))
thresh0 <- unlist(threshMat[r1,jT], use.names=FALSE)

#### Gluing together 

L <- list.files(paste0(data_wd, subfold, REG), pattern="eventflow_HT.*\\.RDS", full.names=T)
Z <- vector("list", length(L))
for(i in seq_len(length(L))){
  Z[[i]] <- readRDS(L[i])
  Z[[i]] <- Z[[i]][,-(1:4)]
}

simfull <- do.call(cbind, Z)

W <- colnames(simfull)
W <- stringr::str_remove(W,pattern="\\.[0-9]+")
W1 <- sapply(1:length(W), function(i){paste0(W[i],".",sum(W[1:i]==W[i]))})
colnames(simfull) <- W1

simfull <- cbind(rn[r1,], simfull)

readr::write_csv(simfull,
                 paste0(data_wd, subfold, REG, "/eventflow_HT_",REG,"_",
                        thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
rm(Z, simfull)

## User friendly coeff object

COEFFS <- readRDS(paste0(data_wd, subfold, REG,  "/coefficients_", REG, "_",
                            thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))

dx <- 1437
dz <- 6

COEFFSnew <- array(0, dim=c(dz,dx,dx))
for(i in 1:dx){
COEFFSnew[,seq_len(i-1),i] <-  COEFFS[,seq_len(i-1),i]
COEFFSnew[,i,i] <- array(NA, dim=c(1,dz))
COEFFSnew[,i+seq_len(dx-i),i] <- COEFFS[,(i-1)+seq_len(dx-i),i]
}
dimnames(COEFFSnew) <- list(c("a","b","c","d","m","s"), paste0("L",1:dx), paste0("L",1:dx))

COEFFSnew <- matrix(aperm(COEFFSnew, c(2,1,3)), dz*dx, dx)
COEFFSnew <- data.frame("Coefficient"=rep(c("a","b","c","d","m","s"), each=dx), rep(paste0("L",1:dx), times=6), COEFFSnew)
colnames(COEFFSnew) <- c("Coefficient", "Site1", paste0("L",1:dx))

readr::write_csv(COEFFSnew, paste0(data_wd, subfold, REG, "/coefficients_", REG, "_",
                            thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
