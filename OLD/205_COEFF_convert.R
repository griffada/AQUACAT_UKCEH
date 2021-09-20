
if(interactive()){
  commandArgs <- function(...){c("12","present","NW")}
}
  #args <- c("01","present","NW")
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}
REG <- "NW"
COEFFS <- readRDS(paste0(data_wd, subfold, "/", REG, "/coefficients_", REG, "_",
                            thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))



dx <- 1437
dz <- 6

COEFFSnew <- array(0, dim=c(dz,dx,dx))
for(i in 1:dx){
  if(i %% 100 == 0){print(i)}
COEFFSnew[,seq_len(i-1),i] <-  COEFFS[,seq_len(i-1),i]
COEFFSnew[,i,i] <- array(NA, dim=c(1,dz))
COEFFSnew[,i+seq_len(dx-i),i] <- COEFFS[,(i-1)+seq_len(dx-i),i]
}
