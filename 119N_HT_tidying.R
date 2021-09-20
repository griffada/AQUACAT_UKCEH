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

## User friendly coeff object

COEFFS <- readRDS(paste0(data_wd, subfold, REG,  "/coefficients_", REG, "_",
                            thresh1,"_", ws1, "_RCM", RCM, suffix, ".rds"))

dx <- 1437
dz <- 6

dx_dim <- ncdim_def(name="v1", units="", longname="Location1", vals=(1L:dx))
dy_dim <- ncdim_def(name="v2", units="", longname="Location2", vals=(1L:dx))
coeff_dim <- ncdim_def(name="coeffx", units="", longname="coefficient",
                       vals=c("a","b","c","d","m","s"))
coeff_var <- ncvar_def(name="coeff", units="",
                       list(dx_dim, dy_dim, coeff_dim), -9999,
                       "Coefficient Value", prec="double", compression=9)

nc.file <- nc_create(paste0(data_wd, subfold, REG, "/coefficients_", REG, "_",
                            thresh1,"_", ws1, "_RCM", RCM, suffix, ".nc"), force_v4=T)

for(i in 1:dx){
  ncvar_put(nc.file, "coeff", COEFFS[,seq_len(i-1),i], start=c(1,1,i), count=c(-1,(i-1),1))
  ncvar_put(nc.file, "coeff", array(NA, dim=c(1,dz)), start=c(1,i,i), count=c(-1,1,1))
  ncvar_put(nc.file, "coeff", COEFFS[,(i-1)+seq_len(dx-i),i], start(1,(i+1),i), count=c(-1,(dx-i),1))
}

nc_close(nc.file)
