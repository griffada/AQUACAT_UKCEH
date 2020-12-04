if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

COEFFS <- readRDS(paste0(wd_id, "slimline/COEFFS_D.rds"))
dim(COEFFS)

N <- 200
COEFFS <- COEFFS[,1:(N-1),1:N]
COEFFS_wide <- array(NA, dim=c(6,N,N))
for(i in 1:N){
  COEFFS_wide[,-i,i] <- COEFFS[,,i]
}

dim(COEFFS_wide) <- c(6*N,N)

conam <- c("a","b","c","d","m","s")
OU <- unlist(outer(conam, paste0(".L",1:N), FUN="paste0"))
OU <- as.vector(outer(conam, paste0(".L",1:N), FUN="paste0"))
dimnames(COEFFS_wide) <- list(OU,
                                   paste0("LL",1:N))

write.csv(COEFFS_wide, paste0(wd_id, "COEFFS_tall.csv"))
