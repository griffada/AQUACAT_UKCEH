if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE", "SEV", "SSC", "SW", "THA", "TRE", "WAL")

if(length(args)==3){
  RCM <- sprintf("%02d", args[1])
  period <- args[2]
  if(period=="present"){
    suffix <- "_198012_201011"
  }else if (period=="future"){
    suffix <- "_205012_201011"
  }
  REG <- args[3]
  if(!any(REG == regions)){
    stop(paste("incorrect call: Rscript 109_HeffTawn_Modelling.R gcm period region \n",
               "- Region must be one of: ANG, ESC, NE, NSC, NW, SE, SEV, SSC, SW, THA, TRE, WAL."))
  }
}

ws1 <- "pc05"
thresh1 <- "POT2"
jT <- which(threshName==thresh1)
jW <- which(wsName == ws1)
print(paste("Running for threshold", thresh1, "at ", ws1, "minimum spread."))

step6_test1 <- readRDS(file=paste0(wd_id, "slimline/step6_D.rds"))
NEW_HT <- step6_test1$MC_sample

PRES_REG <- readr::read_csv(paste0(wd_id, "regionalEvents_exampleNW_0.csv"))

eventDF <- readr::read_csv(paste0(data_wd,"TestData/eventdf_POT2_pc05.csv"))
r1 <- which(rn_regions$REGION == REG) # length = 1437
event_region <- eventDF[r1,]

PRES_REG <- readr::read_csv(paste0(data_wd,subfold, "eventdf_region_",
                                   REG,"_",thresh1,"_", ws1, "_RCM", RCM, suffix, ".csv"))
dimnames(PRES_REG) <- list(paste0("LR", 1:ncol(PRES_REG)),
                          paste0("E", 1:nrow(PRES_REG)))

COEFF_tall <- data.frame(from=numeric(),
                         to=numeric(),
                         a=numeric(),
                         b=numeric(),
                         c=numeric(),
                         d=numeric(),
                         m=numeric(),
                         s=numeric())

for(i in 1:dim(COEFFS)[3]){
  if(i %% 25 == 0){print(i)}
  for(j in 1:dim(COEFFS)[2]){
    row <- COEFFS[,j,i]
    row <- c(i, ifelse(i <= j, j+1, j), row)
    COEFF_tall[nrow(COEFF_tall)+1,] <- row
  }
}

MARGINALS <- t(sapply(MODELS, function(x){c(x$par, x$threshold)}))
colnames(MARGINALS) <- c("scale", "shape", "threshold")


rl_db_file <- paste0(data_wd,subfold, "/SQLite/eventdb_HT_region_", REG, "_RCM",
                     RCM, suffix, ".sqlite")

rl_db <- src_sqlite(rl_db_file, create=TRUE)

copy_to(rl_db,
        PRES_REG,
        temporary=FALSE,
        overwrite=TRUE)

copy_to(rl_db,
        COEFF_tall,
        temporary=FALSE,
        overwrite=TRUE)

copy_to(rl_db,
        NEW_HT,
        temporary=FALSE,
        overwrite=TRUE)

copy_to(rl_db,
        MARGINALS,
        temporary=FALSE,
        overwrite=TRUE) 