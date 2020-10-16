
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

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



reRL <- reshape2::recast(inRL, formula=as.formula(eventNo ~ gpp),
                         id.var=c("loc","Easting", "Northing"))

pivotRL <- tidyr::pivot_wider(inRL, names_from=eventNo, values_from=gpp,
                              names_prefix="E")
