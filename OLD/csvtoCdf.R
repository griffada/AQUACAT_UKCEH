if(interactive()){commandArgs <- function(...){c("05","present", "NW")}}
REG <- NULL
#### SETUP ####----------------------
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

csvToCdf <- function(RCM, period, method, region=NULL, threshold="POT2", cover="pc01"){

  if(is.null(region)){
    rg <- ""
    rg_suff <- paste0("_",threshold,"_",cover)
    nrow <- 19914L
  }else{
    rg <- "NW/"
    rg_suff <- paste0("_region_",region)
    nrow <- 1437L
  }
  ape_path <- paste0(data_wd, subfold, rg, "eventape_", method, rg_suff, "_RCM", RCM, suffix, ".csv")
  dpe_path <- paste0(data_wd, subfold, rg, "eventdpe_", method, rg_suff, "_RCM", RCM, suffix, ".csv") 
  flow_path <- paste0(data_wd, subfold, rg, "eventflow_", method, rg_suff, "_RCM", RCM, suffix, ".csv")      

  dp_top <- length(readdf(dpe_path, nrows=1)) - 4
  
  loc_dim <- ncdim_def(name="loc", units="Num", longname="Location",
                       vals=(1L:nrow))
  event_dim <- ncdim_def("event", units="Num", vals=1:dp_top, 
                         unlim=T, longname="Event")
  dpe_var <- ncvar_def("dpe",
                       "ProbOfExc",
                       list(loc_dim, event_dim),
                       -9999,
                       "Daily Probability of Exceedance",
                       prec="float", compression=9)
  ape_var <- ncvar_def("ape",
                       "ProbOfExc",
                       list(loc_dim, event_dim),
                       -9999,
                       "Annual Probability of Exceedance",
                       prec="float", compression=9)
  flow_var <- ncvar_def("flow",
                        "cumecs",
                        list(loc_dim, event_dim),
                        -9999,
                        "Peak Flow",
                        prec="float", compression=9)
  row_var <- ncvar_def("row",
                       "Num",
                       list(loc_dim),
                       -9999,
                       "Row",
                       prec="integer")
  col_var <- ncvar_def("col",
                       "Num",
                       list(loc_dim),
                       -9999,
                       "Column",
                       prec="integer")
  north_var <- ncvar_def("northing",
                         "metres",
                         list(loc_dim),
                         -9999,
                         "Northing",
                         prec="integer")
  east_var <- ncvar_def("easting",
                         "metres",
                         list(loc_dim),
                        -9999,
                         "Easting",
                         prec="integer")
  event_var <- ncvar_def("eventno",
                         "",
                         list(event_dim),
                         -9999,
                         prec="integer")
  
  savepath <- paste0(data_wd, subfold, rg, "event", method, "_RCM", RCM, suffix, ".nc")
  nc.file <- nc_create(savepath,
    list(row_var, col_var, north_var, east_var, flow_var, ape_var, dpe_var, event_var), force_v4 = T)
  
  print("nc created")
  
  ape <- readdf(ape_path)

  ncvar_put(nc.file, row_var, unlist(ape[,1]))
  ncvar_put(nc.file, col_var, unlist(ape[,2]))
  ncvar_put(nc.file, north_var, unlist(ape[,3]))
  ncvar_put(nc.file, east_var, unlist(ape[,4]))
  ncvar_put(nc.file, ape_var, signif(unlist(ape[,-(1:4)]),5))
  rm(ape)
  dpe <- readdf(dpe_path)
  ncvar_put(nc.file, dpe_var, signif(unlist(dpe[,-(1:4)]),5))
  rm(dpe)
  flow <- readdf(flow_path)
  ncvar_put(nc.file, flow_var, signif(unlist(flow[,-(1:4)]), 4))
  ncvar_put(nc.file, event_var, as.numeric(substr(colnames(flow)[-(1:4)],2,6)))
  rm(flow)
  
  ncatt_put(nc.file,0, "RCM", RCM)
  ncatt_put(nc.file,0, "period", period)
  ncatt_put(nc.file,0, "event threshold", thresh1)
  ncatt_put(nc.file,0, "area lower limit", ws1)
  ncatt_put(nc.file,0, "Region", REG)
  print(nc.file)
  nc_close(nc.file)
  return(savepath)
}

#c1 <- csvToCdf(RCM=RCM, period=period, method="OBS", region=NULL, threshold=thresh1, cover=ws1)
c2 <- csvToCdf(RCM=RCM, period=period, method="EC", region=NULL, threshold=thresh1, cover=ws1)
# c3 <- csvToCdf(RCM=RCM, period=period, method="OBS", region=REG, threshold=thresh1, cover=ws1)
# c4 <- csvToCdf(RCM=RCM, period=period, method="EC", region=REG, threshold=thresh1, cover=ws1)
# c5 <- csvToCdf(RCM=RCM, period=period, method="HT", region=REG, threshold=thresh1, cover=ws1)
# c6 <- csvToCdf(RCM=RCM, period=period, method="EC2", region=NULL, threshold=thresh1, cover=ws1)

