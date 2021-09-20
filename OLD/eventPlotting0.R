library(fields)
library(ncdf4)
library(raster)
# library(knitr)
# library(kableExtra)
library(dplyr)

##### SETUP #####--------------------------------------------------------------
if(interactive()){commandArgs <- function(...){c("08","present", "NW")}}

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}
regions <- c("ANG", "ESC", "NE", "NSC", "NW", "SE",
             "SEV", "SSC", "SW", "THA", "TRE", "WAL")
if(length(args)==3){
  RCM <- sprintf("%02d", as.numeric(args[1]))
  period <- args[2]
  if(period=="present"){
    suffix <- "_198012_201011"
  }else if (period=="future"){
    suffix <- "_205012_208011"
  }
  REG <- args[3]
  if(!any(REG == regions)){
    stop(paste("incorrect call: Rscript 110_HT_PoEEstimation.R gcm period region \n",
    "- Region must be one of: ANG, ESC, NE, NSC, NW, SE, SEV, SSC, SW, THA, TRE, WAL."))
  }
}
print(paste("Running for threshold", thresh1, "at", ws1, "minimum spread."))
subfold <- paste0("RCM", RCM, suffix, "/NW/")

suffix_pres <- "_198012_201011"
subfold_pres <- paste0("RCM", RCM, suffix_pres, "/")
ncpres <- paste0(g2g_wd, "dmflow_RCM", RCM, suffix_pres, "_out.nc")

ncin <- nc_open(ncpres)
map1 <- ncvar_get(ncin, start=c(1,1,1), count=c(-1,-1,1))

#### DATA #####----------------------------------------------------------------

# For the RCM/region/period
newEventApe <- readdf(paste0(data_wd,subfold, "eventape_OBS_region_NW",
                             "_RCM",RCM,suffix,".csv"))[,-(1:4)]

newEventApe <- newEventApe[,!is.na(newEventApe[1,])]

#map1 <- raster(paste0(wd_id, "RiverNetwork.asc"))


threshMat <- readRDS(paste0(data_wd, "RCM", RCM, suffix,"/threshMat_RCM", 
                             RCM, suffix,".rds"))

r1 <- which(rn_regions$REGION == REG)
NREG <- length(r1)
thresh0 <- unlist(threshMat[r1,jT], use.names=FALSE)

partable <- readdf(paste0(data_wd,"RCM", RCM, suffix, 
                           "/paramtableG_",thresh1, "_RCM", RCM, suffix, ".csv"),
                            col_types=cols(.default= col_double()))
partable <- partable[r1,]
rn_reg <- as.data.frame(rn[r1,])

eventSumm <- readdf(paste0(data_wd,"RCM", RCM, suffix, 
                           "/NW/eventSumm_OBS_NW_",thresh1, "_pc01_RCM", RCM, suffix, ".csv"),
                            col_types=cols("season"=col_character(),
                                           .default= col_double()))



# shapefiles for all NW files

hyd_areas <- readOGR(paste0(data_wd, "hydrometricAreas/hyd_areas.shp"))
HANW <- hyd_areas[which(hyd_areas$HA_NUM %in% c(68:77, 103)),] %>%
          buffer(dissolve=T)
#plot(HANW)
rnm <- readRDS(paste0(wd_id, "rivernetworkmatrix.rds"))
rnm[rnm > 0] <- 1
rnm_raster <- raster(rnm, xmn=0, xmx=1e6, ymn=0, ymx=7e5)
rnm_NW <- crop(t(rnm_raster), HANW)
rnm_NW <- mask(rnm_NW, HANW)

#plot(rnm_NW, col=heat.colors(20), asp=1)
#lines(HANW)
#### FUNCTIONS ####--------------------------------------------------------

eventPlotterNoMap_Regional <- function(flow_tab, threshold_tab, base_map,
                              filename=NULL, loc_only=FALSE){
  
  NH <- nrow(flow_tab) # number of river network cells
  if(nrow(flow_tab) != nrow(threshold_tab)){
    stop("flow cells should match threshold cells")
  }
  
  rn <- flow_tab[,1:2] # Eastings and Northings
  ex <- extent(base_map)
  base_map <- t(as.matrix(base_map))
  base_map[!is.na(base_map)] <- -0.03
  
  inunLocs <- which(flow_tab[,3] < threshold_tab[,3], arr.ind=T)
  inunLevel <- round((length(inunLocs))/NH, 2) * 100 # % inundation
  
  below_t <- (flow_tab[,3] <= threshold_tab[,3])
  
  if (loc_only){
    flow_tab[!below_t, 3] <- 0.5
    flow_tab[below_t, 3] <- 1
    colpal <- c("grey90", "grey80", "red", "blue")
    brks <- c(-0.04, -0.02, -0.01, 0.99, 1.99)
    MF <- 1
  }else{
    flow_tab[!below_t, 3] <- -0.01
    MF <- ceiling(max(flow_tab[,3]))
    colpal <- c("grey90", "grey80",
                viridis(n=19, alpha=1, begin=0, end=0.75, direction=-1, option="plasma"))
    brks <- c(-0.04, -0.02, seq(0, MF, length.out=20))
  }
  
  for(n in 1:NH){
    base_map[cellFromXY(rnm_NW, rn_reg[n,3:4])] <- flow_tab[n,3]
  }  
  
  flow2 <- as.matrix(base_map)
  
  fgh <- function(){
    par(mar=c(3,3,1,1), mgp=c(2,1,0))
    image.plot(flow2,
             x=seq(ex[1],ex[2],by=1000),
             y=seq(ex[3],ex[4],by=1000),
             zlim=c(-2, MF),
             col=colpal,
             breaks=brks, asp=1,
             xlab="Easting", ylab="Northing",
             axis.args = list(tick=F, labels=F))
    image.plot(flow2,
             z=c(0,MF),
             col=colpal[-(1:2)],
             breaks=brks[-(1:2)],
             legend.only=T)
    pu <- par('usr')
    text(pu[1] + (pu[2] - pu[1])*0.1, pu[3] + (pu[4] - pu[3])*0.9,
         paste0("POT2 threshold\n", inunLevel, "% inundation"),
         adj=c(0,0.5))
  }
  if (!is.null(filename)) {
    png(filename, res=300, width=70, height=70, units='mm', pointsize=9)
    fgh()
    dev.off()
  }
  fgh()
  
  return(flow2)
}
debug(eventPlotterNoMap_Regional)

which(eventSumm$nclusters > 2)


##### PLOTTING #####-------------------------------------------------------

num <- 103
ev <- eventPlotterNoMap_Regional(cbind(rn[r1,1:2], newEventApe[,num]),
                  cbind(rn[r1,1:2], rep(1-0.135, length(r1))),
                  base_map=rnm_NW,
                  filename="./Plots/BiggestArea_OBS_NW.png",
                  loc_only=F)


num <- 9
ev <- eventPlotterNoMap_Regional(cbind(rn[r1,1:2], newEventApe[,num]),
                  cbind(rn[r1,1:2], rep(1-0.135, length(r1))),
                  rnm_NW,
                  filename="./Plots/Peakiest_OBS_NW.png",
                  loc_only=F)

num <- 168
ev <- eventPlotterNoMap_Regional(cbind(rn[r1,1:2], newEventApe[,num]),
                  cbind(rn[r1,1:2], rep(1-0.135, length(r1))),
                  rnm_NW,
                  filename="./Plots/BiggestAPOE_OBS_NW.png",
                  loc_only=F)

num <- 31
ev <- eventPlotterNoMap_Regional(cbind(rn[r1,1:2], newEventApe[,num]),
                  cbind(rn[r1,1:2], rep(1-0.135, length(r1))),
                  rnm_NW,
                  filename="./Plots/Multicluster_OBS_NW_08pres.png",
                  loc_only=F)






