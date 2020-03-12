#####
# Adam Griffin, 2020-01-27
#
# Preliminary script investigating structure and use of daily flow model 
# outputs.
# 
# Created ABG 2020-01-27
#
#####

##### SETUP ##### --------------------------------------------------------

#### Packages -------------------------------------------------------
library(ncdf4)
library(fields)


#### Data -----------------------------------------------------------
ncname <- "K:/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
ncin <- nc_open(ncname)
print(ncin)

istart <- 390
Di <- 25
jstart <- 400
Dj <- 25


# get flow data from ncin.
ncflow <- ncvar_get(ncin, varid="dmflow", start=c(istart,jstart,1), count=c(Di,Dj,100))

for(i in 1:5){
  print(range(ncflow[,,i], na.rm=T))
}


##### NEW STRUCTURE ##### --------------------------------------------------

# Create a new version of ncin, replacing flow with probability of 
# non-exceedence, based on the period of record.

# subset to test code on.
istart <- 390
Di <- 25
jstart <- 400
Dj <- 25


D <- ncin$dim
D$Northing$len <- Di
D$Easting$len <- Dj

D$Northing$vals <- D$Northing$vals[0:(Di-1) + istart]
D$Easting$vals <- D$Easting$vals[0:(Dj-1) + jstart]


var_poe <- ncvar_def("poe", units="", dim=D, missval=-9999,
                     longname="Probability of exceedence", prec="float")
nc_poe <- nc_create("dmpoe_test.nc", vars=var_poe)
ncatt_get(nc_poe, "poe")
ncatt_put(nc_poe, "poe", attname="missing_value", attval=-9999)

NN <- sapply(1:3, function(n){nc_poe$dim[[n]]$len})

tslice <- ncvar_get(ncin, varid="dmflow", start=c(1,1,500), count=c(-1, -1, 1))
tslice2 <- tslice[390:415, 425:400]
image.plot(tslice2)

for(i in 1:Di){
  for(j in 1:Dj){
    xyslice <- ncvar_get(ncin, varid="dmflow",
                         start=c(i+istart-1, j+jstart-1, 1),
                         count=c(1,1,-1))
    vals <- xyslice
    vals[xyslice == -1] <- -1
    vals[is.na(xyslice)] <- NA
    if(any(!is.na(xyslice) & xyslice > -1)){
      print(paste0("(",i, ", ", j, ")"))
      poe <- ecdf(xyslice[!is.na(xyslice) & xyslice > -1])
      vals[!is.na(xyslice) & xyslice > -1] <- 1 - poe(xyslice[xyslice > -1])
    }
    ncvar_put(nc_poe, "poe", vals,
            start=c(i,j,1),
            count=c(1,1,-1))
  }
}

ticker <- 0
for(k in 1:NN[3]){
  tslice <- ncvar_get(nc_poe, varid="poe",
                      start=c(1,1,k),
                      count=c(-1,-1,1))
  if(any(tslice > 0.5, na.rm=T)){
    ticker <- ticker + 1
    if(ticker < 20){
    print(k)
    image.plot((1-tslice[,Dj:1])/tslice[,Dj:1],
               main=paste0("Time = ", k), zlim=c(-1,1))
    }
    # }else{
    #   break
    # }
  }
}
