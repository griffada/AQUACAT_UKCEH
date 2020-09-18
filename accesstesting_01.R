library(ncdf4)

ncin <- nc_open("/prj/aquacat/run_hmfg2g/run_baseline_ukcp18rcm12km/dmflow_out.nc", readunlim=FALSE) 
I <- c(300:304, 600:604)
J <- c(300:304, 600:604)
XT <- matrix(NA,ncol=5,nrow=10)
YT <- matrix(NA,ncol=5,nrow=10)
for(i in 1:10){
  print("timeslice")
  print(XT[i,] <- system.time({
  XX <- ncvar_get(ncin, "dmflow",
                  start=c(1, 1, I[i]),
                  count=c(-1, -1, 1))
  }))
  print("spaceslice")
  print(YT[i,] <- system.time({
    YY <- ncvar_get(ncin, "dmflow",
                    start=c(I[i], J[i], 1),
                    count=c(1, 1, -1))
  }))
}

apply(XT,2,mean)
apply(YT,2,mean)