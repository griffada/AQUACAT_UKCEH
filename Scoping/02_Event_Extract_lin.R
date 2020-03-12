######
# Adam Griffin, 2020-02-21
#
# Event extraction from daily flow. Practiced on one 30-year period of data.
# Rewritten for linux to reduce read-write times for netCDF.

# For aquaCAT, Project 07441.
# 
# Created ABG 2020-02-21
#
#####

## SETUP ## ----------------------
library(ncdf4)

## DATA ## -----------------------
ncname_lin <-"/prj/ResilRiskInds/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
# This file is ~36GB on the linux server.
ncin <- nc_open(ncname_lin)
print(ncin)


## EVENT EXTRACT ##
# Get grid of river network
tstart <- 1
Dt <- 5
ncwide <- ncvar_get(ncin, "dmflow", start=c(1,1,tstart), count=c(-1, -1, Dt))
HasData <- which(apply(ncwide, c(1,2),function(v){sum(v[!is.na(v)] > -1) == Dt}), arr.ind=T)
NH <- nrow(HasData)

qpot5_grid <- ncvar_get(ncin, "dmflow", start=c(1,1,tstart), count=c(-1, -1, 1))
# Reset grid to easily spot things
qpot5_grid[!is.na(qpot5_grid) & qpot5_grid > -1] <- 0

print("loop start")
for(n in 1:30){
  if(n %% 5 == 0){print(paste(n, Sys.time()))}
  i <- HasData[n,1]
  j <- HasData[n,2]
  tslice <- ncvar_get(ncin, varid="dmflow",
                       start=c(i, j,  1),
                       count=c(1, 1, -1))
  # find quantile for POT5
  qpot5 <- quantile(as.vector(tslice), prob=c(1 - (5/365)))
  qpot5_grid[i,j] <- qpot5
}

save(qpot5_grid, file="qpot5_grid.RDa")
