# Adam Griffin, 2020-06-02
# 
# NetCDF speed testing: new vs old versions of files.
#
# 

library(ncdf4)

#oldnc <- "/prj/ResilRiskInds/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
oldnc <- "V:/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
#newnc <- "/prj/aquacat/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc"
#newnc <- "S:/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc"
newnc <- "S:/Data/RCM01_198012_201011/dmflow_copy_RCM01_198012_201011.nc"
#extranc <- "/prj/aquacat/Data/RCM01_198012_201011/dmflow__RCM01_198012_201011.nc"
extranc <- "S:/Data/RCM01_198012_201011/dmflow_rechunked_RCM01_198012_201011.nc"

ST <-  proc.time()
ncin <- nc_open(oldnc)  # this is a huge file, do not open without reason.
print(paste("OLD nc_open:", proc.time() - ST)[3])
#print(ncin)

ptm <- proc.time()
ncwide <- ncvar_get(ncin, "dmflow",
                    start=c(1,1,1),
                    count=c(-1, -1, 1))
print(paste("OLD wide read:", proc.time() - ptm)[3])

wm <- which(!is.na(ncwide) & ncwide > -1, arr.ind=T)

ptm <- proc.time()
nc_deep <-  ncvar_get(ncin, "dmflow",
                      start=c(wm[1,1],wm[1,2],1),
                      count=c(1, 1, -1))
print(paste("OLD deep read:", proc.time() - ptm)[3])

nc_close(ncin)

###

ST <-  proc.time()
ncin <- nc_open(newnc)  # this is a huge file, do not open without reason.
print(paste("NEW nc_open:", proc.time() - ST)[3])
#print(ncin)

ptm <- proc.time()
ncwide <- ncvar_get(ncin, "dmflow",
                    start=c(1,1,1),
                    count=c(-1, -1, 1))
print(paste("NEW wide read:", proc.time() - ptm)[3])

wm <- which(!is.na(ncwide) & ncwide > -1, arr.ind=T)

ptm <- proc.time()
nc_deep <-  ncvar_get(ncin, "dmflow",
                      start=c(wm[1,1],wm[1,2],1),
                      count=c(1, 1, -1))
print(paste("NEW deep read:", proc.time() - ptm)[3])

nc_close(ncin)

###

ST <-  proc.time()
ncin <- nc_open(extranc)  # this is a huge file, do not open without reason.
print(paste("EXTRA nc_open:", proc.time() - ST)[3])
#print(ncin)

ptm <- proc.time()
ncwide <- ncvar_get(ncin, "dmflow",
                    start=c(1,1,1),
                    count=c(-1, -1, 1))
print(paste("EXTRA wide read:", proc.time() - ptm)[3])

wm <- which(!is.na(ncwide) & ncwide > -1, arr.ind=T)

ptm <- proc.time()
nc_deep <-  ncvar_get(ncin, "dmflow",
                      start=c(wm[1,1],wm[1,2],1),
                      count=c(1, 1, -1))
print(paste("EXTRA deep read:", proc.time() - ptm)[3])

nc_close(ncin)






