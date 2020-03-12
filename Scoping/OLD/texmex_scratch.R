library(texmex)
library(ncdf4)
library(fields)


#### Data -----------------------------------------------------------
ncname <- "K:/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"  # ~ 36GB
ncin <- nc_open(ncname)
print(ncin)

istart <- 390
Di <- 25
jstart <- 400
Dj <- 25
tstart <- 101
Dt <- 2000

# get flow data from ncin.
ncflow <- ncvar_get(ncin, varid="dmflow",
                    start=c(istart,jstart,101), count=c(Di,Dj,Dt))
# for Di = 25, Dj = 25, Dt = 2000, ncflow is ~10MB

for(i in 1:5){
  print(range(ncflow[,,i], na.rm=T))
}

mymex <- mex(winter, mqu = .7, penalty="none", dqu=.7, which = "NO")
plot(mymex)
# Only do 10 replicates to keep CRAN checks happy. Do many more in any
# real application
myboot <- bootmex(mymex, R=10)
plot(myboot)
mypred <- predict(myboot,  pqu=.95)
summary(mypred , probs = c( .025, .5, .975 ))

ncflat <- rbind(ncflow[1,,], ncflow[2,,])
apply(ncflat, 1, function(v){sum(v+1 < 1e-8)})

HasData <- which(apply(ncflow, c(1,2),
                       function(v){sum(v[which(!is.na(v))] > -1) == 2000}), arr.ind=T)
ncData <- matrix(NA, ncol=58, nrow=2000)
for(i in 1:58){
  ncData[,i] <- ncflow[HasData[i,1], HasData[i,2],]
}

mmex <- mex(ncData, mqu=.7, penalty="none", dqu=.7, which=2)
myboot <- bootmex(mmex, R=10)
plot(myboot)
mypred <- predict(myboot,  pqu=.75)
summary(mypred , probs = c( .025, .5, .975 ))


MGPD <- migpd(ncData[,1:6], mqu=.7, penalty="none", maxit=1000)

mexdep <- mexDependence(MGPD, which=1, dqu=.7, margins="gumbel")
mall <- mexAll(ncData[,1:6], mqu=0.7, dqu=0.7)

mmc <- mexMonteCarlo(10, mall)
