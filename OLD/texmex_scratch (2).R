## Initial exploration of data in subset. --------
plot(nctime, ncData[,1], type='l')

# vector prealloc.
qPOT5 <- rep(NA, NH)
ncPOT5 <- vector("list", NH)
EP <- vector("list", NH)
dfPOT5 <- vector("list", NH)
Nday <- 20 # days over threshold per year on average
probThres <- c(1 - (Nday/365)) 

gpafit <- vector("list", NH)
for(i in 1:NH){
  qPOT5[i] <- quantile(ncData[,i], prob=probThres)
  ncPOT5[[i]] <- (ncData[,i] >= qPOT5[i]) # exceedence flag
  df <- data.frame(time=nctime, flow=ncData[,i])
  EP[[i]] <- decluster(x=df, threshold = qPOT5[i], which.cols=2)
  EP[[i]][EP[[i]] <= qPOT5[i]] <- NA  # hide values below threshold
  df[,2] <- EP[[i]]
  dfPOT5[[i]] <- df
  
  # Fit a Gumbel distribution for return period (to match marginals in HT.)
  dfp <- df[!is.na(df[,2]),]
  gpafit[[i]] <- fevd(x=flow, data=dfp, threshold=qPOT5[i], type="Gumbel")$results
}

K <- 4
qPOT5[K] <- quantile(ncData[,K], prob=probThres)
declus_POT <- decluster(x=ncData[,K], threshold = qPOT5[K])
plot(nctime/365, ncData[,K], type='l')
points(dfPOT5[[K]]$time/365, dfPOT5[[K]]$flow, col=2, pch=1)
declus_times<- declus_POT>qPOT5[K]
points(nctime[declus_times]/365, declus_POT[declus_times], col=3, pch=4)

library(ilaprosUtils)
X <- unlist(subset(dfPOT5[[K]]['flow'], !is.na(dfPOT5[[K]]['flow'])))
Y <- glo.fit(ncData[ncData[,K]>qPOT5[K],K])
par <- glo.fit(X)$mle
F1 <- pglo(X, loc=par[1], scale = par[2], sh=par[3])
(1/(1-F1))
retPlot(Y)
fevd(X, threshold=qPOT5[K], type='GP')

### Issues - Gumbel give Really extreme return periods. GLO gives reasonable ones.

head(dfPOT5[[K]]['flow'])
F0 <-pevd(X, loc=gpafit[[K]]$par[1], scale=gpafit[[K]]$par[2],
          threshold=qPOT5[K],
          type='Gumbel')
1/(1-F0)



probThres <- c(1 - (25/365))
for(i in 5:10){
  qPOT5[K] <- quantile(ncData[,i], prob=probThres)
  declus_POT <- decluster(x=ncData[,i], threshold = qPOT5[i])
  
  q_up <- quantile(ncData[,i], prob=exp(c(-8,-3)))
  
  threshrange.plot(declus_POT, r=q_up, type='GP', nint=31)
  text(0,0, paste0("site",i))
  M <- mrlplot(declus_POT)
  text(0,0, paste0("site",i))
  
}
fevd(ncData[, K], threshold=qPOT5[K], type='GP')

threshrange.plot(declus_POT, r=quantile, type='GP', nint=31, verbose=TRUE)
M <- mrlplot(declus_POT)

qPOT5 <- apply(ncData, 2, FUN=function(x){quantile(x, prob=probThres)})
declus_POT <- decluster(x=ncData[,6:10], threshold = qPOT5[6:10])
