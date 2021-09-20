OMAX <- apply(obs_apoe, 2, min)
dOMAX <- density(log(OMAX))
plot(dOMAX)
lines(dJMAX, col=3)


eventSet <- obs_apoe  
NE <- ncol(obs_apoe)
JMAX <- rep(NA,  length(OMAX))
threshVal0 <- 2/360
threshAnn <- (1-exp(-2))
for(i in 1:length(JMAX)){
  U <- sample.int(NE, 1)
  #print(U)
  V <- c(rep(U,20),sample.int(NE, 2, replace=T))
  eventSubset <-  eventSet[,U]
  su <- sum(eventSubset < threshAnn)
  su <- max(0, min(NH, floor(su*rnorm(1, 1, 0.03))))
  tg <- tailsGen(n=NH, nabove=su, threshVal=threshAnn, eventSet=eventSet[,U])
  JMAX[i] <- min(tg)
}
dJMAX <- density(log(JMAX))
plot(dJMAX)
lines(dOMAX, col=2)

#plot(density(log(tg[tg < threshAnn])), ylim=c(0,1))
lu <- log(eventSet[,U])
lv <- unlist(log(eventSet[,V]))
# lines(density(lu[lu<log(threshAnn)]), col=2)
# lines(density(lv[lv<log(threshAnn)]), col=3)

ht <- hist(log(tg[tg<threshAnn]), breaks=seq(-6,0,by=0.25))
hu <- hist(lu[lu<log(threshAnn)], breaks=seq(-6,0,by=0.25)+0.02)  

hv <- hist(lv[lv<log(threshAnn)], breaks=seq(-6,0,by=0.25)+0.01) 
plot(ht)
plot(hv, border=3)
plot(hu, add=T, border=2)
plot(ht, add=T, border=1)

nexc <- apply(eventSet, 2, function(x){sum(x < threshAnn)})



PPP <- pobs(t(obs_apoe))

CCC <- empCopula(PPP, smoothing="beta")
CCQ <- empCopula(t(obs_apoe), smoothing="beta")

rr <- rCopula(1000,CCC)
rrQ <- rCopula(1000,CCQ)


