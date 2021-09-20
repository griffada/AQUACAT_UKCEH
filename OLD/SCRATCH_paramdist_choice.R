commandArgs <- function(...){c("01","present")}
  
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

library(ilaprosUtils)
library(extRemes)
library(lmomco)

threshDayExcList <- lapply(seq_len(length(threshVal)),
                           function(i){vector("list", NH)})


ST <-  Sys.time()
ncin <- nc_open(ncoriginal, readunlim=FALSE)
  # this is a huge file, do not open without reason.
print(Sys.time() - ST)
print(ncin)


##### THRESHOLD EXTRACT #####------------------------------------
## Get grid of river network ##-------------------
tStart <- 1
deltaT <- 1
###### Only uncomment if needed #####
# ncwide <- ncvar_get(ncin, "dmflow",
# 					start=c(1, 1, tStart),
# 					count=c(-1, -1, deltaT))
# 
# rn <- which(apply(ncwide, c(1,2),
#                  function(v){sum(v[!is.na(v)] > -1) == deltaT}), arr.ind=T)
# 
# rn <- data.frame(row=rn[,1], col=rn[,2])
# east <- ncvar_get(ncin, "Easting")
# north <- ncvar_get(ncin, "Northing")
# rn$east <- east[rn[,1]]
# rn$nor <- north[rn[,2]]
################################### 

medianThresh <- function(vals, threshVal=(2/360)){
  vq <- c()
  for(i in 1:30){
    v <- vals[1:(1*360) + (i-1)*360]
    vq[i] <- quantile(v, probs=1 - threshVal)
  }
  return(median(vq))
}

NH <- nrow(rn)
## Initialise dfs and lists ##-------------------------
threshMat <- matrix(NA, ncol=NT, nrow=NH)

partable <- data.frame(meanint=numeric(),
                        threshold=numeric(),
                        loc=numeric(),
                        scale=numeric(),
                        shape=numeric(),
                        threshquan=numeric(),
                        pot2=numeric())

##### Get quantiles for thresholds #####--------------------------------
print("loop start")
ST <- Sys.time()
print(ST)
ST0 <- ST
print(paste("NH =", NH))

h <- 19001

i <- rn$row[h]
j <- rn$col[h]
tSlice <- as.vector(ncvar_get(ncin, varid="dmflow",
                     start=c(i, j,  1),
                     count=c(1, 1, -1)))

# find quantile for threshold
thresh <- quantile(tSlice, prob=c(1 - threshVal[jT]), na.rm=T) #vec
mt <- medianThresh(tSlice, threshVal=(2/360))
threshMat[h,jT] <- mt

NV <- 40

meanint <- 30/NV



k <- jT
#for(k in 1:NT){
# save which days cell n was exceeded.
threshDayExcList[[k]][[h]] <- which(tSlice > mt)
  
# fit GPA to peaks above the threshold
#if(k == jT){ # No point doing it for all 5 thresholds.
    
thresh1b <- ecdf(tSlice)(mt)
    
o1 <- try({
  ep1 <- (ilaprosUtils::extractPeaks(vecObs=tSlice, mintimeDiff=7) == 1)

  min_ep <- ecdf(tSlice)(
    sort(tSlice[which(ep1)], decreasing=T)[min(NV, sum(ep1))])
  #print(min_ep)
})

if(inherits(o1, "try-error") | min_ep < 0.8){
  print(paste("min_ep =", min_ep))
  print(paste(h, " Using UKFE2 and medianThresh"))

  df <- data.frame(Date=1:10800, Flow=tSlice)
  
  o2 <- try({ep1 <- POT_extract_UKFE2(df,
                      div=quantile(tSlice[tSlice < mt], probs=0.67),
                      thresh=thresh1b, Plot=FALSE)
  })
  if(inherits(o2, "try-error")){
    print(n)
    stop("This still doesn't work.")
  }
  ep1 <- 1:10800 %in% ep1$Date
}
thresh0 <- 0.9999 * 
  sort(tSlice[which(ep1)], decreasing=T)[min(NV, sum(ep1))]

ep1 <- ep1 & (tSlice >= thresh0)
  
#}
peak_vals <- tSlice[ep1]
meanint <- 30/sum(ep1)
  # fit a GPA distribution
o <- try({
    at_site_gpa <- pargpa(lmoms(peak_vals), xi=thresh0)
    invisible({cdfgpa(peak_vals, at_site_gpa)})
    at_site_gpa <- at_site_gpa$para
    #lmomco uses opposite shape to fevd
})

f <- function(x){
 q = (1 - ecdf(tSlice)(thresh0))*(1 - x)
 q1 =  1 - exp(-360*q)
 log(1/q1 - 1)
}

g <- function(y){
  q1 =  1 - exp(-360*(1-y))
  log(1/q1 - 1)
}

LM <- lmoms(peak_vals) 
LM$lambdas[1] <- median(peak_vals)
LM$ratios[2] <- LM$lambdas[2]/LM$lambdas[1]
at_site_gpa <- pargpa(lmoms(peak_vals), xi=thresh0)
pygpa = quagpa(prev, at_site_gpa)
cgpa = cdfgpa(peak_vals, at_site_gpa)
pxgpa = (1 - ecdf(tSlice)(thresh0))*(1 - cgpa)
Tgpa = 1/(1-cgpa)

at_site_gpa2 <- pargpa(lmoms(peak_vals))
pygpa = quagpa(prev, at_site_gpa)
cgpa = cdfgpa(peak_vals, at_site_gpa)
pxgpa = (1 - ecdf(tSlice)(thresh0))*(1 - cgpa)
Tgpa = 1/(1-cgpa)

at_site_gpa3 <- pargpa(LM, xi=thresh0)
pygpa3 = quagpa(prev, at_site_gpa3)
cgpa3 = cdfgpa(peak_vals, at_site_gpa3)
pxgpa3 = (1 - ecdf(tSlice)(thresh0))*(1 - cgpa3)


atsiteparglo <- parglo(lmoms(peak_vals))
cglo <- cdfglo(peak_vals, parglo(lmoms(peak_vals)))
pxglo = (1 - ecdf(tSlice)(thresh0))*(1 - cglo)

LM <- lmoms(peak_vals) 
LM$lambdas[1] <- median(peak_vals)
LM$ratios[2] <- LM$lambdas[2]/LM$lambdas[1]
atsiteparglo3 <- parglo(LM)
pyglo3 = quaglo(prev, atsiteparglo3)
#atsiteparglo3$para['xi'] <- median(peak_vals)
cglo3 <- cdfglo(peak_vals, atsiteparglo3)
pxglo3 = (1 - ecdf(tSlice)(thresh0))*(1 - cglo3)


mixedpar <- vec2par(c(thresh0, atsiteparglo$para[2:3]), type='gpa')
cmixed <- cdfgpa(peak_vals, mixedpar)
Tmixed <- 1/(1-cmixed)
pxmixed = (1 - ecdf(tSlice)(thresh0))*(1 - cmixed)

Tglo = 1/(1-cglo)

Yglo <- 1/(1 - exp(-360*pxglo))

Yglo3 <- 1/(1 - exp(-360*pxglo3))

Ygpa <- 1/(1 - exp(-360*pxgpa))

Ygpa3 <- 1/(1 - exp(-360*pxgpa3))

Zglo3 <- 1/(1- exp(-20*pxglo3))


Ymixed <- 1 / (1 - exp(-360*pxmixed))

plot(Yglo, peak_vals)
points(Yglo3, peak_vals, col="blue")
points(Ygpa, peak_vals, col="red")
points(Ygpa3, peak_vals, col="orange")
points(Ymixed, peak_vals, col="green")
points(Zglo3, peak_vals, col="blue", pch=4)



yrev = seq(-8,8,0.1)
trev = 1 + exp(yrev)
prev = 1 - 1/trev
pyglo = quaglo(prev, atsiteparglo)
anrev = 1/(1 - exp(-360*(1-prev)))

pymixed = quagpa(prev, mixedpar)



gringorten <- function(v){
  ((length(v) + 1 - rank(v)) - 0.44)/(length(v) + 0.12)
}
gg = gringorten(peak_vals)

gg2 = (1 - ecdf(tSlice)(thresh0))*(1 - gg)

png("./FFC3.png", width=200, height=200, units="mm", res=150, pointsize=11)
plot(f(prev), pyglo, type='l', xlab="Reduced Variate", ylab="Flow",
     main="Present Return period")
lines(f(prev), pygpa, col=2)
lines(f(prev), pymixed, col=3)
points(f(1-gg), peak_vals)
lines(f(prev), pyglo3, col="blue", lwd=2)

lines(f(prev), pygpa3, col="orange", lwd=2)

#lines(g(prev), pyglo, col=1, lty=2)
axis(3, at=f(prev)[0:8*20+1], labels=round((1 + exp(f(prev)))[0:8*20+1],2))
legend("topleft", legend=c("True GLO", "True GPA", "Mixed distribution", "GL Lmedian", "GPA Lmedian"),
       col=c("black", "red","green","blue", "orange"), lwd=c(1,1,1,2,2))
dev.off()


plot(peak_vals, Ymixed, log='xy', col=3)
points(peak_vals, Ygpa, col=2)
points(peak_vals, Yglo, col=1)



plot(peak_vals, Tgpa, log='xy')
points(peak_vals, Tglo, col=2)
points(peak_vals, Tmixed, col=3)



threshquan <- ecdf(tSlice)(thresh0)

cd40 = cdfgpa(peak_vals, at_site_gpa)

px60 = (1 - ecdf(tSlice)(thresh0))*(1 - cglo)


ax60 = 1 - exp(-360*px60)


allvals_40 = tSlice[tSlice > mt] #n=112
allgpa40 = pargpa(lmoms(allvals_40), xi=mt)
allcdfgpa40 = cdfgpa(allvals_40, allgpa40)
allpg40 = (1 - ecdf(tSlice)(mt))*(1-allcdfgpa40)
allag40 = 1 - exp(-360*allpg40)

ub = quantile(tSlice, prob=0.99)
allvals_40 = tSlice[tSlice > ub] #n=112

allvalsx = sort(allvals_40, decreasing=T)
allgpa40 = pargpa(lmoms(allvalsx), xi=ub)
allcdfgpa40 = cdfgpa(allvals_40, allgpa40)
allpg40 = (1 - ecdf(tSlice)(ub))*(1-allcdfgpa40)
allag40 = 1 - exp(-360*allpg40)
length(allvals_40)
c(sum(allag40 < 0.5), sum(allag40 < 0.1), sum(allag40 < 0.01))

cd60 = cdfgpa(peak_vals, at_site_gpa)
px40 = (1 - ecdf(tSlice)(thresh0))*(1 - cd40)
ax40 = 1 - exp(-360*px40)

allvals_60 = tSlice[tSlice > mt] #n=182

atglo40 = parglo(lmoms(peak_vals40))
cdfglo40 = cdfglo(peak_vals40, atglo40)
atglo60 = parglo(lmoms(peak_vals60))
cdfglo60 = cdfglo(peak_vals60, atglo60)
pg60 = (1 - ecdf(tSlice)(4.899))*(1-cdfglo60)
pg40 = (1 - ecdf(tSlice)(5.66))*(1-cdfglo40)
p4060 = 1 - ecdf(tSlice)(peak_vals60[peak_vals60 < 5.66])

ag60 = 1 - exp(-360*pg60)
ag40 = 1 - exp(-360*pg40)

#gloALL = glod.fit(tSlice)
atgloALL = parglo(lmoms(tSlice))
cdfgloALL = cdfglo(tSlice, atgloALL)
pgALL = (1-cdfgloALL)

agALL = 1 - exp(-360*pgALL)


plot(peak_vals60, 1/ax60, ylim=c(0,60))
points(peak_vals40, 1/ax40, col=2)
points(peak_vals60, 1/ag60, pch=4)
points(peak_vals40, 1/ag40, col=2, pch=4)
points(tSlice[tSlice > 5], 1/agALL[tSlice > 5], col=3)
points(peak_vals60, 1/(1 - ecdf(tSlice)(peak_vals60)), col=4)
points(allvals_40, 1/allag40, col=8)

c(sum(allag40 < 0.5), sum(allag40 < 0.1), sum(allag40 < 0.01))

plot(tSlice, col=1+ep1)





