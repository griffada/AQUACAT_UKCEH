ai <- c(1, "a", 2)
f <- function(x){ifelse(x==2, stop("Hi."), x)}

for(i in -3:3){
 o <- try({k <- f(i);-10})
 if(inherits(o, "try-error")){
   print(i)
   print("oops")
   warning("Help!")
 }else{
   print(k)
   print(o)
 }
}

vals <- tSlice
vq <- c()
for(i in 1:30){
 v <- vals[1:(1*360) + (i-1)*360]
 vq[i] <- quantile(v, probs=1 - (2/360))
}
vm <- c()
vc <- c()
vz <- c()
for(i in 1:30){
   v <- vals[1:360 + (i-1)*360]
   vm[i] <- sum(v > mean(vq))
   vc[i] <- sum(v > median(vq))
   vz[i] <- sum(v > thresh[jT])
}
plot(vals, log='y')
lines(rep(vq,each=1*360), col=2)