library(lmomco)

pargev_in <- vec2par(c(99, 29, -0.4), type='gev')
peak_vals <- rlmomco(60, pargev_in)

gringorten <- function(v){
  ((length(v) + 1 - rank(v)) - 0.44)/(length(v) + 0.12)
}
gg = gringorten(peak_vals)

ygg <- log(1/gg - 1)


yrev = seq(-8,8,0.1)
trev = 1 + exp(yrev)
prev = 1 - 1/trev

LM <- lmoms(peak_vals)
pg2 <- pargpa(LM, xi=min(peak_vals))
pygpa2 = quagpa(prev, pg2)
pyglo2 <- quaglo(prev, parglo(LM))

pg3 <- pg2
pg3$para['xi'] <- min(peak_vals)*0.99
pygpa3 = quagpa(prev, pg3)

pygev <- pargev(LM)
qv <- quagev(prev, pygev)


LM$lambdas[1] <- (median(peak_vals) + mean(peak_vals))/2
LM$ratios[2] <- LM$lambdas[2]/LM$lambdas[1]
pg <- pargpa(LM, xi=min(peak_vals))
pygpa = quagpa(prev, pg)
pyglo <- quaglo(prev, parglo(LM))

PV2 <- sort(peak_vals, decreasing=T)[-2]
LM2 <- lmoms(PV2)
LM2$lambdas[1] <- (median(PV2) + mean(PV2))/2
LM2$ratios[2] <- LM2$lambdas[2]/LM2$lambdas[1]
pyglo3 <- quaglo(prev, parglo(LM2))

plot(ygg, peak_vals)
lines(yrev, pygpa)
lines(yrev, pygpa2, col=2)
lines(yrev, pyglo, col=3)
lines(yrev, pyglo2, col=4)
lines(yrev, pygpa3, col="purple")
lines(yrev, pyglo3, col="orange")
lines(yrev, qv, col=1, lwd=2)
abline(v=yrev[which.min(abs(trev - 5))])
abline(v=yrev[which.min(abs(trev - 10))])
abline(v=yrev[which.min(abs(trev - 50))])
abline(v=yrev[which.min(abs(trev - 100))])
abline(v=yrev[which.min(abs(trev - 2))])


