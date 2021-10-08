
x <- 0:500/50
f <- pnorm(x, 5, 2)
g <- pnorm(x, 5, 1)
h <- sqrt(f*g)
plot(x,f,type='l', ylim=c(0,1))
lines(x,g,col=2)
lines(x,h,lty=2)
abline(h=0.64)


u <- qnorm(0.68, 5, 2)
l <- qnorm(0.64, 5, 1)
fg <- approx(c(l,u), c(pnorm(l,5,1), pnorm(u,5,2)), xout=seq(u,l,length.out=20))
lines(fg$x, fg$y, col=3)


library(lmomco)
x <- 0:2000
parsGP <- vec2par(c(90,12,-0.1), type='gpa')
parsGL <- vec2par(c(110,12,-0.3), type='glo')

obs_events <- sort(rlmomco(500, vec2par(c(92,11.5,-0.15), type='gpa')))
parsGP <- pargpa(lmoms(obs_events))
parsGL <- parglo(lmoms(obs_events))

ygpa <- 1-cdfgpa(x, gpa_par)
yglo <- 1-cdfglo(x, glo_par)

hm <- sqrt(ygpa*yglo)
rplim <- -log(1-(0.02))/360
u0 <- quaglo(1-rplim, glo_par)
l0 <- quagpa(1-rplim, gpa_par)
if(u0 <= l0){
  l <- 1.5*quagpa(1-rplim, gpa_par)
  u <- u0
}else{
  u <- 1.5*quaglo(1-rplim, glo_par) 
  l <- l0
}
signif(c(u0,u,l0,l,1-cdfgpa(l, gpa_par), 1-cdfglo(u, glo_par)),4)

fg <- approx(c(l,u),
             c(1-cdfgpa(l, gpa_par), 1-cdfglo(u, glo_par)),
             xout=x)
fg1 <- fg
fg1$y[fg1$x < l] <- -1
fg1$y[fg1$x > u] <- 2



plot(obs_events, glo_poe0, type='p', log='y', ylim=c(1e-6, 1))
points(obs_events, gpa_poe0, col=2)
points(obs_events, gpa_geom, col=3)
points(lint$x, lint$y, col=4)
points(obs_events, pmin(lint$y, gpa_geom), pch=3, lty=2)
points(obs_events, pmin(glo_poe0, gpa_poe0), pch=2, col=2)
points(obs_events, gpa_li, col=5, pch=4)
abline(v=glim[1])
abline(h=mid_changept, col=2)
abline(v=glim[2])

max(1/glo_poe0)




