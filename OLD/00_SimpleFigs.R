# Adam Griffin, 2020-02-10
# Script to generate simple figures for reports and presentations


##### SETUP ##### -------------------------------------------------------------
library(mvtnorm)
library(here)
here::here()

library(lattice)
library(ggplot2)
library(fields)


##### SIMPLE BELL CURVES ##### ------------------------------------------------

x <- seq(-30, 30, by=0.1)
y1 <- curve(30*exp(-0.5*(x^2)/20), from=-30, to=30)
y2 <- curve(12*exp(-0.5*((x/8)^4)/60), from=-30, to=30)

png(here::here("DataVis/Plots/bellcurves2.png"), width=120, height=80, units='mm', 
               pointsize=10, res=300)
  par(mfrow=c(2,1), mar=c(3.1,3.1,1,1), mgp=c(2,1,0))
  plot(NA, NA, xlim=c(-30,30), ylim=c(0,30), yaxt='n',
       xlab="Distance from centre (km)", 
       ylab="Return Period (yrs)",frame.plot=F, cex.lab=0.8)
  polygon(c(-30, y1$x,30), c(0, y1$y, 0),col="#9ebdf0", border=NA)
  axis(2, at=c(0,10,20,30), labels=c(0,25,50,75))
  text(30,5,"Small area with very extreme RP at centre", pos=2, cex=0.8)
  plot(NA, NA, xlim=c(-30,30), ylim=c(0,30), yaxt='n',
       xlab="Distance from centre (km)", 
       ylab="Return Period (yrs)", frame.plot=F, cex.lab=0.8)
  axis(2, at=c(0,10,20,30), labels=c(0,25,50,75))
  polygon(c(-32,y2$x,32), c(0,y2$y, 0), col="#86a86d", border=NA)
  text(30,5,"Wide area with less extreme, consistent RP across area",
       pos=2, cex=0.8)
dev.off()
