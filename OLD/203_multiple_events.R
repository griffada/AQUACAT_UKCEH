#~~~~~~~~~~~~~~~~~~~~~~~
# Adam Griffin
#
# Evaluating event connectedness from observed and simulated events.
#
# For aquaCAT, Project 07441.
#
# OUTPUTS: NewEventHT_***.csv,
#          coefficients.rds, 
#          zscores.rds, 
#          depStruct.rds
# 
# Created ABG 2021-04-09
#
#~~~~~~~~~~~~~~~~~~~~~~~

library(raster)
library(sp)
library(rgdal)
library(texmex)
library(ggplot2)
library(patchwork)


if(interactive()){
  commandArgs <- function(...){c("01","present","NW")}
}
  #args <- c("01","present","NW")
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}

ukshape <- readOGR("/prj/aquacat/Other Data/uk_outline_1000m.shp")


##### Extremal Dependence #####

# This is an investigation on whether the extremal dependence of peak flows changes 
# between now and the future.

obs_events  <- as.data.frame(readr::read_csv(paste0(
  "/prj/aquacat/Data/RCM01_198012_201011/eventflow_OBS_POT2_pc01_RCM01_198012_201011.csv"),
  
 # data_wd,subfold, "eventdpe_OBS_",
  #                        thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"),
                          col_types=cols(.default = col_double())))[,-(1:4)]
obs_events2  <- as.data.frame(readr::read_csv(paste0(
  "/prj/aquacat/Data/RCM01_205012_208011/eventflow_OBS_POT2_pc01_RCM01_205012_208011.csv"),
  
 # data_wd,subfold, "eventdpe_OBS_",
  #                        thresh1,"_", ws1,"_RCM", RCM, suffix, ".csv"),
                          col_types=cols(.default = col_double())))[,-(1:4)]


#obs_events <- obs_events[,-(1:4)]
# M <- 100
# dM <- 120
# SQ <- which(rn_regions$REGION=="NW")[seq(1,1400,by=5)]
# M <- length(SQ)
# rn_sub <- rn[SQ,]
# rn_reg_sub <-rn_regions[SQ,]
# rndist <- as.matrix(dist(rn_sub[,3:4]))/1000
# 
# #max(dist(rn[,3:4]))/1000
# 
# chiapply <- matrix(1, nrow=M, ncol=M)
# chibarapply <- matrix(1, nrow=M, ncol=M)
# chilong <- matrix(NA, ncol=7, nrow=(M^2 - M)/2)
# 
# 
# 
# k <- 0
# for(i in 1:(M-1)){
#   if(i %% 10 == 0){print(i)}
#  for(j in (i+1):M){
#    k <- k+1
#    chi_init <- chi(matrix(c(unlist(obs_events[SQ[i],]), unlist(obs_events[SQ[j],])), ncol=2), nq=100)
#    chibarmean <- apply(tail(chi_init$chibar),2,mean,na.rm=T)
#    chimean <- apply(tail(chi_init$chi),2,mean,na.rm=T)
#    
#    chibarapply[j,i] <- chibarapply[i,j]<- chibarmean[2]
#    
#    if(chibarmean[3] > 0.99){
#      chiapply[j,i] <- chiapply[i,j]<- chimean[2]
#    }else{
#      chiapply[j,i] <- chiapply[i,j] <- 0
#      chimean[2] <- 0
#    }
#   chilong[k,] <- c(chibarmean, chimean, rndist[j,i])
#  }
# 
# }
# 
# chiapply <- matrix(1, nrow=M, ncol=M)
# chibarapply <- matrix(1, nrow=M, ncol=M)
# chilong2 <- matrix(NA, ncol=7, nrow=(M^2 - M)/2)
# k <- 0
# for(i in 1:(M-1)){
#   if(i %% 10 == 0){print(i)}
#  for(j in (i+1):M){
#    k <- k+1
#    chi_init <- chi(matrix(c(unlist(obs_events2[SQ[i],]), unlist(obs_events2[SQ[j],])), ncol=2), nq=100)
#    chibarmean <- apply(tail(chi_init$chibar),2,mean,na.rm=T)
#    chimean <- apply(tail(chi_init$chi),2,mean,na.rm=T)
#    
#    chibarapply[j,i] <- chibarapply[i,j]<- chibarmean[2]
#    
#    if(chibarmean[3] > 0.99){
#      chiapply[j,i] <- chiapply[i,j]<- chimean[2]
#    }else{
#      chiapply[j,i] <- chiapply[i,j] <- 0
#      chimean[2] <- 0
#    }
#   chilong2[k,] <- c(chibarmean, chimean, rndist[j,i])
#  }
# 
# }
# colnames(chilong) <- colnames(chilong2) <- c("chibarlb","chibar","chibarub","chilb","chi", "chiub","dist")
# # plot(chilong2[,3:2], xlab="dist", ylab="chi")
# # points(chilong2[,c(3,1)], col=ifelse(chilong[,2]<0.01,2,3), pch=4)
# chilong <- data.frame(chilong)
# chilong2 <- data.frame(chilong2)
# chilong$period <- "present"
# chilong2$period <- "future"
# chilong$asydep <- (chilong$chi < 0.01)
# chilong2$asydep <- (chilong2$chi < 0.01)
# chilongALL <- rbind(chilong, chilong2)

# par(mfrow=c(2,1))
# plot(chilong[,c(7,5)], xlab="dist", ylab="chi", ylim=c(0,1), col=NA)
# #points(chilong[,c(7,2)], col=ifelse(chilong[,5]<0.01,2,3), pch=4)
# ord <- order(chilong[,7])
# smooth_chilong <- rollapply(chilong[ord,], width=20, FUN=mean, by.column=T)
# polygon(c(smooth_chilong[,7],rev(smooth_chilong[,7])),
#         c(smooth_chilong[,1],rev(smooth_chilong[,3])),col="grey80", border=NA)
# points(chilong[,c(7,2)], col=ifelse(chilong[,5]<0.01,2,3), pch=4, lwd=2)
# points(chilong[,c(7,5)], col=1, pch=1, lwd=2)
# abline(v=max(chilong[(chilong[,5]>0.01),7]))
# 
# plot(chilong2[,c(7,5)], xlab="dist", ylab="chi", ylim=c(0,1), col=NA)
# #points(chilong[,c(7,2)], col=ifelse(chilong[,5]<0.01,2,3), pch=4)
# ord <- order(chilong2[,7])
# smooth_chilong2 <- rollapply(chilong2[ord,], width=20, FUN=mean, by.column=T)
# polygon(c(smooth_chilong2[,7],rev(smooth_chilong2[,7])),
#         c(smooth_chilong2[,1],rev(smooth_chilong2[,3])),col="grey80", border=NA)
# points(chilong2[,c(7,2)], col=ifelse(chilong2[,5]<0.01,2,3), pch=4, lwd=2)
# points(chilong2[,c(7,5)], col=1, pch=1, lwd=2)
# abline(v=max(chilong2[(chilong2[,5]>0.01),7]))
# 
# plot(chilong$dist, chilong$chi, xlab="dist", ylab="chi", ylim=c(0,1), col=ifelse(chilong$chi < 0.01, NA, 1), pch=16)
# points(chilong2$dist, chilong2$chi, col=ifelse(chilong2$chi < 0.01, NA, 2), pch=17)
# points(chilong$dist, chilong$chibar, col=ifelse(chilong$chi < 0.01, 1, NA), pch=1)
# points(chilong2$dist, chilong2$chibar, col=ifelse(chilong2$chi < 0.01, 2, NA), pch=2)

# png(paste0(wd_id,"chibarplot1.png"), type='cairo',width=150, height=150, units='mm', res=200, pointsize=11)
#   p1 <- ggplot(chilongALL, aes(x=dist, y=chi)) + 
#     ylim(0.05,1) +
#     xlim(0,800) + 
#     geom_bin2d(binwidth=c(25,0.025)) + 
#     #geom_point(size=0.01) + 
#     facet_grid(.~period)
#   p2 <- ggplot(chilongALL, aes(x=dist, y=chibar)) + 
#     ylim(-1,1) +
#     xlim(0,800) + 
#     geom_bin2d(binwidth=c(25,0.05)) + 
#     #geom_point(size=0.01) + 
#     facet_grid(.~period)
#   #library(patchwork)
#   p1/p2
# dev.off()

EDENNO <- 6390
rndistEDEN <- apply(rn, 1, function(x){sqrt(sum((x[3:4] - rn[EDENNO,3:4])^2))})

chiapply <- matrix(1, nrow=NH, ncol=1)
chibarapply <- matrix(1, nrow=NH, ncol=1)
chilong <- matrix(NA, ncol=7, nrow=NH)

for(j in 1:NH){
  if(j == EDENNO){next}
  if((j < 10) | (j %% 1000 == 0)){print(j)}

   chi_init <- chi(matrix(c(unlist(obs_events[EDENNO,]), unlist(obs_events[j,])),
                          ncol=2), nq=100)
   chibarmean <- apply(tail(chi_init$chibar),2,mean,na.rm=T)
   chimean <- apply(tail(chi_init$chi),2,mean,na.rm=T)
   
   chibarapply[j] <- chibarapply[j]<- chibarmean[2]
   
   if(chibarmean[3] > 0.99){
     chiapply[j] <- chiapply[j] <- chimean[2]
   }else{
     chiapply[j] <- chiapply[j] <- 0
     chimean[2] <- 0
   }
  chilong[j,] <- c(chibarmean, chimean, rndistEDEN[j])
}

chilong2 <- matrix(NA, ncol=7, nrow=NH)

for(j in 1:NH){
  if(j == EDENNO){next}
  if((j < 10) | (j %% 1000 == 0)){print(j)}

   chi_init <- chi(matrix(c(unlist(obs_events2[EDENNO,]), unlist(obs_events2[j,])),
                          ncol=2), nq=100)
   chibarmean <- apply(tail(chi_init$chibar),2,mean,na.rm=T)
   chimean <- apply(tail(chi_init$chi),2,mean,na.rm=T)
   
   chibarapply[j] <- chibarapply[j]<- chibarmean[2]
   
   if(chibarmean[3] > 0.99){
     chiapply[j] <- chiapply[j] <- chimean[2]
   }else{
     chiapply[j] <- chiapply[j] <- 0
     chimean[2] <- 0
   }
  chilong2[j,] <- c(chibarmean, chimean, rndistEDEN[j])
}

colnames(chilong) <- colnames(chilong2) <- c("chibarlb","chibar","chibarub","chilb","chi", "chiub","dist")
chilong <- data.frame(chilong)
chilong2 <- data.frame(chilong2)
chilong$period <- "present"
chilong2$period <- "future"
chilong$asydep <- (chilong$chi < 0.01)
chilong2$asydep <- (chilong2$chi < 0.01)
chilongALL <- rbind(chilong, chilong2)

print(getwd())
png("./chilongplots2.png", type="cairo")
  p1 <- ggplot(chilongALL, aes(x=dist, y=chi)) + 
    ylim(0.05,1) +
    xlim(0,800) + 
    geom_bin2d(binwidth=c(25,0.025)) + 
    #geom_point(size=0.01) + 
    facet_grid(.~period)
  p2 <- ggplot(chilongALL, aes(x=dist, y=chibar)) + 
    ylim(-1,1) +
    xlim(0,800) + 
    geom_bin2d(binwidth=c(25,0.05)) + 
    #geom_point(size=0.01) + 
    facet_grid(.~period)
  #library(patchwork)
  p1/p2
dev.off()