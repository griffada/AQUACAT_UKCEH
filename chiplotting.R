library(raster)
library(sp)
library(rgdal)
library(texmex)
library(ggplot2)
library(patchwork)
library(viridis)

if(interactive()){
  commandArgs <- function(...){c("01","present","NW")}
}
  #args <- c("01","present","NW")
if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else{
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}


obs_events <- nc_open(paste0(data_wd,subfold, "eventOBSMar_",
                      thresh1,"_", ws1,"_RCM", RCM, suffix, ".nc"))

present_obsflow <- ncvar_get(obs_events, "flow")

suffix <- "_205012_208011"
subfold <- paste0("RCM", RCM, suffix, "/")

obs_events <- nc_open(paste0(data_wd,subfold, "eventOBSMar_",
                      thresh1,"_", ws1,"_RCM", RCM, suffix, ".nc"))
future_obsflow <- ncvar_get(obs_events, "flow")


M <- 100000

chitable <- data.frame(matrix(NA, ncol=8, nrow=0))
colnames(chitable) <-
  c("chibarlb","chibar","chibarub","chilb","chi", "chiub","dist", "period")

for(m in 1:M){
  if(m < 10 | m %% 1000 == 0){print(m)}
  ij <- sample.int(NH, 2)
  i <- ij[1]; j <- ij[2]
  dist <- sqrt(sum((rn[i,3:4] - rn[j,3:4])^2))/1000
  chi_pres <- chi(present_obsflow[,c(i,j)], nq=100)
  chi_futu <- chi(future_obsflow[,c(i,j)], nq=100)
  chibarmean_p <- apply(tail(chi_pres$chibar),2,mean,na.rm=T)
  chimean_p <- apply(tail(chi_pres$chi),2,mean,na.rm=T)
  chibarmean_f <- apply(tail(chi_futu$chibar),2,mean,na.rm=T)
  chimean_f <- apply(tail(chi_futu$chi),2,mean,na.rm=T)
  chitable[2*(m-1)+1,] <- c(chibarmean_p, chimean_p, dist, "1980-2010")
  chitable[2*m,] <- c(chibarmean_f, chimean_f, dist, "2050-2080")
  if(chitable[2*(m-1)+1,3] <= 0.99){
    chitable[2*(m-1)+1,5] <- 0
  }
  if(chitable[2*m,3] <= 0.99){
    chitable[2*m,5] <- 0
  }
}
for(i in 1:7){
chitable[,i] <- as.numeric(unlist(chitable[,i]))
}

print(getwd())
png("./Plots/chilongplots2.png", type="cairo")
  p1 <- ggplot(chitable, aes(x=dist, y=chi)) + 
    ylim(0.05,1) +
    xlim(0,900) + 
    geom_bin2d(binwidth=c(25,0.025)) + 
    scale_fill_viridis(option = "D") +
    theme_bw() + 
    labs(x="Distance (km)",
         y=expression(paste(chi, phantom(10), "(Dependent)")))+
    facet_grid(.~period)
  p2 <- ggplot(chitable, aes(x=dist, y=chibar)) + 
    ylim(-1,1) +
    xlim(0,900) + 
    geom_bin2d(binwidth=c(25,0.05)) + 
    scale_fill_viridis(option = "D") +
    theme_bw() + 
    labs(x="Distance (km)",
         y=expression(paste(bar(chi), phantom(10), "(Independent)"))) +
    facet_grid(.~period)
  p1/p2
ggsave("./Plots/chilongplots2.png", width=120, height=120, units='mm')
dev.off()