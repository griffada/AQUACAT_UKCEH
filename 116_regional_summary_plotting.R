######
# Adam Griffin, 2021-02-16
#
# Summarising size and time of observed extreme events.
# 
# For aquaCAT, Project 07441.
# 
# Created ABG 2021-02-16
#####

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) %in% c("Cen","Fed")){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

library(dplyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
#display.brewer.all()

# thresh1 <- "POT2"
# ws1 <- "pc05"
# # Only using POT2 and 2% inundation minimums.
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)

# Summarylist <- list()
# rcms <- c("01","04","05","06","07","08","09","10","11","12","13","15")
# for (i in rcms) {
# 
#   for (j in c("present","future")) { 
#     
#     for (k in c("EC","OBS","HT")) {
#       suffix <- ifelse(j=="present", "_198012_201011", "_205012_208011")
#       subfold <- paste0("RCM", i, suffix, "/NW/")
#       print(paste(i,j,k))
#       
#       if(file.exists(paste0(
#           data_wd,subfold, "eventSumm_",k,"_NW_", thresh1,
#           "_",ws1,"_RCM",i,suffix,".csv"))){
#         
#         S <- readr::read_csv(
#           paste0(
#             data_wd,subfold, "eventSumm_",k,"_NW_", thresh1,
#             "_",ws1,"_RCM",i,suffix,".csv"),
#           col_types = cols(
#             eventNumber = col_double(),
#             eventDay = col_double(),
#             eventLength = col_double(),
#             area = col_double(),
#             peakA = col_double(),
#             peakD = col_double(),
#             season = col_character(),
#             nclusters = col_double(),
#             peakyness = col_double()
#         )
#       )
#       S$rcm <- as.numeric(i)
#       S$period <- j
#       S$method <- k
#       #head(S)
#       Summarylist[[length(Summarylist) + 1]] <- S
#       }
#     }
#   }
# }
# FullSumm <- do.call(rbind, Summarylist)
# #FullSumm <- FullSumm0[FullSumm0$area>99,]
# par(mar=c(3,3,1,1), mgp=c(2,1,0))
# 
# #FullSumm <-FullSumm %>%
# #  dplyr::filter(rcm==4)# %>%
#   #group_by(eventLength, season) %>%
#   #summarise(areabar=median(area))
# 
FullSumm$season <- factor(x=FullSumm$season, levels=c("DJF","MAM","JJA","SON"), ordered=T)
FullSumm$period <- factor(x=FullSumm$period, levels=c("present","future"), ordered=T)
FullSumm$method <- factor(x=FullSumm$method, levels=c("OBS", "EC", "HT"), ordered=T)
FullSumm$RPA <- 1/FullSumm$peakA
FullSumm <- FullSumm[!is.na(FullSumm$area),]
# 
# readr::write_csv(FullSumm, "./Data/FullSumm_plotting.csv")
FullSumm <- readr::read_csv("./Data/FullSumm_plotting.csv")

FullSumm$weight <- sapply(levels(FullSumm$method)[FullSumm$method],
                          function(x){switch(x, "OBS"=(1/6319), "EC"=(1/50436), "HT"=(1/9408))})

ggplot(FullSumm, aes(x=RPA, y=area+1)) +
  geom_bin2d(aes(fill=after_stat(ncount))) +
  xlim(0,1000) +
  theme_bw() + 
  scale_y_log10() +
  scale_fill_viridis(option = "D",  trans = 'log', breaks= c(0, 0.001, 0.01, 0.1, 1)) + 
  facet_grid(period~method) +
  labs(x="Return Period", y="Area (km\U00B2)", fill = "Density")
ggsave(paste0("S:/Plots/allfactors3ec_region.png"),
       type='cairo-png', width=160, height=120, units='mm')

ggplot(FullSumm, aes(x=eventLength, y=area+1)) +
  geom_bin2d(binwidth=c(1,0.1)) +
    theme_bw() + 
  scale_y_log10() +
  scale_fill_viridis(option = "D",  trans = 'log', breaks= c(0.5,5,50,500)) + 
  facet_grid(period~season) +
  labs(x="Duration of Event (days)", y="Area (km\U00B2)")
ggsave(paste0("S:/Plots/area_period2_region.png"),
       type='cairo-png', width=160, height=120, units='mm')

ggplot(FullSumm, aes(x=season)) +   # Needs to be proportions
  geom_bar(position="dodge",aes(fill=period, weight=weight)) + 
    theme_bw() + 
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  facet_grid(~method) +
  labs(x="Season", y="Number of events", fill="Period")
ggsave(paste0("S:/Plots/events_season2_region.png"),
       type='cairo-png', width=160, height=120, units='mm')

ggplot(FullSumm, aes(x=eventLength)) +
  geom_histogram(binwidth=2, aes(fill=season)) + 
    theme_bw() + 
  xlim(0, 13) +
  facet_grid(period~season) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Duration(days)", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/duration_period_season2_region.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm, aes(x=area)) +
  geom_histogram(aes(fill=method, weight=weight)) + 
    theme_bw() + 
  scale_x_log10() +
  #coord_cartesian(xlim = c(0, 5)) +
  facet_grid(period~method) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Area (km\U00B2)", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/area_period_season3_region.png"),
       width=120, height=120, units='mm')

FullSumm2 <- FullSumm[FullSumm$peakA > 0,]
#FullSumm2$peak[FullSumm2$peak == 0] <- 1e-4

ggplot(FullSumm2, aes(x=RPA)) +
  geom_density(aes(fill=season)) + 
    theme_bw() + 
  scale_x_log10() +
  #coord_cartesian(xlim = c(0, 500)) +
  facet_grid(period~season) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Area (km\U00B2)", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/RPA2_region.png"),
       width=120, height=120, units='mm')


ggplot(FullSumm, aes(x=season)) +
  geom_bar(position="dodge",aes(fill=period)) + 
    theme_bw() + 
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Season", y="Number of events", fill="Period")
ggsave(paste0("S:/Plots/season period_region.png"),
       type='cairo-png', width=160, height=120, units='mm')

ggplot(FullSumm, aes(x=nclusters)) +
  geom_histogram(aes(fill=method), binwidth=1) + 
    theme_bw() + 
  #scale_x_log10() +
  #coord_cartesian(xlim = c(0, 5)) +
  facet_grid(period~method) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Clusters", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/clusters2_region.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm, aes(x=peakyness, after_stat(density))) +
  geom_histogram(aes(fill=method) , binwidth=0.1) + 
    theme_bw() + 
  #scale_x_log10() +
  xlim(0,0.99)+
  #coord_cartesian(xlim = c(0, 5)) +
  facet_grid(period~method) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Peakyness", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/peaky_period_season2_region.png"),
       width=120, height=120, units='mm')


ggplot(FullSumm, aes(x=peakyness, y=area)) +
  geom_bin2d(aes(fill=after_stat(ncount))) +
  theme_bw() + 
  #scale_y_log10() +
  scale_fill_viridis(option = "D",  trans = 'log', breaks= c(0,0.0001,0.001,0.01,0.1,1)) + 
  facet_grid(period~method) +
  labs(x="Peakyness", y="Area (km\U00B2)")
ggsave(paste0("S:/Plots/peaky_period_area_method3_region.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm, aes(x=area, fill=factor(nclusters))) +
  geom_histogram(position="stack") + 
  theme_bw() + 
  scale_x_log10() +
  #coord_cartesian(xlim = c(0, 5)) +
  facet_grid(period~season) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  labs(x="Area", y="Propn of events", fill="N Clusters")
ggsave(paste0("S:/Plots/area2_clusters_region.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm2, aes(x=peakA, y=peakD)) +
  geom_bin2d() +
  theme_bw() + 
  #coord_cartesian(ylim = c(0, 25000)) +
  scale_y_log10() +
  scale_x_log10() +
  scale_fill_viridis(option = "D", trans = 'log', breaks= c(0.5,5,50,500,5000)) + 
  facet_grid(method~season) +
  labs(x="Peak Annual PoE", y="Peak Daily Annual PoE")
ggsave(paste0("S:/Plots/annVdaily_region.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm, aes(x=RPA, after_stat(density))) +
  geom_histogram(aes(fill=method)) + 
  theme_bw() + 
  scale_x_continuous(trans="log10") +
  #coord_cartesian(xlim = c(0, 1000)) +
  facet_grid(period~method) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Return Period (years)", y="Event Density")
ggsave(paste0("S:/Plots/peakA_method3_region.png"),
       width=120, height=120, units='mm')