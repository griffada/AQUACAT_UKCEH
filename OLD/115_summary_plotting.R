######
# Adam Griffin, 2021-02-16
#
# Summarising size and time of observed extreme events.
# 
# For aquaCAT, Project 07441.
# 
# Created ABG 2021-02-16
#
#
#
#####

if(substr(osVersion,1,3) == "Win"){
  source("S:/CodeABG/setup_script_00.R")
}else if (substr(osVersion,1,3) == "Fed"){
  source("/prj/aquacat/CodeABG/setup_script_00.R")
}else{
  source("~/AQUACAT/CodeABG/setup_script_00.R")
}

library(dplyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
display.brewer.all()


# thresh1 <- "POT2"
# ws1 <- "pc05"
# # Only using POT2 and 2% inundation minimums.
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)

Summarylist <- list()
rcms <- c("01","04","05","06","07","08","09","10","11","12","13","15")
for (i in rcms) {

  for (j in c("present","future")) { 
    suffix <- ifelse(j=="present", "_198012_201011", "_205012_208011")
    subfold <- paste0("RCM", i, suffix, "/")
    print(paste(i,j))
    S <- readr::read_csv(paste0(data_wd,subfold, "eventSumm_OBS_",
                      thresh1,"_", ws1,"_RCM", i, suffix, ".csv"),
                      col_types=cols(
                        eventNumber = col_double(),
                        eventDay = col_double(),
                        eventLength = col_double(),
                        area = col_double(),
                        peak = col_double(),
                        season = col_character()))
    S$rcm <- as.numeric(i)
    S$period <- j
    Summarylist[[length(Summarylist)+1]] <- S
  }
}
FullSumm0 <- do.call(rbind, Summarylist)
FullSumm <- FullSumm0[FullSumm0$area>99,]
par(mar=c(3,3,1,1), mgp=c(2,1,0))

#FullSumm <-FullSumm %>%
#  dplyr::filter(rcm==4)# %>%
  #group_by(eventLength, season) %>%
  #summarise(areabar=median(area))

FullSumm$season <- factor(x=FullSumm$season, levels=c("DJF","MAM","JJA","SON"), ordered=T)
FullSumm$period <- factor(x=FullSumm$period, levels=c("present","future"), ordered=T)

ggplot(FullSumm, aes(x=eventLength, y=area)) +
  geom_bin2d() +
  scale_y_log10() +
  scale_fill_viridis(option = "D",  trans = 'log', breaks= c(0.5,5,50,500)) + 
  facet_grid(period~season) +
  labs(x="Duration of Event (days)", y="Area (km\U00B2)")
ggsave(paste0("S:/Plots/allfactors.png"), type='cairo-png', width=160, height=120, units='mm')

ggplot(FullSumm, aes(x=eventLength, y=area)) +
  geom_bin2d() +
  scale_y_log10() +
  scale_fill_viridis(option = "D",  trans = 'log', breaks= c(0.5,5,50,500)) + 
  facet_grid(~period) +
  labs(x="Duration of Event (days)", y="Area (km\U00B2)")
ggsave(paste0("S:/Plots/area_period.png"), type='cairo-png', width=160, height=120, units='mm')

ggplot(FullSumm, aes(x=season)) +
  geom_bar(position="dodge",aes(fill=period)) + 
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Season", y="Number of events", fill="Period")
ggsave(paste0("S:/Plots/events_season.png"), type='cairo-png', width=160, height=120, units='mm')

ggplot(FullSumm, aes(x=eventLength)) +
  geom_histogram(binwidth=2, aes(fill=season)) + 
  coord_cartesian(xlim = c(0, 30)) +
  facet_grid(period~season) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Duration(days)", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/duration_period_season.png"), width=120, height=120, units='mm')

ggplot(FullSumm, aes(x=area)) +
  geom_histogram(aes(fill=season)) + 
  scale_x_log10() +
  #coord_cartesian(xlim = c(0, 5)) +
  facet_grid(period~season) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Area (km\U00B2)", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/area_period_season.png"), width=120, height=120, units='mm')

FullSumm2 <- FullSumm[FullSumm$peak > 0,]
#FullSumm2$peak[FullSumm2$peak == 0] <- 1e-4

ggplot(FullSumm2, aes(x=1/peak)) +
  geom_density(aes(fill=season)) + 
  scale_x_log10() +
  #coord_cartesian(xlim = c(0, 500)) +
  facet_grid(period~season) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Area (km\U00B2)", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/peak_period_season.png"), width=120, height=120, units='mm')