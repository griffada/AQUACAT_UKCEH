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

suppressPackageStartupMessages({
library(dplyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
})
#display.brewer.all()

# thresh1 <- "POT2"
# ws1 <- "pc05"
# # Only using POT2 and 2% inundation minimums.
# jT <- which(threshName==thresh1)
# jW <- which(wsName == ws1)

Summarylist <- list()
rcms <- c("01","04","05","06","07","08","09","10","11","12","13","15")
for (i in rcms) {

  for (j in c("present","future")) {

    for (k in c("OBS", "EC")) {
      suffix <- ifelse(j=="present", "_198012_201011", "_205012_208011")
      subfold <- paste0("RCM", i, suffix)
      print(paste(i,j,k))

      if(file.exists(paste0(
          data_wd,subfold, "/eventSumm_",k,"_", thresh1,
          "_",ws1,"_RCM",i,suffix,".csv"))){

        S <- readdf(
          paste0(
            data_wd,subfold, "/eventSumm_",k,"_", thresh1,
            "_",ws1,"_RCM",i,suffix,".csv"))
      S$rcm <- as.numeric(i)
      S$period <- j
      S$method <- k
      #head(S)
      Summarylist[[length(Summarylist) + 1]] <- S
      }
    }
  }
}
# for(i in 1:length(Summarylist)){
#  if(ncol(Summarylist[[i]]) == 13) 
#  Summarylist[[i]] <- Summarylist[[i]][,-6]
# }
FullSumm <- do.call(rbind, Summarylist)

HTsumm <- readdf("S:/Data/eventSumm_HT_NW_POT2_pc01_ALL.csv")
HTsumm$method <- "HT"
HTsumm <- HTsumm[,-12]

FullSumm <- rbind(FullSumm, HTsumm)

#FullSumm <- FullSumm0[FullSumm0$area>99,]
par(mar=c(3,3,1,1), mgp=c(2,1,0))

#FullSumm <-FullSumm %>%
#  dplyr::filter(rcm==4)# %>%
  #group_by(eventLength, season) %>%
  #summarise(areabar=median(area))

FullSumm$season <- factor(x=FullSumm$season, levels=c("DJF","MAM","JJA","SON"), ordered=T)
FullSumm$period <- factor(x=FullSumm$period, levels=c("present","future"), ordered=T)
FullSumm$method <- factor(x=FullSumm$method, levels=c("OBS", "EC", "HT"), ordered=T)
FullSumm$RPA <- 1/FullSumm$peakA_mid
# FullSumm <- FullSumm[!is.na(FullSumm$area),]
# 
# readr::write_csv(FullSumm, "./Data/FullSumm_plotting_national.csv")
FullSumm <- readr::read_csv("./Data/FullSumm_plotting_national.csv")

# FullSumm$weight <- sapply(levels(FullSumm$method)[FullSumm$method],
#               function(x){switch(x, "OBS"=(1/14778), "EC"=(1/120000), "HT"=1)})

FullSumm$weight <- 0 
FullSumm$weight[FullSumm$method == "OBS" & FullSumm$period=="present"] <- 
  1/sum(FullSumm$method == "OBS" & FullSumm$period == "present")
FullSumm$weight[FullSumm$method == "OBS" & FullSumm$period=="future"] <- 
  1/sum(FullSumm$method == "OBS" & FullSumm$period == "future")
FullSumm$weight[FullSumm$method == "EC" & FullSumm$period == "present"] <-
                  1/sum(FullSumm$method == "EC" & FullSumm$period == "present")
FullSumm$weight[FullSumm$method == "EC" & FullSumm$period == "future"] <-
                  1/sum(FullSumm$method == "EC" & FullSumm$period == "future")
FullSumm$weight[FullSumm$method == "HT" & FullSumm$period == "present"] <-
                  1/sum(FullSumm$method == "HT" & FullSumm$period == "present")
FullSumm$weight[FullSumm$method == "HT" & FullSumm$period == "future"] <-
                  1/sum(FullSumm$method == "HT" & FullSumm$period == "future")

FullSumm$peakA_mid[FullSumm$peakA_mid > 1000] <- 1100

al <- as_labeller(c("OBS"="G2G", "EC"="EC", "HT"="HT"))

pl <- as_labeller(c("present"="1980-2010", "future"="2050-2080"))

ggplot(FullSumm[FullSumm$method != "HT",], aes(x=1/peakA_mid, y=area+1)) +
  geom_bin2d() +
  theme_bw() + 
  scale_y_log10() +
  scale_x_log10(limits=c(1,4900)) +
  scale_fill_viridis(option = "D",  trans = 'log', breaks= c(0, 1, 10, 100, 1000)) + 
  facet_grid(period~method, labeller=labeller(period=pl, method=al)) +
  labs(x="Return Period", y="Area (km\U00B2)", colour="Number of events")
ggsave(paste0("S:/Plots/area_rpa_period5.png"),
       type='cairo-png', width=160, height=120, units='mm')

ggplot(FullSumm[FullSumm$method == "OBS",], aes(x=eventLength, y=1/peakA)) +
  geom_bin2d(binwidth=c(1,0.1)) +
    theme_bw() + 
  scale_y_log10(limits=c(1,300)) +
  xlim(0,14) +
  scale_fill_viridis(option = "D",  trans = 'log', breaks= c(1,10,100)) + 
  facet_grid(period~season, labeller=labeller(period=pl)) +
  labs(x="Duration of Event (days)", y="Return Period (years)")
ggsave(paste0("S:/Plots/area_period5.png"),
       type='cairo-png', width=160, height=120, units='mm')


ggplot(FullSumm, aes(x=season)) +   # Needs to be proportions
  geom_bar(position="dodge",aes(fill=period, weight=weight)) + 
    theme_bw() + 
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  facet_grid(~method, labeller=al) +
  labs(x="Season", y="Proportion of events", fill="Period")
ggsave(paste0("S:/Plots/events_season4.png"),
       type='cairo-png', width=160, height=120, units='mm')

ggplot(FullSumm, aes(x=eventLength)) +
  geom_histogram(binwidth=1, aes(fill=season, weight=weight)) + 
    theme_bw() + 
  xlim(0, 15) +
  facet_grid(period~season, labeller=labeller(period=pl)) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Duration(days)", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/duration_period_season3.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm[FullSumm$method != "HT" & abs(FullSumm$RPA - 50)>1,], aes(x=RPA)) +
  geom_density(aes(fill=period, weight=weight), adjust=1.55) + 
    theme_bw() + 
  scale_x_log10(limits=c(1,300)) +
  #coord_cartesian(xlim = c(0, 5)) +
  facet_grid(period~method, labeller=labeller(method=al,period=pl)) +
  guides(fill="none") +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Return Period (years)", y="Density of events", fill="Season")
ggsave(paste0("S:/Plots/area_period_season4.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm[FullSumm$method != "HT",], aes(x=area)) +
  geom_histogram(aes(fill=season, weight=weight)) + 
  theme_bw() + 
  scale_x_log10() +
  #coord_cartesian(xlim = c(0, 5)) +
  facet_grid(method~season, labeller=labeller(method=al,period=pl)) +
  guides(fill="none") +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Area (km\U00B2)", y="Density of events", fill="Season")
ggsave(paste0("S:/Plots/area_season4.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm[FullSumm$period == "present" & FullSumm$method == "OBS",], aes(x=area)) +
  geom_histogram() + 
    theme_bw() + 
  scale_x_log10() +
  #coord_cartesian(xlim = c(0, 5)) +
  facet_wrap(~rcm) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Area (km\U00B2)", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/rcm_variability_present_area3.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm[FullSumm$period == "present" &
                  FullSumm$rcm %in% c(1,7,15) & 
                  FullSumm$method != "HT",], aes(x=1/peakA_mid)) +
  geom_density() + 
    theme_bw() + 
  scale_x_log10() +
  #coord_cartesian(xlim = c(0, 5)) +
  facet_wrap(method~rcm) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Return Period (years)", y="Probability Density", fill="Season")
ggsave(paste0("RPA_rcm_method3.png"),
       width=120, height=120, units='mm')

FullSumm2 <- FullSumm[FullSumm$peakA > 1e-3,]
#FullSumm2$peak[FullSumm2$peak == 0] <- 1e-4

ggplot(FullSumm2, aes(x=RPA)) +
  geom_density(aes(fill=season)) + 
    theme_bw() + 
  scale_x_log10() +
  #coord_cartesian(xlim=c(1,750)) +
  facet_grid(period~season) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Return Period (years)", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/RPA_period3.png"),
       width=120, height=120, units='mm')


ggplot(FullSumm, aes(x=season)) +
  geom_bar(position="dodge",aes(fill=period)) + 
    theme_bw() + 
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Season", y="Number of events", fill="Period")
ggsave(paste0("S:/Plots/season period3.png"),
       type='cairo-png', width=160, height=120, units='mm')

ggplot(FullSumm, aes(x=nclusters)) +
  geom_histogram(aes(fill=method, weight=weight), binwidth=1) + 
    theme_bw() + 
  #scale_x_log10() +
  #coord_cartesian(xlim = c(0, 5)) +
  facet_grid(period~method) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Clusters", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/clusters3.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm, aes(x=peakyness, after_stat(density))) +
  geom_histogram(aes(fill=method) , binwidth=0.1) + 
    theme_bw() + 
  #scale_x_log10() +
  #xlim(0,0.99)+
  #coord_cartesian(xlim = c(0, 5)) +
  facet_grid(cols=vars(period)) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Peakyness", y="Number of events", fill="Season")
ggsave(paste0("S:/Plots/peaky_period_season3.png"),
       width=120, height=120, units='mm')


ggplot(FullSumm, aes(x=peakyness, y=area)) +
  geom_bin2d(aes(fill=after_stat(ncount))) +
  theme_bw() + 
  #scale_y_log10() +
  scale_fill_viridis(option = "D",  trans = 'log', breaks= c(0,0.0001,0.001,0.01,0.1,1)) + 
  facet_grid(period~method) +
  labs(x="Peakyness", y="Area (km\U00B2)")
ggsave(paste0("S:/Plots/peaky_period_area_method3.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm, aes(x=area, fill=factor(nclusters))) +
  geom_histogram(position="stack") + 
  theme_bw() + 
  scale_x_log10() +
  #coord_cartesian(xlim = c(0, 5)) +
  facet_grid(period~season) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  labs(x="Area (km\U00B2)", y="Propn of events", fill="N Clusters")
ggsave(paste0("S:/Plots/area2_clusters3.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm, aes(x=peakA, y=peakD)) +
  geom_bin2d(aes(fill=after_stat(ncount))) +
  theme_bw() + 
  #coord_cartesian(ylim = c(0, 25000)) +
  scale_y_log10() +
  scale_x_log10() +
  scale_fill_viridis(option = "D", trans = 'log', breaks= c(0.001,0.01,0.1,1)) + 
  facet_grid(method~season)# +
  labs(x="Peak Annual PoE", y="Peak Daily Annual PoE")
ggsave(paste0("S:/Plots/annVdaily3.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm[FullSumm$method != "HT",], aes(x=1/peakA_mid)) +
  geom_histogram(aes(weight=weight)) + 
  theme_bw() + 
  scale_x_continuous(trans="log10", limits=c(1,4900)) +
  facet_grid(period~method) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Return Period (years)", y="Probability Density")
ggsave(paste0("S:/Plots/peakA_hist4.png"),
       width=120, height=120, units='mm')

ggplot(FullSumm, aes(x=area)) +
  geom_histogram(aes(weight=weight)) + 
  theme_bw() + 
  scale_x_continuous(trans="log10") +
  coord_cartesian(xlim = c(10, 12000)) +
  facet_grid(period~method) +
  guides(fill=FALSE) +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Area (km\U00B2)", y="Probability Density")
ggsave(paste0("S:/Plots/Area_hist3.png"),
       width=120, height=120, units='mm')



FullSumm$RPA2big <- "All"
FullSumm$RPA2big[FullSumm$RPA > 2] <- "Q2"
FullSumm$RPA2big[FullSumm$RPA > 5] <- "Q5"
FullSumm$RPA2big <- factor(FullSumm$RPA2big, levels=c("All","Q2","Q5"), ordered=T)

ggplot(FullSumm, aes(x=area, fill=RPA2big, alpha=0.5)) +
  geom_density() + 
    theme_bw() + 
  scale_x_log10() +
  #coord_cartesian(xlim = c(0, 5)) +
  facet_grid(cols=vars(method)) +
  guides(fill="legend", alpha="none") +
  scale_fill_brewer(palette = "Dark2")  +
  theme_bw() +
  labs(x="Area (km\U00B2)", y="Probability Density", fill="Season")
ggsave(paste0("RPA_rcm_method3.png"),
       width=120, height=120, units='mm')