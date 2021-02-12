#####
# Adam Griffin, 2021-01-19
#
# Rejigging coastal extremes; 46 sites
#
#######

library(readr)
library(tidyverse)
library(reshape2)
library(lubridate)

raw_ext <- readr::read_csv("S:/Data/CoastalSumm/extremes.txt", skip=1,
                           col_types="ccdcdcc")
colnames(raw_ext) <- c("Site", "Day_min", "MMIN", "Day_max",
                       "MMAX", "Datum", "Channel")

raw_ext <- raw_ext %>%
  dplyr::select(Site, Day_max, MMAX) %>%
  dplyr::mutate(Date = lubridate::dmy_hm(Day_max),
                month = month(Date),
                year = year(Date))
raw_ext$Day <- floor(time_length(raw_ext$Date - raw_ext$Date[1], "days")) + 1

raw_counts <- raw_ext %>%
  group_by(month, year) %>%
  summarise(eventsum = n())

raw_thresh <- raw_ext %>%
  group_by(Site) %>%
  mutate(thresh = quantile(MMAX, probs=0.75)) %>%
  ungroup() %>%
  dplyr::filter(thresh < MMAX)

event_days <- sort(unique(floor_date(raw_thresh$Date, "day")))
event_conc <- time_length(event_days - event_days[1], "days")

VV <- split(event_days, cumsum(diff(c(-Inf, event_conc)) != 1))
VV2 <- sapply(VV, length)
VV <- VV[VV2 > 1]

df <- raw_ext

df_1day <- raw_ext %>%
  left_join(raw_ext %>% group_by(Day) %>% summarise(Count = n()), by = c('Day'))
# Count columns shows the total number of basins with concurrent AMAX in a 
# given window length (L)

# Order dataframe by Day and Count columns in decreasing order
df_1day <- df_1day[with(df_1day, order(Day, decreasing = TRUE)), ]
df_1day <- df_1day[with(df_1day, order(Count, decreasing = TRUE)), ]
row.names(df_1day) <- NULL

# Start from L = 1-day episodes
df_Ndays <- df_1day

# Search for any other AMAX within n days prior dmax and add this counts to dmax
# dmax dates represent L = 1-day dates

count <- function(this.day, lag) {
  cond1 <- (df_Ndays$Day == this.day - lag)
  cond1[which(df_Ndays$Day == this.day - lag) <
          max(which(df_Ndays$Day == this.day))] <- FALSE
  if (any(cond1)) {
    df_Ndays$Count[df_Ndays$Day == this.day] <<- 
      mean(df_Ndays$Count[condition1]) + df_Ndays$Count[df_Ndays$Day == this.day]
    df_Ndays$Count[cond1] <<- 0     
    # Assign a value of Count = 0 to the AMAX found in the other days that
    # are not dmax dates
  }
}

lag <- seq_len(5)   
# Change lag accordingly to the time window (L) under investigation 
# e.g. lag = seq_len(1) for L = 2-days; lag = seq_len(5) for L = 6-days; etc..
for (i in unique(df_Ndays$Day)) {
  temp <- df_Ndays$Count[df_Ndays$Date == i]
  if (any(temp > 0)) {
    print(temp)
    sapply(lag, count, this.day = i)
  }
}

# Order dataframe by Day and Count_1 columns in decreasing order
df_Ndays <- df_Ndays[with(df_Ndays, order(Day, decreasing = TRUE)), ]
df_Ndays <- df_Ndays[with(df_Ndays, order(Count, decreasing = TRUE)), ]


