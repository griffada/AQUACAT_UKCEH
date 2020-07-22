#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 Adam Griffin, 2020-02-21

 Event extraction from daily flow. Practiced on one 30-year period of data.
 Rewritten for linux to reduce read-write times for netCDF.

 For aquaCAT, Project 07441.
 
Created ABG 2020-02-21
Updated ABG 2020-03-05 - continued development
"""

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import numpy.ma as ma
import platform

if platform.system() == "Windows":
    mainfolder = "S:"
else:
    mainfolder = "/prj/aquacat"

nc_f = mainfolder + "/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc"


ncin = Dataset(nc_f, 'r')
ncshape = ncin.variables['dmflow'].shape  #(10957, 1250, 700)
flow_wide = ncin.variables['dmflow'][0,:,:]

qPOT5 = np.array(flow_wide.data)
qPOT5[:,:] = -999
dx = 50
events_per_year = 5
flow_thres = 1 - (events_per_year/365)
inun_thres = 0.01


extractPOT = function(series, datetime=Null, threshold=0, timeOfRise=0):
    thrConst = 2/3
    if (!is.data.frame(series)) {
    if (datetime is Null){
      datetime_real <- seq(series.shape[0])
    }
    datetime_real <- datetime
    series <- pd.DataFrame({"datetime":datetime, "obs":series}})
  } else {
    datetime_real <- sort(series[, 1])
    series$datetime <- order(series[, 1])
  }

  mintimeDiff <- 3*timeOfRise

  #setup
  tt <- series
  names(tt) <- c("time", "obs")
  tt <- tt[order(tt$time), ]

  NR <- nrow(tt)

  #find all the peaks and troughs
  cc <- pastecs::turnpoints(tt$obs)  # cc$pos is timepoints without tied values
  NP <- length(cc$pos)
  sub <- data.frame(time = tt$time[cc$pos],
                    flow = tt$obs[cc$pos],
                    peak = cc$peaks,  # peak = 1, nonpeak = 0
                    pit = cc$pits)  # pits = troughs
  sub$peak[sub$flow < threshold] <- 0  # remove peaks below threshold
  keep <-     rep(0, NP)
  ObsEA <-    (1:NR)[-cc$pos]
  obsPeaks <- which(sub$peak == 1)
  obsPits <-  which(sub$pit  == 1)
  keep[obsPeaks[1]] <- 1  # keep first peak to start with

  #peak validation
  for (i in 2:length(obsPeaks)) {
    now <- obsPeaks[i]
    prev <- max(intersect(obsPeaks[1:(i-1)], which(keep==1)))
    # which still remaining obsPeak occured most recently
    keep[now] <- 1
    tp <- sub$time[prev]
    fp <- sub$flow[prev]
    tn <- sub$time[now]
    fn <- sub$flow[now]
    if (fn == fp) {
      smallest <- prev
    }else{
      smallest <- ifelse(sub$flow[prev] < sub$flow[now], prev, now)
    }
    if ((tn - tp) < mintimeDiff) {
      keep[smallest] <- 0
    }else{
      minThrough <- min(sub$flow[(prev+1):(now-1)])
      if (minThrough > (min(fp, fn) * thrConst)) {
        keep[smallest] <- 0
      }
    }
  }

  # formatting for output
  isPeak <- rep(0, NR)
  isPeak[cc$pos] <- keep
  isPeak[ObsEA] <- 0

  pot <- tt[isPeak == 1, ]
  pot[, 1] <- datetime_real[isPeak == 1]

  return(list(is_peak=isPeak, pot=pot))
}



for i in range(0, ncshape[1], dx):
    for j in range(0, ncshape[2], dx):
        if flow_wide.mask[i:(i+dx), j:(j+dx)].all():
            # mask is true if masking -9999, skip the chunk if all masked
            continue
        else:
            print(str(i) + " " + str(j))
            # extract dmflow in chunks to reduce read-write time
            flow_sub = ncin.variables['dmflow'][0:ncshape[0],
                                                i:min(i+dx, ncshape[1]),
                                                j:min(j+dx, ncshape[2])]
            # compute quantile
            q5 = np.quantile(flow_sub.data, q=flow_thres, axis=0)
            qPOT5[i:(i+dx), j:(j+dx)] = q5
            # preserve "non-river" value of -1
            #qPOT5[np.where(flow_wide[i:(i+dx),
            #                j:(j+dx)] == -1)] = -1

# save with stored Events/yr
pd.DataFrame(qPOT5).to_csv(
    "{}/Data/qPOT_{}epy_sans_mask2.csv".format(mainfolder, events_per_year),
    index=False, header=None)

pd.DataFrame(flow_wide.mask).to_csv(
    "{}/Data/qPOT_mask2.csv".format(mainfolder, events_per_year),
    index=False, header=None)

# == dates of extreme events == 
extremes = []
qPOT6 = np.array(qPOT5)

for tt in range(0, ncshape[0]):
    if tt % 100 == 0:
        print("tt = " + str(tt))
    flow_slice = ncin.variables['dmflow'][tt,:,:]
    inun = np.sum((flow_slice > qPOT6) & (qPOT6 > -1))
    if inun > int(inun_thres*np.sum(qPOT6 > -1)):
        #print(str(tt) + ": " + str(inun) + " / " + str(np.sum(qPOT6 > -1)))
        extremes.append([tt, inun])
        
extremes = pd.DataFrame(extremes, columns=['Day', 'Inun'])
        
pd.DataFrame(extremes).to_csv(
    "{}/Data/extremedates_{}inun_{}epy2.csv".format(
        mainfolder, str(int(inun_thres*100)), str(events_per_year)),
    index=False, header=True)

ncin.close()