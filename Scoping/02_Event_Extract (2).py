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

nc_f = "/prj/ResilRiskInds/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
mainfolder = "/users/sgsys/adagri/aquacat"

ncin = Dataset(nc_f, 'r')
ncshape = ncin.variables['dmflow'].shape  #(10957, 1250, 700)
flow_wide = ncin.variables['dmflow'][0,:,:]

qPOT5 = np.array(flow_wide.data)
qPOT5[:,:] = -999
dx = 50
events_per_year = 5
flow_thres = 1 - (events_per_year/365)
inun_thres = 0.01

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