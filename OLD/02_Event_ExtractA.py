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

#%%
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import numpy.ma as ma
import platform
from datetime import datetime as dt
from scipy.signal import find_peaks, peak_prominences

if platform.system() == "Windows":
    mainfolder = "S:"
else:
    mainfolder = "/prj/aquacat"

nc_f = mainfolder + "/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc"

ncin = Dataset(nc_f, 'r')
ncshape = ncin.variables['dmflow'].shape  #(10800, 1000, 700)
flow_wide = ncin.variables['dmflow'][0,:,:]

qPOT5 = np.array(flow_wide.data)
qPOT5[:,:] = -999
dx = 50
events_per_year = 5
flow_thres = 1 - (events_per_year/365)
inun_thres = 0.01

rn = np.empty((0, 2))
k = 0
for i in range(0, ncshape[1]):
    for j in range(0, ncshape[2]):
        if not (flow_wide.data[i,j] == -1. or flow_wide.mask[i,j] == True):
            #print([i,j])
            rn = np.append(rn, [[i,j]], axis=0)
            k = k + 1

#%%
            
q5 = []

for k in range(rn.shape[1]):
    # extract dmflow in chunks to reduce read-write time    
    tic = time.perf_counter()
    flow_sub = ncin.variables['dmflow'][0:ncshape[0], rn[k,0], rn[k,1]]
    toc = time.perf_counter()
    print(f"Downloaded the tutorial in {toc - tic:0.4f} seconds")
    
    # compute quantile
    q5.append(np.quantile(flow_sub.data, q=flow_thres, axis=0))
    # preserve "non-river" value of -1
    #qPOT5[np.where(flow_wide[i:(i+dx),
    #                j:(j+dx)] == -1)] = -1

rn = np.column_stack((rn, q5))

pd.DataFrame(rn).to_csv("/prj/aquacat/CodeABG/InterimData/things.csv")

ncin.close()