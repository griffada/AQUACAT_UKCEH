#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 10:10:42 2020

@author: adagri
"""

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

nc_f = "/prj/ResilRiskInds/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
mainfolder = "/users/sgsys/adagri/aquacat"

ncin = Dataset(nc_f, 'r')
ncshape = ncin.variables['dmflow'].shape  #(10957, 1250, 700)
flow_wide = ncin.variables['dmflow'][0,:,:]

qPOT5 = flow_wide
dx = 50
events_per_year = 5
flow_thres = 1 - (events_per_year/365)
inun_thres = 0.01


flow_long = ncin.variables['dmflow'][:,1,1]
plt.scatter(range(0,ncshape[0]), flow_long.mask)

# == Data ==

QPOT5 = pd.read_csv(
    mainfolder + "/Data/qPOT_" + str(events_per_year) + "epy_sans_mask2.csv", header=None)

Qmask = pd.read_csv(
    mainfolder + "/Data/qPOT_mask2.csv", header=None)

extremes = pd.read_csv(
    "{}/Data/extremedates_{}inun_{}epy2.csv".format(
        mainfolder, str(int(inun_thres*100)), str(events_per_year)), header=1)

# == Plot of QPOT5 ==

grey50 = np.array([220/256, 220/256, 220/256, 1])
cmap = plt.cm.viridis
cmap.set_under(color=grey50)

fig = plt.figure(dpi=300, figsize=(8,11))
ax=fig.add_subplot(111)
#im = ax.imshow(flow_wide, cmap=cmap, vmin=0.1)
im = ax.imshow(QPOT5, cmap=cmap,
              norm=LogNorm(vmin=0.01, vmax=1000))
cbar = plt.colorbar(im, ax=ax)
plt.tight_layout()
plt.show()
# plt.savefig(mainfolder + "/Figures/qPOT5plot_20200303.png")


# == plots of Example widespread events ==
plt.hist(extremes.iloc[:,1], bins=30)
plt.show()
ncin.close()
np.max(extremes.iloc[:,1])
tpos = np.where(extremes.iloc[:,1] == np.max(extremes.iloc[:,1]))[0]
tt = extremes.iloc[1896,0]
flow_slice = ncin.variables['dmflow'][tt,:,:]
qPOT7 = 1*(flow_slice > qPOT5)
qPOT5b = ma.array(qPOT7, mask=flow_wide.mask)
fig = plt.figure(dpi=300, figsize=(8,11))
ax=fig.add_subplot(111)

grey50 = np.array([220/256, 220/256, 220/256, 1])
cmap = plt.cm.viridis
cmap.set_under(color=grey50)
im = ax.imshow(qPOT5b[200:1250, 75:680], cmap=cmap, vmin=0.1)
im = ax.imshow(qPOT5b[200:1250, 75:680], cmap=cmap,
               norm=LogNorm(vmin=0.2, vmax=1000))
cbar = plt.colorbar(im, ax=ax)
plt.tight_layout()
plt.show()