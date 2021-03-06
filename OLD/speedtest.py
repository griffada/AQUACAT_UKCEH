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
import time

nc_old = "/prj/ResilRiskInds/hmfg2g_runs/baseline_qx_setting/dmflow_out.nc"
nc_new = "/prj/aquacat/run_hmfg2g/outputs/dmflow_RCM01_198012_201011_out.nc"

tic = time.perf_counter()
ncin = Dataset(nc_new, 'r')
toc = time.perf_counter()
print(f"NEW read: {toc - tic:0.4f} seconds")

tic = time.perf_counter()
flow_wide = ncin.variables['dmflow'][0,:,:]
toc = time.perf_counter()
print(f"NEW wide read: {toc - tic:0.4f} seconds")

for j in range(0, ncshape[2]):
    if not (flow_wide.data[100,j] == -1. or flow_wide.mask[100,j] == True):
        rn = j
        break

tic = time.perf_counter()
flow_sub = ncin.variables['dmflow'][0:ncshape[0], 100, rn]
toc = time.perf_counter()
print(f"NEW deep read: {toc - tic:0.4f} seconds")

tic = time.perf_counter()
ncin = Dataset(nc_old, 'r')
toc = time.perf_counter()
print(f"OLD read: {toc - tic:0.4f} seconds")

ncshape = ncin.variables['dmflow'].shape  #(10800, 1000, 700)

tic = time.perf_counter()
flow_wide = ncin.variables['dmflow'][0,:,:]
toc = time.perf_counter()
print(f"OLD wide read: {toc - tic:0.4f} seconds")

for j in range(0, ncshape[2]):
    if not (flow_wide.data[100,j] == -1. or flow_wide.mask[100,j] == True):
        rn = j
        break

tic = time.perf_counter()
flow_sub = ncin.variables['dmflow'][0:ncshape[0], 100, rn]
toc = time.perf_counter()
print(f"OLD deep read: {toc - tic:0.4f} seconds")