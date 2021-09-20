# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 11:17:32 2021
Author: Adam Griffin, UKCEH
Project: AQUACAT
Script to convert csv files to netcdf for March datasets.
"""
import os
import netCDF4 as nc
import pandas as pd
import numpy as np
import re

RCMS = ["01","04","05","06","07","08","09","10","11","12","13","15"]
periods = ["198012_201011","205012_208011"]

# CHANGE THIS TO THE TOP LEVEL OF THE FOLDER THE CSVs ARE IN
toplevel = r'C:/Users/adagri/Dropbox/2021 Mar 12 - from CEH'

# CHANGE THIS TO THE TOP LEVEL OF THE FOLDER THE NETCDFs ARE IN
outlevel = toplevel #r'S:/Data'

# CHANGE THIS TO WHERE THE hasData files are, they should exist in the toplevel folder.
rn = pd.read_csv("".join([toplevel,"/hasData_primary.csv"]))
rnreg = pd.read_csv("".join([toplevel,"/hasData_regions.csv"]))

method="OBS" # "EC" OR "OBS", "HT" AS NEEDED

regional = False

if regional:
    subfold='/NW'
    fileinfix = 'region_NW'
    rn = rn[rnreg.REGION=="NW"]
else:
    subfold=''
    fileinfix = 'POT2_pc01'

for rcm in RCMS:
    for period in periods:
        
        print(rcm)
        print(period)
        
        foldername = "".join([toplevel,r'/RCM', rcm, '_', period, subfold])
        
        flowfile = "".join([foldername,r'/eventflow_',method,'_',fileinfix,
                            '_RCM',  rcm,'_',period,'.csv'])
        flow = pd.read_csv(flowfile).iloc[:,4:]
        
        NE = flow.shape[1]
        NH = flow.shape[0]
        
        apefile = "".join([foldername,r'/eventape_', method, '_',fileinfix,
                           '_RCM',rcm,'_',period,'.csv'])
        ape = pd.read_csv(apefile).iloc[:,4:]
        
        dpefile = "".join([foldername,r'/eventdpe_', method, '_',fileinfix,
                           '_RCM',rcm, '_',period,'.csv'])
        dpe = pd.read_csv(dpefile).iloc[:,4:]
        
        print('data downloaded')
        
        ncpath = "".join([outlevel, r'/RCM', rcm, '_', period, subfold,
                          r'/event', method, 'Mar_POT2_pc01_RCM',rcm,
                          '_',period, '.nc'])
        ncfile = nc.Dataset(ncpath, mode='w')
        
        loc_dim = ncfile.createDimension("loc", size=NH)
        
        event_dim = ncfile.createDimension("event", size=None)
        
        dpe_var = ncfile.createVariable("dpe", np.float32,
                                        dimensions=('loc', 'event'),
                                        complevel=4, chunksizes=[4,250])
        dpe_var.units="PoE"
        dpe_var.long_name="Daily Probability of Exceedance"
        ape_var = ncfile.createVariable("ape", np.float32,
                                        dimensions=('loc', 'event'),
                                        complevel=4, chunksizes=[4,250])
        ape_var.units="PoE"
        ape_var.long_name="Annual Probability of Exceedance"
        flow_var = ncfile.createVariable("flow", np.float32,
                                        dimensions=('loc', 'event'),
                                        complevel=4, chunksizes=[4,250])
        flow_var.units="cumecs"
        flow_var.long_name="Peak Flow"
        
        event_dim_var = ncfile.createVariable('event', np.int32, ('event',))
        event_dim_var.long_name="Event"
        loc_dim_var = ncfile.createVariable('loc', np.int32, ('event',))
        loc_dim_var.long_name="Location"
        
        event_var = ncfile.createVariable("eventNo", np.int32, dimensions=('event',))
        event_var.long_name = "Matching OBS event"
        
        row_var = ncfile.createVariable("row", np.int32, dimensions=('loc',))
        row_var.long_name="Row"
        
        col_var = ncfile.createVariable("col", np.int32, dimensions=('loc',))
        col_var.long_name="Column"
        
        north_var = ncfile.createVariable("northing", np.int32, dimensions=('loc',))
        north_var.long_name = "Northing"
        north_var.units="m" 
        
        east_var = ncfile.createVariable("easting", np.int32, ('loc',))
        east_var.long_name = "Easting"
        east_var.units="m" 
        
        print('variables defined')
        
        row_var[:] = rn.row
        col_var[:] = rn.col
        east_var[:] = rn.east
        north_var[:] = rn.nor
        
        event_var[:] = [int(re.split('[E\.]',y)[1]) for y in flow.columns]
        
        flow_var[:,:] = flow
        ape_var[:,:] = ape
        dpe_var[:,:] = dpe
        
        ncfile.RCM = rcm
        ncfile.period = period
        ncfile.event_threshold = "POT2"
        ncfile.area_lower_limit = "pc01"
        ncfile.method='OBS'
        
        print('data added')
        
        ncfile.close()
        print('ncfile closed')