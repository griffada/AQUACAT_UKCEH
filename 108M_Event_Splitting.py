# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 16:54:00 2021

@author: adagri
"""

import sys
import os
import netCDF4 as nc
import pandas as pd
import numpy as np
import re
import yaml
import time

#rcm="04"
#period="198012_201011"

rcm = sys.argv[1]

if sys.argv[2] == "present":
    period = "198012_201011"
else:
    period = "205012_208011"
    
RCMS = ["01","04","05","06","07","08","09","10","11","12","13","15"]
periods = ["198012_201011","205012_208011"]

# CHANGE THIS TO THE TOP LEVEL OF THE FOLDER THE CSVs ARE IN
toplevel = r'S:/Data' #r'/prj/aquacat/Data' #r'C:/Users/adagri/Dropbox/2021 Mar 12 - from CEH'

# CHANGE THIS TO THE TOP LEVEL OF THE FOLDER THE NETCDFs ARE IN
outlevel = toplevel #r'S:/Data'
    
# yaml_path = f"{outlevel}/RCM{rcm}_{period}/settings.yaml"
# with open(yaml_path) as ym:
    # list_doc = yaml.safe_load(ym)

# if list_doc['HTsplit']:
    # print("Split already done. Closing 108.")
# else:
    # CHANGE THIS TO WHERE THE hasData files are, they should exist in the toplevel folder.
rn = pd.read_csv("".join([toplevel,"/hasData_primary.csv"]))
rnreg = pd.read_csv("".join([toplevel,"/hasData_Regions.csv"]))
rnNW = np.where(rnreg.REGION == "NW")[0].tolist()
rn1 = rn[rnreg.REGION == "NW"]

thresh_path = f"{outlevel}/RCM{rcm}_{period}/threshMat_RCM{rcm}_{period}.csv"
threshvec = pd.read_csv(thresh_path).iloc[rnNW,1]

subfold=''
fileinfix = 'POT2_pc01'

print("Finished Setup")

for method in ["OBS", "EC"]:
    #method="OBS"
    print(f"{method} national events open.")
    
    ncpath = f"{outlevel}/RCM{rcm}_{period}{subfold}/NW/event{method}_region_NW_RCM{rcm}_{period}.nc"
    
    ec_events = nc.Dataset(ncpath, mode='r')
    vvec_all = ec_events.variables['flow'][:,:]
    
    eventNo = ec_events.variables['eventNo'][:]
    
    vvec_reg = vvec_all[:,:]
    avec_reg = ec_events.variables['ape'][:,:]
    avec_mid_reg = ec_events.variables['ape_mid'][:,:]
    dvec_reg = ec_events.variables['dpe'][:,:]
    event_reg = eventNo[:]   
    ec_events.close()
    
    print("national events closed.")
    
    ncpath_reg = f"{outlevel}/RCM{rcm}_{period}{subfold}/NW/event{method}_region_NW_RCM{rcm}_{period}B.nc"
    ncfile = nc.Dataset(ncpath_reg, mode='w')
    
    loc_dim = ncfile.createDimension("loc", size=1437)
    event_dim = ncfile.createDimension("event", size=None)
    
    dpe_var = ncfile.createVariable("dpe", np.float32,
                                    dimensions=('event', 'loc'),
                                    complevel=4, chunksizes=[250,4])
    dpe_var.units="PoE"
    dpe_var.long_name="Daily Probability of Exceedance"
    
    ape_var = ncfile.createVariable("ape", np.float32,
                                    dimensions=('event', 'loc'),
                                    complevel=4, chunksizes=[250,4])
    ape_var.units="PoE"
    ape_var.long_name="Annual Probability of Exceedance"
    
    ape_varm = ncfile.createVariable("ape_mid", np.float32,
                                    dimensions=('event', 'loc'),
                                    complevel=4, chunksizes=[250,4])
    ape_varm.units="PoE"
    ape_varm.long_name="Annual Probability of Exceedance (mid risk)"
    
    flow_var = ncfile.createVariable("flow", np.float32,
                                    dimensions=('event', 'loc'),
                                    complevel=4, chunksizes=[250,4])
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
    
    rn2 = np.array(rn1)
    
    row_var[:] = rn2[:,0]
    col_var[:] = rn2[:,1]
    east_var[:] = rn2[:,2]
    north_var[:] = rn2[:,3]
    
    ncfile.RCM = rcm
    ncfile.period = period
    ncfile.event_threshold = "POT2"
    ncfile.area_lower_limit = "pc01"
    ncfile.method='OBS'
    
    ncfile.variables['flow'][:,:] = vvec_reg.T
    ncfile.variables['ape'][:,:] = avec_reg.T
    ncfile.variables['ape_mid'][:,:] = avec_mid_reg.T
    ncfile.variables['dpe'][:,:] = dvec_reg.T
    ncfile.variables['eventNo'][:] = event_reg.data
    ncfile.region="NW"
    ncfile.close()
    
    print(f"nc closed, {method}")
    
yaml_path = f"{outlevel}/RCM{rcm}_{period}/settings.yaml"
with open(yaml_path) as ym:
    list_doc = yaml.safe_load(ym)

list_doc['HTsplit'] = True

with open(yaml_path, 'w') as ym:
    yaml.dump(list_doc, ym)
    
print(time.strftime("%Y-%m-%d %H:%M:%S"))
print(f"108 finished, {rcm}, {period}")