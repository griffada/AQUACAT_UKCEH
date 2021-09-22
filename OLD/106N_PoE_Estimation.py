# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 13:04:10 2021

@author: adagri
"""
import os
import netCDF4 as nc
import pandas as pd
import numpy as np
import re
import sys
import yaml
import time
from scipy.stats import genpareto
from statsmodels.distributions.empirical_distribution import ECDF


def dpeApeCalculator(h, vals, obs_events, pars, thresh_val):
    events_per_year = 2
    rare_limit = 1e-3
    vals0 = list(vals.data)
    vals0.extend([-1, 1e8])
    ecd = ECDF(vals0)
    gpa_poe = 1 - genpareto.cdf(obs_events,
                                loc=pars[0],
                                scale=pars[1],
                                c=pars[2])
    ecd_poe = 1 - ecd(obs_events)
    gpa_poe[np.isnan(gpa_poe)] = 1
    whext = np.array(obs_events) > thresh_val
    gpa_poe[whext] = gpa_poe[whext]*(1-ecd(thresh_val))
    gpa_poe[(~whext)] = ecd_poe[(~whext)]
    
    valsape = 1 - np.exp(-360*gpa_poe)
    
    valsape[valsape < 1./5000.] = 1./5000.
    gpa_poe[gpa_poe < 5e-7] = 5e-7
    
    return([valsape, gpa_poe])


#rcm="04"
#period="205012_208011"

rcm = sys.argv[1]

if sys.argv[2] == "present":
    period = "198012_201011"
else:
    period = "205012_208011"

print(f"Running RCM {rcm} for {period}.")

# CHANGE THIS TO THE TOP LEVEL OF THE FOLDER THE CSVs ARE IN
toplevel = r"/prj/aquacat/Data" #r'S:/Data' #
# CHANGE THIS TO THE TOP LEVEL OF THE FOLDER THE NETCDFs ARE IN
outlevel = toplevel 

# CHANGE THIS TO WHERE THE hasData files are, they should exist in the toplevel folder.
rn = pd.read_csv(f"{toplevel}/hasData_primary.csv")
rnreg = pd.read_csv(f"{toplevel}/hasData_Regions.csv")

method="OBS"  # DON'T CHANGE!!!

regional = False

if regional:
    subfold='/NW'
    fileinfix = 'NW_POT2_pc01'
    rn = rn[rnreg.REGION=="NW"]
    NH = 1437
else:
    subfold=''
    fileinfix = 'POT2_pc01'
    NH = 19914

ncpath0 = (f"{outlevel}/RCM{rcm}_{period}{subfold}/event{method}"
           f"_POT2_pc01_RCM{rcm}_{period}.nc")

ec_events = nc.Dataset(ncpath0, mode='r')
vvec_all = ec_events.variables['flow'][:,:]
avec_all = np.zeros(shape=vvec_all.shape)
dvec_all = np.zeros(shape=vvec_all.shape)

eventNo = list(ec_events.variables["eventNo"][:])
NE = np.sum([i > 0 for i in eventNo])

ec_events.close()

g2g_wd = "/prj/aquacat/run_hmfg2g/outputs/"
suffix_pres = "_198012_201011"
subfold_pres = f"RCM{rcm}{suffix_pres}/"
ncpres = f"{g2g_wd}dmflow_RCM{rcm}{suffix_pres}_out.nc"
ncin_pres = nc.Dataset(ncpres, mode='r')

param_path = (f"{outlevel}/RCM{rcm}_{period}/paramtableG"
              f"_POT2_RCM{rcm}_{period}.csv")
param_table = pd.read_csv(param_path)

thresh_path = f"{outlevel}/RCM{rcm}_{period}/threshMat_RCM{rcm}_{period}.csv"
threshvec = pd.read_csv(thresh_path).iloc[:,1]

print("Setup complete")
print("loop start")
start_time = time.time()
for h in range(NH):
  #print(h)
    if (h < 10) | (h % 1000 == 0): # time recording
        print(h)
        print("--- %s seconds ---" % (time.time() - start_time))
        start_time = time.time()
    vals = ncin_pres.variables["dmflow"][:,rn.col[h]-1, rn.row[h]-1]
    
    avec_all[:,h], dvec_all[:,h] = dpeApeCalculator(h=h,
                                              vals=vals,
                                              obs_events=list(vvec_all[:,h]),
                                              pars=param_table.iloc[h,2:5],
                                              thresh_val=threshvec[h])

print("Saving outputs")

ncin_pres.close()

ncpath = (f"{outlevel}/RCM{rcm}_{period}{subfold}/event{method}"
          f"_POT2_pc01_RCM{rcm}_{period}.nc")
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

event_var[:] = eventNo

flow_var[:,:] = vvec_all.T
ape_var[:,:] = avec_all.T
dpe_var[:,:] = dvec_all.T

ncfile.RCM = rcm
ncfile.period = period
ncfile.event_threshold = "POT2"
ncfile.area_lower_limit = "pc01"
ncfile.method=method

ncfile.close()

yaml_path = f"{outlevel}/RCM{rcm}_{period}/settings.yaml"
with open(yaml_path) as ym:
    list_doc = yaml.safe_load(ym)

list_doc['EC2ape'] = True
list_doc['EC2dpe'] = True

with open(yaml_path, 'w') as ym:
    yaml.dump(list_doc, ym)
    
print("done 117N")