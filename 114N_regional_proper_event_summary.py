# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 11:17:32 2021
Author: Adam Griffin, UKCEH
Project: AQUACAT
Script to summarise EC events.
"""
import os
import netCDF4 as nc
import pandas as pd
import numpy as np
import re
import sys
import yaml

rcm = sys.argv[1]

if sys.argv[2] == "present":
    period = "198012_201011"
else:
    period = "205012_208011"

print("Running RCM {} for {}.".format(rcm, period))

RCMS = ["01","04","05","06","07","08","09","10","11","12","13","15"]
periods = ["198012_201011","205012_208011"]

# CHANGE THIS TO THE TOP LEVEL OF THE FOLDER THE CSVs ARE IN
toplevel = r"/prj/aquacat/Data"

# CHANGE THIS TO THE TOP LEVEL OF THE FOLDER THE NETCDFs ARE IN
outlevel = toplevel #r'S:/Data'

# CHANGE THIS TO WHERE THE hasData files are, they should exist in the toplevel folder.
rn = pd.read_csv("".join([toplevel,"/hasData_primary.csv"]))
rnreg = pd.read_csv("".join([toplevel,"/hasData_Regions.csv"]))

regional = True

if regional:
    subfold='/NW'
    fileinfix = 'NW_POT2_pc01'
    rn = rn[rnreg.REGION=="NW"]
    NH = 1437
else:
    subfold=''
    fileinfix = 'POT2_pc01'
    NH = 19914

#for method in ["OBS", "EC", "HT"]:
method="OBS" # "EC" OR "OBS", "HT" AS NEEDED

ncpath = (f"{outlevel}/RCM{rcm}_{period}{subfold}/event{method}_{fileinfix}"
          f"_RCM{rcm}_{period}.nc")
ncfile = nc.Dataset(ncpath, mode='r')

param_path = (f"{outlevel}/RCM{rcm}_{period}{subfold}/paramtableG"
              f"_POT2_RCM{rcm}_{period}.csv")

param_table = pd.read_csv(param_path)

thresh_path = f"{outlevel}/RCM{rcm}_{period}{subfold}/threshMat_RCM{rcm}_{period}.csv"

threshvec = pd.read_csv(thresh_path).iloc[:,1]
if regional:
    threshvec = threshvec[rnreg.REGION == "NW"]

init_path = (f"{outlevel}/RCM{rcm}_{period}/initialSummary_RCM{rcm}_{period}.csv")
init_table = pd.read_csv(init_path)

summ_path = (f"{outlevel}/RCM{rcm}_{period}{subfold}/eventSumm"
             f"_OBS_POT2_pc01_RCM{rcm}_{period}.csv")


threshvec = pd.read_csv(thresh_path).iloc[:,1]
if regional:
    threshvec = threshvec[rnreg.REGION == "NW"]

summtable_out = pd.DataFrame(columns=["eventNumber", "eventDay", "eventLength",
                                      "area","peakA", "peakD", "season",
                                      "nclusters","peakyness"])

eventNo = list(ncfile.variables["eventNo"][:])
NE = np.sum([i > 0 for i in eventNo])

avec_all = ncfile.variables['ape'][:,:]
avec_mid = ncfile.variables['ape_mid'][:,:]
vvec_all = ncfile.variables['flow'][:,:]
dvec_all = ncfile.variables['dpe'][:,:]

ncfile.close()

print("Setup complete")

import time

start_time = time.time()
for i in range(NE):
    if (i < 10) or (i % 200) == 0:
        print(i)
        print("--- %s seconds ---" % (time.time() - start_time))
        start_time = time.time()
    ni = eventNo[i]
    vvec = sum(vvec_all[i,:] > threshvec)
    avec = min(avec_all[i,:])
    amid = min(avec_mid[i,:])
    dvec = min(dvec_all[i,:])
    D = summtable.iloc[ni-1,1:3]
    summtable_out.loc[i] = [ni, D[0], D[1], vvec, avec, dvec, summtable.iloc[ni-1,6], 0, 0]

print("--- %s seconds ---" % (time.time() - start_time))
print("done")

yaml_path = "".join([outlevel, '/RCM', rcm, '_', period, '/settings.yaml'])

summtable_out.to_csv("".join([outlevel, r'/RCM', rcm, '_', period, subfold,
                  r'/eventSumm_', method, 'B_', fileinfix,'_RCM',rcm,
                  '_',period, '.csv'])) 

with open(yaml_path) as ym:
    list_doc = yaml.safe_load(ym)

list_doc['HTsumm'] = True
list_doc['propsumm'] = "114N.py"

with open(yaml_path, 'w') as ym:
    yaml.dump(list_doc, ym)
    
    print(time.strftime("%Y-%m-%d %H:%M:%S"))
print("Files saved and YAML updated. End.")