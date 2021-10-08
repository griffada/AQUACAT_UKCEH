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
import time
import gc

def season(x):
  return ["DJF","MAM","JJA","SON"][((x // 90) % 4)]

#rcm = "01"
#period="205012_208011"

rcm = sys.argv[1]

if sys.argv[2] == "present":
    period = "198012_201011"
else:
    period = "205012_208011"

print(f"Running RCM {rcm} for {period}.")

RCMS = ["01","04","05","06","07","08","09","10","11","12","13","15"]
periods = ["198012_201011","205012_208011"]

# CHANGE THIS TO THE TOP LEVEL OF THE FOLDER THE CSVs ARE IN
toplevel = r"/prj/aquacat/Data" #r'S:/Data'

# CHANGE THIS TO THE TOP LEVEL OF THE FOLDER THE NETCDFs ARE IN
outlevel = toplevel #r'S:/Data'

# CHANGE THIS TO WHERE THE hasData files are, they should exist in the toplevel folder.
rn = pd.read_csv(f"{toplevel}/hasData_primary.csv")
rnreg = pd.read_csv(f"{toplevel}/hasData_Regions.csv")

method="EC"

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

if len(sys.argv) > 3:
    if sys.argv[3] == "FF":
        subfold = '_FF'

ncpath = (f"{outlevel}/RCM{rcm}_{period}{subfold}/event{method}"
          f"_POT2_pc01_RCM{rcm}_{period}.nc")
ncfile = nc.Dataset(ncpath, mode='r')

param_path = (f"{outlevel}/RCM{rcm}_{period}{subfold}/paramtableG_POT2_RCM{rcm}"
              f"_{period}.csv")

param_table = pd.read_csv(param_path)

thresh_path = f"{outlevel}/RCM{rcm}_{period}{subfold}/threshMat_RCM{rcm}_{period}.csv"

threshvec = pd.read_csv(thresh_path).iloc[:,1]
if regional:
    threshvec = threshvec[rnreg.REGION == "NW"]

summ_path = (f"{outlevel}/RCM{rcm}_{period}{subfold}/eventSumm_OBS"
             f"_POT2_pc01_RCM{rcm}_{period}.csv")
summtable = pd.read_csv(summ_path)

init_path = (f"{outlevel}/RCM{rcm}_{period}{subfold}/initialSummary"
             f"_RCM{rcm}_{period}.csv")
init_table = pd.read_csv(init_path)

summtable_out = pd.DataFrame(columns=["eventNumber", "eventDay", "eventLength",
                                      "area","peakA", "peakA_mid", "peakD", "season",
                                      "nclusters","peakyness"])

eventNo = list(ncfile.variables["eventNo"][:])
NE = np.sum([i > 0 for i in eventNo])

avec_all = ncfile.variables['ape'][:,:]
avec_mid = ncfile.variables['ape_mid'][:,:]

print("Setup complete")

start_time = time.time()
for i in range(NE):
    if (i < 10) or (i % 1000) == 0:
        print(i)
        print("--- %s seconds ---" % (time.time() - start_time))
        start_time = time.time()
    ni = eventNo[i]
    vvec = 0
    avec = min(avec_all[i,:])
    amid = min(avec_mid[i,:])
    dvec = 0
    D = init_table.iloc[ni-1,:] # Done in R afterwards
    #seas = init_table.iloc[ni-1, 3]) # Done in R afterwards
    summtable_out.loc[i] = [ni, D[0], D[1],
                            vvec, avec, amid, dvec,
                            D[3], 0, 0]

print("--- %s seconds ---" % (time.time() - start_time))
print("avec amid done")

del avec_all
del avec_mid
gc.collect()

vvec_all = ncfile.variables['flow'][:,:]
dvec_all = ncfile.variables['dpe'][:,:]

for i in range(NE):
    if (i < 10) or (i % 1000) == 0:
        print(i)
        print(">>> %s seconds >>>" % (time.time() - start_time))
        start_time = time.time()
    ni = eventNo[i]
    summtable_out.iloc[i,3] = sum(vvec_all[i,:] > threshvec)
    summtable_out.iloc[i,6] = min(dvec_all[i,:])



print(">>> %s seconds >>>" % (time.time() - start_time))
print("vvec dvec done")


ncfile.close()

yaml_path = f"{outlevel}/RCM{rcm}_{period}/settings.yaml"

summpath_out = (f"{outlevel}/RCM{rcm}_{period}{subfold}/eventSumm_"
                f"{method}_POT2_pc01_RCM{rcm}_{period}.csv")

summtable_out.to_csv(summpath_out, index=False) 

with open(yaml_path) as ym:
    list_doc = yaml.safe_load(ym)

list_doc['EC2summ'] = True
list_doc['ECpropsumm'] = "113bN.py"

with open(yaml_path, 'w') as ym:
    yaml.dump(list_doc, ym)

print(time.strftime("%Y-%m-%d %H:%M:%S"))
print("Files saved and YAML updated. End.")
