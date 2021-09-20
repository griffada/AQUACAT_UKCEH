# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 10:08:32 2021

@author: adagri
"""

from netCDF4 import Dataset
import numpy as np
import pandas as pd

dpe = pd.read_csv("S:/Data/RCM01_198012_201011/eventdpe_EC_POT2_pc01_RCM01_198012_201011.csv")
ape = pd.read_csv("S:/Data/RCM01_198012_201011/eventdpe_EC_POT2_pc01_RCM01_198012_201011.csv")
flow = pd.read_csv("S:/Data/RCM01_198012_201011/eventflow_EC_POT2_pc01_RCM01_198012_201011.csv")

ncfile = Dataset("S:/Data/RCM01_198012_201011/nctest2.nc", mode='w')
loc_dim = ncfile.createDimension('loc', 19914)
event_dim = ncfile.createDimension('event', None)
dpe_var = ncfile.createVariable('dpe', np.float32, ('event','loc'))
dpe_var.units = 'PoE'
dpe_var.standard_name = "Daily_PoE"
dpe_var[:,:] = dpe.iloc[:,4:].T

ape_var = ncfile.createVariable('ape', np.float32, ('event','loc'))
ape_var.units = 'PoE'
ape_var.standard_name = "Annual_PoE"
ape_var[:,:] = ape.iloc[:,4:].T

flow_var = ncfile.createVariable('flow', np.float32, ('event','loc'))
flow_var.units = 'cumec'
flow_var.standard_name = "Peak_Flow"
flow_var[:,:] = flow.iloc[:,4:].T

row_var = ncfile.createVariable('row', np.int32, ('loc'))
row_var[:] = ape.iloc[:,0]

col_var = ncfile.createVariable('col', np.int32, ('loc'))
col_var[:] = ape.iloc[:,1]

north_var = ncfile.createVariable('northing', np.int32, ('loc'))
north_var.units = 'm'
north_var.standard_name = 'Northing'
north_var[:] = ape.iloc[:,2]

east_var = ncfile.createVariable('easting', np.int32, ('loc'))
east_var.units = 'm'
east_var.standard_name = 'Easting'
east_var[:] = ape.iloc[:,3]

ncfile.rcm = "RCM01"
ncfile.period = "present"
ncfile.threshold = "POT2"
ncfile.inundation = "pc01"

print(ncfile)
ncfile.close()
