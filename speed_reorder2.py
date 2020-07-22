# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 14:17:25 2020

@author: adagri
"""
#%%
from netCDF4 import Dataset
import numpy as np
import time
from datetime import datetime
nc_trial = "/prj/aquacat/CodeABG/InterimData/dmflow_trial.nc"

ncin = Dataset(nc_trial, 'r', format='NETCDF4')

# check netCDF file format
print(ncin.file_format)

# get axis data
print(ncin.dimensions.keys())
print(ncin.dimensions['Time'])
tin = ncin.variables['Time']
northing = ncin.variables['Northing']
easting = ncin.variables['Easting']

# get length of axis data
ntime = len(tin)
nn = len(northing)
ne = len(easting)
nt = ncin.variables['dmflow'].shape[0]

#flow_wide = ncin.variables['dmflow'][10000,:,:]
# print axis
#print(tin[:])
#print(latitude[:])
#print(longitude[:])

# get variables
#print(ncin.variables.keys())
#print(ncin.variables['t2m'])

# read data
#vin = ncin.variables['dmflow'][0:1000,:,:]
#vin2 = ncin.variables['dmflow'][0:1000,100,211]
#print(vin.long_name)
#print(vin.units)

#------------------
# write netCDF file
#------------------

# open a netCDF file to write
ncout = Dataset('/prj/aquacat/CodeABG/InterimData/trial2.nc', 'w', format='NETCDF4')

nt=10800
nn=1000
ne=700

# define axis size
ncout.createDimension('Time', nt)  # unlimited
ncout.createDimension('Northing', nn)
ncout.createDimension('Easting', ne)

# create time axis
thyme = ncout.createVariable('Time', np.float32, ('Time'))
thyme.long_name = 'Time'
thyme.units = 'days since 1961-1-1 0:0:0'
thyme.calendar = 'standard'
thyme.axis = 'T'

# create latitude axis
nor = ncout.createVariable('Northing', np.float32, ('Northing'))
nor.standard_name = 'Northing'
nor.long_name = 'Northing'
nor.units = 'GB National Grid'
nor.axis = 'Y'

# create longitude axis
eas = ncout.createVariable('Easting', np.float32, ('Easting'))
eas.standard_name = 'Easting'
eas.long_name = 'Easting'
eas.units = 'GB National Grid'
eas.axis = 'X'

# create variable array
vout = ncout.createVariable('dmflow', np.float32,
                            ('Northing', 'Easting', 'Time'))
vout.long_name = 'Daily mean river flow'
vout.standard_name = 'dmflow'
vout.units = 'm3 s-1'

# copy axis from original dataset
thyme[:] = tin[0:nt]
nor[:] = northing[:]
eas[:] = easting[:]
A = [i for i in range(0,11000,500)] + [10800]
#%%
print("The big reordering.")
print(datetime.now())


for i in range(len(A)-1):
    vout[:,:,A[i]:A[i+1]] = np.moveaxis(
        ncin.variables['dmflow'][A[i]:A[i+1],:,:],0,2)
    print(i)

# close files
    
ncin.close()
ncout.close()


tic = time.perf_counter()
ncini = ncin.variables['dmflow'][0:500,100,211]
toc = time.perf_counter()
print(f"IN read: {toc - tic:0.4f} seconds")

tic = time.perf_counter()
ncouto = ncout.variables['dmflow'][100,211, 0:500]
toc = time.perf_counter()
print(f"OUT read: {toc - tic:0.4f} seconds")
