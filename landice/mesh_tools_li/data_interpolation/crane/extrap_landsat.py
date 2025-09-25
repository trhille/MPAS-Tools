#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 17:23:11 2021
This script uses nearest-neighbor extrapolation to fill in ice-free areas in 
bedmachine dataset for interpolation to MPAS mesh. 
First 
@author: trevorhillebrand
"""
from scipy.interpolate import NearestNDInterpolator
from netCDF4 import Dataset
import numpy as np
import time
from shutil import copyfile


dataFile = 'landsat_18dec2002_20feb2003_hsc15_htc70_inc2_3031_new_hp.nc'
outFile = copyfile(dataFile, dataFile + '_extrap')
data = Dataset(outFile, 'r+')
# Get data and masks for each variable.
# Remember that in this case, a mask value
# of True means that the value is missing there.
vx = data.variables["vx_masked"][:].data.ravel()
vx_mask = np.logical_not(data.variables["vx_masked"][:].mask.ravel())
vy = data.variables["vy_masked"][:].data.ravel()
vy_mask = np.logical_not(data.variables["vy_masked"][:].mask.ravel())
# Mask needs an extra speed threshold of ~1500 m/yr to remove spurious fast cells.
spd_thresh = 4.  # m/day
mask = np.logical_and(np.logical_and(vx_mask, vy_mask), np.sqrt(vx**2. + vy**2.) < spd_thresh)

x1 = data.variables["x1"][:]
y1 = data.variables["y1"][:]
x = np.arange(0,  data.dimensions["x"].size)
y = np.arange(0,  data.dimensions["y"].size)

tic = time.perf_counter()
print('Constructing meshgrid')
xGrid, yGrid = np.meshgrid(x1,y1)
xx = xGrid.ravel()
yy = yGrid.ravel()
toc = time.perf_counter()
print('Constructed meshgrid in {:.2f} seconds'.format(toc-tic))

tic = time.perf_counter()
print('Beginning building vx interpolator')
interp_vx = NearestNDInterpolator(list(zip(xx[mask], yy[mask])), vx[mask])

print('Beginning building vy interpolator')
interp_vy = NearestNDInterpolator(list(zip(xx[mask], yy[mask])), vy[mask])

toc = time.perf_counter()
print('Finished building interpolators in {:.2f} seconds'.format(toc - tic))

tic = time.perf_counter()
print('Beginning vx extrapolation')
vx_extrap = interp_vx(xGrid,yGrid)
toc = time.perf_counter()
print('Extrapolation completed in {:.2f} seconds'.format(toc - tic))
data.variables["vx"][:] = vx_extrap
data.sync()
del vx_extrap

tic = time.perf_counter()
print('Beginning vy extrapolation')
vy_extrap = interp_vy(xGrid, yGrid)
toc = time.perf_counter()
print('Extrapolation completed in {:.2f} seconds'.format(toc - tic))
data.variables["vy"][:] = vy_extrap
data.sync()
del vy_extrap

data.createVariable("mask", datatype="f", dimensions=("y","x"))
data.variables["mask"][:] = np.reshape(mask, (np.shape(y1)[0], np.shape(x1)[0]))
data.sync()

data.close()
