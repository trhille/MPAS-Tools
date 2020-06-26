#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 14:32:20 2020

@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np

sourceFile = '/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_1to10km_betaonly_041720.nc'
destFile = '/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_1to10km_betaonly_041720_bedTopoTrimmed.nc'

sourceData = Dataset(sourceFile, 'r')
sourceData.set_auto_mask(False)

destData = Dataset(destFile, 'w', format="NETCDF4")

xCell = sourceData.variables["xCell"][:]
yCell = sourceData.variables["yCell"][:]
thk = sourceData.variables["thickness"][0,:]
bedTopo = sourceData.variables["bedTopography"][0,:]

# Get domain boundary cells. This snippet taken from mark_domain_boundaries_dirichlet.py
nCells = len(sourceData.dimensions['nCells'])
cONc = sourceData.variables['cellsOnCell'][:]
nEdgesOnCell = sourceData.variables['nEdgesOnCell'][:]

mask = bedTopo*False
for i in range(nCells):
   nE = nEdgesOnCell[i]
   if min(cONc[i, :nE]) == 0:
      mask[i] = 1

mask =  (mask + (thk > 200.) + (bedTopo > 0)) >= 1

bedTopoNew = bedTopo[mask]
xCellNew = xCell[mask]
yCellNew = yCell[mask]

destData.createDimension('nCells', len(xCellNew))
destData.createDimension('Time', None)
destData.createVariable('xCell', 'f', ('nCells'))
destData.createVariable('yCell', 'f', ('nCells'))
destData.createVariable('bedTopography', 'f', ('Time', 'nCells'))

destData.variables["xCell"][:] = xCellNew
destData.variables["yCell"][:] = yCellNew
destData.variables["bedTopography"][0,:] = bedTopoNew

sourceData.close()
destData.close()