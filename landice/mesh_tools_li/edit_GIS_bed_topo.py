#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 14:32:20 2020

@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np

sourceFile = '/Users/trevorhillebrand/Documents/Greenland/greenland_1km_2020_04_20.epsg3413.icesheetonly.nc'
destFile = '/Users/trevorhillebrand/Documents/Greenland/greenland_1km_2020_04_20.epsg3413.icesheetonly.nobathymetry.nc'

sourceData = Dataset(sourceFile, 'r')

destData = Dataset(destFile, 'w', format="NETCDF4")

x = sourceData.variables["x1"][:]
y = sourceData.variables["y1"][:]
thk = sourceData.variables["thk"][0,:]
bedTopo = sourceData.variables["topg"][0,:]

xGrid, yGrid = np.meshgrid(x, y)

xmin = -4.5e5 #lefthand x extent of area to cut out
xmax = -3.5e5 #righthand x extent of area to cut out
ymin = -1.125e6 #bottom y of area to cut out
ymax = -1.015e6 #top y " "

# The mask here has many conditions, but np.logical_and only takes two arguments
mask = np.logical_and(np.logical_or(thk<50., bedTopo<0), xmin<xGrid)
mask = np.logical_and(mask, xGrid<xmax)
mask = np.logical_and(mask, ymin<yGrid)
mask = np.logical_and(mask, yGrid<ymax)
mask = mask==False

bedTopo = bedTopo[mask]


xCell = xGrid[mask]
yCell = yGrid[mask]

destData.createDimension('Time', None)
destData.createDimension('nCells', len(xCell))
destData.createVariable('xCell', 'f', ('nCells'))
destData.createVariable('yCell', 'f', ('nCells'))
destData.createVariable('bedTopography', 'f', ('Time', 'nCells'))

destData.variables["xCell"][:] = xCell
destData.variables["yCell"][:] = yCell
destData.variables["bedTopography"][0,:] = bedTopo

sourceData.close()
destData.close()