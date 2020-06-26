#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Trevor Hillebrand
"""
from netCDF4 import Dataset
import numpy as np
import h5py

filepath = '/Users/trevorhillebrand/Documents/Greenland/tmp.nc'

data = Dataset(filepath, 'r+')
data.set_auto_mask(False)

# load in velocity and uncertainty fields
x = data.variables["x1"][:]
y = data.variables["y1"][:]
vx = data.variables["vx"][0, :, :]
vy = data.variables["vy"][0, :, :]
vxErr = data.variables["ex"][0, :, :]
vyErr = data.variables["ey"][0, :, :]

# velocity and uncertainty fields have very large missing_values. replace those with nans
vx[vx == data.variables["vx"].missing_value] = np.nan # x-component of surface velocity
vy[vy == data.variables["vy"].missing_value] = np.nan # y-component of surface velocity
vxErr[vxErr == data.variables["ex"].missing_value] = np.nan # x-component surface velocity error
vyErr[vyErr == data.variables["ey"].missing_value] = np.nan # y-component surface velocity error

surfaceSpeed = np.sqrt(vx**2 + vy**2)
surfaceSpeed[surfaceSpeed == 0] = np.nan #enable division by surfaceSpeed when value is zero

vErr = 0.5 * np.sqrt( (2 * vxErr * vx)**2 + (2 * vyErr * vy)**2 ) / surfaceSpeed
vErr[np.isnan(vErr)] = data.variables["ex"].missing_value

#assert "vErr" not in data.variables.keys()
#
#data.createVariable("vErr", 'f4', dimensions=("time", "y1", "x1"))
    
data.variables["vErr"][0, :, :] = vErr

data.close()