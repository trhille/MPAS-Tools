#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 13:14:43 2020
Removes floating ice from a mpas mesh file
@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np

input_path = ('/Users/trevorhillebrand/Documents/mpas/TESTS/landice/greenland/Humboldt_3to30km/standard_configuration/create_mesh/Humboldt_3to30km.nc')

data = Dataset(input_path, 'r+')
data.set_auto_mask(False)

bed = data.variables["bedTopography"][0, :]
thk = data.variables["thickness"][0, :]
x = data.variables["xCell"][:]

floating_ice_mask = ((thk*910/1028+bed)<0)*(thk>0)

Humboldt_boundary_x = -287873.
Humboldt_x = (x < Humboldt_boundary_x)

thk[Humboldt_x] = thk[Humboldt_x] * (1-floating_ice_mask)[Humboldt_x]

data.variables["thickness"][0, :] = thk

data.close()