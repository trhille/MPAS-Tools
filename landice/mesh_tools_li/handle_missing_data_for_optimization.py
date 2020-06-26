#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 10:47:39 2020
Handle missing data for Albany optimization. To be used after interpolate_to_mpasli_grid.py
@author: trevorhillebrand
"""
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import sys
import numpy as np
import os
from netCDF4 import Dataset
from optparse import OptionParser
import geopandas as gpd
import math
from collections import OrderedDict
from shapely.geometry import Point, Polygon
import scipy.spatial
import time
from datetime import datetime

print("== Gathering information.  (Invoke with --help for more details. All arguments are optional)\n")
parser = OptionParser()
parser.description = __doc__
parser.add_option("-f", "--file", dest="inputFile", help="name of source (input) file in MPASLI format.", default="landice_grid.nc", metavar="FILENAME")
parser.add_option("-g", "--geometric_features", dest="geoFeaturesPath", help="Path to geometric_features directory", default="landice_grid.nc", metavar="FILENAME")

for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

data = Dataset(options.inputFile, 'r+')

# load data
thk = data.variables["thickness"][:]

xSpeed = data.variables["observedSurfaceVelocityX"][:] #units of m/s
ySpeed = data.variables["observedSurfaceVelocityY"][:]
speed = np.sqrt(xSpeed**2 + ySpeed**2)

speedErr = data.variables["observedSurfaceVelocityUncertainty"][:]

dHdt = data.variables["observedThicknessTendency"][:] #units of m/s
dHdtErr = data.variables["observedThicknessTendencyUncertainty"][:]

smb = data.variables["sfcMassBal"][:] # units of kg / m^2 / s
smbErr = data.variables["sfcMassBalUncertainty"][:]

#
latCell = data.variables["latCell"][:] * 180/np.pi
lonCell = data.variables["lonCell"][:] * 180/np.pi
nCells = data.dimensions["nCells"].size

# find missing values. Missing values in ORNL dataset are large
missing_value = 100.
missing_xSpeed = (xSpeed > missing_value)
missing_ySpeed = (ySpeed > missing_value)
missing_speed = missing_xSpeed + missing_ySpeed

missing_dHdt = (dHdt > missing_value)

missing_smb = (smb > missing_value)

# set 0 values missing data in velocity, SMB, and  dHdt fields
xSpeed[missing_xSpeed] = 0.
ySpeed[missing_ySpeed] = 0.
smb[missing_smb] = 0.
dHdt[missing_dHdt] = 0.

# set large (1e20) uncertainty for missing data in velocity and dHdt fields
missing_uncertainty_value = 1e20
speedErr[missing_speed] = 1e20
dHdtErr[missing_dHdt] = 1e20
smbErr[thk==0] = 8e-5

# set dHdt error based on surface velocity
# 50 m/yr speed threshold from Csatho et al. (2014) for dHdt uncertainties:
# <50 m/yr --> 0.01 m/yr dHdt uncertainty; >=50 m/yr --> 50 cm/yr dHdt uncertainty
s_per_year = 60. * 60. * 24. * 365.
speedThreshold = 50. / s_per_year # convert to m/s
dHdtErr_slow = 0.01 / s_per_year
dHdtErr_fast = 0.5 /s_per_year

slow = (speed<speedThreshold)
fast = (speed>=speedThreshold)
dHdtErr[slow] = dHdtErr_slow
dHdtErr[fast] = dHdtErr_fast



# set dHdt error for fast-flowing by basin 1 m·y−1for northwest and 
# southwest drainage basins; and 2 m·y−1for Jak (westCentralGreenland) and southeast drainage basins
#geojsonPath = options.geoFeaturesPath + '/geometric_data/landice/region/'
#

# loop through directories and read data from Greenland regions
#for directory in os.listdir(options.geoFeaturesPath):
#    if 'westCentralGreenland' in directory or 'southEastGreenland' in directory \
#        or 'northWestGreenland' in directory or 'southWestGreenland' in directory:
#       
#        basin = gpd.read_file(options.geoFeaturesPath + directory + '/region.geojson')
#        poly = basin['geometry'].iloc[0]
#        print(str(basin['name']))
#        nCellsInBasin = 0
#
#        for cell in np.arange(0, nCells):
#            
#            pt = Point(lonCell[cell], latCell[cell])
#            if cell % (nCells // 10) == 0.:
#                print("Cell {} of {}".format(cell, nCells))
#            if pt.within(poly):
#                nCellsInBasin += 1
#                if 'westCentralGreenland' in directory or 'southEastGreenland' in \
#                    directory:
#                        dHdtErr[0, cell] = 2.0 / s_per_year
#                elif 'northWestGreenland'  in directory or 'southWestGreenland' in \
#                    directory:
#                        dHdtErr[0, cell] = 1.0 / s_per_year
#        print('Found {} cells in basin {}'.format(nCellsInBasin, directory))
            
dHdtErr[slow] = dHdtErr_slow            
            
data.variables["observedSurfaceVelocityX"][:]= xSpeed #units of m/s
data.variables["observedSurfaceVelocityY"][:] = ySpeed

data.variables["observedSurfaceVelocityUncertainty"][:] = speedErr

data.variables["observedThicknessTendency"][:] = dHdt
data.variables["observedThicknessTendencyUncertainty"][:] = dHdtErr

data.variables["sfcMassBal"][:] = smb
data.variables["sfcMassBalUncertainty"][:] = smbErr

            
data.close()
        
