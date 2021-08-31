#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 15:16:04 2020
Plot  bed and surface topography from MPAS netCDF along a transect defined in csv
@author: trevorhillebrand
"""

import numpy as np
import csv
from netCDF4 import Dataset
from optparse import OptionParser
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

parser = OptionParser(description='Plot transect from MPAS netCDF')
parser.add_option("-d", "--data", dest="data_file", help="the MPAS netCDF file")
parser.add_option("-c", "--coords", dest="coords_file", help="csv file defining transect")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

dataset = Dataset(options.data_file, 'r')
dataset.set_always_mask(False)
xCell = dataset.variables["xCell"][:]
yCell = dataset.variables["yCell"][:]
areaCell = dataset.variables["areaCell"][:]
# only take thickness and speed of dynamic ice
if "cellMask_dynamicIce" in dataset.variables.keys():
    thk = dataset.variables["thickness"][:] * dataset.variables["cellMask_dynamicIce"][:]
    speed = dataset.variables["surfaceSpeed"][:] * dataset.variables["cellMask_dynamicIce"][:] * 3600. * 24. * 365.
else:
    print('Variable cellMask_dynamicIce was not found. Speed interpolation may be incorrect.' +
          'Please run convert_landice_bitmasks.py and try again')
    thk = dataset.variables["thickness"][:] * ( dataset.variables["cellMask"][:] & 2 )//2
    speed = dataset.variables["surfaceSpeed"][:] * ( dataset.variables["cellMask"][:] & 2 )//2 * 3600. * 24. * 365.
    
bedTopo = dataset.variables["bedTopography"][0,:]
startYear = 2000
time = dataset.variables["daysSinceStart"][:] / 365. + startYear
time=time[0:30]
thk = thk*(thk<5e3)

coords = []
x = []
y = []

with open(options.coords_file, newline='') as csvfile:
     reader = csv.reader(csvfile, delimiter=',')
     
     for row in reader:
         coords.append(row)
         x.append(float(row[0]))
         y.append(float(row[1]))

xArray = np.array(x) * 1000.0
yArray= np.array(y) * 1000.0

d_distance = np.zeros(len(xArray))
for ii in np.arange(1, len(xArray)):
    d_distance[ii] = np.sqrt( (xArray[ii] - xArray[ii-1])**2 + (yArray[ii] - yArray[ii-1])**2 )

distance = np.cumsum(d_distance)

transectFig, transectAx = plt.subplots(2,1)
thickAx = transectAx[0]
thickAx.grid()
speedAx = transectAx[1]
speedAx.grid()
timeColors = cm.plasma(np.linspace(0,1,len(time)))

plt.rcParams.update({'font.size': 16})

bed_interpolant = LinearNDInterpolator(np.vstack((xCell, yCell)).T, bedTopo)
bed_transect = bed_interpolant(np.vstack((xArray, yArray)).T)

for timeSlice in np.arange(0, len(time)):

    thk_interpolant = LinearNDInterpolator(np.vstack((xCell, yCell)).T, thk[timeSlice,:])
    thk_transect = thk_interpolant(np.vstack((xArray, yArray)).T)
    thk_transect[(thk_transect + bed_transect) < 1.] = 0.0
    thickAx.plot( (distance[-1]-distance)/1000., thk_transect + bed_transect, color=timeColors[timeSlice])
    
    speed_interpolant = LinearNDInterpolator(np.vstack((xCell, yCell)).T, speed[timeSlice,:])
    speed_transect = speed_interpolant(np.vstack((xArray, yArray)).T)

    speed_transect[speed_transect==0] = np.nan
    speed_transect[thk_transect == 0.] = np.nan
    speedAx.plot( (distance[-1]-distance)/1000. , speed_transect, color=timeColors[timeSlice])
    
thickAx.plot( (distance[-1]-distance)/1000., bed_transect, color='black')

speedAx.set_xlabel('Distance (km)')
speedAx.set_ylabel('Surface\nspeed (m/yr)')
speedAx.set_xlim((0,35))
thickAx.set_ylabel('Elevation\n(m asl)')
thickAx.set_xlim((0,35))
transectFig.subplots_adjust(hspace=0.4)

cbar = plt.colorbar(cm.ScalarMappable(cmap='plasma'), ax=transectAx)
cbar.set_ticks(np.linspace(0,1,6))
cbar.set_label('Year', size='large')
years = str(time).lstrip('[').rstrip(']').split()
cbar.set_ticklabels(years[0::len(years)//5])
cbar.ax.tick_params(labelsize=12)

plt.show()

dataset.close()