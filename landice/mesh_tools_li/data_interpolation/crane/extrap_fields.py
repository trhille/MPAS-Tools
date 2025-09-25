#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 17:23:11 2021
This script uses nearest-neighbor extrapolation to fill in ice-free areas in 
bedmachine dataset for interpolation to MPAS mesh. 
First 
@author: trevorhillebrand
"""
from argparse import ArgumentParser
from scipy.interpolate import NearestNDInterpolator
from netCDF4 import Dataset
import numpy as np
import time
from shutil import copyfile


print("Gathering Information ... ")
parser = ArgumentParser(
        prog='extrap_Measures.py',
        description="Extrapolates MEaSUREs velocity data into regions of missing data.")
parser.add_argument("-i", dest="dataFile",
                    help="Original data file", metavar="FILENAME")
parser.add_argument("-o", dest="outFile",
                    help="Output file with extrapolated fields", metavar="FILENAME")
parser.add_argument("-v", dest="variables",
                    help="Space-separated list of variables to extrapolate",
                    metavar="VARNAME", nargs='+',)
args = parser.parse_args()

copyfile(args.dataFile, args.outFile)
data = Dataset(args.outFile, 'r+')
x1 = data.variables["x"][:]
y1 = data.variables["y"][:]
x = np.arange(0,  data.dimensions["x"].size)
y = np.arange(0,  data.dimensions["y"].size)

tic = time.perf_counter()
print('Constructing meshgrid')
xGrid, yGrid = np.meshgrid(x1,y1)
xx = xGrid.ravel()
yy = yGrid.ravel()
toc = time.perf_counter()
print(f'Constructed meshgrid in {toc - tic:.2f} seconds')
# Get data and masks for each variable.
# Remember that in this case, a mask value
# of True means that the value is missing there.
for var in args.variables:
    var_array = data.variables[var][:].data.ravel()
    var_mask = np.logical_not(data.variables[var][:].mask.ravel())

    tic = time.perf_counter()
    print(f'Beginning building {var} interpolator')
    interp_var = NearestNDInterpolator(list(zip(xx[var_mask], yy[var_mask])), var_array[var_mask])
    toc = time.perf_counter()
    print(f'Finished building {var} interpolator in {toc - tic:.2f} seconds')

    tic = time.perf_counter()
    print(f'Beginning {var} extrapolation')
    extrap_var = interp_var(xGrid,yGrid)
    toc = time.perf_counter()
    print(f'{var} extrapolation completed in {toc - tic:.2f} seconds')
    data.variables[var][:] = extrap_var
    data.sync()
    del extrap_var

data.close()
