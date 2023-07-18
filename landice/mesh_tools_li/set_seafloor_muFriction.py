#!/usr/bin/env python
"""
Script to extrapolate arbitrary variable

Created on Mon Feb1 2021

@author: Matt Hoffman, Trevor Hillebrand
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import sys
from datetime import datetime


parser = OptionParser(description='')

parser.add_option("-f", "--file", dest="nc_file", 
                  help="the mpas file to write to")
parser.add_option("-r", "--regions_file", dest="regions_file", 
                  help=".nc file containing region masks")

parser.add_option("-p", "--percentile", dest="percentile", default=None, 
                  help="percentile to define regional seafloor muFriction")
parser.add_option("-s", "--speed_threshold", dest="speed_threshold", 
                  default=None, 
                  help="grounded ice speed above which values are used to" +
                  " define muFriction. units of m/yr")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()

dataset = Dataset(options.nc_file, 'r+')
regions = Dataset(options.regions_file, 'r')
speed_threshold = float(options.speed_threshold)
percentile = float(options.percentile)
muFriction = dataset.variables['muFriction'][0, :]
muFrictionNew = muFriction.copy()
region_masks = regions.variables['regionCellMasks'][:]
# Extrapolation
nCells = len(dataset.dimensions['nCells'])



thickness = dataset.variables['thickness'][0,:]
bed = dataset.variables["bedTopography"][0,:]
cellsOnCell = dataset.variables['cellsOnCell'][:]
nEdgesOnCell = dataset.variables['nEdgesOnCell'][:]
xCell = dataset.variables["yCell"][:]
yCell = dataset.variables["xCell"][:]
spd = np.sqrt(dataset.variables["uReconstructX"][0, :, 0]**2 + 
              dataset.variables["uReconstructY"][0, :, 0]**2)


# Define region of good data to extrapolate from.  Different methods for different variables
groundedMask = (thickness > (-1028.0 / 910.0 * bed))
keepCellMask = np.copy(groundedMask)

for iCell in range(nCells):
    for n in range(nEdgesOnCell[iCell]):
        jCell = cellsOnCell[iCell, n] - 1
        if (groundedMask[jCell] == 1):
            keepCellMask[iCell] = 1
            continue
keepCellMask *= (muFriction > 0)  # ensure zero muFriction does not get extrapolated

keepCellMaskNew = np.copy(keepCellMask)  # make a copy to edit that will be used later
keepCellMaskOrig = np.copy(keepCellMask)  # make a copy to edit that can be edited without changing the original

# loop over regions and determine appropriate muFriction value
for region in np.arange(regions.dimensions['nRegions'].size):
    region_mask = region_masks[:, region] == 1
    speed_mask = (spd * region_mask * keepCellMask * 3.154e7) > speed_threshold
    seafloor_mu = np.quantile(muFriction[speed_mask], 0.05)
    replace_mu_mask = np.where( np.logical_and(
                                    np.logical_and((muFriction > seafloor_mu), 
                                    (keepCellMask == 0)), region_mask))
    muFrictionNew[replace_mu_mask] = seafloor_mu 

dataset.variables['muFriction'][0,:] = muFrictionNew # Put updated array back into file.
# === Clean-up =============

# Update history attribute of netCDF file
thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + " ".join(sys.argv[:])
if hasattr(dataset, 'history'):
   newhist = '\n'.join([thiscommand, getattr(dataset, 'history')])
else:
   newhist = thiscommand
setattr(dataset, 'history', newhist)


dataset.close()
