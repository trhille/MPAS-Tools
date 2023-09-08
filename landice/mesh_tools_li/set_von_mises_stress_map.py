#!/usr/bin/env python
'''
Simple script for creating an input file with region-by-region values for
the von Mises stress threshold.  This can be used to assign optimal
regional values identified through a tuning process.

A region mask file is required as input.  Values to assign are hardcoded
below.  There is not error checking that the number of values match the
number of regions, so use with care.

The script outputs a file called von_mises_calving_parameters.nc with the
assigned regional values.

Matt Hoffman, 9/19/2022
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt

parser = OptionParser(description=__doc__)
parser.add_option("-n", dest="fileRegions", help="region file name.", metavar="FILENAME")
options, args = parser.parse_args()

f = Dataset(options.fileRegions, 'r')
regionCellMasks = f.variables['regionCellMasks'][:]
nRegions = len(f.dimensions['nRegions'])
nCells = len(f.dimensions['nCells'])

fout = Dataset("von_mises_calving_parameters.nc", 'w')
fout.createDimension('nCells', nCells)
fout.createDimension('Time', None)
grdVM = fout.createVariable('groundedVonMisesThresholdStress', 'd', ('Time', 'nCells',))
fltVM = fout.createVariable('floatingVonMisesThresholdStress', 'd', ('Time', 'nCells',))

values=[
        150.0,  # DML
        200.0,  # Enderby
        150.0,  # Amery
        150.0,  # Phillipi, Denman

        200.0,  # Totten
        200.0,  # Mertz
        200.0,  # Victoria Land
        100.0,  # Ross

        150.0,  # Getz
        200.0,  # Thwaites/PIG
        200.0,  # Bellingshausen
        100.0,  # George VI

        100.0,  # Larsen A-C
        100.0,  # Larsen E
        150.0,  # FRIS
        100.0,  # Brunt-Stancomb
        ]

grdVM[:]=100.0e3
fltVM[:]=100.0e3
for r in range(nRegions):
    mask = np.nonzero(regionCellMasks[:,r] == 1)[0]
    grdVM[0, mask] = values[r] * 1000.0
    fltVM[0, mask] = values[r] * 1000.0

fout.close()
f.close()
