#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 19:15:45 2020
This creates a mass budget (melting, calving, SMB) from model output.nc file. 
Output fields must include thickness, calvingThickness, faceMeltRateApplied, sfcMassBalApplied,
groundedMarineMarginMask

@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-f", dest="filename", help="filename for plotting", metavar="FILENAME")
options, args = parser.parse_args()

filenames = options.filename.split(',') #create list of filenames to loop over
nFiles = len(filenames)
rhoi = 910.
s_per_day = 86400.

if nFiles > 3:
    nRows = 2
    nCols = int(nFiles / 2 + nFiles % 2)
else:
    nRows = 1
    nCols = nFiles

fig, axs = plt.subplots(nrows=nRows, ncols=nCols, sharex=True, sharey=True)

axs = axs.ravel() #easier to index with flattened axs array; remove last subplot if nFiles is odd
if nFiles % 2 == 1:
    fig.delaxes(axs[-1])

for filename,ax in zip(filenames,axs):
    f = Dataset(filename, 'r')
    f.set_auto_mask(False)

    deltat = np.gradient(f.variables["daysSinceStart"][:]) * s_per_day
    yr = f.variables["daysSinceStart"][:] / 365.

    thk = f.variables["thickness"][:]
    sfcMassBal = f.variables["sfcMassBalApplied"][:]
    faceMeltingThickness = f.variables["faceMeltingThickness"][:] #m
    calvingThickness = f.variables["calvingThickness"][:]
    xCell = f.variables["xCell"][:]
    areaCell = f.variables["areaCell"][:]

    cellAreaArray = np.tile(areaCell, (np.shape(calvingThickness)[0],1))

    totalVol = np.sum(thk * cellAreaArray, axis=1)
    calvingVolFlux = np.sum(calvingThickness * cellAreaArray,axis=1) #m^3
    faceMeltVolFlux = np.sum(faceMeltingThickness * cellAreaArray,axis=1) # m^3
    sfcMassBalVolFlux = np.sum(sfcMassBal * cellAreaArray, axis=1) / 910. * deltat

    massBudget = sfcMassBalVolFlux - faceMeltVolFlux - calvingVolFlux

    budgetSumPlot, = ax.plot(yr, np.cumsum(massBudget) - massBudget[0], c='tab:blue');
    faceMeltPlot, = ax.plot(yr, np.cumsum(-faceMeltVolFlux), c='tab:purple')
    sfcMassBalPlot, = ax.plot(yr, np.cumsum(sfcMassBalVolFlux), c='tab:pink')
    calvingPlot, = ax.plot(yr, np.cumsum(-calvingVolFlux), c='tab:green')
    totalVolChangePlot, = ax.plot(yr, totalVol - totalVol[0], c='tab:orange', linestyle='dotted'); 
    ax.set_xlabel('yrs')
    ax.set_ylabel('volume change (m^3)')
    ax.grid()

axs[0].legend([budgetSumPlot, faceMeltPlot, sfcMassBalPlot,  calvingPlot, totalVolChangePlot],
               ['total budget', 'undercutting', 'SMB', 'calving', 'total volume change'])

plt.show()
