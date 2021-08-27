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
nRows = nFiles * 2
nCols = nFiles

fig, axs = plt.subplots(nrows=nRows, ncols=nCols, sharex=True, sharey=False)

if nRows > 1 or nCols > 1:
    axs = axs.ravel() #easier to index with flattened axs array; remove last subplot if nFiles is odd

floatValue = 4 # value in cellMask denoting floating ice
dynamicValue = 2 # value in cellMask denoting dynamic ice.
for ii, filename in enumerate(filenames):
    f = Dataset(filename, 'r')
    g = Dataset(filename.replace('output_all_timesteps', 'globalStats'))
    f.set_auto_mask(False)
    deltat = np.gradient(f.variables["daysSinceStart"][:]) * s_per_day
    yr = f.variables["daysSinceStart"][:] / 365.
    nCells = f.dimensions["nCells"].size
    cellsOnCell = f.variables["cellsOnCell"][:]

    GLflux = g.variables["groundingLineFlux"][:]
    GLflux = np.interp(yr, g.variables["daysSinceStart"][:] / 365., GLflux)
    thk = f.variables["thickness"][:]
    sfcMassBal = f.variables["sfcMassBalApplied"][:]
    basalMassBal = f.variables["basalMassBal"][:]
    faceMeltingThickness = f.variables["faceMeltingThickness"][:] #m
    calvingThickness = f.variables["calvingThickness"][:]
    xCell = f.variables["xCell"][:]
    areaCell = f.variables["areaCell"][:]
    
    cellMask = f.variables["cellMask"][:]
    floatMask = (cellMask & floatValue) // floatValue
    dynamicMask = (cellMask & dynamicValue) // dynamicValue
    groundedMask = (thk > 1) - floatMask
    masks = [groundedMask, floatMask]
    
    # add non-dynamic cells fringing grounded ice to groundedMask
    print(nCells)
    print('Adding non-dynamic cells neighboring grounded ice to groundedMask')
    for iTime in np.arange(0, len(deltat)):
        for iCell in np.arange(0, nCells):
            if (dynamicMask[iTime, iCell] == 0) and (np.sum(groundedMask[iTime, cellsOnCell[iCell]-1]) >= 1):
                groundedMask[iTime, iCell] = 1
                floatMask[iTime, iCell] = 0

    print('Finished adding non-dynamic cells to groundedMask')

    # loop over floating and grounded masks
    for jj, mask in enumerate(masks):
        cellAreaArray = np.tile(areaCell, (np.shape(calvingThickness)[0],1))

        totalVol = np.sum(thk * mask * cellAreaArray, axis=1)
        calvingVolFlux = np.sum(calvingThickness * mask * cellAreaArray,axis=1) #m^3
        faceMeltVolFlux = np.sum(faceMeltingThickness * mask * cellAreaArray,axis=1) # m^3
        sfcMassBalVolFlux = np.sum(sfcMassBal * mask * cellAreaArray, axis=1) / 910. * deltat
        basalMassBalVolFlux = np.sum(basalMassBal * mask * cellAreaArray, axis=1) / 910. * deltat
        GLvolFlux = GLflux * deltat / 3.154e7 / 910. #m^3

        if mask is groundedMask:
            title = 'Grounded Ice'
            GLvolFlux = GLvolFlux * -1. #mass loss for grounded ice
            basalMassBalVolFlux *= 0. # I don't know why, but there are huge negative spikes for grounded ice, likely related to surges?
        elif mask is floatMask:
            title = 'Floating Ice'
            GLvolFlux = 0. * GLvolFlux #np.cumsum(totalVol - totalVol[0] - sfcMassBalVolFlux + faceMeltVolFlux + calvingVolFlux)

        massBudget = sfcMassBalVolFlux + basalMassBalVolFlux - faceMeltVolFlux - calvingVolFlux + GLvolFlux

        if mask is groundedMask:
            GLfluxPlot, = axs[jj].plot(yr, np.cumsum(GLvolFlux), c='black')
        elif mask is floatMask:
            GLfluxPlot, = axs[jj].plot(yr, totalVol - totalVol[0] - np.cumsum(massBudget) + massBudget[0], c='black')
        #budgetSumPlot, = axs[jj].plot(yr, np.cumsum(massBudget) - massBudget[0], c='tab:blue');
        faceMeltPlot, = axs[jj].plot(yr, np.cumsum(-faceMeltVolFlux), c='tab:purple')
        basalMassBalPlot, = axs[jj].plot(yr, np.cumsum(basalMassBalVolFlux), c='tab:cyan')
        sfcMassBalPlot, = axs[jj].plot(yr, np.cumsum(sfcMassBalVolFlux), c='tab:pink')
        calvingPlot, = axs[jj].plot(yr, np.cumsum(-calvingVolFlux), c='tab:green')
        totalVolChangePlot, = axs[jj].plot(yr, totalVol - totalVol[0], c='tab:orange', linestyle='dotted'); 
        axs[jj].set_xlabel('yrs')
        axs[jj].set_ylabel('volume change (m^3)')
        axs[jj].grid()
        axs[jj].set_title(title)

axs[0].legend([GLfluxPlot, faceMeltPlot, sfcMassBalPlot,  basalMassBalPlot, calvingPlot, totalVolChangePlot],
               ['GL flux', 'undercutting', 'SMB', 'BMB', 'calving', 'total volume change'])

plt.show()
