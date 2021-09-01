#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 21:10:36 2020
Plot ensemble volume and sea level. Plots a desired variable from globalStats.nc
files within a given set of directories. Runs in the same directory are plotted
in different colors (up to 9; change colormap if more are needed) with the same
linestyle. Additional directories will loop through the same colors with a 
different linestyle. Thus, this script is useful for separating out runs by a 
single characteristic (linestlye). For example, RCP2.6 runs with dashed lines,
RCP8.5 runs with solid lines. Also plots globa mean sea level equivalent 
for volumeAboveFloatation.
@author: trevorhillebrand
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import os
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt
import matplotlib as mpl

rhoi = 910.0
rhosw = 1028.0
print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-d", dest="ensembleDirs", help="directory containing ensemble members (strings separated by commas; no spaces)", metavar="FILENAME")
parser.add_option("-b", dest="boundsDirs", help="directory containing ensemble members (strings separated by commas; no spaces)", metavar="FILENAME")
parser.add_option("-v", dest="variableName", help="variable(s) to plot, separated by commas", default = "volumeAboveFloatation", metavar="FILENAME")
parser.add_option("-c", dest="controlFiles", help="comma-separated paths to control run(s) to subtract from ensemble members", metavar="FILENAME")
parser.add_option("--plotChange", dest="plotChange", help="Use this flag to plot the change since start", default=False)

parser.add_option("-u", dest="units", help="units for mass/volume: m3, kg, Gt", default="m3", metavar="FILENAME")
options, args = parser.parse_args()

if options.plotChange == 'True':
    plotChange = True
elif options.plotChange == 'False':
    plotChange = False
else:
    print('Invalid value for plotChange. Should either be True or False. Setting to False.')
    
ensembleDirs = options.ensembleDirs.split(',') # split ensemble directories into list
variableName = options.variableName.split(',')

if options.controlFiles:
    controlFiles = options.controlFiles.split(',') # split control files into list
else:
    controlFiles = None
    
if options.boundsDirs:
    boundsDirs = options.boundsDirs.split(',') # split control files into list
else:
    boundsDirs = None

print("Using ice density of {} kg/m3".format(rhoi))
print("Using seawater density of {} kg/m3".format(rhosw))

# create axes to plot into
nCols = 3
varFig, varAx = plt.subplots(len(variableName), nCols , sharey='row', sharex=True)
varAx = varAx.ravel()

for ax in varAx:
    ax.grid()

plotLines = [] #empty list to fill with lines
plotLineNames = [] #empty list to fill with filenames
plotBounds = []
plotBoundNames = []

def plotEnsembleBounds(boundsDir, controlFile=None):
    boundsFiles = sorted(os.listdir(boundsDir)) # get filenames in directory
    
    for boundsFile in boundsFiles:
        if 'globalStats' not in boundsFile:
            boundsFiles.remove(boundsFile)
            
            
    f1 = Dataset(boundsDir + boundsFiles[0], 'r')
    yr = f1.variables['daysSinceStart'][:] / 365.0
    f2 = Dataset(boundsDir + boundsFiles[1], 'r')
    if controlFile:
        #interpolate control run onto ensemble member time vector
        controlData = Dataset(controlFile, 'r')
        controlInterp = np.interp(yr, controlData.variables['daysSinceStart'][:]/365.0, 
                       controlData.variables[variable][:])
    
    var2plot1 = f1.variables[variable][:] \
                 - f1.variables[variable][0] \
                 - controlInterp + controlInterp[0]
                
    var2plot2 = np.interp(yr, f2.variables['daysSinceStart'][:]/365.0,
                          f2.variables[variable][:] 
                          - f2.variables[variable][0]) \
                          - controlInterp + controlInterp[0]
    if 'CNRM' in boundsFiles[0]:
        plotAx = varAx[row*nCols + 2]
    elif 'HadGEM2' in boundsFiles[0]:
        plotAx = varAx[row*nCols + 1]
    elif 'MIROC5' in boundsFiles[0]:
        plotAx = varAx[row*nCols]
        
    tmpFill = plotAx.fill_between(yr+2007., var2plot1, var2plot2, facecolor='tab:grey', alpha = 0.6)
    plotBounds.append(tmpFill)
    plotBoundNames.append(boundsFiles[0])
    
    return plotBounds, plotBoundNames

def plotEnsemble(ensDir, controlFile=None):
    ensembleFiles = sorted(os.listdir(ensDir)) # get filenames in directory
    for ensembleMember in ensembleFiles:
        if 'globalStats' in ensembleMember:
            
            f = Dataset(ensDir+ensembleMember,'r')
            yr = f.variables['daysSinceStart'][:]/365.0
            
            initVol = f.variables["volumeAboveFloatation"][0]
                
            if plotChange is True:
                var2plot = f.variables[variable][:] \
                         - f.variables[variable][0]
                initVol = 0.
            else: 
                var2plot = f.variables[variable][:]
                         
            # subtract off variables from control run
            if controlFile:
                #interpolate control run onto ensemble member time vector
                controlData = Dataset(controlFile, 'r')
                controlInterp = np.interp(yr, controlData.variables['daysSinceStart'][:]/365.0, 
                               controlData.variables[variable][:])

                var2plot = var2plot - controlInterp + controlInterp[0]

            if 'CNRM' in ensembleMember:
                plotAx = varAx[row*nCols + 2]
            elif 'HadGEM2' in ensembleMember:
                plotAx = varAx[row*nCols + 1]
            elif 'MIROC5' in ensembleMember:
                plotAx = varAx[row*nCols]
                
            tmpLine, = plotAx.plot(yr+2007., var2plot, 
                                   label=ensembleMember)
            plotLines.append(tmpLine)
            plotLineNames.append(ensembleMember)

            f.close()
    return plotLines, plotLineNames, initVol

# The following functions for plotting sea-level equiv axis
def VAF2seaLevel(vol):
    return -(vol-initVol) / 3.62e14 * rhoi / rhosw * 1000.

def seaLevel2VAF(vol):
    return -(vol-initVol) * 3.62e14 * rhosw / rhoi / 1000.

def addSeaLevAx(axName):
    seaLevAx = axName.secondary_yaxis('right', functions=(VAF2seaLevel, seaLevel2VAF))
    seaLevAx.set_ylabel('Sea-level\ncontribution (mm)', fontsize=16)

f = Dataset(controlFiles[0], 'r')
for row, variable in enumerate(variableName):
    controlIndex=0
    boundsIndex=0
    controlFile=None
    for directory in ensembleDirs:
        print("Ensemble {}".format(directory))
        if controlFiles:
            controlFile=controlFiles[controlIndex]
        if boundsDirs and boundsIndex <= len(boundsDirs)-1:
            plotEnsembleBounds(boundsDirs[boundsIndex], controlFile)
        plotLines, plotLineNames, initVol = plotEnsemble(directory, controlFile)
        controlIndex += 1
        boundsIndex += 1
    
    units=f.variables[variable].units
    if variable == "volumeAboveFloatation":
        addSeaLevAx(varAx[nCols - 1])
    if variable == 'volumeAboveFloatation' and plotChange is True:
        varAx[row * nCols].set_ylabel('$\Delta$ volume above\nfloatation (m$^3$)', fontsize=16)
    elif variable == 'volumeAboveFloatation' and plotChange is  False:
        varAx[row * nCols].set_ylabel('volume above\nfloatation (m$^3$)', fontsize=16)
    elif variable != 'volumeAboveFloatation' and plotChange is True:
        varAx[row * nCols].set_ylabel('$\Delta$ {} (${}$)'.format(variable, units), fontsize=16)
    else:
        varAx[row * nCols].set_ylabel('{} (${}$)'.format(variable, units), fontsize=16)

f.close()

# Special plotting for humboldt ensemble:        
for ax in varAx[-nCols::]:
    ax.set_xlabel('Year', fontsize=16)

varAx[0].set_title('MIROC5')
varAx[1].set_title('HadGEM2')
varAx[2].set_title('CNRM')

varFig.tight_layout()
varAx[0].set_xlim(left=2007, right=2100.)
plt.rcParams.update({'font.size': 16})

for line, lineName in zip(plotLines, plotLineNames):

    if 'shelfMelt10myr' in lineName:
        line.set_linewidth(1)
    elif 'shelfMelt20myr' in lineName:
        line.set_linewidth(2)
    elif 'shelfMelt30myr' in lineName:
        line.set_linewidth(3)
    if 'm5_' in lineName or 'm10_' in lineName:
        lowCalving = 'VM180'
        medCalving = 'VM170'
        highCalving = 'VM160'
    elif 'm7_' in lineName or 'm25_' in lineName:
        if 'HadGEM2' in lineName or 'CNRM' in lineName:
            lowCalving = 'VM190'
            medCalving = 'VM180'
            highCalving = 'VM170'
        elif 'MIROC5' in lineName:
            lowCalving = 'VM180'
            medCalving = 'VM170'
            highCalving = 'VM160'
    elif 'm1_' in lineName:
        lowCalving = 'VM180'
        medCalving = 'VM170'
        highCalving = 'VM150'
    if 'smb_only' in lineName:
        line.set_color('tab:pink')
    elif '2017calvingFront' in lineName:
        line.set_color('tab:grey')
    elif 'calvingVelocityData' in lineName:
        line.set_color('tab:grey')
    elif highCalving in lineName:
        line.set_color('tab:purple')
    elif medCalving in lineName:
        line.set_color('tab:blue')
    elif lowCalving in lineName:
        line.set_color('tab:cyan')
    if 'm1_' in lineName:
        line.set_color('grey')
    if 'm3_' in lineName:
        line.set_linestyle('none')
    if 'm5_' in lineName or 'm10_' in lineName:
        line.set_linestyle('solid')
    elif 'm7_' in lineName or 'm25_' in lineName:
        line.set_linestyle('dashed')
    
    if 'm10_' in lineName or 'm25_' in lineName:
        line.set_color('tab:pink')

for bound, boundName in zip(plotBounds, plotBoundNames):
    if 'm3_' in boundName:
        bound.set_color('tab:grey')
        bound.set_edgecolor('none')
        bound.set_alpha(0.4)
        bound.set_zorder(0)
    elif 'm1_' in boundName:
        bound.set_edgecolor('tab:grey')
        bound.set_facecolor('tab:grey')
        #bound.set_hatch('xxxxxx')
        bound.set_alpha(0.6)
        
varFig.set_size_inches(nCols * 5, len(varAx) * 2)
varFig.subplots_adjust(wspace=0.15, hspace=0.15)

plt.show()