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

parser.add_option("-u", dest="units", help="units for mass/volume: m3, kg, Gt", default="m3", metavar="FILENAME")
options, args = parser.parse_args()

ensembleDirs = options.ensembleDirs.split(',') # split ensemble directories into list
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

#set colormap for plots
colormap = mpl.colors.TABLEAU_COLORS
#colormap = mpl.colors.XKCD_COLORS
colorlist = list(colormap.items())
#Get units
f = Dataset(controlFiles[0], 'r')
units=f.variables[options.variableName].units
f.close()
#set linestyles to loop through (e.g., to separate out RCP scenarios)
linestyleList = ['solid', 'dashed', 'dotted', 'dashdot']
linestyleIndex = 0 # initialize for loop

# create axes to plot into
varFig, varAx = plt.subplots(1,2, sharey=True, sharex=True)
#ratioFig, ratioAx = plt.subplots(1,1)
varAx[0].grid()
varAx[1].grid()
#ratioAx.grid()
plotLines = [] #empty list to fill with lines
plotLineNames = [] #empty list to fill with filenames
plotBounds = []
plotBoundNames = []

def VAF2seaLevel(vol):
    return -vol / 3.62e14 * rhoi / rhosw * 1000.

def seaLevel2VAF(vol):
    return -vol * 3.62e14 * rhosw / rhoi / 1000.

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
                       controlData.variables[options.variableName][:])
    
    var2plot1 = f1.variables[options.variableName][:] \
                 - f1.variables[options.variableName][0] \
                 - controlInterp + controlInterp[0]
                
    var2plot2 = np.interp(yr, f2.variables['daysSinceStart'][:]/365.0,
                          f2.variables[options.variableName][:] 
                          - f2.variables[options.variableName][0]) \
                          - controlInterp + controlInterp[0]
    if 'HadGEM2' in boundsFiles[0]:
        plotAx = varAx[1]
    elif 'MIROC5' in boundsFiles[0]:
        plotAx = varAx[0]
        
    tmpFill = plotAx.fill_between(yr+2007., var2plot1, var2plot2, facecolor='tab:grey', alpha = 0.6)
    plotBounds.append(tmpFill)
    plotBoundNames.append(boundsFiles[0])
    
    return plotBounds, plotBoundNames

def plotEnsemble(ensDir, controlFile=None):
#    print("Reading and plotting file: {}".format(fname))
    colorIndex = 0 #initialize index to loop through color list
    ensembleFiles = sorted(os.listdir(ensDir)) # get filenames in directory
    for ensembleMember in ensembleFiles:
        if 'globalStats' in ensembleMember:
            
            f = Dataset(ensDir+ensembleMember,'r')
            yr = f.variables['daysSinceStart'][:]/365.0
            
            # get units
            try: 
                units = f.variables[options.variableName].units
            except:
                units = 'm^3'
        
            var2plot = f.variables[options.variableName][:] \
                         - f.variables[options.variableName][0]
                         
            # subtract off variables from control run
            if controlFile:
                #interpolate control run onto ensemble member time vector
                controlData = Dataset(controlFile, 'r')
                controlInterp = np.interp(yr, controlData.variables['daysSinceStart'][:]/365.0, 
                               controlData.variables[options.variableName][:])
                
                var2plot = var2plot - controlInterp + controlInterp[0]

            if 'HadGEM2' in ensembleMember:
                plotAx = varAx[1]
            elif 'MIROC5' in ensembleMember:
                plotAx = varAx[0]
                
            tmpLine, = plotAx.plot(yr+2007., var2plot, 
                                   label=ensembleMember)
            plotLines.append(tmpLine)
            plotLineNames.append(ensembleMember)
            

            colorIndex += 1 # go to next color
            f.close()
    return plotLines, plotLineNames
    

def addSeaLevAx(axName):
    seaLevAx = axName.secondary_yaxis('right', functions=(VAF2seaLevel, seaLevel2VAF))
    seaLevAx.set_ylabel('Sea-level\nequivalent (mm)', fontsize=16)

controlIndex=0
boundsIndex=0
controlFile=None



for directory in ensembleDirs:
    print("Ensemble {}".format(directory))
    if controlFiles:
        controlFile=controlFiles[controlIndex]
    if boundsDirs and boundsIndex <= len(boundsDirs)-1:
        plotEnsembleBounds(boundsDirs[boundsIndex], controlFile)
    plotLines, plotLineNames = plotEnsemble(directory, controlFile)
    controlIndex += 1
    linestyleIndex += 1
    boundsIndex += 1
    
if options.variableName == "volumeAboveFloatation":
    addSeaLevAx(varAx[-1])

varAx[0].set_xlabel('Year', fontsize=16)
varAx[1].set_xlabel('Year', fontsize=16)
varAx[0].set_ylabel('$\Delta$ {} (${}$)'.format(options.variableName, units), fontsize=16)

#varAx.legend()
varFig.tight_layout()
#varAx.set_ylim(bottom=-7e12, top=0)
varAx[0].set_xlim(left=2000, right=2100.)
varAx[1].set_xlim(left=2000, right=2100.)
#ratioAx.set_xlim(left=0, right=100.)
#set a reasonable fontsize
plt.rcParams.update({'font.size': 16})

# Now plot ratio of change high-res:low-res

#def plotRatios():
#    colorIndex=0
#    ensembleFiles = sorted(os.listdir(ensembleDirs[0]))
#    for ensembleMember in ensembleFiles:
#        if 'globalStats.nc' in ensembleMember:
#            
#            f1 = Dataset(ensembleDirs[0]+ensembleMember,'r')
#            f2 = Dataset(ensembleDirs[1]+ensembleMember,'r')
#            f1Var = f1.variables[options.variableName][:]
#            f2Var = f2.variables[options.variableName][:]
#            yr = f1.variables['daysSinceStart'][:]/365.0
#            
#            # interpolate f2 onto f1 time
#            f2Varinterp = np.interp(yr, f2.variables['daysSinceStart'][:]/365.0, 
#                               f2Var[:])
#            
#           # get units
#           # units = f1.variables[options.variableName].units
#
#           # subtract off variables from control run
#           
#            if controlFiles:
#                #interpolate control run onto ensemble member time vector
#                control1Data = Dataset(controlFiles[0], 'r')
#                control1Interp = np.interp(yr, control1Data.variables['daysSinceStart'][:]/365.0, 
#                               control1Data.variables[options.variableName][:])
#                control2Data = Dataset(controlFiles[1], 'r')
#                control2Interp = np.interp(yr, control2Data.variables['daysSinceStart'][:]/365.0, 
#                               control2Data.variables[options.variableName][:])
#                
#                var2plot = (f2Varinterp - f2Varinterp[0] - (control2Interp - control2Interp[0])) / \
#                            (f1Var - f1Var[0] - (control1Interp - control1Interp[0]))
#            else:
#               var2plot = (f2Varinterp - f2Varinterp[0]) / (f1Var - f1Var[0]) 
#                        
#                        
#            if 'red' in colorlist[colorIndex][0]: #skip red for colorblind safety
#                colorIndex += 1
#                
#            ratioAx.plot(yr, var2plot, color=colorlist[colorIndex][0], 
#                    linestyle=linestyleList[0], label=ensembleMember)
#            colorIndex += 1 # go to next color
#            f1.close()
#            f2.close()
#            
#    return f2Varinterp
#  
#
#plotRatios()
#
#ratioAx.set_xlabel('Year')
#ratioAx.set_ylabel('SLR ratio \n(high-res/low-res)')
#ratioFig.tight_layout()
##set a reasonable fontsize
#plt.rcParams.update({'font.size': 16})
# Special plotting for humboldt ensemble:
for line, lineName in zip(plotLines, plotLineNames):
    if 'smb_only' in lineName:
        line.set_color('tab:pink')
    elif 'draftCalving' in lineName:
        line.set_color('tab:green')
    elif 'VM160' in lineName:
        line.set_color('tab:purple')
    elif 'VM170' in lineName:
        line.set_color('tab:blue')
    elif 'VM180' in lineName:
        line.set_color('tab:cyan')
    if 'm1_' in lineName:
        line.set_linestyle('none')
    if 'm3_' in lineName:
        line.set_linestyle('none')
    if 'm10_' in lineName:
        line.set_linestyle('solid')
    elif 'm25_' in lineName:
        line.set_linestyle('dashed')
        
for bound, boundName in zip(plotBounds, plotBoundNames):
    if 'm3_' in boundName:
        bound.set_color('tab:grey')
        bound.set_edgecolor('none')
        bound.set_alpha(0.4)
        bound.set_zorder(0)
    elif 'm1_' in boundName:
        bound.set_edgecolor('tab:grey')
        bound.set_facecolor('none')
        bound.set_hatch('xxxxxx')
        bound.set_alpha(0.6)
        
varFig.set_size_inches(15, 5)
varFig.subplots_adjust(wspace=0.15)

plt.show()
#varFig.savefig('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_melt_calv_ensemble/Humboldt_only_runs/projections/SL_contributions.png', dpi=300)
