#!/usr/bin/env python
'''
Script to plot common time-series from one or more landice regionalStats files.
Currently only useful for whole-AIS simulations.


Matt Hoffman, 8/23/2022
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import matplotlib.pyplot as plt

rhoi = 910.0


print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-1", dest="file1inName", help="input filename", default="globalStats.nc", metavar="FILENAME")
parser.add_option("-2", dest="file2inName", help="input filename", metavar="FILENAME")
parser.add_option("-3", dest="file3inName", help="input filename", metavar="FILENAME")
parser.add_option("-4", dest="file4inName", help="input filename", metavar="FILENAME")
parser.add_option("-5", dest="file5inName", help="input filename", metavar="FILENAME")
parser.add_option("-6", dest="file6inName", help="input filename", metavar="FILENAME")
parser.add_option("-u", dest="units", help="units for mass/volume: m3, kg, Gt", default="Gt", metavar="FILENAME")
parser.add_option("-n", dest="fileRegionNames", help="region name filename.  If not specified, will attempt to read region names from file 1.", metavar="FILENAME")
parser.add_option("-r", dest="plot_regions", help="indices of regions to plot. comma-separated", default=None)  # Ross, ASE, FRIS = 7,9,14
parser.add_option("-v", dest="plot_var", help="name of variable to plot", default="regionalVolumeAboveFloatation")
options, args = parser.parse_args()

print("Using ice density of {} kg/m3 if required for unit conversions".format(rhoi))

if options.units == "m3":
   massUnit = "m$^3$"
elif options.units == "kg":
   massUnit = "kg"
elif options.units == "Gt":
   massUnit = "Gt"
else:
   sys.exit("Unknown mass/volume units")
print("Using volume/mass units of: ", massUnit)

# Get nRegions and yr from first file
f = Dataset(options.file1inName, 'r')
nRegions = len(f.dimensions['nRegions'])
yr = f.variables['daysSinceStart'][:]/365.0
if options.plot_regions is None:
    plot_regions = range(nRegions)
else:
    plot_regions = [int(i) for i in options.plot_regions.split(',')]

# Get region names from file
if options.fileRegionNames:
   fn = Dataset(options.fileRegionNames, 'r')
   rNamesIn = fn.variables['regionNames'][:]
else:
   rNamesIn = f.variables['regionNames'][:]
# Process region names
rNamesOrig = list()
for r in range(nRegions):
    thisString = rNamesIn[r, :].tobytes().decode('utf-8').strip()  # convert from char array to string
    rNamesOrig.append(''.join(filter(str.isalnum, thisString)))  # this bit removes non-alphanumeric chars

# Antarctic data from:
# Rignot, E., Bamber, J., van den Broeke, M. et al. Recent Antarctic ice mass loss from radar interferometry
# and regional climate modelling. Nature Geosci 1, 106-110 (2008). https://doi.org/10.1038/ngeo102
# Table 1: Mass balance of Antarctica in gigatonnes (10^12 kg) per year by sector for the year 2000
# https://www.nature.com/articles/ngeo102/tables/1
# and
# Rignot, E., S. Jacobs, J. Mouginot, and B. Scheuchl. 2013. Ice-Shelf Melting Around Antarctica. Science 341 (6143): 266-70. https://doi.org/10.1126/science.1235798.
# Note: May want to switch to input+, net+
# Note: Some ISMIP6 basins combine multiple Rignot basins.  May want to separate if we update our regions.
ISMIP6basinInfo = {
        'ISMIP6BasinAAp': {'name': 'Dronning Maud Land', 'input': [60,9], 'outflow': [60,7], 'net': [0, 11], 'shelfMelt': [57.5]},
        'ISMIP6BasinApB': {'name': 'Enderby Land', 'input': [39,5], 'outflow': [40,2], 'net': [-1,5], 'shelfMelt': [24.6]},
        'ISMIP6BasinBC': {'name': 'Amery-Lambert', 'input': [73, 10], 'outflow': [77,4], 'net': [-4, 11], 'shelfMelt': [35.5]},
        'ISMIP6BasinCCp': {'name': 'Phillipi, Denman', 'input': [81, 13], 'outflow': [87,7], 'net':[-7,15], 'shelfMelt': [107.9]},
        'ISMIP6BasinCpD': {'name': 'Totten', 'input': [198,37], 'outflow': [207,13], 'net': [-8,39], 'shelfMelt': [102.3]},
        'ISMIP6BasinDDp': {'name': 'Mertz', 'input': [93,14], 'outflow': [94,6], 'net': [-2,16], 'shelfMelt': [22.8]},
        'ISMIP6BasinDpE': {'name': 'Victoria Land', 'input': [20,1], 'outflow': [22,3], 'net': [-2,4], 'shelfMelt': [22.9]},
        'ISMIP6BasinEF': {'name': 'Ross', 'input': [61+110,(10**2+7**2)**0.5], 'outflow': [49+80,(4**2+2^2)**0.5], 'net': [11+31,(11*2+7**2)**0.5], 'shelfMelt': [70.3]},
        'ISMIP6BasinFG': {'name': 'Getz', 'input': [108,28], 'outflow': [128,18], 'net': [-19,33], 'shelfMelt': [152.9]},
        'ISMIP6BasinGH': {'name': 'Thwaites/PIG', 'input': [177,25], 'outflow': [237,4], 'net': [-61,26], 'shelfMelt': [290.9]},
        'ISMIP6BasinHHp': {'name': 'Bellingshausen', 'input': [51,16], 'outflow': [86,10], 'net': [-35,19], 'shelfMelt': [76.3]},
        'ISMIP6BasinHpI': {'name': 'George VI', 'input': [71,21], 'outflow': [78,7], 'net': [-7,23], 'shelfMelt': [152.3]},
        'ISMIP6BasinIIpp': {'name': 'Larsen A-C', 'input': [15,5], 'outflow': [20,3], 'net': [-5,6], 'shelfMelt': [32.9]},
        'ISMIP6BasinIppJ': {'name': 'Larsen E', 'input': [8,4], 'outflow': [9,2], 'net': [-1,4], 'shelfMelt': [4.3]},
        'ISMIP6BasinJK': {'name': 'FRIS', 'input': [93+142, (8**2+11**2)**0.5], 'outflow': [75+145,(4**2+7**2)**0.5], 'net': [18-4,(9**2+13**2)**0.5], 'shelfMelt': [155.4]},
        'ISMIP6BasinKA': {'name': 'Brunt-Stancomb', 'input': [42+26,(8**2+7**2)**0.5], 'outflow': [45+28,(4**2+2**2)**0.5], 'net':[-3-1,(9**2+8**2)**0.5], 'shelfMelt': [10.4]}
        }

# Parse region names to more usable names, if available
rNames = [None]*nRegions
for r in range(nRegions):
    if rNamesOrig[r] in ISMIP6basinInfo:
        rNames[r] = ISMIP6basinInfo[rNamesOrig[r]]['name']
    else:
        rNames[r] = rNamesOrig[r]

#print(rNames)

ncol = 2
nrow = 1

# Set up Figure 1: volume stats overview
fig1, axs1 = plt.subplots(nrow, ncol, figsize=(8, 6), num=1, sharex=True, sharey=True)
for ax in axs1:
    ax.grid()
    ax.set_xlabel('Year')
axs1[0].set_ylabel(f'Cumulative mass change ({options.units})')

# Set up unit conversion factors to be used when reading variables
if options.units == "m3":
    volUnitFactor = 1.0
    massUnitFactor = 1.0 / rhoi
elif options.units == "kg":
    volUnitFactor = rhoi
    massUnitFactor = 1.0
elif options.units == "Gt":
    volUnitFactor = rhoi / 1.0e12
    massUnitFactor = 1.0 / 1.0e12
else:
    sys.exit("ERROR: Unknown unit specified")


def plotStat(fname, addToLegend=False):
    if fname is None:
        return
    print("Reading and plotting file: {}".format(fname))

    name = fname

    if 'AE02' in fname:
       ax = axs1[0]
    elif 'AE03' in fname:
       ax = axs1[1]

    if 'gamma21000' in fname or 'm10/' in fname or 'no_thermal' in fname:
       linestyle = 'dashed'
    elif 'gamma9620' in fname or 'm1/' in fname:
       linestyle = 'dotted'
    else:
       linestyle = 'solid'

    colors = ['tab:green', 'tab:purple', 'tab:cyan']

    f = Dataset(fname,'r')
    yr = f.variables['daysSinceStart'][:]/365.0
    dt = f.variables['deltat'][:]/(3600.0*24.0*365.0) # in yr
    #yr = yr-yr[0]  # uncomment to align all start dates
    dtnR = np.tile(dt.reshape(len(dt),1), (1,nRegions))  # repeated per region with dim of nt,nRegions
    nRegionsLocal = len(f.dimensions['nRegions'])
    if nRegionsLocal != nRegions:
        sys.exit(f"ERROR: Number of regions in file {fname} does not match number of regions in first input file!")

    # Fig 1: summary plot
    if "Volume" in options.plot_var:
        plot_var = f.variables[options.plot_var][:] * volUnitFactor
    else:
        plot_var = f.variables[options.plot_var][:]
    plot_var = plot_var[:,:] - plot_var[0,:]
    for r, color in zip(plot_regions, colors):
        ax.plot(yr, plot_var[:,r], linestyle=linestyle, color=color, linewidth=2)

    f.close()

for in_file in [options.file1inName, options.file2inName, options.file3inName, options.file4inName, options.file5inName, options.file6inName]:
    plotStat(in_file)

print("Generating plot.")
fig1.tight_layout()
plt.show()

