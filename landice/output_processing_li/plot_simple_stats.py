#!/usr/bin/env python
'''
Script to plot desired time-series from landice global or regional stats files.
Currently only useful for whole-AIS simulations.

Trevor Hillebrand 11/2023
'''

import sys
import os
import re
import numpy as np
from netCDF4 import Dataset
from argparse import ArgumentParser
import matplotlib.pyplot as plt

rhoi = 910.0


print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = ArgumentParser(description=__doc__)
parser.add_argument("-f", dest="files", help="input filenames, space delimited", type=str, nargs='*', default="globalStats.nc", metavar="FILENAME")
parser.add_argument("-u", dest="units", help="units for mass/volume: m3, kg, Gt", default="Gt", metavar="FILENAME")
parser.add_argument("-n", dest="fileRegionNames", help="region name filename.  If not specified, will attempt to read region names from file 1.", metavar="FILENAME")
parser.add_argument("-r", dest="plot_regions", help=("indices of regions to plot. comma-separated.",
                                                   "If not specified, will plot all available regions."),
                  default=None)  # Ross, ASE, FRIS = 7,9,14
parser.add_argument("-v", dest="plot_var", help="name of variable to plot", default="regionalVolumeAboveFloatation")
parser.add_argument("--regional", dest="regional", help="Whether to plot regional stats", action='store_true')
args = parser.parse_args()

print("Using ice density of {} kg/m3 if required for unit conversions".format(rhoi))

if args.units == "m3":
   massUnit = "m$^3$"
elif args.units == "kg":
   massUnit = "kg"
elif args.units == "Gt":
   massUnit = "Gt"
else:
   sys.exit("Unknown mass/volume units")
print("Using volume/mass units of: ", massUnit)

in_files = args.files
# Get nRegions and yr from first file
f = Dataset(in_files[0], 'r')
yr = f.variables['daysSinceStart'][:]/365.0

# Create a dictionary with entries run_file : {hist_file}
# We will append plot setting to this later.
run_dict = {}
for file in in_files:
    if file is None:
        continue
    if 'hist' in file:
        hist_file = None
    else:
        if "ctrl" in file:
            hist_file = re.sub("ctrlAE", "hist", file)
        else:
            hist_file = re.sub("expAE..", 'hist', file)
        if 'cleaned' in hist_file and not os.path.isfile(hist_file):
            tmpname = hist_file.removesuffix('.cleaned')
            if os.path.isfile(tmpname):
                print(f'{hist_file} does not exist, plotting {tmpname} instead.')
                hist_file = tmpname
            else:
                print(f'{fname} and {tmpname} do not exist. Skipping.')
                hist_file = None
    run_dict[file] = {}
    run_dict[file]['hist_file'] = hist_file

# Antarctic data from:
# Rignot, E., Bamber, J., van den Broeke, M. et al. Recent Antarctic ice mass loss from radar interferometry
# and regional climate modelling. Nature Geosci 1, 106-110 (2008). https://doi.org/10.1038/ngeo102
# Table 1: Mass balance of Antarctica in gigatonnes (10^12 kg) per year by sector for the year 2000
# https://www.nature.com/articles/ngeo102/tables/1
# and
# Rignot, E., S. Jacobs, J. Mouginot, and B. Scheuchl. 2013. Ice-Shelf Melting Around Antarctica. Science 341 (6143): 266-70. https://doi.org/10.1126/science.1235798.
# Note: May want to switch to input+, net+
# Note: Some ISMIP6 basins combine multiple Rignot basins.  May want to separate if we update our regions.
# Also assign each region a color for plotting. These were chosen such that the three major WAIS basins would
# look good when plotted together. The rest are fairly random, with an attempt to avoid colors that are too similar.
ISMIP6basinInfo = {
        'ISMIP6BasinAAp': {'name': 'Dronning Maud Land', 'color': 'tab:blue', 'input': [60,9], 'outflow': [60,7], 'net': [0, 11], 'shelfMelt': [57.5]},
        'ISMIP6BasinApB': {'name': 'Enderby Land', 'color': 'tab:orange', 'input': [39,5], 'outflow': [40,2], 'net': [-1,5], 'shelfMelt': [24.6]},
        'ISMIP6BasinBC': {'name': 'Amery-Lambert', 'color': 'tab:red', 'input': [73, 10], 'outflow': [77,4], 'net': [-4, 11], 'shelfMelt': [35.5]},
        'ISMIP6BasinCCp': {'name': 'Phillipi, Denman', 'color': 'tab:brown', 'input': [81, 13], 'outflow': [87,7], 'net':[-7,15], 'shelfMelt': [107.9]},
        'ISMIP6BasinCpD': {'name': 'Totten', 'color': 'tab:pink', 'input': [198,37], 'outflow': [207,13], 'net': [-8,39], 'shelfMelt': [102.3]},
        'ISMIP6BasinDDp': {'name': 'Mertz', 'color': 'tab:gray', 'input': [93,14], 'outflow': [94,6], 'net': [-2,16], 'shelfMelt': [22.8]},
        'ISMIP6BasinDpE': {'name': 'Victoria Land', 'color': 'tab:olive', 'input': [20,1], 'outflow': [22,3], 'net': [-2,4], 'shelfMelt': [22.9]},
        'ISMIP6BasinEF': {'name': 'Ross', 'color': 'tab:green', 'input': [61+110,(10**2+7**2)**0.5], 'outflow': [49+80,(4**2+2^2)**0.5], 'net': [11+31,(11*2+7**2)**0.5], 'shelfMelt': [70.3]},
        'ISMIP6BasinFG': {'name': 'Getz', 'color': 'lightcoral', 'input': [108,28], 'outflow': [128,18], 'net': [-19,33], 'shelfMelt': [152.9]},
        'ISMIP6BasinGH': {'name': 'Thwaites/PIG', 'color': 'tab:purple', 'input': [177,25], 'outflow': [237,4], 'net': [-61,26], 'shelfMelt': [290.9]},
        'ISMIP6BasinHHp': {'name': 'Bellingshausen', 'color': 'teal', 'input': [51,16], 'outflow': [86,10], 'net': [-35,19], 'shelfMelt': [76.3]},
        'ISMIP6BasinHpI': {'name': 'George VI', 'color': 'maroon', 'input': [71,21], 'outflow': [78,7], 'net': [-7,23], 'shelfMelt': [152.3]},
        'ISMIP6BasinIIpp': {'name': 'Larsen A-C', 'color': 'violet', 'input': [15,5], 'outflow': [20,3], 'net': [-5,6], 'shelfMelt': [32.9]},
        'ISMIP6BasinIppJ': {'name': 'Larsen E', 'color': 'blueviolet', 'input': [8,4], 'outflow': [9,2], 'net': [-1,4], 'shelfMelt': [4.3]},
        'ISMIP6BasinJK': {'name': 'FRIS', 'color': 'tab:cyan', 'input': [93+142, (8**2+11**2)**0.5], 'outflow': [75+145,(4**2+7**2)**0.5], 'net': [18-4,(9**2+13**2)**0.5], 'shelfMelt': [155.4]},
        'ISMIP6BasinKA': {'name': 'Brunt-Stancomb', 'color': 'goldenrod', 'input': [42+26,(8**2+7**2)**0.5], 'outflow': [45+28,(4**2+2**2)**0.5], 'net':[-3-1,(9**2+8**2)**0.5], 'shelfMelt': [10.4]}
        }

if args.regional:
    nRegions = len(f.dimensions['nRegions'])
    if args.plot_regions is None:
        plot_regions = range(nRegions)
    else:
        plot_regions = [int(i) for i in args.plot_regions.split(',')]

    # Get region names from file
    if args.fileRegionNames:
       fn = Dataset(args.fileRegionNames, 'r')
       rNamesIn = fn.variables['regionNames'][:]
    else:
       rNamesIn = f.variables['regionNames'][:]
    # Process region names
    rNamesOrig = list()
    for r in range(nRegions):
        thisString = rNamesIn[r, :].tobytes().decode('utf-8').strip()  # convert from char array to string
        rNamesOrig.append(''.join(filter(str.isalnum, thisString)))  # this bit removes non-alphanumeric chars

    # Parse region names to more usable names, if available.
    # Assign colors to each region, using dictionary above.
    rNames = [None]*nRegions
    region_colors = [None]*nRegions
    for r in range(nRegions):
        if rNamesOrig[r] in ISMIP6basinInfo:
            rNames[r] = ISMIP6basinInfo[rNamesOrig[r]]['name']
            region_colors[r] = ISMIP6basinInfo[rNamesOrig[r]]['color']
        else:
            rNames[r] = rNamesOrig[r]

else:
    # Hard-code plot settings for ISMIP6 runs. These will only apply if plotting global stats,
    # as regional stats will need to use different colors for different regions.
    # Color corresponds to climate model, linestyle to scenario (solid for 8.5, dashed for 2.6),
    # and alpha (opacity) corresponds to Forcing (opaque for 2300, semi-transparent for Repeat).
    # Hydrofracture runs will all have their own panel, so no need for special plot settings for those.
    global_plot_settings = {
            'hist': {'color': 'black', 'linestyle': 'solid', 'alpha': 1., 'tier': '1'},
            'ctrlAE': {'color': 'black', 'linestyle': 'solid', 'alpha': 1., 'tier': '1'},
            'expAE01': {'color': 'tab:gray', 'linestyle': 'dashed', 'alpha': 0.6, 'tier': '1'},
            'expAE02': {'color': 'tab:blue', 'linestyle': 'solid', 'alpha': 1., 'tier': '1'},
            'expAE03': {'color': 'tab:orange', 'linestyle': 'solid', 'alpha': 1., 'tier': '1'},
            'expAE04': {'color': 'tab:purple', 'linestyle': 'solid', 'alpha': 1., 'tier': '1'},
            'expAE05': {'color': 'tab:pink', 'linestyle': 'solid', 'alpha': 1., 'tier': '1'},
            'expAE06': {'color': 'tab:pink', 'linestyle': 'dashed', 'alpha': 1., 'tier': '1'},
            'expAE07': {'color': 'tab:gray', 'linestyle': 'solid', 'alpha': 0.6, 'tier': '2a'},
            'expAE08': {'color': 'tab:orange', 'linestyle': 'solid', 'alpha': 0.6, 'tier': '2a'},
            'expAE09': {'color': 'tab:purple', 'linestyle': 'solid', 'alpha': 0.6, 'tier': '2a'},
            'expAE10': {'color': 'tab:pink', 'linestyle': 'dashed', 'alpha': 1, 'tier': '2a'},
            'expAE11': {'color': 'tab:blue', 'linestyle': 'solid', 'alpha': 1., 'tier': '2b'},
            'expAE12': {'color': 'tab:orange', 'linestyle': 'solid', 'alpha': 1., 'tier': '2b'},
            'expAE13': {'color': 'tab:purple', 'linestyle': 'solid', 'alpha': 1., 'tier': '2b'},
            'expAE14': {'color': 'tab:pink', 'linestyle': 'solid', 'alpha': 1., 'tier': '2b'}
            }
    # Determine plot settings for each run based on the above dict.
    for run_item in run_dict:
        for scen_item in global_plot_settings:
            if scen_item in run_item:  #i.e., if the substring of the experiment name occurs in the run name
                run_dict[run_item]['color'] = global_plot_settings[scen_item]['color']
                run_dict[run_item]['linestyle'] = global_plot_settings[scen_item]['linestyle']
                run_dict[run_item]['alpha'] = global_plot_settings[scen_item]['alpha']
                run_dict[run_item]['tier'] = global_plot_settings[scen_item]['tier']

if args.regional:
    ncol = 2
    nrow = 1
else:
    ncol = 3
    nrow = 1

# Set up Figure 1: volume stats overview
fig1, axs1 = plt.subplots(nrow, ncol, figsize=(8, 6), num=1, sharex=True, sharey=True)
for ax in axs1:
    ax.grid()
    ax.set_xlabel('Year')
axs1[0].set_ylabel(f'Cumulative mass change ({args.units})')

# Set up unit conversion factors to be used when reading variables
if args.units == "m3":
    volUnitFactor = 1.0
    massUnitFactor = 1.0 / rhoi
elif args.units == "kg":
    volUnitFactor = rhoi
    massUnitFactor = 1.0
elif args.units == "Gt":
    volUnitFactor = rhoi / 1.0e12
    massUnitFactor = 1.0 / 1.0e12
else:
    sys.exit("ERROR: Unknown unit specified")


def plotStat(fname):
    if fname is None:
        return
    # If the cleaned file doesn't exist, plot the unclean version.
    # This is useful for hist files, which were added to the list
    # of files dynamically, but usually aren't cleaned.
    if 'cleaned' in fname and not os.path.isfile(fname):
        tmpname = fname.removesuffix('.cleaned')
        if os.path.isfile(tmpname):
            print(f'{fname} does not exist, plotting {tmpname} instead.')
            fname = tmpname
        else:
            print(f'{fname} and {tmpname} do not exist. Skipping.')
            return

    print("Reading and plotting file: {}".format(fname))

    name = fname

    # Hard-code some settings for ISMIP6 sensitivity tests.
    if args.regional:
        if 'AE02' in fname:
            axs = [axs1[0]]  # for list comprehension plotting, below
        elif 'AE03' in fname:
            axs = [axs1[1]]  # for list comprehension plotting, below
        elif 'hist' in fname or 'ctrlAE' in fname:
            axs = axs1
        if 'gamma21000' in fname or 'm10/' in fname or 'no_thermal' in fname:
           linestyle = 'dashed'
        elif 'gamma9620' in fname or 'm1/' in fname:
           linestyle = 'dotted'
        else:
           linestyle = 'solid'
    else:
        if run_dict[fname]['tier'] == '1':
           axs = [axs1[0]]
        if run_dict[fname]['tier'] == '2':
           axs = [axs1[1]]
        if run_dict[fname]['tier'] == '3':
           axs = [axs1[2]]

    f = Dataset(fname,'r')
    yr = f.variables['daysSinceStart'][:]/365.0
    dt = f.variables['deltat'][:]/(3600.0*24.0*365.0) # in yr
    plot_var = f.variables[args.plot_var][:]

    if run_dict[fname]['hist_file'] is not None:
        hist = Dataset(run_dict[fname]['hist_file'], 'r')
        hist_yr = hist.variables['daysSinceStart'][:]/365.0
        hist_dt = hist.variables['deltat'][:]/(3600.0*24.0*365.0) # in yr
        hist_plot_var = hist.variables[args.plot_var][:]

        # Concatenate hist data with exp data
        yr = np.concatenate((hist_yr, yr))
        dt = np.concatenate((hist_dt, dt))
        plot_var = np.concatenate((hist_plot_var, plot_var))

    # Fig 1: summary plot
    if "Volume" in args.plot_var:
        plot_var *= volUnitFactor

    plot_var = plot_var[:] - plot_var[0]
    if args.regional:
        dtnR = np.tile(dt.reshape(len(dt),1), (1,nRegions))  # repeated per region with dim of nt,nRegions
        nRegionsLocal = len(f.dimensions['nRegions'])
        if nRegionsLocal != nRegions:
            sys.exit(f"ERROR: Number of regions in file {fname} does not match number of regions in first input file!")
        for r in plot_regions:
            lines, = [ax.plot(yr, plot_var[:,r], linestyle=linestyle,
                          color=region_colors[r], linewidth=2,
                          label=rNames[r]) for ax in axs]
        if linestyle == 'solid':  # Only add solid curves to legend to avoid repetition
            axs1[0].legend()

    else:
        [ax.plot(yr, plot_var, color=run_dict[fname]['color'],
                 linestyle=run_dict[fname]['linestyle'],
                 alpha=run_dict[fname]['alpha'], linewidth=2) for ax in axs]

    f.close()

for file in in_files:
    plotStat(file)

print("Generating plot.")
fig1.tight_layout()
plt.show()

