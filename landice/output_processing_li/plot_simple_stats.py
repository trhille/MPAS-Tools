#!/usr/bin/env python
'''
Script to plot desired time-series from landice global or regional stats files.
Currently only useful for whole-AIS simulations.

Trevor Hillebrand 11/2023
'''

import sys
import os
import re
import string
import numpy as np
from netCDF4 import Dataset
from argparse import ArgumentParser
import matplotlib.pyplot as plt

rhoi = 910.0
rhosw = 1028.

print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = ArgumentParser(description=__doc__)
parser.add_argument("-d", dest="directories", help="input directories; space delimited", type=str, nargs='*', default=None, metavar="FILENAME")
parser.add_argument("-u", dest="units", help="units for mass/volume: m3, kg, Gt", default="Gt", metavar="FILENAME")
parser.add_argument("-n", dest="fileRegionNames", help="region name filename.  If not specified, will attempt to read region names from file 1.", metavar="FILENAME")
parser.add_argument("-r", dest="plot_regions", type=int, nargs='*',
                    help=("indices of regions to plot. Space delimited.",
                          "If not specified, will plot all available regions."),
                    default=None)  # Ross, ASE, FRIS = 7,9,14
parser.add_argument("-v", dest="plot_var", help="name of variable to plot", default="volumeAboveFloatation")
parser.add_argument('-s', dest='save_filename', default=None, help='Path to save .png to, if desired.')
parser.add_argument('--start_year', dest='start_year', default=0., type=float, help='Year value to assign at beginning of time series.')
parser.add_argument("--regional", dest="regional", help="Whether to plot regional stats", action='store_true')
parser.add_argument("--sensitivity", dest="sensitivity",
                    help="Whether this is a sensitivity test vs a core ISMIP6 experiment",
                    action='store_true')
parser.add_argument("--plot_obs", dest="plot_obs", help="Whether to plot observations", action='store_true')
parser.add_argument("--hist", dest="hist", help="Whether this is a historical run only", action='store_true')
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

in_directories = args.directories
# Get nRegions and yr from first file
f = Dataset(in_directories[0] + "/globalStats.nc", 'r')
yr = f.variables['daysSinceStart'][:]/365.0
f.close()

# Create a dictionary with entries run_file : {hist_directory}
# We will append plot setting to this later.
run_dict = {}
for directory in in_directories:
    if directory is None:
        continue
    if 'hist' in directory:
        hist_directory = None
    else:
        if "ctrl" in directory:
            hist_directory = re.sub("ctrlAE", "hist", directory)
        else:
            hist_directory = re.sub("expAE..", 'hist', directory)
    run_dict[directory] = {}
    run_dict[directory]['hist_directory'] = hist_directory

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
        'ISMIP6BasinBC': {'name': 'Amery-Lambert', 'color': 'tab:brown', 'input': [73, 10], 'outflow': [77,4], 'net': [-4, 11], 'shelfMelt': [35.5]},
        'ISMIP6BasinCCp': {'name': 'Phillipi, Denman', 'color': 'tab:red', 'input': [81, 13], 'outflow': [87,7], 'net':[-7,15], 'shelfMelt': [107.9]},
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
    fn = Dataset(args.fileRegionNames, 'r')
    nRegions = len(fn.dimensions['nRegions'])
    rNamesIn = fn.variables['regionNames'][:]
    fn.close()

    if args.plot_regions is None:
        plot_regions = range(nRegions)
    else:
        #plot_regions = [int(i) for i in args.plot_regions.split(',')]
        plot_regions = args.plot_regions

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


# Hard-code plot settings for ISMIP6 runs. These will only apply if plotting global stats,
# as regional stats will need to use different colors for different regions.
# Color corresponds to climate model, linestyle to scenario (solid for 8.5, dashed for 2.6),
# and alpha (opacity) corresponds to Forcing (opaque for 2300, semi-transparent for Repeat).
# Hydrofracture runs will all have their own panel, so no need for special plot settings for those.
global_plot_settings = {
        'hist': {'color': 'black', 'linestyle': 'solid', 'alpha': 1., 'tier': '1', 'esm': None},
        'ctrlAE': {'color': 'black', 'linestyle': 'solid', 'alpha': 1., 'tier': '1', 'esm': "control"},
        'expAE01': {'color': 'tab:gray', 'linestyle': 'dashed', 'alpha': 0.6, 'tier': '1', 'esm': "NorESM1-M-RCP2.6-repeat"},
        'expAE02': {'color': 'tab:blue', 'linestyle': 'solid', 'alpha': 1., 'tier': '1', 'esm': "CCSM4-RCP8.5-2300"},
        'expAE03': {'color': 'tab:orange', 'linestyle': 'solid', 'alpha': 1., 'tier': '1', 'esm': "HadGEM2-RCP8.5-2300"},
        'expAE04': {'color': 'tab:red', 'linestyle': 'solid', 'alpha': 1., 'tier': '1', 'esm': "CESM2-ssp5-85-2300"},
        'expAE05': {'color': 'tab:pink', 'linestyle': 'solid', 'alpha': 1., 'tier': '1', 'esm': "UKESM-ssp5-85-2300"},
        'expAE06': {'color': 'tab:pink', 'linestyle': 'solid', 'alpha': 0.6, 'tier': '1', 'esm': "UKESM-ssp5-85-repeat"},
        'expAE07': {'color': 'tab:gray', 'linestyle': 'solid', 'alpha': 0.6, 'tier': '2a', 'esm': "NorESM1-M-RCP8.5-repeat"},
        'expAE08': {'color': 'tab:orange', 'linestyle': 'solid', 'alpha': 0.6, 'tier': '2a', 'esm': "HadGEM2-RCP8.5-repeat"},
        'expAE09': {'color': 'tab:red', 'linestyle': 'solid', 'alpha': 0.6, 'tier': '2a', 'esm': "CESM2-ssp5-85-repeat"},
        'expAE10': {'color': 'tab:pink', 'linestyle': 'dashed', 'alpha': 1, 'tier': '2a', 'esm': "UKESM-ssp1-26-2300"},
        'expAE11': {'color': 'tab:blue', 'linestyle': 'solid', 'alpha': 1., 'tier': '2b', 'esm': "CCSM4-RCP8.5-2300-h"},
        'expAE12': {'color': 'tab:orange', 'linestyle': 'solid', 'alpha': 1., 'tier': '2b', 'esm': "HadGEM2-RCP8.5-2300-h"},
        'expAE13': {'color': 'tab:red', 'linestyle': 'solid', 'alpha': 1., 'tier': '2b', 'esm': "CESM2-ssp5-85-2300-h"},
        'expAE14': {'color': 'tab:pink', 'linestyle': 'solid', 'alpha': 1., 'tier': '2b', 'esm': "UKESM-ssp5-85-h"}
        }
# Determine plot settings for each run based on the above dict.
for run_item in run_dict:
    for scen_item in global_plot_settings:
        if scen_item in run_item:  #i.e., if the substring of the experiment name occurs in the run name
            run_dict[run_item]['color'] = global_plot_settings[scen_item]['color']
            run_dict[run_item]['linestyle'] = global_plot_settings[scen_item]['linestyle']
            run_dict[run_item]['alpha'] = global_plot_settings[scen_item]['alpha']
            run_dict[run_item]['tier'] = global_plot_settings[scen_item]['tier']
            run_dict[run_item]['label'] = global_plot_settings[scen_item]['esm']

if args.sensitivity and not args.regional:
    ncol = 1
    nrow = 1
    col_scale = 6
    row_scale = 4
    sharey = [False]
elif args.hist and args.regional:
    ncol = 2
    nrow = 1
    col_scale = 6
    row_scale = 4
    sharey = [False, False] 
elif args.regional:
    ncol = 3
    nrow = 1
    col_scale = 4
    row_scale = 4
    sharey = [False, False, True]
else:
    ncol = 3
    nrow = 1
    col_scale = 4
    row_scale = 4
    sharey = [False, True, True]

# Set up Figure 1: volume stats overview
fig1 = plt.figure(figsize=(col_scale*ncol, row_scale*nrow), num=1)
axs1 = []
for ind, share in zip(range(1, ncol * nrow + 1), sharey):
    if share:
        share_ax = axs1[-1]  # share with previous axis. Not very general
    else:
        share_ax = None
    axs1.append(fig1.add_subplot(nrow, ncol, ind, sharey=share_ax))

alphabet = list(string.ascii_lowercase)  # For plot labels
for ii, ax in enumerate(axs1):
    ax.grid()
    ax.set_xlabel('Year', fontsize=13)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.text(0.01, 1.01, f'({alphabet[ii]})', transform=ax.transAxes, fontsize=12)

# Set up unit conversion factors to be used when reading variables
if args.units == "m3":
    volUnitFactor = 1.0
    massUnitFactor = 1.0 / rhoi
elif args.units == "kg":
    volUnitFactor = rhoi
    massUnitFactor = 1.0
elif args.units == "Gt":
    if "olumeAboveFloatation" in args.plot_var:
        # more convenient to display as 10^6 Gt
        args.units = "10$^6$ Gt"
        volUnitFactor = rhoi / 1.e6 / 1.0e12
        massUnitFactor = 1.0 / 1.e6 / 1.0e12
        axs1[0].set_ylabel(f'Cumulative change in mass\nabove floatation ({args.units})')
    else:
        volUnitFactor = rhoi / 1.0e12
        massUnitFactor = 1.0 / 1.0e12
        if "roundedIceVolume" in args.plot_var:
            axs1[0].set_ylabel(f'Grounded ice mass change ({args.units})', fontsize=13)
else:
    sys.exit("ERROR: Unknown unit specified")

def VAF_to_sea_level(vol):
    return vol / volUnitFactor / 3.62e14 * rhoi / rhosw * 1000.

def sea_level_to_VAF(vol):
    return vol * volUnitFactor * 3.62e14 * rhosw / rhoi / 1000. 

def add_sea_lev_ax(axName):
    sea_lev_ax = axName.secondary_yaxis('right', functions=(VAF_to_sea_level, sea_level_to_VAF))
    return sea_lev_ax
def plot_stats(directory, file):
    if directory is None:
        return
    # If a cleaned file exists, plot that. Conversely, if a cleaned
    # filed is specified but does not exist, plot the unclean version.
    # This is useful for hist files, which were added to the list
    # of files dynamically, but usually aren't cleaned.
    fname = directory + "/" + file + "Stats.nc"

    if "regionalStats" in fname:
        regional = True
    else:
        regional = False

    if os.path.isfile(fname + '.cleaned'):
        pltname = fname + '.cleaned'
        print(f'Found cleaned file {pltname}. Plotting from that.')
    elif 'cleaned' in fname and not os.path.isfile(fname):
        tmpname = fname.removesuffix('.cleaned')
        if os.path.isfile(tmpname):
            print(f'{fname} does not exist, plotting {tmpname} instead.')
            pltname = tmpname
        else:
            print(f'{fname} and {tmpname} do not exist. Skipping.')
            return
    else:
        pltname = fname

    print("Reading and plotting file: {}".format(pltname))

    # Hard-code some settings for ISMIP6 sensitivity tests.
    if args.sensitivity:
        if 'gamma21000' in fname or 'm10' in fname or 'enthalpy' in fname or '2km' in fname:
            linestyle = 'dashed'
        elif 'gamma9620' in fname or 'm1/' in fname or 'no_thermal' in fname:
            linestyle = 'dotted'
        elif "m3" in fname:
            linestyle = 'dashdot'
        else:
            linestyle = 'solid'
        if not args.regional:
             axs = [axs1[0]]

    if "depth_integrated" in fname or "3rd_order" in fname:
        linewidth = 1
    else:
        linewidth = 2

    if regional:
        if 'AE02' in fname or 'AE11' in fname:
            axs = [axs1[1]]  # for list comprehension plotting, below
        elif 'AE03' in fname or 'AE12' in fname:
            axs = [axs1[2]]  # for list comprehension plotting, below
        elif 'hist' in fname or 'ctrlAE' in fname:
            if args.hist:
                axs = [axs1[1]]
            else:
                axs = axs1

    if not args.sensitivity and not args.regional:
        if run_dict[directory]['tier'] == '1':
           axs = [axs1[0]]
        if run_dict[directory]['tier'] == '2a':
           axs = [axs1[1]]
        if run_dict[directory]['tier'] == '2b':
           axs = [axs1[2]]

    f = Dataset(pltname,'r')
    yr = f.variables['daysSinceStart'][:]/365.0
    dt = f.variables['deltat'][:]/(3600.0*24.0*365.0) # in yr

    # Deal with different variable names for global and regional stats
    if regional:
        plot_var_name = "regional" + args.plot_var[0].capitalize() + args.plot_var[1:]
    else:
        plot_var_name = args.plot_var

    plot_var = f.variables[plot_var_name][:]

    if run_dict[directory]['hist_directory'] is not None:
        hist = Dataset(run_dict[directory]['hist_directory'] + "/" + file + "Stats.nc", 'r')
        hist_yr = hist.variables['daysSinceStart'][:]/365.0
        hist_dt = hist.variables['deltat'][:]/(3600.0*24.0*365.0) # in yr
        hist_plot_var = hist.variables[plot_var_name][:]

        # Concatenate hist data with exp data
        yr = np.concatenate((hist_yr, yr))
        dt = np.concatenate((hist_dt, dt))
        plot_var = np.concatenate((hist_plot_var, plot_var))

    # Fig 1: summary plot
    if "olume" in args.plot_var:  # leave out the "v" to make it case insensitive
        plot_var *= volUnitFactor

    plot_var = plot_var[:] - plot_var[0]
    if args.regional and regional:
        dtnR = np.tile(dt.reshape(len(dt),1), (1,nRegions))  # repeated per region with dim of nt,nRegions
        nRegionsLocal = len(f.dimensions['nRegions'])
        if nRegionsLocal != nRegions:
            sys.exit(f"ERROR: Number of regions in file {fname} does not match number of regions in first input file!")
        for r in plot_regions:
            if linestyle == 'solid' and linewidth == 2:  # Only add solid curves to legend to avoid repetition
                label = rNames[r]
            else:
                label = None
            lines, = [ax.plot(args.start_year + yr, plot_var[:,r], linestyle=linestyle,
                          color=region_colors[r], linewidth=linewidth,
                          label=label) for ax in axs]
            if args.plot_obs:
                [mn, sig] = ISMIP6basinInfo[rNamesOrig[r]]['net']
                [ax.fill_between(args.start_year + yr, yr*(mn-sig), yr*(mn+sig),
                                 color=region_colors[r], alpha=0.2) for ax in axs]
            axs1[1].legend()


    else:
        if args.sensitivity:
            lstyle = linestyle
            label = None
        else:
            lstyle = run_dict[directory]['linestyle']
            label = run_dict[directory]['label']

        if args.regional and not regional:
            axs = [axs1[0]]
        [ax.plot(args.start_year + yr, plot_var, color=run_dict[directory]['color'],
                 linestyle=lstyle,
                 alpha=run_dict[directory]['alpha'],
                 linewidth=linewidth,
                 label=label) for ax in axs]
#        [ax.legend(fontsize='small') for ax in axs]
        if args.plot_obs:
            mnTot=0.0
            sigTot = 0.0
            for reg in ISMIP6basinInfo:
                [mn, sig] = ISMIP6basinInfo[reg]['net']
                mnTot += mn
                sigTot += sig**2
            sigTot = sigTot**0.5
            [ax.fill_between(args.start_year + yr, yr*(mnTot-sigTot),
                             yr*(mnTot+sigTot), color='k', alpha=0.2) for ax in axs]

    f.close()

if args.regional:
    files = ["global", "regional"]
else:
    files = ["global"]

for directory in in_directories:
    for file in files:
        plot_stats(directory, file)

if "olumeAboveFloatation" in args.plot_var:
    for ax in axs1:
        sea_lev_ax = add_sea_lev_ax(ax)
    sea_lev_ax.set_ylabel('Sea-level\nequivalent (mm)')

[ax.legend(fontsize='small') for ax in axs1]

print("Generating plot.")
fig1.tight_layout()
if args.save_filename is not None:
    fig1.savefig(args.save_filename, format="pdf", bbox_inches='tight')
plt.show()

