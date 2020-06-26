#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 10:54:15 2020

@author: trevorhillebrand
This scriot compares MALI spinup with observed velocities, thickness, 
and dHdt using 2- and Inf-norms of residuals through time
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import sys
import numpy as np
from netCDF4 import Dataset
from optparse import OptionParser
import math
from collections import OrderedDict
import scipy.spatial
import time
from datetime import datetime
import matplotlib.pyplot as plt

# Parse options. Just need the MALI output .nc file and the .nc file containing observations
print("== Gathering information.  (Invoke with --help for more details. All arguments are optional)\n")
parser = OptionParser()
parser.description = __doc__
parser.add_option("-o", "--output", dest="modelOutputFile", help="MALI output.nc file.", default="output.nc", metavar="FILENAME")
parser.add_option("-d", "--data", dest="observationsFile", help="File containing observations", default="landice_grid.nc", metavar="FILENAME")
parser.add_option("-s", "--save", dest="saveFigs", help="Save figures or not?", default="False", metavar="FILENAME")
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = parser.parse_args()



mod = Dataset(options.modelOutputFile, 'r')
mod.set_auto_mask(False)

obs = Dataset(options.observationsFile, 'r')
obs.set_auto_mask(False)


# Check that files are on the same grid. If not, interpolate model onto observations grid.
x_mod = mod.variables["xCell"][:]
y_mod = mod.variables["yCell"][:]
x_obs = obs.variables["xCell"][:]
y_obs = obs.variables["yCell"][:]


#For now, just give a warning. In future, could add function to automatically do the interpolation without overwriting files
if np.sum(x_mod - x_obs) != 0 or np.sum(y_mod - y_obs) != 0:
    print("Data and model output are not on the same grid. Use interpolate_to_mpasli_grid.py")
else:
    x = x_mod 
    y = y_mod 
    
meshDensity = mod.variables["meshDensity"][:]
s_per_yr = 3600. * 24. * 365.

#Load pertinent fields for comparison
obsThicknessOrig = obs.variables["thickness"][0,:]
obsSurfaceSpeedOrig = s_per_yr * np.sqrt(obs.variables["observedSurfaceVelocityX"][0,:]**2 + \
                       obs.variables["observedSurfaceVelocityX"][0,:]**2)
obs_dHdtOrig = s_per_yr * obs.variables["observedThicknessTendency"][0,:]  # convert dHdt from m/s to m/yr

#Remove unphysical values from observations and replace with NaNs
obsThickness = obsThicknessOrig
obsSurfaceSpeed = obsSurfaceSpeedOrig
obs_dHdt = obs_dHdtOrig

obsThickness[np.abs(obsThicknessOrig) > 5e3] = np.nan
obsSurfaceSpeed[np.abs(obsSurfaceSpeedOrig) > 30e3] = np.nan
obs_dHdt[np.abs(obs_dHdtOrig) > 100.] = np.nan

modThicknessOrig = mod.variables["thickness"][:]
modSurfaceSpeedOrig = s_per_yr * mod.variables["surfaceSpeed"][:]
mod_dHdtOrig = mod.variables["dHdt"][:]

modTime = mod.variables["daysSinceStart"][:]/365 #Time in yrs for plotting
nTime = len(modTime)

# Close netCDFs 
obs.close()
mod.close()


# Remove cells where there is no data. Save a copy to preserve size for mapping
modThickness = modThicknessOrig[:, np.isfinite(obsThicknessOrig)]
modSurfaceSpeed = modSurfaceSpeedOrig[:, np.isfinite(obsSurfaceSpeedOrig)]
mod_dHdt = mod_dHdtOrig[:, np.isfinite(obs_dHdtOrig)]

obsThickness = obsThickness[np.isfinite(obsThicknessOrig)]
obsSurfaceSpeed = obsSurfaceSpeed[np.isfinite(obsSurfaceSpeedOrig)]
obs_dHdt = obs_dHdt[np.isfinite(obs_dHdtOrig)]

# calculate residuals
thicknessRes = modThickness - np.tile(obsThickness, (nTime,1))
surfaceSpeedRes = modSurfaceSpeed - np.tile(obsSurfaceSpeed, (nTime, 1))
dHdtRes = mod_dHdt - np.tile(obs_dHdt, (nTime, 1))

# calculate L-2 and L-infinity norms for residuals (model - obs)
thicknessRes_2norm = np.linalg.norm(thicknessRes, ord=2, axis=1)
thicknessRes_InfNorm = np.linalg.norm(thicknessRes, ord=np.inf, axis=1)

surfaceSpeedRes_2norm = np.linalg.norm(surfaceSpeedRes, ord=2, axis=1)
surfaceSpeedRes_InfNorm = np.linalg.norm(surfaceSpeedRes, ord=np.inf, axis=1)

dHdtRes_2norm = np.linalg.norm(dHdtRes, ord=2, axis=1)
dHdtRes_InfNorm = np.linalg.norm(dHdtRes, ord=np.inf, axis=1)

# calculate stats
thicknessRes_mean = np.mean(thicknessRes, axis=1)
thicknessRes_stdev = np.std(thicknessRes, axis=1)
thicknessRMSE = np.sqrt(np.mean(thicknessRes**2, axis=1))
surfaceSpeedRes_mean = np.mean(surfaceSpeedRes, axis=1)
surfaceSpeedRes_stdev = np.std(surfaceSpeedRes, axis=1)
surfaceSpeedRMSE = np.sqrt(np.mean(surfaceSpeedRes**2, axis=1))
dHdtRes_mean = np.mean(dHdtRes, axis=1)
dHdtRes_stdev = np.std(dHdtRes, axis=1)
dHdtRMSE = np.sqrt(np.mean(dHdtRes**2, axis=1))

# Plot up norms of residuals over time
plt.ioff() # turn interactive mode off to plot multiple figures
save_dpi = 250 # dpi at which to save figures

thick_fig, thick_axs = plt.subplots(3)
thick_fig.suptitle('Thickness residuals (m)')
thick_axs[0].plot(modTime, thicknessRMSE)
thick_axs[0].set_ylabel('RMSE')
thick_axs[1].plot(modTime, thicknessRes_InfNorm)
thick_axs[1].set_ylabel('Inf norm')
thick_axs[2].plot(modTime, np.abs(thicknessRes_mean))
thick_axs[2].set_xlabel('Time (yrs)')
thick_axs[2].set_ylabel('Mean abs. error')
plt.subplots_adjust(hspace=0.5)

#plot surfaceSpeed residual stats
speed_fig, speed_axs = plt.subplots(3)
speed_fig.suptitle('surfaceSpeed residuals (m/yr)')
speed_axs[0].plot(modTime, (surfaceSpeedRMSE))
speed_axs[0].set_ylabel('RMSE')
speed_axs[1].plot(modTime, (surfaceSpeedRes_InfNorm))
speed_axs[1].set_ylabel('Inf norm')
speed_axs[2].plot(modTime, np.abs(surfaceSpeedRes_mean))
speed_axs[2].set_xlabel('Time (yrs)')
speed_axs[2].set_ylabel('Mean abs. error')
plt.subplots_adjust(hspace=0.5)

#Plot dHdt residuals
dHdt_fig, dHdt_axs = plt.subplots(3)
dHdt_fig.suptitle('dHdt residuals (m/yr)')
dHdt_axs[0].plot(modTime, dHdtRMSE)
dHdt_axs[0].set_ylabel('RMSE')
dHdt_axs[1].plot(modTime, dHdtRes_InfNorm)
dHdt_axs[1].set_ylabel('Inf norm')
dHdt_axs[2].plot(modTime, np.abs(dHdtRes_mean))
dHdt_axs[2].set_xlabel('Time (yrs)')
dHdt_axs[2].set_ylabel('Mean abs. error')
plt.subplots_adjust(hspace=0.5)


#marker size based on meshDensity for plotting maps
mark_s = 3 * (np.abs(np.log10(meshDensity)) + 0.1)

#map of thickness residuals
thickmap_fig = plt.figure(4)
thickmap = plt.scatter(x/1000., y/1000., s=mark_s, c=modThicknessOrig[-1,:] - obsThicknessOrig, \
                       vmin = -200, vmax = 200, cmap = 'seismic')
thickmap_cbar = plt.colorbar(thickmap)
thickmap_cbar.set_label('thickness residuals (m)', rotation=270, labelpad = 12)


#map of surface speed residuals
speedmap_fig = plt.figure(5)
speedmap = plt.scatter(x/1000., y/1000.,s=mark_s,c=np.log10(np.abs(modSurfaceSpeedOrig[-1,:] - \
                        obsSurfaceSpeedOrig)), cmap='magma', vmin=-2.5, vmax=2.5)
speedmap_cbar = plt.colorbar(speedmap)
speedmap_cbar.set_label('log10 abs surfaceSpeed residuals ($10^x$ m/yr)', rotation=270, labelpad = 12)

#map of dHdt residuals
dHdt_map_fig = plt.figure(6)
dHdt_map = plt.scatter(x/1000., y/1000., s=mark_s, \
                       c=np.log10(np.abs(mod_dHdtOrig[-1,:] - obs_dHdtOrig)),
                       cmap='BuPu', vmin=-2, vmax=2)
dHdt_map_cbar = plt.colorbar(dHdt_map)
dHdt_map_cbar.set_label('log10 dHdt residuals ($10^x$ m/yr)', rotation=270, labelpad = 12)

#Save figures
if options.saveFigs=='True':
    thick_fig.savefig('thick_tmp.png', bbox_inches='tight', dpi=save_dpi)
    speed_fig.savefig('speed_tmp.png', bbox_inches='tight', dpi=save_dpi)
    dHdt_fig.savefig('dHdt_tmp.png', bbox_inches='tight', dpi=save_dpi)
    thickmap_fig.savefig('thickmap_tmp.png', bbox_inches='tight', dpi=save_dpi)
    speedmap_fig.savefig('speedmap_tmp.png', bbox_inches='tight', dpi=save_dpi)
    dHdt_map_fig.savefig('dHdtmap_tmp.png', bbox_inches='tight', dpi=save_dpi)

plt.show()  # show figures



    