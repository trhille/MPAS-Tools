#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 08:38:40 2020
Plot forcing timeseries
@author: trevorhillebrand
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

runoffData26 = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/ISMIP6_forcings/Humboldt_1to10km_forcings/MIROC5-rcp26/Humboldt_1to10km_MIROC5-rcp26_Runoff.nc', 'r')
runoffData85 = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/ISMIP6_forcings/Humboldt_1to10km_forcings/MIROC5-rcp85/Humboldt_1to10km_MIROC5-rcp85_Runoff.nc', 'r')

SMBData26 = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/ISMIP6_forcings/Humboldt_1to10km_forcings/MIROC5-rcp26/Humboldt_1to10km_MIROC5-rcp26_sfcMassBal.nc','r')
SMBData85 = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/ISMIP6_forcings/Humboldt_1to10km_forcings/MIROC5-rcp85/Humboldt_1to10km_MIROC5-rcp85_sfcMassBal.nc','r')

TFData26 = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/ISMIP6_forcings/Humboldt_1to10km_forcings/MIROC5-rcp26/Humboldt_1to10km_MIROC5-rcp26_oceanThermalForcing.nc','r')
TFData85 = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/ISMIP6_forcings/Humboldt_1to10km_forcings/MIROC5-rcp85/Humboldt_1to10km_MIROC5-rcp85_oceanThermalForcing.nc','r')

thkData = Dataset('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/Humboldt_1to10km_betaonly_041720.nc','r')

iceMask = thkData.variables["thickness"][:]

middleCellID = 3191
xCell = runoffData26.variables["xCell"][:]
yCell = runoffData26.variables["yCell"][:]
areaCell = runoffData26.variables["areaCell"][:]
nTime = runoffData26.dimensions["Time"].size
time = np.arange(1950, 1950 + nTime)
areaCellArray = np.tile(areaCell, (nTime,1))

humboldtMask = (xCell < -285704)
humboldtMaskArray = np.tile(humboldtMask, (nTime,1))
iceMaskArray = np.tile(iceMask, (nTime,1))

runoff_26 = runoffData26.variables["ismip6Runoff"][:] * 3600 * 24. * 365
runoff_85 = runoffData85.variables["ismip6Runoff"][:] * 3600 * 24. * 365
total_runoff_26 = np.sum(runoff_26 * iceMaskArray ,  axis=1)
total_runoff_85 = np.sum(runoff_85 * iceMaskArray , axis=1)

SMB26 = SMBData26.variables["sfcMassBal"][:] * 3600 * 24. * 365
total_SMB26 = np.sum(SMB26 * areaCellArray * iceMaskArray, axis=1)

SMB85 = SMBData85.variables["sfcMassBal"][:] * 3600 * 24. * 365
total_SMB85 = np.sum(SMB85 * areaCellArray * iceMaskArray, axis=1)

TF26 = TFData26.variables["ismip6_2dThermalForcing"][:]
TF85 = TFData85.variables["ismip6_2dThermalForcing"][:]
mean_TF26 = np.mean(TF26 * humboldtMaskArray * (TF26 > 0), axis=1)
mean_TF85 = np.mean(TF85 * humboldtMaskArray * (TF85 > 0), axis=1)

plt.rcParams.update({'font.size': 20})

fig, ax = plt.subplots(3,1, figsize=(8,10))
fig.tight_layout()
fig.subplots_adjust(hspace=0.4)
ax[0].plot(time, runoff_26[:, middleCellID])
ax[0].plot(time, runoff_85[:, middleCellID])
ax[0].set_xlim(left=2000, right=2100)
ax[0].set_ylabel('Runoff\n(kg m$^{-2}$ yr$^{-1}$)')
ax[0].grid()

ax[1].plot(time, total_SMB26)
ax[1].plot(time, total_SMB85)
ax[1].set_ylim(bottom=(1.05*np.min(total_SMB85)), top=1.05*np.max(total_SMB26))
ax[1].set_xlim(left=2000, right=2100)
ax[1].set_ylabel('total SMB\n(kg/yr)')
ax[1].grid()

ax[2].plot(time, mean_TF26)
ax[2].plot(time, mean_TF85)
ax[2].set_xlim(left=2000, right=2100)
ax[2].set_xlabel('Year')
ax[2].set_ylabel('Mean thermal\nforcing (°C)')
ax[2].grid()


