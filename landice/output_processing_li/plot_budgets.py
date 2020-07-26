#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 19:15:45 2020
This creates a mass budget (melting, calving, SMB) from model output.nc file. 
Output fields must include calvingThickness, faceMeltRateApplied, sfcMassBalApplied,
and may include basalMassBal

@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

filename = ('/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_1to10km/calving_'
            + 'melting_test/RCP_ensembles/rcp26/calving_melting_smb.nc')

f = Dataset(filename, 'r')
f.set_auto_mask(False)
rhoi = 910.
s_per_day = 86400.

deltat = np.gradient(f.variables["daysSinceStart"][:]) * s_per_day
yr = f.variables["daysSinceStart"][:] / 365.

thkAnnual = f.variables["thickness"][:]
sfcMassBal = f.variables["sfcMassBalApplied"][:]
faceMeltRateApplied = f.variables["faceMeltRateApplied"][:] #m/s
calvingThickness = f.variables["calvingThickness"][:] # m
xCell = f.variables["xCell"][:]
areaCell = f.variables["areaCell"][:]

#humboldtMask = (xCell < -285704)
#humboldtMaskArray = np.tile(humboldtMask, (np.shape(calvingThickness)[0],1))
cellAreaArray = np.tile(areaCell, (np.shape(calvingThickness)[0],1))

HumboldtVol = np.sum(thkAnnual * cellAreaArray, axis=1)
calvingVolFlux = np.sum(calvingThickness * cellAreaArray,axis=1) #m^3
faceMeltVolFlux = np.sum(faceMeltRateApplied, axis=1) * deltat # m^3
sfcMassBalVolFlux = np.sum(sfcMassBal * cellAreaArray, axis=1) / 910. * deltat

massBudget = sfcMassBalVolFlux - faceMeltVolFlux - calvingVolFlux

plt.plot(yr, np.cumsum(massBudget)); 
plt.plot(yr, HumboldtVol - HumboldtVol[0]); plt.xlabel('yrs'); 
plt.ylabel('volume change (m^3)'); plt.show()

plt.show()
