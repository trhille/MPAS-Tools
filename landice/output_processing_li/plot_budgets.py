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

filename = '/lustre/scratch4/turquoise/trhille/Humboldt_1to10km_r02_20210112/m20/floatKill_TFminus0.5/output_all_timesteps.nc' 
         

f = Dataset(filename, 'r')
f.set_auto_mask(False)
rhoi = 910.
s_per_day = 86400.

deltat = np.gradient(f.variables["daysSinceStart"][:]) * s_per_day
yr = f.variables["daysSinceStart"][:] / 365.

thk = f.variables["thickness"][:]
sfcMassBal = f.variables["sfcMassBalApplied"][:]
faceMeltRateApplied = f.variables["faceMeltRateApplied"][:] #m/s
calvingThickness = f.variables["calvingThickness"][:]
groundedCalvingThickness = calvingThickness * f.variables["groundedMarineMarginMask"][:]# m
floatingCalvingThickness = calvingThickness * (1 - f.variables["groundedMarineMarginMask"][:])
xCell = f.variables["xCell"][:]
areaCell = f.variables["areaCell"][:]

cellAreaArray = np.tile(areaCell, (np.shape(calvingThickness)[0],1))

totalVol = np.sum(thk * cellAreaArray, axis=1)
calvingVolFlux = np.sum(calvingThickness * cellAreaArray,axis=1) #m^3i
groundedCalvingVolFlux = np.sum(groundedCalvingThickness * cellAreaArray,axis=1) #m^3
floatingCalvingVolFlux = np.sum(floatingCalvingThickness * cellAreaArray,axis=1) #m^3
faceMeltVolFlux = np.sum(faceMeltRateApplied, axis=1) * deltat # m^3
sfcMassBalVolFlux = np.sum(sfcMassBal * cellAreaArray, axis=1) / 910. * deltat

massBudget = sfcMassBalVolFlux - faceMeltVolFlux - calvingVolFlux

budgetSumPlot, = plt.plot(yr, np.cumsum(massBudget) - massBudget[0], c='tab:blue');
faceMeltPlot, = plt.plot(yr, np.cumsum(-faceMeltVolFlux), c='tab:purple')
sfcMassBalPlot, = plt.plot(yr, np.cumsum(sfcMassBalVolFlux), c='tab:pink')
groundedCalvingPlot, = plt.plot(yr, np.cumsum(-groundedCalvingVolFlux), c='tab:green', linestyle='dashed')
floatingCalvingPlot, = plt.plot(yr, np.cumsum(-floatingCalvingVolFlux), c='tab:green', linestyle='dotted')
calvingPlot, = plt.plot(yr, np.cumsum(-calvingVolFlux), c='tab:green')
totalVolChangePlot, = plt.plot(yr, totalVol - totalVol[0], c='tab:orange', linestyle='dotted'); 
plt.xlabel('yrs')
plt.ylabel('volume change (m^3)')
plt.legend([budgetSumPlot, faceMeltPlot, sfcMassBalPlot, groundedCalvingPlot, floatingCalvingPlot, calvingPlot, totalVolChangePlot],
           ['total budget', 'faceMelt', 'sfcMassBal', 'grounded calving', 'floating calving', 'total calving', 'total volume change'])
plt.grid()

plt.show()
