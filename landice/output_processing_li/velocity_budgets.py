#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 13:37:25 2020

@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

filepath = '/Users/trevorhillebrand/Documents/mpas/MALI_output/Humboldt_3to30km/' + \
            'calving_melting_test/RCP_ensembles/relaxation_tests/faceMelt_vonMises_calving.nc' 
            
data = Dataset(filepath, 'r')

rhoi = 910.

xCell = data.variables["xCell"][:]
yCell = data.variables["yCell"][:]
yrs = data.variables["daysSinceStart"][:] / 365.
deltat = data.variables["deltat"][:]
deltatCellArray = np.tile(deltat, (np.shape(xCell)[0],1)).T
areaCell = data.variables["areaCell"][:]
areaCellTimeArray = np.tile(areaCell, (np.shape(yrs)[0],1))

# Define a mask that cuts out Petermann. 
humboldtMask = (xCell < -285704)
humboldtMaskTimeArray = np.tile(humboldtMask, (np.shape(yrs)[0],1))

thk = data.variables["thickness"][:]
surfaceSpeed = data.variables["surfaceSpeed"][:]
calvingVelocity = data.variables["calvingVelocity"][:]
calvingThickness = data.variables["calvingThickness"][:]
faceMeltSpeed = data.variables["faceMeltSpeed"][:]
faceMeltRateApplied = data.variables["faceMeltRateApplied"][:]
meltingThickness = data.variables["meltingThickness"][:]
groundedMarginMask = data.variables["groundedMarineMarginMask"][:]
dHdt = data.variables["dHdt"][:] / (60. * 60. * 24. * 365.)
dynamicThickening = data.variables["dynamicThickening"][:] / (60. * 60. * 24. * 365.)
sfcMassBal = data.variables["sfcMassBalApplied"][:]


#Create mask of calving and melting for humboldt
calvingMask = (calvingThickness > 0.)
meltingCalvingMask = (groundedMarginMask + calvingMask) * humboldtMaskTimeArray

# First calculate volume change along front from calving and melting. Sum along time axis
volCalved = np.sum(calvingThickness * areaCellTimeArray * meltingCalvingMask, axis=1)
volFaceMelted = np.sum(meltingThickness * areaCellTimeArray * meltingCalvingMask, axis=1)
volAdvected = np.sum(dynamicThickening * deltatCellArray * areaCellTimeArray * meltingCalvingMask, axis=1)
volSMB = np.sum(sfcMassBal / rhoi * deltatCellArray * areaCellTimeArray * meltingCalvingMask, axis=1)

voldHdt = np.sum(dHdt * deltatCellArray * areaCellTimeArray * meltingCalvingMask, axis=1)

# Now balance melting, calving, and ice velocities at ice front
velocityBalance = (4/5*surfaceSpeed - calvingVelocity - faceMeltSpeed) * meltingCalvingMask
frontMigration = velocityBalance * deltatCellArray
frontVelFig, frontVelAx = plt.subplots(2,1)
frontVelAx[0].grid()
velBalPlot = frontVelAx[0].scatter(xCell/1000,yCell/1000, s=15, c=velocityBalance[1,:] * 60 * 60 * 24 * 365)
frontVelAx[0].set_xlim(left=-4.5e2, right=-3e2)
frontVelAx[0].set_ylim(top=-1000, bottom=-1125)
velBal_cbar = plt.colorbar(velBalPlot, ax=frontVelAx[0])
velBal_cbar.set_label('front velocity (m/yr)')

maxFrontMigrationPlot = frontVelAx[1].plot(yrs, np.min(frontMigration, axis=1))
frontVelAx[1].set_xlabel('years')
frontVelAx[1].set_ylabel('Maximum retreat\nduring timestep (m)')

frontVelFig.subplots_adjust(hspace=0.7)

# Plot budgets
plt.rcParams.update({'font.size': 14}) #Set a nice font size
massBudgetFig, massBudgetAx = plt.subplots(1,1)
massBudgetAx.grid()
massBudgetAx.plot(yrs, np.cumsum(volAdvected + volSMB - volCalved - volFaceMelted),
        linewidth=3, color='tab:blue', label='Advected + SMB\n- Calved - Melted')
massBudgetAx.plot(yrs, np.cumsum(voldHdt), color='tab:orange', linestyle='dashdot', label='dHdt')
massBudgetAx.plot(yrs, np.cumsum(volSMB), color='tab:purple', label='SMB')
massBudgetAx.plot(yrs, np.cumsum(-volFaceMelted), color='tab:pink', label='Melted')
massBudgetAx.plot(yrs, np.cumsum(volAdvected), color='tab:green', label='Advected')
massBudgetAx.plot(yrs, np.cumsum(-volCalved), color='black', label='Calved')

massBudgetAx.legend()

massBudgetAx.set_ylabel('Cumulative $\Delta$ Volume (m$^3$)')
massBudgetAx.set_xlabel('years')