#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 18:06:25 2020
This script adds necessary variables and time dimension to an MPAS mesh
@author: trevorhillebrand
"""
from netCDF4 import Dataset
import numpy as np

meshPath = '/Users/trevorhillebrand/Documents/mpas/MALI_output/Thwaites/thwaites_1km_meshonly.nc'

mesh = Dataset(meshPath, 'r+')

# create dimensions and variables
#mesh.createDimension('Time')
#mesh.createDimension('nVertLevels', size=10)
#
#mesh.createVariable('thickness', 'f', ('Time', 'nCells'))
#mesh.createVariable('bedTopography', 'f', ('Time', 'nCells'))
#mesh.createVariable('sfcMassBal', 'f', ('Time', 'nCells'))
mesh.createVariable('basalHeatFlux', 'f', ('Time', 'nCells'))
#mesh.createVariable('floatingBasalMassBal', 'f', ('Time', 'nCells'))
#mesh.createVariable('temperature', 'f', ('Time', 'nCells', 'nVertLevels'))
#mesh.createVariable('surfaceAirTemperature', 'f', ('Time', 'nCells'))
#mesh.createVariable('beta', 'f', ('Time', 'nCells'))
#mesh.createVariable('observedSurfaceVelocityX', 'f', ('Time', 'nCells'))
#mesh.createVariable('observedSurfaceVelocityY', 'f', ('Time', 'nCells'))
#mesh.createVariable('observedSurfaceVelocityUncertainty', 'f', ('Time', 'nCells'))
#mesh.createVariable('observedThicknessTendency', 'f', ('Time', 'nCells'))
#mesh.createVariable('observedThicknessTendencyUncertainty', 'f', ('Time', 'nCells'))
#mesh.createVariable('thicknessUncertainty', 'f', ('Time', 'nCells'))
#mesh.createVariable('ismip6shelfMelt_basin', 'f', ('Time', 'nCells'))
#mesh.createVariable('ismip6shelfMelt_deltaT', 'f', ('Time', 'nCells'))
#


mesh.close()