#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mark ice-edge cells with Dirichlet boundary conditions within an arbitrary region defined by a shapefile.
If no shapefile is specified, the script will mark all ice-edge cells with Dirichlet boundary conditions.

Authors: Trevor Hillebrand, Matthew Hoffman
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
from optparse import OptionParser
from netCDF4 import Dataset
import geopandas as gpd
from matplotlib.path import Path

print("** Gathering information.")
parser = OptionParser()
parser.add_option("-f", "--file", dest="file", help="grid file to modify; default: landice_grid.nc", metavar="FILE")
parser.add_option("-s", "--shapefile", dest="shapefile", default=None, help="optional shapefile mask to search within", metavar="FILE")
options, args = parser.parse_args()

f = Dataset(options.file,'r+')
f.set_auto_mask(False)

# Parse data from netCDF file
xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]
nVertInterfaces = len(f.dimensions['nVertInterfaces'])
dirichletVelocityMask = f.variables['dirichletVelocityMask'][0,:,:]
thickness = f.variables['thickness'][0, :]
cellsOnCell = f.variables['cellsOnCell'][:]
nCells = len(f.dimensions['nCells'])
nEdgesOnCell = f.variables['nEdgesOnCell'][:]
keepCellMask = (thickness.copy() * 0) # define cell mask by thickness
keepCellMask[thickness > 0.0] = 1
print('Num of cells with ice: {}'.format(sum(keepCellMask)))

#create a copy of dirichletVelocityMask to re-mark edges at the end
dirichletEdges = dirichletVelocityMask.copy()

#initialize mask in case shapefile is not used
shpMask = np.ones(np.shape(xCell))

def mask2shapefile(options, xCell, yCell):
    #read in shapefile as a geopandas geodataframe
    shp = gpd.read_file(options.shapefile)
    shp.to_crs(epsg=3413, inplace=True)
    
    #convert xCell and yCell to a list of tuples for comparison with shapefile
    wholeDomainCoordinates = []
    for ii in np.arange(0, len(xCell)):
        tmptup = (xCell[ii], yCell[ii])
        wholeDomainCoordinates.append(tmptup)
    
    #get exterior coordinates from shapefile polygon 
    shpExtCoordinates = []
    shpX, shpY = shp.loc[0, 'geometry'].exterior.coords.xy
    for ii in np.arange(0, len(shpX)):
        tmptup = [shpX[ii], shpY[ii]]
        shpExtCoordinates.append(tmptup)
    
    # Convert shapefile coordinates to polygon and determine which 
    # of the cells in the model domain lie within the prescribed polgyon
    shpExtPolygon = Path(shpExtCoordinates)
    shpMask = shpExtPolygon.contains_points(wholeDomainCoordinates)
   
    return(shpMask)
    
    
def findMarginCells(keepCellMask, nCells, cellsOnCells, nEdgesOnCell):
    # find list of margin cells
    iceCells = np.nonzero(keepCellMask == 1)[0]
    for i in range(len(iceCells)):
        iCell = iceCells[i]
        for neighbor in cellsOnCell[iCell,:nEdgesOnCell[iCell]]-1:  # the -1 converts from the fortran indexing in the variable to python indexing
            if thickness[neighbor] == 0.0:
                dirichletVelocityMask[iCell,:] = 1
    return(dirichletVelocityMask)
    

dirichletVelocityMask = findMarginCells(keepCellMask, nCells, cellsOnCell, nEdgesOnCell)

if options.shapefile != None:
    shpMask = mask2shapefile(options, xCell, yCell)
                
# apply the shapefile mask to all layers of dirichletVelocityMask
dirichletVelocityMask = (dirichletVelocityMask * np.tile(shpMask, (nVertInterfaces,1)).T      
                          + dirichletEdges)

# In case, you've run this more than once and added old dirichlet masks together, reset all values to 0 and 1
dirichletVelocityMask[dirichletVelocityMask>0] = 1

f.variables['dirichletVelocityMask'][0,:] = dirichletVelocityMask

f.close()