#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
import netCDF4
import jigsawpy
from matplotlib import pyplot as plt

path2nc = 'DEM_2002.85300894_shifted_H_V_EPSG3031.nc'
# get needed fields from dataset
f = netCDF4.Dataset(path2nc,'r+')
f.set_auto_mask(False) # disable masked arrays

x1 = f.variables['x'][:]
y1 = f.variables['y'][:]
sfc = f.variables['surfaceBerthier'][:]

# subset data - optional
step=1
x1=x1[::step]
y1=y1[::step]
sfc=sfc[::step, ::step]

dx = x1[1] - x1[0] # assumed constant and equal in x and y
nx = len(x1)
ny = len(y1)

sz = sfc.shape

# define extent of region to mesh
xx0=x1.min(); xx1=x1.max();
yy0=y1.min(); yy1=y1.max();
geom_points = np.array([ # list of xy "node" coordinates
       ((xx0, yy0), 0),
       ((xx1, yy0), 0),
       ((xx1, yy1), 0),
       ((xx0, yy1), 0)],
       dtype=jigsawpy.jigsaw_msh_t.VERT2_t)

geom_edges = np.array([    # list of "edges" between nodes
       ((0, 1), 0),
       ((1, 2), 0),
       ((2, 3), 0),
       ((3, 0), 0)],
       dtype=jigsawpy.jigsaw_msh_t.EDGE2_t)

# flood fill to remove island, icecaps, etc.
searchedMask = np.zeros(sz)
floodMask = np.zeros(sz)
iStart=sz[0]//2
jStart=sz[1]//2
floodMask[iStart,jStart]=1

neighbors=np.array([[1,0], [-1,0], [0,1], [0,-1]])

lastSearchList = np.ravel_multi_index([[iStart],[jStart]], sz, order='F')

sfc[sfc<=5.0] = 0.0
   # flood fill -------------------
cnt = 0
print("Beginning flood fill")
while len(lastSearchList) > 0:
       #print(cnt)
       cnt += 1
       newSearchList = np.array([], dtype='i')

       for iii in range(len(lastSearchList)):
           [i, j] = np.unravel_index(lastSearchList[iii], sz, order='F')
           # search neighbors
           for n in neighbors:
               ii=i+n[0]; jj=j+n[1];  # subscripts to neighbor
               if ii<sz[0] and jj<sz[1] and searchedMask[ii,jj] == 0:  # only consider unsearched neighbors
                   searchedMask[ii,jj] = 1  # mark as searched

                   if sfc[ii, jj] > 0.0:
                       floodMask[ii,jj] = 1  # mark as ice
                       newSearchList = np.append(newSearchList, np.ravel_multi_index([[ii],[jj]], sz, order='F')[0]) # add to list of newly found  cells
       lastSearchList = newSearchList
       # optional - plot flood mask
       #plt.pcolor(floodMask)
       #plt.show()


# apply flood fill
sfc[floodMask==0] = np.nan
print("Flood fill complete")
    #   # make masks -------------------
neighbors=np.array([[1,0], [-1,0], [0,1], [0,-1], [1,1], [-1,1], [1,-1], [-1,-1]])
    
iceMask = sfc>0.0
marginMask = np.zeros(sz, dtype='i')
for n in neighbors:
    marginMask = np.logical_or(marginMask, np.logical_not(np.roll(iceMask, n, axis=[0,1])))
marginMask = np.logical_and(marginMask, iceMask) # where ice exists and neighbors non-ice locations
    
f.variables["surfaceBerthier"][:]=np.ma.MaskedArray(sfc)
print("done")
f.close()

