#!/usr/bin/env python3

import datetime
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


ds = xr.open_dataset("ANT_G1920V01_IceShelfMelt.nc")

# We want the 2002 - 2003 average to be consistent with SMB, velocity, 
# dH/dt, etc. This file doesn't have a numerical time variable, which
# is annoying, but Mar 2002 is index 40 and Dec 2003 is index 47.
mask = ds.ID
nx = ds.sizes["x"]
ny = ds.sizes["y"]

t_start = 40
t_end = 48
n_time = t_end - t_start
melt = ds.melt[t_start:t_end]
time = ds.time[t_start:t_end]

# Now we want to extract the remnant Larsen B (ID value 51)
# and average over cells adjacent to the grounding line (ID value 255)
larsen_b_mask_value = 51
land_mask_value = 255
ocean_mask_value = 0

larsen_b_mask = mask == larsen_b_mask_value
larsen_b_mask_tile = np.tile(larsen_b_mask, (n_time, 1, 1))

melt.values[larsen_b_mask_tile]

# Warning: Slow for-loop! We only need to loop over the Larsen B, though, so let's create
# a little box for it. If you need to run this over large regions, it will take a long
# time, so might want to figure out how to vectorize, but it's quick for Larsen B alone
y_ind, x_ind = np.where(larsen_b_mask)
buffer_size = 2
x_ind_min, x_ind_max = np.min(x_ind) - buffer_size, np.max(x_ind) + buffer_size + 1
y_ind_min, y_ind_max = np.min(y_ind) - buffer_size, np.max(y_ind) + buffer_size + 1

has_land_neighbor, has_neighbor_with_land_neighbor = np.zeros_like(mask), np.zeros_like(mask)
has_ocean_neighbor, has_neighbor_with_ocean_neighbor = np.zeros_like(mask), np.zeros_like(mask)

for y in range(y_ind_min, y_ind_max):
    for x in range (x_ind_min, x_ind_max):
        has_land_neighbor[y, x] = np.any(mask[y-1:y+2, x-1:x+2] == land_mask_value)
        has_ocean_neighbor[y, x] = np.any(mask[y-1:y+2, x-1:x+2] == ocean_mask_value)
        
for y in range(y_ind_min, y_ind_max):
    for x in range (x_ind_min, x_ind_max):       
        has_neighbor_with_land_neighbor[y, x] = np.any(has_land_neighbor[y-1:y+2, x-1:x+2])
        has_neighbor_with_ocean_neighbor[y, x] = np.any(has_ocean_neighbor[y-1:y+2, x-1:x+2])

# We want all cells that are within two cells of land and >2 cells from open ocean
ocean_neighbor_mask = np.logical_or(has_ocean_neighbor, has_neighbor_with_ocean_neighbor)
land_neighbor_mask = np.logical_or(has_land_neighbor, has_neighbor_with_land_neighbor)
tmp_mask = np.logical_and(1 - ocean_neighbor_mask, land_neighbor_mask)

mask_final = np.logical_and(tmp_mask, larsen_b_mask)
mask_final_tile = np.tile(mask_final, (n_time, 1, 1))

melt_mean = np.nanmean(melt.values[mask_final_tile])
melt_std = np.nanstd(melt.values[mask_final_tile])
print (f"{melt_mean:.2f} ± {melt_std:.2f} m/yr")

# Convert from m/yr ice equiv with rhoi=917 kg/m3 to kg / m2 / s
rhoi = 917.  # use density from Paolo et al. dataset
s_per_yr = 24. * 60. * 60. * 365.
melt_mean_kgm2s = melt_mean * rhoi / s_per_yr
melt_std_kgm2s = melt_std * rhoi / s_per_yr 

mesh = xr.open_dataset("Crane.nc", mode="a", engine="netcdf4")
mesh["floatingBasalMassBal"][:] = melt_mean_kgm2s * (mesh.bedTopography.values < 0.)
# Need to create this variable because it's not automatically added to the mesh.
mesh["floatingBasalMassBalUncertianty"] = (
    ("Time", "nCells"),
    melt_std_kgm2s * (mesh.bedTopography.data<0.) )
mesh.to_netcdf("Crane.nc", mode="a")

ds.close()
mesh.close()
