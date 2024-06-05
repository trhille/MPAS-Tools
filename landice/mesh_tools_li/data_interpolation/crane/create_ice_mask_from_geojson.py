#!/usr/bin/env python3

import geopandas as gpd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from shapely.geometry import MultiPoint, Point


mesh_path = "Crane.nc"
mesh = Dataset(mesh_path, "r+")
mesh.set_auto_mask(False)

shp_path = "Crane_ice_extent.geojson"
shp = gpd.read_file(shp_path)

x_cell = mesh.variables["xCell"][:]
y_cell = mesh.variables["yCell"][:]

shp_xy = shp.to_crs('EPSG:3031')
mesh_xy = gpd.GeoDataFrame(
    index=range(len(x_cell)), crs="EPSG:3031",
    geometry=[Point([x, y]) for x, y in zip(x_cell, y_cell)])

ice_mask = np.zeros_like(x_cell, dtype='int')
for ii, (x, y) in enumerate(zip(x_cell, y_cell)):
    ice_mask[ii] = shp_xy.contains(Point([x, y]))[0]

print(f"contains {np.sum(ice_mask)} out of {len(x_cell)} cells")

if "iceMask_2002" not in mesh.variables.keys():
    mesh.createVariable("iceMask_2002", "i", dimensions=("nCells"))
mesh.variables["iceMask_2002"][:] = ice_mask
mesh.close()
