#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.warp import reproject, Resampling, transform_bounds, calculate_default_transform
from rasterio.windows import from_bounds
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator
from scipy.ndimage import distance_transform_edt
import subprocess


high_quality_path = "/global/cfs/cdirs/m4288/users/trhille/data/Crane/Rebesco_bathymetry/LarsenB_bathy.tif"
low_quality_path = "/global/cfs/cdirs/m4288/users/trhille/data/Crane/Huss_Farinotti_tc-8-1261-2014-supplement/SOM/bedrock/bedrock.tif"
N = 20  # width of the blending zone in pixels
output_path = f"Rebesco_HF14_BedMachinev3_merged_buffered_feathering_N{N}_extended"

bedmachine_path = "/global/cfs/cdirs/m4288/users/trhille/data/BedMachineAntarctica-v3/"
bedmachine_file = "BedMachineAntarctica-v3"  # leave out .nc

# Buffered feathering approach:
# 1. Create N pixel-side buffer around high quality DEM
# 2. NN interpolate into the buffer
# 3. In that "gutter" region around the high quality DEM, 
#    have blending weight vary from 0 at the outer edge to 1
#    at the inner edge (i.e., at the actual high quality DEM boundary)
# 4. Blend the two DEMs across that N pixel region using that 
#    weighting scheme, so that we get a fully smooth transition 
#    over N pixels, and we use the high quality data exactly
#    as-is for the areas where it exists.

# It would be much better to trim the DEMs to a region
# around Crane Glacier, but for some reason that has been 
# very difficult to make work correctly. So, for now we do
# all operations on the full data set and then trim it using
# ncks at the end. Not ideal, not fast, but it works.
minx, miny, maxx, maxy = -2637500, 1205700, -2384100, 1285900
target_bounds_orig = (minx, miny, maxx, maxy)
low_crs = CRS.from_epsg(3031)

# Load lower-res DEM:
with rasterio.open(low_quality_path) as src_low:
    low_subset = src_low.read(1)
    low_transform = src_low.transform
    low_bounds = src_low.bounds

# Load high-quality DEM
with rasterio.open(high_quality_path) as src_high:
    fill = -9999
    high_crs = src_high.crs
    target_bounds = transform_bounds(low_crs, high_crs, *low_bounds)
    window = from_bounds(*target_bounds, transform=src_high.transform)
    high_dem = src_high.read(1, window=window,
                             boundless=True,
                             masked=False,
                             fill_value=fill)
    high_profile = src_high.profile
    high_transform = src_high.window_transform(window)
    #high_transform = src_high.transform
    high_bounds = src_high.bounds

# Reproject low-quality subset
reprojected_low = np.zeros_like(high_dem, dtype=np.float32)
reproject(
    source=low_subset,
    destination=reprojected_low,
    src_transform=low_transform,
    src_crs=low_crs,
    dst_transform=high_transform,
    dst_crs=high_crs,
    resampling=Resampling.bilinear,
    src_nodata=-9999,
    dst_nodata=np.nan
)

# Clean and correct
high_dem[high_dem == -9999] = np.nan
reprojected_low[reprojected_low == -9999] = np.nan

high_mask = ~np.isnan(high_dem)

# Nearest neighbor extrapolation
_, nearest_idx = distance_transform_edt(~high_mask, return_indices=True)
high_extrapolated = high_dem[tuple(nearest_idx)]

# Distance to nearest high-quality pixel
dist_to_high = distance_transform_edt(~high_mask)

# Gutter mask: only in a ring N pixels wide
gutter_mask = (dist_to_high > 0) & (dist_to_high <= N)

# Linear blending weight: 1 at edge of high_mask, 0 at outer edge of gutter
blend_weight = np.zeros_like(high_dem, dtype=np.float32)
blend_weight[gutter_mask] = 1 - (dist_to_high[gutter_mask] / N)

# Blending only in gutter region
blended = blend_weight * high_extrapolated + (1 - blend_weight) * reprojected_low

# Final merged surface
merged = np.full_like(high_dem, np.nan)
merged[high_mask] = high_dem[high_mask]  # use original high-quality DEM
merged[gutter_mask] = blended[gutter_mask]  # blend in gutter
merged[~high_mask & (dist_to_high > N)] = reprojected_low[~high_mask & (dist_to_high > N)]  # far from high DEM
merged[np.logical_and(reprojected_low > 0., ~high_mask)] = reprojected_low[np.logical_and(reprojected_low > 0., ~high_mask)]  # only apply merging below sea level.

subprocess.run(['ncks', '-O', '-d', 'x,1750,1950', '-d', 'y,4045,4280', f"{bedmachine_path}"f"{bedmachine_file}.nc", f"{bedmachine_file}_Crane.nc"])

subdataset_path = f"NETCDF:{bedmachine_file}_Crane.nc:bed"

with rasterio.open(subdataset_path) as src_bed:
    bed_data = src_bed.read(1).astype(np.float32)
    bed_data[bed_data <= -9998] = np.nan  # Clean BedMachine nodata
    bed_transform = src_bed.transform
    bed_crs = src_bed.crs

filled_bedmachine = np.full_like(merged, np.nan, dtype=np.float32)

reproject(
    source=bed_data,
    destination=filled_bedmachine,
    src_transform=bed_transform,
    src_crs=bed_crs,
    dst_transform=high_transform,
    dst_crs=high_crs,
    resampling=Resampling.bilinear,
    src_nodata=np.nan,
    dst_nodata=np.nan
)

filled_merged = merged.copy()
gap_mask = np.isnan(filled_merged)
filled_merged[gap_mask] = filled_bedmachine[gap_mask]

# Now let's put this all back into EPSG 3031
bounds = rasterio.transform.array_bounds(
    merged.shape[0], merged.shape[1], high_transform
)
target_crs = CRS.from_epsg(3031)
transform_3031, width_3031, height_3031 = calculate_default_transform(
    high_crs,         # CRS of your current DEM
    target_crs,
    filled_merged.shape[1],
    filled_merged.shape[0],
    *bounds   # affine transform of original
)

filled_merged_reprojected = np.empty((height_3031, width_3031), dtype=np.float32)

reproject(
    source=filled_merged,
    destination=filled_merged_reprojected,
    src_transform=high_transform,
    src_crs=high_crs,
    dst_transform=transform_3031,
    dst_crs=target_crs,
    resampling=Resampling.bilinear,  # or Resampling.nearest if categorical
    src_nodata=np.nan,
    dst_nodata=np.nan
)

profile = high_profile.copy()
profile.update({
    'crs': target_crs,
    'transform': transform_3031,
    'width': width_3031,
    'height': height_3031,
    'dtype': 'float32',
    'nodata': np.nan
})

with rasterio.open(f'{output_path}.tif', "w", **profile) as dst:
    dst.write(filled_merged_reprojected, 1)

subprocess.run(['gdal_translate', '-of', 'netCDF', f'{output_path}.tif', f'{output_path}.nc'])

subprocess.run(['ncks', '-O', '-d', 'x,18000,22300', '-d', 'y,25000,29000', f'{output_path}.nc', f'{output_path}_small.nc'])
#!rm Rebesco_HF14_BedMachinev3.5_merged_buffered_feathering_N{N}_extended.nc

subprocess.run(['ncrename', '-v', 'Band1,topg', f'{output_path}_small.nc'])

subprocess.run(['ncap2', '-A', '-s', '"x1=x; y1=y"', f'{output_path}_small.nc'])
