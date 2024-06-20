#!/bin/bash
e3sm_conda
mpas_tools=/global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools/landice/

cd $mpas_tools
git checkout crane_interpolation
cd -

echo $PWD
cp ../landice/crane/mesh_gen/mesh/Crane.nc .

# Do all edits on local copy of landsat velocity file
cp /global/cfs/cdirs/fanssie/users/trhille/data/Crane/velocities_for_trevor/landsat_both/18dec2002_20feb2003/landsat_18dec2002_20feb2003_hsc15_htc70_inc2_3031_new_hp.nc .

ncap2 -A -s "x1 = x; y1 = y" landsat_18dec2002_20feb2003_hsc15_htc70_inc2_3031_new_hp.nc

python /global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools/mesh_tools/create_SCRIP_files/create_SCRIP_file_from_MPAS_mesh.py \
    -m Crane.nc -s Crane_scrip.nc

python /global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools/mesh_tools/create_SCRIP_files/create_SCRIP_file_from_planar_rectangular_grid.py \
    -i landsat_18dec2002_20feb2003_hsc15_htc70_inc2_3031_new_hp.nc \
    -s landsat_18dec2002_20feb2003_hsc15_htc70_inc2_3031_new_hp_scrip.nc -p ais-bedmap2 -r 2

# Might need an interactive node to use ESMF_RegridWeightGen
ESMF_RegridWeightGen --source landsat_18dec2002_20feb2003_hsc15_htc70_inc2_3031_new_hp_scrip.nc \
    -d Crane_scrip.nc --weight landsat_to_MPAS_weights.nc --method conserve --netcdf4 \
    --dst_regional --src_regional --ignore_unmapped

ncap2 -A -s "observedSurfaceVelocityXLandSat = observedSurfaceVelocityX * 0.0;
             observedSurfaceVelocityYLandSat = observedSurfaceVelocityY * 0.0;
             landSatMask = observedSurfaceVelocityX * 0.0" Crane.nc

python extrap_landsat.py

python /global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools/landice/mesh_tools_li/interpolate_to_mpasli_grid.py \
    -s landsat_18dec2002_20feb2003_hsc15_htc70_inc2_3031_new_hp.nc_extrap -d Crane.nc -m e -w landsat_to_MPAS_weights.nc \
    -v observedSurfaceVelocityXLandSat observedSurfaceVelocityYLandSat landSatMask

# Cut out a subdomain around Crane for fast interpolation
ncks -O -d x,1750,1950 -d y,4045,4280 \
    /global/cfs/cdirs/fanssie/users/trhille/data/BedMachineAntarctica-v3.5/BedMachineAntarctica-v3.5_edits_floodFill_extrap.nc \
    BedMachineAntarctica-v3.5_edits_floodFill_extrap_Crane.nc

python /global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools//mesh_tools/create_SCRIP_files/create_SCRIP_file_from_planar_rectangular_grid.py \
    -i BedMachineAntarctica-v3.5_edits_floodFill_extrap_Crane.nc \
    -s BedMachineAntarctica-v3.5_edits_floodFill_extrap_Crane_scrip.nc -p ais-bedmap2 -r 2

ESMF_RegridWeightGen --source BedMachineAntarctica-v3.5_edits_floodFill_extrap_Crane_scrip.nc \
    -d Crane_scrip.nc --weight BedMachine_to_MPAS_weights.nc --method conserve \
    --netcdf4 --dst_regional --src_regional --ignore_unmapped

ncap2 -A -s "surfaceBedMachine = bedTopography * 0.0" Crane.nc

python /global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools/landice/mesh_tools_li/interpolate_to_mpasli_grid.py \
    -s BedMachineAntarctica-v3.5_edits_floodFill_extrap_Crane.nc -d Crane.nc -m e -w BedMachine_to_MPAS_weights.nc

ncap2 -A -s "where(iceMask < 0.5) thickness = 0.0" Crane.nc

# Reproject Berthier DEM to EPSG 3031
# TODO: Berthier DEM uses elevation wrt EGM96 geoid, while bedmachine uses EIGEN-6C4
conda activate gdal

gdalwarp -t_srs EPSG:3031 /global/cfs/cdirs/fanssie/users/trhille/data/Crane/Berthier_DEMs/DEM_2002.85300894_shifted_H_V.tif \
    DEM_2002.85300894_shifted_H_V_EPSG3031.tif

gdal_translate -of netCDF DEM_2002.85300894_shifted_H_V_EPSG3031.tif DEM_2002.85300894_shifted_H_V_EPSG3031.nc

# Trim and reproject ASTER DEM to EPSG 3031
gdalwarp -t_srs EPSG:3031 /global/cfs/cdirs/fanssie/users/trhille/data/Crane/APDEM100m.tif APDEM100m_EPSG3031.tif

gdal_translate -of netCDF APDEM100m_EPSG3031.tif APDEM100m_EPSG3031.nc

ncks -O -d x,1000,2000 -d y,3800,5000 APDEM100m_EPSG3031.nc APDEM100m_EPSG3031_Crane.nc

# Create the 2002 ice mask using a pre-defined geojson file. This could maybe be improved,
# but using the gappy datasets to define the ice extent is also problematic.
python create_ice_mask_from_geojson.py

conda deactivate  # no longer need gdal

# Remap higher-resolution Berthier DEM onto Aster DEM grid
python /global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools/mesh_tools/create_SCRIP_files/create_SCRIP_file_from_planar_rectangular_grid.py \
    -i DEM_2002.85300894_shifted_H_V_EPSG3031.nc \
    -s DEM_2002.85300894_shifted_H_V_EPSG3031_scrip.nc -p ais-bedmap2 -r 2

python /global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools//mesh_tools/create_SCRIP_files/create_SCRIP_file_from_planar_rectangular_grid.py \
    -i APDEM100m_EPSG3031_Crane.nc \
    -s APDEM100m_EPSG3031_Crane.scrip.nc -p ais-bedmap2 -r 2

ESMF_RegridWeightGen --source DEM_2002.85300894_shifted_H_V_EPSG3031_scrip.nc \
    -d APDEM100m_EPSG3031_Crane.scrip.nc --weight Berthier_to_ASTER_weights.nc \
    --method conserve --netcdf4 --dst_regional --src_regional --ignore_unmapped

ncrename -v Band1,surfaceBerthier DEM_2002.85300894_shifted_H_V_EPSG3031.nc
ncrename -v Band1,surfaceAster APDEM100m_EPSG3031_Crane.nc

# Flood-fill Berthier DEM to remove noise and icebergs? In theory
# it seems like a good idea, but it ends up with a really bad surface
# on the MALI mesh.
# python dem_flood_fill.py

# Now remap Berthier onto ASTER mesh (finer onto coarser) before
# combining and interpolating to MALI mesh.
ncremap -m Berthier_to_ASTER_weights.nc DEM_2002.85300894_shifted_H_V_EPSG3031.nc Berthier_ASTER_composite_DEM.nc
ncks -A -v surfaceAster APDEM100m_EPSG3031_Crane.nc Berthier_ASTER_composite_DEM.nc

# Use Berthier data where valid, elswhere use ASTER
ncap2 -A -s "surfaceBerthierAster = surfaceBerthier" Berthier_ASTER_composite_DEM.nc
ncap2 -A -s "where(surfaceBerthier > 0.0 && surfaceBerthier < 5.0e3) surfaceBerthierAster = surfaceBerthier;
             elsewhere surfaceBerthierAster = surfaceAster;
             where(surfaceBerthierAster < 0.0) surfaceBerthierAster = 0.0" Berthier_ASTER_composite_DEM.nc
ncap2 -A -s "x1 = x; y1 = y"  Berthier_ASTER_composite_DEM.nc

# Now we want to interpolate the composite DEM onto the MALI mesh
ESMF_RegridWeightGen --source APDEM100m_EPSG3031_Crane.scrip.nc \
    -d Crane_scrip.nc --weight DEM_to_MPAS_weights.nc --method conserve \
    --netcdf4 --dst_regional --src_regional --ignore_unmapped

ncap2 -A -s "surfaceBerthierAster = bedTopography * 0.0" Crane.nc

# Interpolate composite DEM surface 
# TODO: Determine if these are relative to ellipsoid or geoid
python /global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools/landice/mesh_tools_li/interpolate_to_mpasli_grid.py \
    -s Berthier_ASTER_composite_DEM.nc -d Crane.nc -m e -w DEM_to_MPAS_weights.nc -v surfaceBerthierAster

# Use MEaSUREs 1995-2001 and 2007-2009 velocities to place bounds on velocities upstream.
# We assume an upper bound of the 2007-2009 velocity + std and a lower bound of the 1995-2001
# velocity - std. The resulting "observed" velocity is the cell-wise average of these, and the
# "uncertainty" is half the spread.
cp /global/cfs/cdirs/piscees/trhille/AIS_datasets/NSIDC_0761_measures_velocities/antarctica_ice_velocity_1995-2001_450m_v01.0_Crane.nc .
ncrename -v VX,VX2001 -v VY,VY2001 -v STDX,STDX2001 -v STDY,STDY2001 -v x,x1 -v y,y1 antarctica_ice_velocity_1995-2001_450m_v01.0_Crane.nc
ncks -A -v VX2001,VY2001,STDX2001,STDY2001,x1,y1 antarctica_ice_velocity_1995-2001_450m_v01.0_Crane.nc measures_upstream_velocities.nc

cp /global/cfs/cdirs/piscees/trhille/AIS_datasets/NSIDC_0761_measures_velocities/antarctica_ice_velocity_2007-2009_450m_v01.0_Crane.nc .
ncrename -v VX,VX2009 -v VY,VY2009 -v STDX,STDX2009 -v STDY,STDY2009 -v x,x1 -v y,y1 antarctica_ice_velocity_2007-2009_450m_v01.0_Crane.nc
ncks -A -v VX2009,VY2009,STDX2009,STDY2009 antarctica_ice_velocity_2007-2009_450m_v01.0_Crane.nc measures_upstream_velocities.nc

ncap2 -A -s "lower_vel_bound_x = VX2001 - STDX2001; lower_vel_bound_y = VY2001 - STDY2001;
             upper_vel_bound_x = VX2009 + STDX2009; upper_vel_bound_y = VY2009 + STDY2009;
             vx = (VX2001 + VX2009) / 2.0;
             vy = (VY2001 + VY2009) / 2.0;
             vErr = (sqrt(upper_vel_bound_x^2.0 + upper_vel_bound_y^2.0) -
                     sqrt(lower_vel_bound_x^2.0 + lower_vel_bound_y^2.0)).abs() / 2.0" \
             measures_upstream_velocities.nc

python /global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools/mesh_tools/create_SCRIP_files/create_SCRIP_file_from_planar_rectangular_grid.py \
    -i measures_upstream_velocities.nc \
    -s measures_upstream_velocities.scrip.nc \
    -p ais-bedmap2 -r 2

ESMF_RegridWeightGen --source measures_upstream_velocities.scrip.nc \
    -d Crane_scrip.nc --weight measures_to_MPAS_weights.nc --method conserve --netcdf4 \
    --dst_regional --src_regional --ignore_unmapped

ncap2 -A -s "observedSurfaceVelocityXMeasures = observedSurfaceVelocityX * 0.0;
             observedSurfaceVelocityYMeasures = observedSurfaceVelocityY * 0.0;
             observedSurfaceVelocityUncertaintyMeasures = observedSurfaceVelocityUncertainty * 0.0" Crane.nc

python /global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools/landice/mesh_tools_li/interpolate_to_mpasli_grid.py \
    -s measures_upstream_velocities.nc \
    -d Crane.nc -m e -w measures_to_MPAS_weights.nc -v observedSurfaceVelocityXMeasures observedSurfaceVelocityYMeasures observedSurfaceVelocityUncertaintyMeasures

# Where LandSat velocities are good, use those (and assume 20% uncertainty).
# Everywhere else, use MEaSUREs data.
ncap2 -A -s "where(landSatMask>0.99999) {
                 observedSurfaceVelocityX = observedSurfaceVelocityXLandSat;
                 observedSurfaceVelocityY = observedSurfaceVelocityYLandSat;
                 observedSurfaceVelocityUncertainty = 0.2 * sqrt(observedSurfaceVelocityXLandSat^2.0 + observedSurfaceVelocityYLandSat^2.0);
              } elsewhere{
                 observedSurfaceVelocityX = observedSurfaceVelocityXMeasures;
                 observedSurfaceVelocityY = observedSurfaceVelocityYMeasures;
                 observedSurfaceVelocityUncertainty = observedSurfaceVelocityUncertaintyMeasures;
              }" Crane.nc

ncap2 -A -s "where(observedSurfaceVelocityUncertainty == 0.0) observedSurfaceVelocityUncertainty = 1.0" Crane.nc
ncap2 -A -s "where(sqrt((observedSurfaceVelocityX^2.0 + observedSurfaceVelocityY^2.0) == 0.0) observedSurfaceVelocityUncertainty = 1.0" Crane.nc
# Get subset of ITS_LIVE dHdt data:
ncks -O -d x,175,225 -d y,775,850 /global/cfs/cdirs/fanssie/users/trhille/data/ITS_LIVE_dHdt/ITS_LIVE_dHdt.nc ITS_LIVE_dHdt_Crane.nc
# ITS_LIVE data are scaled, so we need to unpack:
ncap2 -A -s "dh_unpacked = float(dh) * dh@scale_factor"  ITS_LIVE_dHdt_Crane.nc
# Berthier DEM is from Nov 2002, velocities are from Dec 2002 to Feb 2003, so let's use dH/dt from Oct 2002 to March 2003.
# Determining those dates using ncview rather than something more reproducible gives time indices 210 to 215.
t_start=210
t_end=215
ncap2 -A -s "dt = time * 0.0; dt(1:-1) = time(1:-1) - time(0:-2)" ITS_LIVE_dHdt_Crane.nc
ncap2 -A -s "dH = dh_unpacked * 0.0" ITS_LIVE_dHdt_Crane.nc
ncap2 -A -s "dH(1:-1,:,:) = dh_unpacked(1:-1,:,:) - dh_unpacked(0:-2,:,:)" ITS_LIVE_dHdt_Crane.nc
ncatted -a units,dt,o,c,"days since previous time" ITS_LIVE_dHdt_Crane.nc
ncatted -a standard_name,dH,o,c,"height anomaly w.r.t. previous time" ITS_LIVE_dHdt_Crane.nc
ncatted -a units,dt,o,c,"days" ITS_LIVE_dHdt_Crane.nc
ncap2 -A -s "dHdt = dH / dt * 365.0" ITS_LIVE_dHdt_Crane.nc
ncatted -a standard_name,dHdt,o,c,"surface height change rate" ITS_LIVE_dHdt_Crane.nc
ncatted -a units,dHdt,o,c,"m/yr" ITS_LIVE_dHdt_Crane.nc
ncks -O -d time,$t_start,$t_end ITS_LIVE_dHdt_Crane.nc ITS_LIVE_dHdt_Crane_Oct02_March03_mean.nc
ncap2 -A -s 'dHdt_mean = dH.sum($time) / dt.sum($time) * 365.0' ITS_LIVE_dHdt_Crane_Oct02_March03_mean.nc
ncatted -a standard_name,dHdt_mean,o,c,"mean surface height change rate over time=$t_start:$t_end" ITS_LIVE_dHdt_Crane_Oct02_March03_mean.nc
ncatted -a units,dHdt_mean,o,c,"m/yr" ITS_LIVE_dHdt_Crane_Oct02_March03_mean.nc
ncap2 -A -s 'dHdtErr = dHdt.rmssdn($time)' ITS_LIVE_dHdt_Crane_Oct02_March03_mean.nc
ncatted -a standard_name,dHdtErr,o,c,"standard deviation of dHdt over averaging period" ITS_LIVE_dHdt_Crane_Oct02_March03_mean.nc
ncatted -a units,dHdtErr,o,c,"m/yr" ITS_LIVE_dHdt_Crane_Oct02_March03_mean.nc
ncap2 -A -s "x1 = x; y1 = y" ITS_LIVE_dHdt_Crane_Oct02_March03_mean.nc
ncatted -a _FillValue,dHdt_mean,o,f,0.0 ITS_LIVE_dHdt_Crane_Oct02_March03_mean.nc
ncatted -a _FillValue,dHdtErr,o,f,3.154e7 ITS_LIVE_dHdt_Crane_Oct02_March03_mean.nc

# Interpolate dHdt to MALI mesh
# TODO: How to treat missing data over lower glacier?
python /global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools/mesh_tools/create_SCRIP_files/create_SCRIP_file_from_planar_rectangular_grid.py \
    -i ITS_LIVE_dHdt_Crane_Oct02_March03_mean.nc \
    -s ITS_LIVE_dHdt_Crane_Oct02_March03_mean.scrip.nc \
    -p ais-bedmap2 -r 2

ESMF_RegridWeightGen --source ITS_LIVE_dHdt_Crane_Oct02_March03_mean.scrip.nc \
    -d Crane_scrip.nc --weight ITS_LIVE_to_MPAS_weights.nc \
    --method bilinear --extrap_method nearestd --netcdf4 --dst_regional --src_regional --ignore_unmapped

python /global/cfs/cdirs/fanssie/users/trhille/MPAS-Tools/landice/mesh_tools_li/interpolate_to_mpasli_grid.py \
    -s ITS_LIVE_dHdt_Crane_Oct02_March03_mean.nc \
    -d Crane.nc -m e -w ITS_LIVE_to_MPAS_weights.nc \
    -v observedThicknessTendency observedThicknessTendencyUncertainty

# Remap SMB, using 2002-2003 average to be consistent with surface, velocity, and dH/dt observations.
# Testing shows that there is not a large difference compared with using the 1995-2017 average.
# This is a hacky use of the COMPASS workflow, but I can't think of a better way to use the existing tools.
cd /global/cfs/cdirs/fanssie/users/trhille/compass/
git checkout trhille/landice/crane_smb
source /global/cfs/cdirs/fanssie/users/trhille/compass/load_dev_compass_1.4.0-alpha.2_pm-cpu_gnu_mpich_albany.sh
mkdir process_RACMO
sed -i -e "s/17,39/21,26/g" \
    /global/cfs/cdirs/fanssie/users/trhille/compass/compass/landice/tests/ismip6_forcing/atmosphere/process_smb_racmo.py
sed -i -e "s/1995.2017/2000-2005/g" \
    /global/cfs/cdirs/fanssie/users/trhille/compass/compass/landice/tests/ismip6_forcing/atmosphere/process_smb_racmo.py
compass setup -t landice/ismip6_forcing/atmosphere -w process_RACMO -f ismip6_forcing.cfg
cd process_RACMO/landice/ismip6_forcing/atmosphere/process_smb_racmo/
compass run
cd -
ncks -A -v sfcMassBal,sfcMassBalUncertainty \
    process_RACMO/landice/ismip6_forcing/atmosphere/process_smb_racmo/Crane_RACMO2.3p2_ANT27_smb_climatology_2000-2005.nc \
    Crane.nc

# Calculate ice thickness from surface and bed topography
ncap2 -O -s "rhoi=910.0; rhosw=1028.0; 
             float_surface = (-rhosw / rhoi * bedTopography) *
                 ((-rhosw / rhoi * bedTopography) > 0.0) + bedTopography" Crane.nc tmp.nc
ncap2 -A -s "float_mask = (surfaceBerthierAster < float_surface) * (bedTopography < 0)" tmp.nc
ncap2 -A -s "thickness = float_mask * ( 1.0 / (1.0 - rhoi/rhosw) *
                 surfaceBerthierAster) + (1.0 - float_mask) *
                 (surfaceBerthierAster - bedTopography)" tmp.nc
ncks -A -v thickness tmp.nc Crane.nc

# Now trim ice extent to match edge of velocities dataset, which was manually defined with a geojson file.
# This was calculated on the MALI mesh above.
ncap2 -A -s "thickness = thickness * iceMask_2002" Crane.nc

cp /global/cfs/cdirs/fanssie/users/trhille/data/Paolo_2023_melt_rates/ANT_G1920V01_IceShelfMelt.nc .
python calculate_floating_basal_mass_bal.py 
# To-do: 1) Correct thickness for firn air content?;
