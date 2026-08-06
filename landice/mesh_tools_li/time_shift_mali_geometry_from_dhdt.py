#!/usr/bin/env python
'''
Time-shift the ice thickness of a MALI mesh between two years using observed
surface-elevation-change datasets from the NASA MEaSUREs ITS_LIVE project.
The mesh can be shifted backward in time (rewound) or forward in time
(fast-forwarded) depending on the choice of start and end years.

Two complementary datasets are used, both provided on the same EPSG:3031
Antarctic Polar Stereographic grid (2916 x 2916, 1920 m posting):

  * Grounded ice (NSIDC-0782, ANT_G1920_GroundedIceHeight_v01.nc)
      variable ``height_change(time, y, x)`` = change in ice-sheet surface
      elevation [m] relative to a fixed reference date.  For grounded ice the
      bed is assumed fixed, so surface-elevation change equals thickness
      change.

  * Floating ice (NSIDC-0792, NSIDC-0792_19920317-20171216_V01.0.nc)
      variable ``height_change(time, y, x)`` = change in ice-shelf surface
      (freeboard) elevation [m].  Assuming hydrostatic equilibrium, the
      corresponding thickness change is the surface change multiplied by the
      flotation factor rhosw / (rhosw - rhoi).

The thickness change applied to the mesh is the cumulative height change at the
requested ``--endYear`` minus that at ``--startYear`` (each sampled at the
nearest available time step in the respective dataset):

    dThickness = height_change(endYear) - height_change(startYear)

To shift a mesh that represents ``startYear`` (e.g. 2015) back to an earlier
``endYear`` (e.g. 2000), this difference is negative where the ice thinned over
that interval, so the resulting mesh is correspondingly thicker in the past.
The same tool can also advance a mesh forward in time simply by choosing
``endYear > startYear``.

Because the grounded and floating datasets do not cover exactly the same
footprint, small gaps (notably near grounding lines) can occur between them.
These gaps are filled by linear interpolation from surrounding valid mesh
cells.  By default the thickness change is applied only to cells that already
contain ice; ice-free cells are left unchanged so that interpolation artifacts
do not appear in the open ocean or on bare land.  Pass ``--allowExtentChange``
to let the ice extent change between the two years: new ice may then appear,
but only where a dataset directly measured a change (the observed ice
footprint) and where the resulting thickness is positive, so the result is
clipped to the ice extent at the target time and artifacts cannot spread to
the mesh edge.  Mesh cells that fall entirely outside the combined data
coverage are left unchanged.  Resulting thicknesses are clamped to be
non-negative.

The input mesh is copied to a new file (``*_shifted_<endYear>.nc`` by default)
before its ``thickness`` field is overwritten; the original mesh is not
modified.  A diagnostic integer field ``dThicknessSource`` is also written to
the output, flagging for each cell how its thickness change was obtained
(0 = unchanged, 1 = directly measured, 2 = gap-filled by interpolation).

NOTE: the MALI mesh must use the same EPSG:3031 projected coordinates as the
ITS_LIVE datasets (xCell / yCell in metres).  A warning is printed if the mesh
cell centres fall largely outside the dataset extent.

Trevor Hillebrand, 2026
'''

import shutil
import numpy as np
import xarray as xr
import cftime
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator, griddata
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser(
        prog='time_shift_mali_geometry_from_dhdt.py',
        description=__doc__)
    parser.add_argument("-m", "--mesh", dest="meshFile", required=True,
                        metavar="FILENAME",
                        help="MALI mesh file to time-shift (EPSG:3031 "
                             "coordinates)")
    parser.add_argument("-g", "--grounded", dest="groundedFile", required=True,
                        metavar="FILENAME",
                        help="Grounded-ice height-change dataset "
                             "(ANT_G1920_GroundedIceHeight_v01.nc)")
    parser.add_argument("-f", "--floating", dest="floatingFile", required=True,
                        metavar="FILENAME",
                        help="Floating-ice height-change dataset "
                             "(NSIDC-0792_19920317-20171216_V01.0.nc)")
    parser.add_argument("-s", "--startYear", dest="startYear", type=float,
                        required=True,
                        help="Nominal year the mesh currently represents "
                             "(e.g. 2015)")
    parser.add_argument("-e", "--endYear", dest="endYear", type=float,
                        required=True,
                        help="Year to shift the mesh to (e.g. 2000)")
    parser.add_argument("-o", "--outFile", dest="outFile", default=None,
                        metavar="FILENAME",
                        help="Output mesh filename. Default: input name with "
                             "'_shifted_<endYear>.nc' suffix")
    parser.add_argument("--rhoi", dest="rhoi", type=float, default=910.0,
                        help="Ice density [kg m-3] used for the floating "
                             "flotation factor (default: 910)")
    parser.add_argument("--rhosw", dest="rhosw", type=float, default=1028.0,
                        help="Sea-water density [kg m-3] used for the floating "
                             "flotation factor (default: 1028)")
    parser.add_argument("--allowExtentChange", dest="allowExtentChange",
                        action="store_true",
                        help="Allow the ice extent to change between the two "
                             "years. New ice may appear only where a dataset "
                             "directly measured a change and where the "
                             "resulting thickness is positive, clipping the "
                             "result to the ice extent at the target time. By "
                             "default only cells that currently contain ice "
                             "are modified.")
    return parser.parse_args()


def nearest_time_index(ds, year):
    '''Return the index of the time step nearest to a (fractional) calendar
    year, along with the actual fractional year of that step.'''
    units = ds['time'].attrs.get('units', 'days since 1950-01-01')
    calendar = ds['time'].attrs.get('calendar', 'standard')
    dates = cftime.num2date(ds['time'].values, units, calendar=calendar)
    years = np.array([d.year + (d.dayofyr - 1) / 365.0 for d in dates])
    idx = int(np.argmin(np.abs(years - year)))
    return idx, years[idx]


def height_change_difference(ds, startYear, endYear):
    '''Cumulative height-change difference height_change(endYear) -
    height_change(startYear), sampled at the nearest time step to each year.
    Returned as a 2D array (y, x) with NaN where either slice is missing.'''
    iStart, yStart = nearest_time_index(ds, startYear)
    iEnd, yEnd = nearest_time_index(ds, endYear)
    hcStart = ds['height_change'].isel(time=iStart).values.astype(np.float64)
    hcEnd = ds['height_change'].isel(time=iEnd).values.astype(np.float64)
    print(f"    startYear {startYear:g} -> nearest step {iStart} "
          f"({yStart:.2f}); endYear {endYear:g} -> nearest step {iEnd} "
          f"({yEnd:.2f})")
    return hcEnd - hcStart


def interp_grid_to_mesh(field2d, xGrid, yGrid, xCell, yCell):
    '''Bilinearly interpolate a 2D grid field (dims y, x) onto mesh cell
    centres.  Handles a descending y coordinate.  Returns NaN for mesh cells
    outside the grid or adjacent to missing (NaN) grid values.'''
    # RegularGridInterpolator requires strictly increasing coordinates.
    if yGrid[0] > yGrid[-1]:
        yGrid = yGrid[::-1]
        field2d = field2d[::-1, :]
    if xGrid[0] > xGrid[-1]:
        xGrid = xGrid[::-1]
        field2d = field2d[:, ::-1]
    interp = RegularGridInterpolator(
        (yGrid, xGrid), field2d,
        method='linear', bounds_error=False, fill_value=np.nan)
    return interp(np.column_stack((yCell, xCell)))


def main():
    args = parse_args()
    flotationFactor = args.rhosw / (args.rhosw - args.rhoi)

    outFile = args.outFile
    if outFile is None:
        outFile = f"{args.meshFile.rsplit('.nc', 1)[0]}" \
                  f"_shifted_{args.endYear:g}.nc"
    print(f"Copying '{args.meshFile}' to '{outFile}'")
    shutil.copy(args.meshFile, outFile)

    # --- Read mesh geometry ---
    mesh = Dataset(outFile, 'r+')
    xCell = mesh.variables['xCell'][:]
    yCell = mesh.variables['yCell'][:]
    thickness = mesh.variables['thickness'][0, :].astype(np.float64)
    nCells = len(xCell)

    # --- Grounded dataset: thickness change = surface-elevation change ---
    print("Reading grounded-ice dataset:")
    with xr.open_dataset(args.groundedFile, decode_times=False) as gds:
        xGrid = gds['x'].values
        yGrid = gds['y'].values
        dThkGrounded = height_change_difference(
            gds, args.startYear, args.endYear)
    dGrounded = interp_grid_to_mesh(dThkGrounded, xGrid, yGrid, xCell, yCell)

    # Sanity check that the mesh overlaps the dataset extent.
    inExtent = ((xCell >= xGrid.min()) & (xCell <= xGrid.max()) &
                (yCell >= yGrid.min()) & (yCell <= yGrid.max()))
    if inExtent.mean() < 0.5:
        print("WARNING: fewer than half of the mesh cell centres fall within "
              "the ITS_LIVE grid extent. Check that the mesh uses EPSG:3031 "
              "(Antarctic Polar Stereographic) coordinates in metres.")

    # --- Floating dataset: apply flotation factor to surface change ---
    print("Reading floating-ice dataset:")
    with xr.open_dataset(args.floatingFile, decode_times=False) as fds:
        xGridF = fds['x'].values
        yGridF = fds['y'].values
        dSurfFloating = height_change_difference(
            fds, args.startYear, args.endYear)
    dThkFloating = dSurfFloating * flotationFactor
    print(f"    applied flotation factor rhosw/(rhosw-rhoi) = "
          f"{flotationFactor:.3f}")
    dFloating = interp_grid_to_mesh(dThkFloating, xGridF, yGridF, xCell, yCell)

    # --- Combine (grounded takes precedence; footprints do not overlap) ---
    # directMask flags cells where at least one dataset *directly* provides a
    # height change (i.e. the cell falls within the observed ice footprint,
    # not an extrapolated value).
    directMask = np.isfinite(dGrounded) | np.isfinite(dFloating)
    dThk = np.where(np.isfinite(dGrounded), dGrounded, dFloating)

    nCovered = directMask.sum()
    print(f"Thickness change defined on {nCovered} of {nCells} mesh cells "
          f"before gap filling")

    # --- Determine which cells may be modified ---
    # By default only cells that currently contain ice are changed, which
    # avoids interpolation artifacts in ice-free ocean/bare-land cells.
    # With --allowExtentChange the ice extent is permitted to change, but new
    # ice may appear only where a dataset directly measured a change (the
    # observed ice footprint); combined with the non-negative clamp on
    # (thickness + dThk) below, this clips the result to the ice extent at the
    # target time and prevents artifacts reaching the mesh edge.
    iceMask = thickness > 0.0
    if args.allowExtentChange:
        modifyMask = inExtent & (iceMask | directMask)
    else:
        modifyMask = inExtent & iceMask

    # --- Fill gaps between the two datasets by linear interpolation from
    #     surrounding directly-measured mesh cells (within the convex hull
    #     only). Only gaps inside modifyMask are filled, so extrapolation
    #     cannot spread beyond the observed/existing ice. ---
    gapFilled = np.zeros(nCells, dtype=bool)
    gaps = (~directMask) & modifyMask
    if gaps.any() and directMask.sum() >= 3:
        filled = griddata(
            np.column_stack((xCell[directMask], yCell[directMask])),
            dThk[directMask],
            np.column_stack((xCell[gaps], yCell[gaps])), method='linear')
        dThk[gaps] = filled
        gapFilled[gaps] = np.isfinite(filled)
        nFilled = gapFilled.sum()
        print(f"Filled {nFilled} gap cells by linear interpolation")

    # Cells without a value (outside coverage) are left unchanged.
    dThk = np.where(np.isfinite(dThk), dThk, 0.0)

    # Restrict the change to the cells we are allowed to modify (removes
    # interpolation artifacts outside the ice).
    dThk = np.where(modifyMask, dThk, 0.0)

    # --- Classify the source of each cell's thickness change (diagnostic) ---
    # 0 = unchanged, 1 = directly measured, 2 = gap-filled by interpolation.
    dThkSource = np.zeros(nCells, dtype=np.int32)
    dThkSource[directMask & modifyMask] = 1
    dThkSource[gapFilled & modifyMask] = 2

    # --- Apply and write ---
    newThickness = np.maximum(thickness + dThk, 0.0)
    nChanged = np.count_nonzero(dThk != 0.0)
    print(f"Updated thickness on {nChanged} cells "
          f"(min dThk {dThk.min():.2f} m, max dThk {dThk.max():.2f} m)")

    mesh.variables['thickness'][0, :] = newThickness

    # Write the diagnostic source field, creating it if necessary.
    if 'dThicknessSource' not in mesh.variables:
        src = mesh.createVariable('dThicknessSource', 'i4', ('Time', 'nCells'))
        src.long_name = ('source of the applied thickness change from '
                         'time_shift_mali_geometry_from_dhdt.py')
        src.flag_values = np.array([0, 1, 2], dtype=np.int32)
        src.flag_meanings = 'unchanged directly_measured gap_filled'
    mesh.variables['dThicknessSource'][0, :] = dThkSource
    print(f"Wrote dThicknessSource diagnostic field "
          f"(directly_measured {np.count_nonzero(dThkSource == 1)}, "
          f"gap_filled {np.count_nonzero(dThkSource == 2)})")

    mesh.close()
    print(f"Wrote time-shifted mesh to '{outFile}'")


if __name__ == '__main__':
    main()
