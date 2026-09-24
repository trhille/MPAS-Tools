#!/usr/bin/env python3
"""Map regional fractional basal-melt uncertainties to a MALI mesh."""

import argparse
import os

import numpy as np

RHO_ICE = 910.0
RHO_OCEAN = 1028.0


def cell_field(dataset, name, time_index):
    """Read an nCells field, selecting time_index if it has a Time axis."""
    variable = dataset.variables[name]
    index = [slice(None)] * variable.ndim
    if "Time" in variable.dimensions:
        index[variable.dimensions.index("Time")] = time_index
    return np.ma.filled(variable[tuple(index)], np.nan).astype(float).reshape(-1)


def decode_region_names(variable):
    characters = np.ma.filled(variable[:], b" ")
    return [b"".join(row).decode().rstrip(" \x00") for row in characters]


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Convert regional basal-melt flux uncertainties to cellwise "
            "floatingBasalMassBal uncertainties."
        )
    )
    parser.add_argument("-f", "--file", required=True, help="MALI mesh/output file")
    parser.add_argument("-r", "--regions", required=True, help="Regional masks file")
    parser.add_argument("-c", "--csv", required=True, help="Regional melt CSV")
    parser.add_argument("-o", "--output", required=True, help="Output NetCDF file")
    parser.add_argument("--time-index", type=int, default=-1)
    args = parser.parse_args(argv)

    if os.path.exists(args.output):
        parser.error(f"output file already exists: {args.output}")

    try:
        from netCDF4 import Dataset
    except ImportError:
        parser.error("netCDF4 is required (for example: conda install netcdf4)")

    with Dataset(args.file) as mesh:
        thickness = cell_field(mesh, "thickness", args.time_index)
        bed = cell_field(mesh, "bedTopography", args.time_index)
        basal_melt = cell_field(mesh, "floatingBasalMassBal", args.time_index)
        units = getattr(mesh.variables["floatingBasalMassBal"], "units", "m s-1")

    with Dataset(args.regions) as regions:
        masks = np.asarray(regions.variables["regionCellMasks"][:], dtype=bool)
        names = decode_region_names(regions.variables["regionNames"])

    table = np.loadtxt(args.csv, delimiter=",", skiprows=1)
    region_indices = table[:, 0].astype(int)
    regional_melt = table[:, 1]
    regional_uncertainty = table[:, -1]

    n_cells = thickness.size
    n_regions = masks.shape[1]
    if masks.shape[0] != n_cells:
        parser.error("regionCellMasks and the MALI file have different nCells")
    if len(names) != n_regions or not np.array_equal(region_indices, np.arange(n_regions)):
        parser.error("CSV rows must be ordered region indices 0..nRegions-1")
    if np.any(~np.isfinite(regional_melt)) or np.any(regional_melt == 0.0):
        parser.error("regional melt values must be finite and nonzero")
    if np.any(~np.isfinite(regional_uncertainty)):
        parser.error("regional uncertainties must be finite")

    fractional_uncertainty = regional_uncertainty / np.abs(regional_melt)
    water_depth = np.maximum(-bed, 0.0)
    floating = (thickness > 1.0) & (RHO_ICE * thickness <= RHO_OCEAN * water_depth)
    memberships = masks.sum(axis=1)
    if np.any(floating & (memberships != 1)):
        parser.error("each floating cell must belong to exactly one region")
    if np.any(floating & ~np.isfinite(basal_melt)):
        parser.error("floatingBasalMassBal contains non-finite values on floating ice")

    cell_uncertainty = np.zeros(n_cells)
    for region in range(n_regions):
        cells = floating & masks[:, region]
        cell_uncertainty[cells] = (
            fractional_uncertainty[region] * np.abs(basal_melt[cells])
        )

    with Dataset(args.output, "w", format="NETCDF3_64BIT_OFFSET") as output:
        output.createDimension("nCells", n_cells)
        variable = output.createVariable(
            "floatingBasalMassBalUncertainty", "f8", ("nCells",)
        )
        variable[:] = cell_uncertainty
        variable.units = units
        variable.long_name = "cellwise uncertainty in floating basal mass balance"
        output.source_csv = os.path.abspath(args.csv)

    print(f"Wrote {args.output}")
    print(f"Mapped {floating.sum()} floating cells")
    print("Region  Fractional uncertainty")
    for index, name, fraction in zip(region_indices, names, fractional_uncertainty):
        print(f"{index:6d}  {fraction:22.6g}  {name}")


if __name__ == "__main__":
    main()
