#!/usr/bin/env python3
"""Estimate grounding-line ice discharge from a MALI/MPAS-Landice mesh.

The script can use either observed surface velocity as depth-uniform plug flow,
or the modeled 3-D uReconstructX/Y fields.  Modeled interface velocities are
averaged into layers and integrated with layerThicknessFractions.

For every interior edge separating an ice-covered grounded cell from an
ice-covered floating cell, the script evaluates

    q_e = rho_i * H_e * (u_e dot n_grounded_to_floating) * dvEdge_e

using arithmetic cell-to-edge averages for H and u.  Positive flux is from
grounded ice into the floating shelf.

Geometry variables (names can be overridden on the command line):
    thickness, bedTopography, cellsOnEdge, dvEdge, xCell, yCell

Velocity variables:
    observed: observedSurfaceVelocityX, observedSurfaceVelocityY
    modeled:  uReconstructX, uReconstructY, layerThicknessFractions

Optional:
    observedSurfaceVelocityUncertainty

Examples:
    python estimate_gl_flux.py -f landice_grid.nc
    python estimate_gl_flux.py -f output.nc --velocity-source modeled
    python estimate_gl_flux.py -f output.nc --velocity-source both
    python estimate_gl_flux.py -f output.nc --time-index -1 --velocity-units m/yr
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np


SECONDS_PER_YEAR = 365.0 * 24.0 * 60.0 * 60.0
KG_PER_GT = 1.0e12


@dataclass
class FluxResult:
    n_grounded_cells: int
    n_floating_cells: int
    n_gl_edges: int
    signed_gt_per_year: float
    outward_gt_per_year: float
    inward_gt_per_year: float
    absolute_gt_per_year: float
    uncertainty_gt_per_year: Optional[float]


def _as_array(variable: Any, time_index: int, expected_dim: str) -> np.ndarray:
    """Read a 1-D mesh field, selecting Time if present."""
    dims = tuple(variable.dimensions)
    index = []
    for dim in dims:
        if dim.lower() == "time":
            index.append(time_index)
        elif dim == expected_dim:
            index.append(slice(None))
        else:
            # Permit singleton dimensions, but reject ambiguous multi-level data.
            size = variable.shape[len(index)]
            if size != 1:
                raise ValueError(
                    f"Variable {variable.name!r} has unsupported dimension "
                    f"{dim!r} of size {size}; expected only Time and {expected_dim}."
                )
            index.append(0)

    if expected_dim not in dims:
        raise ValueError(
            f"Variable {variable.name!r} does not have expected dimension "
            f"{expected_dim!r}; dimensions are {dims}."
        )
    data = variable[tuple(index)]
    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)
    return np.asarray(data, dtype=float).reshape(-1)


def _as_cell_interface_array(variable: Any, time_index: int) -> np.ndarray:
    """Read a cell-by-vertical velocity field, regardless of dimension order."""
    dims = tuple(variable.dimensions)
    if "nCells" not in dims or "nVertInterfaces" not in dims:
        raise ValueError(
            f"Variable {variable.name!r} must contain nCells and "
            f"nVertInterfaces; dimensions are {dims}."
        )

    index = []
    retained_dims = []
    for axis, dim in enumerate(dims):
        if dim.lower() == "time":
            index.append(time_index)
        elif dim in ("nCells", "nVertInterfaces"):
            index.append(slice(None))
            retained_dims.append(dim)
        elif variable.shape[axis] == 1:
            index.append(0)
        else:
            raise ValueError(
                f"Variable {variable.name!r} has unsupported dimension {dim!r} "
                f"of size {variable.shape[axis]}."
            )

    data = variable[tuple(index)]
    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)
    data = np.asarray(data, dtype=float)
    cell_axis = retained_dims.index("nCells")
    vertical_axis = retained_dims.index("nVertInterfaces")
    data = np.moveaxis(data, (cell_axis, vertical_axis), (0, 1))
    if data.ndim != 2:
        raise ValueError(f"Could not reduce {variable.name!r} to a 2-D field; shape is {data.shape}.")
    return data


def depth_average_velocity(
    velocity: np.ndarray,
    layer_thickness_fractions: np.ndarray,
) -> np.ndarray:
    """Average interface velocities into layers, then return the vertical mean."""
    velocity = np.asarray(velocity, dtype=float)
    fractions = np.asarray(layer_thickness_fractions, dtype=float).reshape(-1)
    if velocity.ndim != 2:
        raise ValueError("Modeled velocity must have shape (nCells, nVertical).")
    if fractions.size == 0 or np.any(~np.isfinite(fractions)) or np.any(fractions < 0.0):
        raise ValueError("layerThicknessFractions must be finite, nonnegative, and nonempty.")
    fraction_sum = float(fractions.sum())
    if fraction_sum <= 0.0:
        raise ValueError("layerThicknessFractions must have a positive sum.")
    if not np.isclose(fraction_sum, 1.0, rtol=1.0e-8, atol=1.0e-12):
        raise ValueError(
            f"layerThicknessFractions sum to {fraction_sum:.16g}, not 1.0."
        )

    if velocity.shape[1] != fractions.size + 1:
        raise ValueError(
            f"Modeled velocity has {velocity.shape[1]} interfaces but "
            f"layerThicknessFractions has {fractions.size} layers."
        )
    layer_velocity = 0.5 * (velocity[:, :-1] + velocity[:, 1:])

    return layer_velocity @ fractions


def _read_cells_on_edge(variable: Any) -> np.ndarray:
    """Read cellsOnEdge and return shape (nEdges, 2), retaining 1-based IDs."""
    dims = tuple(variable.dimensions)
    data = variable[:]
    if np.ma.isMaskedArray(data):
        data = data.filled(0)
    data = np.asarray(data, dtype=np.int64)

    if data.ndim != 2:
        raise ValueError(f"{variable.name!r} must be 2-D; got shape {data.shape}.")
    if data.shape[1] == 2:
        return data
    if data.shape[0] == 2:
        return data.T
    raise ValueError(
        f"Could not identify TWO dimension in {variable.name!r}; shape is "
        f"{data.shape}, dimensions are {dims}."
    )


def _units_text(variable: Any) -> str:
    value = getattr(variable, "units", "")
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    return str(value).strip()


def velocity_factor_to_m_per_year(units_option: str, metadata_units: str) -> tuple[float, str]:
    """Return multiplier converting input velocity to m/yr and a description."""
    if units_option != "auto":
        chosen = units_option
    else:
        normalized = metadata_units.lower().replace(" ", "").replace("**", "^")
        per_second = ("s-1" in normalized or "s^-1" in normalized or "/s" in normalized)
        per_year = any(token in normalized for token in ("yr-1", "yr^-1", "/yr", "year-1", "a-1"))
        if per_second and not per_year:
            chosen = "m/s"
        elif per_year and not per_second:
            chosen = "m/yr"
        else:
            raise ValueError(
                "Could not infer velocity units from NetCDF metadata "
                f"{metadata_units!r}. Pass --velocity-units m/yr or m/s."
            )

    if chosen == "m/s":
        return SECONDS_PER_YEAR, "m/s (converted to m/yr using 365 days/yr)"
    return 1.0, "m/yr"


def compute_flux(
    thickness: np.ndarray,
    bed: np.ndarray,
    cells_on_edge_one_based: np.ndarray,
    edge_length: np.ndarray,
    x_cell: np.ndarray,
    y_cell: np.ndarray,
    velocity_x_m_per_year: np.ndarray,
    velocity_y_m_per_year: np.ndarray,
    velocity_uncertainty_m_per_year: Optional[np.ndarray] = None,
    *,
    ice_density: float = 910.0,
    water_density: float = 1028.0,
    sea_level: float = 0.0,
    min_thickness: float = 1.0,
) -> FluxResult:
    """Compute plug-flow grounding-line discharge.

    Edge-normal directions are calculated from the vector joining the two
    cell centers, which is normal to their shared edge on a planar MPAS Voronoi
    mesh.  The uncertainty calculation treats the supplied value as an
    independent, isotropic 1-sigma uncertainty for each cell's x and y
    velocity components.  Shared cells are accounted for.
    """
    thickness = np.asarray(thickness, dtype=float)
    bed = np.asarray(bed, dtype=float)
    ux = np.asarray(velocity_x_m_per_year, dtype=float)
    uy = np.asarray(velocity_y_m_per_year, dtype=float)
    x_cell = np.asarray(x_cell, dtype=float)
    y_cell = np.asarray(y_cell, dtype=float)
    cells = np.asarray(cells_on_edge_one_based, dtype=np.int64)
    edge_length = np.asarray(edge_length, dtype=float)

    n_cells = thickness.size
    if any(a.shape != (n_cells,) for a in (bed, ux, uy, x_cell, y_cell)):
        raise ValueError(
            "thickness, bed, xCell, yCell, velocity X, and velocity Y must "
            "have the same 1-D shape."
        )
    n_edges = cells.shape[0]
    if cells.shape != (n_edges, 2):
        raise ValueError("cellsOnEdge must have shape (nEdges, 2).")
    if edge_length.shape != (n_edges,):
        raise ValueError("dvEdge must have length nEdges.")
    if ice_density <= 0.0 or water_density <= 0.0:
        raise ValueError("Densities must be positive.")

    finite_geometry = np.isfinite(thickness) & np.isfinite(bed)
    ice = finite_geometry & (thickness > min_thickness)
    water_depth = np.maximum(sea_level - bed, 0.0)
    # Positive flotation residual means ice overburden exceeds ocean pressure.
    flotation_residual = ice_density * thickness - water_density * water_depth
    grounded = ice & (flotation_residual > 0.0)
    floating = ice & ~grounded

    # MPAS connectivity is 1-based; zero denotes no neighboring cell.
    interior = (cells[:, 0] > 0) & (cells[:, 1] > 0)
    c0_raw = cells[:, 0] - 1
    c1_raw = cells[:, 1] - 1
    valid_ids = (
        (c0_raw >= 0) & (c0_raw < n_cells)
        & (c1_raw >= 0) & (c1_raw < n_cells)
    )
    # Clipping makes indexing safe before the interior/valid-ID mask is applied.
    c0 = np.clip(c0_raw, 0, max(n_cells - 1, 0))
    c1 = np.clip(c1_raw, 0, max(n_cells - 1, 0))

    gf = interior & valid_ids & grounded[c0] & floating[c1]
    fg = interior & valid_ids & floating[c0] & grounded[c1]
    gl = gf | fg
    edge_ids = np.flatnonzero(gl)

    if edge_ids.size == 0:
        return FluxResult(
            int(grounded.sum()), int(floating.sum()), 0,
            0.0, 0.0, 0.0, 0.0,
            0.0 if velocity_uncertainty_m_per_year is not None else None,
        )

    i = c0[edge_ids]
    j = c1[edge_ids]
    # The vector from c0 to c1 is normal to their shared Voronoi edge. Reverse
    # it where c1 is grounded so every normal points grounded -> floating.
    orientation = np.where(gf[edge_ids], 1.0, -1.0)
    dx = x_cell[j] - x_cell[i]
    dy = y_cell[j] - y_cell[i]
    center_distance = np.hypot(dx, dy)
    if np.any(~np.isfinite(center_distance)) or np.any(center_distance <= 0.0):
        bad = int((~np.isfinite(center_distance) | (center_distance <= 0.0)).sum())
        raise ValueError(f"Found {bad} grounding-line edges with invalid cell-center geometry.")
    nx = orientation * dx / center_distance
    ny = orientation * dy / center_distance

    h_edge = 0.5 * (thickness[i] + thickness[j])
    ux_edge = 0.5 * (ux[i] + ux[j])
    uy_edge = 0.5 * (uy[i] + uy[j])
    normal_velocity = ux_edge * nx + uy_edge * ny
    mass_flux = ice_density * h_edge * normal_velocity * edge_length[edge_ids]  # kg/yr

    finite_edge = (
        np.isfinite(mass_flux)
        & np.isfinite(edge_length[edge_ids])
        & (edge_length[edge_ids] > 0.0)
    )
    if not np.all(finite_edge):
        bad = int((~finite_edge).sum())
        raise ValueError(f"Found {bad} grounding-line edges with invalid flux or nonpositive dvEdge.")

    signed = float(mass_flux.sum())
    outward = float(mass_flux[mass_flux > 0.0].sum())
    inward = float(mass_flux[mass_flux < 0.0].sum())
    absolute = float(np.abs(mass_flux).sum())

    uncertainty_gt_per_year = None
    if velocity_uncertainty_m_per_year is not None:
        sigma = np.asarray(velocity_uncertainty_m_per_year, dtype=float)
        if sigma.shape != (n_cells,):
            raise ValueError("Velocity uncertainty must have the same shape as thickness.")
        gl_cells = np.unique(np.concatenate((i, j)))
        if np.any(~np.isfinite(sigma[gl_cells])) or np.any(sigma[gl_cells] < 0.0):
            raise ValueError("Velocity uncertainty is missing, nonfinite, or negative on grounding-line cells.")

        # Q is linear in cell velocity. Aggregate coefficients first so that
        # correlations caused by one cell participating in several edges are
        # represented exactly under the independent-cell assumption.
        k = ice_density * h_edge * edge_length[edge_ids]
        coeff_x = np.zeros(n_cells, dtype=float)
        coeff_y = np.zeros(n_cells, dtype=float)
        np.add.at(coeff_x, i, 0.5 * k * nx)
        np.add.at(coeff_x, j, 0.5 * k * nx)
        np.add.at(coeff_y, i, 0.5 * k * ny)
        np.add.at(coeff_y, j, 0.5 * k * ny)
        variance = np.sum((sigma * coeff_x) ** 2 + (sigma * coeff_y) ** 2)
        uncertainty_gt_per_year = math.sqrt(float(variance)) / KG_PER_GT

    return FluxResult(
        n_grounded_cells=int(grounded.sum()),
        n_floating_cells=int(floating.sum()),
        n_gl_edges=int(edge_ids.size),
        signed_gt_per_year=signed / KG_PER_GT,
        outward_gt_per_year=outward / KG_PER_GT,
        inward_gt_per_year=inward / KG_PER_GT,
        absolute_gt_per_year=absolute / KG_PER_GT,
        uncertainty_gt_per_year=uncertainty_gt_per_year,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Estimate observed and/or modeled grounding-line discharge on a MALI mesh.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-f", "--file", required=True, help="MALI NetCDF input file")
    parser.add_argument("--time-index", type=int, default=-1, help="Time record (Python indexing)")
    parser.add_argument(
        "--velocity-source", choices=("observed", "modeled", "both"), default="observed",
        help="Velocity field(s) used for the flux calculation",
    )
    parser.add_argument("--ice-density", type=float, default=910.0, help="Ice density (kg m-3)")
    parser.add_argument("--water-density", type=float, default=1028.0, help="Ocean-water density (kg m-3)")
    parser.add_argument("--sea-level", type=float, default=0.0, help="Sea-surface elevation (m)")
    parser.add_argument(
        "--min-thickness", type=float, default=1.0,
        help="Cells at or below this thickness are treated as ice-free (m)",
    )
    parser.add_argument(
        "--velocity-units", choices=("auto", "m/yr", "m/s"), default="auto",
        help="Velocity units; auto reads each NetCDF units attribute",
    )

    names = parser.add_argument_group("NetCDF variable names")
    names.add_argument("--thickness-var", default="thickness")
    names.add_argument("--bed-var", default="bedTopography")
    names.add_argument("--cells-on-edge-var", default="cellsOnEdge")
    names.add_argument("--edge-length-var", default="dvEdge")
    names.add_argument("--x-cell-var", default="xCell")
    names.add_argument("--y-cell-var", default="yCell")
    names.add_argument("--velocity-x-var", default="observedSurfaceVelocityX")
    names.add_argument("--velocity-y-var", default="observedSurfaceVelocityY")
    names.add_argument("--uncertainty-var", default="observedSurfaceVelocityUncertainty")
    names.add_argument("--modeled-velocity-x-var", default="uReconstructX")
    names.add_argument("--modeled-velocity-y-var", default="uReconstructY")
    names.add_argument("--layer-fractions-var", default="layerThicknessFractions")
    names.add_argument(
        "--no-uncertainty", action="store_true",
        help="Do not read or propagate the observed-velocity uncertainty",
    )
    return parser


def _require_variable(dataset: Any, name: str) -> Any:
    if name not in dataset.variables:
        available = ", ".join(sorted(dataset.variables))
        raise KeyError(f"Required variable {name!r} not found. Available variables: {available}")
    return dataset.variables[name]


def _convert_velocity_pair_to_m_per_year(
    ux: np.ndarray,
    uy: np.ndarray,
    ux_variable: Any,
    uy_variable: Any,
    units_option: str,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Convert a pair of velocity components using their metadata or an override."""
    x_units = _units_text(ux_variable)
    y_units = _units_text(uy_variable)
    x_factor, description = velocity_factor_to_m_per_year(units_option, x_units)
    y_factor, _ = velocity_factor_to_m_per_year(units_option, y_units)
    if y_factor != x_factor:
        raise ValueError(
            f"Velocity component units are incompatible: X={x_units!r}, Y={y_units!r}."
        )
    return ux * x_factor, uy * y_factor, description


def _print_result(label: str, result: FluxResult, units_description: str) -> None:
    print(f"\n{label} velocity")
    print(f"  Input velocity units:         {units_description}")
    print(f"  Grounded ice cells:           {result.n_grounded_cells}")
    print(f"  Floating ice cells:           {result.n_floating_cells}")
    print(f"  Grounding-line edges:         {result.n_gl_edges}")
    print(f"  Signed GL discharge:          {result.signed_gt_per_year:.6g} Gt/yr")
    print(f"  Outward-only contribution:    {result.outward_gt_per_year:.6g} Gt/yr")
    print(f"  Inward-only contribution:     {result.inward_gt_per_year:.6g} Gt/yr")
    print(f"  Sum of absolute contributions:{result.absolute_gt_per_year: .6g} Gt/yr")
    if result.uncertainty_gt_per_year is not None:
        print(f"  Approx. 1-sigma uncertainty:  {result.uncertainty_gt_per_year:.6g} Gt/yr")
        print("    (Assumes independent cells and isotropic component uncertainties.)")


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        from netCDF4 import Dataset
    except ImportError:
        print(
            "ERROR: This script requires netCDF4 (for example, `conda install netcdf4` "
            "or `python -m pip install netCDF4`).",
            file=sys.stderr,
        )
        return 2

    try:
        with Dataset(args.file, "r") as ds:
            thickness = _as_array(_require_variable(ds, args.thickness_var), args.time_index, "nCells")
            bed = _as_array(_require_variable(ds, args.bed_var), args.time_index, "nCells")
            cells_on_edge = _read_cells_on_edge(_require_variable(ds, args.cells_on_edge_var))
            edge_length = _as_array(
                _require_variable(ds, args.edge_length_var), args.time_index, "nEdges"
            )
            x_cell = _as_array(_require_variable(ds, args.x_cell_var), args.time_index, "nCells")
            y_cell = _as_array(_require_variable(ds, args.y_cell_var), args.time_index, "nCells")

            results: dict[str, FluxResult] = {}
            units_descriptions: dict[str, str] = {}

            if args.velocity_source in ("observed", "both"):
                observed_ux_var = _require_variable(ds, args.velocity_x_var)
                observed_uy_var = _require_variable(ds, args.velocity_y_var)
                observed_ux = _as_array(observed_ux_var, args.time_index, "nCells")
                observed_uy = _as_array(observed_uy_var, args.time_index, "nCells")
                observed_ux, observed_uy, units_descriptions["observed"] = (
                    _convert_velocity_pair_to_m_per_year(
                        observed_ux, observed_uy, observed_ux_var, observed_uy_var,
                        args.velocity_units,
                    )
                )

                uncertainty = None
                if not args.no_uncertainty:
                    if args.uncertainty_var in ds.variables:
                        uncertainty_var = ds.variables[args.uncertainty_var]
                        uncertainty = _as_array(uncertainty_var, args.time_index, "nCells")
                        uncertainty_factor, _ = velocity_factor_to_m_per_year(
                            args.velocity_units, _units_text(uncertainty_var)
                        )
                        uncertainty *= uncertainty_factor
                    else:
                        print(
                            f"WARNING: {args.uncertainty_var!r} is absent; "
                            "continuing without observed-velocity uncertainty.",
                            file=sys.stderr,
                        )

                results["observed"] = compute_flux(
                    thickness=thickness,
                    bed=bed,
                    cells_on_edge_one_based=cells_on_edge,
                    edge_length=edge_length,
                    x_cell=x_cell,
                    y_cell=y_cell,
                    velocity_x_m_per_year=observed_ux,
                    velocity_y_m_per_year=observed_uy,
                    velocity_uncertainty_m_per_year=uncertainty,
                    ice_density=args.ice_density,
                    water_density=args.water_density,
                    sea_level=args.sea_level,
                    min_thickness=args.min_thickness,
                )

            if args.velocity_source in ("modeled", "both"):
                modeled_ux_var = _require_variable(ds, args.modeled_velocity_x_var)
                modeled_uy_var = _require_variable(ds, args.modeled_velocity_y_var)
                modeled_ux_interfaces = _as_cell_interface_array(
                    modeled_ux_var, args.time_index
                )
                modeled_uy_interfaces = _as_cell_interface_array(
                    modeled_uy_var, args.time_index
                )
                if modeled_ux_interfaces.shape != modeled_uy_interfaces.shape:
                    raise ValueError(
                        "Modeled X and Y velocity fields have different shapes: "
                        f"{modeled_ux_interfaces.shape} and {modeled_uy_interfaces.shape}."
                    )
                layer_fractions = _as_array(
                    _require_variable(ds, args.layer_fractions_var),
                    args.time_index,
                    "nVertLevels",
                )
                modeled_ux = depth_average_velocity(modeled_ux_interfaces, layer_fractions)
                modeled_uy = depth_average_velocity(modeled_uy_interfaces, layer_fractions)
                modeled_ux, modeled_uy, units_descriptions["modeled"] = (
                    _convert_velocity_pair_to_m_per_year(
                        modeled_ux, modeled_uy, modeled_ux_var, modeled_uy_var,
                        args.velocity_units,
                    )
                )
                results["modeled"] = compute_flux(
                    thickness=thickness,
                    bed=bed,
                    cells_on_edge_one_based=cells_on_edge,
                    edge_length=edge_length,
                    x_cell=x_cell,
                    y_cell=y_cell,
                    velocity_x_m_per_year=modeled_ux,
                    velocity_y_m_per_year=modeled_uy,
                    ice_density=args.ice_density,
                    water_density=args.water_density,
                    sea_level=args.sea_level,
                    min_thickness=args.min_thickness,
                )
    except (KeyError, ValueError, IndexError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print(f"Input file: {args.file}")
    print(f"Time index: {args.time_index}")
    if "observed" in results:
        _print_result("Observed surface (plug-flow)", results["observed"], units_descriptions["observed"])
    if "modeled" in results:
        _print_result("Modeled depth-averaged", results["modeled"], units_descriptions["modeled"])
    if "observed" in results and "modeled" in results:
        observed_flux = results["observed"].signed_gt_per_year
        modeled_flux = results["modeled"].signed_gt_per_year
        difference = modeled_flux - observed_flux
        print("\nModeled minus observed")
        print(f"  Signed-flux difference:      {difference:.6g} Gt/yr")
        if observed_flux != 0.0:
            print(f"  Relative difference:         {100.0 * difference / observed_flux:.6g}%")
    print("Positive signed discharge is from grounded ice toward floating ice.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
