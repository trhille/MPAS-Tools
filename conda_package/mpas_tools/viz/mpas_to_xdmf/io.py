import glob
import importlib.resources
import itertools
import os

import h5py
import numpy as np
import xarray as xr
from jinja2 import Template
from tqdm import tqdm

from mpas_tools.viz.mpas_to_xdmf.geometry import (
    _build_cell_geometry,
    _build_edge_geometry,
    _build_vertex_geometry,
)
from mpas_tools.viz.mpas_to_xdmf.mesh import _get_ds_mesh
from mpas_tools.viz.mpas_to_xdmf.time import _set_time

# Fields with extra dimensions (e.g. nVertLevels) are unwrapped into one
# 2D field per index.  Reading each of those fields individually is very
# slow, because a single vertical level is strided across the whole variable
# on disk, so the entire variable has to be read once per level.  Instead we
# read as many indices at a time as we can, capped by this many bytes so that
# memory use stays bounded on very large meshes.
_DEFAULT_MAX_READ_BYTES = 2 * 1024**3


def _load_dataset(mesh_filename, time_series_filenames, variables, xtime_var):
    """
    Load the MPAS mesh file and optionally combine it with time series
    files into a mesh dataset and a dataset to convert.

    Parameters
    ----------
    mesh_filename : str
        Path to the MPAS mesh file.
    time_series_filenames : list of str or str
        List of filenames or a wildcard string for time series files.
    variables : list of str
        List of variables to convert. Special keys include:
        - 'allOnCells': Load all variables with dimension `nCells`.
        - 'allOnEdges': Load all variables with dimension `nEdges`.
        - 'allOnVertices': Load all variables with dimension `nVertices`.
    xtime_var : str
        Name of the variable containing time information.

    Returns
    -------
    ds_mesh : xarray.Dataset
        The mesh dataset.
    ds : xarray.Dataset
        The dataset to convert, which may include the mesh variables.
    """
    # Load the mesh file
    ds_mesh = xr.open_dataset(mesh_filename)

    if time_series_filenames is None:
        ds = ds_mesh
    else:
        if isinstance(time_series_filenames, str):
            # Deterministic ordering for reproducibility
            file_list = sorted(glob.glob(time_series_filenames))
        else:
            file_list = time_series_filenames

        if len(file_list) == 0:
            raise ValueError(
                f'No time series files found matching '
                f"'{time_series_filenames}'"
            )

        if len(file_list) == 1:
            # If only one file, open it directly
            ds = xr.open_dataset(file_list[0])
        else:
            ds = xr.open_mfdataset(
                file_list,
                combine='nested',
                concat_dim='Time',
                data_vars='minimal',
                coords='minimal',
                compat='override',
                decode_times=False,
                decode_timedelta=False,
            )

    ds_mesh = _get_ds_mesh(ds_mesh)
    if variables is not None:
        selected_vars = set()
        for var in variables:
            if var == 'allOnCells':
                selected_vars.update(
                    [v for v in ds.data_vars if 'nCells' in ds[v].dims]
                )
            elif var == 'allOnEdges':
                selected_vars.update(
                    [v for v in ds.data_vars if 'nEdges' in ds[v].dims]
                )
            elif var == 'allOnVertices':
                selected_vars.update(
                    [v for v in ds.data_vars if 'nVertices' in ds[v].dims]
                )
            else:
                selected_vars.add(var)
        if xtime_var is not None:
            selected_vars.add(xtime_var)
        ds = ds[list(selected_vars)]

    _set_time(ds=ds, xtime_var=xtime_var)

    return ds_mesh, ds


def _convert_to_xdmf(
    ds,
    ds_mesh,
    out_dir,
    extra_dims=None,
    quiet=False,
    float32=False,
    max_read_bytes=_DEFAULT_MAX_READ_BYTES,
):
    """
    Convert an xarray Dataset to XDMF + HDF5 format.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset to convert.
    ds_mesh : xarray.Dataset
        The mesh dataset.
    out_dir : str
        Directory where XDMF and HDF5 files will be saved.
    extra_dims : dict, optional
        Dictionary mapping extra dimensions to the indices to write.
    quiet : bool, optional
        If True, suppress progress output. Default is False.
    float32 : bool, optional
        If True, write floating-point fields in single precision.
    max_read_bytes : int, optional
        Approximate limit on the number of bytes read from ``ds`` at a time.
    """
    os.makedirs(out_dir, exist_ok=True)

    kwargs = dict(
        out_dir=out_dir,
        extra_dims=extra_dims,
        quiet=quiet,
        float32=float32,
        max_read_bytes=max_read_bytes,
    )

    if 'nCells' in ds.dims:
        _convert_cells_to_xdmf(ds, ds_mesh, **kwargs)
    if 'nEdges' in ds.dims:
        _convert_edges_to_xdmf(ds, ds_mesh, **kwargs)
    if 'nVertices' in ds.dims:
        _convert_vertices_to_xdmf(ds, ds_mesh, **kwargs)


def _convert_cells_to_xdmf(ds, ds_mesh, **kwargs):
    """
    Convert cell-centered data to XDMF + HDF5 format.
    """
    ds_cell_geom = _build_cell_geometry(ds_mesh)
    cell_vars = [var for var in ds.data_vars if 'nCells' in ds[var].dims]
    ds_cells = ds[cell_vars]
    _write_xdmf(
        ds_cell_geom, ds_cells, suffix='Cells', topo_dim='nCells', **kwargs
    )


def _convert_edges_to_xdmf(ds, ds_mesh, **kwargs):
    """
    Convert edge-centered data to XDMF + HDF5 format.
    """
    ds_edge_geom = _build_edge_geometry(ds_mesh)
    edge_vars = [var for var in ds.data_vars if 'nEdges' in ds[var].dims]
    ds_edges = ds[edge_vars]
    _write_xdmf(
        ds_edge_geom, ds_edges, suffix='Edges', topo_dim='nEdges', **kwargs
    )


def _convert_vertices_to_xdmf(ds, ds_mesh, **kwargs):
    """
    Convert vertex-centered data to XDMF + HDF5 format.
    """
    ds_vertex_geom = _build_vertex_geometry(ds_mesh)
    vertex_vars = [var for var in ds.data_vars if 'nVertices' in ds[var].dims]
    ds_vertices = ds[vertex_vars]
    # each vertex is repeated once per kite; the map is applied to each field
    # after it has been read rather than as a (much slower) indexed read
    vert_to_kite_map = ds_vertex_geom['vert_to_kite_map'].values
    _write_xdmf(
        ds_vertex_geom,
        ds_vertices,
        suffix='Vertices',
        topo_dim='nVertices',
        topo_map=vert_to_kite_map,
        **kwargs,
    )


def _write_xdmf(
    ds_geom,
    ds_data,
    out_dir,
    suffix,
    topo_dim,
    extra_dims=None,
    topo_map=None,
    quiet=False,
    float32=False,
    max_read_bytes=_DEFAULT_MAX_READ_BYTES,
):
    """
    Write data to HDF5 and metadata to XDMF format.

    Parameters
    ----------
    ds_geom : xarray.Dataset
        Dataset containing geometry information (e.g., points, connectivity).
    ds_data : xarray.Dataset
        Dataset containing time-varying data to write.
    out_dir : str
        Directory where XDMF and HDF5 files will be saved.
    suffix : str
        Suffix to append to output filenames (e.g., 'Cells', 'Edges').
    topo_dim : str
        The dimension the fields are defined on (e.g. ``'nCells'``).
    extra_dims : dict, optional
        Dictionary mapping extra dimensions to the indices to write.  Each
        variable is unwrapped into one field per combination of indices.
    topo_map : numpy.ndarray, optional
        Indices along ``topo_dim`` to map each field onto the geometry.
    quiet : bool, optional
        If True, suppress progress output. Default is False.
    float32 : bool, optional
        If True, write floating-point fields in single precision.
    max_read_bytes : int, optional
        Approximate limit on the number of bytes read from ``ds_data`` at a
        time.
    """
    h5_basename = f'fieldsOn{suffix}.h5'
    h5_filename = os.path.join(out_dir, h5_basename)
    xdmf_filename = os.path.join(out_dir, f'fieldsOn{suffix}.xdmf')

    variables_metadata = _field_metadata(ds_data, extra_dims)

    # Write HDF5 file
    with h5py.File(h5_filename, 'w') as h5_file:
        # Write geometry
        h5_file.create_dataset('Points', data=ds_geom['points'].values)
        h5_file.create_dataset('Cells', data=_as_index_type(ds_geom['cells']))

        # Calculate total progress steps
        total_steps = sum(
            ds_data.sizes['Time'] if field['has_time'] else 1
            for field in variables_metadata
        )

        # Write time-varying and static data with progress bar
        if quiet:
            iterator = None
        else:
            iterator = tqdm(total=total_steps, desc=f'Writing {suffix}')
        for var_name in ds_data.data_vars:
            if iterator is not None:
                iterator.set_description(f'Processing {var_name}')
            _write_variable(
                h5_file=h5_file,
                da=ds_data[var_name],
                var_name=var_name,
                extra_dims=extra_dims,
                topo_dim=topo_dim,
                topo_map=topo_map,
                float32=float32,
                max_read_bytes=max_read_bytes,
                iterator=iterator,
            )
        if iterator is not None:
            iterator.close()

    # Load XDMF template from external file
    package = 'mpas_tools.viz.mpas_to_xdmf.templates'
    filename = 'xdmf_template.xml'
    with importlib.resources.open_text(package, filename) as template_file:
        xdmf_template = Template(template_file.read())

    # Render XDMF file
    cells = ds_geom['cells']
    times = ds_data['Time'].values if 'Time' in ds_data.dims else []

    xdmf_content = xdmf_template.render(
        times=times,
        num_elements=cells.shape[0],
        num_verts=cells.shape[1],
        num_points=ds_geom['points'].shape[0],
        variables=variables_metadata,
        suffix=suffix,
        h5_basename=h5_basename,
    )

    with open(xdmf_filename, 'w') as xdmf_file:
        xdmf_file.write(xdmf_content)


def _as_index_type(da):
    """
    Return the connectivity array as 32-bit integers if the indices fit,
    halving the size of the largest static array in the HDF5 file.
    """
    values = da.values
    if (
        values.dtype.itemsize > 4
        and values.size > 0
        and values.max() < np.iinfo(np.int32).max
    ):
        values = values.astype(np.int32)
    return values


def _expand_extra_dims(var_name, da, extra_dims):
    """
    Determine the fields a variable is unwrapped into, one per combination of
    extra-dimension indices.

    Parameters
    ----------
    var_name : str
        The name of the variable.
    da : xarray.DataArray
        The variable to unwrap.
    extra_dims : dict or None
        Dictionary mapping extra dimensions to the indices to write.

    Returns
    -------
    dims : list of str
        The extra dimensions present on ``da``, in the order given by
        ``extra_dims``.
    fields : list of tuple
        A ``(name, indices)`` pair for each field, where ``indices`` gives the
        index along each dimension in ``dims``.
    """
    dims = [dim for dim in (extra_dims or {}) if dim in da.dims]
    index_lists = [extra_dims[dim] for dim in dims]
    fields = [
        (var_name + ''.join(f'_{index}' for index in indices), indices)
        for indices in itertools.product(*index_lists)
    ]
    return dims, fields


def _field_metadata(ds_data, extra_dims):
    """
    Build the list of fields to write, in the order they are written, for use
    in the XDMF template.
    """
    metadata = []
    for var_name in ds_data.data_vars:
        da = ds_data[var_name]
        _, fields = _expand_extra_dims(var_name, da, extra_dims)
        has_time = 'Time' in da.dims
        metadata.extend(
            {'name': name, 'has_time': has_time} for name, _ in fields
        )
    return metadata


def _read_group_size(da, dims, extra_dims, topo_dim, max_read_bytes):
    """
    Determine how many indices along the first extra dimension can be read at
    once without exceeding ``max_read_bytes``.  At least one index is always
    read.
    """
    per_index = da.sizes[topo_dim] * da.dtype.itemsize
    for dim in dims[1:]:
        indices = extra_dims[dim]
        per_index *= max(indices) - min(indices) + 1
    return max(1, int(max_read_bytes // max(per_index, 1)))


def _write_variable(
    h5_file,
    da,
    var_name,
    extra_dims,
    topo_dim,
    topo_map,
    float32,
    max_read_bytes,
    iterator,
):
    """
    Unwrap a variable into one HDF5 dataset per combination of extra-dimension
    indices (and per time index), reading as many indices at a time as the
    memory budget allows.
    """
    dims, _ = _expand_extra_dims(var_name, da, extra_dims)
    has_time = 'Time' in da.dims
    time_indices = range(da.sizes['Time']) if has_time else [None]

    if dims:
        group_size = _read_group_size(
            da, dims, extra_dims, topo_dim, max_read_bytes
        )
        first_indices = extra_dims[dims[0]]
        groups = [
            first_indices[start : start + group_size]
            for start in range(0, len(first_indices), group_size)
        ]
    else:
        groups = [None]

    for t_index in time_indices:
        da_t = da.isel(Time=t_index) if has_time else da
        time_suffix = f'_t{t_index}' if has_time else ''
        for group in groups:
            block, offsets, index_lists = _read_block(
                da_t, dims, extra_dims, group
            )
            _write_block(
                h5_file=h5_file,
                block=block,
                dims=dims,
                offsets=offsets,
                index_lists=index_lists,
                var_name=var_name,
                time_suffix=time_suffix,
                topo_map=topo_map,
                float32=float32,
                iterator=iterator,
            )
            # release the block before the next one is read, so that only one
            # is ever in memory at a time
            del block


def _read_block(da, dims, extra_dims, group):
    """
    Read the smallest contiguous span of ``da`` that covers ``group`` (a set of
    indices along ``dims[0]``) together with all the requested indices of the
    remaining extra dimensions, in a single pass over the variable on disk.

    Returns the block along with the offset of each extra dimension within it
    and the indices that the block was read for.
    """
    if group is None:
        return da.compute(), {}, []

    offsets = {dims[0]: min(group)}
    selection = {dims[0]: slice(min(group), max(group) + 1)}
    for dim in dims[1:]:
        indices = extra_dims[dim]
        offsets[dim] = min(indices)
        selection[dim] = slice(min(indices), max(indices) + 1)
    index_lists = [group] + [extra_dims[dim] for dim in dims[1:]]
    return da.isel(selection).compute(), offsets, index_lists


def _write_block(
    h5_file,
    block,
    dims,
    offsets,
    index_lists,
    var_name,
    time_suffix,
    topo_map,
    float32,
    iterator,
):
    """
    Write one HDF5 dataset per combination of extra-dimension indices in an
    in-memory block.
    """
    for indices in itertools.product(*index_lists):
        name = (
            var_name + ''.join(f'_{index}' for index in indices) + time_suffix
        )
        field = block.isel(
            {
                dim: index - offsets[dim]
                for dim, index in zip(dims, indices, strict=True)
            }
        )
        values = field.values
        if topo_map is not None:
            values = values[topo_map]
        if float32 and values.dtype.kind == 'f':
            values = values.astype(np.float32)
        h5_file.create_dataset(name, data=values)
        if iterator is not None:
            iterator.update(1)


def _parse_extra_dims(dimension_list, ds):
    """
    Parse and prompt for indices of extra dimensions.

    Parameters
    ----------
    dimension_list : list of str
        List of dimensions and indices in the format <dimension>=<indices>.
    ds : xarray.Dataset
        Dataset to extract dimensions from.

    Returns
    -------
    extra_dims : dict
        Dictionary mapping dimensions to their selected indices.
    """
    extra_dims = {}
    unspecified_dims = []

    for dim in ds.dims:
        if dim in ['Time', 'nCells', 'nEdges', 'nVertices']:
            continue

        # Check if the dimension is specified in the dimension_list
        specified = False
        if dimension_list:
            for dim_spec in dimension_list:
                if dim_spec.startswith(f'{dim}='):
                    indices = dim_spec.split('=')[1]
                    extra_dims[dim] = _parse_indices(indices, ds.sizes[dim])
                    specified = True
                    break

        if not specified:
            unspecified_dims.append(dim)

    # If there are unspecified dimensions, display a detailed prompt
    if unspecified_dims:
        print('\nThe following dimensions require indices to be specified:')
        print('You can enter indices in one of the following formats:')
        print("  - A single index (e.g., '0')")
        print("  - A comma-separated list of indices (e.g., '0,1,2')")
        print("  - A range with optional stride (e.g., ':' or '0:10:2')")
        print('  - Leave blank to skip fields with this dimension.')
        print()

        for dim in unspecified_dims:
            print(f"Dimension '{dim}' has size {ds.sizes[dim]}.")
            indices = input(f"Enter indices to keep for '{dim}': ")
            extra_dims[dim] = _parse_indices(indices, ds.sizes[dim])

    return extra_dims


def _parse_indices(index_string, dim_size):
    """
    Parse an index string into a list of indices.

    Parameters
    ----------
    index_string : str
        Index string (e.g., "0,1,2", "0:10:2", ":").
    dim_size : int
        Size of the dimension.

    Returns
    -------
    indices : list of int
        Parsed indices.
    """
    if not index_string:
        return []
    if ':' in index_string:
        # Support slice notation like ':', '0:10', '0:10:2', etc.
        parts = index_string.split(':')
        # Validate that parts has at most 3 elements
        if len(parts) > 3:
            raise ValueError(
                f"Invalid index string '{index_string}': too many colons. "
                'Expected at most two colons.'
            )
        # Pad parts to length 3 with empty strings if needed
        while len(parts) < 3:
            parts.append('')
        # Convert to int or None
        start = int(parts[0]) if parts[0] else 0
        stop = int(parts[1]) if parts[1] else dim_size
        step = int(parts[2]) if parts[2] else 1
        return list(range(start, stop, step))
    return [int(i) for i in index_string.split(',')]


def _process_extra_dims(ds, extra_dims):
    """
    Process extra dimensions in the dataset by ensuring all are covered and
    dropping variables with extra dimensions for which no indices were
    selected.

    Variables that are kept retain their extra dimensions here.  They are
    unwrapped into one field per combination of indices (with the suffix
    ``_<index>`` for each extra dimension) as they are written, so that each
    variable only needs to be read from disk once.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to process for extra dimensions.
    extra_dims : dict
        Dictionary mapping dimensions to their selected indices.

    Returns
    -------
    ds : xarray.Dataset
        The processed dataset with extra dimensions handled.

    Raises
    ------
    ValueError
        If any extra dimensions are not covered by the extra_dims dictionary.
    """
    basic_dims = ['Time', 'nCells', 'nEdges', 'nVertices']
    for dim in ds.dims:
        if dim not in basic_dims and dim not in (extra_dims or {}):
            raise ValueError(
                f"Dimension '{dim}' is not covered by the extra_dims "
                f'dictionary.'
            )

    if extra_dims:
        # Drop variables with a dimension for which no indices were selected
        vars_to_drop = {
            var: None
            for dim, indices in extra_dims.items()
            if not indices
            for var in ds.data_vars
            if dim in ds[var].dims
        }
        if vars_to_drop:
            ds = ds.drop_vars(list(vars_to_drop))

    return ds
