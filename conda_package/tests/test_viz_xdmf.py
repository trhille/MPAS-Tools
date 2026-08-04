import os
import sys

import h5py
import numpy as np
import pytest
import xarray as xr

from mpas_tools.io import write_netcdf
from mpas_tools.viz.mpas_to_xdmf.io import (
    _load_dataset,
    _parse_indices,
    _process_extra_dims,
)
from mpas_tools.viz.mpas_to_xdmf.mpas_to_xdmf import MpasToXdmf, main
from mpas_tools.viz.mpas_to_xdmf.time import _set_time

from .util import get_test_data_file

TEST_MESH = get_test_data_file('mesh.QU.1920km.151026.nc')

NTIME = 3
NVERT_LEVELS = 5
NTRACERS = 2


def _make_test_dataset(ds_mesh):
    """
    Build a dataset of fields on cells, edges and vertices with an extra
    ``nVertLevels`` dimension (and an ``nTracers`` dimension on one field) for
    testing how extra dimensions are unwrapped.
    """
    ncells = ds_mesh.sizes['nCells']
    nedges = ds_mesh.sizes['nEdges']
    nvertices = ds_mesh.sizes['nVertices']
    rng = np.random.default_rng(seed=42)

    ds = xr.Dataset()
    ds['temperature'] = (
        ('Time', 'nCells', 'nVertLevels'),
        rng.random((NTIME, ncells, NVERT_LEVELS)),
    )
    ds['tracers'] = (
        ('Time', 'nTracers', 'nCells', 'nVertLevels'),
        rng.random((NTIME, NTRACERS, ncells, NVERT_LEVELS)),
    )
    ds['bottomDepth'] = ('nCells', rng.random(ncells))
    ds['maxLevelCell'] = (
        'nCells',
        rng.integers(1, NVERT_LEVELS, ncells).astype(np.int32),
    )
    ds['normalVelocity'] = (
        ('Time', 'nEdges', 'nVertLevels'),
        rng.random((NTIME, nedges, NVERT_LEVELS)),
    )
    ds['vorticity'] = (
        ('Time', 'nVertices', 'nVertLevels'),
        rng.random((NTIME, nvertices, NVERT_LEVELS)),
    )
    return ds


def _convert_test_dataset(out_dir, **kwargs):
    """
    Convert the test dataset, keeping a non-contiguous subset of the vertical
    levels so that the read grouping has to handle gaps.
    """
    ds_mesh = xr.open_dataset(TEST_MESH)
    ds = _make_test_dataset(ds_mesh)
    extra_dims = {'nVertLevels': [0, 2, 4], 'nTracers': [0, 1]}
    converter = MpasToXdmf(ds=ds, ds_mesh=ds_mesh)
    converter.convert_to_xdmf(
        out_dir=str(out_dir), extra_dims=extra_dims, quiet=True, **kwargs
    )
    return ds


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_load_mesh_only():
    converter = MpasToXdmf()
    converter.load(mesh_filename=TEST_MESH)
    assert isinstance(converter.ds, xr.Dataset)
    assert isinstance(converter.ds_mesh, xr.Dataset)
    # Should have mesh dimensions
    assert 'nCells' in converter.ds.dims


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_set_time_with_no_xtime():
    converter = MpasToXdmf()
    converter.load(mesh_filename=TEST_MESH)
    # Should create a 'Time' variable if 'Time' in dims
    if 'Time' in converter.ds.dims:
        assert 'Time' in converter.ds
        arr = converter.ds['Time'].values
        assert np.all(arr == np.arange(converter.ds.sizes['Time']))


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_convert_to_xdmf(tmp_path):
    converter = MpasToXdmf()
    variables = ['xCell', 'areaCell', 'cellsOnCell']
    extra_dims = {'maxEdges': [0]}
    converter.load(mesh_filename=TEST_MESH, variables=variables)
    out_dir = tmp_path / 'out'
    converter.convert_to_xdmf(str(out_dir), extra_dims=extra_dims)
    # Check that output files exist for cells
    assert (out_dir / 'fieldsOnCells.h5').exists()
    assert (out_dir / 'fieldsOnCells.xdmf').exists()


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_extra_dims(tmp_path):
    converter = MpasToXdmf()
    converter.load(mesh_filename=TEST_MESH)
    # Simulate an extra dimension if present
    extra_dims = {}
    for dim in converter.ds.dims:
        if dim not in ['Time', 'nCells', 'nEdges', 'nVertices']:
            extra_dims[dim] = [0]
    out_dir = tmp_path / 'out_extra'
    converter.convert_to_xdmf(str(out_dir), extra_dims=extra_dims)
    assert (out_dir / 'fieldsOnCells.h5').exists()


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_load_with_time_series_and_variables(tmp_path):
    ts1 = tmp_path / 'ts1.nc'
    ts2 = tmp_path / 'ts2.nc'

    # Simulate a time series by adding xtime and area variables
    ds = xr.open_dataset(TEST_MESH)
    ds['xtime'] = ('Time', ['0001-01-01_00:00:00'])
    ds['area'] = (('Time', 'nCells'), ds.areaCell.values[None, :])
    write_netcdf(ds, ts1)
    ds['xtime'] = ('Time', ['0001-01-02_00:00:00'])
    write_netcdf(ds, ts2)

    variables = ['areaCell', 'area']

    converter = MpasToXdmf()
    converter.load(
        mesh_filename=TEST_MESH,
        time_series_filenames=[str(ts1), str(ts2)],
        variables=variables,
    )
    print(converter.ds)
    for var in variables:
        assert var in converter.ds.data_vars, (
            f'Variable {var} not found in dataset'
        )
    assert converter.ds.sizes['Time'] == 2


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_process_extra_dims_drop(tmp_path):
    converter = MpasToXdmf()
    converter.load(mesh_filename=TEST_MESH)

    # drop all variables with extra dimensions
    extra_dims = {
        'maxEdges': [],
        'maxEdges2': [],
        'TWO': [],
        'vertexDegree': [],
    }

    ds = _process_extra_dims(converter.ds, extra_dims=extra_dims)
    for dim in extra_dims:
        assert dim not in ds.dims, f'Dimension {dim} should be dropped'


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_extra_dims_unwrapped_values(tmp_path):
    """
    Each combination of extra-dimension indices becomes its own field, named
    with one ``_<index>`` suffix per extra dimension, and holds the values of
    the corresponding slice of the source variable.
    """
    out_dir = tmp_path / 'out_unwrap'
    ds = _convert_test_dataset(out_dir)

    with h5py.File(out_dir / 'fieldsOnCells.h5', 'r') as h5_file:
        names = set(h5_file.keys())
        for t_index in range(NTIME):
            for level in [0, 2, 4]:
                name = f'temperature_{level}_t{t_index}'
                expected = ds.temperature.isel(
                    Time=t_index, nVertLevels=level
                ).values
                assert np.array_equal(np.asarray(h5_file[name]), expected)
                for tracer in range(NTRACERS):
                    name = f'tracers_{level}_{tracer}_t{t_index}'
                    expected = ds.tracers.isel(
                        Time=t_index, nVertLevels=level, nTracers=tracer
                    ).values
                    assert np.array_equal(np.asarray(h5_file[name]), expected)
        # fields without extra dimensions keep their original name
        assert np.array_equal(
            np.asarray(h5_file['bottomDepth']), ds.bottomDepth.values
        )

    # levels that were not requested are not written
    assert 'temperature_1_t0' not in names
    assert 'temperature_3_t0' not in names


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_vertex_fields_mapped_to_kites(tmp_path):
    """
    Vertex-centered fields are repeated once per kite of the dual mesh, so
    they are longer than the ``nVertices`` dimension of the source data.
    """
    out_dir = tmp_path / 'out_vertices'
    ds = _convert_test_dataset(out_dir)

    with h5py.File(out_dir / 'fieldsOnVertices.h5', 'r') as h5_file:
        cells = np.asarray(h5_file['Cells'])
        field = np.asarray(h5_file['vorticity_2_t1'])
        assert field.shape == (cells.shape[0],)
        # every written value comes from the source field
        source = ds.vorticity.isel(Time=1, nVertLevels=2).values
        assert np.all(np.isin(field, source))


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_max_read_bytes_does_not_change_results(tmp_path):
    """
    ``max_read_bytes`` only controls how many indices are read at a time, so a
    limit small enough to force one read per index must give the same output as
    a limit large enough to read every index at once.
    """
    one_pass = tmp_path / 'out_one_pass'
    many_passes = tmp_path / 'out_many_passes'
    _convert_test_dataset(one_pass, max_read_bytes=1 << 30)
    # 1 byte forces the minimum of one index of the first extra dimension per
    # read
    _convert_test_dataset(many_passes, max_read_bytes=1)

    for basename in ['fieldsOnCells.h5', 'fieldsOnEdges.h5']:
        with (
            h5py.File(one_pass / basename, 'r') as expected,
            h5py.File(many_passes / basename, 'r') as actual,
        ):
            assert set(expected.keys()) == set(actual.keys())
            for name in expected:
                assert np.array_equal(
                    np.asarray(expected[name]), np.asarray(actual[name])
                ), f'{basename}:/{name} differs'

    for basename in ['fieldsOnCells.xdmf', 'fieldsOnEdges.xdmf']:
        assert (one_pass / basename).read_text() == (
            many_passes / basename
        ).read_text()


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_float32(tmp_path):
    """
    ``float32=True`` writes floating-point fields in single precision but
    leaves integer fields and the mesh geometry alone.
    """
    out_dir = tmp_path / 'out_float32'
    ds = _convert_test_dataset(out_dir, float32=True)

    with h5py.File(out_dir / 'fieldsOnCells.h5', 'r') as h5_file:
        assert h5_file['temperature_0_t0'].dtype == np.float32
        assert h5_file['bottomDepth'].dtype == np.float32
        assert h5_file['maxLevelCell'].dtype == np.int32
        # the geometry stays in double precision so cell shapes are unaffected
        assert h5_file['Points'].dtype == np.float64
        expected = ds.temperature.isel(Time=0, nVertLevels=0).values
        assert np.array_equal(
            np.asarray(h5_file['temperature_0_t0']),
            expected.astype(np.float32),
        )


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_process_extra_dims_keeps_selected_dims(tmp_path):
    """
    Variables with selected indices keep their extra dimensions, which are
    unwrapped as the variables are written rather than up front.
    """
    converter = MpasToXdmf()
    converter.load(mesh_filename=TEST_MESH)

    extra_dims = {
        'maxEdges': [0, 1],
        'maxEdges2': [],
        'TWO': [],
        'vertexDegree': [],
    }
    ds = _process_extra_dims(converter.ds, extra_dims=extra_dims)

    assert 'maxEdges' in ds.dims
    assert 'verticesOnCell' in ds.data_vars
    assert 'maxEdges' in ds.verticesOnCell.dims
    for dim in ['maxEdges2', 'TWO', 'vertexDegree']:
        assert dim not in ds.dims, f'Dimension {dim} should be dropped'


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_process_extra_dims_uncovered_dim():
    """An extra dimension that is not listed in ``extra_dims`` is an error."""
    converter = MpasToXdmf()
    converter.load(mesh_filename=TEST_MESH)
    with pytest.raises(ValueError, match='is not covered'):
        _process_extra_dims(converter.ds, extra_dims={'maxEdges': [0]})


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_set_time_invalid_xtime(tmp_path):
    ts1 = tmp_path / 'ts1.nc'
    # Simulate a time-depndent variable and add xtime
    ds = xr.open_dataset(TEST_MESH)
    ds['xtime'] = ('Time', ['0001-01-01_00:00:00'])
    ds['area'] = (('Time', 'nCells'), ds.areaCell.values[None, :])
    write_netcdf(ds, ts1)

    converter = MpasToXdmf()
    converter.load(mesh_filename=TEST_MESH, time_series_filenames=[str(ts1)])
    # Should raise ValueError if xtime_var is not present
    with pytest.raises(ValueError):
        _set_time(ds=converter.ds, xtime_var='not_a_var')


def test_parse_indices_invalid_cases():
    # Should raise on mixed slice/list
    with pytest.raises(ValueError):
        _parse_indices('1:3,5', 5)
    # Should raise on invalid string
    with pytest.raises(ValueError):
        _parse_indices('foo', 5)


def test_parse_indices_valid_cases():
    # Empty list
    assert _parse_indices('', 5) == []
    # Single index
    assert _parse_indices('0', 5) == [0]
    # Comma-separated list
    assert _parse_indices('1,2,3', 5) == [1, 2, 3]
    # Slice notation
    assert _parse_indices('0:3', 5) == [0, 1, 2]
    # Slice with stride
    assert _parse_indices('0:5:2', 5) == [0, 2, 4]
    # Full slice
    assert _parse_indices(':', 4) == [0, 1, 2, 3]


def test_main_cli(monkeypatch, tmp_path):
    # Test CLI entry point with minimal arguments
    mesh = TEST_MESH
    if not os.path.exists(mesh):
        pytest.skip('Test mesh not available')
    out_dir = tmp_path / 'cli_out'
    sys_argv = ['prog', '-m', mesh, '-o', str(out_dir), '-v', 'areaCell']
    monkeypatch.setattr(sys, 'argv', sys_argv)
    # Patch input to always return blank (skip extra dims)
    monkeypatch.setattr('builtins.input', lambda _: '')
    main()
    assert (out_dir / 'fieldsOnCells.h5').exists()


@pytest.mark.skipif(
    not os.path.exists(TEST_MESH), reason='Test mesh not available'
)
def test_load_dataset_missing_variable():
    # Should not raise if variable is missing in mesh, but should raise if not
    # present at all
    with pytest.raises(KeyError):
        _load_dataset(TEST_MESH, None, ['not_a_var'], None)
