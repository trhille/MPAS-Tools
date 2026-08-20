#!/usr/bin/env python

import numpy as np
import xarray as xr

from mpas_tools.vector.reconstruct import reconstruct_variable


def test_reconstruct_planar():
    """
    On a planar mesh, the zonal and meridional components are the x and y
    components, with no rotation applied.
    """
    ds_mesh, coeffs, var_on_edges = _get_mesh(on_a_sphere='NO')

    ds_out = xr.Dataset()
    reconstruct_variable(
        'velocity', var_on_edges, ds_mesh, coeffs, ds_out, quiet=True
    )

    assert np.allclose(ds_out.velocityZonal, ds_out.velocityX)
    assert np.allclose(ds_out.velocityMeridional, ds_out.velocityY)


def test_reconstruct_planar_without_lat_lon():
    """
    A planar mesh need not have ``latCell`` and ``lonCell`` at all.
    """
    ds_mesh, coeffs, var_on_edges = _get_mesh(on_a_sphere='NO')
    ds_mesh = ds_mesh.drop_vars(['latCell', 'lonCell'])

    ds_out = xr.Dataset()
    reconstruct_variable(
        'velocity', var_on_edges, ds_mesh, coeffs, ds_out, quiet=True
    )

    assert np.allclose(ds_out.velocityZonal, ds_out.velocityX)
    assert np.allclose(ds_out.velocityMeridional, ds_out.velocityY)


def test_reconstruct_spherical():
    """
    On a spherical mesh, the Cartesian components are rotated into the local
    geographic frame at each cell center.
    """
    ds_mesh, coeffs, var_on_edges = _get_mesh(on_a_sphere='YES')

    ds_out = xr.Dataset()
    reconstruct_variable(
        'velocity', var_on_edges, ds_mesh, coeffs, ds_out, quiet=True
    )

    lat = ds_mesh.latCell
    lon = ds_mesh.lonCell
    u_x = ds_out.velocityX
    u_y = ds_out.velocityY
    u_z = ds_out.velocityZ
    zonal = -u_x * np.sin(lon) + u_y * np.cos(lon)
    merid = -(u_x * np.cos(lon) + u_y * np.sin(lon)) * np.sin(
        lat
    ) + u_z * np.cos(lat)

    assert np.allclose(ds_out.velocityZonal, zonal)
    assert np.allclose(ds_out.velocityMeridional, merid)


def test_reconstruct_missing_on_a_sphere():
    """
    A mesh without the ``on_a_sphere`` attribute is treated as spherical, as
    it was before planar meshes were handled.
    """
    ds_mesh, coeffs, var_on_edges = _get_mesh(on_a_sphere='YES')
    ds_with_attr = ds_mesh
    ds_without_attr = ds_mesh.copy()
    del ds_without_attr.attrs['on_a_sphere']

    ds_ref = xr.Dataset()
    reconstruct_variable(
        'velocity', var_on_edges, ds_with_attr, coeffs, ds_ref, quiet=True
    )
    ds_out = xr.Dataset()
    reconstruct_variable(
        'velocity', var_on_edges, ds_without_attr, coeffs, ds_out, quiet=True
    )

    assert np.allclose(ds_out.velocityZonal, ds_ref.velocityZonal)
    assert np.allclose(ds_out.velocityMeridional, ds_ref.velocityMeridional)


def test_reconstruct_on_a_sphere_padded():
    """
    MPAS files often pad the ``on_a_sphere`` attribute with spaces.
    """
    ds_mesh, coeffs, var_on_edges = _get_mesh(on_a_sphere='NO              ')

    ds_out = xr.Dataset()
    reconstruct_variable(
        'velocity', var_on_edges, ds_mesh, coeffs, ds_out, quiet=True
    )

    assert np.allclose(ds_out.velocityZonal, ds_out.velocityX)
    assert np.allclose(ds_out.velocityMeridional, ds_out.velocityY)


def _get_mesh(on_a_sphere):
    """
    A tiny synthetic mesh, along with reconstruction coefficients and a
    variable on edges, for testing ``reconstruct_variable()``
    """
    n_cells = 4
    max_edges = 6
    n_edges = n_cells * max_edges
    rng = np.random.default_rng(seed=0)

    ds_mesh = xr.Dataset()
    ds_mesh['edgesOnCell'] = (
        ('nCells', 'maxEdges'),
        np.arange(1, n_edges + 1).reshape(n_cells, max_edges),
    )
    if on_a_sphere.strip() == 'NO':
        # as in meshes from ``mpas_tools.planar_hex``
        lat = np.zeros(n_cells)
        lon = np.zeros(n_cells)
    else:
        lat = np.deg2rad(np.array([-80.0, -20.0, 15.0, 70.0]))
        lon = np.deg2rad(np.array([10.0, 135.0, 200.0, 350.0]))
    ds_mesh['latCell'] = (('nCells',), lat)
    ds_mesh['lonCell'] = (('nCells',), lon)
    ds_mesh.attrs['on_a_sphere'] = on_a_sphere

    coeffs = xr.DataArray(
        rng.random((n_cells, max_edges, 3)),
        dims=('nCells', 'maxEdges', 'R3'),
    )
    var_on_edges = xr.DataArray(rng.random(n_edges), dims=('nEdges',))

    return ds_mesh, coeffs, var_on_edges
