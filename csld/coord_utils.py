#!/usr/bin/env python3

"""
Utilities for manipulating coordinates or list of coordinates, under periodic
boundary conditions or otherwise. Many of these are heavily vectorized in
numpy for performance.
"""

# adapted from original version in pymatgen version from pymatgen

import numpy as np
import math

def small_fractional(x):
    return x-np.round(x)

def find_in_coord_list(coord_list, coord, atol=1e-8):
    """
    Find the indices of matches of a particular coord in a coord_list.

    Args:
        coord_list: List of coords to test
        coord: Specific coordinates
        atol: Absolute tolerance. Defaults to 1e-8. Accepts both scalar and
            array.

    Returns:
        Indices of matches, e.g., [0, 1, 2, 3]. Empty list if not found.
    """
    if len(coord_list) == 0:
        return []
    diff = np.array(coord_list) - np.array(coord)[None, :]
    return np.where(np.all(np.abs(diff) < atol, axis=1))[0]


def in_coord_list(coord_list, coord, atol=1e-8):
    """
    Tests if a particular coord is within a coord_list.

    Args:
        coord_list: List of coords to test
        coord: Specific coordinates
        atol: Absolute tolerance. Defaults to 1e-8. Accepts both scalar and
            array.

    Returns:
        True if coord is in the coord list.
    """
    return len(find_in_coord_list(coord_list, coord, atol=atol)) > 0


def is_coord_subset(subset, superset, atol=1e-8):
    """
    Tests if all coords in subset are contained in superset.
    Doesn't use periodic boundary conditions

    Args:
        subset, superset: List of coords

    Returns:
        True if all of subset is in superset.
    """
    c1 = np.array(subset)
    c2 = np.array(superset)
    is_close = np.all(np.abs(c1[:, None, :] - c2[None, :, :]) < atol, axis=-1)
    any_close = np.any(is_close, axis=-1)
    return np.all(any_close)


def coord_list_mapping(subset, superset):
    """
    Gives the index mapping from a subset to a superset.
    Subset and superset cannot contain duplicate rows

    Args:
        subset, superset: List of coords

    Returns:
        list of indices such that superset[indices] = subset
    """
    c1 = np.array(subset)
    c2 = np.array(superset)
    inds = np.where(np.all(np.isclose(c1[:, None, :], c2[None, :, :]),
                           axis=2))[1]
    result = c2[inds]
    if not np.allclose(c1, result):
        if not is_coord_subset(subset, superset):
            raise ValueError("subset is not a subset of superset")
    if not result.shape == c1.shape:
        raise ValueError("Something wrong with the inputs, likely duplicates "
                         "in superset")
    return inds


def coord_list_mapping_pbc(subset, superset, atol=1e-8):
    """
    Gives the index mapping from a subset to a superset.
    Subset and superset cannot contain duplicate rows

    Args:
        subset, superset: List of frac_coords

    Returns:
        list of indices such that superset[indices] = subset
    """
    c1 = np.array(subset)
    c2 = np.array(superset)

    diff = c1[:, None, :] - c2[None, :, :]
    diff -= np.round(diff)
    inds = np.where(np.all(np.abs(diff) < atol, axis = 2))[1]

    #verify result (its easier to check validity of the result than
    #the validity of inputs)
    test = c2[inds] - c1
    test -= np.round(test)
    if not np.allclose(test, 0):
        if not is_coord_subset_pbc(subset, superset):
            raise ValueError("subset is not a subset of superset")
    if not test.shape == c1.shape:
        raise ValueError("Something wrong with the inputs, likely duplicates "
                         "in superset")
    return inds


def get_linear_interpolated_value(x_values, y_values, x):
    """
    Returns an interpolated value by linear interpolation between two values.
    This method is written to avoid dependency on scipy, which causes issues on
    threading servers.

    Args:
        x_values: Sequence of x values.
        y_values: Corresponding sequence of y values
        x: Get value at particular x

    Returns:
        Value at x.
    """
    a = np.array(sorted(zip(x_values, y_values), key=lambda d: d[0]))

    ind = np.where(a[:, 0] >= x)[0]

    if len(ind) == 0 or ind[0] == 0:
        raise ValueError("x is out of range of provided x_values")

    i = ind[0]
    x1, x2 = a[i - 1][0], a[i][0]
    y1, y2 = a[i - 1][1], a[i][1]

    return y1 + (y2 - y1) / (x2 - x1) * (x - x1)


def all_distances(coords1, coords2):
    """
    Returns the distances between two lists of coordinates

    Args:
        coords1: First set of cartesian coordinates.
        coords2: Second set of cartesian coordinates.

    Returns:
        2d array of cartesian distances. E.g the distance between
        coords1[i] and coords2[j] is distances[i,j]
    """
    c1 = np.array(coords1)
    c2 = np.array(coords2)
    z = (c1[:, None, :] - c2[None, :, :]) ** 2
    return np.sum(z, axis=-1) ** 0.5


def pbc_diff(fcoords1, fcoords2):
    """
    Returns the 'fractional distance' between two coordinates taking into
    account periodic boundary conditions.

    Args:
        fcoords1: First set of fractional coordinates. e.g., [0.5, 0.6,
            0.7] or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
            coord or any array of coords.
        fcoords2: Second set of fractional coordinates.

    Returns:
        Fractional distance. Each coordinate must have the property that
        abs(a) <= 0.5. Examples:
        pbc_diff([0.1, 0.1, 0.1], [0.3, 0.5, 0.9]) = [-0.2, -0.4, 0.2]
        pbc_diff([0.9, 0.1, 1.01], [0.3, 0.5, 0.9]) = [-0.4, -0.4, 0.11]
    """
    fdist = np.subtract(fcoords1, fcoords2)
    return fdist - np.round(fdist)

def pbc_images(nogamma=False):
    r = np.arange(-1, 2)
    arange = r[:, None] * np.array([1, 0, 0])[None, :]
    brange = r[:, None] * np.array([0, 1, 0])[None, :]
    crange = r[:, None] * np.array([0, 0, 1])[None, :]
    images = arange[:, None, None] + brange[None, :, None] + \
        crange[None, None, :]
    images = images.reshape((27,3))
    return np.delete(images, 13, 0) if nogamma else images



def pbc_shortest_vectors(lattice, fcoords1, fcoords2):
    """
    Returns the shortest vectors between two lists of coordinates taking into
    account periodic boundary conditions and the lattice.

    Args:
        lattice: lattice to use
        fcoords1: First set of fractional coordinates. e.g., [0.5, 0.6, 0.7]
            or [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]]. It can be a single
            coord or any array of coords.
        fcoords2: Second set of fractional coordinates.

    Returns:
        array of displacement vectors from fcoords1 to fcoords2
        first index is fcoords1 index, second is fcoords2 index
    """
    #ensure correct shape
    fcoords1, fcoords2 = np.atleast_2d(fcoords1, fcoords2)

    #ensure that all points are in the unit cell
    fcoords1 = np.mod(fcoords1, 1)
    fcoords2 = np.mod(fcoords2, 1)

    images = pbc_images()

    #create images of f2
    shifted_f2 = fcoords2[:, None, :] + images[None, :, :]

    cart_f1 = lattice.get_cartesian_coords(fcoords1)
    cart_f2 = lattice.get_cartesian_coords(shifted_f2)

    #all vectors from f1 to f2
    vectors = cart_f2[None, :, :, :] - cart_f1[:, None, None, :]

    d_2 = np.sum(vectors ** 2, axis=3)
    a, b = np.indices([len(fcoords1), len(fcoords2)])
    return vectors[a, b, np.argmin(d_2, axis=2)]


def find_in_coord_list_pbc(fcoord_list, fcoord, atol=1e-8):
    """
    Get the indices of all points in a fractional coord list that are
    equal to a fractional coord (with a tolerance), taking into account
    periodic boundary conditions.

    Args:
        fcoord_list: List of fractional coords
        fcoord: A specific fractional coord to test.
        atol: Absolute tolerance. Defaults to 1e-8.

    Returns:
        Indices of matches, e.g., [0, 1, 2, 3]. Empty list if not found.
    """
    if len(fcoord_list) == 0:
        return []
    fcoords = np.tile(fcoord, (len(fcoord_list), 1))
    fdist = fcoord_list - fcoords
    fdist -= np.round(fdist)
    return np.where(np.all(np.abs(fdist) < atol, axis=1))[0]


def in_coord_list_pbc(fcoord_list, fcoord, atol=1e-8):
    """
    Tests if a particular fractional coord is within a fractional coord_list.

    Args:
        fcoord_list: List of fractional coords to test
        fcoord: A specific fractional coord to test.
        atol: Absolute tolerance. Defaults to 1e-8.

    Returns:
        True if coord is in the coord list.
    """
    return len(find_in_coord_list_pbc(fcoord_list, fcoord, atol=atol)) > 0


def is_coord_subset_pbc(subset, superset, atol=1e-8):
    """
    Tests if all fractional coords in subset are contained in superset.

    Args:
        subset, superset: List of fractional coords

    Returns:
        True if all of subset is in superset.
    """
    c1 = np.array(subset)
    c2 = np.array(superset)
    dist = c1[:, None, :] - c2[None, :, :]
    dist -= np.round(dist)
    is_close = np.all(np.abs(dist) < atol, axis=-1)
    any_close = np.any(is_close, axis=-1)
    return np.all(any_close)


# def lattice_points_in_supercell_ijk(sc_mat)
#     """
#     sc_mat: supercell defining matrix
#     return lattice points in integer [i j k], i.e. fractional coordinates relative to primitive cell
#     """
#     dim = len(SCmat)
#     lat=Flatten[Outer[List, Sequence@@(Range[Min[#], Max[#]]&/@Transpose[Join[Total[SCmat[[#]]]&/@DeleteCases[Subsets[Range[dim]], {}], {ConstantArray[0,{dim}]}]])],dim-1];
# QSCmat=Inverse[Transpose[SCmat]];
# lat=Select[lat, (And@@(Function[z, (z>=0)&&(z<1)]/@(QSCmat.#)))&];
# If[Length[lat]!= Abs@Det[SCmat], Print["ERROR: found ", Length[lat], " lattice points in supercell ", SCmat]];
# lat
# ];



def lattice_points_in_supercell(supercell_matrix, relative_to_prim=False):
    """
    Returns the list of points on the original lattice contained in the
    supercell in fractional coordinates (with the supercell basis).
    e.g. [[2,0,0],[0,1,0],[0,0,1]] returns [[0,0,0],[0.5,0,0]]

    Args:
        supercell_matrix: 3x3 matrix describing the supercell

    Returns:
        numpy array of the fractional coordinates
    """
    diagonals = np.array(
        [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1],
         [1, 1, 0], [1, 1, 1]])
    d_points = np.dot(diagonals, supercell_matrix)

    minimax = np.array([np.min(d_points, axis=0), np.max(d_points, axis=0) + 1])

    ar = np.arange(minimax[0, 0], minimax[1, 0])[:, None] * \
         np.array([1, 0, 0])[None, :]
    br = np.arange(minimax[0, 1], minimax[1, 1])[:, None] * \
         np.array([0, 1, 0])[None, :]
    cr = np.arange(minimax[0, 2], minimax[1, 2])[:, None] * \
         np.array([0, 0, 1])[None, :]

    all_points = ar[:, None, None] + br[None, :, None] + cr[None, None, :]
    all_points = all_points.reshape((-1, 3))

    frac_points = np.dot(all_points, inv_sc = np.linalg.inv(supercell_matrix))

    tvects = frac_points[np.where(np.all(frac_points < 1 - 1e-10, axis=1)
                                  & np.all(frac_points >= -1e-10, axis=1))]
    assert len(tvects) == np.round(np.abs(np.linalg.det(supercell_matrix)))
    if relative_to_prim:
        return tvects
    else:
        return np.around(np.dot(tvects), supercell_matrix)


def barycentric_coords(coords, simplex):
    """
    Converts a list of coordinates to barycentric coordinates, given a
    simplex with d+1 points. Only works for d >= 2.

    Args:
        coords: list of n coords to transform, shape should be (n,d)
        simplex: list of coordinates that form the simplex, shape should be
            (d+1, d)

    Returns:
        a LIST of barycentric coordinates (even if the original input was 1d)
    """
    coords = np.atleast_2d(coords)

    t = np.transpose(simplex[:-1, :]) - np.transpose(simplex[-1, :])[:, None]
    all_but_one = np.transpose(
        np.linalg.solve(t, np.transpose(coords - simplex[-1])))
    last_coord = 1 - np.sum(all_but_one, axis=-1)[:, None]
    return np.append(all_but_one, last_coord, axis=-1)


def get_angle(v1, v2, units="degrees"):
    """
    Calculates the angle between two vectors.

    Args:
        v1: Vector 1
        v2: Vector 2
        units: "degrees" or "radians". Defaults to "degrees".

    Returns:
        Angle between them in degrees.
    """
    d = np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)
    d = min(d, 1)
    d = max(d, -1)
    angle = math.acos(d)
    if units == "degrees":
        return math.degrees(angle)
    elif units == "radians":
        return angle
    else:
        raise ValueError("Invalid units {}".format(units))


def supercell_latticepoints(sc_matrix):
    """
    return supercell lattice points

    Args:
        sc_matrix: A scaling matrix for transforming the lattice
            vectors. Has to be all integers.
    """
    def range_vec(i):
        low = 0
        high = 0
        # print("debug ", i, sc_matrix[:, i])
        for z in sc_matrix[:, i]:
            if z > 0:
                high += z
            else:
                low += z
        return np.arange(low, high + 1)

    arange = range_vec(0)[:, None] * np.array([1, 0, 0])[None, :]
    brange = range_vec(1)[:, None] * np.array([0, 1, 0])[None, :]
    crange = range_vec(2)[:, None] * np.array([0, 0, 1])[None, :]
    all_points = arange[:, None, None] + brange[None, :, None] + \
                 crange[None, None, :]
    all_points = all_points.reshape((-1, 3))

    # find the translation vectors (in terms of the initial lattice vectors)
    #that are inside the unit cell defined by the scale matrix
    #we're using a slightly offset interval from 0 to 1 to avoid numerical
    #precision issues
    frac_points = np.dot(all_points, np.linalg.inv(sc_matrix))
    tvects = all_points[np.where(np.all(frac_points < 1 - 1e-10, axis=1)
                                 & np.all(frac_points >= -1e-10, axis=1))]
    assert len(tvects) == np.round(abs(np.linalg.det(sc_matrix)))
    return tvects

def match_p1p0(p1, p0):
    """

    :param p1:
    :param p0:
    :return:
    """
    dp = np.mod(p1 - p0, 1)
    for i in range(p0.shape[0]):
        for j in range(p0.shape[1]):
            if dp[i, j] >= 0.49:
                dp[i, j] -= 1
            elif dp[i, j] <= -0.49:
                dp[[i, j]] += 1

    return p0 + dp


def ReadPBC2Cart(fn_p1, p0, tol=1e-13):
    """

    :param fn_p1: filename of p1
    :param p0: ideal coordinates
     R: lattice vectors
    :return:
    """
    from .interface_vasp import Poscar
    from .util.mathtool import mychop
    str1 = Poscar.from_file(fn_p1).structure
    p1 = str1.frac_coords
    return mychop(np.dot(match_p1p0(p1, p0)-p0, str1.lattice.matrix), tol)

