#!/usr/bin/env python
"""Collection of simple functions useful in computational chemistry scripting.

Many of the following functions are used to make operations on xyz coordinates
of molecular structure. When refering to ``xyz_data`` bellow, the following
structures (also used in :py:mod:`~comp_chem_utils.molecule_data`) is assumed::

    atom 1 label  and corresponding xyz coordinate
    atom 2 label  and corresponding xyz coordinate
    : : :
    atom N label  and corresponding xyz coordinate

For example the ``xyz_data`` of a Hydrogen molecule along the z-axis
should be passed as::

    >>> xyz_data
    [['H', 0.0, 0.0, 0.0], ['H', 0.0, 0.0, 1.0]]

"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"


import os
import shutil
import numpy as np
from scipy.spatial.distance import pdist

from comp_chem_utils.periodic import element


def vel_auto_corr(vel, max_corr_index, tstep):
    """Calculate velocity autocorrelation function.

    Args:
        vel (list): The velocities along the trajectory are given as
            a list of ``np.array()`` of dimensions N_atoms . 3.

        max_corr_index (int): Maximum number of steps to consider
            for the auto correlation function. In other words, it
            corresponds to the number of points in the output function.
    """

    max_step_index = len(vel) - max_corr_index
    natoms = vel[0].shape[0]

    G = np.zeros((max_corr_index), dtype='float64')


    for itau in range(max_corr_index):
        for it in range(max_step_index):
            #for i in range(natoms):
            #    G[itau] += np.dot(vel[it][i,:], vel[it+itau][i,:])

            G[itau] += np.trace(np.dot(vel[it],np.transpose(vel[it+itau])))

    G = G / (natoms * max_corr_index * tstep)
    xpts = np.arange(max_corr_index)
    xpts = xpts * tstep

    return xpts, G


def get_lmax_from_atomic_charge(charge):
    """Return the maximum angular momentum based on the input atomic charge.

    This function is designed to return LMAX in a CPMD input file.

    Args:
        charge (int): Atomic charge.

    Returns:
        'S' for H or He; 'P' for second row elements; and 'D' for heavier elements.
    """

    if charge <= 2:
        # H or He
        lmax = 'S'
    elif charge <= 10:
        lmax = 'P'
    else:
        lmax = 'D'

    return lmax


def get_file_as_list(filename, raw=False):
    """Read a file and return it as a list of lines (str).

    By default comments (i.e. lines starting with #)
    and empty lines are ommited. This can be changed
    by setting ``raw=True``

    Args:
        filename (str): Name of the file to read.
        raw (bool, optional): To return the file as it is,
            i.e. with comments and blank lines. Default
            is ``raw=False``.

    Returns:
        A list of lines (str)."""

    lines = []
    with open(filename,'r') as myfile:
        for line in myfile:
            if raw:
                lines.append(line)
            else:
                # remove empty lines
                if line.strip():
                    # remove comments
                    if line.strip()[0] != '#':
                        lines.append(line)

    return lines


def make_new_dir(dirn):
    """Make new empty directory.

    If the directory already exists it is erased and replaced.

    Args:
        dirn (str): Name for the new directory (can include path).
    """

    if not os.path.exists(dirn):
        os.makedirs(dirn)
    else:
        try:
            os.removedirs(dirn)
        except(OSError):
            print("WARNING: erasing (not empty) directory! {}".format(dirn))
            shutil.rmtree(dirn)

        os.makedirs(dirn)


def center_of_mass(xyz_data):
    """Calculate center of mass of a molecular structure based on xyz coordinates.

    Args:
        xyz_data (list): xyz atomic coordinates arranged as described above.

    Returns:
        3-dimensional ``np.array()`` containing the xyz coordinates of the
        center of mass of the molecular structure. The unit of the center of mass
        matches the xyz input units (usually Angstroms).
    """

    totM = 0.0
    COM = np.zeros((3))
    for line in xyz_data:
        symbol = line[0]
        coord = np.array(line[1:])

        mass = element(symbol).mass
        totM += mass
        COM += coord*mass

    COM = COM/totM
    return COM


def change_vector_norm(fix, mob, R):
    """Scale a 3-D vector defined by two points in space to have a new norm R.

    The input vector is defined by a fix position in 3-D space ``fix``,
    and a mobile position ``mob``. The function returns a new mobile
    position such that the new vector has the norm R.

    Args:
        fix (np.array): xyz coordinates of the fix point.
        mob (np.array): Original xyz coordinates of the mobile point.
        R (float): Desired norm for the new vector.

    Returns:
        The new mobile position as an ``np.array()`` of dimenssion 3.
    """

    unit = mob - fix
    unit = unit/np.linalg.norm(unit)

    # return new position
    return fix + R * unit


def change_vector_norm_sym(pos1, pos2, R):
    """Symmetric version of change_vector_norm function.

    In other word both positions are modified symmetrically.
    """

    unit = pos2 - pos1
    norm = np.linalg.norm(unit)
    unit = unit/norm

    shift = (R-norm)/2.0

    # return new position
    new_pos1 = pos1 - unit * shift
    new_pos2 = pos2 + unit * shift
    return new_pos1, new_pos2


def get_rmsd(xyz_data1, xyz_data2):
    """Calculate RMSD between two sets of coordinates.

    The Root-mean-square deviation of atomic positions is calculated as

    .. math::
        RMSD = \\sqrt{ \\frac{1}{N} \\sum_{i=1}^N \\delta_{i}^{2} }

    Where ``\delta_i`` is the distance between atom i in ``xyz_data1`` and in
    ``xyz_data2``.

    Args:
        xyz_data1 (list): List of atomic coordinates for the first structure
            arranged as described above for xyz_data.
        xyz_data2 (list): Like ``xyz_data1`` but for the second structure.

    Returns:
        The RMSD (float).
    """

    rmsd = 0
    for c1, c2 in zip(xyz_data1, xyz_data2):

        d1 = np.array([c1[x] for x in range(1,4)])
        d2 = np.array([c2[x] for x in range(1,4)])
        vector = d2 - d1
        rmsd += np.dot(vector, vector)

    rmsd = rmsd/(len(coord1))
    return np.sqrt(rmsd)


def get_distance(xyz_data, atoms, box_size=None):
    """Calculate distance between two atoms in xyz_data.

    Args:
        xyz_data (list): xyz atomic coordinates arranged as described above.
        atoms (list): list of two indices matching the two rows of the
            xyz_data for wich the distance should be calculated.

    Returns:
        Distance between the two atoms in the list as a ``float``, the
        unit will match the input unit in the ``xyz_data``.

    """
    coord1 = np.array([xyz_data[atoms[0]][x] for x in range(1,4)])
    coord2 = np.array([xyz_data[atoms[1]][x] for x in range(1,4)])
    vector = coord2 - coord1

    if box_size:
       for i,x in enumerate(vector):
          if abs(x) > box_size/2.0:
             vector[i] = box_size - abs(x)

    return np.linalg.norm(vector)


def get_distance_matrix(xyz_data, box_size=None):

   # make np.array
   natoms = len(xyz_data)
   coord = np.zeros((natoms,3), dtype='float')
   for i,line in enumerate(xyz_data):
      coord[i,:] = line[1:]

   if box_size:
      npairs = natoms * (natoms - 1)
      matrix = np.zeros((npairs,3), dtype='float')

      for i in range(natoms):
         for j in range(i+1,natoms):
            ij = i + (j-1)*natoms

            matrix[ij,:] = coord[i,:] - coord[j,:]

      # find out which element to shift:
      # basically shift is an array of same shape as matrix
      # with zeros every where except where the elements of
      # matrix are larger than box_size/2.0
      # in that case shift as the value box_size
      shift = box_size * (matrix > box_size/2.0).astype(int)

      # we can now shift the matrix as follows:
      matrix = abs(shift - matrix)

      # and get the distances...
      matrix = np.linalg.norm(matrix, axis=1)

   else:
      matrix = pdist(coord)

   return matrix

def get_distance_matrix_2(xyz_data1, xyz_data2, box_size=None):
   # repeat as above for 2 different sets of coordinates
   nat1 = len(xyz_data1)
   coord1 = np.zeros((nat1,3), dtype='float')
   for i,line in enumerate(xyz_data1):
      coord1[i,:] = line[1:]

   nat2 = len(xyz_data2)
   coord2 = np.zeros((nat2,3), dtype='float')
   for i,line in enumerate(xyz_data2):
      coord2[i,:] = line[1:]

   if box_size:
      matrix = np.zeros((nat1,nat2,3), dtype='float')

      for i in range(nat1):
         for j in range(nat2):
            matrix[i,j,:] = coord[i,:] - coord[j,:]

      # find out which element to shift:
      # basically shift is an array of same shape as matrix
      # with zeros every where except where the elements of
      # matrix are larger than box_size/2.0
      # in that case shift as the value box_size
      shift = box_size * (matrix > box_size/2.0).astype(int)

      # we can now shift the matrix as follows:
      matrix = abs(shift - matrix)

      # and get the distances...
      matrix = np.linalg.norm(matrix, axis=2)

   else:
      matrix = cdist(coord1, coord2)

   return matrix

def get_angle(xyz_data, atoms):
    """Calculate angle between three atoms in xyz_data.

    Args:
        xyz_data (list): xyz atomic coordinates arranged as described above.
        atoms (list): list of three indices matching the rows of the
            xyz_data for wich the angle should be calculated.

    Returns:
        Angle between the three atoms in the list as a ``float`` in degrees.

    """

    coord1 = np.array([xyz_data[atoms[0]][x] for x in range(1,4)])
    coord2 = np.array([xyz_data[atoms[1]][x] for x in range(1,4)])
    coord3 = np.array([xyz_data[atoms[2]][x] for x in range(1,4)])
    vec1 = coord1-coord2
    vec2 = coord3-coord2

    return np.degrees( np.arccos( np.dot(vec1,vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)) ) )

def get_dihedral_angle(table, atoms):
    """Calculate dihedral angle defined by 4 atoms in xyz_data.

    It relies on the praxeolitic formula (1 sqrt, 1 cross product).

    Args:
        xyz_data (list): xyz atomic coordinates arranged as described above.
        atoms (list): list of 4 indices matching the rows of the
            xyz_data for wich the dihedral angle should be calculated.

    Returns:
        Dihedral angle defined by the 4 atoms in the list as a ``float`` in degrees.

    """

    p0 = np.array([table[atoms[0]][x] for x in range(1,4)])
    p1 = np.array([table[atoms[1]][x] for x in range(1,4)])
    p2 = np.array([table[atoms[2]][x] for x in range(1,4)])
    p3 = np.array([table[atoms[3]][x] for x in range(1,4)])

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    return np.degrees(np.arctan2(y, x))


