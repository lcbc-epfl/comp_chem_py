#!/usr/bin/env python
"""Collection of simple functions useful in computational chemistry scripting."""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"


import os
import shutil
import numpy as np

from periodic import element


def get_lmax_from_atomic_charge(charge):
    if charge <= 2:
        # H or He
        lmax = 'S'
    elif charge <= 10:
        lmax = 'P'
    else:
        lmax = 'D'

    return lmax

def get_file_as_list(filename, raw=False):
    """Read a file and return it as a list.

    By default comments (i.e. lines starting with #)
    and empty lines are ommited. This can be changed
    by setting ``raw=True``

    :param str filename: Name of the file to read.
    :param bool raw: Return file as it is.

    :returns: A list of lines (str)."""

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
    If the directory already exists it is erased and replaced
    
    :param str dirn: Name for the new directory (can include path)."""

    if not os.path.exists(dirn):
        os.makedirs(dirn)
    else:
        try:
            os.removedirs(dirn)
        except(OSError):
            print "WARNING: erasing (not empty) directory! ",dirn
            shutil.rmtree(dirn)

        os.makedirs(dirn)

def center_of_mass(xyz):
    """Calculate center of mass from xyz coordinates (xyz file format)."""
    totM = 0.0
    COM = np.zeros((3))
    for line in xyz:
        symbol = line[0]
        coord = np.array(line[1:])

        mass = element(symbol).mass
        totM += mass
        COM += coord*mass

    COM = COM/totM
    return COM

def change_vector_norm(fix, mob, R):
    """Given a fix position: fix(x,y,z), 
    and a mobile position: mob(x,y,z),
    and a new vector magnitude: R,
    we output new coordinates for the mobile position
    such that the distance fix-mob is R without changing the direction."""
    
    # convert to np arrays
    unit = mob - fix
    unit = unit/np.linalg.norm(unit)

    # return new position
    return fix + R * unit


def get_rmsd(coord1, coord2):
    """Calculate RMSD between two sets of coordinates."""

    rmsd = 0
    for c1, c2 in zip(coord1, coord2):

        d1 = np.array([c1[x] for x in range(1,4)])
        d2 = np.array([c2[x] for x in range(1,4)])
        vector = d2 - d1
        rmsd += np.dot(vector, vector)

    rmsd = rmsd/(len(coord1))
    return np.sqrt(rmsd)


def get_distance(table, atoms):
    """Calculate distance between two atoms in table of xyz type."""
    coord1 = np.array([table[atoms[0]][x] for x in range(1,4)])
    coord2 = np.array([table[atoms[1]][x] for x in range(1,4)])
    vector = coord2 - coord1

    return np.linalg.norm(vector)

def get_angle(table, atoms):
    """Calculate angle between three atoms in table of xyz type."""
    coord1 = np.array([table[atoms[0]][x] for x in range(1,4)])
    coord2 = np.array([table[atoms[1]][x] for x in range(1,4)])
    coord3 = np.array([table[atoms[2]][x] for x in range(1,4)])
    vec1 = coord1-coord2
    vec2 = coord3-coord2

    return np.degrees( np.arccos( np.dot(vec1,vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)) ) )

def get_dihedral_angle(table, atoms):
    """Calculate angle between three atoms in table of xyz type.

    Praxeolitic formula (1 sqrt, 1 cross product)."""
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


