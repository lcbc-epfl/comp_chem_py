#!/usr/bin/env python
"""
Set of routine to deal with xyz data of QMMM calculations.

In particular for solute/solvent system it is possible to extract 
a QM region based on a distance criterion from the solvent.
"""

import numpy as np
from comp_chem_utils.utils import center_of_mass

def read_xyz_structure(xyz, nat_solute, nat_solvent):
    """decompose xyz data as solute + group of solvent molecules"""
    xyz_struc = []
    xyz_struc.append(xyz[0:nat_solute])

    nlines = len(xyz)
    xyz_struc += [ xyz[x:x+nat_solvent] for x in xrange(nat_solute, nlines, nat_solvent) ]

    return xyz_struc


def get_distance(xyz_struc, ref):
    """For each block in xyz_struc, calculate distance to ref"""

    distance = [0.0] # solute has distance zero by definition
    for block in xyz_struc[1:]:
        # calculate distance based on Oxygen atom in water molecule
        # i.e. first atom
        coord = np.array(block[0][1:])
        distance.append(np.linalg.norm(ref-coord))

    return distance


def get_QM_region(xyz, nat_slu, nat_sol, nmol_sol):
    """Extract QM information from the full set of coordinates.
    

    The following structure is assumed:

    First ``nat_slu`` atoms are the atoms of the solute molecules 
    The rest are the solvent molecules arranged in blocks. 
    For example, for a water solvent (``nat_sol = 3``)::

        O ...
        H ...
        H ...
        O ...
        H ...
        H ...
        :

    Args:
        xyz (list): xyz_data (list of list) for the whole system.
            Each sublist contains an atomic symbol and 3 (x,y,z) 
            coordinates.

        nat_slu (int): number of atoms of the solute molecule.

        nat_sol (int): number of atoms of a single solvent molecule.

        nmol_sol (int): total number of solvent molecules that
            should be included in the QM region.

    Returns:
        xyz_QM, atom_types, max_dist

        xyz_QM (list): xyz_data for the QM region. 

        atom_types (dic): Each key is an atomic symbol and the 
            corresponding value is a list of indices of the atom 
            of that type present in the QM region. The indices 
            refer to the ordering of the original xyz_data.

        max_dist (float): Distance (in AA) between the center of
            mass of the solute and the center of mass of the farthest
            solvent molecule included in the QM region
    """

    # the distance is measured with respect to the COM of the solute
    ref = center_of_mass(xyz[0:nat_slu])

    # convert xyz data in block structure (solute and solvent)
    xyz_struc = read_xyz_structure(xyz, nat_slu, nat_sol)

    # associate a distance to each block
    dist = get_distance(xyz_struc, ref)

    # get the distance for the farthest molecule to be included in the QM region
    max_dist = sorted(dist)[nmol_sol]

    # Build a collection of atom types.
    # each atom type is characterized by its symbol (e.g. Cl)
    # for each type a list of indices is kept
    idx = 1
    atom_types = {}
    xyz_QM = []
    for d, b in zip(dist, xyz_struc):
        # only include atoms "close" to the solute
        if d<=max_dist:
            for atom in b:
                xyz_QM.append(atom)
                symb = atom[0]
                if symb in atom_types:
                    # update the list of indices
                    atom_types[symb].append( idx )
                else:
                    # start new atom type list
                    atom_types[symb] = [ idx ]
                idx +=1
        else:
            idx += len(b)

    return xyz_QM, atom_types, max_dist




