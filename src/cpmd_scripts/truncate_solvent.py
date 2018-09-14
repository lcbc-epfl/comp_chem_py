#!/usr/bin/env python
description="""
Reduce the size of the solvent around a solute.

Take CPMD TRAJEC.xyz file as input.
The molecular structure is assumed to be composed of a solute part and a solvent part.
The purpose of this script is to generate a new trajectory file
with a reduced number of solvent molecules around the solute.
"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import numpy as np
import os
import sys
import argparse

from comp_chem_utils.utils import center_of_mass
from comp_chem_utils.cpmd_utils import read_TRAJEC_xyz, write_TRAJEC_xyz


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

# ---------------------------------------------------------------------
# for compatibility with autodoc in sphinx
if __name__ == '__main__':

    # Tell the users about the assumptions on the xyz structure
    print("""
        The following structure is assumed:
            First X atoms are the atoms of the solute
            The rest are the solvent molecules (e.g. water molecules)
            arranged in blocks. For example:
          O ...
          H ...
          H ...
          O ...
          H ...
          H ...
          :
          """)
    
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("NATQM", type=int, help="Number of atoms in QM region or in solute")
    parser.add_argument("--xyz_inp", type=str, help="Input name for trajectory file in xyz format", default='TRAJEC.xyz')
    parser.add_argument("--xyz_out", type=str, help="Output name for trajectory file in xyz format")
    parser.add_argument("--NATMM", type=int, help="Number of atoms per solvent molecules (e.g. 3 for water)", default=3)
    parser.add_argument("--NMOLMM", type=int, help="Number of solvent molecules to keep in the output structure", default=0)
    args = parser.parse_args()
    if not args.xyz_out:
        args.xyz_out = '{}_{}.xyz'.format(args.xyz_inp[:-4], args.NMOLMM)
    
    # print info
    print 'Input trajectory read from  : {}'.format(args.xyz_inp)
    print 'Output trajectory printed to: {}\n'.format(args.xyz_out)
    print 'Number of atoms in QM region or in solute: {}'.format(args.NATQM)
    print 'Number of atoms per solvent molecules    : {}'.format(args.NATMM)
    print 'Number of solvent molecules to keep      : {}\n'.format(args.NMOLMM)
    
    # read trajectory file
    steps, traj_xyz = read_TRAJEC_xyz(args.xyz_inp)
    
    print 'Number of structures read from input     : {}\n'.format(len(traj_xyz))
    
    # get truncated xyz data for each xyz geometry in trajectory
    new_traj = []
    super_max_dist = 0.0
    for xyz in traj_xyz:
        # reference point is taken from COM of solute
        ref = center_of_mass(xyz[0:args.NATQM]) 
    
        # convert xyz data in block structure (solute and solvent)
        xyz_struc = read_xyz_structure(xyz, args.NATQM, args.NATMM)
    
        # associate a distance to each block
        dist = get_distance(xyz_struc, ref)
    
        # get the distance for the farthest molecule to be included
        max_dist = sorted(dist)[args.NMOLMM] 
        if max_dist > super_max_dist:
            super_max_dist = max_dist
    
        # remove all blocks with distance larger than max_dist
        new_struc = []
        for d, b in zip(dist, xyz_struc):
            if d<=max_dist:
                new_struc.append(b)
    
        # remove block structure
        new_traj.append( [atom for block in new_struc for atom in block] )
        
    
    # write new truncated trajectory file
    write_TRAJEC_xyz(steps, new_traj, output=args.xyz_out)
    
    print "Maximum distance from solute center of mass: {} AA\n".format(max_dist)
    print "New trajectory printed to {} file".format(args.xyz_out)


