#!/usr/bin/env python
"""Reduce the size of the solvent around a solute.

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

from comp_chem_utils.qmmm_utils import get_QM_region
from comp_chem_utils.cpmd_utils import read_TRAJEC_xyz, write_TRAJEC_xyz


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
    
    parser = argparse.ArgumentParser(description=__doc__)
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

        new_xyz, atom_types, max_dist = get_QM_region(xyz, args.NATQM, args.NATMM, args.NMOLMM)

        if max_dist > super_max_dist:
            super_max_dist = max_dist
    
        new_traj.append( new_xyz )
        
    # write new truncated trajectory file
    write_TRAJEC_xyz(steps, new_traj, output=args.xyz_out)

    print "Maximum distance from solute center of mass: {} AA\n".format(max_dist)
    print "New trajectory printed to {} file".format(args.xyz_out)


