#!/usr/bin/env python
"""Converts a CPMD TRAJECTORY file into a TRAJEC.xyz file."""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import argparse

from comp_chem_utils.conversions import BOHR_TO_ANG
from comp_chem_utils.cpmd_utils import read_TRAJECTORY, read_GEOMETRY_xyz, write_TRAJEC_xyz


# ---------------------------------------------------------------------
# for compatibility with autodoc in sphinx
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xyz", type=str, help="xyz file for atom type structure (e.g. GEOMETRY.xyz)")
    parser.add_argument("traj", type=str, help="CPMD TRAJECTORY file with xyz coordinates.")
    args = parser.parse_args()

    # get list of atom symbols
    geo_xyz = read_GEOMETRY_xyz(args.xyz)
    symb = [line[0] for line in geo_xyz]
    
    # read trajectory file
    stp, xyz = read_TRAJECTORY(args.traj)[:2]
    
    # convert to angstroms and xyz format
    new_traj = []
    for x in xyz:
    
        new_xyz = []
        for idx, line in enumerate(x):
            new_xyz.append( [symb[idx]] + list(line * BOHR_TO_ANG) )
    
        new_traj.append(new_xyz)
    
    
    # write output file to disk
    out = raw_input('Output filename for TRAJEC.xyz file [NEW_TRAJEC.xyz]:\n')
    if out == '':
        out = 'NEW_TRAJEC.xyz'
    
    write_TRAJEC_xyz(stp, new_traj, out)
