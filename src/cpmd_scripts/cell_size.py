#!/usr/bin/env python
description="""
Determine the size of the CELL for a CPMD calculation 

Two options are available:
    reading from xyz format
    reading from gromos format

In the case of gromos format, a list of residues 
in the QM region has to be specified.

xyz data can be in Aangstroms or in Bohrs
"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import sys
import argparse
import numpy as np

from comp_chem_utils.utils import get_file_as_list
from comp_chem_utils.conversions import BOHR_TO_ANG
from comp_chem_utils.molecule_data import read_xyz_table

# for compatibility with autodoc in sphinx
if __name__ == '__main__':

   parser = argparse.ArgumentParser(description=description)
   parser.add_argument("inpfile", type=str, help="input file with xyz coordinates (.xyz file or gromos format)")
   
   # READ INPUT FILE
   inp = get_file_as_list(parser.parse_args().inpfile)
   
   
   # GET INFO FROM USER
   print "Format of coordinates:"
   print "   0: xyz format"
   print "   1: gromos format"
   frmt = input("")
   
   print "Units of coordinates:"
   print "   0: Aangstroms (AA)"
   print "   1: Bohr radius (a.u.)"
   unit = input("")
   
   
   # READ COODRINATES
   if frmt==0:
       # read xyz file
       table = read_xyz_table(inp[2:])
       xyz = []
       for line in table:
           xyz.append( [ float(x) for x in line[1:4] ])
   
   elif frmt==1:
   
       # get list of residue names from user
       rn_list = raw_input("List of residues in the QM region:\n").split()
       
       # read gromos.crd file
       
       # extract coordinates of the atoms in the QM region
       xyz = []
       for line in inp:
           ml = line.split()
           try:
               rn = ml[1]
               if rn in rn_list:
                   # for some reason in gromos.crd there is a factor 10 
                   # inserted in the coordinates.......
                   xyz.append( [10.0*float(x) for x in ml[4:7]] )
       
           except(IndexError):
               # line is not appropriate
               continue
   
   coord = np.array(xyz)
   
   if unit==0:
       # convert data to Bohr
       coord = coord/BOHR_TO_ANG
   
   
   # CALCULATE CELL SIZE AND PRINT
   minxyz = np.amin(coord, axis=0) # get min xyz values
   maxxyz = np.amax(coord, axis=0) # get max xyz values
   
   diff = np.subtract(maxxyz, minxyz)
   print "                     {:^12} {:^12} {:^12} ".format('x','y','z')
   print "Raw max diff (AA)  : {:12.6f} {:12.6f} {:12.6f} ".format(*diff)
   print "Raw max diff (a.u.): {:12.6f} {:12.6f} {:12.6f} ".format(*diff/BOHR_TO_ANG)
   print ""
   print "Cell size    (AA)  : {:12.6f} {:12.6f} {:12.6f} ".format(*diff+7)
   print "Cell size    (a.u.): {:12.6f} {:12.6f} {:12.6f} ".format(*(diff+7)/BOHR_TO_ANG)
   print ""
   diff = (diff+7)
   dim = np.amax(diff)/BOHR_TO_ANG
   print "Recommended (cubic) cell dimmension = {:9.3f} (AA)".format(np.amax(diff))
   print "Recommended (cubic) cell dimmension = {:9.3f} (a.u.)".format(dim)
   print ""
   print "   CELL"
   print "   {:.3f} 1. 1. 0 0 0".format(dim)

