#!/usr/bin/env python
"""PDOS.py is a program which automatically extracts information from Gaussian 
(http://www.gaussian.com/) output and calculates partial density of states 
(PDOS) based upon Lowdin or Mulliken orbital analysis."""


__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import argparse

from PDOS.modules.PDOS_module import calculate_and_plot_pdos

# for compatibility with autodoc in sphinx
if __name__ == '__main__':

    # read user arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gaussian_file", type=str, help="Gaussian log (output) file containing relevant info")
    parser.add_argument("group_file", type=str, help="file with the orbital groups for which you want to plot PDOSs")
    parser.add_argument("--npts", type=int, help="number of points to use to plot the PDOS", default=0)
    parser.add_argument("--fwhm", type=float, help="Full width at half maximum [eV]", default=1.0)
    parser.add_argument("--xmin", type=int, help="Minimum value on x axis", default=1)
    parser.add_argument("--xmax", type=int, help="Maximum value on x axis", default=-1)
    parser.add_argument("--plot_dos", type=bool, help="To plot or not the total DOS", default=True)
    parser.add_argument("--mulliken", type=bool, help="To use mulliken pop analysis instead of Lowdin", default=False)
    args = parser.parse_args()
    
    # plot PDOS
    calculate_and_plot_pdos(args.gaussian_file, args.group_file, 
            npts=args.npts, 
            fwhm=args.fwhm, 
            xmin=args.xmin, 
            xmax=args.xmax, 
            plot_dos=args.plot_dos, 
            mulliken=args.mulliken)


