#!/usr/bin/env python
"""
Simple script to plot absorption spectra from CPMD TDDFT output files.

Usage::

    python cpmd_spec.py cpmd1.out cpmd2.out ...

Where cpmdx.out denote a CPMD output file for a TDDFT calculation.
At least one output file should be given.
"""

import sys
import seaborn as sns
import matplotlib.pyplot as plt

import comp_chem_utils.spectrum as sp


if __name__=="__main__":
    # TRY TO READ LIST OF OUTPUT FILES
    try:
        flist = sys.argv[1:]
        if not flist:
            raise Exception('We need some parameters...')
    except:
        print(__doc__)
        sys.exit()


    # READ THE OUTPUT FILES AND EXTRACT SPECTRUM INFORMATION
    exc = []
    osc = []
    for filename in flist:
        e, f = sp.read_spectrum(filename,"cpmd")
        exc.append(e)
        osc.append(f)

    # FOR EACH SET OF DATA:
    #   - PLOT THE ABS. CROSS SECTION
    #   - PLOT A STICK SPECTRUM
    #lab=['2.62', '2.66', '2.70','2.74','2.78']
    #lab=['Vacuum 2.62 AA','Vacuum 2.70 AA','Original QM, 3 eth.', 'New QM, closest eth.']
    #lab=['Vacuum 2.70 AA','snapshot 1', 'snapshot 2', 'snapshot 3']
    for i,(e,f) in enumerate(zip(exc, osc)):
        x, y = sp.plot_spectrum(e,f,plot=False)
        tmp = plt.plot(x, y, label=flist[i])
        #tmp = plt.plot(x, y, label=lab[i])
        plt.vlines(e, -0.1, f, color=tmp[0].get_color())
    
    # ADD BLACK LINE AT f=0
    plt.axhline(color='k')

    # SET UP THE PLOT AND SHOW IT
    plt.ylabel('Cross section [$\\AA^2 \\cdot $ molecules$^{-1}$]')
    plt.xlabel('Transition energy [eV]')
    plt.legend()
    plt.show()


