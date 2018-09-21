#!/usr/bin/env python3

from comp_chem_utils.spectrum import read_spectrum, plot_spectrum

output = 'newtonX_iodine_g09.out'
    
# get excitation energies [eV] and oscillator strength
exc, osc = read_spectrum( output, 'newton-x', verbose=True)

plot_spectrum(exc, osc, unit_in='ENERGY: eV',
        nconf=100, fwhm=0.1, temp=0.0, refraction=1.0, 
        ctype='lorentzian', kind='CROSS_SECTION', plot=True)



