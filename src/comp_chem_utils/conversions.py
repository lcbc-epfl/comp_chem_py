#!/usr/bin/env python
"""Utilities to convert quantities from and to several units"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import numpy as np

import comp_chem_utils.physcon as const


# electron-volt: by definition, it is the amount of energy gained (or lost) by the charge of a 
# single electron moving across an electric potential difference of one volt.
# 1 volt (1 joule per coulomb, 1 J/C) multiplied by the elementary charge:
EV_TO_JOULES = const.value('charge-e')

# Hartrees: from the definition we have:
AU_TO_JOULES = 2.0*const.value('rydberg')*const.value('planck')*const.value('lightvel')

# Time conversion
AU_TO_S = const.value('dirac') / AU_TO_JOULES
AU_TO_FS = 1.0e15 * AU_TO_S
AU_TO_PS = 1.0e12 * AU_TO_S

# Space conversion
BOHR_TO_ANG = 1.0e10 * const.value('bohrradius')

# ------------------------------------------------------------------------
# SPECIFIC TO THE CALCULATION OF ELECTRONIC SPECTRA
#
# In order to arrive at the extinction or absorption coefficient
# in the conventionally used units [M-1 . cm-1],
# we follow the following procedure:
#
# The spectral function S(w) is the sum over states of products 
# of a unitless oscillator strength and a line shape function.
# The line shape function is expressed in inverse frequency units (seconds).
# Or more generally in the reciprocal units of the excitation energies.
#
# The absorption cross section sigma(w) can be expressed in [Angstroms^2] as
# sigma(w) = SPEC_TO_SIGMA * S(w)
SPEC_TO_SIGMA = 1.0e20 * np.pi * const.value('charge-e') * const.value('charge-e') / ( 2.0 * const.value('mass-e') * const.value('lightvel') * const.value('elec-const'))
# assuming S(w) is expressed in reciprocal 'ANG. FREQ: s-1', i.e. seconds!
# if S(w) is expressed in reciprocal "another_unit" it can be converted to reciprocal 'ANG. FREQ: s-1' as
# S(w) = S(w) / convert("another_unit", 'ANG. FREQ: s-1')
#
# The extinction coefficient is then obtained in [M-1 . cm-1] as
# eps(w) = SIGMA_TO_EPS * sigma(w)
SIGMA_TO_EPS = 1.0e-16 * 1.0e-3 * const.value('avogadro') / np.log(10.0)
# assuming sigma(w) is expressed in 'Angstroms^2'
#
# SPEC_TO_SIGMA = [Angrstom^2 * s-1]
# SIGMA_TO_EPS = [mol-1] * 1.0e-16 * 1.0e-3
# SPEC_TO_SIGMA * SIGMA_TO_EPS = [Angrstom^2 * s-1] * [mol-1] * 1.0e-16 * 1.0e-3
#                              = 1.0e-3 [cm^2 * s-1] * [mol-1 * L] * [L-1]
#                              = [cm^2 * s-1] * [M-1] * [cm-3]
#                              = [cm-1 * s-1] * [M-1]
#                              = [M-1 * cm-1] * [s-1]
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# allowed units for excitation energies or frequencies
convert_to_joules = {
        'ENERGY: J': 1.0,
        'ENERGY: eV': EV_TO_JOULES,
        'ENERGY: a.u.': AU_TO_JOULES,
        'FREQ: s-1': const.value('planck'),
        'ANG. FREQ: s-1': const.value('dirac'), # hbar
        'WAVE NUMBER: cm-1': 1.0e2*const.value('lightvel')*const.value('planck'),
        'WAVE LENGTH: nm': 1.0e9*const.value('lightvel')*const.value('planck') # CAREFUL WITH THAT ONE!!!
        }

def convert(X, from_u, to_u):
    """ returns X converted from 'from_u' unit to 'to_u' unit"""
    # test input
    if from_u not in convert_to_joules:
        print("unit {} not part of convert_to_joules dictionary".format(from_u))
        sys.exit("Wrong unit in energy conversion function")
    if to_u not in convert_to_joules:
        print("unit {} not part of convert_to_joules dictionary".format(to_u))
        sys.exit("Wrong unit in energy conversion function")

    # Wavelength is inversly proportional to the eneergy
    # and has to be treated with special care!

    # CONVERSION TO JOULES
    if from_u == 'WAVE LENGTH: nm':
        # inversely proportional conversion from_u --> Joules
        X_Joules = convert_to_joules[from_u] / X
    else:
        # proportional conversion from_u --> Joules
        X_Joules = X * convert_to_joules[from_u]

    # CONVERSION TO TARGET UNIT
    if to_u == 'WAVE LENGTH: nm':
        # inversely proportional conversion Joules --> to_u
        newX = convert_to_joules[to_u] / X_Joules
    else:
        # proportional conversion Joules --> to_u
        newX = X_Joules / convert_to_joules[to_u]

    return newX

# Conversion tests
def test_conversion(value):
    print("\nTesting all conversions with reference value = {}\n".format(value))
    for from_u in convert_to_joules:
        for to_u in convert_to_joules:
            print(' {:15.8g} {:20} = {:15.8g} {:20}'.format(value, from_u, convert(value, from_u, to_u), to_u))

def print_constants():
    print('\nConversion constants:\n')
    print('EV_TO_JOULES  = {}'.format(EV_TO_JOULES))
    print('AU_TO_JOULES  = {}'.format(AU_TO_JOULES))
    print('SPEC_TO_SIGMA = {}'.format(SPEC_TO_SIGMA))
    print('SIGMA_TO_EPS  = {}'.format(SIGMA_TO_EPS))

if __name__=="__main__":

    print_constants()

    test_conversion(1.0)


