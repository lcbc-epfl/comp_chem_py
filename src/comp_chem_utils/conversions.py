#!/usr/bin/env python
"""Functions and constants to convert quantities from and to different units.

When executed directly, calls the two test functions:
    
* :py:func:`~comp_chem_utils.print_constants`
* :py:func:`~comp_chem_utils.test_conversion`
"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import numpy as np

import comp_chem_utils.physcon as const


# Energy conversions
EV_TO_JOULES = const.value('charge-e')
"""Energy conversion from eV to J.

By definition an electron-volt is the amount of energy gained (or lost) 
by the charge of a single electron moving across an electric potential difference 
of one volt. So 1 eV is 1 volt (1 joule per coulomb, 1 J/C) multiplied by the 
elementary charge.
"""

AU_TO_JOULES = 2.0*const.value('rydberg')*const.value('planck')*const.value('lightvel')
"""Energy conversion from Hartree (a.u.) to J.

1 Hartree (a.u.) is defined as

.. math::
    E_{h} = 2 R_\\infty h c
"""


# Time conversion
AU_TO_S = const.value('dirac') / AU_TO_JOULES
"""Time conversions from a.u. to seconds."""

AU_TO_FS = 1.0e15 * AU_TO_S
"""Time conversions from a.u. to femto-seconds."""

AU_TO_PS = 1.0e12 * AU_TO_S
"""Time conversions from a.u. to pico-seconds."""

# Space conversion
BOHR_TO_ANG = 1.0e10 * const.value('bohrradius')


# ------------------------------------------------------------------------
# SPECIFIC TO THE CALCULATION OF ELECTRONIC SPECTRA
#
# TODO: move this documentation to the spectrum documentation
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
"""Conversion constant  between a spectral function ``S`` expressed in seconds (reciprocal 
angular frequency unit) and the absorption cross section ``\sigma`` expressed in ``Angstroms^2``.

.. math::
    \sigma(\omega) = \\text{SPEC\_TO\_SIGMA} \cdot S(\omega)

.. math::
    \\text{SPEC\_TO\_SIGMA} = 10^{20} \\frac{\\pi e^{2} }{2 m_{e} c \\epsilon_{0}} \cdot S(\\omega)

For more details see the documentation of the :py:mod:`~comp_chem_utils.spectrum` module.
"""
# if S(w) is expressed in reciprocal "another_unit" it can be converted to reciprocal 'ANG. FREQ: s-1' as
# S(w) = S(w) / convert("another_unit", 'ANG. FREQ: s-1')
#
# The extinction coefficient is then obtained in [M-1 . cm-1] as
# eps(w) = SIGMA_TO_EPS * sigma(w)
SIGMA_TO_EPS = 1.0e-16 * 1.0e-3 * const.value('avogadro') / np.log(10.0)
"""Conversion constant between an absorption cross section ``\sigma`` expressed in ``Angstroms^2`` 
and the extinction coefficient expressed in ``M^{-1} . cm^{-1}``.

.. math::
    \epsilon(\omega) = \\text{SIGMA\_TO\_EPS} \cdot \\sigma(\omega)

.. math::
    \\text{SIGMA\_TO\_EPS} = 10^{-16} \\frac{N_{A}}{10^{3} \\ln{10}} \cdot S(\\omega)

For more details see the documentation of the :py:mod:`~comp_chem_utils.spectrum` module.
"""
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
    """Convert the energy value ``X`` from the unit ``from_u`` to the unit ``to_u``.
    
    Args:
        X (float): Energy value to convert.
        from_u (str): Unit of the input energy value ``X`` given as one 
            of the keys of the dictionary ``convert_to_joules``.
        to_u (str): Unit of the output energy value ``X`` given as one 
            of the keys of the dictionary ``convert_to_joules``.
            
    Returns:
        The ``X`` value is returned expressed in the ``to_u`` unit.

    Example:
        >>> from comp_chem_utils import conversions as c
        >>> c.convert(1.0, 'WAVE LENGTH: nm', 'ENERGY: eV')
        1239.8419292004205
    """

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
def test_conversion(value=1.0):
    """Use all possible conversions for a single value and printout the results.

    This is intended as a test.

    Args:
        value (float, optional): hypothetical energy value. Default is 1.0.

    """

    print("\nTesting all conversions with reference value = {}\n".format(value))
    for from_u in convert_to_joules:
        for to_u in convert_to_joules:
            print(' {:15.8g} {:20} = {:15.8g} {:20}'.format(value, from_u, convert(value, from_u, to_u), to_u))

def print_constants():
    """Print all conversion constants defined in this module."""
    print('\nConversion constants:\n')
    print('EV_TO_JOULES  = {}'.format(EV_TO_JOULES))
    print('AU_TO_JOULES  = {}'.format(AU_TO_JOULES))
    print('AU_TO_S       = {}'.format(AU_TO_S))
    print('AU_TO_FS      = {}'.format(AU_TO_FS))
    print('AU_TO_PS      = {}'.format(AU_TO_PS))
    print('SPEC_TO_SIGMA = {}'.format(SPEC_TO_SIGMA))
    print('SIGMA_TO_EPS  = {}'.format(SIGMA_TO_EPS))


if __name__=="__main__":

    print_constants()
    test_conversion()


