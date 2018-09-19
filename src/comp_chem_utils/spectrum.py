#!/usr/bin/env python
"""Collection of methods to read calculate and plot electronic spectra"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import sys
import matplotlib.pyplot as plt
import numpy as np

from comp_chem_utils.utils import get_file_as_list
import comp_chem_utils.physcon as const
import comp_chem_utils.conversions as conv


class Search(object):
    """
    Class for defining how to search spectrum information
    in an output file. Typically from a computational chemistry package.
    """

    def __init__(self, kind, search_str, col_id, offset=0, trim=(0,None), cfac=1.0):
        self.kind = kind
        self.search_str = search_str
        self.col_id = col_id
        self.offset = offset
        self.trim = trim
        self.cfac = cfac

    def get(self, lines, idx):

        try:
            # split data line as a list
            ml = lines[idx + self.offset ].split()
             
            # get element of interest in the line
            data = ml[ self.col_id ][self.trim[0]:self.trim[1]]
             
            # return value converted to proper units
            return float( data ) * self.cfac

        except Exception as e:
            print('ERROR: class info:')
            print('     kind       :{}'.format(self.kind))
            print('     search_str :{}'.format(self.search_str))
            print('     col_id     :{}'.format(self.col_id))
            print('     offset     :{}'.format(self.offset))
            print('     trim       :{}'.format(self.trim))
            print('     cfac       :{}'.format(self.cfac))
            print('ERROR: output line:{}'.format(lines[idx + self.offset ]))
            print('ERROR: {}'.format(e))
            sys.exit('ERROR(read_spectrum): failed to read data from output')


def search_exc(kind):
    """Get excitation energies from an output file"""

    if kind == 'gaussian':
        return Search(kind, ' Excited State ', 4)

    if kind == 'cpmd':
        return Search(kind, 'TD_OS_BERRY|            dE=', -4)

    if kind == 'newton-x':
        return Search(kind, ' Vertical excitation (eV):', 3)

    if kind == 'lsdalton':
        AU_TO_EV = conv.convert(1.0, 'ENERGY: a.u.', 'ENERGY: eV')
        return Search(kind, 'LOFEX: excitation energy   =', -1, cfac=AU_TO_EV)

    if kind == 'turbo_cc2':
        return Search(kind, '  |    frequency : ', 5)

    if kind == 'turbo_tddft':
        return Search(kind, ' Excitation energy / eV:', -1)

def search_osc(kind):
    """Get oscillator strengths from an output file"""

    if kind == 'gaussian':
        return Search(kind, ' Excited State ', 8, trim=(2,None) )

    if kind == 'cpmd':
        return Search(kind, 'TD_OS_BERRY|            dE=', -1)

    if kind == 'newton-x':
        return Search(kind, ' Oscillator strength:', -1)

    if kind == 'lsdalton':
        return Search(kind, 'LOFEX: oscillator strength =', -1)

    if kind == 'turbo_cc2':
        return Search(kind, 'oscillator strength (length gauge)   :', -1)

    if kind == 'turbo_tddft':
        return Search(kind, ' Oscillator strength:', -1, offset=4)


def print_spectrum(exc_ener,strength,output):
    """Print table summary of Excitation energies and oscillator strengths."""

    print('\n   Energy (eV)   Strength  ')
    print(  '  ------------------------ ')
    for e, f in zip(exc_ener, strength):
        print('    {:7.3f}    {:10.5f}   '.format(e,f))
    print('\nSpectrum information succesfully read from {} \n'.format(output))


def read_spectrum(output, kind, verbose=False):
    """Read file output and parse it to find excitation energies and oscillator strengths"""
    out = get_file_as_list(output, raw=True)
    
    exc_ener = []
    strength = []
    # to read one number we need
    #   output type
    #   search string
    #   line offset
    #   column index
    s_exc = search_exc(kind)
    s_osc = search_osc(kind)

    for idx, line in enumerate(out):
    
        if s_exc.search_str in line:
            exc_ener.append( s_exc.get(out, idx) )

        if s_osc.search_str in line:
            strength.append( s_osc.get(out, idx) )

    if verbose:
        print_spectrum(exc_ener,strength,output)

    return np.asarray(exc_ener), np.asarray(strength)


def read_table_spectrum(output, search_str, skip=0, pos_e=1, pos_f=2, verbose=False):
    """Read spectrum data arranged as a table in the output file"""

    out = get_file_as_list(output, raw=True)
    
    exc_ener = []
    strength = []
    for idx, line in enumerate(out):
        if search_str in line:
            keep_reading = True
            i=0
            while(keep_reading):
                try:
                    nl = out[i+skip+idx].split()
                    exc_ener.append( float(nl[pos_e]) )
                    strength.append( float(nl[pos_f]) )
                    i+=1
                except:
                    keep_reading = False

    if verbose:
        print_spectrum(exc_ener,strength,output)

    return np.asarray(exc_ener), np.asarray(strength)


def spectral_function(
        exc, osc, # Excitation frequencies/energies and oscillator strengths
        nconf=1.0, # Number of structures used for conformational sampling
        unit_in='ENERGY: eV', # Unit of the excitations in input
        fwhm=None, # Full Width at Half Maximum for the convolution function, in [unit_in]
        ctype='lorentzian', # Convolution function | options [lorentzian, gaussian]
        x_range=None, # Excitation range, in [unit_in]
        x_reso=None # Number of grid points per unit_in 
        ):
    """
    Calculate the spectral function from theoretical data, 
    excitation energies (exc) and oscillator strength (osc)

    S(w) = 1/nconf sum_{R=1,nconf} sum_{i=1,nstates} osc[i;R] g(w - exc[i;R], fwhm)

    g(w) is the spectral line shape function, (lorentzian or gaussian), 
    expressed in reciprocal angular frequency units, i.e. seconds (per molecules or structure)
    """

    # set default value for fwhm to 0.1 eV
    if not fwhm:
        fwhm = conv.convert(0.1, 'ENERGY: eV', unit_in)

    # set default value for resolution to 100 pts per eV
    if not x_reso:
        x_reso = 100.0 / conv.convert(1.0, 'ENERGY: eV', unit_in)

    # copy input excitation energies to 'ANG. FREQ: s-1' units
    ang_freq = conv.convert(exc, unit_in, 'ANG. FREQ: s-1')
    fwhm_freq = conv.convert(fwhm, unit_in, 'ANG. FREQ: s-1')

    # get x-axis range
    if x_range:
        xmin= x_range[0] 
        xmax= x_range[1]
    else:
        xmax = max(exc) + 4.0*fwhm
        xmin = min(exc) - 4.0*fwhm

    # set number of points from range and resolution
    npts = int((xmax - xmin)*x_reso)

    # set x-values 
    xpts = np.linspace(xmin, xmax, npts)
    ypts = np.zeros(npts)

    # get function parameters from FWHM:
    if ctype=='lorentzian':
        delta = fwhm_freq
        norm = delta/(2.0 * np.pi * nconf)
    elif ctype=='gaussian':
        delta = fwhm_freq/np.sqrt(2.0 * np.log(2.0))
        norm = np.sqrt( 2.0/(delta * delta * np.pi) )/nconf

    print("npts = {}".format(npts))
    print("xmax = {} [{}]".format(xmax, unit_in))
    print("xmin = {} [{}]".format(xmin, unit_in))
    print("Convolution with a {} function and FWHM = {} [{}]".format(ctype, fwhm, unit_in))
    print("Number of structures/conformations = {}".format( nconf))

    # make convolution of spectrum
    for i, f in enumerate(osc):
        tmp = (ang_freq[i] - conv.convert(xpts, unit_in, 'ANG. FREQ: s-1')) 
        if ctype=='lorentzian':
            ypts += f*norm/( tmp*tmp + (delta/2.0)**2.0 )

        elif ctype=='gaussian':
            ypts += f*norm*np.exp( - 2.0 * (tmp/delta)**2.0 )

        else:
            sys.exit('Unknown convolution type: {}'.format(ctype))

    # return spectral function as data arrays in [unit_in] unit
    return xpts, ypts


def temperature_effect(E, unit_in='ENERGY: eV', temp=0.0):
    """ Calculate tenperature factor for the absorption cross section"""
    if temp < 1.0e-4:
        factor = 0.0
    else:
        # temperature is in Kelvin
        # Boltzmann constant is in J.K-1
        # So E must be in joules
        E_Joules = conv.convert(E, unit_in, 'ENERGY: J')
        factor = np.exp( - E_Joules / (const.value('boltzmann') * temp) )
        
    return 1.0 - factor


def cross_section(
        freq, # Grid data of the x-axis (energy/frequencies) corresponding to spectral_func
        spectral_func, # Grid data of the spectral function in reciprocal angular frequency units (seconds)
        unit_in='ENERGY: eV', # Input energy units for the freq array
        temp=0.0, # Temperature in Kelvin
        refraction=1.0 # Refraction index
        ):
    """
    Calculate the absorption cross-section from the spectral function.
    The cross section is just a scaled version of the spectral function
    that might depend on temperature and the refraction index of the medium.

    The spectral function is calculated in [Angrstom^2 * molecule-1]
    """
    
    sigma = []
    for E, S in zip(freq, spectral_func):
        # SPEC_TO_SIGMA is the main conversion constant defined in the conversions module
        sigma.append( S * temperature_effect(E, unit_in=unit_in, temp=temp) * conv.SPEC_TO_SIGMA / refraction )

    return np.asarray(sigma)


def plot_stick_spectrum(exc, osc, color=None, label=None):

    for e, f in zip(exc, osc):
        myh = plt.vlines(e, -0.1, f, color=color, label=label)

    return myh


spectra_kinds={
        'STICKS': 'Oscillator strength [Arbitrary units]', 
        'SPECTRAL FUNC': 'Spectral function [s. $\\cdot $ molecules$^{-1}$]', 
        'CROSS_SECTION': 'Cross section [$\\AA^2 \\cdot $ molecules$^{-1}$]', 
        'EXPERIMENTAL': 'Molar absorptivity [M$^{-1}  \\cdot $cm$^{-1}$]'
        }

def plot_spectrum(exc, osc, unit_in='ENERGY: eV',
        nconf=1, fwhm=0.1, temp=0.0, refraction=1.0, 
        ctype='lorentzian', kind='CROSS_SECTION', plot=True):
    """
    Plot a spectrum based on theoretical data points.

    Excitation energies and oscillator strengths.
    Several kinds of representation can be chosen.
    The x and y points are returned
    """

    if kind=='STICKS':
        plot_stick_spectrum(exc, osc, color='b')
        if not plot:
            print('That was quite stupid...')
            return exc, osc

    else:
        # get spectral function in seconds per molecule
        xpts, ypts = spectral_function(exc, osc, nconf=nconf, unit_in=unit_in, fwhm=fwhm, ctype=ctype)

        if kind!='SPECTRAL FUNC':
            # get abs. cross section in Angstrom^2 per molecule
            ypts = cross_section(xpts, ypts, unit_in=unit_in, temp=temp, refraction=refraction)

            if kind!='CROSS_SECTION':
                # get Molar absorptivity in M^-1 cm-1
                ypts = ypts * conv.SIGMA_TO_EPS

        if plot:
            # plot grid points
            plt.plot(xpts, ypts)

    if plot:
        # horizontal black line at 0
        plt.axhline(0, color='k')

        # axis labels
        plt.xlabel(unit_in)
        plt.ylabel(spectra_kinds[kind])

        # show plot
        plt.show()

    # return grid points
    return xpts, ypts



