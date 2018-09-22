#!/usr/bin/env python
"""Collection of functions to read, calculate, and plot electronic spectra."""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import sys
import matplotlib.pyplot as plt
import numpy as np

from comp_chem_utils.utils import get_file_as_list
import comp_chem_utils.physcon as const
import comp_chem_utils.conversions as conv


class Search(object):
    """Class for defining how to search spectrum information in an computational chemistry output file.

    The information needed are excitation energies and oscillator strength.
    The assumption here is that in the output file this information can be located
    from a string ``search_str``.

    The data to be read can be shifted from the ``search_str`` by an ``offset`` number
    of lines. It is then found in a specific ``col_id`` inside the line which is splitted
    based on blank characters. It can be trimmed by ``trim`` and converted by ``cfac``.

    Example:
        Possible output file (named ``'my_output'``)::

            Excited State 
            E=1.0000cm-1 f=0.00000

        In that case to extract the energy in eV we would have to do::

            >>> from comp_chem_utils.spectrum import Search
            >>> my_search_obj = Search('my_output', 'Excited State', 0, offset=1, trim=(2,7), cfac=1.2398e-4)
            >>> my_search_obj.get_all()[0]
            0.00012398

    Args:
        fname (str): Name of the output file containing the spectrum information. Can include path to file.

        search_str (str): String that will be search in the output file
            to locate the data to be extracted.

        col_id (int): The line matching the search string is splitted into a list with :py:meth:`str.split()` 
            and the element ``col_id`` of the list is extracted.

        offset (int, optional): In case the line of interest is not the line with ``search_str`` the
            ``offset`` can be used to shift to a line above (negative offset) or bellow the matching line.
            Default value is 0.

        trim (tuple): In case the element ``col_id`` extracted from the line contains more than
            the desired data. A subset of the string can be extracted with trim. See example.
            Default is ``(0,None)``, i.e. it takes the whole element.

        cfac (float, optional): In case the extracted data is not in the desired unit, it can be
            multiplied by ``cfac`` to convert it. Default value is 1.0.
    """

    def __init__(self, fname, search_str, col_id, offset=0, trim=(0,None), cfac=1.0):
        self.fname = fname
        self.search_str = search_str
        self.offset = offset
        self.col_id = col_id
        self.trim = trim
        self.cfac = cfac

    def get_all(self):
        """Extract all the relevant data depending on attribute values.
        
        Return:
            The extracted data is returned as floats in a 1-D ``np.array()``.
        """

        try:
            out = get_file_as_list(self.fname, raw=True)
    
            all_data = []

            for idx, line in enumerate(out):
    
                if self.search_str in line:
                    all_data.append( self.get_single(out, idx) )

        except Exception as e:
            print('ERROR: class info:')
            print('     fname      :{}'.format(self.fname))
            print('     search_str :{}'.format(self.search_str))
            print('     col_id     :{}'.format(self.col_id))
            print('     offset     :{}'.format(self.offset))
            print('     trim       :{}'.format(self.trim))
            print('     cfac       :{}'.format(self.cfac))
            print('ERROR: output line:{}'.format(lines[idx + self.offset ]))
            print('ERROR: {}'.format(e))
            sys.exit('ERROR(read_spectrum): failed to read data from output')

        return np.asarray(all_data)

    def get_single(self, lines, idx):
        """Once a matching line has been located this method returns the associated data."""

        # split data line as a list
        ml = lines[idx + self.offset ].split()
         
        # get element of interest in the line
        data = ml[ self.col_id ][self.trim[0]:self.trim[1]]
         
        # return value converted to proper units
        return float( data ) * self.cfac


def search_exc(kind, output):
    """Get excitation energies from an output file.
    
    The ``kind`` of output file is used to determine the parameters to be used 
    to set the :py:class:`~comp_chem_utils.spectrum.Search` class.
    """

    if kind == 'gaussian':
        mys = Search(output, ' Excited State ', 4)

    if kind == 'cpmd':
        mys = Search(output, 'TD_OS_BERRY|            dE=', -4)

    if kind == 'newton-x':
        mys = Search(output, ' Vertical excitation (eV):', 3)

    if kind == 'lsdalton':
        AU_TO_EV = conv.convert(1.0, 'ENERGY: a.u.', 'ENERGY: eV')
        mys = Search(output, 'LOFEX: excitation energy   =', -1, cfac=AU_TO_EV)

    if kind == 'turbo_cc2':
        mys = Search(output, '  |    frequency : ', 5)

    if kind == 'turbo_tddft':
        mys = Search(output, ' Excitation energy / eV:', -1)

    return mys.get_all()


def search_osc(kind, output):
    """Get oscillator strengts from an output file.
    
    The ``kind`` of output file is used to determine the parameters to be used 
    to set the :py:class:`~comp_chem_utils.spectrum.Search` class.
    """

    if kind == 'gaussian':
        mys = Search(output, ' Excited State ', 8, trim=(2,None) )

    if kind == 'cpmd':
        mys = Search(output, 'TD_OS_BERRY|            dE=', -1)

    if kind == 'newton-x':
        mys = Search(output, ' Oscillator strength:', -1)

    if kind == 'lsdalton':
        mys = Search(output, 'LOFEX: oscillator strength =', -1)

    if kind == 'turbo_cc2':
        mys = Search(output, 'oscillator strength (length gauge)   :', -1)

    if kind == 'turbo_tddft':
        mys = Search(output, ' Oscillator strength:', -1, offset=4)

    return mys.get_all()


def print_spectrum(exc_ener,strength,output):
    """Print table summary of Excitation energies and oscillator strengths."""

    print('\n   Energy (eV)   Strength  ')
    print(  '  ------------------------ ')
    for e, f in zip(exc_ener, strength):
        print('    {:7.3f}    {:10.5f}   '.format(e,f))
    print('\nSpectrum information succesfully read from {} \n'.format(output))


def read_spectrum(output, kind, verbose=False):
    """Read an output file and parse it to find excitation energies and oscillator strengths.

    Note: 
        If the spectrum data is written as a table, i.e. in the following format::

            THIS IS A TABLE OF SPECTRUM DATA:
            E (eV)     f
            1.000    0.000
            1.000    0.000
            :        :

        Then the :py:func:`~comp_chem_utils.spectrum.read_table_spectrum` function
        should be used instead.
    
    Args: 
        output (str): Name of the output file containing the spectrum information. 
            Can include path to file.

        kind (str): Kind of output, e.g. 'gaussian' or 'lsdalton'. This is used in the
            :py:func:`~comp_chem_utils.spectrum.search_exc` and
            :py:func:`~comp_chem_utils.spectrum.search_osc` 
            to set the parameters of the :py:class:`~comp_chem_utils.spectrum.Search`
            class. New ``kind`` have to be implemented in those functions.

        verbose (bool, optional): If ``True`` the extracted data will be printed out.
            Default is ``False``.

    Return:
        exc, osc

        Excitation energies and oscillator strengths are returned in the form of the two
        1-D ``np.array()``, ``exc`` and ``osc``.
    """

    out = get_file_as_list(output, raw=True)
    
    exc = search_exc(kind, output)
    osc = search_osc(kind, output)

    if verbose:
        print_spectrum(exc,osc,output)

    return exc, osc


def read_table_spectrum(output, search_str, offset=0, pos_e=0, pos_f=1, verbose=False):
    """Read spectrum data arranged as a table in the output file.
    
    Args:
        output (str): Name of the output file containing the spectrum information. 
            Can include path to file.

        search_str (str): String that will be search in the output file
            to locate the data to be extracted.

        offset (int, optional): Number of lines to skip after the matching line before
            the table starts. Default is 0, which means that the table is assumed to start
            right after the matching line.

        pos_e (int, optional): Column index for the excitation energies. 
            Default is 0.

        pos_f (int, optional): Column index for the oscillator strengths.
            Default is 1.

        verbose (bool, optional): If ``True`` the extracted data will be printed out.
            Default is ``False``.

    Return:
        exc, osc

        Excitation energies and oscillator strengths are returned in the form of the two
        1-D ``np.array()``, ``exc`` and ``osc``.
    """

    out = get_file_as_list(output, raw=True)
    
    exc_ener = []
    strength = []
    for idx, line in enumerate(out):
        if search_str in line:
            keep_reading = True
            i=0
            while(keep_reading):
                try:
                    nl = out[i+offset+idx].split()
                    exc_ener.append( float(nl[pos_e]) )
                    strength.append( float(nl[pos_f]) )
                    i+=1
                except:
                    keep_reading = False

    if verbose:
        print_spectrum(exc_ener,strength,output)

    return np.asarray(exc_ener), np.asarray(strength)


def spectral_function(exc, osc, unit_in='ENERGY: eV', nconf=1, fwhm=None, ctype='lorentzian',
        x_range=None, x_reso=None):
    """Calculate the spectral function from theoretical data (excitation energies and oscillator strengths).

    Note:
        It is not recommended to use this function directly. Instead the 
        :py:func:`~comp_chem_utils.spectrum.plot_spectrum` function should 
        be used which serves as a wrapper and gives more flexibility on the output data.

    The spectral function is calculated as

    .. math::
        S(\omega) = \\frac{1}{N_\\text{conf}} \sum_{R=1}^{N_\\text{conf}} \sum_{i=1}^{N_\\text{states}}
        f_{i}(R) \cdot g( \omega - \omega_{i}(R), \\delta) 

    This expression is very well descrined in, e.g. [Barbatti2010a]_.
    (f_i, \omega_i) is the pair of input excitation energies and oscillator strengths, while 
    \omega is the incident frequency. The spectral line shape function is g which depends on the 
    Full Width at Half Maximum \delta, it is expressed in reciprocal angular frequency units, 
    i.e. seconds (per molecules or structure).

    Args:
        exc: Input excitation energies given in a one dimenssional ``np.array()``.

        osc: Input oscillator strengths given in a one dimenssional ``np.array()``.

        unit_in (str, optional): String describing the unit used for the input 
            excitation energies. The string must correspond to one of the keys of 
            the dictionary ``convert_to_joules`` in the :py:mod:`~comp_chem_utils.conversions`
            module. Default is ``'ENERGY: eV'``.

        nconf (int, optional): Total number of conformations used in the input data.
            This is used for normalization to a single structure. Default is 1.

        fwhm (float, optional): Full width at Half Maximum used in the convolution function. 
            It must be given in the same units as ``unit_in``. The Default value is ``None``,
            which will be latter changed to correspond to 0.1 eV.

        ctype (str, optional): Defines the type of convolution. Either ``'lorentzian'`` which is default 
            or ``'gaussian'``.

        x_range (list, optional): This is a two-value list defining the range of energy data for which
            the spectral function has to be calculated. It should be given in the same units 
            as ``unit_in``. It is default to ``None`` which will be latter changed to 
            appropriate values related to the FWHM.

        x_reso (int, optional): Resolution of the spectral function given as the number of grid points
            per energy unit (``unit_in``). The default value is ``None``,  which will be latter changed 
            to correspond to 100 pts per eV.

    Returns:
        xpts, ypts

        Those are two ``np.array()`` containing the grid points required to plot the spectral function.
        No matter the unit of the input energies, (``unit_in``). The spectral function (in ypts) is 
        expressed in reciprocal angular frequency (seconds). While the xpts values are expressed in ``unit_in``
        energy unit.
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


def temperature_effect(E, unit_in='ENERGY: eV', temp=None):
    """Calculate tenperature factor for the absorption cross section.

    The temperature effect is calculated as

    .. math::
        1.0 - \exp [ - E / (k_{B} T) ]

    where E is the input energy which will be converted to Joules.

    Args:
        E (float): Input energy or frequency expressed in ``unit_in``.

        unit_in (str, optional): String describing the unit used for the input 
            excitation energies. The string must correspond to one of the keys of 
            the dictionary ``convert_to_joules`` in the :py:mod:`~comp_chem_utils.conversions`
            module. Default is ``'ENERGY: eV'``.

        temp (float, optional): Working temperature in Kelvin. The default value is ``None``, 
            which means that no temperature effect will be calculated.

    Returns:
        The temperature effect as a float. If ``temp=None``, the value returned is 1.0.
    """
    if not temp:
        factor = 0.0
    else:
        # temperature is in Kelvin
        # Boltzmann constant is in J.K-1
        # So E must be in joules
        E_Joules = conv.convert(E, unit_in, 'ENERGY: J')
        factor = np.exp( - E_Joules / (const.value('boltzmann') * temp) )
        
    return 1.0 - factor


def cross_section(xpts, ypts, unit_in='ENERGY: eV', temp=None, refraction=1.0):
    """Calculate the absorption cross-section from the spectral function.

    The cross section is just a scaled version of the spectral function
    that might depend on temperature and the refraction index of the medium.

    Note:
        It is not recommended to use this function directly. Instead the 
        :py:func:`~comp_chem_utils.spectrum.plot_spectrum` function should 
        be used which serves as a wrapper and gives more flexibility on the output data.

        The ``xpts``, and ``ypts`` input data are expected to come directly
        from a call to :py:func:`~comp_chem_utils.spectrum.spectral_function`.

    Args:
        xpts (np.array): Grid points corresponding to the x-axis of the spectral
            function. They are expressed in ``unit_in``.

        ypts (np.array): Grid points corresponding to the y-axis of the spectral
            function. They have to be given in reciprocal angular frequency units 
            (seconds).

        unit_in (str, optional): String describing the unit used for the input 
            excitation energies. The string must correspond to one of the keys of 
            the dictionary ``convert_to_joules`` in the :py:mod:`~comp_chem_utils.conversions`
            module. Default is ``'ENERGY: eV'``.

        temp (float, optional): Working temperature in Kelvin. The default value is ``None``, 
            which means that no temperature effect will be calculated.

        refraction (float, optional): The refraction index of the medium. Default value is 1.0.

    Return:
        The absorption cross section is returned in ``Angstrom^2 . molecule^{-1}``.
        It is calculated from the spectral function as

        .. math::
            \\sigma(\omega) = S (\omega) \cdot \\text{temp\_effect} 
            \cdot \\text{SPEC\_TO\_SIGMA} n.

        where the temp_effect comes from :py:func:`~comp_chem_utils.spectrum.temperature_effect`,
        ``SPEC_TO_SIGMA`` is the main conversion constant from 
        :py:mod:`~comp_chem_utils.spectrum.conversions`, and n is the refraction index.

        The absorption cross section is returned as grid points in an ``np.array()``.
    """
    
    sigma = []
    for E, S in zip(xpts, ypts):
        # SPEC_TO_SIGMA is the main conversion constant defined in the conversions module
        sigma.append( S * temperature_effect(E, unit_in=unit_in, temp=temp) * conv.SPEC_TO_SIGMA / refraction )

    return np.asarray(sigma)


def plot_stick_spectrum(exc, osc, color=None, label=None):
    """Plot a stick spectrum from theoretical data.

    The raw excitation energies and oscillator strength are
    used to make a stick spectrum.

    The color and label arguments should be specified in order
    to ensure a uniform color for all the sticks.

    Args:
        exc: Input excitation energies given in a one dimenssional ``np.array()``.

        osc: Input oscillator strengths given in a one dimenssional ``np.array()``.

        color (optional): color code used in matplotlib.

        label (optional): label to describe the data.

    Return:
        Handle object comming out of the plt call that can be used for 
        the legend.
    """

    for e, f in zip(exc, osc):
        myh = plt.vlines(e, -0.1, f, color=color, label=label)

    return myh


spectra_kinds={
        'STICKS': 'Oscillator strength [Arbitrary units]', 
        'SPECTRAL_FUNC': 'Spectral function [s. $\\cdot $ molecules$^{-1}$]', 
        'CROSS_SECTION': 'Cross section [$\\AA^2 \\cdot $ molecules$^{-1}$]', 
        'EXPERIMENTAL': 'Molar absorptivity [M$^{-1}  \\cdot $cm$^{-1}$]'
        }

def plot_spectrum(exc, osc, unit_in='ENERGY: eV',
        nconf=1, fwhm=0.1, temp=0.0, refraction=1.0, 
        ctype='lorentzian', kind='CROSS_SECTION', plot=True):
    """Plot a spectrum based on theoretical data points.

    This is the main function of the module that should be used
    to calculate and plot electronic spectra.

    Args:
        exc: Input excitation energies given in a one dimenssional ``np.array()``.

        osc: Input oscillator strengths given in a one dimenssional ``np.array()``.

        unit_in (str, optional): String describing the unit used for the input 
            excitation energies. The string must correspond to one of the keys of 
            the dictionary ``convert_to_joules`` in the :py:mod:`~comp_chem_utils.conversions`
            module. Default is ``'ENERGY: eV'``.

        nconf (int, optional): Total number of conformations used in the input data.
            This is used for normalization to a single structure. Default is 1.

        fwhm (float, optional): Full width at Half Maximum used in the convolution function. 
            It must be given in the same units as ``unit_in``. The Default value is ``None``,
            which will be latter changed to correspond to 0.1 eV.

        temp (float, optional): Working temperature in Kelvin. The default value is ``None``, 
            which means that no temperature effect will be calculated.

        refraction (float, optional): The refraction index of the medium. Default value is 1.0.

        ctype (str, optional): Defines the type of convolution. Either ``'lorentzian'`` which is default 
            or ``'gaussian'``.

        kind (str, optional): String describing the type of spectrum that should be calculated.
            It has to be one of the following::

                'STICKS': 'Oscillator strength [Arbitrary units]', 
                'SPECTRAL_FUNC': 'Spectral function [s. $\\cdot $ molecules$^{-1}$]', 
                'CROSS_SECTION': 'Cross section [$\\AA^2 \\cdot $ molecules$^{-1}$]', 
                'EXPERIMENTAL': 'Molar absorptivity [M$^{-1}  \\cdot $cm$^{-1}$]'

            The default is ``kind='CROSS_SECTION'``.

        plot (bool, optional): To control wether to plot or not the spectrum. The default
            is to plot it.

    Return:
        xpts, ypts

        Those are two ``np.array()`` containing the grid points required to plot the desired spectrum.
        The xpts array contains the x-axis values in ``unit_in`` unit, while the ypts array contains
        the spectrum intensity with units depending on the chosen ``kind``.

        By default or if ``plot=True``, the function will plot the desired spectrum and show it
        using the ``matplotlib`` tool. 
    """

    if kind=='STICKS':
        plot_stick_spectrum(exc, osc, color='b')
        if not plot:
            print('That was quite stupid...')
            return exc, osc

    else:
        # get spectral function in seconds per molecule
        xpts, ypts = spectral_function(exc, osc, unit_in=unit_in, nconf=nconf, fwhm=fwhm, ctype=ctype)

        if kind!='SPECTRAL_FUNC':
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



