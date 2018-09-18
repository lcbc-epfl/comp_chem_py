PDOS
====

..
   OLD README FILE STUFF
   [![Documentation Status](https://readthedocs.org/projects/pdos/badge/?version=latest)](http://pdos.readthedocs.io/en/latest/?badge=latest)
   # PDOS
   Clone the repository, 
   ```
   $ git clone https://gitlab.com/pablobaudin/comp_chem_py/PDOS.git  
   ```
   - [Documentation](http://pdos.readthedocs.io/)
   - Version: 2.0
   - Licensed under [LGPLv2.1](./LICENSE)


PDOS stands for Partial Density of States.
PDOS.py is a program which automatically extracts information from Gaussian 
(http://www.gaussian.com/) output and calculates PDOS based upon Lowdin or 
Mulliken orbital analysis. 

The first version of the program was written in January 2011 as part of 
a bachelor project under the supervision of Prof. Mark E. Casida at the
Joseph Fourier University (Grenoble I), France.

Another program with this functionality is the Python program GaussSum 
(http://gausssum.sourceforge.net/) and PDOS.py has been checked against GaussSum. 

An advantage of PDOS.py is that you can plot multiple PDOS as well as the 
total density of states (DOS) on the same graph.  We needed this for our
project and it does not seem to be very easy to do with GaussSum.

GaussSum and PDOS.py differ in their definitions of the gaussian convolution 
which is done.  In GaussSum, the gaussians always have unit HEIGHT.  In
PDOS.py, the gaussians always have unit AREA.  This latter choice seems more
logical to us.  This means that the ratio of peak heights calculated with
GaussSum to that of PDOS.py is, 

.. math::
   \frac{GaussSum}{PDOS.py} = \frac{1}{2} \sqrt{ \frac{\pi}{\log(2)} } \cdot FWHM

where FWHM is the full width at half maximum.


.. include:: ../src/PDOS/doc/pdos.rst
.. include:: ../src/PDOS/doc/gaussian.rst
.. include:: ../src/PDOS/doc/running.rst
.. include:: ../src/PDOS/doc/example.rst
.. include:: ../src/PDOS/doc/display.rst
