.. PDOS documentation master file, created by
   sphinx-quickstart on Sun Oct  9 18:00:01 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

About the PDOS program
======================

.. toctree::
   :maxdepth: 2

   pdos
   gaussian
   running
   example
   display

PDOS.py is a program which automatically extracts information
from Gaussian (http://www.gaussian.com/) output and calculates partial density 
of states (PDOS) based upon Lowdin or Mulliken orbital analysis.  
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


