PDOS
====

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
For example of applications see [Wawire2013a]_ and [Magero2017]_.

.. include:: ../src/PDOS/doc/pdos.rst
.. include:: ../src/PDOS/doc/gaussian.rst
.. include:: ../src/PDOS/doc/running.rst
.. include:: ../src/PDOS/doc/example.rst
.. include:: ../src/PDOS/doc/display.rst
..
.. include:: ../src/PDOS/doc/code.rst

References
----------

.. [Wawire2013a] Muhavini Wawire, C. *et al.*, (2014). 
   Density-functional study of luminescence in polypyridine ruthenium complexes. 
   *Journal of Photochemistry and Photobiology A: Chemistry*, **276**, 8–15. 
   http://doi.org/10.1016/j.jphotochem.2013.10.018.

.. [Magero2017] Magero, D. *et al.*, (2017). 
   Partial density of states ligand field theory (PDOS-LFT): Recovering a LFT-like picture 
   and application to photoproperties of ruthenium(II) polypyridine complexes. 
   *Journal of Photochemistry and Photobiology A: Chemistry*, **348**, 305–325. 
   http://doi.org/10.1016/j.jphotochem.2017.07.037


