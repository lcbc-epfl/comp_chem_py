Information from quantum chemistry packages
-------------------------------------------

The central function of the ``PDOS`` program, :py:func:`~PDOS.PDOS_modules.calculate_and_plot_pdos`, 
requires information from a quantum chemistry calculation.
In the current version, this information is assumed to be comming from
the Gaussian package. This information is read through the function
:py:func:`~comp_chem_utils.read_gaussian.get_gaussian_info`.

This function is going through the output file generated by Gaussian.
It is imperative that Gaussian has been run with the right options so that all
needed information can be extracted from output file.
We thus recommend the following options::

   pop=full iop(3/33=1,3/36=-1)

The output file should then include,

#. The number of basis functions and electrons in the system (closed shell is assumed), ::

     Nbasis = xxxx
     xxx alpha electrons

#. The overlap matrix preceded by the line, ::

      *** Overlap *** 

#. The eigenvalues and the molecular orbital coefficients. preceded by the lines, ::

     Molecular Orbital Coefficients
     EIGENVALUES

Alternatively, a new function (similar to :py:func:`~comp_chem_utils.read_gaussian.get_gaussian_info`)
can be added to the code in order to extract the relevant information
from another quantum chemistry package. 


