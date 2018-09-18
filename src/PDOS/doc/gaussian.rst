Running Gaussian
================

It is imparative that Gaussian be run with the right options so that all
needed information can be extracted from the Gaussian (.log) output file.
We recommend the options ::

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

Alternatively, the code can easily be modified to read the overlap matrix
and the MO energies and coefficients from other sources.


