Installing and running PDOS.py 
------------------------------

To run the ``PDOS`` script, you must first create a text file with the orbital groups 
for which you want to plot PDOS.  This file has the following format, ::

   s orbitals
   1,2,6,10,11,15
   p orbitals
   3-5,7-9,12-14,16-18

Here the lines ``s orbitals`` and ``p orbitals`` are labels for each PDOS curve.  What
follows is a index list of the AOs whose individual PDOS will be summed to make the curve.
The AOs are in the same order as in the Gaussian output file.
For convenience, you can create this file with the Python script ``group_of_AOs``.
Just type, ::

   group_of_AOs gaussian_file

and follow the instructions.
You should then be ready to use the PDOS script as follows, ::

   PDOS gaussian_file group_of_AOs.inp

The use of default parameters is recommended.
The graphical interface can also be used, (it might be more convenient to change
the default settings such as the energy range or the FWHM.). In order to 
access the GUI just type in your terminal, ::

   PDOS_GUI

then there will be a pop-up graphical interface with boxes to fill in,
only the two input files are required, the other field can be left 
intact in order to use the default parameters.
Click on the button marked ``PLOT`` and wait for it.


