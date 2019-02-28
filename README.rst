|Documentation Status|
|MIT license| 
.. |PyPI status|
.. |PyPI version shields.io|
|DOI|

.. |Documentation Status| image:: https://readthedocs.org/projects/comp_chem_py/badge/?version=latest
   :target: http://comp_chem_py.readthedocs.io/?badge=latest
   
.. |MIT license| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: LICENSE
   
.. |PyPI status| image:: https://img.shields.io/pypi/status/comp_chem_py.svg
   :target: https://pypi.python.org/pypi/comp_chem_py/

.. |PyPI version shields.io| image:: https://img.shields.io/pypi/v/comp_chem_py.svg
   :target: https://pypi.python.org/pypi/comp_chem_py/

.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2580170.svg
   :target: https://doi.org/10.5281/zenodo.2580170

comp_chem_py: python library for computational chemistry
========================================================

The comp_chem_py package is a collections of python modules and scripts 
that can be usefull in computational chemistry.
The package evolves with my needs in the field.
Feel free to re-use and modify as you wish and :ref:`contact` me if you have any question or comment.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   comp_chem_utils
   cpmd_scripts
   PDOS

Citation
========

If you use this library in a program or publication, please add the following reference:

* *comp_chem_py, a python library for computational chemistry*, Pablo Baudin,
  (Version vX.Y), (2019). http://doi.org/10.5281/zenodo.2580170


Set up and installation
=======================

Dependencies
------------

The comp_chem_py library depends on the following modules:

* scipy_
* numpy_
* matplotlib_
* MySQLdb_

When installing the comp_chem_py library with pip_ those modules will be installed along if needed.

Using git and Linux/UNIX
------------------------

In order to use the library, first clone it::

   git clone --recursive https://gitlab.com/pablobaudin/comp_chem_py.git

Then export the PATH in your ``~/.bashrc``::

   COMP_CHEM_PATH=/path/to/comp_chem_py/root/dir
   export PATH=${COMP_CHEM_PATH}/bin:$PATH
   export PYTHONPATH=${COMP_CHEM_PATH}/external:$PYTHONPATH
   export PYTHONPATH=${COMP_CHEM_PATH}/src:$PYTHONPATH

Using pip
---------

pip_ is the standard Python package manager. 
To install ``comp_chem_py`` with pip_:: 

   pip install comp_chem_py

This will install the ``comp_chem_py`` package with all available dependencies as regular pip controlled packages.

Documentation
=============

Please take a look at the code `documentation <http://comp_chem_py.readthedocs.io/>`_ for more details. 


.. todo::
   #. Add tests. Maybe with doctest. They should also be included in the setup.py. See,
         https://python-packaging.readthedocs.io/en/latest/testing.html

   #. Consider making comp_chem_utils an external submodule.

   #. test ./setup.py script with pip.

   #. Find out how to deal with the dependency on MySQLdb package.

   #. Add the ebsel package into setup.py using the git repo as,
         dependency_links=['http://github.com/user/repo/tarball/master#egg=package-1.0']

.. _contact:

Contact
=======

* Pablo Baudin
* Scientific collatoborator at `LCBC EPFL <https://lcbc.epfl.ch>`_.
* `pablo.baudin@epfl.ch <pablo.baudin@epfl.ch>`_

