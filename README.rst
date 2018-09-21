|Documentation Status|

.. |Documentation Status| image:: https://readthedocs.org/projects/comp_chem_py/badge/?version=latest
   :target: http://comp_chem_py.readthedocs.io/?badge=latest
   
.. [![License](https://img.shields.io/badge/license-%20BSD--3-blue.svg)](../master/LICENSE)

Set up and installation
=======================

In order to use the library, first clone it::

   git clone --recursive https://gitlab.com/pablobaudin/comp_chem_py.git

Then export the PATH in your ``~/.bashrc``::

   COMP_CHEM_PATH=/path/to/comp_chem_py/root/dir
   export PATH=${COMP_CHEM_PATH}/bin:$PATH
   export PYTHONPATH=${COMP_CHEM_PATH}/external:$PYTHONPATH
   export PYTHONPATH=${COMP_CHEM_PATH}/src:$PYTHONPATH


Documentation
=============

* [Master branch](https://comp-chem-py.readthedocs.io/en/latest/)


.. todo::
   #. Keep improving documentation.
   #. Remove old comp_chem_py libraries and rename the new one.
   #. Make library open source.
   #. Link documentation to read the docs.
   #. Add link to documentation in this README file
   #. Add tests. Maybe with doctest.
   #. Consider making comp_chem_utils an external submodule.


