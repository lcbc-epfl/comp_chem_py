|Documentation Status|
|MIT license| 

.. |Documentation Status| image:: https://readthedocs.org/projects/comp_chem_py/badge/?version=latest
   :target: http://comp_chem_py.readthedocs.io/?badge=latest
   
.. |MIT license| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: LICENSE
   

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

* Please take a look at the code `documentation <http://comp_chem_py.readthedocs.io/en/latest/>`_ for more details. 


.. todo::
   #. Keep improving documentation.
   #. Remove old_comp_chem_py repository/group and all subrepositories or archive it.
   #. Add tests. Maybe with doctest.
   #. Consider making comp_chem_utils an external submodule.
   #. Add section about requirements and dependencies on documentation.


