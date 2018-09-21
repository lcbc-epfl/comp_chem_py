#!/usr/bin/env python

from distutils.core import setup

packages=[
'comp_chem_utils',
'cpmd_scripts',
'PDOS',
'PDOS.modules',
'simp_sec_quant',
]

package_dir={
        '':'src'
        }

setup(name='comp_chem_py',
      version='1.0',
      description='Computational chemistry library in python',
      author='Pablo Baudin',
      author_email='pablo.baudin@epfl.ch',
      package_dir=package_dir,
      packages=packages,
     )

