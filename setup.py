#!/usr/bin/env python

from distutils.core import setup

def readme():
    with open('README.rst') as f:
        return f.read()

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

scripts=[
        'bin/basis_count',
        'bin/cell_size',
        'bin/group_of_AOs',
        'bin/MSDB',
        'bin/PDOS',
        'bin/PDOS_GUI',
        'bin/plot_PES',
        'bin/pyrun',
        'bin/select_state',
        'bin/truncate_solvent',
    ]

dependencies=[
        'matplotlib',
        'numpy',
        'scipy',
        'MySQLdb',
        #
        # Other imported modules that I probably don't need to specify here...
        #'itertools',
        #'tkinter',
        #'Tkinter',
        #'argparse',
        #'sys',
        #'shutil',
        #'os',
        #'math',
        #'subprocess',
    ]

classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Unix Shell',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ]

setup(name='comp_chem_py',
        version='1.0',
        description='Computational chemistry library in python',
        long_description=readme(),
        license='MIT',
        author='Pablo Baudin',
        author_email='pablo.baudin@epfl.ch',
        url='https://comp-chem-py.readthedocs.io/en/latest/',
        package_dir=package_dir,
        packages=packages,
        scripts=scripts,
        install_requires=dependencies,
        classifiers=classifiers,
        # this is to include files in MANIFEST.in file
        include_package_data=True,
    )



