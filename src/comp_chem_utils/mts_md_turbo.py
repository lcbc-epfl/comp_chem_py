#!/usr/bin/env python
"""Multiple time step algorithm for Molecular Dynamics with Turbomole.

Default usage
-------------

The following files are required:
    - GEOMETRY: a CPMD file for starting positions and velocities.
    - define_high.inp: Input to TURBOMOLE's define script to set up
        the high level electronic structure calculations.
    - define_low.inp: same as define_high.inp but for low level.

From a python script or interpreter::

    # setup the system, the list of atoms needs to be the same as in the GEOMETRY file.
    system = system_settings(['O','H','H'])

    # Run MD
    data = md_driver(system)

The output dictionary ``data``, contains the positions, velocities, forces, and energies
along the MD trajectory.

The MD and TURBOMOLE parameters can be modified by setting the ``md_settings`` and
``turbo_settings`` objects and passing them to the ``md_driver``.
"""

import sys
import os
from subprocess import check_call

import numpy as np

from comp_chem_utils.utils import get_file_as_list
from comp_chem_utils.cpmd_utils import read_GEOMETRY
from comp_chem_utils.conversions import BOHR_TO_ANG, ATOMIC_MASS_AU, AU_TO_KELVIN
from comp_chem_utils.periodic import element
from cpmd_scripts.MTS_scripts.turbo_call import prepare_and_run_turbo

try:
    # enable management of modules through python
    execfile('/software/modules/3.2.10/Modules/3.2.10/init/python.py')
except:
    print('WARNING: Failed to enable module environment.')
    print('WARNING: TURBOMOLE cannot be loaded!')

class system_settings(object):

    def __init__(self, atoms, linear=False):
        self.atoms = atoms
        self.natoms = len(atoms)
        self.mass = np.array([element(x).mass * ATOMIC_MASS_AU for x in atoms])
        # degrees of freedom
        if self.natoms == 1:
            self.NDOF = 0
            sys.exit('Cannot do MD with 1 atom...')
        elif linear:
            print("Warning: Linear Molecule, you sure?!")
            self.NDOF = 3 * self.natoms - 5
        else:
            self.NDOF = 3 * self.natoms - 6

class md_settings(object):

    def __init__(self, max_iter=10000, time_step=10.0, mts_factor=1):
        self.time_step = time_step
        self.max_iter = max_iter
        self.mts_factor = mts_factor

    def output(self):
        print("RUNNING MTS BO MOLECULAR DYNAMICS\n")
        print("     time step  = {} a.u.".format(self.time_step))
        print("     Max iter   = {}".format(self.max_iter))
        print("     MTS factor = {}".format(self.mts_factor))

class turbo_settings(object):

    def __init__(self, high_input='define_high.inp', low_input='define_low.inp', nnodes=24):
        self.high_input = high_input
        self.low_input = low_input
        self.nnodes = nnodes

    def load(self):
        """Load the environment necessary to run TURBOMOLE."""
        module('load','turbomole/7.1.1')
        os.environ["PARA_ARCH"] = "SMP"
        os.environ["PARNODES"] = str(self.nnodes)
        #module('purge')
        #os.environ["PATH"] += os.environ['TURBODIR'] + '/bin/em64t-unknown-linux-gnu_smp'


def md_driver(sys_in, md_in=md_settings(), turbo_in=turbo_settings()):
    """Run a molecular dynamic using the RESPA algorithm of Tuckerman."""

    global system, md_parm, turbo_parm
    system = sys_in
    md_parm = md_in
    turbo_parm = turbo_in

    md_parm.output()

    # PREPARE TURBOMOLE
    turbo_parm.load()

    nmts = md_parm.mts_factor
    data = initialization(nmts=nmts)

    # current position and velocities
    x = data['pos'][-1]
    v = data['vel'][-1]

    # curent forces
    f_low = data['f_low'][-1]
    f_high = data['f_high'][-1]
    f_tot = f_high * nmts - f_low * (nmts-1)


    # START MD LOOP
    for idx in range(1,md_parm.max_iter):

        # UPDATE VELOCITIES
        v = vel_update(v,f_tot)

        # UPDATE POSITIONS
        x = x + v * md_parm.time_step

        # GET NEW FORCES AND ENERGIES
        if idx%nmts==0:
            # step index is a multiple of MTS factor
            # this a "large" step, we compute High level forces
            f_low, e_low = get_forces_and_energy(x, 'low')
            f_high, e_high = get_forces_and_energy(x, 'high')
            f_tot = f_high * nmts - f_low * (nmts-1)
        else:
            # Just get low level forces as effective forces
            f_low, e_low = get_forces_and_energy(x, 'low')
            f_high, e_high = np.zeros_like(f_low), 0.0
            f_tot = f_low

        # FINAL VELOCITY UPDATE
        v = vel_update(v,f_tot)

        # UPDATE DATA WITH CURRENT VALUES
        data = data_update(data, x, v, f_low, f_high, e_low, e_high, idx, nmts)

    return data



def initialization(fname='GEOMETRY', nmts=1):
    """
    Initialize positions, velocities, forces, and energies.

    The initial positions and velocities are read from a CPMD GEOMETRY file
    in atomic units. (Bohr not Angstrom!)

    The Low and High level forces and energies are then calculated by
    calling the Turbomole program.
    """

    # read GEOMETRY file
    x, v = read_GEOMETRY(fname)

    # Get forces and energies
    f_low, e_low = get_forces_and_energy(x, 'low')
    f_high, e_high = get_forces_and_energy(x, 'high')

    # store everything in the data dictionary
    data = {
            'pos':[],
            'vel':[],
            'f_low':[],
            'f_high':[],
            'e_low':[],
            'e_high':[],
            }
    data = data_update(data, x, v, f_low, f_high, e_low, e_high, idx=0, nmts=nmts)

    return data


def get_forces_and_energy(x, level):
    """Calculate forces and energies at a given level by calling TURBOMOLE."""

    def write_position(fname, x):
        """Write current geometry to disk in Angstroms."""
        with open(fname,'w') as xyz_file:
           xyz_file.write(str(system.natoms)+'\n')
           xyz_file.write(' GEOMETRY from CPMD\n')
           for symb, at in zip(system.atoms, x):
               at = at * BOHR_TO_ANG
               xyz_file.write('{:<2} {:13.6f} {:13.6f} {:13.6f}\n'.format(symb, *at))

    def read_forces(fname):
        """Read forces generated by TURBOMOLE. Assumed to be in a.u."""
        lines = get_file_as_list(fname, raw=True)
        lines = lines[2:]
        f = np.zeros((system.natoms,3))
        for i in range(system.natoms):
            f[i,:] = [float(x) for x in lines[i].split()[1:] ]
        return f

    def read_energy():
        """Read turbomole SCF energy from the energy file generated in the last calculation."""
        lines = get_file_as_list('TURBOMOLE/energy')
        return float(lines[1].split()[1])


    # SELECT TURBOMOLE INPUT BASED ON LEVEL
    if level=='high':
        define_inp = turbo_parm.high_input

    elif level=='low':
        define_inp = turbo_parm.low_input

    else:
        sys.exit('Wrong level in get_forces_and_energy: {}'.format(level))

    geometry = 'geometry.xyz'
    forces = 'forces.xyz'

    # WRITE POSITIONS TO DISK
    write_position(geometry, x)

    # CALL TURBOMOLE
    prepare_and_run_turbo(os.getcwd(), 'DFT', '', geometry, define_inp, forces)

    # READ FORCES AND ENERGY FROM DISK
    f = read_forces(forces)
    e = read_energy()

    return f, e



def vel_update(v,f):
    """Velocity Verlet update of Velocities."""
    for i in range(system.natoms):
        v[i,:] += md_parm.time_step * f[i,:] / (2.0 * system.mass[i])
    return v


def data_update(data, x, v, f_low, f_high, e_low, e_high, idx, nmts):
    """Update the data dictionnary with the current iteration values.

    The data are also written to disk."""

    data['pos'].append( x )
    data['vel'].append( v )
    data['f_low'].append( f_low )
    data['e_low'].append( e_low )
    data['f_high'].append( f_high )
    data['e_high'].append( e_high )

    write_xyz_data('POSITIONS.xyz', data['pos'], conv=BOHR_TO_ANG, nmts=nmts, idx=[idx])
    write_xyz_data('VELOCITIES.xyz', data['vel'], nmts=nmts, idx=[idx])
    write_xyz_data('FORCES_LOW.xyz', data['f_low'], nmts=nmts, idx=[idx])
    write_energies("ENERGIES_LOW", data['e_low'], data['vel'], nmts=nmts, idx=[idx])

    if idx%nmts==0:
        write_energies("ENERGIES_HIGH", data['e_high'], data['vel'], nmts=nmts, idx=[idx])
        write_xyz_data('FORCES_HIGH.xyz', data['f_high'], nmts=nmts, idx=[idx])

    return data

def write_xyz_data(fname, data, conv=1.0, nmts=1, idx=[]):
    """Write or update XYZ type data to disk (positions, velocities, forces...)."""
    if not idx:
        idx = np.array( range(0, len(data), nmts) )

    if 0 in idx:
        open_kind = 'w'
    else:
        open_kind = 'a'

    with open(fname,open_kind) as xyz_file:
        for i,xyz in enumerate(data):
            if i in idx:
                xyz_file.write(str(system.natoms)+'\n')
                xyz_file.write(' STEP: {} \n'.format(i))
                for symb, at in zip(system.atoms, xyz):
                    at = at * conv
                    xyz_file.write('{:<2} {:13.6f} {:13.6f} {:13.6f}\n'.format(symb, *at))

def write_energies(fname, E_pot, v, nmts=1, idx=[]):
    """ print to fname:

        Step: MD step index (depends on nmts)
        E_pot: potential energy in a.u.
        E_kin: kinetic energy in a.u.
        E_tot: Constant energy (E_pot + E_kin) in a.u.
        Temp: Kinetic temperature in Kelvin
    """

    def get_kinetic_energy(v_all):
        """Calculate kinetic energy from velocities."""
        e_kin = np.zeros(len(v_all))
        for idx, v in enumerate(v_all):
            for i in range(system.natoms):
                const = 0.5 * system.mass[i]
                e_kin[idx] += const * (v[i,0]*v[i,0] + v[i,1]*v[i,1] + v[i,2]*v[i,2] )
        return e_kin

    def get_temperature(e_kin):
        """Calculate kinetic temperature from kinetic energy."""
        return e_kin * 2.0 * AU_TO_KELVIN / system.NDOF

    # define arrays to print
    Step = np.array( range(len(E_pot)) )
    E_pot = np.array(E_pot)
    E_kin = get_kinetic_energy(v)
    E_tot = E_pot + E_kin
    Temp = get_temperature(E_kin)

    # Check which step is to ne printed
    if not idx:
        idx = np.array( range(0, len(E_pot), nmts) )

    if 0 in idx:
        open_kind = 'w'
    else:
        open_kind = 'a'

    # PRINT TO FILE
    header = ['Step index', 'E_pot a.u.', 'E_kin a.u.', 'E_tot a.u.', 'Temp Kelv.']
    with open(fname, open_kind) as energy:
        if 0 in idx:
            energy.write('{:^15}   {:^15}   {:^15}   {:^15}   {:^15}\n'.format(*header) )
        for s,ep,ek,et,tp in zip(Step,E_pot,E_kin,E_tot,Temp):
            if s in idx:
                energy.write('{:^15d}   {:15.8e}   {:15.8e}   {:15.8e}   {:15.5f}\n'.format(s,ep,ek,et,tp) )

if __name__=='__main__':

    # INPUT SETTINGS
    system = system_settings(['O','O','H','H','H','H'])
    md_parm = md_settings(high_input='define_high.inp', low_input='define_low.inp', nnodes=24)
    turbo_parm = turbo_settings(max_iter=10000, time_step=10.0, mts_factor=1)

    data = md_driver(system, md_parm=md_parm, turbo_parm=turbo_parm)

