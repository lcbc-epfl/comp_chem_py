#!/usr/bin/env python
"""Collection of functions to read and manipulate CPMD output files"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import sys
import os
import numpy as np

from comp_chem_utils.molecule_data import read_xyz_table, xyz_file
from comp_chem_utils.utils import get_file_as_list, get_lmax_from_atomic_charge
from comp_chem_utils.conversions import AU_TO_PS
from comp_chem_utils.periodic import element


def read_standard_file(fn):
    """Read standard trajectory file and export it.

    Here a trajectory file is understood as a file arranged
    as columns in which the first column contains the
    step indices and the others any type of informations.

    For example the ENERGIES or SH_STATE.dat files enter
    in this category.

    Args:
        fn (str): Name of the trajectory file.

    Returns:
        steps, info

        steps is as list of step indices (int) as read from
        the first column of the trajectory file, while info
        is an array (``np.array()``) of floats containing the
        rest of the information.
    """

    lines = get_file_as_list(fn)

    steps = [int(l.split()[0]) for l in lines]
    ninfo = len(lines[0].split()) - 1

    info = np.zeros( (len(steps), ninfo) )
    for i, line in enumerate(lines):
        info[i,:] = [float(line.split()[j+1]) for j in range(ninfo)]

    return steps, info

def read_TRAJEC_xyz(fn):
    """Read TRAJEC.xyz file and export it.

    An xyz type information is read for each step in the trajectory file,
    where xyz info is assumed to have the following format::

        number of atoms
        title line = step index
        atom 1 label  and corresponding xyz coordinate
        atom 2 label  and corresponding xyz coordinate
        : : :
        atom N label  and corresponding xyz coordinate

    Args:
        fn (str): Name of the TRAJEC.xyz file.

    Returns:
        steps, traj_xyz

        steps (list): List of step indices (int) as read from
            the title line of each xyz_data block.

        traj_xyz (list): List of xyz_data blocks. Each block
            is itself a list of the lines containing the coordinates
            in the xyz format. Each line is a 4 item list with first
            index the atom symbol and the last 3 items are the xyz
            coodinates.
    """

    lines = get_file_as_list(fn)
    natoms = int(lines[0])

    # 1 xyz_data consist of natoms lines plus the two first lines (natoms + step index)
    traj_xyz = []
    steps = []
    for xyz_data in (lines[x:x+natoms+2] for x in range(0,len(lines),natoms+2)):
        # test that the xyz_data is not empty
        if (xyz_data[0].strip()):
            # read_xyz_table return a list of length natoms
            # each line is a list with the atom's symbol and the 3 coordinates
            traj_xyz.append( read_xyz_table(xyz_data[2:natoms+2]) )
            steps.append( int(xyz_data[1].split()[1]) )

    return steps, traj_xyz

def read_GEOMETRY_xyz(fn):
    """Like read_TRAJEC_xyz for a single geometry and without the step index."""

    lines = get_file_as_list(fn)
    natoms = int(lines[0])
    return read_xyz_table(lines[2:natoms+2])

def write_TRAJEC_xyz(steps, traj_xyz, output):
    """Write a TRAJEC.xyz file (CPMD style) to the output file.

    Args:
        steps (list): Step indices. See read_TRAJEC_xyz() function.

        traj_xyz (list): xyz data. See read_TRAJEC_xyz() function.

        output (str): Name (and path) fo the file in which the
            information will be written.
    """

    with open(output, 'w') as myf:
        for s, x in zip(steps, traj_xyz):
            myf.write('{:12d}\n'.format(len(x)))
            myf.write(' STEP:{:12d}\n'.format(s))
            for line in x:
                #myf.write('{:>2} {:13.6f} {:13.6f} {:13.6f}\n'.format(*line))
                myf.write('{0[0]:<2} {0[1]:13.6f} {0[2]:13.6f} {0[3]:13.6f}\n'.format(line))


def split_TRAJEC_data(steps, traj_xyz, verbose=True, interactif=True, write_xyz=False,
        start_i=1, nstep=100, delta=None, dt=5, name=''):
    """Select equidistant xyz data snapshot from trajectory data.

    From a starting index, a number of snapshots and a number of
    steps between each snapshot. A new trajectory data is generated
    with only a subset of steps.

    Args:
        steps (list): Step indices. See read_TRAJEC_xyz() function.

        traj_xyz (list): xyz data. See read_TRAJEC_xyz() function.

        verbose (bool, optional): Print more information while running.
            Default ``is True``.

        interactif (bool, optional): Let the user chose the parameters
            interactively. Default ``is True``.

        write_xyz (bool, optional): Write an xyz file to disk for every step.
            Default ``is False``.

        start_i (int, optional): Step index for the first snapshot.
            Default is 1

        nstep (int, optional): Total number of final snapshots.
            Default is 100.

        delta (int, optional): Number of steps between two snapshots.
            Default is ``None``. It will be changed to the maximum possible
            value depending on other paramters.

        dt (int, optional): Time step used in MD [in a.u.] (to provide
            print out the actual time between the snapshots). Default
            is 5 a.u.

        name (str, optional): String/Title to be used in xyz filename
            in case ``write_xyz = True``.

    Return:
        new_steps, new_traj_xyz

        Just a subset of ``steps`` and ``traj_xyz``.
    """

    if verbose:
        print("Total number of steps    : {}".format(len(steps)))
        print("First and last step index: {} - {}".format(steps[0], steps[-1]))
        print("Delta between two steps  : {}\n".format(steps[1]-steps[0]))

    if interactif:
        # ask user to chose the relevant snapshots
        start_i = int( raw_input("Step index for first snapshot (between 0 and {}): \n".format(len(steps)-1)))
        nstep = int( raw_input("Total number of snapshot (between 1 and {}): \n".format(len(steps)-start_i)))
        delta = (len(steps) - start_i)/(nstep-1) - 1
        inp = raw_input("Delta between two step indices (max and default = {}): \n".format(delta))
        if inp:
            delta = int(inp)

        dt = float( raw_input("Time step used in MD [in a.u.] (to provide time between snapshots): \n"))
        write_xyz = raw_input("Write xyz file ? [y/n]\n").strip().lower() == 'y'

        if write_xyz:
            name = raw_input("String/Title to be used in xyz filename:\n")

    else:
        # use input or default value
        if not delta:
            delta = (len(steps) - start_i)/(nstep-1) - 1

    # time between snapshots is given by d_steps * dt [a.u.]
    # where d_steps is the actual number of steps between two snapshots:
    d_steps = steps[start_i + delta] - steps[start_i]

    # stoping index in steps list
    stop_i = start_i + delta*nstep

    # sanity check
    if stop_i>len(steps):
        stop_i = len(steps)
        if verbose:
            print("WARNING new number of steps: {}".format( len(range(start_i, stop_i, delta)) ))

    if verbose:
        print("Starting index:      {}".format( start_i ))
        print("Last index:          {}".format( stop_i ))
        print("Number of snapshots: {}\n".format( nstep ))

        print("# steps between two snapshots: {}".format(d_steps))
        print("Time between snapshots [a.u.]: {}".format(d_steps * dt))
        print("Time between snapshots [ps]:   {}\n".format(d_steps * dt * AU_TO_PS))


    # Loop over chosen configurations
    new_traj_xyz = []
    new_steps = []
    for i in range(start_i, stop_i, delta):

        new_steps.append( steps[i] )
        new_traj_xyz.append( traj_xyz[i] )

        # create xyz file
        if write_xyz:
            fname = '{:010d}_{}.xyz'.format(new_steps[-1], name)

            xyzf = xyz_file()
            xyzf.read_from_table(new_traj_xyz[-1], title=name)
            xyzf.out_to_file(fname=fname)

    return new_steps, new_traj_xyz


def read_GEOMETRY(fn):
    """Simplified version of read_FTRAJECTORY for a single structure."""
    data = np.loadtxt(fn)
    xyz = data[:,0:3]
    vel = data[:,3:]
    return xyz, vel


def read_TRAJECTORY(fn, verbose=False):
    """Just a wrapper to read_FTRAJECTORY."""
    return read_FTRAJECTORY(fn, forces=False, verbose=verbose)


def read_FTRAJECTORY(fn, forces=True, verbose=False):
    """Read the FTRAJECTORY file from a CPMD run and export it.

    Three different formats can be read.

    #. The TRAJECTORY file (with ``forces=False``)::

        Column 0: step index
        Column 1-3: xyz coordinates
        Column 4-6: velocities

    #. The FTRAJECTORY file (with ``forces=True``)::

        Column 0: step index
        Column 1-3: xyz coordinates
        Column 4-6: velocities
        Column 7-9: forces

    #. The FTRAJECTORYMTS file (with ``forces=True``, this file is not produced anymore)::

        Column 0: step index
        Column 1-3: low level forces
        Column 4-6: high level forces
        Column 7-9: xyz coordinates

    """

    if verbose:
        print('INFO: entering read_FTRAJECTORY')

    stp = []
    xyz = []
    vel = []
    if forces:
        fce = []

    lines = get_file_as_list(fn)
    if verbose:
        print('INFO: file read')

    # read total number of steps from last line
    nstep = int(lines[-1].split()[0])

    # get number of atoms as (# lines) / (# steps)
    nlines = len(lines)
    natoms = get_natoms(lines)

    if verbose:
        print('INFO: number of steps: {}'.format(nstep))
        print('INFO: number of atoms: {}'.format(natoms))

    blocks = [ lines[x:x+natoms] for x in range(0, nlines, natoms) ]

    if verbose:
        print('INFO: blocks read')


    for blk in blocks:

        # step index from first column of first line in block
        stp.append( int(blk[0].split()[0]) )

        if verbose:
            print('INFO: read step: {:10d} of {:10d}'.format(stp[-1], nstep))

        myxyz = np.zeros( (natoms, 3) )
        myvel = np.zeros( (natoms, 3) )
        if forces:
            myfce = np.zeros( (natoms, 3) )
        for i, ml in enumerate(blk):

            l = ml.split()
            myxyz[i,:] = [float(l[x]) for x in range(1,4)]
            myvel[i,:] = [float(l[x]) for x in range(4,7)]
            if forces:
                myfce[i,:] = [float(l[x]) for x in range(7,10)]

        xyz.append( myxyz )
        vel.append( myvel )
        if forces:
            fce.append( myfce )

    if forces:
        return stp, xyz, vel, fce
    else:
        return stp, xyz, vel


ref_code = {
        'steps' :-1,
        'E_kel' : 0,
        'Temp'  : 1,
        'E_KS'  : 2,
        'E_cla' : 3,
        'E_ham' : 4,
        'RMS'   : 5,
        'CPU_t' : 6
        }

def read_ENERGIES(fn, code, factor=1, HIGH=True):
    """Read ENERGIES file from CPMD and return data based on code:

    Args:
        fn (str): filename of the ENERGIES file (can include path to the file).

        code: List of codes written as strings describing which information
            should be extracted from the file. The available codes are::

                'steps' : Step indices
                'E_kel' : Electronic kinetic energy (only for CPMD)
                'Temp'  : Temperature [K]
                'E_KS'  : Kohn-Sham energy [a.u.]
                'E_cla' : Classical energy, E_KS + E_kin (constant for BOMD)
                'E_ham' : 0 for BOMD
                'RMS'   : Nuclear displacement wrt initial position (?)
                'CPU_t' : CPU time

        factor (int): Integer factor used to skip some information. E.g. for
            an MTS calculation the high level information can be extracted by
            setting factor equals to the MTS factor used in the calculation.
            Default value is 1 (every time step info is extracted).

        HIGH (bool): define wether to extract the information every ``factor``
            step (this is default, `HIGH=True`) or the negative counter part
            meaning the info is extracted for every step except every factor
            steps `HIGH=False`.

    Return:
        The function returns a dictionary with keys the input codes and
        with values an array (``np.array``) containing the corresonding
        information as ``floats``. The exception is the value for ``'steps'``
        which is a simple list of integers.

        For example if the input codes are ``code=['steps','E_cla']``, the
        output dictionary will have the form::

            >>> read_ENERGIES('ENERGIES', ['steps','E_cla'])
            {'E_cla': array([-546.99550862, -546.99770079, ..., -546.96549996]),
            'steps': [1, 2, ..., 10000]}

    """

    steps, info = read_standard_file(fn)
    myrange = range(factor-1,len(steps),factor)

    to_return = {}
    for i in code:
        if HIGH:
            # extract information every factor steps
            if i=='steps':
                to_return[i] = [steps[j] for j in myrange]
            else:
                l = info[:,ref_code[i]]
                to_return[i] = [l[j] for j in myrange]
        else:
            # extract negative info (i.e. all except every factor steps)
            if i=='steps':
                to_return[i] = [steps[j] for j in range(len(steps)) if j not in myrange ]
            else:
                l = info[:,ref_code[i]]
                to_return[i] = [l[j] for j in range(len(l)) if j not in myrange]

    return to_return


def read_SH_ENERG(fn, nstates=None, factor=1):
    """Read SH_ENERG.dat file and exctract info in dictionary"""
    steps, info = read_standard_file(fn)

    sh_data = {}
    sh_data['steps'] = [ x*factor for x in steps ]

    maxstates = len(info[0,:]) - 1
    if not nstates:
        nstates = maxstates

    if nstates > maxstates:
        nstates = maxstates

    for i in range(nstates):
        sh_data[ 'State {}'.format(i) ] = info[:,i]

    sh_data[ 'Driving State' ] = info[:,-1]
    return sh_data


def read_MTS_EXC_ENERG(fn, nstates, MTS_FACTOR, HIGH):
    """Read MTS_EXC_ENERG.dat file and exctract info in dictionary.

    Either the HIGH or the LOW level info will be exctracted.
    """

    steps, info = read_standard_file(fn)

    myrange = range(1, len(info[:,0]), (MTS_FACTOR+1))

    if HIGH:
        SUB='HIGH'
        sub_info = np.array([info[x,:] for x in myrange] )
        sub_steps = np.array([steps[x] for x in myrange] )

    else:
        SUB='LOW'
        sub_info = np.array([info[x,:] for x in range(len(info[:,0])) if x not in myrange] )
        sub_steps = np.array([steps[x] for x in range(len(info[:,0])) if x not in myrange] )

    mts_data = {}
    mts_data['steps'] = sub_steps

    maxstates = len(sub_info[0,:]) - 1
    if not nstates:
        nstates = maxstates

    if nstates > maxstates:
        nstates = maxstates

    for i in range(nstates):
        mts_data[ 'E({}) State {}'.format(SUB, i) ] = sub_info[:,i]

    mts_data[ 'E Driving' ] = sub_info[:,-1]
    return mts_data


def get_time_info(fn):
    """Read CPMD input file and extract time step info."""

    # set defaults
    TIMESTEP = 5
    MAXSTEP = 10000
    USE_MTS = False
    MTS_FACTOR = 1

    # read input file
    lines = get_file_as_list(fn)

    # parse it
    for j, line in enumerate(lines):
        if 'TIMESTEP' in line and 'FACTOR' not in line:
            TIMESTEP = float(lines[j+1].split()[0])
        elif 'MAXSTEP' in line:
            MAXSTEP = int(lines[j+1].split()[0])
        elif 'BOMD_FORCES MTS' in line:
            USE_MTS = True
        elif 'TIMESTEP_FACTOR' in line:
            MTS_FACTOR = int(lines[j+1].split()[0])

    return TIMESTEP, MAXSTEP, USE_MTS, MTS_FACTOR

def get_natoms(lines):
    """Get the number of atoms from a file like TRAJECTORY."""
    istep = lines[0].split()[0]
    natoms=0
    for line in lines:
        new_i = line.split()[0]
        if new_i != istep:
            return natoms
        else:
            natoms+=1

def xyz_to_cpmd_atoms(xyz_data=None, xyz_filename=None, PPs=None):
    """Return list of strings as in CPMD ATOMS section."""

    xyzf = xyz_file()
    if xyz_data:
        xyzf.read_from_table(xyz_data)
    elif xyz_filename:
        xyzf.read_xyz(xyz_filename)
    else:
        sys.exit('Wrong input!!')

    lines = []
    lines.append('&ATOMS')
    # print atom types
    for atype in xyzf.atom_types:

        lines.append('*{}_{}'.format(atype.symb, PPs))

        lmax = get_lmax_from_atomic_charge(atype.charge)

        lines.append('   LMAX={}'.format(lmax) )
        lines.append('   {}'.format(atype.natoms) )

        for x,y,z in zip(atype.xvals, atype.yvals, atype.zvals):
            line = "   {:20.10f} {:20.10f} {:20.10f}".format( x, y, z )
            lines.append( line )

    lines.append('&END')

    return lines

def qmmm_cpmd_atoms(atom_types, PPs=None):
    """Return list of strings as in CPMD ATOMS section for QMMM calculation."""

    lines = []
    lines.append('&ATOMS')
    # print atom types
    for atype, idx in atom_types.items():

        lines.append('*{}_{}'.format(atype, PPs))

        atomic_charge = element(atype).atomic
        lmax = get_lmax_from_atomic_charge(atomic_charge)

        lines.append('   LMAX={}'.format(lmax) )
        lines.append('   {}'.format( len(idx)) )

        idx_str = ' '.join( [str(x) for x in idx ] )
        lines.append( '   ' +idx_str )

    lines.append('&END')

    return lines

def xyz_data_to_np(xyz_data):
    """Convert xyz_data information to a tuple (atoms, coord)

    Where atoms is a list of the atomic symbols and coord is
    a matrix: np.array[natoms, 3]"""

    atoms = [ l[0] for l in xyz_data]
    coord = [ l[1:] for l in xyz_data]
    return (atoms, np.array(coord))

if __name__ == "__main__":
    # TODO: actually test the result of the functions

    print("Testing {} module\n".format(os.path.basename(__file__)))

    fold = './tests/CPMD_files/'

    fn = '{}SH_ENERG.dat'.format(fold)
    read_standard_file(fn)
    print("Test of function: read_standard_file: OK")

    # read ./test/TRAJEC.xyz file and print it to new file
    fn = "{}TRAJEC.xyz".format(fold)
    steps, traj_xyz = read_TRAJEC_xyz(fn)
    print("Test of function: read_TRAJEC_xyz   : OK")

    output = '{}TRAJEC_1.xyz'.format(fold)
    write_TRAJEC_xyz(steps, traj_xyz, output)
    print("Test of function: write_TRAJEC_xyz  : OK")
    #os.rm(output)

    fn = '{}FTRAJECTORY'.format(fold)
    read_FTRAJECTORY(fn)
    print("Test of function: read_FTRAJECTORY  : OK")

    fn = '{}ENERGIES'.format(fold)
    read_ENERGIES(fn, ref_code)
    print("Test of function: read_ENERGIES     : OK")
    print("\nAll tests were successfully executed! :D")


