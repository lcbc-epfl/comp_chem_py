#!/usr/bin/env python
"""
Script used to obtain nuclear forces from turbomole based.

The script is made for used in combination with CPMD.
I.e. The nuclear positions are produced along a CPMD trajectory
and Turbomole provide the corresponding nuclear forces.
CPMD is taking care of the propagation.

The nuclear positions are read from:
    - GEO_XYZ [default is "geometry.xyz"]

The nulcear forces are written to:
    - FORCES_XYZ [default is "FORCES.xyz]

In addition a input for Turbomole is required:
    - DEF_INP [default is "define.inp"]

If Turbomole is crashing, the file EXT.error is produced 
and should be detected by CPMD.
"""

import os
import sys
import argparse
import numpy as np
from shutil import copy, rmtree
from subprocess import check_call


# GLOBAL VARIABLE FOR JOB TYPE
RIMP2 = 'RIMP2'
RICCS = 'RICCS'
RICC2 = 'RICC2'
TDDFT = 'TDDFT'
DFT   = 'DFT'
JOBS = [RIMP2, RICCS, RICC2, TDDFT, DFT]
    

def quit_cpmd(statement):
    """write EXT.error file that will make CPMD quit"""
    print(statement)
    with open(os.path.join(CPDIR, 'EXT.error'), 'w') as error:
        error.write(statement)
    sys.exit(statement)


def call_bash(command):
    """execute bash command"""
    return check_call(['bash','-c',command])


def run_turbo(prog, out, inp=''):
    """
    Run TURBOMOLE program 'prog' and put results in 'out'.

    Eventually use 'inp' as input to the program.
    """

    # RUN
    if inp:
        call_bash('{} < {} &> {}'.format(prog, inp, out ))
    else:
        call_bash('{} &> {}'.format(prog, out ))


    # CHECK FOR PROPER TERMINATION
    success = False
    with open(out, 'r') as outf:
        #success = 'all done' in  outf.read().split('\n')[-6]
        success = 'ended normally' in  outf.read().split('\n')[-2]

    if not success:
        quit_cpmd('TURBOMOLE Error in {} see {} file'.format(prog, out))

    return success


def get_natom_from_xyz(XYZ):
    """Read number of atoms from the first line of an XYZ file"""

    with open(XYZ, 'r') as f:
        try:
            natom = int( f.read().split('\n')[0] )
        except:
            quit_cpmd("Problem while reading xyz file: {}".format(XYZ))

    return natom


def get_energies(JOB, TITLE, grad_file, efile):
    """Create a file with energies from Turbomole for analysis purposes"""

    # extract energies from TURBOMOLE output file
    energies = []
    diag = 0.0
    with open(grad_file, 'r') as exstates:
        for line in exstates:
            if JOB in [RICC2, RICCS]:
                if 'Final CC2 energy' in line:
                    energies.append( float(line.split()[-2]) )
                elif 'Energy of reference wave function' in line and JOB == RICCS:
                    energies.append( float(line.split()[-1]) )
                elif 'Total energy of excited state:' in line:
                    energies.append( float(line.split()[-1]) )
                elif 'D1 diagnostic' in line:
                    diag = float(line.split()[-2]) 

            elif JOB == TDDFT:
                if 'Total energy:' in line:
                    energies.append( float(line.split()[-1]) )

            else:
                quit_cpmd('Undefined JOB for TURBOMOLE 1: {}'.format(JOB))

    # print energies to file
    with open(efile, 'a') as out:
        ml = '{:>80} Diagnostic: {:20.10e} '.format(TITLE, diag)
        for e in energies:
            ml += ' {:20.10e} '.format(e)
        out.write( ml+'\n' )


def get_forces(natom, TITLE, grad_file, forcefile):
    """Read forces from turbomole and write them to a file for CPMD"""

    # read forces from main output file
    with open(grad_file, 'r') as grad_f:
        grad = grad_f.read().split('\n')
    forces = read_forces(grad) 

    # print dictionary to forcefile
    with open(forcefile, 'w') as out:
        out.write( ' {:10d} \n'.format(natom) )
        out.write( TITLE + '\n' )
        for sym, force in forces:
            out.write( '{:>5}   {:20.10e}   {:20.10e}   {:20.10e}\n'.format(sym, *force) )


def read_forces(grad):
    """Read forces in Turbomole's format"""

    # forces is a list of tuple (key, value)
    # each key represents an atom
    # each value is a 3D vector containing the gradient
    forces = []
    for i, line in enumerate(grad):
        if ' ATOM ' in line:
            # typical line: ATOM      1 o           2 o           3 h
            # so n = (7-1)/2 = 3
            n = (len(line.split())-1)/2
            for j in range(n):
                # get symbol
                sym = line.split()[(j+1)*2].upper()
         
                # get force
                force = np.zeros(3)
                for k in range(3):
                    force[k] = -float( grad[i+k+1].split()[j+1].replace("D", "E") ) 
         
                # save info
                forces.append( (sym, force) ) 

    return forces


def prepare_and_run_turbo(CPDIR_IN, JOB, TITLE, GEO_XYZ, DEF_INP, FORCES_XYZ):
    """
    Obtain forces from a turbomole calculation
    
    This is the main routine here:
        1) Prepare calculation (directory, input files...)
        2) Run Turbomole calculation of nuclear forces
        3) Read Turbomole's output file and create new file with nuclear forces
    """

    global CPDIR
    CPDIR = CPDIR_IN

    # CREATE DIRECTORY WHERE TURBOMOLE WILL RUN
    TBDIR = os.path.join(CPDIR, 'TURBOMOLE')
    if os.path.isdir(TBDIR):
        # rm previous dir if it exists
        rmtree(TBDIR)
    # create empty directory
    os.mkdir(TBDIR)

    # PREPARE TURBOMOLE INPUT FILES
    os.system( 'x2t {} > {}/coord'.format(GEO_XYZ, TBDIR) )
    copy( DEF_INP, TBDIR )
    os.chdir(TBDIR)
    run_turbo('define', 'define.out', inp=DEF_INP)

    # WE MIGHT WANT TO MODIFY THE CONTROL FILE
    #update_control_file()
    

    # RUN TURBOMOLE
    if JOB == DFT:
        grad_file = 'dft_grad.out'
        run_turbo('dscf', 'dft.out')
        run_turbo('grad', grad_file)
    
    elif JOB == TDDFT:
        grad_file = 'tddft_grad.out'
        exc_file = 'tddft_exc.out'
        run_turbo('dscf', 'dft.out')
        run_turbo('escf', exc_file)
        run_turbo('egrad', grad_file)
    
    elif JOB == RICCS:
        grad_file = 'riccs.out'
        run_turbo('dscf', 'hf.out')
        run_turbo('ricc2', grad_file)
    
    elif JOB == RICC2:
        grad_file = 'ricc2.out'
        run_turbo('dscf', 'hf.out')
        run_turbo('ricc2', grad_file)
    
    elif JOB == RIMP2:
        grad_file = 'rimp2.out'
        run_turbo('dscf', 'hf.out')
        run_turbo('ricc2', grad_file)
    
    else:
        quit_cpmd('Undefined JOB for TURBOMOLE 3: {}'.format(JOB))
    
    
    # CREATE FORCES_XYZ FILE
    natom = get_natom_from_xyz(os.path.join(CPDIR, GEO_XYZ))
    forcefile = os.path.join(CPDIR, FORCES_XYZ)
    get_forces(natom, TITLE, grad_file, forcefile)
    
    # EVENTUALLY CREATE FILE WITH EXCITATION ENERGIES
    if JOB in [TDDFT, RICCS, RICC2]:
        energyfile = os.path.join(CPDIR, 'exc_energies')
        get_energies(JOB, TITLE, grad_file, energyfile)

    # get back to CPMD directory
    os.chdir(CPDIR)


if __name__ == '__main__':

    # READ INPUT ARGUMENTS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("JOB", 
            type=str, 
            choices=JOBS,
            help="String used to decide which programs to use when running Turbomole")
    parser.add_argument("TITLE", 
            type=str, 
            help="Title line describing the JOB (functional, basis set...)")
    parser.add_argument("--GEO_XYZ", 
            type=str, 
            default='geometry.xyz',
            help="Filename for the input geometry (from CPMD)")
    parser.add_argument("--DEF_INP", 
            type=str, 
            default='define.inp',
            help="Input file for 'define' used to create 'control' file")
    parser.add_argument("--FORCES_XYZ", 
            type=str, 
            default='forces.xyz',
            help="Output file containing the forces produced by Turbomole")
    args = parser.parse_args()

    # GET CURRENT WORKING DIRECTORY
    CPDIR = os.getcwd()

    # MAKE SURE ALL THE REQUIRED FILES ARE PRESENT
    # i.e. DEF_INP and GEO_XYZ
    files = [args.GEO_XYZ, args.DEF_INP]
    for f in files:
        if not os.path.exists(f):
            quit_cpmd('File: {} missing!'.format(f))
    
    # PREPARE EVERYTHING AND RUN TURBOMOLE
    prepare_and_run_turbo(CPDIR, args.JOB, args.TITLE, args.GEO_XYZ, args.DEF_INP, args.FORCES_XYZ)


