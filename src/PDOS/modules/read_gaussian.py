#!/usr/bin/env python
"""This file contains a set of functions usefull to extract some information from Gaussian output files."""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"


import numpy as np

from comp_chem_utils.periodic import element
from comp_chem_utils.utils import get_file_as_list
from comp_chem_utils.conversions import convert


class gauss_atom:
    """Contain basis information about an atom."""

    def __init__(self):
        self.label = "X"
        self.index = 0
        self.basis = {}

    def output(self):
        print("label = {}".format(self.label))
        print("index = {}".format(self.index))
        print("basis = {}".format(self.basis))


def get_gaussian_info(filename):
    """Extract relevant information from a Gaussian output file.
    
    Goes through a gaussian output file looking for the information
    required to calculate PDOS.

    Args:
        filename (str): Name of the Gaussian output file, it may
            also include the path to the file.

    Returns:
        nocc (int): Number of occupied orbitals.

        nbas (int): Number of basis functions (atomic orbitals).

        overlap (np.array): Squared AO overlap matrix. 

        epsilon (list): Molecular orbital energies.

        coef (np.array): Matrix of the molecular orbital coefficients.
    """

    # Store the whole gaussian file in output
    # ---------------------------------------
    output = get_file_as_list(filename)
    
    
    # Get Number of basis functions
    # -----------------------------
    nline = get_line(output,"   NBasis = ")
    nbas = int(output[nline].split()[2])

    # get number of occupied orbitals (assume closed shell)
    # -----------------------------------------------------
    nline = get_line(output,"alpha electrons")
    nocc = int(output[nline].split()[0])

    print("Number of basis functions   = {}".format(nbas))
    print("Number of occupied orbitals = {}\n".format(nocc))
    print("Reading Overlap and MO transformation matrices...")
    
    
    # Read overlap matrix from gaussian output
    # ----------------------------------------
    nline = get_line(output," *** Overlap ***") + 1
    overlap = read_triangular_matrix(nline,nbas,output)
    
    
    # Read Orbital energies and Orbital coefficients
    # ----------------------------------------------
    nline = get_line(output,"Molecular Orbital Coefficients") + 1
    
    # The MO coefficients are stored in the coef matrix.
    # we assume the number of molecular orbitals equal to
    # the number of atomic basis functions.
    coef=np.zeros((nbas,nbas))
    epsilon=[] # orbital energies
    
    iline = nline + 1 # skip column headers
    # first read complete blocks
    for iblock in range(nbas/5): 
        iline = read_orbitals_and_coefs(iblock,nbas,iline,5,epsilon,coef,output)
    
    # If nbas is not a multiple of five the next routine 
    # reads the last data of the matrix
    rest=nbas%5 
    if (rest!=0):
        iblock = nbas/5
        iline = read_orbitals_and_coefs(iblock,nbas,iline,rest,epsilon,coef,output)

    print("Gaussian output file parsed")
    return nocc,nbas,overlap,epsilon,coef


def get_line(output,string):
    """ Return the line number in 'output' containing the first occurence of 'string' """ 

    for i,line in enumerate(output):

        if (string in line):
            return i


def read_triangular_matrix(start,dim,output):

    # The GAUSSIAN's matrix are made by N/5 block of 5 columns and N lines each one.
    # (caution: there is just the half matrix in the GAUSSIAN file so we have to rebuild the other part)
    
    mat = np.zeros((dim,dim))
    iline = start
    # Loop over each block of 5 columns
    for iblock in range(dim/5):  
    
        iline += 1 # increment index in output (empty line between blocks)
    
        # Reads each line of a block (the length of the block decrease with the index dim)
        iline = read_triangular_block(iblock,dim,mat,iline,output)

    # If dim is not a multiple of five the next routine reads the last data of the matrix
    rest=dim%5 
    if (rest!=0):
    
        iline += 1 # increment index in output (empty line between blocks)
    
        iblock = dim/5
        iline = read_triangular_block(iblock,dim,mat,iline,output)

    return mat


def read_triangular_block(iblock,dim,mat,iline,output):

    start_ind = 5*iblock # start row and column index depending on block
    for irow in range(start_ind,dim): 

        # read line and replace double precision fortran notation with standard exponent
        test=output[iline].replace("D","e").split()[1:]

        # read the (maximum) 5 numbers in a row
        for i,f in enumerate(test):
            # i goes from 0 to maximum 4, we add start_ind depending on block
            icol = start_ind + i
            mat[irow,icol] = float(f)
            mat[icol,irow] = float(f) # symmetrize
    
        iline += 1 # increment index in output (read next line)

    return iline


def read_orbitals_and_coefs(iblock,nbas,iline,ncol,epsilon,coef,output):

    # read orbital energies
    # ---------------------
    iline += 1 
    test = output[iline]
    for i in range (ncol):
        # Conversion: [Hartree] --> [eV]
        e = float(test.split()[i+2])
        epsilon.append( convert( e, 'ENERGY: a.u.', 'ENERGY: eV') )
    
    # Read each line of a block  
    # -------------------------
    start_ind = iblock*5
    for irow in range(0,nbas): 
        iline += 1 # go to next line
        test = output[iline].split()[-ncol:] # get ncol last values in line

        for i,f in enumerate(test):
            icol = start_ind + i
            coef[irow,icol] = float(f)
    
    # Reads the empty lines  
    iline += 2 # skip lines between blocks
    return iline



def read_basis(iline, nbas, output):

    basis = {}
    # read first basis set for given atom
    # the gaussian line look like:
    # 1 1   Ru 1S           .00000    .00000    .00000    .00000    .00000
    #
    # in that case label = 1S and we store the index 1.
    label = output[iline].split()[3]
    basis[label] = int(output[iline].split()[0])

    iline += 1
    if (iline >= nbas):
        return iline, basis

    # now the gaussian line should look like:
    # 2        2S           .00000    .00000    .00000    .00000    .00000
    #
    # check if next line is a next atom
    next_atom = element(output[iline].split()[2])
    while(not next_atom):

        label = output[iline].split()[1]
        basis[label] = int(output[iline].split()[0])
        iline += 1
        if (iline >= nbas):
            break
        next_atom = element(output[iline].split()[2])

    # the atom's basis has been read, we are done
    return iline, basis


def read_atom_and_basis(gaussian_file):

    # Store the whole gaussian file in output
    # ---------------------------------------
    myfile = open(gaussian_file,"r")
    output = []
    for line in myfile:
        output.append(line)
    myfile.close()
        
    
    # Parse Gaussian output file
    # make a dictionary of atoms and basis functions
    # ----------------------------------------------
    nline = get_line(output,"   NBasis = ")
    nbas = int(output[nline].split()[2])
    nline = get_line(output,"EIGENVALUES") + 1
    output = output[nline:nline + nbas]
    gauss_list = []
    
    # read all lines in output
    iline = 0
    while (iline < nbas):
    
        # read output line
        line = output[iline]
    
        # Check for atomic symbol in second position
        if element(line.split()[2]):
            atom = gauss_atom()
            atom.index = int(line.split()[1])
            atom.label = line.split()[2]
        else:
            print("Something wrong while reading gaussian file")
            exit()
    
        # read basis set for atom
        iline, atom.basis = read_basis(iline, nbas, output)
        gauss_list.append(atom)

    return gauss_list


def read_boolean(string):
    """ read boolean from users (yes/no question where yes is default)"""

    # read capitalized input
    inp = raw_input(string).upper() 

    # set logical based on input
    if (inp == "N" or inp == "NO"):
        logical = False
    else:
        logical = True

    return logical



