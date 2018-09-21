#!/usr/bin/env python
"""This script can be used to generate the group_of_AOs file required by PDOS.""" 


__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import sys
import argparse

from PDOS.modules.read_gaussian import read_atom_and_basis, read_boolean


def get_list_of_basis(atom, pos_bas):
    """Read set of basis functions from user
    
    It makes sure to correct for upper/lower cases 
    and generalizations such as P = 1PX, 1PY, 1PZ, 2PX, 2PY..."""

    print "Insert list of atomic orbitals corresponding to atom ", atom
    print "that you which to consider forthe group of AOs such as, 1S, 2PZ..."
    inp_basis = raw_input("You can also just choose all 's' and 'p' orbitals by wrinting, S, P\n")
    # split line by comma and remove blanks
    grp_basis = [l.strip().upper() for l in inp_basis.split(",")]

    # make detailed list in case user specified partial orbitals such as S or P
    list_of_basis = []
    for bas in grp_basis:
        if (bas in pos_bas):
            # add it
            list_of_basis.append(bas)
        else:
            # check if part of it is
            for pb in pos_bas:
                if (bas in pb):
                    list_of_basis.append(pb)

    return list_of_basis


# for compatibility with autodoc in sphinx
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gaussian_file", type=str, help="Gaussian log (output) file containing relevant info")
    args = parser.parse_args()
    
    gauss_list = read_atom_and_basis(args.gaussian_file)
    
    # count and print number of different atoms
    # -----------------------------------------
    atoms = []
    prev = ""
    for atom in gauss_list:
    
        if (atom.label != prev):
            atoms.append(atom.label)
            prev = atom.label
    
    print "Found ", len(atoms), " type(s) of atom(s)"
    print atoms,"\n"
    
    
    
    print "Let's build the groups of AOs"
    print "-----------------------------\n"
    groups = {}
    get_new_group = True
    while(get_new_group):
        list_of_AOs = []
        list_of_atoms = []
    
        # Make groups of atomic orbitals by interacting with user
        # -------------------------------------------------------
        print "Insert the list of atoms to consider for the group of AOs as a comma separated list such as,"
        inp_atoms = raw_input("H, He, Li \n")
        
        # split line by comma and remove blanks
        list_of_atoms = [l.strip() for l in inp_atoms.split(",")]
        
        # For each atom in the group add user defined AOs
        # -----------------------------------------------
        for atom in list_of_atoms:
    
            # sanity check
            if (atom not in atoms):
                print "The atoms you have chosen are not part of the Gaussian output file..."
                exit()
        
            # find atom in dictionary to print its basis functions 
            # it assumes that all atoms of a given type have the 
            # same set of basis functions
            for gatom in gauss_list:
                found_it = False
                if (atom == gatom.label):
    
                    print "Possible basis functions for ", atom
                    pos_bas = sorted(gatom.basis.keys())
                    print pos_bas, "\n"
                    found_it = True
    
                if (found_it):
                    # exit loop when the basis functions of atom have been printed
                    break
    
            # get_ list of basis functions to include from user
            list_of_basis = get_list_of_basis(atom, pos_bas)
    
            # add basis function to current group
            grp_label = raw_input("Insert the label to associate to this group\n")
            groups[grp_label] = []
    
            # for each atom of he current type add the basis set indices to the group
            for gatom in gauss_list:
    
                if (atom == gatom.label):
    
                    # loop over basis functions of that atom
                    for bas_label in gatom.basis:
                        if (bas_label in list_of_basis):
                            groups[grp_label].append(gatom.basis[bas_label])
    
    
        get_new_group = read_boolean("Create another group? y/n [y]:\n")
    
    
    # print group to file
    # -------------------
    output = open("group_of_AOs.inp","w")
    for grp_label in groups:
        output.write(grp_label+"\n")
    
        bas = groups[grp_label]
        output.write(",".join([str(b) for b in bas]))
        output.write("\n")
    
    output.close()

