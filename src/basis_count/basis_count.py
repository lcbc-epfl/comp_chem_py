#!/usr/bin/env python
"""Calculate number of basis function from molecular structure and basis set"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import sys
import argparse

from ebsel.src.EMSL_local import EMSL_local
from ebsel.src.conversion import Converter

from comp_chem_utils.molecule_data import mol_file


shells = ["S", "P", "D", "F", "G", "H", "I", "K", "L", "M"]
conv = Converter().parse_multi_g94


def degenerate_from_symbol(symbol):
    """
    return the number of degenerate levels for a given orbital angular momentum L
    
    L is given by its letter symbol: S, P, D...
    """
    if symbol=="SP":
        symbol="P"
    return 2 * shells.index(symbol) + 1


def nbasis(bse):
    """Calculate number of basis functions from a BasisSetEntry object"""
    nbasis = 0
    for symb, nfunc in bse.functions_per_shell.items():
        nbasis += nfunc * degenerate_from_symbol(symb)

    return nbasis


def nbasis_from_molbas_data(molbas_data):
    """returns the number of basis function for a molecular structure"""

    # get basis set database in gaussian format
    bs_g94 = EMSL_local(fmt="g94")
   
    ntot = 0 
    new_molbas_data = {}
    for bas, el_dic in molbas_data.items():
   
        # get basis set for specific elements as a string
        basis = "\n".join(bs_g94.get_basis( bas, elements=el_dic.keys() ))
   
        # convert basis set to BasisSetEntry object
        bse = conv(basis, "db/Gaussian94.db")
   
        new_el_dic = {}
        for el_bas in bse:
            # get number of basis function for a single specific atom
            nbas = nbasis( el_bas )
            # multiply by occurance of that specific atom in molbas_data
            oc = el_dic[el_bas.symbol]
            ntot += nbas * oc

            # create new dictionary with basis set info
            new_el_dic[el_bas.symbol] = [oc, nbas]

        new_molbas_data[bas] = new_el_dic
   
    return new_molbas_data, ntot


def print_molbas_data(molbas_data, ntot=0, basis=False):
    print "\n\nPRINT MOLECULAR DATA:"
    print "====================="
    for bas, el_dic in molbas_data.items():
        print "\n\nBasis set: {}".format(bas)
        if basis:
            print "    ------------------------------------"
            print "      Symbol   Num.      Nbas     Ntot"
            print "    ------------------------------------"
            for el, l in el_dic.items():
                print "      {:3}:    {:5d}    {:6d}   {:6d}".format(el, l[0], l[1], l[0]*l[1])
        else:
            print "    ------------------"
            print "      Symbol   Num."
            print "    ------------------"
            for el, oc in el_dic.items():
                print "      {:3}:    {:5d}".format(el, oc)
    if basis:
        print "\n\nTOTAL NUMBER OF BASIS FUNC: {:6d}".format(ntot)


def molbas_data_from_user():
    molbas_data = {}
    while(True):
        basis=raw_input("\nBasis set, e.g. cc-pVDZ (q=quit): ")
        if not basis or basis=='q':
            break
        el_dic={}
        while(True):
            atype=raw_input("Atom type, e.g. He (q=quit): ")
            if not atype or atype=='q':
                break
            Natype=int(raw_input("Number of atoms of type {} with basis {}: ".format(atype, basis)))
            el_dic[atype] = Natype

        molbas_data[basis] = el_dic

    return molbas_data


def molbas_data_from_dalton(fn):
    """Build molecular data from a MOLECULE.INP file (dalton format)"""

    # create mol_file object from MOLECULE.INP:
    mymol = mol_file()
    mymol.read_mol(fn)

    # create molbas_data dictionary from mymol info
    molbas_data = {}
    if mymol.atom_basis:
        # several basis set for different atom types
        for atype in mymol.atom_types:

            # get potentially existing dictionary
            el_dic = molbas_data.get(atype.basis)
            if not el_dic:
                # create empty dic. 
                el_dic = {}

            # add occurence of atom type
            el_dic[atype.symb] = atype.natoms

            # add (updated) dictionary to molbas_data
            molbas_data[atype.basis] = el_dic

    else:
        # only one basis set:
        el_dic = {}
        for atype in mymol.atom_types:
            # add occurence for each atom type
            el_dic[atype.symb] = atype.natoms

        molbas_data[mymol.basis] = el_dic

    return molbas_data



if __name__ == "__main__":

    # read optional arguments from user
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--mol", type=str, help="Molecular structure in Dalton format (MOLECULE.INP)")
    args = parser.parse_args()

    if not args.mol:
        # tell user that he might also give a MOLECULE.INP as argument
        parser.print_help()

        # no MOLECULE.INP file get mobas_data from interface with user
        molbas_data = molbas_data_from_user()

    else:
        # read mobas_data from MOLECULE.INP file
        molbas_data = molbas_data_from_dalton(args.mol)

    # Test for water in the aug-cc-pVDZ' basis 
    #molbas_data = {"cc-pVDZ" : {"H":2}, "aug-cc-pVDZ" : {"O":1} }

    new_molbas_data, ntot = nbasis_from_molbas_data(molbas_data)
    print_molbas_data(new_molbas_data, ntot=ntot, basis=True)

