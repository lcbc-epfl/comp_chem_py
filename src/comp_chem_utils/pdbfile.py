#!/usr/bin/env python
"""Suppose to deal with pdb file format"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import numpy as np

from comp_chem_utils.periodic import element
from comp_chem_utils.utils import *


# remoteness indicator code sequence
remote_code = ['A','B','G','D','E','Z','H']

class pdb_file(object):
    """for now we assume a very simple structure:
    CRYST1
    ATOM
    ATOM
    ...
    END
    ATOM
    ATOM
    ...
    """

    def __init__(self):
        self.name = ''
        self.first_line = ''
        self.structures = []

    def nstr(self):
        return len(self.structures)

    def change_atom(self, new_name, new_symb, old_name, old_symb):
        """change the atom name and symbol for a specific type of atom in the pdbfile"""
        new_str = []
        for struc in self.structures:
            struc.change_atom(new_name, new_symb, old_name, old_symb)
            new_str.append(struc)
        self.structures = new_str

    def change_bond_length(self, atom_fix, atom_mob, R):
        """move atom_mob such that the distance to atom_fix is R"""
        new_str = []
        for struc in self.structures:
            struc.change_bond_length(atom_fix, atom_mob, R)
            new_str.append(struc)
        self.structures = new_str

    def read(self, filen):
        """Read pdb file"""
        lines = get_file_as_list(filen)

        self.name = filen.split("/")[-1][:-4]
        self.first_line = lines[0].rstrip()

        new_structure = True
        iline = 1
        while new_structure:
            # check if current line is an atom:
            if iline >= len(lines):
                # end of file
                new_structure = False

            elif not lines[iline]:
                # end of file
                new_structure = False

            elif lines[iline][0:4] == 'ATOM':

                # read structure
                struc = structure()
                name = self.name + '_' + str(self.nstr()+1)
                iline = struc.read(lines, iline, name)
                self.structures.append(struc)

            elif lines[iline][0:3] == 'END':
                iline +=1
                new_structure = True

    def output(self):

        lines = []
        lines.append(self.first_line)
        for struc in self.structures:
            lines.extend(struc.output())
            lines.append('END')

        return lines

    def out_pdbfile(self, fn=''):

        if not fn:
            fn=self.name

        with open(fn,'w') as myfile:
            for line in self.output():
                myfile.write(line + '\n')


    def keep_only(self, list_of_res):
        """delete all residues from the structures that dont belong to the list"""
        newstruc = []
        for struc in self.structures:
            struc.keep_only(list_of_res)
            newstruc.append(struc)

        self.structures = newstruc


class structure(object):
    """for now we assume a structure is just a list of residues"""

    def __init__(self):
        self.name = ''
        self.residues = []

    def nres(self):
        return len(self.residues)

    def change_atom(self, new_name, new_symb, old_name, old_symb):
        """change the atom name and symbol for a specific type of atom in the structure"""
        new_res = []
        for res in self.residues:
            res.change_atom(new_name, new_symb, old_name, old_symb)
            new_res.append(res)
        self.residues = new_res

    def change_bond_length(self, atom_fix, atom_mob, R):
        """move atom_mob such that the distance to atom_fix is R"""
        new_res = []
        for res in self.residues:
            res.change_bond_length(atom_fix, atom_mob, R)
            new_res.append(res)
        self.residues = new_res

    def read(self, lines, iline, name):

        self.name = name
        new_residue = True
        while new_residue:
            
            # check if current line is an atom:
            if iline >= len(lines):
                # end of file
                new_residue = False

            elif lines[iline][0:4] == 'ATOM':

                # read residue
                res = residue()
                iline = res.read(lines, iline)
                 
                # add residue to structure
                self.residues.append(res)

            else:
                new_residue = False

        return iline

    def output(self):

        lines = []
        for res in self.residues:
            lines.extend(res.output())

        return lines

    def keep_only(self, list_of_res):
        """delete all residues from the structure that dont belong to the list"""
        newres = []
        for res in self.residues:
            if res.number in list_of_res:
                newres.append(res)

        self.residues = newres

    def get_list_of_res(self):
        """return list of residues in structure"""
        list_of_res = []
        for res in self.residues:
            list_of_res.append(res.number)

        return list_of_res
    
    def export_xyz(self):
        """export xyz data as a 4 columns and natoms lines table"""
        table = []
        for res in self.residues:
            for myat in res.atoms:
                line = []
                line.append(myat.symbol.strip())
                line.extend(myat.coord.tolist())
                table.append(line)

        return table


class residue(object):

    def __init__(self):
        self.name = ''
        self.number = 0
        self.atoms = []
        self.nelec = 0

    def add_atom(self,atom):
    
        success = False # Has the atom been added succesfuly?

        # if list of atoms is empty the residue info is coming from atom
        if not self.atoms:
            self.name = atom.residue_name
            self.number = atom.residue_number
            self.atoms.append(atom)
            success = True
            self.nelec += element(atom.symbol.strip()).atomic 

        # else we check for consistency
        elif (self.name == atom.residue_name) and (self.number == atom.residue_number):
            self.atoms.append(atom)
            success = True
            self.nelec += element(atom.symbol.strip()).atomic 
      
        return success

    def change_atom(self, new_name, new_symb, old_name, old_symb):
        """change the atom name and symbol for a specific atom in the residue"""

        new_atoms = []
        for myat in self.atoms:
            if (myat.name.out() == old_name) and (myat.symbol == old_symb):
                myat.name.read(new_name)
                myat.symbol = new_symb

            new_atoms.append(myat)

        self.atoms = new_atoms
        # check and correct for atom names
        self.rename_h_atoms()

    def rename_h_atoms(self):

        for rc in remote_code:
            ibr = 1
            new_atoms = []
            for myat in self.atoms:
                if (myat.name.remote==rc) and (myat.name.symbol==' H'):
                    myat.name.branch=str(ibr)
                    ibr += 1
                new_atoms.append(myat)
            self.atoms = new_atoms


    def change_bond_length(self, atom_fix, atom_mob, R):
        """move atom_mob such that the distance to atom_fix is R"""

        # get coordinates for the two atoms
        new_atoms = []
        ifix = -1
        imob = -1
        for i,myat in enumerate(self.atoms):
            if myat.name.out() == atom_fix:
                ifix = i
            if myat.name.out() == atom_mob:
                imob = i

            new_atoms.append(myat)

        if (ifix == -1) or (imob == -1):
            print("ERROR(residue: change_bond_length): atoms not found!")
            return

        # get new coordinates
        new_atoms[imob].coord = change_vector_norm(new_atoms[ifix].coord, new_atoms[imob].coord, R)

        # update coordinates for mobile atom line
        self.atoms = new_atoms


    def read(self, lines, iline):

        new_atom = True
        while new_atom:

            # check if current line is an atom:
            if iline >= len(lines):
                # end of file
                new_atom = False

            elif lines[iline][0:4] == 'ATOM':

                # read atom line
                myat = atom()
                myat.read(lines[iline])
                 
                # try adding atom to residue
                new_atom = self.add_atom(myat)
                if new_atom:
                    iline += 1

            else:
                new_atom = False

        return iline


    def output(self):

        lines = []
        for myat in self.atoms:
            lines.append(myat.output())

        return lines


class atom(object):

    def __init__(self):
        self.serial_number = 0
        self.name = atom_name()
        self.loc_ind = ''
        self.residue_name = ''
        self.chain_id = ''
        self.residue_number = 0
        self.code_insert_residue = ''
        self.coord = np.zeros(3)
        self.occupancy = 0.0
        self.temp_factor = 0.0
        self.segment_id = ''
        self.symbol = ''

    def read(self,line):
        """read line from pdb file that starts with ATOM and init atom type"""
        self.serial_number = int(line[6:11])
        self.name.read(line[12:16])
        self.loc_ind = line[16]
        self.residue_name = line[17:20]
        self.chain_id = line[21]
        self.residue_number = get_res_num(line)
        self.code_insert_residue = line[26]
        for i,x in enumerate( line[30:54].split() ):
            self.coord[i] = float(x)
        self.occupancy = float(line[54:60])
        self.temp_factor = float(line[60:66])
        self.segment_id = line[72:76]
        self.symbol = line[76:78]
        self.line = line


    def update(self):
        """reconstruct self.line from data in atom"""
        self.line = 'ATOM  '
        self.line += "{:5d}".format(self.serial_number) + ' '
        self.line += self.name.out()
        self.line += self.loc_ind 
        self.line += self.residue_name + ' '
        self.line += self.chain_id 
        self.line += "{:4d}".format(self.residue_number)
        self.line += self.code_insert_residue + '   '
        self.line += "{:8.3f}{:8.3f}{:8.3f}".format(*self.coord.tolist())
        self.line += "{:6.2f}".format(self.occupancy)
        self.line += "{:6.2f}".format(self.temp_factor) + '      '
        self.line += self.segment_id
        self.line += self.symbol

    def output(self):
        """return string to be printed as a line in a pdb file"""
        self.update()
        return self.line

class atom_name(object):

    def __init__(self):
        self.symbol = '  '
        self.remote = ' ' # remoteness indicator code
        self.branch = ' ' # branch indicator code

    def read(self, name):
        self.symbol = name[0:2]
        self.remote = name[2]
        self.branch = name[3]

    def out(self):
        return self.symbol + self.remote + self.branch


def get_res_num(line):
    """get residue number from line starting with ATOM"""
    return int(line[22:26])


