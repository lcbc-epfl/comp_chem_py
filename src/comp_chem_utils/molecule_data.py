#!/usr/bin/env python
"""Set of classes to deal with molecular data and formats"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"


from comp_chem_utils.conversions import BOHR_TO_ANG
from comp_chem_utils.mysql_tables import mol_info_table, mol_xyz_table
from comp_chem_utils.periodic import element

hostname = 'localhost'
username = 'root'
password = ''
database = 'molecules'


class mol_data(object):
    """ All information regarding a given molecule. """

    # INITIALIZTION FUNCTIONS
    # -----------------------
    def __init__(self):
        self.idd = 0
        self.name = ""
        self.chem_name = ""
        self.note = ""
        self.tab_name = ""
        self.charge = 0
        self.natoms = 0
        self.natom_types = 0
        self.atom_types = []

    def init_from_mol(self, mol_file):
        """ update mol_data with data in mol (DALTON) format"""
        self.name = mol_file.fname.split("/")[-1][:-4]
        self.note = mol_file.title1+" "+mol_file.title2
        self.charge = mol_file.charge
        self.natoms = mol_file.natoms
        self.natom_types = mol_file.natom_types
        self.atom_types = mol_file.to_angstrom().atom_types
        self.get_chem_name()
        return self.check() # check validity of mol_data"

    def init_from_xyz(self, xyz_file):
        """ update mol_data with data in xyz format"""
        self.name = xyz_file.fname.split("/")[-1][:-4]
        self.note = xyz_file.title
        #print("WARNING: init_from_xyz setting default charge to zero")
        self.charge = 0
        self.natoms = xyz_file.natoms
        self.natom_types = xyz_file.natom_types
        self.atom_types = xyz_file.atom_types
        self.get_chem_name()
        return self.check() # check validity of mol_data"

    def init_from_db(self, cur, idd):
        """ update mol_data with data from database"""

        # get main information from mol_info table
        main_tab = mol_info_table()
        cur.execute( main_tab.get_row(idd) )
        for (idd, name, chem_name, note, charge, natoms, natom_types) in cur:
            self.idd = idd
            self.name = name
            self.chem_name = chem_name
            self.note = note
            self.charge = charge
            self.natoms = natoms
            self.natom_types = natom_types

        # get xyz data from molecule specific table
        xyz_tab = mol_xyz_table( idd )
        cur.execute( xyz_tab.get_table() )
        old_charge = 0
        for (idd, symb, charge, x, y, z) in cur:

            if charge!=old_charge:
                if (old_charge>0):
                    # add old atype to atom_types list
                    self.atom_types.append(atype)
                # new atom type
                atype = atom_type()
                old_charge = charge

            atype.charge = charge
            atype.symb = symb
            atype.xvals.append(x)
            atype.yvals.append(y)
            atype.zvals.append(z)
            atype.natoms += 1

        # add last atype to atom_types list
        self.atom_types.append(atype)
        self.get_chem_name()
        return self.check() # check validity of mol_data"
            


    # OUTPUT FUNCTIONS
    # ----------------
    def output(self):
        """ return list of strings containing mol_data"""
        lines = []
        lines.append("ID         : "+str(self.idd))
        lines.append("NAME       : "+self.name)
        lines.append("CHEM NAME  : "+self.chem_name)
        lines.append("NOTE       : "+self.note)
        lines.append("TABLE NAME : "+self.tab_name)
        lines.append("CHARGE     : "+str(self.charge))
        lines.append("NATOMS     : "+str(self.natoms))
        lines.append("NATOM TYPES: "+str(self.natom_types))
        # print atom types
        for atype in self.atom_types:
            lines.extend(atype.output())

        return lines

    def out_to_mol(self):
        """ generate mol_file data from mol_data """
        mol = mol_file()
        mol.fname = self.name.replace(' ','_')+".mol"
        mol.charge = self.charge
        mol.natoms = self.natoms
        mol.bohr = False
        if len(self.note.split())>8:
            mol.title1 = self.note[:6].replace('\n',' ')
            mol.title2 = self.note[6:].replace('\n',' ')
        else:
            mol.title1 = self.note.replace('\n',' ')

        mol.natom_types = self.natom_types
        mol.atom_types = self.atom_types

        return mol

    def out_to_xyz(self):
        """ generate xyz_file data from mol_data """
        xyz = xyz_file()
        xyz.fname = self.name.replace(' ','_')+".xyz"
        xyz.natoms = self.natoms
        xyz.title = self.note.replace('\n',' ')
        xyz.natom_types = self.natom_types
        xyz.atom_types = self.atom_types

        return xyz

    def out_to_db(self, cur):
        """ add mol_data as a new entry to the molecules database"""
        if self.check():  # check validity of mol_data before saving to database"
            return self.check()

        main_tab = mol_info_table()
        cur.execute(main_tab.add_row(), (self.name, self.chem_name, self.note, self.charge, self.natoms, self.natom_types) )
        # get mol_id for mol_info
        self.idd = cur.lastrowid

        # create xyz table:
        xyz_tab = mol_xyz_table( self.idd )
        cur.execute( xyz_tab.create_table() )
        for atype in self.atom_types:
            atype.add_db_entry(cur, xyz_tab.add_row())

    def get_chem_name(self):
        """ get/set chemical formula based on mol_data"""
        self.chem_name = ""
        for atype in self.atom_types:
            self.chem_name += "{0.symb}_{0.natoms} ".format(atype)

        # add total charge
        self.chem_name += " ({:+d})".format(self.charge)

    def check(self):
        """ check consistency of mol_data"""
        # sum of atoms in atom_types must be equal to natoms
        nat = 0
        for atype in self.atom_types:
            nat += atype.natoms

        if nat!=self.natoms:
            print("ERROR: Inconsistency in mol_data (natoms)".format(nat,self.natoms))
            return "Inconsistency in mol_data (natoms)"



#======================================================================================
class atom_type(object):
    """ atom type as used in mol files (DALTON format)"""

    def __init__(self):
        self.symb = "X"
        self.charge = 0
        self.natoms = 0
        self.bohr = False
        self.basis = ""
        self.aux_basis = ""
        self.xvals = []
        self.yvals = []
        self.zvals = []

    def to_angstrom(self):
        """ convert coordinates from bohr to angstrom"""
        if self.bohr:
            self.xvals = [x*BOHR_TO_ANG for x in self.xvals]
            self.yvals = [y*BOHR_TO_ANG for y in self.yvals]
            self.zvals = [z*BOHR_TO_ANG for z in self.zvals]

        self.bohr = False
        return self

    def read_from_mol(self, lines, iline, bohr):
        self.bohr = bohr
        line = lines[iline].split()
        for elmt in line:
            el = elmt.split("=")
            test = el[0].lower()
            if test == "atoms":
                self.natoms = int(float(el[1]))
            elif test == "charge":
                self.charge =int(float(el[1]))
            elif test == "basis" or el[0] == "bas":
                self.basis = el[1]
            elif test == "aux":
                self.aux_basis = el[1]

        self.symb = element(self.charge).symbol
        # read all atoms in type
        for nato in xrange(self.natoms):
            iline += 1
            line = lines[iline].split()
            self.xvals.append(float(line[1]))
            self.yvals.append(float(line[2]))
            self.zvals.append(float(line[3]))

    def output(self,mol=False):
        """ return list of strings as in mol file (Atomtypes section) """
        lines = []
        if mol:
            # when call for mol file as in DALTON format
            line = "Charge="+str(self.charge)+" Atoms="+str(self.natoms)
            if self.basis!='':
                line += " Basis="+self.basis
                if self.aux_basis!='':
                    line += " Aux="+self.aux_basis
            lines.append(line)

        for x,y,z in zip(self.xvals, self.yvals, self.zvals):
            line = "{:<4} {:20.10f} {:20.10f} {:20.10f}".format( self.symb, x, y, z )
            lines.append( line )

        return lines

    def add_db_entry(self, cur, add_mol_xyz):
        """ add xyz data as a new entry to the molecules database"""

        for x, y, z in zip(self.xvals, self.yvals, self.zvals):
            cur.execute(add_mol_xyz, (self.symb, self.charge, x, y, z) )

    def set_from_GUI(self, atype_sec):
        """ set atom_type information from GUI """
        self.symb = atype_sec.symb
        self.charge = atype_sec.charge
        self.natoms = atype_sec.natoms

        for x, y, z in  zip(atype_sec.xvals, atype_sec.yvals, atype_sec.zvals):
            self.xvals.append( x.get() )
            self.yvals.append( y.get() )
            self.zvals.append( z.get() )



#======================================================================================
class mol_file(object):
    """ structure and data of mol files (DALTON format)"""

    def __init__(self):
        self.fname = ""
        self.atom_basis = False
        self.basis = "cc-pVDZ"
        self.aux_basis = "cc-pVDZ-RI"
        self.title1 = ""
        self.title2 = ""
        self.natoms = 0
        self.natom_types = 0
        self.symmetry = False
        self.bohr = True
        self.charge = 0
        self.atom_types = []

    def check(self):
        """ check consistency of mol_file"""
        # sum of atoms in atom_types must be equal to natoms
        nat = 0
        for atype in self.atom_types:
            nat += atype.natoms

        if nat!=self.natoms:
            print("ERROR: Inconsistency in mol_file data (natoms)".format(nat,self.natoms))
            return "Inconsistency in mol_file data (natoms)"

    def read_mol(self, filename):
        """ update mol_file data from mol file"""
        self.fname = filename
        lines = []
        with open(self.fname, "r") as xyzfile:
            for line in xyzfile:
                lines.append(line.rstrip('\n'))

        self.atom_basis = lines[0].strip().lower() == "atombasis"
        if self.atom_basis:
            offset = 1
        else:
            offset = 2
            # read basis
            line = lines[1].split()
            self.basis = line[0]
            if len(line)>1:
                # read aux basis
                if line[1][0:3].lower() == "aux":
                    self.aux_basis = line[1][4:]

        # read title lines:
        self.title1 = lines[offset]
        self.title2 = lines[offset+1]

        # read main line:
        line = lines[offset+2].split()
        for elmt in line:
            el = elmt.split("=")
            test = el[0].lower()
            if test == "atomtypes":
                self.natom_types=int(float(el[1]))
            elif test == "charge":
                self.charge=int(float(el[1]))
            elif test == "nosymmetry":
                self.symmetry = False
            elif test == "angstrom":
                self.bohr = False 

        # read all atom types
        iline = offset + 2
        self.atom_types = [] 
        for natp in xrange(self.natom_types):
            iline += 1
            new_atype = atom_type()
            new_atype.read_from_mol(lines, iline, self.bohr)
            self.atom_types.append( new_atype )
            self.natoms += self.atom_types[-1].natoms
            iline += self.atom_types[-1].natoms

        self.natom_types = len(self.atom_types)
        return self.check() # check validity of data


    def to_angstrom(self):
        """ convert coordinates from bohr to angstrom"""
        if self.bohr:
            # convert data
            ang_types = []
            for atype in self.atom_types:
                ang_types.append(atype.to_angstrom())
            self.atom_types = ang_types

        self.bohr = False
        return self
            

    def output(self):
        """ return list of strings as in mol file """
        lines = []
        # print basis info
        if (self.atom_basis):
            lines.append("ATOMBASIS")
        else:
            lines.append("BASIS")
            line = self.basis
            if (self.aux_basis!=''):
                line += " Aux="+self.aux_basis
            lines.append(line)

        # print title lines
        lines.append(self.title1)
        lines.append(self.title2)

        # print main line
        line = "Atomtypes="+str(self.natom_types)
        line += " Charge="+str(self.charge)
        if not self.symmetry:
            line += " Nosymmetry"
        if not self.bohr:
            line += " Angstrom"
        lines.append(line)

        # print atom types
        for atype in self.atom_types:
            lines.extend(atype.output(True))

        return lines

    def out_to_file(self, fname=''):
        if self.check():
            return self.check() # check validity of data"

        if not fname:
            fname=self.fname

        with open(fname, 'w') as outfile:
            for line in self.output():
                outfile.write(line + "\n")


        
#======================================================================================
class xyz_file(object):
    """ structure and data of xyz files"""

    def __init__(self):
        self.fname = ""
        self.natoms = 0
        self.title = ""
        self.natom_types = 0
        self.atom_types = []

    def check(self):
        """ check consistency of xyz_file"""
        # sum of atoms in atom_types must be equal to natoms
        nat = 0
        for atype in self.atom_types:
            nat += atype.natoms

        if nat!=self.natoms:
            print("ERROR: Inconsistency in xyz_file data (natoms)".format(nat,self.natoms))
            return "Inconsistency in xyz_file data (natoms)"

    def read_xyz(self, filename):
        """ update xyz_file data from xyz file"""
        self.fname = filename
        lines = []
        with open(self.fname, "r") as xyzfile:
            for line in xyzfile:
                lines.append(line.rstrip('\n'))
        
        # Get number of atoms:
        try:
            self.natoms = int(lines[0])
        except:
            return 'file: '+self.fname+' does not fit the xyz format (cannot read natoms)'

        # check that file contains enough lines
        # i.e. at least 2 + natoms
        if len(lines) < 2 +self.natoms:
            return 'file: '+self.fname+' does not fit the xyz format (too few lines)'
        
        # get title line
        self.title=str(lines[1])

        # make table:
        table = read_xyz_table(lines[2:self.natoms+2])

        # update self:
        self.read_from_table(table)
        

    def read_from_table(self, table, fname='', title=''):
        """xyz data is just a table with 4 columns and natoms lines"""
        if fname:
            self.fname = fname
        if title:
            self.title = title
        self.natoms = len(table)

        # read labels and xyz coordinates
        charges = []
        xvals = []
        yvals = []
        zvals = []
        for line in table:
            # get atomic number from symbol
            charge = element(line[0]).atomic
            charges.append(charge) 
            xvals.append(float(line[1]))
            yvals.append(float(line[2]))
            zvals.append(float(line[3]))

        self.atom_types = get_clean_atom_types(charges, xvals, yvals, zvals)
        self.natom_types = len(self.atom_types)
        return self.check() # check validity of data"
            
        
    def output(self):
        """ return list of strings as in xyz file """
        lines = []
        lines.append(str(self.natoms))
        lines.append(self.title)
        # print atom types
        for atype in self.atom_types:
            lines.extend(atype.output())

        return lines

    def out_to_file(self, fname=''):
        if self.check():
            return self.check() # check validity of data"

        if not fname:
            fname=self.fname

        with open(fname, 'w') as outfile:
            for line in self.output():
                outfile.write(line + "\n")

def read_xyz_table(lines):
    table = []
    for row in lines:
        line = []
        line.append(row.split()[0])
        line.append(float(row.split()[1]))
        line.append(float(row.split()[2]))
        line.append(float(row.split()[3]))
        table.append(line)
    return table 

def get_clean_atom_types(charges, xvals, yvals, zvals, bohr=False):
    """ order and merge atom types"""
    # get non redundant charges in decreasing order
    charge_types = sorted(set(charges), reverse=True)
    atypes = []
    # build atom types
    for charge in charge_types:
        atype = atom_type()
        # get atomic symbol
        atype.symb = element(charge).symbol
        atype.charge = charge
        atype.bohr = bohr
        for c, x, y, z in zip(charges, xvals, yvals, zvals):
            if c == charge:
                atype.xvals.append(x)
                atype.yvals.append(y)
                atype.zvals.append(z)
                atype.natoms += 1
        
        # add atype to atom_type list
        atypes.append(atype)

    return atypes


