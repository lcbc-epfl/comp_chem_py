#!/usr/bin/env python
"""
Definitions of classes to work with molecular data.

The mol_info_table class is used to store and manipulated 
information about a molecule.

The mol_xyz_table class is used to store and manipulated
the xyz coodinate of a molecule.

The functions defined here are basically wrappers to
MySQL execution lines.
Each function returns a MySQL command that should be
executed by the calling routine.
"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

mol_info = "mol_info"
mol_xyz = "mol_xyz"
sql_idd = 'int unsigned NOT NULL auto_increment'
sql_int = 'int NOT NULL'
sql_str = 'varchar(4) NOT NULL'
sql_flt = 'double(40,20) NOT NULL'
sql_txt = 'text NOT NULL'


# functions in common between the two types of tables
# ---------------------------------------------------
def mysql_add_row(table):
    headers = table.headers[1][0]
    code = '%s'
    for (lab, typ) in table.headers[2:]:
        headers += ', {0}'.format(lab)
        code += ', %s'
    
    return "INSERT INTO {0.name} ({1}) VALUES ({2}) ".format(table, headers, code)


def mysql_create_table(table):
    line = "CREATE TABLE {0.name} ( ".format(table)
    for (lab, typ) in table.headers:
        line += "{0} {1}, ".format(lab, typ) 

    # add primary key:
    line += 'PRIMARY KEY (id) ) '
    return line


# ---------------------------------------------------
class mol_info_table(object):
    """ handle the main database table mol_info"""

    def __init__(self):
        self.name = mol_info
        self.headers = [
                ('id', sql_idd),
                ('name', sql_txt),
                ('chem_name', sql_txt),
                ('note', sql_txt),
                ('charge', sql_int),
                ('natoms', sql_int),
                ('natom_types', sql_int)
                ]
        self.ncol = len(self.headers)

    def create_table(self):
        return mysql_create_table(self)

    def find_duplicates(self, chem_name):
        return 'SELECT id, name FROM {0.name} WHERE chem_name = "{1}"'.format(self, chem_name)

    def get_col(self, headers):
        s = ', '.join(headers) 
        return "SELECT {0} FROM {1.name}".format(s, self)

    def get_row(self, idd):
        return "SELECT * FROM {0.name} WHERE id = {1}".format(self, idd)

    def add_row(self):
        return mysql_add_row(self)
        
    def delete_row(self, idd):
        return "DELETE FROM {0.name} WHERE id={1}".format(self, idd)

    def update(self, col, new, idd):
        return 'UPDATE {0} SET {1}="{2}" WHERE id={3}'.format(self.name, col, new, idd) 



class mol_xyz_table(object):
    """ handle the coordinate database table mol_xyz"""

    def __init__(self, idd):
        self.idd = idd
        self.name = "_".join( [mol_xyz, str(self.idd)] )
        self.headers = [
                ('id', sql_idd),
                ('labels', sql_txt),
                ('charge', sql_int),
                ('xvals', sql_flt),
                ('yvals', sql_flt),
                ('zvals', sql_flt),
                ]
        self.ncol = len(self.headers)

    def create_table(self):
        return mysql_create_table(self)

    def get_table(self):
        return "SELECT * FROM {0.name} ORDER BY charge DESC".format(self)

    def delete_table(self):
        return "DROP TABLE IF EXISTS {0.name}".format(self)

    def add_row(self):
        return mysql_add_row(self)

