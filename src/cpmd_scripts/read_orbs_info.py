#!/usr/bin/env python
"""
Read CPMD TDDFT output file and extract informations about the orbitals involved in the electronic transitions.

A CPMD input string is printed which will enable to generate grid files to plot all the
orbital involved in the electronic transitions.
"""

import sys
from comp_chem_utils.utils import get_file_as_list

def get_orbital_energies(lines):
    eocc = []
    evir = []
    for i,line in enumerate(lines):

        # READ ORBITAL ENERGIES
        if 'EIGENVALUES(EV) AND OCCUPATION' in line:

            idx = i + 1
            myl = lines[idx].split()
            while 'CHEMICAL' not in myl:

                # read first orbital energy in the line
                if int(float(myl[2])) == 2.0:
                    eocc.append(float(myl[1]))
                elif int(float(myl[2])) == 0.0:
                    evir.append(float(myl[1]))

                # check for second orbital energy
                if len(myl) == 6:

                    # read second orbital energy in the line
                    if int(float(myl[-1])) == 2.0:
                        eocc.append(float(myl[-2]))
                    elif int(float(myl[-1])) == 0.0:
                        evir.append(float(myl[-2]))

                idx += 1
                myl = lines[idx].split()

    return eocc, evir


def get_transition_orbitals(lines):
    homo = []
    lumo = []
    for i,line in enumerate(lines):

        # READ ORBITALS INVOLVED IN TRANSITIONS
        if all((word in line for word in ['TRANSITION', 'HOMO', 'LUMO'])):
            homo.append( int(line.split()[-5]) )
            lumo.append( int(line.split()[-1]) )


    # convert list of orbitals involved in transitions to indices
    homo = list(set(homo))
    homo.sort(reverse=True)
    lumo = list(set(lumo))
    lumo.sort()
    return homo, lumo

error = """
Please provide a CPMD TDDFT output file as argument:

python read_orbs_info.py cpmd_tddft.out"""

if __name__=="__main__":

    # read CPMD output file
    try:
        cp = get_file_as_list(sys.argv[1])
    except:
        sys.exit(error)

    eocc, evir = get_orbital_energies(cp)

    # homo and lumo contain a list of indices of the orbital involved in the electronic transisions
    # the indices are with respect to the homo and lumo orbitals
    # For example homo = [0,3,7] means that the HOMO, HOMO-3 and HOMO-7 occupied orbitals are involved.
    # lumo = [0,2,4] means that the LUMO, LUMO+2 and LUMO+4 virtual orbitals are involved.
    homo, lumo = get_transition_orbitals(cp)

    norbs = len(homo) + len(lumo)
    nocc = len(eocc)
    nvir = len(evir)

    # format strings for CPMD input that plot orbitals
    str_orb = ''
    for orb in homo:
        str_orb += '-{} '.format( nocc - orb )
    for orb in lumo:
        str_orb += '-{} '.format( nocc + 1 + orb )


    print("""
    RESTART WAVEFUNCTION COORDINATES LINRES LATEST
    KOHN-SHAM ENERGIES
        {}
    RHOOUT BANDS
        {}
        {}
    """.format(nvir, norbs, str_orb)
    )

