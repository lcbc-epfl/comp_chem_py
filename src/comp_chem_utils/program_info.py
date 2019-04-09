#!/usr/bin/env python
__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import sys

class executable(object):
    """This class is used to deal with available programs and their environments"""

    def __init__(self, prog, version=None):
        self.prog = prog
        self.version = version

    def get_exe(self):
        exe=''
        if self.prog==LSDALTON:

            if self.version=='OMP 2018':
                exe = '/home/pbaudin/Work/softwares/lsdalton_2018/build_omp_intel17_new/lsdalton.x'

            elif self.version=='MPI CPSD':
                exe = '/home/pbaudin/Work/softwares/lsdalton_mp_cpsd/build/lsdalton.x'

        elif self.prog==ADF:
            exe = 'adf'

        elif self.prog==TURBO:
            # for turbmole executables are job specific
            exe = ''

        elif self.prog==CPMD:

            if self.version=='ELISA MTS':
                exe = '/home/pbaudin/Work/MTS_CPMD/cpmd_elisa/bin/cpmd.x'

            elif self.version=='PABLO REFERENCE MTS':
                exe = '/home/pbaudin/Work/MTS_CPMD/cpmd_reference/bin/cpmd.x'

            elif self.version=='MARTIN SVN':
                exe = '/home/pbaudin/Work/softwares/cpmd_martin/bin/cpmd.x'

            elif self.version=='PABLO PORTING MTS':
                exe = '/home/pbaudin/Work/MTS_CPMD/cpmd_porting/CPMD/bin_elisa/bin/cpmd.x'

            elif self.version=='PABLO PORTING MASTER':
                exe = '/home/pbaudin/Work/MTS_CPMD/cpmd_porting/CPMD/bin_master/bin/cpmd.x'

            elif self.version=='PABLO PORTING TDMTS':
                exe = '/home/pbaudin/Work/MTS_CPMD/cpmd_porting/CPMD/bin_tdmts/bin/cpmd.x'

            elif self.version=='FIDIS TDMTS':
                exe = '/home/pbaudin/cpmd_td_mts/CPMD/bin_intel/bin/cpmd.x'

            elif self.version=='ML PLAYGROUND':
                exe = '/home/pbaudin/Work/MTS_CPMD/cpmd_porting/CPMD/bin_ml_playground/bin/cpmd.x'

            else:
                exe = 'cpmd.x'

        else:
            sys.exit('Unknowm program: {}'.format(self.prog) )

        return exe


    def get_env(self):

        env = ['source /software/ENV/modules']
        env.append( "EXE='{}'".format( self.get_exe() ) )

        if self.prog==ADF:
            if self.version:
                env.append( 'module load adf/{}'.format(self.version) )
            else:
                env.append( 'module load adf' )

        elif self.prog==LSDALTON:
            if self.version=='OMP 2018':
                env.append( 'module load intel/17.0.4' )

            elif self.version=='MPI CPSD':
                env.append( 'module load intelmpi/17.0.4' )

        elif self.prog==TURBO:
            if self.version:
                env.append( 'module load turbomole/{}'.format(self.version) )
            else:
                env.append( 'module load turbomole' )
            env.append( 'export PARA_ARCH=SMP' )
            env.append( 'export PATH=$TURBODIR/bin/em64t-unknown-linux-gnu_smp:$PATH' )
            env.append( 'export PARNODES=$np' )


        elif self.prog==CPMD:
            env.append( "export PP_LIBRARY_PATH='/software/cpmd/PPs/'")

            if not self.version:
                env.append( 'module load cpmd' )

            elif self.version=='RELEASE 4.1':
                env.append( 'module load cpmd/4.1-QMMM/intel-15.0.3' )

            elif self.version=='POLYDEFKIS QMMM':
                env.append( 'module load cpmd/4.1-QMMM/intel-15.0.3-polydefkis' )

            elif self.version=='FIDIS TDMTS':
                # On FIDIS there is no need to source the modules package
                # so I just remove it here
                env.pop(0)
                env.pop(1)
                env.append( "export PP_LIBRARY_PATH='/home/pbaudin/cpmd_PPs/'" )
                env.append( 'module load intel intel-mpi intel-mkl' )

            else:
                env.append( 'module load intel/15.0.3' )
                env.append( 'module load intelmpi/5.1.1' )
                env.append( 'export I_MPI_WAIT_MODE=1' )
                env.append( 'export I_MPI_PIN_DOMAIN=omp' )

        else:
            sys.exit('Unknowm program: {}'.format(prog) )

        return env


    def get_exec_line(self):

        exec_line = ''
        if self.prog==ADF:
            exec_line = '$EXE -n $np < ${job}.inp  > ${job}.out'

        elif self.prog==LSDALTON:
            if self.version=='OMP 2018':
                exec_line = '$EXE'

            elif self.version=='MPI CPSD':
                exec_line = 'mpirun -np $np $EXE'

        elif self.prog==TURBO:
            exec_line = 'dscf > ${job}.out'

        elif self.prog==CPMD:
            if self.version=='FIDIS TDMTS':
                exec_line = 'srun $EXE ${job}.inp > ${job}.out'
            else:
                exec_line = 'mpirun -np $np $EXE ${job}.inp  > ${job}.out'

        return exec_line


# an executable is defined by a program name and a version

# a program is defined by a name and has a set of versions available as a list

# a version is defined by a string

# to each executable is associated an environment and a execution line

# all programs and versions
ADF='adf'
LSDALTON='lsdalton'
CPMD='cpmd'
TURBO='turbo'

programs = {
        ADF: ['2013.01b',  '2016.101',  '2016.107',  '2017.113'],
        LSDALTON: ['OMP 2018','MPI CPSD'],
        TURBO: ['6.5','7.1.1'],
        CPMD: [
    'ELISA MTS',
    'PABLO REFERENCE MTS',
    'MARTIN SVN',
    'RELEASE 4.1',
    'POLYDEFKIS QMMM',
    'PABLO PORTING MTS',
    'PABLO PORTING MASTER',
    'PABLO PORTING TDMTS',
    'FIDIS TDMTS',
    'ML PLAYGROUND',
    ]
    }


def avail_progs():
    print('')
    print('AVAILABLE PROGRAMS: VERSIONS:')
    print('=============================')
    for name in programs:
        print('')
        for vrs in programs[name]:
            print('    {}: {}'.format(name, vrs))
    print('')

def choose_exec():
    avail_progs()

    prog = raw_input('ENTER SELECTED PROGRAM:\n')
    if prog not in programs:
        sys.exit('ERROR: program not available!')

    print('')
    print('AVAILABLE VERSIONS:')
    print('==================')
    for vrs in programs[prog]:
        print('    {}'.format(vrs) )

    print('')
    vrs = raw_input('ENTER SELECTED VERSION:\n')
    if vrs not in programs[prog]:
        print('WARNING: default version selected!')

    # return executable
    return executable(prog, vrs)



