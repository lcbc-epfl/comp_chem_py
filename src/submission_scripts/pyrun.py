#!/usr/bin/env python
description="""Writes a bash scripts that will submit a calculation to scratch

All the input arguments are global variables in this script!
"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import sys
import os
import argparse

from comp_chem_utils.utils import get_file_as_list
from comp_chem_utils.program_info import executable, choose_exec, avail_progs
from comp_chem_utils.sbatch_info import sbatch_info, sbatch_line


# change default scratch directory depending on computer/cluster
scr='/scratch/pbaudin'

# assuming the submit_job.sh file is in the same dir as the current file
script_dir = os.path.dirname(os.path.abspath(__file__))
script_name= '{}/submit_job.sh'.format(script_dir)

def print_input(verbose=False):
    for key, arg in args.items():
        print arg.output(verbose=verbose)

def read_command_line():
    """Read the submission line and split the arguments.

    Every key is followed by its corresponding value(s).
    Possibilities include:

    key value
    key=value
    key "value1 value2 ..."
    key="value1 value2 ..."

    a list = [key1, value1, key2, value1.1, value1.2, key3, ...]
    is returned
    """

    # split key=value into key, value
    l = []
    for x in sys.argv[1:]:
        l += x.split('=')

    # merge back lists of values between quotations
    lst = []
    s=''
    open_q = False
    for x in l:
        if open_q:
            if x[-1]=='"':
                # finish merging
                open_q = False
                s+=x.split('"')[0]
                lst.append(s)
                s=''
            else:
                # just add value to s
                s+=' '+x
        else:
            if x[0]=='"':
                # start merging
                open_q = True
                s+=x.split('"')[1]
            else:
                # standard value
                lst.append(x)
    return lst


def show_usage():
    print "Usage: {} [-h] [--avail] [--inp] [-f input_file] [args**]\n".format(os.path.basename(__file__))

    print description

    print "\nOptional argunents:"
    print "  -h, --help     Show this help message and exit"
    print "  --avail        List available programs and exit"
    print "  --inp          Print short input file and exit"
    print "  --INP          Print commented input file and exit"
    print "  -f INPUT       Use INPUT as input file to read arguments"
    print "  args**         list of arguments written as shown below:"
    print "                 KEY=VALUE or KEY VALUE"
    sys.exit()

class argument(object):
    def __init__(self, name, help=None, default=None, booltype=False):
        self.name = name
        self.help = help
        self.value = default
        self.default = default
        self.booltype = booltype

    def read_file(self, lines):
        for line in lines:
            key = line.split()[0]
            if key==self.name:
                if self.booltype:
                    self.value = True
                else:
                    self.value = line[4:].strip()

        if not self.value:
            self.value = self.default

        return self.value

    def read_line(self, lst):
        for i, el in enumerate(lst):
            if el == self.name:
                if self.booltype:
                    self.value = True
                else:
                    self.value = lst[i+1]

        if not self.value:
            self.value = self.default

        return self.value

    def output(self, verbose=False):
        if verbose:
            return "# Purpose: {}\n# Default: {}\n{}\n".format(self.help, self.default, self.name)
        else:
            return self.name

args={
        'PRG': argument('PRG', 'Name of the software to be used'),
        'JOB': argument('JOB', 'Name describing the current job'),
        'VRS': argument('VRS', 'Version for the selected software', default=None),
        'INP': argument('INP', 'List of input files or folders (to be copied to SCR)', default =None),
        'OUT': argument('OUT', 'List of output files (to be copied from SCR to WRK)', default = '*'),
        'WRK': argument('WRK', 'Working directory', default='$PWD'),
        'SCR': argument('SCR', 'Scratch directory', default=None),
        'BSH': argument('BSH', 'Name of bash script', default=None),
        'NOX': argument('NOX', 'Do not execute submission batch script', default=False, booltype=True),
        'NOV': argument('NOV', 'No verbose: turn off output printings', default=False, booltype=True),
        'MAIL': argument('MAIL', 'Send a mail after the calculation is done', default=False, booltype=True),

        'SLUR': argument('SLUR', 'Script submitted on a slurm system (with sbatch)', default=False, booltype=True),
        'time': argument('time', 'Max wall time, format: days-hours:minutes', default='0-1:00'),
        'nodes': argument('nodes', 'Number of nodes to request', default='1'),
        'partition': argument('partition', 'Which partition to use (debug or parallel)', default='debug')
        }

def read_input():
    """Try to read the input file or the submission line

    Based on input data, some options or actions are setup
    in case it fails, the usage documentation is printed.
    """
    global PRG,JOB,VRS,INP,OUT,WRK,SCR,BSH,NOX,NOV,MAIL,SLUR,time,nodes,partition
    try:
        if sys.argv[1] in ['-h','--help']:
            show_usage()

        elif sys.argv[1] == ('--avail'):
            avail_progs()
            sys.exit()

        elif sys.argv[1] == ('--inp'):
            print_input()
            sys.exit()

        elif sys.argv[1] == ('--INP'):
            print_input(verbose=True)
            sys.exit()

        elif sys.argv[1] == ('-f'):
            # read from file
            pyrun_inp = get_file_as_list(sys.argv[2])

            PRG = args['PRG'].read_file(pyrun_inp)
            JOB = args['JOB'].read_file(pyrun_inp)
            VRS = args['VRS'].read_file(pyrun_inp)
            INP = args['INP'].read_file(pyrun_inp)
            OUT = args['OUT'].read_file(pyrun_inp)
            WRK = args['WRK'].read_file(pyrun_inp)
            SCR = args['SCR'].read_file(pyrun_inp)
            BSH = args['BSH'].read_file(pyrun_inp)
            NOX = args['NOX'].read_file(pyrun_inp)
            NOV = args['NOV'].read_file(pyrun_inp)
            MAIL = args['MAIL'].read_file(pyrun_inp)

            SLUR = args['SLUR'].read_file(pyrun_inp)
            time = args['time'].read_file(pyrun_inp)
            nodes = args['nodes'].read_file(pyrun_inp)
            partition = args['partition'].read_file(pyrun_inp)

        else:
            # read from command line
            lst = read_command_line()

            PRG = args['PRG'].read_line(lst)
            JOB = args['JOB'].read_line(lst)
            VRS = args['VRS'].read_line(lst)
            INP = args['INP'].read_line(lst)
            OUT = args['OUT'].read_line(lst)
            WRK = args['WRK'].read_line(lst)
            SCR = args['SCR'].read_line(lst)
            BSH = args['BSH'].read_line(lst)
            NOX = args['NOX'].read_line(lst)
            NOV = args['NOV'].read_line(lst)
            MAIL = args['MAIL'].read_line(lst)

            SLUR = args['SLUR'].read_line(lst)
            time = args['time'].read_line(lst)
            nodes = args['nodes'].read_line(lst)
            partition = args['partition'].read_line(lst)

    except Exception as e:
        print e
        show_usage()



def check_input():
    """Some checks on input arguments"""
    global PRG,SCR,BSH

    if not PRG:
        print """Compulsory keywords PRG not found in input file"""
        show_usage()

    if not JOB:
        print """Compulsory keywords JOB not found in input file"""
        show_usage()

    if not SCR:
        SCR = '{}/{}'.format(scr, JOB)

    if not BSH:
        BSH = './submit_{}.sh'.format(JOB)

    if PRG=='select':
        PRG = choose_exec()
    else:
        PRG = executable(PRG, version=VRS)


def print_param():
    """Print out input parameters to terminal"""

    print "BASH SCRIPT PARAMETERS:\n"

    print "{:<26}: {}\n".format("Name of bash script", BSH)

    print "{:<26}: {}/{}".format("Programme used", PRG.prog, PRG.version)
    print "{:<26}: {}".format("Jobname", JOB)
    print "{:<26}: {}".format("Working directory", WRK)
    print "{:<26}: {}\n".format("Scratch directory", SCR)

    print "{:<26}: {}".format("Files copied to scratch", INP)
    print "{:<26}: {}\n".format("Files copied from scratch", OUT)

    print "{:<26}: {}".format("Exectute bash script", not NOX)
    print "{:<26}: {}".format("Output printings", not NOV)
    print "{:<26}: {}".format("Send e-mail", MAIL)
    print "{:<26}: {}".format("Running on SLURM", SLUR)
    if SLUR:
        print "   {:<23}: {}".format("Number of nodes", nodes)
        print "   {:<23}: {}".format("Time limit", time)
        print "   {:<23}: {}".format("Partition", partition)


def create_submission_script():
    """Read, update, and return a template for a bash submission script"""

    template = get_file_as_list(script_name, raw=True)

    nlines = []
    for line in template:

        # add current line to new bash script
        nlines.append( line.strip() )

        # add new info to bash script
        if '# JOBNAME' in line:
            nlines.append( 'job={}'.format(JOB) )

        elif '# PROG ENVIRONMENT' in line:
            nlines.extend ( PRG.get_env() )

        elif '# WORK' in line:
            nlines.append( 'wrk={}'.format(WRK) )

        elif '# SCRATCH' in line:
            nlines.append( 'scr={}'.format(SCR) )

        elif ('# CP FILES TO SCR' in line):
            if INP:
                nlines.append( 'cp -r {} $scr'.format( INP ) )

        elif ('# RUN PROG' in line):
            nlines.append( 'cd $scr' )
            nlines.append( PRG.get_exec_line() )

        elif ('# SAVE FILES FROM SCR' in line):
            if OUT:
                nlines.append( 'out=$wrk/OUT_{}'.format(JOB) )
                nlines.append( 'rm -r $out' )
                nlines.append( 'mkdir -p $out' )
                nlines.append( 'cd $scr' )
                nlines.append( 'cp -r {} $out'.format( OUT ) )
                nlines.append( 'cd -' )

        elif ('# SEND EMAIL' in line):
            if MAIL:
                nlines.append( 'mail -s "$job is done on $hst" pablo.baudin@epfl.ch <<< $text' )

    return nlines


if __name__ == "__main__":

    # read input line or file
    read_input()

    # sanity check on input
    check_input()

    verbose = not NOV

    # print parameters
    if verbose:
        print_param()

    # create submission script from input data
    nlines = create_submission_script()

    # setup submission line
    if SLUR:
       sbatch_info['job-name'].set( JOB )
       sbatch_info['time'].set( time )
       sbatch_info['nodes'].set( nodes )
       sbatch_info['partition'].set( partition )
       sbatch = sbatch_line( BSH )

    else:
       sbatch = 'nohup {0} > {1}.stdo 2> {1}.stde &'.format(BSH, JOB)

    if verbose:
        print '\nSubmission line:'
        print sbatch

    # write script to disk
    with open(BSH,'w') as script:
        for line in nlines:
            script.write('{}\n'.format(line))

        # add submission line as comment at the end of the file
        script.write("# To submit this batch script use the following command:\n")
        script.write("# {}\n".format(sbatch))

    # eventually execute script
    os.system('chmod u+x '+BSH)
    if not NOX:
        os.system(sbatch)


