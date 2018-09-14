#!/usr/bin/env python
description="""plot potential energy sufaces from CPMD calculations

The main task here is to extract ENERGIES and TIMES information
and eventually plot it"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt

from comp_chem_utils.cpmd_utils import read_ENERGIES, get_time_info, read_SH_ENERG, read_MTS_EXC_ENERG
from comp_chem_utils.conversions import AU_TO_FS


def get_nstates():
    try:
        nstates = int(raw_input("\nNumber of states to read (empty=all):\n"))
    except(ValueError):
        nstates = None

    return nstates


def plot_data_in_dict(data):

    tstep = data['time_step']
    times = AU_TO_FS * tstep * np.asarray(data['steps'])

    for lab in data:
        if lab in ['steps', 'time_step']:
            continue
        elif 'driving' in lab.lower():
            plt.plot(times, data[lab], 'k', label=lab, linewidth=3)
        else:
            plt.plot(times, data[lab], label=lab)


# for compatibility with autodoc in sphinx
if __name__ == '__main__':

   files_data = {}
   
   nfiles = len(sys.argv[1:])
   iplot=1
   for ifiles, cpmd_inp in enumerate(sys.argv[1:]):
   
       jobname = cpmd_inp.split('/')[-1][:-4]
       folder = '/'.join( cpmd_inp.split('/')[:-1] )
       if folder=='':
           folder='./'
       else:
           folder+='/'
       
       energy_files = sys.argv[2:]
       
       # READ TIME INFO
       TIMESTEP, MAXSTEP, USE_MTS, MTS_FACTOR, MTS_TSH = get_time_info(cpmd_inp)
       print """
       TIMESTEP: {} [a.u.]
       MAXSTEP: {}
       USE_MTS: {}
       MTS FACTOR: {}
       MTS TSH LEVEL: {}
       """.format(TIMESTEP, MAXSTEP, USE_MTS, MTS_FACTOR, MTS_TSH)
       
       data = {}
       
       # check for energy files in folder
       print """
       Looking for files with energy data:
            ENERGIES
            SH_ENERG.dat
            MTS_EXC_ENERG.dat
       """
       
       energy_files = []
       if os.path.isfile(folder+'ENERGIES'):
           if raw_input('Use ENERGIES file? [y/n]\n').strip().lower() == 'y':
               energy_files.append(folder+'ENERGIES')
       
       if os.path.isfile(folder+'SH_ENERG.dat'):
           if raw_input('Use SH_ENERG.dat file? [y/n]\n').strip().lower() == 'y':
               energy_files.append(folder+'SH_ENERG.dat')
       
       if os.path.isfile(folder+'MTS_EXC_ENERG.dat'):
           if raw_input('Use MTS_EXC_ENERG.dat file? [y/n]\n').strip().lower() == 'y':
               energy_files.append(folder+'MTS_EXC_ENERG.dat')
           
       if not energy_files:
           sys.exit('None of the relevant energy files have been found in: {}'.format(folder))
       
       for energies in energy_files:
       
           print("\n READING: {}".format(energies) )
           if 'ENERGIES' in energies:
               print("\nChoose info to read (space separated list)")
               print("   E_kel:   Electronic kinetic energy (only for CPMD)")
               print("   Temp:    Temperature [K]")
               print("   E_KS:    Kohn-Sham energy [a.u.]")
               print("   E_cla:   Classical energy, E_KS + E_kin (constant for BOMD)")
               print("   E_ham:   0 for BOMD")
               print("   RMS:     Nuclear displacement wrt initial position (?)")
               print("   CPU_t:   CPU time")
               e_codes = ['steps']
               e_codes.extend( raw_input("").split() )
           
               data['ENERGIES'] = read_ENERGIES(energies, e_codes)
           
               # add time step info
               data['ENERGIES']['time_step'] = TIMESTEP
           
           
           elif 'SH_ENERG.dat' in energies:
           
               nstates = get_nstates()
           
               if USE_MTS and MTS_TSH == 'HIGH':
                   data['SH_ENERG'] = read_SH_ENERG(energies, nstates, MTS_FACTOR) 
               else:
                   data['SH_ENERG'] = read_SH_ENERG(energies, nstates) 
           
               # add time step info
               data['SH_ENERG']['time_step'] = TIMESTEP
           
           
           elif 'MTS_EXC_ENERG.dat' in energies:
               if not USE_MTS:
                   sys.exit('Incoherent input: read MTS_EXC_ENERG.dat and not USE_MTS!!')
           
               read_high = raw_input("\nRead MTS high level energies? [y/n]\n").strip().lower() == 'y'
               if read_high:
                   nstates = get_nstates()
           
                   data['HIGH_MTS_EXC_ENERG'] = read_MTS_EXC_ENERG(energies, nstates, MTS_FACTOR, HIGH=True)
           
                   # add time step info
                   data['HIGH_MTS_EXC_ENERG']['time_step'] = TIMESTEP
           
           
               read_low  = raw_input("\nRead MTS low level energies? [y/n] \n").strip().lower() == 'y'
               if read_low:
                   nstates = get_nstates()
                   data['LOW_MTS_EXC_ENERG'] = read_MTS_EXC_ENERG(energies, nstates, MTS_FACTOR, HIGH=False)
           
                   # add time step info
                   data['LOW_MTS_EXC_ENERG']['time_step'] = TIMESTEP
   
           files_data[jobname] = data
   
   
       # Now we shoud have a dictionary of data ready to plot:
       # each data set is also organized as a dictionary
       # there should be at least one entry with key 'steps'
       # and one entry with key time_step that will be used to 
       # determine the time info.
       # The rest of the entries will be printed to file or plotted.
       
       # Print out data ready for plotting
       for title, info in data.items():
   
           plt.subplot(nfiles, len(data), iplot)
           plot_data_in_dict( info )
           iplot +=1
   
           with open('{}{}_{}.dat'.format(folder, jobname, title), 'w') as out:
            
               time_step = info['time_step']
               steps = info['steps']
            
               del info['time_step']
               del info['steps']
            
               # make list of labels:
               labels = ['steps (fs)']
               labels.extend( list(info.keys()) )
            
               # print it as header
               header = '#' + len(labels) * ' {:^20}' + '\n'
               out.write( header.format(*labels) )
            
               content = ' '+len(labels)*' {:^20.12g}'+'\n'
               for i,step in enumerate(steps):
                   line = [ time_step * AU_TO_FS * step ]
                   for lab in labels[1:]:
                       line.append( info[lab][i] )
            
                   out.write( content.format( *line ) ) 
   
   
   plt.xlabel('Time [fs]')
   #plt.ylabel('Energy [a.u.]')
   plt.legend()
   plt.show()

