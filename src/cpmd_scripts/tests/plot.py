#!/usr/bin/env python
"""
Example of a practical use for the cpmd_utils.py module:
Read information from CPMD Trajectory.
Eventually calculate statistical quantities.
Plot data.
"""

import sys
import os
import matplotlib.pyplot as plt

from comp_chem_py.cpmd_utils import read_standard_file, read_ENERGIES


def get_max_traj(f):

    fmax = 0.0
    natoms = f.shape[0]
    for i in range(natoms):
        # get norm of force vector
        ftmp = np.linalg.norm(f[i,:])

        if ftmp > fmax: 
            fmax = ftmp

    return fmax

def get_mean_traj(f):

    fmean = 0.0
    natoms = f.shape[0]
    for i in range(natoms):
        # get norm of force vector
        fmean += np.linalg.norm(f[i,:])

    return fmean/natoms


# calculation directory
dr = sys.argv[1]+'/'
fig, axarr = plt.subplots(2, sharex=True)


try:
    ts = float(sys.argv[2])
    steps = np.arange(0.0, ts*100000, ts)
    print 'steps from input:',len(steps),steps[-1]
except:
    # extract time information from cpmd.out
    fn = sys.argv[1]+'.out'
    steps = [0.0]
    tstp = [0.0]
    with open(dr+fn, 'r') as cpmd:
        for line in cpmd:
            if 'New value for dt ions' in line:
                tstp.append( float( line.split()[-1] ) )
                steps.append( steps[-1] + tstp[-1] )
    print 'steps from file:',len(steps),steps[-1]
    axarr[1].plot(steps, tstp, label='Time step' )




# ENERGIES
dic = read_ENERGIES(dr+'ENERGIES', ['E_KS', 'E_cla'])

n = len(dic['E_cla'])

axarr[0].axhline(dic['E_cla'][0], color='k')
axarr[0].plot(steps[:n], dic['E_cla'], label='Total (conserved) energy')
#axarr[0].plot(dic['E_KS'], label='Potential energy')
axarr[0].set_ylabel('Energies [a.u.]')
axarr[0].legend()

info = read_standard_file(dr+'SH_ENERG.dat')[1]
axarr[0].plot(steps[:n], info[:n,0], label='S_0')
axarr[0].plot(steps[:n], info[:n,1], '+', label='S_1')
axarr[0].plot(steps[:n], info[:n,2], label='S_2')
axarr[0].plot(steps[:n], info[:n,3], label='S_3')


## VELOCITIES AND FORCES
#veloc, forces = read_FTRAJECTORY(dr+'FTRAJECTORY')[2:4]
#
#fmax = []
#fmean = []
#vmax = []
#vmean = []
#for v, f in zip(veloc, forces):
#    fmax.append(get_max_traj(f))
#    fmean.append(get_mean_traj(f))
#
#    vmax.append(get_max_traj(v))
#    vmean.append(get_mean_traj(v))
#
#axarr[1].plot(fmax, label='max force')
#axarr[1].plot(fmean, label='mean force')
#axarr[1].set_ylabel('Forces [a.u.]')
#axarr[1].legend()
#
#axarr[2].plot(vmax, label='max velocity')
#axarr[2].plot(vmean, label='mean velocity')
#axarr[2].set_ylabel('Velocity [a.u.]')
#axarr[2].legend()
#
#
## TRANSITION PROBABILITY
#probs = read_standard_file(dr+'SH_PROBS.dat')[1]
#
#pmax = probs.max(1)
#pmean = probs.mean(1)
#
#axarr[3].plot(pmax, label='max proba')
#axarr[3].plot(pmean, label='mean proba')
#axarr[3].set_ylabel('Transition probability')
#axarr[3].legend()



fig.subplots_adjust(hspace=0.0)
plt.show()

