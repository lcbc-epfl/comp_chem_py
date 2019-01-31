#!/usr/bin/env python
"""
Stochastic selection of excited state based on oscillator strength

This simple scripts reads an output file from CPMD.
In particular the oscillator strengths from TDDFT calculation.
And stochastically select a state based on the different strengths.
"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import sys
from numpy.random import choice

from comp_chem_utils.spectrum import read_spectrum

# for compatibility with autodoc in sphinx
if __name__ == '__main__':

    verbose = False
    
    name = 'OUT_{0:02d}_tddft_pbe0_10x/{0:02d}_tddft_pbe0_10x.out'
    name = 'OUT_{0:02d}_tddft_pbe_10x/{0:02d}_tddft_pbe_10x.out'
    n = 1 # int(raw_input('\nSize of sample:\n'))
    nstate = 6
    
    full = []
    # get oscillator strength from cpmd 
    for idx in range(50):
    #for fn in sys.argv[1:]:
    
        fn = name.format( idx+1 )
        f = read_spectrum(fn, 'cpmd', verbose)[1]
        f = f[0:nstate]
        #nstate = len(f)
        
        # get probability from f
        tot = sum(f)
        p = [x/tot for x in f]
        
        if verbose:
            print "State transition probability:"
            for i,(j,k) in enumerate(zip(f,p)):
                print ' {:4d}  {:10.5f}  {:4.2f} '.format(i+1,j,k)
        
        
        states = [x+1 for x in choice(nstate, n, p=p)]
        
        print '{}: state = {}'.format(fn, states[0])
        full.extend(states)
    
        # checking
        #print ''
        #for i in range(nstate):
        #    print '   {:4d}    {:4.2f}'.format(i+1, states.count(i+1)/float(n) )
    
    print ' '.join(map(str, full))
    print ''
    for i in range(nstate):
        print '   {:4d}    {:4.2f}'.format(i+1, full.count(i+1)/float(len(full)) )
