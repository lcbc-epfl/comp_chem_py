#!/usr/bin/env python
"""
Determine optimal CP_GROUP value for CPMD MPI calculations.

This is based on empirical investigation and is not guarented
to give optimal performance at all!

The attempt here is:

* to have a uniform distribution of plane waves on all CP groups.

* to minimize CP_GROUP with that constraint, where CP_GROUP is the
  number of MPI tasks per CP_GROUP, i.e. CP_GROUP size.

* to make the number of MPI tasks a multiple of the number of CPUs on a node.
"""

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"

import pandas as pd


if __name__=="__main__":

    max_nnodes      = int(input("Maximum number of nodes that can be used   :  "))

    # if a node has 16 CPUs this should probably be 16
    ntasks_per_node = int(input("Number of MPI tasks per node               :  "))

    # This can be obtained by submitting a dummy CPMD calculation
    # with your system and plane wave cutoff. Look for "REAL SPACE MESH".
    mesh            = int(input("Size of the REAL SPACE MESH in X direction :  "))


    foundit = False
    solutions = {
        'Number of nodes':[],
        'CP_GROUP (size of)':[],
        'Number of PW per group':[],
    }

    # We just go brute force and try all possibilities
    for N  in range(1,max_nnodes+1):

        # CP_GROUP has to be lower than the number of MPI tasks
        ntasks = N * ntasks_per_node
        for ntasks_per_cpgrp in range(1,ntasks):

            # check that ntasks_per_cpgrp is a divisor of ntasks
            if ntasks%ntasks_per_cpgrp == 0:
                n_cp_groups = ntasks//ntasks_per_cpgrp
            else:
                continue

            # check that the number of CP_GROUPs  is a divisor of the mesh size
            if mesh%n_cp_groups == 0:
                # we have a winner
                foundit = True
                npw_per_cpgrp = mesh//n_cp_groups

                solutions['Number of nodes'].append( N )
                solutions['CP_GROUP (size of)'].append( ntasks_per_cpgrp )
                solutions['Number of PW per group'].append( npw_per_cpgrp )

                # this one is the first value of CP_GROUP that works with
                # the current number of nodes. It is thuse the smalest value
                # which is what we want, so we stop here for that number of nodes
                break

            else:
                continue

    print()
    if foundit:
        # convert solutions to a pandas data frame
        # and sort it based on the number of plane waves
        df = pd.DataFrame(data=solutions)
        df = df.sort_values(by='Number of PW per group')
        print( df.to_string(index=False) )

    else:
        print( "Could not find an optimal solution... :(" )
        print( "Try to increase the number of nodes." )
        print( "Or change the MESH size." )



