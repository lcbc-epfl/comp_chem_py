#!/usr/bin/env python

__author__="Pablo Baudin"
__email__="pablo.baudin@epfl.ch"


import matplotlib.pyplot as plt 
import scipy.linalg as sc
import numpy as np      

from PDOS.modules.read_gaussian import get_gaussian_info


def plot_dos_and_pdos(nocc, eps, xmin, xmax, xpts, groups, pdos, gaussian_file, plot_dos, dos):
    """plot DOS and PDOS: This routine can be modify to alter the final display of the graph"""

    # plot orbital energy sticks
    # --------------------------
    for iorb, energy in enumerate(eps):

        # if the orbital energy is in the x range
        if (energy >= xmin) and (energy <= xmax):
            if (iorb < nocc):
                # plot occupied orbitals in cyan
                col ='c'
            else:
                # plot virtual orbitals in magenta
                col ='m'

            # plot stick at the bottom of the window
            plt.axvline(energy, 0, 0.1, color=col)
     

    # plot the total density of states (DOS)
    # --------------------------------------
    ymax = 0.0
    if (plot_dos):

        plt.plot(xpts, dos, "r:", label="DOS")
        ymax = dos.max()
   

    # plot the partial DOS for each group
    # -----------------------------------
    for igrp,groupname in enumerate(groups):

        plt.plot(xpts, pdos[igrp], label=groupname)
        ymax = max(ymax, pdos[igrp].max())
    

    # final touch plus show it!
    # -------------------------
    ft = 16
    # y-range is defined such that:
    #   y=0.0 is at 10% from the bottom of the window (bellow are the sticks)
    #   there is a buffer space of 10% with the largest (P)DOS value
    ymin = -ymax/8.0
    ymax = -9.0*ymin
    plt.ylim([ymin,ymax])
    plt.xlim([xmin,xmax])

    plt.xlabel("Energy [eV]", fontsize=ft)
    plt.ylabel("(Partial) DOS [eV$^{-1}$]", fontsize=ft)
    plt.title("DOS & PDOS of "+gaussian_file, fontsize=ft)

    plt.tick_params(axis='both', labelsize=ft)
    plt.legend(fontsize=ft)
    plt.tight_layout()

    plt.show()



def calculate_and_plot_pdos(gaussian_file, group_file, npts=0, fwhm=1.0, xmin=1, xmax=-1, plot_dos=True, mulliken=False):
    """Calculate and plot partial density of states (PDOS)
    
    This function read a gaussian output files containing the 
    necessary information for calculating total and partial
    density of states. The DOS and PDOS are then calculated 
    and plotted."""
  
    # get gaussian info:
    #   nocc:    number of occupied orbitals (assuming closed shell)
    #   nbas:    number of basis functions (assuming #AOs = #MOs)
    #   overlap: AO overlap matrix
    #   eps:     molecular canonical orbitals energies
    #   coef:    canonical MO to AO transformation matrix
    nocc, nbas, overlap, eps, coef = get_gaussian_info(gaussian_file)

    if (mulliken):
        print("Calculating Mulliken partial charges...")
        # get partial Mulliken charges
        # q_{alpha,i} = sum_{beta} S_{alpha,beta} C_{beta,i} C_{alpha,i}  
        part_pop = np.dot(overlap, coef)
        part_pop = np.multiply(part_pop, coef)
    else:
        print("Calculating Lowdin partial charges...")
        # get partial Lowdin charges
        # q_{alpha,i} = sum_{beta} S_{alpha,beta}^{1/2} C_{beta,i} C_{alpha,i} S_{alpha,beta}^{1/2} 
        Shalf = sc.sqrtm(overlap)
        part_pop = np.dot(Shalf, coef)**2.0
    
    
    # Calculation of gross orbital populations
    total_pop = part_pop.sum(0) # sum over all AOs (rows=0)


    # make groups of atomic orbitals based on file "group_file"
    # ---------------------------------------------------------
    print("Reading group of AOs file: {}\n".format(group_file))
    inputfile = open(group_file,"r")

    # there are different groups of orbitals, for each group we have a 
    # groupname to which a list of atomic orbital is associated
    groups = {}
    
    for line in inputfile:
        # the first line give the name of the group
        groupname = line.strip()
        orbitals = []
        line = inputfile.next()
        parts = line.split(",")
        for x in parts:
            temp = x.split("-")
            try:
                if len(temp)==1:
                    orbitals.append(int(temp[0]))
                else:
                    orbitals.extend(range(int(temp[0]),int(temp[1])+1))
            except:
                print("Warning: error while reading {}".format(group_file))

        # associate list of orbitals to each group
        groups[groupname] = orbitals
    
    inputfile.close()
    

    # Define the parameters of the gaussian convolution based
    # on user inputs or predefined defaults
    # -------------------------------------------------------
    print("Gaussian convolution of data...")

    # get x-axis range
    if (xmax - xmin) <= 0:
        xmax = (eps[-1] + 5)
        xmin = (eps[0] - 5)
    
    # default is to have 100 points per eV
    if (npts <= 0):
        npts = int((xmax - xmin)*100)

    # set x-values 
    xpts = np.linspace(xmin, xmax, npts)

    # the gaussian curve is defined as g(x) = norm * exp (- alpha * delta_x**2) where,
    alpha = 4.0*np.log(2.0)/(fwhm*fwhm)
    norm = np.sqrt(alpha/np.pi)
    dos = np.zeros(npts)

    print("xmax = {} [eV]".format(xmax))
    print("xmin = {} [eV]".format(xmin))
    print("npts = {}".format(npts))
    print("FWHM = {}".format(fwhm))
    print("norm = {}".format(norm))
    print("alpha = {}".format(alpha))


    # calculate the total density of states (DOS)
    # -------------------------------------------
    if (plot_dos):

        # convolution the total population of each MO:
        for iorb, pop in enumerate(total_pop):

            dos += pop*norm*np.exp(-alpha*(eps[iorb]-xpts)**2.0)
    

    # calculate the partial DOS for each group
    # ----------------------------------------
    pdos = []
    for igrp,groupname in enumerate(groups):

        # convolution the partial population of each MO,
        # In the part_pop array the rows correspond to AOs
        # while the column correspond to MOs. We therefore
        # want to sum up the AOs (the rows) corresponding to
        # each groupA

        # add the population of each AO in the group
        my_part_pop = np.zeros(nbas)
        for elmt in groups[groupname]:

            my_part_pop += part_pop[elmt-1,:]

        pdos.append(np.zeros(npts))
        # convolution the partial population of each MO:
        for iorb, pop in enumerate(my_part_pop):

            pdos[igrp] += pop*norm*np.exp(-alpha*(eps[iorb]-xpts)**2.0)
        
    print("\nConvolution of data done!")
    

    # make actual plots
    # -----------------
    plot_dos_and_pdos(nocc, eps, xmin, xmax, xpts, groups, pdos, gaussian_file, plot_dos, dos)
 

