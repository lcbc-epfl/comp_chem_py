#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 11:49:30 2019

@author: pbaudin
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

import comp_chem_utils.spectrum as sp

flist = sys.argv[1:]

exc = np.empty(0)
osc = np.empty(0)
exc = []
osc = []
# read all files
for filename in flist:
    e, f = sp.read_spectrum(filename,"cpmd")
    #exc = np.concatenate( (exc, e) )
    #osc = np.concatenate( (osc, f) )
    exc.append(e)
    osc.append(f)
    
#sp.plot_spectrum(exc, osc)
xpts =[]
ypts = []
# collect spectra
for i,(e,f) in enumerate(zip(exc, osc)):
    x, y = sp.plot_spectrum(e,f,plot=False)
    xpts.append(x)
    ypts.append(y)
    tmp = plt.plot(x, y, label=flist[i])
    plt.vlines(e, -0.1, f, color=tmp[0].get_color())
    
plt.legend()
plt.show()
    
   
   
