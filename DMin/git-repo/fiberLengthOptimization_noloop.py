# -*- coding: utf-8 -*-
"""
Created on Thu Feb 06 17:03:27 2014

@author: dim1
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import csv

os.chdir('O:\\OFM\\Maser\\PM for OPO and DFG\\Fiber Length Optimization Output')
dirlist = os.listdir('.')

n = 2**15 # n points

dT      = 0.001
dTSeconds = dT * 1e-12
twidth  = dT*(n/2.0)
T = np.linspace(-twidth/2,twidth/2,n) #ps

LITfile = 'LIT_SMF0.55EDF1.3new.csv'
Zfile = 'Z_SMF0.55EDF1.3new.csv'
        
SMFlength = 0.55
EDFlength = 1.30
Z = np.genfromtxt(Zfile,delimiter=',')
LIT = np.genfromtxt(LITfile,delimiter=',')
LIT = 10**(LIT/10)

x, y = np.unravel_index(np.argmax(LIT), np.shape(LIT))

outplots = LIT[x, :]

SMFfinal = Z[x] - SMFlength - EDFlength

lengths = 'SMF_'+str(SMFlength+0.5)+'_EDF_'+str(EDFlength)+'_finalSMF_'+str(SMFfinal)

plt.plot(T,outplots)
        
#csv.writer(open('lengths.csv', 'wb'), delimiter=',').writerows(lengths)
#np.savetxt('maxtraces.csv',outplots,delimiter=',')