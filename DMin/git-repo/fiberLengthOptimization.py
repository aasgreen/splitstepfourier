# -*- coding: utf-8 -*-
"""
Created on Thu Feb 06 17:03:27 2014

@author: dim1
"""
import numpy as np
#import matplotlib.pyplot as plt
import os
import csv

os.chdir('O:\\OFM\\Maser\\PM for OPO and DFG\\Fiber Length Optimization Output')
dirlist = os.listdir('.')

#T = np.genfromtxt('T.csv',delimiter=',') #time grid, in ps

j = 0

for each in dirlist:
    if each.find("LIT_SMF") == 0:
        j += 1

lengths = []

i = 0

for each in dirlist:
    if each.find("Z_SMF") == 0:
        
        SMFlength = np.float(each[each.find("SMF") + 3 : each.find("EDF")])
        EDFlength = np.float(each[each.find("EDF") + 3 : each.find(".csv")])
        Z = np.genfromtxt(each,delimiter=',')
        LITfile = "LIT" + each[each.find("_SMF"): ]
        LIT = np.genfromtxt(LITfile,delimiter=',')
        LIT = 10**(LIT/10)

        if i == 0:
            outplots = np.zeros((j, len(LIT[0, :])))
        
        x, y = np.unravel_index(np.argmax(LIT), np.shape(LIT))
        
        outplots[i, :] = LIT[x, :]

        SMFfinal = Z[x] - SMFlength - EDFlength

        lengths.append('SMF_'+str(SMFlength+0.5)+'_EDF_'+str(EDFlength)+'_finalSMF_'+str(SMFfinal))

        i += 1
        
csv.writer(open('lengths.csv', 'wb'), delimiter=',').writerows(lengths)
np.savetxt('maxtraces.csv',outplots,delimiter=',')