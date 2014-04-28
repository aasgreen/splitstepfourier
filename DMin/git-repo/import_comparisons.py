# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:29:52 2014

@author: dim1
"""
import numpy as np
import matplotlib.pyplot as plt

pyfull = np.loadtxt('pydump.csv',delimiter=',')
pydata = pyfull[0:len(pyfull[:,0])>>1,:] + 1j * pyfull[len(pyfull[:,0])>>1,:]

matfull = np.loadtxt('MATLABdump.csv',delimiter=',')
matfull = np.transpose(matfull)
matdata = matfull[0:len(matfull[:,0])>>1,:] + 1j * matfull[len(matfull[:,0])>>1,:]

def flatten(x):
    return x[0,:]

pydata = flatten(pydata)
matdata = flatten(matdata)

plt.plot(np.abs(pydata))
plt.plot(np.abs(matdata))
plt.show()