# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 15:08:15 2014

@author: dim1
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

data = np.genfromtxt('florian_frog.txt',delimiter=',')
wl = data[:,0]
phase = data[:,1]
f = interp1d(wl,phase)

newwl = np.linspace(wl[0],wl[-1],1e3)
newphase = f(newwl)

plt.plot(newwl,newphase)
plt.show()