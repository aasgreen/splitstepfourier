# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 16:41:04 2014

@author: dim1
"""
import numpy as np
from scipy.interpolate import interp1d

def getGainsForGNLSE(freqs,fiber):
    baseDir = 'O:\\OFM\\Maser\\PM for OPO and DFG\\'
    fiberRoot = str(baseDir)
    if fiber == 'er7pm':
        fiberDataFile = 'er7pm_gain.csv'
    fiberData = np.genfromtxt(str(fiberRoot+fiberDataFile),delimiter=',')
    freqAxis = fiberData[0,:]
    gainInterp = interp1d(freqAxis[::-1], fiberData[1,::-1], kind='linear',bounds_error=False,fill_value=0)    
    fractGain = gainInterp(freqs)
    
    return fractGain