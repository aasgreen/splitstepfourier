# -*- coding: utf-8 -*-
"""
Created on Wed Jan 08 15:39:57 2014

@author: dim1
"""
import numpy as np
import pyfftw

class fftcomputer:
    def __init__(self,gridsize):
        self.gridsize = gridsize
        self.a = pyfftw.n_byte_align_empty(gridsize,16,'complex128')
        self.b = pyfftw.n_byte_align_empty(gridsize,16,'complex128')
        self.fftw = pyfftw.FFTW(self.a,self.b,direction='FFTW_FORWARD')
        
        self.c = pyfftw.n_byte_align_empty(gridsize,16,'complex128')
        self.d = pyfftw.n_byte_align_empty(gridsize,16,'complex128')
        self.ifftw = pyfftw.FFTW(self.c,self.d,direction='FFTW_BACKWARD')
    def fft(self,data, copy = True):
        self.a[:] = data
        if copy:
            return np.copy(self.fftw())
        else:
            return self.fftw()
            
    def ifft(self,data, copy = True):
        self.c[:] = data
        if copy:
            return np.copy(self.ifftw())
        else:
            return self.ifftw()