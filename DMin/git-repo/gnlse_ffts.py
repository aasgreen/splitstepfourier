# -*- coding: utf-8 -*-
"""
Created on Sat Apr 19 12:08:09 2014

@author: Gabe-Local
"""
from numpy import fft

def FFT_t(A,ax=0):
    return fft.ifftshift(fft.ifft(fft.fftshift(A,axes=(ax,)),axis=ax),axes=(ax,))
def IFFT_t(A,ax=0):
    return fft.ifftshift(fft.fft(fft.fftshift(A,axes=(ax,)),axis=ax),axes=(ax,)) 

# these last two are defined in laserFOAM but never used
def FFT_x(self,A):
        return fft.ifftshift(fft.fft(fft.fftshift(A)))
def IFFT_x(self,A):
        return fft.ifftshift(fft.ifft(fft.fftshift(A)))
