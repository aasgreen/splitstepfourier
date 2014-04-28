# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 14:16:59 2014

@author: dim1
"""
import numpy as np
from scipy.misc import factorial
from scipy import fftpack
from fftw_transforms import fftcomputer as fftw
import odespy

def complexin(x):
    a = np.real(x)
    b = np.imag(x)
    return np.hstack((a,b))
    
def complexout(x):
    a = x[0:len(x)>>1]
    b = x[len(x)>>1:]
    return a + 1j*b

def gnlse(T,A,w0,gamma,betas,gain,fr,RT,fiberlength,nsaves):
    # From Supercontinuum Generation in Optical Fibers
    # Modified slightly to include gain
    
    # Time Grid Parameters
    n = len(T)
    dT = T[1] - T[0]
    V = 2*np.pi*np.transpose(np.arange(-n/2,n/2))/(n*dT)      # Frequency grid
    
    fftc = fftw(n) 
    
    # Calculate g from Gain via 10^(G(dB)/10) = exp(gain)
    alpha = np.log(10**(gain/10))
    
    # Dispersion
    B = 0
    for i in range(len(betas)):
        B = B + betas[i]/factorial(i+2)*V**(i+2)
        
    L = 1j*B - alpha/2
    
    eps = np.finfo(np.float).eps #machine epsilon
    
    if abs(w0) > eps: # If w0>0 then include shock
        gamma = gamma/w0
        W = V + w0 # for a shock this is the frequency
    else:
        W = 1 # W=1 for no shock
    
    RW = n*fftc.ifft(fftpack.ifftshift(np.transpose(RT))) # Frequency domain Raman
    L = fftpack.fftshift(L)
    W = fftpack.fftshift(W) # Reorder for fft space
    
    #**** Define equation to impliment propagation equation
    def rhs(AW, z):
        AW = complexout(AW)
        AT = fftc.fft(AW*np.exp(L*z)) # AT = amplitude, time
        IT = abs(AT)**2 # IT = intensity, time
        if (len(RT)==1 or abs(fr) < eps): # No Raman case
            M = fftc.ifft(AT*IT) # response function
        else:
            RS = dT*fr*fftc.fft(fftc.ifft(IT)*RW) # Raman convolution
            M = fftc.ifft(AT*((1-fr)*IT + RS)) # Response function
        R = 1j*gamma*W*M*np.exp(-L*z) # Full RHS of propagation eqn
        R = complexin(R)
        return R
    
    def isempty(x):
        if np.shape(x) == 0:
            val = 1
        else:
            val = 0
        return val    
    
    #**** Define function to print ODE integrator status
    def report(z, y, flag):
        status = 0
        if isempty(flag):
            print round(z/fiberlength*1000) / 10,'% complete'
        return status
    
    # Set up and run ODE integrator
    Z = np.linspace(0, fiberlength, nsaves)       # set up save points
    
    # Set integrator options
    solver = odespy.Dop853(rhs)
    solver.set_initial_condition(complexin(fftpack.ifft(A)))
    
#    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-12, \
#                'NormControl', 'on', 'OutputFcn', report)
    AW, Z = solver.solve(Z) # Run integrator
    # Process output
    AW = AW[:,:n] + 1j*AW[:,n:] 
    AT = np.zeros(np.shape(AW))
    for i in range(len(AW[:,0])):
        AW[i,:] = AW[i,:]*np.exp(np.transpose(L)*Z[i]) # Change variables
        AT[i,:] = fftpack.ifftshift(fftc.fft(AW[i,:])) # Time domain
        AW[i,:] = fftpack.fftshift(AW[i,:])/dT # Scale
    W = V + w0 # Absolute frequency grid
    
    return [Z, AT, AW, W]