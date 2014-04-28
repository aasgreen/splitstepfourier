# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 17:09:10 2014

@author: dim1
"""
import numpy as np
from scipy.misc import factorial
from getGainsForGNLSE import getGainsForGNLSE
from scipy import fftpack
from fftw_transforms import fftcomputer as fftw
import odespy
#import tables
#import matplotlib.pyplot as plt
#from scipy.integrate import ode

def complexin(x):
    y = np.zeros(2*len(x))
    y[0:len(y)>>1] = np.real(x)
    y[len(y)>>1: ] = np.imag(x)
    return y
    
def complexout(x):
    a = x[0:len(x)>>1]
    b = x[len(x)>>1:]
    return a + 1j*b

def gnlsewithgain(T,A,w0,gamma,betas,gain,fr,RT,fiberlength,nsaves,isgain=True):
    # From Supercontinuum Generation in Optical Fibers
    # Modified slightly to include gain   
    
    # Time Grid Parameters
    n = len(T)
    
    dT = np.double(T[1] - T[0])
    V = 2*np.pi*np.transpose(np.arange(-n/2,n/2, dtype=np.double))/(n*dT)      # Frequency grid
    
    fftc = fftw(int(n)) 
    
    # Calculate g from Gain via 10^(G(dB)/10) = exp(gain)
    alpha = -np.log(10**(np.double(gain)/10))
#    print 'gain: %f'%gain
#    print 'alpha: %f'%alpha
    
    # Dispersion
    B = np.zeros(n)
    for i in range(len(betas)):
        B = B + betas[i]/factorial(i+2)*V**(i+2)
    
    # Gain. To first approximation, we can import a gain spectrum

    calcGain = alpha * getGainsForGNLSE( (w0+V) * 1e12/ (2*np.pi), 'er7pm')

    if isgain == True:
        L = 1j*B - calcGain/2
    else:
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
#    print "Total Power Diagnostic"
#    print np.sum(np.abs(fftc.ifft(A)*np.exp(fiberlength*L))**2)
    
#    if isgain == True:
#        np.savetxt('pydump.csv',[getGainsForGNLSE( (w0+V) * 1e12/ (2*np.pi), 'er7pm'), np.abs(A)], delimiter=',')    
    
    temp_store = np.zeros(n*2, dtype=np.double)
    temp_store_cplx = np.zeros(n, dtype=np.complex128)
    IT = np.zeros(n, dtype=np.complex128)
    AT = np.zeros(n, dtype=np.complex128)
    RS = np.zeros(n, dtype=np.complex128)
    R = np.zeros(n, dtype=np.complex128)
    M = np.zeros(n, dtype=np.complex128)    
    def rhs(AW, z):      #rhs(z, AW)  
        temp_store_cplx.real = AW[:n]
        temp_store_cplx.imag = AW[n:]
        AT[:] = fftc.fft(temp_store_cplx*np.exp(L*z), copy=False)
        IT[:] = np.abs(AT)**2 # IT = intensity, time
        if (len(RT)==1 or abs(fr) < eps): # No Raman case
            M[:] = fftc.ifft(AT*IT, copy=False) # response function
        else:
            RS[:] = dT*fr*fftc.fft(fftc.ifft(IT, copy=False)*RW, copy=False) # Raman convolution
            M[:] = fftc.ifft(AT*((1-fr)*IT + RS), copy=False) # Response function
        R[:] = 1j*gamma*W*M*np.exp(-L*z) # Full RHS of propagation eqn        
        #report(z)
        
#        h5file = tables.open_file('pyout.h5',mode='w',title='Base')
#        gcols = h5file.create_group(h5file.root, 'nextdir', 'label')
#        h5file.create_array(gcols,'arrayset',R,'output of rhs function')
#        h5file.close()        
        temp_store[:n] = R.real[:]
        temp_store[n:] = R.imag[:]
        return temp_store
    
    #**** Define equation to impliment propagation equation
    

    # Set up and run ODE integrator
    multfactor = 10
    Z = np.linspace(0, fiberlength, nsaves*multfactor + 1)       # set up save points
    
    # Set integrator options
    solver = odespy.RK4(rhs,rtol=1e-5,atol=1e-12)
    solver.set_initial_condition(complexin(fftc.ifft(A))) #complexin
    
    AWout, Z = solver.solve(Z) # Run integrator    
    
    AW = AWout[::multfactor,:]
    Z = Z[::multfactor]
    
#    solver = ode(rhs)
#    solver.set_integrator("dop853", 
#                          first_step = fiberlength / 100000.0,
#                          max_step   = fiberlength/1000.0, nsteps = -1)
#    solver.set_initial_value(complexin(fftc.ifft(A)))
#    AWout = []    
#    Z = []
#    ncalls = 0
#    while solver.t < fiberlength:
#        print (solver.t/fiberlength)*100,'%'
#        #AWout.append(solver.integrate(solver.t, step = True))
#        solver.integrate(fiberlength, step = True)
#        Z.append(solver.t)
#        ncalls += 1
#        print ncalls
    
#    AW = np.array(AWout)
    
    # Process output
    AW = AW[:,:n] + 1j*AW[:,n:]

    AT = np.zeros(np.shape(AW), dtype=np.complex64)
    for i in range(len(Z)):
        AW[i,:] = AW[i,:]*np.exp(np.transpose(L)*Z[i]) # Change variables
        AT[i,:] = fftpack.ifftshift(fftc.fft(AW[i,:])) # Time domain
        AW[i,:] = fftpack.fftshift(AW[i,:])/dT # Scale
    W = V + w0 # Absolute frequency grid

    return [Z, AT, AW, W]