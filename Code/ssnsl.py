import sys
import os
import numpy as np
import math
import param
import matplotlib.pyplot as plt
'''This program will do it's best to compute the split-step 
fourier transform to model nonlinear optical pulse propogation'''

if __name__ == '__main__':
    
    #Arrays

    tau =param.dtau*np.arange(-param.nt/2,param.nt/2)
    omega = math.pi/param.Tmax*np.linspace(-param.nt/2,param.nt/2, num=param.nt)


    #Input Field Profile

    if param.mshape == 0:
        uu=np.cosh(tau)**-1*np.exp(-0.5j*param.chirp0*tau**2) #soliton
    else:
        uu = np.exp(-0.5*(1+1j*param.chirp0)*tau**(2*param.mshape))

    
    #Plot Spectrum
    spec =np.fft.ifftshift( np.fft.fft(uu))
    freq = np.fft.ifftshift(np.fft.fftfreq(param.nt,d=param.dtau))
    f = plt.figure()
    plt.subplot(2,1,1)
    plt.plot(tau, np.absolute(uu)**2, 'ko-') 
    plt.title('Power of Input Pulse')
    plt.xlabel('Time')
    plt.subplot(2,1,2)
    plt.plot(np.fft.ifftshift(np.fft.fftfreq(param.nt,d=param.dtau)), np.absolute(spec)**2, 'ko-')
    plt.title('Power of Spectrum')
    plt.xlabel('Frequency')
    
    #Store Dispersion
    dispersion = np.exp(1j*0.5*param.beta2*freq**2*param.deltaz)
    hhz = 1j*N**2*param.deltaz

    #use symmetrized split step method

    temp =uu*np.exp(np.absolute(uu)**2*hhz/2)
    for x in range(param.stepNum):
        fTemp= np.fft.fft(temp)*dispersion
        uu = np.fft.ifft(fTemp)
        temp= uu*np.exp(np.absolute(uu)**2*hhz)
    uu = temp*np.exp(-np.absolute(uu)**2*hhz/2) #Final Field
    temp = np.fft.fftshift(np.fft.fft(uu)*(param.nt*dtau)/math.sqrt(2*math.pi)

    figure(2)
    subplot(2,1,1)
    plt.plot(tau, np.absolute(uu)**2, 'ko-')
    subplot(2,1,2)
    plt.plot(freq, np.absolute(temp)**2, 'ko-')
