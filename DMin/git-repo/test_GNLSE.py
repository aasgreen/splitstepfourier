# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 11:01:15 2014

@author: dim1

adapting Gabe's MATLAB code for Python
"""

# test the NLSE propagator

import os

os.chdir('O:\\OFM\\Maser\\PM for OPO and DFG\\')

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from fftw_transforms import fftcomputer as fftw
from scipy import fftpack
from gnslewithgain import gnlsewithgain
from DTabulationToBetas import DTabulationToBetas

plt.close('all')

c = 299792458*1e9/1e12 # c in nm/ps
n = 2**15 # n points

fftc = fftw(n)

# For simulation of the actual comb modes,
# we will want to set up the grid in frequency space,
# fixing the number of comb modes and their seperation

# The GNLSE code starts from the time-domain representation, so
# we need to calculate the time-domain parameters

#dF      = 250 * 1e6 * 1e-12 # 1/ps
#dF      = 50000 * 1e6 * 1e-12 # 1/ps
#dOmega  = 2*pi*dF # rad/ps
dT      = 0.001
dTSeconds = dT * 1e-12
twidth  = dT*(n/2.0)

# Set the average power coupled into the nonlinear fiber
powerIn = 0.033 # Watts
frep = 100e6
epp = powerIn / frep
wavelength = 1560.0
lowWL = 1100.0
highWL = 1900.0

w0 = (2.0*np.pi*c)/wavelength # reference angular frequency (2pi THz)
T = np.linspace(-twidth/2,twidth/2,n) # time grid
V = 2*np.pi*np.transpose(np.arange(-n/2,n/2))/(n*dT) # Frequency grid
W = V+w0 # Absolute frequency grid

print 'The time width is',twidth,'ps'

# Now load up data from the cavity simulation
# Format is frequency (Hz) | Intensityhttp://www.sportsonearth.com/article/67226344 (arb) | phase (radians)

baseDir='O:\\OFM\\Maser\\PM for OPO and DFG'
pulseData = np.genfromtxt(str(baseDir+'\\pulseEnteringEDF.csv'),delimiter=',')

mval = max(pulseData[:,1])
peakInd = np.argmax(pulseData[:,1])
intensityFunc = interp1d(pulseData[:,0], (pulseData[:,0]<3e17/1630)*pulseData[:,1], kind='linear',bounds_error=False,fill_value=0)
phaseFunc = interp1d(pulseData[:,0], pulseData[:,2], kind='linear', bounds_error=False,fill_value=0)

intensityInterp = intensityFunc(V*1e12/(2*np.pi))
phaseInterp = phaseFunc(V*1e12/(2*np.pi))

# Pad pulseData so that we have n points
padN = n-len(pulseData)
print 'Padding extra',padN,'points'

pulseData = np.vstack((np.zeros((np.floor(padN/2),3)), pulseData, np.zeros((np.ceil(padN/2),3))))
 
# Add random noise to the spectrum
noiseLevel = -560.0 # db from peak
freqNoise = (np.sqrt(max(pulseData[:,1])) * 10**(noiseLevel/10)*np.exp(1j*2*np.pi*np.random.rand(n,1)))[:,0]

# Frequency domain
eScale = np.sqrt(epp / ( dTSeconds*sum(abs(fftc.fft(np.sqrt(intensityInterp)))**2)))
print eScale
print '..'

# With real chirp
EfieldIn =  (freqNoise+eScale *np.sqrt(intensityInterp) * np.exp(1j*phaseInterp))
# Without real chirp
#EfieldIn =  (freqNoise+eScale *np.sqrt(intensityInterp))

# Time domain
A = fftpack.ifftshift(fftc.fft(fftpack.fftshift(EfieldIn)))
#A = A * sqrt(epp / sum(dTSeconds .* abs(A).*abs(A)));

plt.figure()
plt.subplot(411)
plt.plot(2*np.pi*c/W, abs(fftpack.fftshift(fftc.ifft(A)))**2,'.')
plt.xlim([1400,1700])

plt.subplot(412)
plt.plot(T,abs(A)**2)
plt.xlim([-1,1])

#plot(T,phase(A),'r')
print dT
print max(abs(fftc.ifft(A))**2)/dT

print 'Max power is',max(abs(A)**2),'W'
print 'Pulse energy is',1e9*sum(dTSeconds * abs(A)**2),'nJ'

plt.subplot(413)

iis = np.logical_and(2*np.pi*c/W>1500,2*np.pi*c/W<1600)

EFA = fftpack.fftshift(fftc.ifft(fftpack.ifftshift(A)))

plt.plot(2*np.pi*c/W[iis],np.angle(EFA[iis])-max(np.angle(EFA[iis])),'bo')
plt.plot(2*np.pi*c/W[iis],np.angle(EfieldIn[iis]),'ro')
plt.xlim([1500,1600])

plt.subplot(414)
plt.plot(2*np.pi*c/W,phaseInterp,'bo')
plt.xlim([1500,1600])

USE_SIMPLE_PULSE=0

c = 299792458*1e9/1e12 # c in nm/ps

#USE_SIMPLE_PULSE=1

if USE_SIMPLE_PULSE == 1:
    n = 2**14 # n points
    twidth = 25.0 # Time window in ps 
    wavelength = 1040.0 # reference wavelength (nm)
    T = np.linspace(-twidth/2,twidth/2,n) # time grid
    lowWL   = 500.0 # low and high wavelengths
    highWL  = 1700.0
    # Input pulse
    power   = 1.0 # peak power (W)
    t0      = 0.030 # input pulse duration (ps)
    A       = np.sqrt(power)/np.cosh(T/t0) # input field (W^0.5)
    
    energyPP = 4e-9 # J
    
    totalEnergy = sum(abs(A)**2)*twidth*1e-12/n # Figure out what the total `intensity' is
    A = A * np.sqrt(energyPP / totalEnergy)


w0 = (2.0*np.pi*c)/wavelength # reference angular frequency (2pi THz)


# Fiber parameters
# Betas = beta2, beta3, ...
#betas   = [-11.830e-3, 8.1038e-5, -9.5205e-8, 2.0737e-10,...
#        -5.3943e-13, 1.3486e-15, -2.5495e-18, 3.0524e-21,...
#        -1.7140e-24];
betas = [32e-3] # beta2 for PM EDF

baseDir='O:\\OFM\\Maser\\PM for OPO and DFG\\'
fiberDispPrefix = str(baseDir)

betasPANDA   = DTabulationToBetas(wavelength, str(fiberDispPrefix+'PANDA1550.csv'),2)
    #betas   = DTabulationToBetas(wavelength, strcat(fiberDispPrefix,'NL1050_ZERO.csv'),2);
    #betas   = DTabulationToBetas(wavelength, strcat(fiberDispPrefix,'SC-3.7-975.csv'),3);
    #betas   = DTabulationToBetas(wavelength, strcat(fiberDispPrefix,'NL-1050-NEG-1.csv'),2);
    #betas   = DTabulationToBetas(wavelength, strcat(fiberDispPrefix,'NL-PM-750.csv'),3);

# Gain fiber parameters
gamma = 5e-3 # 1/W/m
gain = 14.0 # dB/m
fiberlength = 1.0 # fiber length (m)

# Compression fiber parameters
gammaPANDA  = 1.1e-3 # 1/W/m
gainPANDA   = 0.0 # dB/m
fiberlengthPANDA = 0.47

# Raman response
fr     = 0.18 # fractional Raman contribution
tau1   = 0.0122 
tau2   = 0.032
RT     = (tau1**2 + tau2**2) /tau1 / tau2**2*np.exp(-T/tau2)*np.sin(T/tau1)
RT[T<0]=0 # Raman response only after impulse
RT     = RT / np.trapz(T, RT) # Normalize to unit integral
 
# simulation parameters
nsaves     = 50.0
 
print 'Running simulation.'
# Run simulation

# Run simulation for amplifier
#[Z, AT, AW, W] = gnlsewithgain(T,fftpack.fftshift(A),w0,gamma,betas,gain,fr,RT,fiberlength,nsaves,isgain=True)
#AFinal = AT[len(AT[:,0]) - 1,:]

# Run simulation for compression fiber
[ZC, ATC, AWC, WC] = gnlsewithgain(T,fftpack.fftshift(A),w0,gammaPANDA,betasPANDA,gainPANDA,fr,RT,fiberlengthPANDA,nsaves,isgain=False)

print 'Pulse energy before compressor is ',1e9*sum(dTSeconds*abs(A)**2),' nJ'
#print 'Pulse energy after amplifier is ',1e9*sum(dTSeconds*abs(AFinal)**2),' nJ'
print 'Pulse energy after compressor is ',1e9*sum(dTSeconds*\
    abs(ATC[len(ATC[:,0]) - 1,:])**2),' nJ'

Z =  ZC#np.hstack((Z, max(Z)+ZC))
AT =  ATC#np.vstack((AT, ATC))
AW =  AWC#np.vstack((AW, AWC))

# Plot output
IW     = abs(AW)**2 # spectral intensity
lIW     = 10*np.log10(abs(AW)**2) # log scale of spectral intensity
mlIW    = np.max(lIW) # max, for scaling
WL      = 2*np.pi*c/W # wavelength grid
iis     = np.logical_and(WL>lowWL,WL<highWL)

if n<2**19:
    plt.figure()
    plt.subplot(121)
    plt.pcolormesh(WL[iis], Z, lIW[:,iis],vmin=mlIW-40.0, vmax=mlIW)
    plt.autoscale(tight=True)
    #plt.caxis([mlIW-40.0, mlIW])
    plt.xlim([lowWL, highWL])
    #plt.shading interp;
    plt.xlabel('Wavelength / nm')
    plt.ylabel('Distance / m')

    lIT = 10*np.log10(abs(AT)**2) # log scale temporal intensity
    mlIT = np.max(lIT) # max, for scaling

    plt.subplot(122)
    plt.pcolormesh(T, Z, lIT,vmin=mlIT-40.0, vmax=mlIT)
    plt.autoscale(tight=True)
    #plt.caxis([mlIT-40.0, mlIT])
    #shading interp;
    plt.xlabel('Delay / ps')
    plt.ylabel('Distance / m')
    plt.xlim([-5,5])

plt.figure() 

plt.subplot(211)
endIdx = len(lIW[:,1])
plt.plot(WL[iis],IW[endIdx - 1,iis]/sum(IW[endIdx - 1,iis]),'+')
plt.plot(WL[iis],IW[0,iis]/sum(IW[0,iis]),'or')
plt.xlabel('Wavelength')
plt.ylabel('Intensity (dB)')

plt.subplot(212)
plt.plot(WL[iis], np.unwrap(np.angle(AW[0,iis])),'bo')
plt.plot(WL[iis], np.unwrap(np.angle(AW[endIdx - 1,iis])),'ro')
plt.xlabel('Wavelength')
plt.ylabel('Phase (rad)')

# Plot pulse in time domain at z closest to specified z0. For looking at
# the instantaneous pulse shape during compression.
z0 = 1.47 # m
err = min(abs(Z-z0))
index = np.argmin(abs(Z-z0))

plt.figure()
err = max((abs(AT[index,:])**2)[::-1])
peakIdx = np.argmax((abs(AT[index,:])**2)[::-1])

plt.plot(1000*(T - T[peakIdx]), (abs(AT[index, :])**2)[::-1])
plt.xlim([-500,500])
#np.savetxt('.\\GNLSEModeledCompressedPulseNoGainNonlinearity.csv',np.transpose([1000*(T - T[peakIdx - 1]), (abs(AT[index, :]**2)[::-1])]),delimiter=',')
    
#dlmwrite('C:\Users\Gabe.Gabe-PC\SkyDrive\Documents\Analysis\2013\09 September\Unsorted\GNLSEModeledCompressedPulseNoGainNonlinearity.csv',    transpose([1000*(T - T(peakIdx)); wrev(abs(AT(index, :)).^2)]),',')

#Movie

#not fully converted to Python!

#plt.figure()
#plt.ylim([-110,30])
#set(gca,'NextPlot','replacechildren');
## Preallocate the struct array for the struct returned by getframe
#F(len(lIW[:,1])) = struct('cdata',[],'colormap',[]);
## Record the movie
#for j in range(len(lIW[:,1])):   
#    plot(WL(iis),lIW(j,iis),'+')
#    F(j) = getframe;
#    print '.'


#movie(F, 10)
  
# Use results to calculate shift of comb mode-centers due to side-mode amplification.
# Calculation needs to know the index of a single main-mode, which was calculated 
# earlier and stored in peakInd. It als needs the filtering factor
# (fcavity/frep)
  
filterFactor = 30.0 * 4.0 # 30 GHz cavity, frep = 0.25 GHz
idxOffset = filterFactor - np.mod(peakInd, filterFactor)-1
startIdx = 1
stopIdx = int(np.floor( (len(AW[0,:])-idxOffset) / filterFactor) - 3)

IW = abs(AW**2)

fractionalShift = np.zeros((stopIdx-startIdx, 2))

frequencyAxis = (c / WL) * 1e12 # frequency axis, in Hz

endOfFiberIdx  = len(IW[:,1])

for i in range(startIdx,stopIdx):
    j = i - startIdx
    fSubSet = frequencyAxis[idxOffset+i*filterFactor - (filterFactor / 4 - 1) \
    : idxOffset+i*filterFactor + (filterFactor / 4 - 1)]
    iSubSet = IW[endOfFiberIdx - 1, idxOffset+i*filterFactor - (filterFactor / 4 -1 ) \
    : idxOffset+i*filterFactor + (filterFactor / 4 - 1)]
    weightedCenter = sum(fSubSet * np.transpose(iSubSet)) / sum (iSubSet)
    fractionalShift[j,0] = frequencyAxis[i*filterFactor]
    fractionalShift[j,1] = (weightedCenter -frequencyAxis[idxOffset+i*filterFactor]) \
    / frequencyAxis[idxOffset+i*filterFactor]            
#    if i == np.mean([startIdx, stopIdx]):
#        plt.figure()
#        plt.plot(fSubSet, 10*log10(iSubSet),'.')     

#plt.figure() 
#plt.plot(3e17 / fractionalShift[:,0], 10*np.log10(abs(fractionalShift[:,1])))

plt.show()