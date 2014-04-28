# -*- coding: utf-8 -*-
"""
Created on Mon Apr 07 16:41:22 2014

@author: dim1
"""
import numpy as np
import matplotlib.pyplot as plt
from SSFM import SSFM
from fiber import fiber
from pulse import pulse
from fftw_transforms import fftcomputer as fftw
#from scipy import fftpack

plt.close('all')

dz = 1e-3
steps = 100
range1 = np.arange(steps)
range2 = range1 + steps
range3 = range2 + steps
maxpower = 0.4
centerwl = 1560.0


PANDAinitial = 1.05 - 0.5
EDF = 1.30
PANDAfinal = 0.5 # 2.0

init = pulse(power_in = 0.020)

init_in = np.copy(init.A)

evol = SSFM(init, dz = 1e-3)
fftc = fftw(init.n)

y = np.zeros(3 * steps)
AW = np.zeros((init.n, 3 * steps))
AT = np.copy(AW)

fiber1 = fiber(length = PANDAinitial, center_wavelength = init.wavelength,
               fibertype = 'PANDA1550') 
fiber2 = fiber(length = EDF, center_wavelength = init.wavelength,
               fibertype = 'EDF07PM') 
fiber3 = fiber(length = PANDAfinal, center_wavelength = init.wavelength,
               fibertype = 'PANDA1550') 
    
y[range1], AW[:, range1], AT[:, range1], pulse1 = evol.Fiber_propagate(
                                                  pulse_in = init,
                                                  fiber = fiber1,
                                                  n_steps = steps)
y[range2], AW[:, range2], AT[:, range2], pulse2 = evol.Fiber_propagate(
                                                  pulse_in = pulse1,
                                                  fiber = fiber2,
                                                  n_steps = steps,
                                                  output_power = maxpower)
y[range3], AW[:, range3], AT[:, range3], pulse3 = evol.Fiber_propagate(
                                                  pulse_in = pulse2,
                                                  fiber = fiber3,
                                                  n_steps = steps)                                                        

y[range2] = y[range2] + PANDAinitial
y[range3] = y[range3] + PANDAinitial + EDF


def delta_wl(wl,freq,delta_freq):
    return delta_freq / freq * wl

wl = centerwl + delta_wl(centerwl, (init.c*1e3) / (centerwl * 1e-9),
                         init.V / (2 * np.pi) * 1e12)
                         
iis = np.logical_and(wl>init.lowWL,wl<init.highWL)
xW = wl[iis]
xT = init.T#[iis]
zW_in = np.transpose(AW)[:,iis]
zT_in = np.transpose(AT)#[:,iis]
zW = 10*np.log10(np.abs(zW_in)**2)
zT = 10*np.log10(np.abs(zT_in)**2)
mlIW = np.max(zW)
mlIT = np.max(zT)

plt.figure()
plt.subplot(121)
plt.pcolormesh(xW, y, zW, vmin = mlIW - 40.0, vmax = mlIW)
plt.autoscale(tight=True)
plt.xlim([init.lowWL, init.highWL])
plt.xlabel('Wavelength (nm)')
plt.ylabel('Distance (m)')

plt.subplot(122)
plt.pcolormesh(xT, y, zT, vmin = mlIT - 40.0, vmax = mlIT)
plt.autoscale(tight=True)
plt.xlabel('Delay (ps)')
plt.ylabel('Distance (m)')
#plt.xlim([-5,5])

plt.show()