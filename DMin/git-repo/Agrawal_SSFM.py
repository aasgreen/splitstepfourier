# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 15:39:12 2014

@author: dim1
"""

import numpy as np
import matplotlib.pyplot as plt
from SSFM import SSFM
from fiber import Fiber
from pulse import Pulse
#from fftw_transforms import fftcomputer as fftw
#from scipy import fftpack

plt.close('all')

steps = 100

centerwl = 1550.0
gamma = 2e-3



fiber_length = 2500.0
P0 = 10.0
T0 = 0.1

init = Pulse(n = 2**13)
init.gen_sech(P0, T0, centerwl)

fiber1 = Fiber() 
fiber1.generate_fiber(fiber_length,centerwl, [0,0], gamma, 0)

evol = SSFM(disable_Raman = True, disable_self_steepening = True,
            local_error = 0.1, suppress_iteration = True)


y = np.zeros(steps)
AW = np.complex64(np.zeros((init.n, steps)))
AT = np.complex64(np.copy(AW))



y, AW, AT, pulse1, = evol.Fiber_propagate(pulse_in = init, fiber = fiber1, 
                                         n_steps = steps)

m = 1
T = init.T
Leff = fiber_length
LNL = 1/(gamma * P0)
dw_T = (2*m/(T0 * 1e-12)) * (Leff/LNL) * (T/T0)**(2*m-1) * np.exp(-(T/T0)**(2*m))                             
                             
wl = 1e9 * 2 * np.pi * init.c_mks / (init.W * 1e12)

loWL = 1200
hiWL = 2000
                         
iis = np.logical_and(wl>loWL,wl<hiWL)

iisT = np.logical_and(init.T>-1,init.T<5)

xW = wl[iis]
xT = init.T[iisT]
zW_in = np.transpose(AW)[:,iis]
zT_in = np.transpose(AT)[:,iisT]
zW = 10*np.log10(np.abs(zW_in)**2)
zT = 10*np.log10(np.abs(zT_in)**2)
mlIW = np.max(zW)
mlIT = np.max(zT)

plt.figure()
plt.plot(T/T0,dw_T)
plt.xlim(-2,2)

x = (init.W - init.w0) / (2* np.pi) * T0

plt.figure()
#plt.plot(init.V, np.abs(AW[:,-1]))
plt.pcolormesh(x, y[0:-1] / LNL, 10 * np.log10(np.abs(np.transpose(AW))**2), vmin = mlIW - 40.0, vmax = mlIW,
               cmap = plt.cm.gray)
plt.xlim(-8, 8)
plt.xlabel(r'($\nu - \nu_0) \times T_0$')
plt.ylabel(r'Distance ($z/L_{NL})$')

#plt.figure()
#plt.subplot(121)
#plt.pcolormesh(xW, y, zW, vmin = mlIW - 40.0, vmax = mlIW)
#plt.autoscale(tight=True)
#plt.xlim([loWL, hiWL])
#plt.xlabel('Wavelength (nm)')
#plt.ylabel('Distance (m)')
#
#plt.subplot(122) 
#plt.pcolormesh(xT, y, zT, vmin = mlIT - 40.0, vmax = mlIT)
#plt.autoscale(tight=True)
#plt.xlabel('Delay (ps)')
#plt.ylabel('Distance (m)')

plt.show()