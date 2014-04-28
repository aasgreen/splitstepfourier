# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 13:29:08 2014

@author: dim1
"""
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 12:42:30 2014

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

steps = 200

centerwl = 850.0
gamma = 1e-3
fiber_length = 1600.0
T0 = 1
beta = [-1, 0]
P0 = (abs(beta[0] * 1e-3) / gamma / T0**2) * 16

init = Pulse(n = 2**13)
init.gen_sech(P0, T0, centerwl)

fiber1 = Fiber()
fiber1.generate_fiber(fiber_length, centerwl, beta, gamma, 0, "ps^n/km")

evol = SSFM(disable_Raman = True, disable_self_steepening = True,
            local_error = 0.001, suppress_iteration = True)

y = np.zeros(steps)
AW = np.complex64(np.zeros((init.n, steps)))
AT = np.complex64(np.copy(AW))

y, AW, AT, pulse1, = evol.Fiber_propagate(pulse_in = init, fiber = fiber1, 
                                         n_steps = steps)
                                         
wl = 2 * np.pi * init.c / (init.W)

loWL = 820
hiWL = 870
                         
iis = np.logical_and(wl>loWL,wl<hiWL)

iisT = np.logical_and(init.T>-5,init.T<5)

xW = wl[iis]
xT = init.T[iisT]
zW_in = np.transpose(AW)[:,iis]
zT_in = np.transpose(AT)[:,iisT]
zW = 10*np.log10(np.abs(zW_in)**2)
zT = 10*np.log10(np.abs(zT_in)**2)
mlIW = np.max(zW)
mlIT = np.max(zT)

D = fiber1.Beta2_to_D(init)

x = (init.W - init.w0) / (2* np.pi) * T0
x2 = init.T / T0
b2 = beta[0] / 1e3 # in ps^2 / m
LD = T0**2 / abs(b2)
ynew = y / LD

plt.figure()
plt.subplot(121)
plt.pcolormesh(x2, ynew, 10*np.log10(np.abs(np.transpose(AT))**2),
               vmin = mlIW - 20.0, vmax = mlIW, cmap = plt.cm.gray)
plt.autoscale(tight=True)
plt.xlim([-4, 4])
plt.xlabel(r'($T / T_0)$')
plt.ylabel(r'Distance ($z/L_{NL})$')
plt.subplot(122)
plt.pcolormesh(x, ynew, 10*np.log10(np.abs(np.transpose(AW))**2),
               vmin = mlIW - 20.0, vmax = mlIW, cmap = plt.cm.gray)
plt.autoscale(tight=True)
plt.xlim([-4, 4])
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