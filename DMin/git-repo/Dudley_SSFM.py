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


plt.close('all')

dz = 1e-3
steps = 500
range1 = np.arange(steps)

centerwl = 835.0
fiber_length = 0.15

init = Pulse(n = 2**13)
init.gen_sech(1e4, 28.4e-3, centerwl)

fiber1 = Fiber()
fiber1.load_from_db( fiber_length, 'dudley')

evol = SSFM(dz = 1e-6, local_error = 0.001)
y = np.zeros(steps)
AW = np.zeros((init.n, steps))
AT = np.copy(AW)

y, AW, AT, pulse1 = evol.Fiber_propagate(pulse_in = init, fiber = fiber1, 
                                         n_steps = steps)
                           
wl = 2 * np.pi * init.c / (init.W )

loWL = 400
hiWL = 1400
                         
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

D = fiber1.Beta2_to_D(init)
beta = fiber1.Beta2(init)
#
#plt.figure()
#plt.subplot(121)
#plt.plot(wl,D,'x')
#plt.xlim(400,1600)
#plt.ylim(-400,300)
#plt.xlabel('Wavelength (nm)')
#plt.ylabel('D (ps/nm/km)')
#plt.subplot(122)
#plt.plot(wl,beta*1000,'x')
#plt.xlim(400,1600)
#plt.ylim(-350,200)
#plt.xlabel('Wavelength (nm)')
#plt.ylabel(r'$\beta_2$ (ps$^2$/km)')

plt.figure()
plt.subplot(121)
plt.pcolormesh(xW, y, zW, vmin = mlIW - 40.0, vmax = mlIW)
plt.autoscale(tight=True)
plt.xlim([loWL, hiWL])
plt.xlabel('Wavelength (nm)')
plt.ylabel('Distance (m)')

plt.subplot(122)
plt.pcolormesh(xT, y, zT, vmin = mlIT - 40.0, vmax = mlIT)
plt.autoscale(tight=True)
plt.xlabel('Delay (ps)')
plt.ylabel('Distance (m)')

plt.show()
