# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 11:12:22 2014

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

fiber1 = Fiber()

steps = 500

centerwl = 1550.0
gamma = 2e-3
fiber_length = 10.0
T0 = 50e-3
D = 4 # ps / km / nm
beta2 = -2 * np.pi * fiber1.c / centerwl**2 * D
beta3 = 0.1
betas = [beta2, beta3]
P0 = abs(betas[0] * 1e-3) / gamma / T0**2
init = Pulse(n = 2**13)
init.gen_sech(P0, T0, centerwl)



fiber1.generate_fiber(fiber_length, centerwl, betas, gamma, 0, "ps^n/km")

evol = SSFM(disable_Raman = False, disable_self_steepening = False,
            local_error = 0.01, suppress_iteration = True)