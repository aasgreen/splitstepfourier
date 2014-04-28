# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 13:56:17 2014

@author: dim1
"""
import numpy as np
from scipy.misc import factorial
from scipy import constants
import matplotlib.pyplot as plt

def DTabulationToBetas(lambda0, DData, polyOrder, DDataIsFile = True):
    # Read in a tabulation of D vs Lambda
    # Expand about lambda0
    makePlots = 0
    if DDataIsFile:
        DTab = np.genfromtxt(DData,delimiter=',',skiprows=1)
    else:
        DTab = DData[:]
            
    # Units of D are ps/nm/km
    # Convert to s/m/m 
    DTab[:,1] = DTab[:,1] * 1e-12 * 1e9 * 1e-3
    c = constants.speed_of_light
    
    omegaAxis = 2*np.pi*c  / (DTab[:,0]*1e-9) - 2*np.pi*c  /(lambda0 * 1e-9)
    # Convert from D to beta via  beta2 = -D * lambda^2 / (2*pi*c) 
    
    betaTwo = -DTab[:,1] * (DTab[:,0]*1e-9)**2 / (2*np.pi*c) 
    # The units of beta2 for the GNLSE solver are ps^2/m; convert
    betaTwo = betaTwo * 1e24
    # Also convert angular frequency to rad/ps
    omegaAxis = omegaAxis * 1e-12 #  s/ps
    
    # How betas are interpreted in gnlse.m:
    #B=0;
    #for i=1:length(betas)
    #    B = B + betas(i)/factorial(i+1).*V.^(i+1);
    #end
    
    # Fit beta2 with high-order polynomial
    polyFitCo = np.polyfit(omegaAxis, betaTwo, polyOrder)
    
    Betas = polyFitCo[::-1]
    
    polyFit = np.zeros((len(omegaAxis),))   

    for i in range(len(Betas)):
        Betas[i] = Betas[i] * factorial(i)
        polyFit = polyFit + Betas[i] / factorial(i)*omegaAxis**i
    
    if makePlots == 1:
#        try:
#            set(0,'CurrentFigure',dispfig);
#        catch ME
#            dispfig = figure('WindowStyle', 'docked');
#        end        
        plt.plot(omegaAxis, betaTwo,'o')
        plt.plot(omegaAxis, polyFit)
        plt.show()
    
    return Betas