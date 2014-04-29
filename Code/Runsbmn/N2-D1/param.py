import math
import numpy as np
distance=1
beta2=1
N=2
mshape=0
chirp0=0
#N=.27
#mshape=0
#chirp0=0

#Simulation Parameters
nt=1024
Tmax=32
stepNum=round(20.*distance*N**2.)
if (stepNum<100):
    stepNum=400
deltaz=distance/stepNum
dtau=(2.*Tmax)/nt

##Arrays
#
#tau =dtau*np.arange(0,nt)
#omega = math.pi/Tmax*np.linspace(-nt/2,nt/2, num=nt)
#
#
##Input Field Profile
#
#if mshape == 0:
#    uu=np.cosh(tau)**-1*np.exp(-0.5j*chirp0*tau**2) #soliton
#else:
#    uu = np.exp(-0.5*(1+1j*chirp0)*tau**(2*mshape))
#
#
