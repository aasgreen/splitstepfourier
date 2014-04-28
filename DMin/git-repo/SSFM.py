# -*- coding: latin1 -*-
#  Copyright (C) 2006-2010 João Luís Silva <jsilva@fc.up.pt>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2, or (at your option)
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#

import numpy as np

import matplotlib.pyplot as plt
from numpy.fft import fftshift, ifftshift
from pulse import Pulse
import pyfftw

#-----------------------------------------------------------------------
# Evolution of the pulse through the material
#-----------------------------------------------------------------------

###     Global variables    ###
#   USE_PYFFTW : True   - > use pyfftw
#                False  - > use numpy fft
#   On my home PC (Intel Q6600) PYFFTW is almost 2x faster
USE_PYFFTW = True


class SSFM:
    METHOD_SSFM,METHOD_RK4IP = range(2)

    
    def __init__(self,  local_error = 0.001, dz = 1e-5,
                 disable_Raman = False, disable_self_steepening = False,
                 suppress_iteration = True):
        self.iter = 0
        self.last_h = -1.0
        self.last_dir = 0.0
        self.eta = 5
        self.local_error = local_error
#        self.method = SSFM.METHOD_SSFM
        self.method = SSFM.METHOD_RK4IP
        self.disable_Raman = disable_Raman
        self.disable_self_steepening = disable_self_steepening
        self.f_R = 0.18
        
        self.tau_1 = 0.0122
        self.tau_2 = 0.0320
        self.dz = dz
        self.dz_min = 1e-12
        self.suppress_iteration = suppress_iteration



    def setup_fftw(self, pulse_in, fiber, output_power, raman_plots = False):
        ''' Call immediately before starting Propagate. This function does two
        things:\n
        1) it sets up byte aligned arrays for fftw\n
        2) it fftshifts betas, omegas, and the Raman response so that no further\n
            shifts are required during integration. This saves lots of time.'''

        
        self.n = pulse_in.n            
        
        self.fft_input    = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.ifft_input   = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.fft_input_2  = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.ifft_input_2 = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.A_I    = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.A2     = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.exp_D  = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.k1     = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.k2     = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.k3     = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.k4     = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.temp   = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')          
        self.Aw     = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')              
        self.A2w    = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.dA     = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.dA2    = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.R_A2   = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.dR_A2  = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')      
        self.omegas = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')        
        self.alpha  = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')                
        self.betas  = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')                
        self.A      = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')
        self.R      = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')
        self.R0     = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')       
        self.Af     = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')
        self.Ac     = pyfftw.n_byte_align_empty((self.n), 16, dtype='complex128')        

        self.omegas[:]      =  pulse_in.V 
        self.betas[:]       =  fiber.get_betas(pulse_in)
        self.alpha[:]       = -np.log(10**(fiber.get_gain(pulse_in, output_power))/10)
        self.gamma          = fiber.gamma
        self.w0             = pulse_in.w0


        


        # To be double sure that there are no problems, also make 2 copies of
        # the FFT objects. This lets us nest ifft_2 around a function using ifft
        # without worrying about potential problems.
        self.ifft   = pyfftw.builders.fft(self.ifft_input)
        self.ifft_2   = pyfftw.builders.fft(self.ifft_input_2)
        self.fft    = pyfftw.builders.ifft(self.fft_input)
        self.fft_2    = pyfftw.builders.ifft(self.fft_input_2)

        if not self.disable_Raman:
            self.CalculateRamanResponseFT(pulse_in)
            if raman_plots:
                plt.subplot(221)
                plt.plot(self.omegas/(2*np.pi), np.abs(self.R-(1-self.f_R)),'b.')                
                plt.title('Abs[R(w)]')
                plt.xlabel('THz')                
                plt.subplot(222)
                plt.plot(self.omegas/(2*np.pi), np.unwrap(np.angle(self.R-(1-self.f_R))),'b.')
                plt.title('Angle[R(w)]')
                plt.xlabel('THz')
                plt.subplot(223)
                plt.plot(pulse_in.T, ifftshift(abs(self.IFFT_t(self.R))))
                plt.title('Abs[R[t]]')
                plt.ylim([0,300])
                plt.subplot(224)
                plt.plot(self.omegas/(2*np.pi), abs(self.FFT_t(self.A)))
                plt.title('Abs[A[w]]')
                plt.xlabel('THz')
                plt.show()
        self.A[:]       = fftshift(pulse_in.A)
        self.omegas[:]  = fftshift(self.omegas)
        self.betas[:]   = fftshift(self.betas)
        self.alpha[:]   = fftshift(self.alpha)
        self.R[:]       = fftshift(self.R)
        self.R0[:]      = fftshift(self.R0)
            
    #-----------------------------------------------------------------------
    # Calculates the Fourier Transform of R(T). See pg 49 of G. P. Agrawal's 
    # "Nonlinear fiber optics"  for details 
    #-----------------------------------------------------------------------
    def CalculateRamanResponseFT(self, pulse):
        ''' Calculate Raman response in frequency domain. Two versions are
            available: the first is the LaserFOAM one, which directly calculates
            R[w]. The second is Dudley-style, which calculates R[t] and then
            FFTs. Note that the use of fftshifts is critical here (FFT_t_shift)
            as is the factor of pulse_width.'''
        # Laserfoam raman function.
        TAU1 = self.tau_1
        TAU2 = self.tau_2
        F_R = self.f_R        
        C = (TAU1**2+TAU2**2)/(TAU1*TAU2**2)        
        for i in range(pulse.n):
            omega = self.omegas[i]
            H_R = C*TAU1*TAU2**2 / \
                  (TAU1**2 + TAU2**2 - 2j*omega*TAU1**2*TAU2 - TAU1**2*TAU2**2*omega**2)
            self.R0[i] = (1.0-F_R) + (F_R * H_R)
        
        # More conventional way of generating this, via Dudley    
        tau1 = self.tau_1
        tau2 = self.tau_2
        T = pulse.T
        RT     = np.zeros(pulse.n, dtype = 'complex128')
        RT     = (tau1**2 + tau2**2) /( tau1 * tau2**2 )*\
                    np.exp(-T/tau2)*np.sin(T/tau1)
        RT[0:pulse.n>>1]=0
        RT[:]     = RT / np.trapz(RT, T)
        #H_R    = pulse.dT*pulse.n*self.FFT_t(fftshift(RT))
        self.R[:]    = ((1.0-F_R) + pulse.twidth*self.FFT_t_shift(F_R * RT))


    #-----------------------------------------------------------------------
    # Advances the current position by delta_z using an adaptive spatial
    # step algorithm.
    # See O.V. Sinkin et al, J. Lightwave Tech. 21, 61 (2003)
    # dir: 1 - Forward propagation
    #     -1 - Inverse propagation
    #-----------------------------------------------------------------------
    def Propagate(self,delta_z, dir=1):        
#        print "Propagate: delta_z",delta_z
        dist = delta_z
        dz = self.dz        

        self.last_h = -1.0  #Force an update of exp_D
        force_last_dz = False
        factor = 2**(1.0/self.eta)

        if (2.0*dz > dist):
            dz = dist/2.0

        while dist>0.0:
            self.Ac[:] = self.A
            self.Af[:] = self.A
            self.Ac[:] = self.Advance(self.Ac,2.0*dz,dir)
            self.Af[:] = self.Advance(self.Af,dz,dir)
            self.Af[:] = self.Advance(self.Af,dz,dir)

            #delta = |Af - Ac| / |Af| 
            delta = self.CalculateLocalError()

            old_dz = dz
            new_dz = dz
            if not self.suppress_iteration:
                print "iteration:",self.iter,"dz:",dz,"distance:", dist, \
                      "local error", delta                


            if delta > 2.0*self.local_error:
                # Discard the solution, decrease step

                new_dz = dz/2.0
                if new_dz >= self.dz_min:
                    dz = new_dz                                     
                    # discard current step
                    continue
                else:
                    # accept step after all
                    pass
            elif (delta >= self.local_error) and (delta<=2.0*self.local_error): 
                # Keep solution, decrease step
                new_dz = dz / factor
                if new_dz >= self.dz_min:
                    dz = new_dz
                else:
                    pass
#                    printf("[%d] limited a step decrease, h = %g, delta = %g\n",self.iter,dz,delta)
            elif (delta >= (0.5*self.local_error)) and (delta<=self.local_error):
                # keep the step
                new_dz = new_dz
            else:     # delta < local_error/2
                # Step too small
                new_dz = dz * factor
                dz = new_dz


            if self.eta==3:
                self.A[:] = (4.0/3.0) * self.Af -(1.0/3.0) * self.Ac
            elif self.eta==5:
                self.A[:] = (16.0/15.0) * self.Af -(1.0/15.0) * self.Ac
            else:
                p = 2**(self.eta-1.0)
                self.A[:] = (p/(p-1.0)) * self.Af -(1.0/(p-1.0)) * self.Ac

            dist -= 2.0*old_dz
            self.iter += 1


#            printf("-> [%d] dz = %g dist = %g (old_dz = %g) (z = %g)\n",n,dz,dist,old_dz,self.z)
            if (2.0*dz > dist) and (dist>2.0*self.dz_min):
                force_last_dz = True
                return_dz = dz
                dz = dist/2.0
#                printf("[%d] dz = %f\n",self.iter,dz)

        if force_last_dz:
            dz = return_dz
        self.dz = dz
        return
 
    def Advance(self,A,dz,dir):
        if self.method == SSFM.METHOD_SSFM:
            if dir==1:
                A = self.LinearStep(A,dz,dir)
                return np.exp(dz*dir*self.NonlinearOperator(A))*A
            else:
                A = np.exp(dz*dir*self.NonlinearOperator(A))*A
                return self.LinearStep(A,dz,dir)
        elif self.method == SSFM.METHOD_RK4IP:
            return self.RK4IP(A,dz,dir)

    def LinearStep(self,A,h,dir):
        if h!=self.last_h or dir!=self.last_dir:
            self.Calculate_expD(h,dir)
            self.last_h = h
            self.last_dir = dir                        
        return self.IFFT_t(self.exp_D * self.FFT_t(A))

    def Deriv(self,Aw):
        """Calculate the temporal derivative using FFT. \n\n MODIFIED from 
        LaserFOAM original code, now input is in frequency space, output is 
        temporal derivative. This should save a few FFTs per iteration."""
        return self.IFFT_t(-1.0j*self.omegas * Aw)

    def NonlinearOperator(self,A):
        if self.disable_Raman:
            if self.disable_self_steepening:
                return 1j*self.gamma*np.abs(A)**2
                
            self.Aw[:] = self.FFT_t(A)
            self.dA[:] = self.Deriv(self.A)
            
            return 1j*self.gamma*np.abs(A)**2 - \
                   (self.gamma/self.w0)*(2.0*self.dA*A.conj() + A*self.dA.conj())
        else:

            self.A2[:]  = np.abs(A)**2   
            self.A2w[:] = self.FFT_t(self.A2)
           
            if self.disable_self_steepening:
                return 1j*self.gamma*self.IFFT_t(self.R*self.A2w)

            self.Aw[:] = self.FFT_t(A)          
            self.R_A2[:]    = self.IFFT_t(self.R*self.A2w)
            self.dA[:]      = self.Deriv(self.Aw)
            self.dA2[:]     = self.Deriv(self.A2w)        
            self.dR_A2[:]   = self.IFFT_t(self.R*self.FFT_t(self.dA2))
            
            return 1j*self.gamma*self.R_A2 - (self.gamma/self.w0)* \
                   (self.dR_A2 + np.where(np.abs(A)>1.0E-5,self.dA*self.R_A2/A,0.0))


    def RK4IP(self,A,h,dir):
        """Fourth-order Runge-Kutta in the interaction picture.
           J. Hult, J. Lightwave Tech. 25, 3770 (2007)."""                   
        self.A_I[:] = self.LinearStep(A,h,dir)  #Side effect: Rely on LinearStep to recalculate self.exp_D for h/2 and direction dir                
        self.k1[:] = self.IFFT_t_2(self.exp_D*self.FFT_t_2(h*dir*self.NonlinearOperator(A)*A))

        self.k2[:] = h * dir * self.NonlinearOperator(self.A_I + self.k1/2.0)*\
                        (self.A_I + self.k1/2.0)
        self.k3[:] = h * dir * self.NonlinearOperator(self.A_I + self.k2/2.0)*\
                        (self.A_I + self.k2/2.0)        
        self.temp[:] = self.IFFT_t_2(self.exp_D*self.FFT_t_2(self.A_I+self.k3))
        self.k4 = h * dir * self.NonlinearOperator(self.temp)*self.temp
        if not self.suppress_iteration:
            print "ks: ",np.sum(abs(self.k1)),np.sum(abs(self.k2)),\
                    np.sum(abs(self.k3)),np.sum(abs(self.k2))
        
        return self.IFFT_t_2(self.exp_D * self.FFT_t_2(self.A_I + self.k1/6.0 +\
                self.k2/3.0 + self.k3/3.0)) + self.k4/6.0

    def Calculate_expD(self,h,dir):        
        self.exp_D[:] = np.exp(dir*h*0.5*(1j*self.betas-self.alpha/2.0))

    def Fiber_propagate(self, pulse_in, fiber, n_steps, output_power = 0):

        n_steps = int(n_steps)
        
        # Copy parameters from pulse and fiber into class-wide variables
                         
        z_positions = np.linspace(0, fiber.length, n_steps + 1)
        delta_z = z_positions[1] - z_positions[0]

        AW = np.complex64(np.zeros((pulse_in.n, n_steps)))
        AT = np.complex64(np.copy(AW))
        
        print "Pulse energy before", fiber.fibertype,":", \
              1e9 * pulse_in.calc_epp(), 'nJ'          

        pulse_out = Pulse()
        pulse_out.clone_pulse(pulse_in)

        self.setup_fftw(pulse_in, fiber, output_power)
        
        for i in range(n_steps):                        
            print "steps:", i, "totaldist:", fiber.length * (1 - np.float(i)/n_steps)
            self.Propagate(delta_z)
            AW[:,i] = ifftshift(self.FFT_t_2(self.A))
            AT[:,i] = ifftshift(self.A)
    
        pulse_out.A[:] = ifftshift(self.A)

        print "Pulse energy after", fiber.fibertype,":", \
              1e9 * pulse_out.calc_epp(), 'nJ'
#        print "alpha out:",self.alpha
        return z_positions, AW, AT, pulse_out
        
    def test_raman(self, pulse_in, fiber, output_power = 0): 
        ''' Function for testing Raman response function. This plots Generates 
            R[w] and makes plots, but does not actually integrate anything.'''
        # Copy parameters from pulse and fiber into class-wide variables                    
        print "Pulse energy before", fiber.fibertype,":", \
              1e9 * pulse_in.calc_epp(), 'nJ'          
        self.setup_fftw(pulse_in, fiber, output_power, raman_plots = True)

    ### Lots of boring FFT code from here on out.
    def FFT_t(self, A):
        if USE_PYFFTW:
            self.fft_input[:] = A
            return self.fft()
        else:
            return np.fft.ifft(A)
    def IFFT_t(self, A):
        if USE_PYFFTW:
            self.ifft_input[:] = A
            return self.ifft()
        else:
            return np.fft.fft(A)
    def FFT_t_shift(self, A):
        if USE_PYFFTW:
            self.fft_input[:] = fftshift(A)
            return ifftshift(self.fft())
        else:
            return np.fft.ifft(A)
    def IFFT_t_shift(self, A):
        if USE_PYFFTW:
            self.ifft_input[:] = fftshift(A)
            return ifftshift(self.ifft())
        else:
            return np.fft.fft(A)            
    def FFT_t_2(self, A):
        if USE_PYFFTW:
            self.fft_input_2[:] = A
            return self.fft_2()
        else:
            return np.fft.ifft(A)
    def IFFT_t_2(self, A):
        if USE_PYFFTW:
            self.ifft_input_2[:] = A
            return self.ifft_2()
        else:
            return np.fft.fft(A)
            

    #-----------------------------------------------------------------------
    # Calculates the relative local error.
    # See O.V. Sinkin et al, J. Lightwave Tech. 21, 61 (2003)
    # Returns |Af - Ac| / |Af|
    #-----------------------------------------------------------------------
    def CalculateLocalError(self):
        denom = np.linalg.norm(self.Af)
        if denom != 0.0:
            return np.linalg.norm(self.Af-self.Ac)/np.linalg.norm(self.Af)
        else:
            return np.linalg.norm(self.Af-self.Ac)
            