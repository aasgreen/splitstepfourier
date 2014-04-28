# -*- coding: utf-8 -*-
"""
Created on Thu Apr 03 09:38:19 2014

@author: dim1
"""
import numpy as np
from scipy.interpolate import interp1d
from scipy import constants 

from gnlse_ffts import FFT_t, IFFT_t

class Pulse:
    """Class which carried all information about the light field. Includes
    functions for generating standard pulse shapes and also loading complex
    spectra from FROG data. Initialize sets number of points in grid.
    Pulse repetition frequency frep is required to convert from pulse energy to
    average power."""
    # Constants
    c = constants.speed_of_light*1e9/1e12 # c in nm/ps
    c_mks = constants.speed_of_light # m/s
    
    def __init__(self, frep = 100e6, n = 2**12):
        self.n = n # n points
        
        # Set the average power coupled into the nonlinear fiber
        self.frep = frep        
        
    def gen_sech(self, peak_power, T0, center_wavelength,
                 time_window = 10.0, GDD = 0, TOD = 0, chirp2 = 0, chirp3 = 0):
        """Generate sech pulse A(t) = sqrt(peak_power [W]) * sech(t/T0 [ps])
        centered at wavelength center_wavelength (nm).
        time_window (ps) sets temporal grid size. Optional GDD and TOD are
        in ps^2 and ps^3."""
        # The free parameters are time_window, wavelength, and grid size.
        self.twidth      = time_window          # Time window in ps
        self.center_wl  = center_wavelength    # reference wavelength (nm)              

        self.setup_grids()
        
        ### Generate pulse
        self.A              = np.sqrt(peak_power)/np.cosh(self.T/T0)
        self.chirp_pulse_W(GDD, TOD)
        self.chirp_pulse_T(chirp2, chirp3, T0)            
        
    def gen_gaussian(self, peak_power, T0, center_wavelength,
                 time_window = 10.0, GDD = 0, TOD = 0, chirp2 = 0, chirp3 = 0):
        """Generate Gaussian pulse A(t) = sqrt(peak_power[W]) * 
            exp( -(t/T0 [ps])^2 / 2 ) centered at wavelength 
            center_wavelength (nm). time_window (ps) sets temporal grid
            size. Optional GDD and TOD are in ps^2 and ps^3."""
        # The free parameters are time_window, wavelength, and grid size.
        self.twidth      = time_window          # Time window in ps
        self.center_wl  = center_wavelength    # reference wavelength (nm)                        

        self.setup_grids()

        self.A           = np.sqrt(peak_power) * np.exp(-self.T**2/(2 * T0**2)) # input field (W^0.5)            
        self.chirp_pulse_W(GDD, TOD)
        self.chirp_pulse_T(chirp2, chirp3, T0)
        
    def gen_frog(self, time_window, center_wavelength, power, 
                 power_is_epp = False,
                 fileloc = '250MHz_comb_Speck.dat', flip_phase = True):
        """Generate pulse from FROG data. Grid is centered at wavelength
        center_wavelength (nm), but pulse properties are loaded from data
        file. If flip_phase is true, all phase is multiplied by -1 [useful
        for correcting direction of time ambiguity]. time_window (ps) sets 
        temporal grid size. 
        
        power sets the pulse energy:
        if power_is_epp is True  then the number is pulse energy [J] 
        if power_is_epp is False then the power is average power [W], and 
        is multiplied by frep to calculate pulse energy"""
        try:
            self.twidth      = time_window          # Time window in ps
            self.center_wl  = center_wavelength    # reference wavelength (nm)            
            self.w0 = (2.0*np.pi*self.c)/self.center_wl # reference angular frequency (2pi THz)                
            
            self.setup_grids()

            if not power_is_epp:                
                power = power * self.frep
            
            # Read in retrieved FROG trace
            frog_data = np.genfromtxt(self.fileloc)                
            
            wavelengths = frog_data[:,0]# (nm)                    
            intensity   = frog_data[:,1]# (arb. units)
            phase       = frog_data[:,2]# (radians)

            if flip_phase:
                phase = -1 * phase
                
            freq_abs    = self.c/wavelengths
            freq_rel    = freq_abs - self.c / self.center_wl
            
            pulse_envelope = interp1d(freq_rel, intensity, kind='linear',
                                      bounds_error=False,fill_value=0)
            phase_envelope = interp1d(freq_rel, phase, kind='linear', 
                                      bounds_error=False,fill_value=0)
                                      
            gridded_intensity   = pulse_envelope(self.V/(2*np.pi))
            gridded_phase       = phase_envelope(self.V/(2*np.pi))
            
            # Calculate time domain complex electric field A
            self.A  =   IFFT_t(gridded_intensity*np.exp(1j*gridded_phase))
            # Calculate normalization factor  to achieve requested 
            # pulse energy
            e_scale = np.sqrt(power / self.calc_epp() )
            self.A      = self.A * e_scale

        except IOError:
            print 'File not found.'
            
    def setup_grids(self):
        ''' Helper function to set up time, frequency, and wavelength grids.
            Requires:   self.twidth, self.n, self.w0
            Generates:  T, dT, dt_seconds, W, V, wl, loWL, hiWL'''
        # Calculate center angualr frequency
        self.w0 = (2.0*np.pi*self.c)/self.center_wl # reference angular frequency (2pi THz)                
        # Create time axis            
        
            
        self.T = np.linspace(-self.twidth/2,self.twidth/2,self.n) # time grid
        self.dT      = self.T[1] - self.T[0]
        self.dt_seconds = self.dT * 1e-12
        
        # Create angular frequency and wavelength axes
        self.V = 2*np.pi*np.transpose(np.arange(-self.n/2,self.n/2))/(self.n*self.dT) # Frequency grid
        self.W = self.V+self.w0 # Absolute frequency grid
        
        # Derive wavelength grids
        self.wl = 2*np.pi*self.c / self.W

    def calc_epp(self):
        ''' Calculate and return energy per pulse via numerical integration
            of A^2 dt'''
        return self.dt_seconds * np.trapz(abs(self.A)**2)
    def chirp_pulse_W(self, GDD, TOD):
        ''' Add GDD and TOD to the pulse.'''
        self.A              = IFFT_t(np.exp(1j*(GDD/2.0)*self.V**2+\
                                1j*(TOD/6.0)*self.V**3)*FFT_t(self.A))
                                
    def chirp_pulse_T(self, chirp2, chirp3, T0):
        self.A = self.A*np.exp(-1j*(chirp2/2.0) * (self.T/T0)**2 + \
                               -1j*(chirp3/3.0) * (self.T/T0)**3)
                                
    def clone_pulse(self, pulse_instance):
        '''Copy all parameters of pulse_instance into this one'''
        p = pulse_instance
        self.twidth     =   p.twidth
        self.center_wl  =   p.center_wl
        self.n          =   p.n
        self.frep       =   p.frep
        self.setup_grids()
        self.A          =   np.copy(p.A)