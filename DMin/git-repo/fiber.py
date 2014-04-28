# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 14:10:46 2014

@author: dim1
"""
import numpy as np
from scipy.interpolate import interp1d
from pymongo import MongoClient
from DTabulationToBetas import DTabulationToBetas
from scipy.misc import factorial
from scipy import constants 

client = MongoClient('ds030827.mongolab.com',30827)
client.gnlse.authenticate('gnlseuser','1qaz2wsx')
fibers = client.gnlse['fibers']

#fiber:
#    "name"
#    "dispersion_format" ["D","GVD"]
#    "nonlinear_parameter"
#    if D:
#        "dispersion_x_units"
#        "dispersion_y_units"
#        "dispersion_data" [2,n]
#    if GVD:
#        "dispersion_gvd_units"
#        "dispersion_gvd_center_wavelength"
#        "dispersion_data"[1,n]
#    "is_gain"
#    if is_gain:
#        "gain_spectrum"
#        "gain_x_units"

class Fiber:
    def __init__(self):
        self.c_mks = constants.speed_of_light
        self.c = constants.speed_of_light * 1e9/1e12 # c in nm/ps
        self.is_simple_fiber = False
        
    def load_from_db(self, length, fibertype, poly_order = 2):
        self.fibertype = fibertype
        self.fiberspecs = fibers.find_one({"name" : fibertype})
        self.length = length
        self.betas = np.array([0])
        self.gamma = self.fiberspecs["gamma"]        
        self.poly_order = poly_order
        self.load_dispersion()
        
    def load_dispersion(self):
        if self.fiberspecs["dispersion_format"] == "D":
            if self.center_wavelength is None:
                print '''Error: you need to specify a center wavelength for
                    dispersion expansion'''
                return None
            self.dispersion_x_units = self.fiberspecs["dispersion_x_units"]
            self.dispersion_y_units = self.fiberspecs["dispersion_y_units"]
            self.x = self.fiberspecs["dispersion_x_data"]
            self.y = self.fiberspecs["dispersion_y_data"]
            return 1
        elif self.fiberspecs["dispersion_format"] == "GVD":
            self.dispersion_gvd_units = self.fiberspecs["dispersion_gvd_units"]
            self.center_wavelength = self.fiberspecs["dispersion_gvd_center_wavelength"]
            # If in km^-1 units, scale to m^-1
            if self.dispersion_gvd_units == 'ps^n/km':
                self.betas = np.array(self.fiberspecs["dispersion_data"]) / 1e3 
            return 1
        else:
            print "Error: no dispersion found."
            return None
    def get_betas(self,pulse):
        B = np.zeros(len(pulse.A))
        if self.fiberspecs["dispersion_format"] == "D":
            self.betas = DTabulationToBetas(pulse.center_wl,
                                            np.transpose(np.vstack((self.x,self.y))),
                                            self.poly_order,
                                            DDataIsFile = False)            
            for i in range(len(self.betas)):
                B = B + self.betas[i]/factorial(i+2)*pulse.V**(i+2)
            return B
        elif self.fiberspecs["dispersion_format"] == "GVD":
            # calculate beta[n]/n! * (w-w0)^n
            # w0 is the center of the Taylor expansion, and is defined by the 
            # fiber. the w's are from the optical spectrum
            fiber_omega0 =  2*np.pi*self.c / self.center_wavelength            
            betas = self.betas
            for i in range(len(betas)):
                betas[i] = betas[i]
                B = B + betas[i] / factorial(i + 2) * (pulse.W-fiber_omega0)**(i + 2)
            return B 
        else:
            return -1
    def get_gain(self,pulse,output_power = 0):
        if self.fiberspecs["is_gain"]:
            if self.is_simple_fiber:
                return self.gain
            else:
                self.gain_x_units = self.fiberspecs["gain_x_units"]
                x = np.array(self.fiberspecs["gain_x_data"])
                y = np.array(self.fiberspecs["gain_y_data"])
                f = interp1d(self.c_mks/x[::-1],y[::-1],kind ='cubic',
                             bounds_error=False,fill_value=0)
                gain_spec = f((pulse.w0+pulse.V) * 1e12/ (2*np.pi))
                self.native_gain =  sum(np.abs(pulse.A * np.exp(gain_spec*self.length/2))**2) / \
                               sum(np.abs(pulse.A)**2)
                print "native gain:", self.native_gain
                self.desired_gain = output_power / (sum( np.abs(pulse.A)**2 * pulse.dTSeconds)*pulse.frep)        
                print "desired gain:", self.desired_gain
                return gain_spec * np.log(self.desired_gain / self.native_gain)
        else:
#            print "no gain"
            return 1
    def Beta2_to_D(self, pulse): # in ps / nm / km
        return -2 * np.pi * self.c / pulse.wl**2 * self.Beta2(pulse) * 1000
    def Beta2(self, pulse):
        dw = pulse.V[1] - pulse.V[0]
        out = np.diff(self.get_betas(pulse), 2) / dw**2
        out = np.append(out[0], out)
        out = np.append(out, out[-1])
        return out
        
    def generate_fiber(self, length, center_wl, betas, gamma, gain = 0, 
                       gvd_units = 'ps^n/m', label = 'Simple Fiber'):                       
        self.length = length
        self.fiberspecs= {}
        self.fiberspecs['dispersion_format'] = 'GVD'
        self.fibertype = label
        if gain == 0:
            self.fiberspecs["is_gain"] = False            
        else:
            self.fiberspecs["is_gain"] = True            
        self.gain = gain

        self.center_wavelength = center_wl
        self.betas = np.copy(np.array(betas))
        self.gamma = gamma
        # If in km^-1 units, scale to m^-1
        if gvd_units == 'ps^n/km':
            self.betas = self.betas * 1.0e-3

def fiberlist(fibercollection = fibers):
    for each in fibercollection.find():
        print each["name"]
    pass