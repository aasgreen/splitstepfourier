# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 14:10:46 2014

@author: dim1
"""
import numpy as np
import os

from pymongo import MongoClient

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

gammas = {"EDF07PM" : 5e-3, "PANDA1550" : 1.1e-3, "HNLF2" : 25e-3, "HNLF3" : 25e-3, "NL1050_ZERO" : 25e-3, "NL-1050-NEG-1" : 25e-3, "NL-PM-750" : 25e-3, "HC1060-02" : 1e-3, "HI1060" : 1e-3, "SC-3.7-975" : 18e-3, "dudley" : 0.11}

for each in os.listdir('.'):
    if each.find('.csv') != -1:
        if each == 'er7pm_gain_formatted.csv':
            pass
        else:
            if each == 'EDF07PM.csv':
                fibertype = each[0:each.find('.csv')]
                fiberspecs = np.genfromtxt(each,delimiter=',')
                gainspecs = np.genfromtxt('er7pm_gain_formatted.csv', delimiter=',')
                fiberstring = {"name" : fibertype, "dispersion_format" : "D", "dispersion_x_units" : "nm", "dispersion_y_units" : "ps/nm/km", "dispersion_x_data" : list(fiberspecs[:,0]), "dispersion_y_data" : list(fiberspecs[:,1]), "gain_x_data" : list(gainspecs[:,0]), "gain_y_data" : list(gainspecs[:,1]), "gain_x_units" : "nm", "is_gain" : True, "gamma" : gammas[fibertype]}
            if each != 'EDF07PM.csv' and each!= 'dudley.csv':
                fibertype = each[0:each.find('.csv')]
                fiberspecs = np.genfromtxt(each,delimiter=',')
                fiberstring = {"name" : fibertype, "dispersion_format" : "D", "dispersion_x_units" : "nm", "dispersion_y_units" : "ps/nm/km", "dispersion_x_data" : list(fiberspecs[:,0]), "dispersion_y_data" : list(fiberspecs[:,1]), "is_gain" : False, "gamma" : gammas[fibertype]}
            if each == 'dudley.csv':
                fibertype = 'dudley'
                fiberspecs = np.genfromtxt(each,delimiter=',')
                fiberstring = {"name" : fibertype, "dispersion_format" : "GVD", "dispersion_gvd_units" : "ps^n/km",  "dispersion_gvd_center_wavelength" : 835.0, "dispersion_data" : list(fiberspecs), "is_gain" : False, "gamma" : gammas[fibertype]}
            try:
                fibers.remove({"name" : fiberstring["name"]})
            except KeyError:
                pass
            fibers.insert(fiberstring)


for each in fibers.find():
    print each["name"]