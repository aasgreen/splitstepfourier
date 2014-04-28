# -*- coding: utf-8 -*-
"""
Created on Wed Jan 08 12:08:45 2014

@author: ycasg
"""

import numpy as np
import tables

filebase = 'C:\\Users\\ycasg\\Downloads\\transmission_100m\\'

molecule  = 'co2_o18'
nFiles = 1

filename0 =  molecule+'_reallynearIR.txt'
filename1 =  molecule+'_nearIR.txt'
filename2 =  molecule+'_midIR.txt'

h5file = tables.open_file(filebase+'spectraDB_100m.h5', mode='a',
                          title = 'Absorption Spectra')

try:
    gcols = h5file.create_group(h5file.root, molecule, molecule)
except tables.exceptions.NodeError:
    gcols = h5file.root.columns
    print "already had col"


if nFiles == 3:
    print('Reading very near IR')
    data0 = np.genfromtxt(filebase+filename0, skip_header=18)
    print('Reading near IR')
    data1 = np.genfromtxt(filebase+filename1, skip_header=18)
    print('Reading mid IR')
    data2 = np.genfromtxt(filebase+filename2, skip_header=18)
    
    print('Generating hdf5 tables')
    h5file.create_array(gcols, 'wavelength', 
                        np.hstack((data0[:,0],
                                   data1[:,0],
                                   data2[:,0])),
                        "Wavelength Axis")
    h5file.create_array(gcols, 't', 
                        np.hstack((data0[:,1],
                                   data1[:,1],
                                   data2[:,1])),
                        "Transmittance Axis")
else:
    if nFiles == 2:
        print('Reading near IR')
        data1 = np.genfromtxt(filebase+filename0, skip_header=18)
        print('Reading mid IR')
        data2 = np.genfromtxt(filebase+filename1, skip_header=18)
        
        print('Generating hdf5 tables')
        h5file.create_array(gcols, 'wavelength', 
                            np.hstack((data1[:,0],
                                       data2[:,0])),
                            "Wavelength Axis")
        h5file.create_array(gcols, 't', 
                            np.hstack((data1[:,1],
                                       data2[:,1])),
                            "Transmittance Axis")    
    if nFiles == 1:
        print('Reading near IR')
        data1 = np.genfromtxt(filebase+filename1, skip_header=18)
        
        print('Generating hdf5 tables')
        h5file.create_array(gcols, 'wavelength',data1[:,0],
                            "Wavelength Axis")
        h5file.create_array(gcols, 't',data1[:,1],
                            "Transmittance Axis")    
print h5file

h5file.close()

