# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 15:38:25 2014

@author: dim1
"""
import numpy as np
import matplotlib.pyplot as plt
import tables

#try:
#    h5file.close()
#except NameError:
#    pass
#
#data = np.genfromtxt('test.csv',delimiter=',')
#
#h5file = tables.open_file('test.h5',mode='w',title='Base')
#
#gcols = h5file.create_group(h5file.root, 'nextdir', 'label')
#
#h5file.create_array(gcols,'arrayset',data,'this is test data')
#
#h5file.close()
#
#try:
#    h5fileload.close()
#except NameError:
#    pass

h5fileload = tables.open_file('pyout.h5', mode='r')
loaddata = h5fileload.root.nextdir.arrayset.read()
h5fileload.close()
