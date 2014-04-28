# -*- coding: latin1 -*-
#  Copyright (C) 2006 João Luís Silva <jsilva@fc.up.pt>
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
from numpy import fft
import string
import math
import time
from scipy.interpolate import interp1d

import units as units_module

class Record:
    pass

#-----------------------------------------------------------------------
# General support routines
#-----------------------------------------------------------------------



def isFloat(s):
    try:
        value = float(s)
        return True
    except:
        return False

def isInt(s):
    try:
        v_int = int(s)
        v_float = float(s)
        if v_int == v_float:
            return True
        else:
            return False
    except:
        return False

def FFT_t(A,ax=0):
    return fft.ifftshift(fft.ifft(fft.fftshift(A,axes=(ax,)),axis=ax),axes=(ax,))

def IFFT_t(A,ax=0):
    return fft.ifftshift(fft.fft(fft.fftshift(A,axes=(ax,)),axis=ax),axes=(ax,))
#    return fft.ifftshift(fft.fft(fft.fftshift(A)))

def FFT_x(A):
    return fft.ifftshift(fft.fft(fft.fftshift(A)))

def IFFT_x(A):
    return fft.ifftshift(fft.ifft(fft.fftshift(A)))

def Factorial(n):
    if n<=1:
        return 1
    i = 2
    fact = 1
    while i<=n:
        fact = fact * i
        i = i + 1

    return fact

######################################################################
# LoadData
######################################################################
def LoadData(filename,ncols=2):
    """Loads a data file, and automatically detects the separator."""
    file = open(filename,"r")

    line = file.readline()
    x = []
    y = []

    data = []
    for i in range(ncols):
        data.append([])

    delim_list = ["\t",","," ",";"]
    found_delim = False
    delim = "," #Default delimiter

    while line:
        line = string.strip(line)

        if not line:
            line = file.readline()
            continue

        if line[0] == '#':
            line = file.readline()
            continue

        if not found_delim:
            for e in delim_list:
                if len(line.split(e))==ncols:
                    delim = e
                    found_delim = True
                    break
            else:
                #Delimiter not found
                print "Invalid data file."
                return []
            
        l = string.split(line,delim)
        for i in range(ncols):
            if (not isFloat(l[i])):
                print "Invalid values found."
                return []

        for i in range(ncols):
            data[i].append(float(l[i]))

        line = file.readline()

    return np.array(data)

def LoadCSV(filename):
    """Loads a data file, and automatically detects the separator."""
    file = open(filename,"r")

    line = file.readline()
    data = []

    delim_list = ["\t",","," ",";"]
    found_delim = False
    delim = "," #Default delimiter

    while line:
        line = string.strip(line)

        if not line:
            line = file.readline()
            continue

        if line[0] == '#':
            line = file.readline()
            continue

        if not found_delim:
            for e in delim_list:
                if len(line.split(e))==2:
                    delim = e
                    found_delim = True
                    break
            else:
                #Delimiter not found
                Error(glb.main.window,"Invalid data file.")
                return []
            
        l = string.split(line,delim)
        if (not isFloat(l[0])) or (not isFloat(l[1])):
            Error(glb.main.window,"Invalid values found.")
            return []

        data.append([float(l[0]),float(l[1])])

        line = file.readline()

    return data

######################################################################
# SaveData
######################################################################
def SaveData(x,y,filename):
    """Saves data to a tab separated file"""
    file = open(filename,"w")
    n = len(x)
    for i in range(n):
        file.write("%g\t%g\n" % (x[i],y[i]))
    file.close()

def l2w(l):
    """Wavelength to angular frequency."""
    return 2.0*math.pi*units.c/l

def w2l(w):
    """Angular frequency to wavelength."""
    return 2.0*math.pi*units.c/w


#def ToProgUnits(value,from_unit):
#    return glb.Units.ToProgUnits(value,from_unit)
#
#def FromProgUnits(value,from_unit):
#    return glb.Units.FromProgUnits(value,from_unit)

def Autozoom(x,y,offset=0.0,level=1000.0,return_ind=False):
    """Return the x limits for the best zoom.
    If return_ind, return the indices instead."""
    m = max(y)
    n = len(x)
    limit = (m-offset)/level
    i = 0
    while i<n:
        if abs(y[i]-offset)>limit:
            break
        i += 1

    j = n-1
    while j>0:
        if abs(y[j]-offset)>limit:
            break
        j -= 1

    if not return_ind:
        if i>j:
            return (x[0],x[-1])
        else:
            return (x[i],x[j])
    else:
        if i>j:
            return (0,-1)
        else:
            return (i,j)

def Phase(A):
    n = len(A)
    Phase = zeros(n,dtype=A.dtype)
    last_phase = 0.0
    offset = 0.0

    i = 0
    while i<n:
        phase = math.atan2(A[i].imag,A[i].real)
        if abs(phase+2.0*math.pi - last_phase) < abs(last_phase-phase):
            offset += 2.0*math.pi
        elif abs(phase-2.0*math.pi - last_phase) < abs(last_phase-phase):
            offset -= 2.0*math.pi

        last_phase = phase
        Phase[i] = phase + offset

        i += 1

    return Phase

def InitialUppercase(s):
    if not s:
        return ""

    return s[0].upper()+s[1:]

def CalculateFWHM_OLD(I):
    """Determines the FWHM of the pulse in the array I. Returns -1 if multiple values are possible."""

    n = len(I)
    max2 = I.max() / 2.0

#Search where does the line crosses max/2
    count = 0
    pos = [0,0]
    intersect = [0,0]
    duration = -1.0
    for j in range(n-1):
        if (I[j]>=max2 and I[j+1]<=max2) or (I[j]<=max2 and I[j+1]>=max2):
            #Between j and j+1
            count = count + 1
            if count > 2:
                return -1
            else:
                pos[count-1] = j
                intersect[count-1] = j + (I[j]-max2)/(I[j]-I[j+1])
    if count == 2:
        return intersect[1]-intersect[0]

    return -1

def CalculateFWHM(I,choose="max"):
    """Determines the FWHM of the pulse in the array I. Returns -1 if multiple values are possible."""
    UP,DOWN = 0,1

    n = len(I)
    max2 = I.max() / 2.0

#Search where does the line crosses max/2
    count = 0
    intersect = []
    dir = []
    for j in range(n-1):
        if (I[j]>=max2 and I[j+1]<=max2) or (I[j]<=max2 and I[j+1]>=max2):
            #Between j and j+1
            count += 1
            intersect.append(j + (I[j]-max2)/(I[j]-I[j+1]))
            if (I[j]>=max2 and I[j+1]<=max2):
                dir.append(DOWN)
            else:
                dir.append(UP)

    duration = []
    n = len(intersect)
#    print "intersect",intersect,"dir",dir
    for i in range(n-1):
        if dir[i]==UP and dir[i+1]==DOWN:
            duration.append(intersect[i+1]-intersect[i])

    if not duration:
        return -1

    if choose.upper()=="MAX":
        return max(duration)
    else:
        return min(duration)


def FreqToWavelength(A,f,freq_min,freq_max):
#    print "FreqToWavelength: A.shape",A.shape,"f.shape",f.shape

    #Remove everything below x THz
    f_ind = FindInd(f,freq_min)
#    print "f_ind",f_ind,"f[f_ind]",f[f_ind]
    if f_ind != -1:
        A = A[f_ind:]
        f = f[f_ind:]

    #Remove everything above x THz
    f_ind = FindInd(f,freq_max)
#    print "f_ind",f_ind,"f[f_ind]",f[f_ind]
    if f_ind != -1:
        A = A[:f_ind]
        f = f[:f_ind]

#    print "FreqToWavelength: A.shape",A.shape,"f.shape",f.shape
    l_min = units.c/f[-1]
    l_max = units.c/f[0]
#    print "l_min",l_min,"l_max",l_max
    l = np.linspace(l_min,l_max,len(f))

    nl = len(l)
    A_l = np.zeros(nl,dtype=A.dtype)

    l_f = units.c/f

    line = 2.0*math.pi*units.c*A[:]/l_f**2
    A_l[:] = interp1d(l_f[::-1],line[::-1])(l)

    A = A_l.copy()
    return A,l

def CalculateOmega(n,T):
    omega = np.zeros(n,dtype=np.float64)

    k = -n//2
    for i in range(n):
        omega[i] = (k/T)*2.0*math.pi
        k += 1

    return omega

