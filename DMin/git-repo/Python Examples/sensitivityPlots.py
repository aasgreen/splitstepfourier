# -*- coding: utf-8 -*-
"""
Created on Wed Jan 08 12:34:00 2014

@author: ycasg
"""

import tables
import matplotlib.pyplot as plt
import numpy

plt.rcParams['axes.linewidth']=1
plt.rcParams['axes.labelsize']=10

#plt.rcParams['axes.titlesize']='

plt.rcParams['xtick.major.size']=0
plt.rcParams['xtick.major.width']=0
plt.rcParams['ytick.major.size']=0
plt.rcParams['ytick.major.width']=0
plt.rcParams['xtick.labelsize']=10
plt.rcParams['ytick.labelsize']=10
plt.rcParams['legend.fontsize']=10
plt.rcParams['ytick.major.pad']=10
#plt.rcParams['font.family'] = 'Myriad Pro'



molecules = ['h2o','n2o','co2_c12','co2_c13']

filebase = 'C:\\Users\\ycasg\\Downloads\\transmission_100m\\'

h5file = tables.open_file(filebase+'spectraDB_100m.h5', mode='r')

BG = h5file.root.h2o_o16
FG = [h5file.root.n2o_o16,
      h5file.root.n2o_o18,
      h5file.root.n2o_n15,
      h5file.root.co2_c12,
      h5file.root.co2_c13,
      h5file.root.h2o_D,
      h5file.root.h2o_o18,
      h5file.root.ch4_c12,
      h5file.root.ch4_c13,
      h5file.root.ch4_D]
FGNames = ['$N_2^{16}O$',
           '$N_2^{18}O$',
           '$^{14}N^{15}N^{16}O$',
           '$^{12}CO_2$',
           '$^{13}CO_2$',
           '$HDO$',
           '$H_2^{18}O$',
           '$CH_4$',
           '$^{13}CH_4$',
           '$CH_3D$']
# Abundances
# N14: 0.99636             N15: 0.00364
# O16: 0.99757             O18: 2e-3
# C12: 0.9893              C13: 0.0107
# H:  0.999885             D: 0.000115

# Air-drying factor
bgScaleFactor = 1

weightFactors = [1, 2e-3, 3.64e-3, # N2O O16, O18, N15
                 0.9893, 0.0107, # CO2 C12 C13
                 bgScaleFactor*0.000115 * 2, bgScaleFactor*4e-3, # H2O H D
                 0.9893, 0.0107,0.000115 * 4] # CH4 C12 H, C13 H, C12 D



colors = ['r','g','b',
          'r','c',
          'm','k',
          'r','r','orange']
pctChange = 1
#molChange = range(len(FGNames))
molChange = [1,2,
             4,
             5,6,
             8,9]

wlAxis0 = BG.wavelength.read()
wlStart = 2.5
wlStop  = 5
wlAxis = wlAxis0[(wlAxis0 > wlStart) * (wlAxis0 < wlStop)]

# Initialize plots arrays
runningTotal = numpy.ones(len(wlAxis)) * 1e-50
prevPlot     = numpy.ones(len(wlAxis)) * 1e-50

fig, ax1 = plt.subplots()
fig.set_size_inches(6.5,4)

ax2 = ax1.twinx()
p = []
labels = []

p.append ( plt.Rectangle((0, 0), 1, 1, fc='k', alpha = 0.3))
labels.append('Total Transmission')  
ax2.set_ylabel('Transmission')
ax2.set_ylim(0,1)


for j in molChange:
    print "Computing for %s"%FGNames[j]   
    runningTotalA = numpy.interp(wlAxis,
                             BG.wavelength.read(),
                             BG.t.read())
    runningTotalB = runningTotalA[:]    

    for i in range(len(FGNames)):
        additionA = numpy.exp(numpy.log(runningTotalA)+\
                    weightFactors[i]*numpy.log(numpy.interp(wlAxis,
                                 FG[i].wavelength.read(),
                                 FG[i].t.read())))
        runningTotalA = additionA[:]
        if j == i:
            additionB = numpy.exp(numpy.log(runningTotalB)+\
                    weightFactors[j]*(1.0-0.01*pctChange)*
                                numpy.log(numpy.interp(wlAxis,
                                 FG[i].wavelength.read(),
                                 FG[i].t.read())))
        else:
            additionB = numpy.exp(numpy.log(runningTotalB)+\
                    weightFactors[i]*
                                numpy.log(numpy.interp(wlAxis,
                                 FG[i].wavelength.read(),
                                 FG[i].t.read())))
    
        runningTotalB = additionB[:]
                     
    tbp  = 10.0*numpy.log10(1e-100+numpy.abs(runningTotalA - runningTotalB))
    print "Plot range %f:%f"%(min(tbp), max(tbp))
    ax1.plot(wlAxis, tbp,
                     color = colors[j], lw=0.5, alpha=0.9)
    prevPlot += numpy.abs(runningTotalA - runningTotalB)
    p.append ( plt.Rectangle((0, 0), 1, 1, fc=colors[j], alpha = 0.9))
    labels.append(FGNames[j])

#ax2.plot(wlAxis, 10*numpy.log10(runningTotalA),
#        'k-', lw = 2.0, alpha = 0.3)

plt.legend(p, labels, loc=2)
print max(prevPlot)                 
                 
plt.xlim(wlStart, wlStop) 
ax1.set_ylim(-100, -30)
ax1.set_ylabel('Change in T per 1 % change\n in concentration (dB %/%)')
ax1.set_xlabel('Wavelength ($\mu$m)')
#ax = plt.gca()
#ax.grid(True)
plt.draw()
#plt.savefig('C:\Users\ycasg\Documents\NRC - Spectral Stuff\\sens_wideband.png',
#            format='png',dpi=600)
#plt.xlim(4.34,4.5)            
#plt.savefig('C:\Users\ycasg\Documents\NRC - Spectral Stuff\\sens_narrowband.png',
#            format='png',dpi=600)
plt.show()
h5file.close()         