# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 09:41:03 2018

@author: J.Wienand
"""
import caesium
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
plt.figure(figsize=(16,10))
lambdals = np.linspace(324.6,1800,num=100)
#print(lambdals)

for level in levels:
    pots = list()
    scatts = list()
    for lambdal in lambdals:
        pot, scatt = caesium.atom.GetFactors(lambdal*10**-9, level, "transitions.csv")
        pots.append(pot)
        scatts.append(scatt)
    plt.subplot(211)
    plt.plot(lambdals, scatts,"-o", marker=",", label=level)
    plt.subplot(212)
    plt.plot(lambdals, pots,"-o", marker=",", label=level)
    #print(pots)
plt.subplot(211)



plt.title("Scattering Rate Prefactor")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Value (m²/Ws)")
plt.legend(loc='upper left');
plt.ylim([0,5])
plt.xlim([200,1800])
plt.axvline(x=532, color="grey", dashes=[6, 2])
plt.axvline(x=755, color="grey", dashes=[6, 2])
plt.axvline(x=880, color="grey", dashes=[6, 2])
plt.axvline(x=1064, color="grey", dashes=[6, 2])
plt.axvline(x=1535, color="grey", dashes=[6, 2])
plt.axhline(y=0, xmin=0, xmax=1800, linewidth=2, color = 'black',dashes=[6, 2])
plt.axvspan(200, 324.6, alpha=0.5, color='grey')

ax = plt.gca()
minorLocator = MultipleLocator(10)
majorLocator = MultipleLocator(100)
majorFormatter = FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)



plt.subplot(212)
plt.title("Optical Potential Prefactor")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Value (ms²)")
plt.ylim([-0.5,0.5])
plt.xlim([200,1800])
plt.legend(loc='upper left');
plt.axvline(x=532, color="grey", dashes=[6, 2])
plt.axvline(x=755, color="grey", dashes=[6, 2])
plt.axvline(x=880, color="grey", dashes=[6, 2])
plt.axvline(x=1064, color="grey", dashes=[6, 2])
plt.axvline(x=1535, color="grey", dashes=[6, 2])
plt.axvspan(200, 324.6, alpha=0.5, color='grey')
plt.tight_layout()
ax = plt.gca()
#plt.axes().yaxis.set_tick_params(which='minor', right = 'on')
ax.get_yaxis().set_tick_params(which='both', direction='in')
ax.get_xaxis().set_tick_params(which='both', direction='in')
ax.get_xaxis().set_tick_params(which='both', right='on',left='on',bottom='on',top='on')

minorLocator = MultipleLocator(10)
majorLocator = MultipleLocator(100)
majorFormatter = FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)

plt.show()