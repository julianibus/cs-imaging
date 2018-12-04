# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 09:41:03 2018

@author: J.Wienand
"""
import caesium
import numpy as np
import matplotlib.pyplot as plt

levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
plt.figure(figsize=(20,10))
lambdals = np.linspace(400.0,1800,num=300)
print(lambdals)

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
    print(pots)
plt.subplot(211)
plt.title("Optical Potential Prefactor")
plt.xlabel("Wavelength (nm)")
plt.legend(loc='upper left');
plt.ylim([0,5])
plt.subplot(212)
plt.title("Scattering Rate Prefactor")
plt.xlabel("Wavelength (nm)")
plt.ylim([-0.5,0.5])
plt.legend(loc='upper left');
plt.show()