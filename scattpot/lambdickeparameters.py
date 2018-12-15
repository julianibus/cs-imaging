# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 15:50:24 2018

@author: J.Wienand
"""


import caesium
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
import pandas as pd
import math

#config
lambdals = [532]#, 755, 880, 1064, 1535]
output = "532.csv"
m = 132.905 * 1.660539* (10**(-27)) #Wikipedkia


c = 2.997*10**8
hbar = 1.054571*10**(-34)
cs = caesium.atom()

#decay rates
g23 = 4.05*10**6
g34 = 6.23*10**6
g35 = 11.4*10**6
g41 = 28.6*10**6
g26 = 0.13*10**6
g27 = 1.10*10**6
g75 = 0.78*10**6
g65 = 0.11*10**6
g51 = 32.8*10**6
g64 = 0.91*10**6
g21 = 1.84*10**6
########


levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
steadystates = []
#decaywleff = {"6s 2 S 1/2": 455.5*10**(-9)*1.84+894.3*10**(-9)*28.6+852.1*10**(-9)*32.8, "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"}
columns = ['laser wavelength','level', 'pot', 'scatt']
newdf = pd.DataFrame(columns=columns)

intensities = 10**(np.linspace(6,12, num=1000)) #1/m^^2

for lambdal in lambdals:
    
    
    for level in levels:
        
        pots = list()
        scatts = list()
        pot, scatt = cs.GetFactors(lambdal*10**-9, level, "transitions_complemented.csv")
        
        row = {'laser wavelength':lambdal,'level': level, 'pot': pot, 'scatt': scatt}
        newdf.loc[len(newdf)] = row
        
        omegas = list()
        lambdickes = list()
        for inten in intensities:
            omega = math.sqrt(pot*inten/m)*2 *math.pi /(lambdal * 10**(-9))
            omegas.append(omega)
        for inten in intensities:
            omega = math.sqrt(pot*inten/m)*2 *math.pi /(lambdal * 10**(-9))
            lambdicke = math.sqrt((hbar*(2*math.pi)/(1000*10**(-9))) /(hbar * omega))
            lambdickes.append(lambdicke)
        plt.loglog(intensities*10**(-1), lambdickes,"-o", marker=",", label=level) ##concvert so that xaxis is in mW/cm^2
        
    plt.title(str(lambdal) + "nm")
    plt.xlabel("Intensity (mW/cm2)")
    plt.ylabel("H.O. Frequency (1/s)")
    plt.legend(loc='upper left');
        
        
plt.show()        
newdf.to_csv(output, sep=";", index=False)
