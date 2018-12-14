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

levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
steadystates = []
#decaywleff = {"6s 2 S 1/2": 455.5*10**(-9)*1.84+894.3*10**(-9)*28.6+852.1*10**(-9)*32.8, "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"}
columns = ['laser wavelength','level', 'pot', 'scatt']
newdf = pd.DataFrame(columns=columns)

intensities = 10**(np.linspace(6,10, num=1000)) #1/m^^2
print (intensities)
exit(0)
for lambdal in lambdals:
    
    
    for level in levels:
        
        pots = list()
        scatts = list()
        pot, scatt = caesium.atom.GetFactors(lambdal*10**-9, level, "transitions_complemented.csv")
        
        row = {'laser wavelength':lambdal,'level': level, 'pot': pot, 'scatt': scatt}
        newdf.loc[len(newdf)] = row
        
        omegas = list()
        lambdickes = list()
        for inten in intensities:
            omega = math.sqrt(pot*inten/m)*2 *math.pi /(lambdal * 10**(-9))
            omegas.append(omega)
        #for inten in intensities:
        #    lambdicke = math.sqrt((hbar*(2*math.pi)/(decaywleff["level"])) /(hbar * omega))
        #    lambdickes.append(lambdicke)
        plt.loglog(intensities*10**(-1), omegas,"-o", marker=",", label=level) ##concvert so that xaxis is in mW/cm^2
        
    plt.title(str(lambdal) + "nm")
    plt.xlabel("Intensity (mW/cm²)")
    plt.ylabel("H.O. Frequency (1/s)")
    plt.legend(loc='upper left');
        
        
        
newdf.to_csv(output, sep=";", index=False)