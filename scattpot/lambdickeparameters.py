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
lambdals = [767, 880, 1064, 1535]
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
l21 = 455.5 * 10**-9 
l23 = 2931.8 * 10**-9 
l35 = 1469.9  * 10**-9 
l41 = 894.3  * 10**-9 
l64 = 3011.1  * 10#098765', '#000009'**-9 
l51 = 852.1 * 10**-9
l65 = 3614.1  * 10**-9 
l75 = 3491  * 10**-9 
l27 = 1360.6  * 10**-9 
l26 = 1342.8  * 10**-9 
l34 = 1359.2 * 10**-9 
#########


#ss = [0.25420717,0.25407444, 0.0583665 , 0.01374444, 0.02891531, 0.03238204, 0.35831011] #STEADY STATES from boch equations @saturation
#in ss: only account for "incoming" transitions= excite up or decay down

#decaywleff = {"6s 2 S 1/2": (l41*g41*ss[0]+l21*g21*ss[1]+l51*g51*ss[4])/(g41*ss[0]+g21*ss[1]+g51*ss[4]), "7p 2 P ?3/2": l21, "7s 2 S 1/2": l23, "6p 2 P ?1/2": (l34*g34*ss[2]+l64*g64*ss[5])/(g34*ss[2]+g64*ss[5]), "6p 2 P ?3/2": (l75*g75*ss[6]+l65*g65*ss[5])/(g75*ss[6]+ g65*ss[5]), "5d 2 D 3/2": l26, "5d 2 D 5/2": l27}

levellabels = ["|1> $6s_{1/2}$", "|2> $7p_{3/2}$", "|3> $7s_{1/2}$", "|4> $6p_{1/2}$", "|5> $6p_{3/2}$", "|6> $5d_{3/2}$", "|7> $5d_{5/2}$"]
levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
columns = ['laser wavelength','level', 'pot', 'scatt']
newdf = pd.DataFrame(columns=columns)

intensities = 10**(np.linspace(6,10, num=30)) #1/m^^2 ##LATTICE INTENSITY
plt.figure(figsize=(10,7))
ax.set_prop_cycle(color=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan'])
c3 = 0
for lambdal in lambdals:
    c3+=1
    plt.subplot(2,2,c3)
    print(lambdal)
    c = -1
    for level in levels:
        c = c + 1
        #098765', '#000009'
        pots = list()
        scatts = list()
        pot, scatt = cs.GetFactors(lambdal*10**-9, level, "transitions_complemented.csv")
        row = {'laser wavelength':lambdal,'level': level, 'pot': pot, 'scatt': scatt}
        newdf.loc[len(newdf)] = row
        
        omegas = list()
        lambdickes = list()
        for inten in intensities:
            omega = math.sqrt(2*math.fabs(pot)*inten/m)*2 *math.pi /(lambdal * 10**(-9))
            omegas.append(omega)
            
        if (pot > 0):
            plt.loglog(intensities, omegas,"-o", marker=",", label=levellabels[c]) ##concert so that xaxis is in mW/cm^2
        else: 
            plt.loglog(intensities, omegas,"-o", marker=",", label=levellabels[c], linestyle = "--") 
            
        
    plt.title(str(lambdal) + "nm")
    plt.xlabel("Lattice Beam Intensity (W/m^2)")
    plt.xlim(10**6,10**10)
    plt.ylabel("H.O. Frequency (1/s)")
    plt.legend(loc=4, prop={'size': 10})
    plt.grid(b=True, which='both', color='0.65', linestyle='-')
   # plt.savefig("HO" + str(lambdal) + ".png")

plt.show()

l21 = 455.5 * 10**-9 
l23 = 2931.8 * 10**-9 
l35 = 1469.9  * 10**-9 
l41 = 894.3  * 10**-9 
l64 = 3011.1  * 10**-9 
l51 = 852.1 * 10**-9
l65 = 3614.1  * 10**-9 
l75 = 3491  * 10**-9 
l27 = 1360.6  * 10**-9 
l26 = 1342.8  * 10**-9 
l34 = 1359.2 * 10**-9 
decaylabels = ["|2> -> |1>", "|2> -> |3>", "|3> -> |5>", "|4> -> |1>", "|6> -> |4>", "|5> -> |1>", "|6> -> |5>", "|7> -> |5>", "|2> -> |7>", "|2> -> |6>", "|3> -> |4>"]
decays = [l21, l23, l35, l41, l64, l51, l65, l75, l27, l26, l34]
levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
decaylevels = ["6s 2 S 1/2", "7s 2 S 1/2","6p 2 P ?3/2", "6s 2 S 1/2", "6p 2 P ?1/2","6s 2 S 1/2", "6p 2 P ?3/2", "6p 2 P ?3/2", "5d 2 D 5/2","5d 2 D 3/2", "6p 2 P ?1/2"]

for lambdal in lambdals:
    plt.clf()
    fig, ax = plt.subplots()
    ax.set_prop_cycle(color=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan'])
    print(lambdal)
    
    for decayi in range(0,len(decays)):

        lambdickes = list()
        for inten in intensities:
            pot, scatt = cs.GetFactors(lambdal*10**-9, decaylevels[decayi], "transitions_complemented.csv")
            omega = math.sqrt(math.fabs(pot)*inten/m)*2 *math.pi /(lambdal * 10**(-9))
            #if (pot < 0): omega = omega*(-1) #dont negate here as sign not important for calculating the lamb dicke parameter
            ELatRec = ((hbar*(2*math.pi)/(decays[decayi]))**2)/(2*m)
            lambdicke = math.sqrt(ELatRec /(hbar * omega))
            lambdickes.append(lambdicke)
            
        plt.loglog(intensities, lambdickes,"-o", marker=",", label=decaylabels[decayi]) ##concert so that xaxis is in mW/cm^2
     
    plt.title(str(lambdal) + "nm")
    plt.xlabel("Lattice Beam Intensity (W/m^2)")
    plt.ylabel("Lamb-Dicke Parameter")
    plt.xlim(10**6,10**10)
    plt.legend(loc=3, prop={'size': 10})
    plt.grid(b=True, which='both', color='0.65', linestyle='-') 
    plt.savefig("LD" + str(lambdal) + ".png")
plt.show()        
newdf.to_csv(output, sep=";", index=False)



