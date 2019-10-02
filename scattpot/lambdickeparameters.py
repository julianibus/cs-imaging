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
from matplotlib.transforms import Transform

#config
lambdals = [535, 767, 872.7, 1064]
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

levellabels = ['$|1\\rangle \; \; 6s_{1/2}$', "$|2\\rangle \; \; 7p_{3/2}$", "$|3\\rangle \; \; 7s_{1/2}$", "$|4\\rangle \; \; 6p_{1/2}$", "$|5\\rangle \; \; 6p_{3/2}$", "$|6\\rangle \; \; 5d_{3/2}$", "$|7\\rangle \; \; 5d_{5/2}$"]
levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
columns = ['laser wavelength','level', 'pot', 'scatt']
newdf = pd.DataFrame(columns=columns)

intensities = 10**(np.linspace(5,10, num=30)) #mW/cm^2 ##LATTICE INTENSITY
plt.figure(figsize=(7,5))
plt.gca().set_prop_cycle(color=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan'])
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
            inten = inten*10
            omega = math.sqrt(2*math.fabs(pot)*inten/m)*2 *math.pi /(lambdal * 10**(-9))
            omegas.append(omega)
            
        if (pot < 0):
            plt.loglog(intensities, omegas,"-o", marker=",", label=levellabels[c]) ##concert so that xaxis is in mW/cm^2
        else: 
            plt.loglog(intensities, omegas,"-o", marker=",", label=levellabels[c], linestyle = "--") 
            
    def topotdepth(x):
        x**2 * (lambdal * 10**(-9))**4 * m**2 /((2*math.pi)**4 * hbar**2)
    def toomega(x):
        math.sqrt(x* ((2*math.pi)**4 * hbar**2)/((lambdal * 10**(-9))**4 * m**2))
    
    #ax = plt.gca()
    #secax = ax.secondary_xaxis('top', functions=(topotdepth, toomega))
    #secax.set_xlabel('angle [rad]')
    plt.gca().text(0.04, 0.85, str(lambdal) + " nm", fontsize=10, transform=plt.gca().transAxes)
    plt.xlabel("$I$ ($W/m^{2}$)")
    plt.xlim(10**5,10**9)
    plt.ylabel("$\omega$ (1/s)")
    plt.legend(loc=4, prop={'size': 6})
    plt.grid(b=True, which='both', color='0.85', linestyle='-')
   # plt.savefig("HO" + str(lambdal) + ".png")
plt.tight_layout()

plt.savefig("omegas.png")
plt.show()
plt.figure(figsize=(7,5))
plt.gca().set_prop_cycle(color=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan'])
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
        pdepths= list()

        for inten in intensities:
            inten = inten*10
            omega = math.sqrt(2*math.fabs(pot)*inten/m)*2 *math.pi /(lambdal * 10**(-9))
            #omegas.append(omega)
            pdepths.append(inten*math.fabs(pot)/(hbar**2/(2*m) *(2 * math.pi/(lambdal * 10**(-9)))**2))
        if (pot < 0):
            plt.loglog(intensities, pdepths,"-o", marker=",", label=levellabels[c]) ##concert so that xaxis is in mW/cm^2
        else: 
            plt.loglog(intensities, pdepths,"-o", marker=",", label=levellabels[c], linestyle = "--") 
            
        
    plt.gca().text(0.04, 0.85, str(lambdal) + " nm", fontsize=10, transform=plt.gca().transAxes)
    plt.xlabel("$I$ ($W/m^{2}$)")
    plt.xlim(10**5,10**9)
    plt.ylabel("$\omega$ (1/s)")
    plt.legend(loc=4, prop={'size': 6})
    plt.grid(b=True, which='both', color='0.85', linestyle='-')
   # plt.savefig("HO" + str(lambdal) + ".png")
plt.tight_layout()

plt.savefig("potdepth.png")
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
decaylabels = ['$|2\\rangle \\rightarrow |1\\rangle $', '$|2\\rangle \\rightarrow |3\\rangle $', '$|3\\rangle \\rightarrow |5\\rangle $', '$|4\\rangle \\rightarrow |1\\rangle $', '$|6\\rangle \\rightarrow |4\\rangle $', '$|5\\rangle \\rightarrow |1\\rangle $', '$|6\\rangle \\rightarrow |5\\rangle $', '$|7\\rangle \\rightarrow |5\\rangle $', '$|2\\rangle \\rightarrow |7\\rangle $', '$|2\\rangle \\rightarrow |6\\rangle $', '$|3\\rangle \\rightarrow |4\\rangle $']
decays = [l21, l23, l35, l41, l64, l51, l65, l75, l27, l26, l34]
levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
decaylevels = ["6s 2 S 1/2", "7s 2 S 1/2","6p 2 P ?3/2", "6s 2 S 1/2", "6p 2 P ?1/2","6s 2 S 1/2", "6p 2 P ?3/2", "6p 2 P ?3/2", "5d 2 D 5/2","5d 2 D 3/2", "6p 2 P ?1/2"]

plt.figure(figsize=(7,5))
c3 = 0
for lambdal in lambdals:
    c3+=1
    plt.subplot(2,2,c3)
    plt.gca().set_prop_cycle(color=["blue", "orange", "green", 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', "black"])

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
     
    plt.gca().text(0.72, 0.85, str(lambdal) + " nm", fontsize=10, transform=plt.gca().transAxes)
    plt.xlabel("$I$ ($mW/cm^{2}$)")
    plt.xlim(10**5,10**9)
    plt.ylim(2*10**-2,10**0)
    plt.ylabel("$\eta$")
    if c3 == 1:    
        plt.ylim(10**-2,10**0)
        plt.legend(prop={'size': 5},loc=8, ncol=4)
    plt.grid(b=True, which='both', color='0.85', linestyle='-')
plt.tight_layout()

plt.savefig("etas.png")
plt.show()        
newdf.to_csv(output, sep=";", index=False)



