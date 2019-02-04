import caesium
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
import pandas as pd
import math

#config
lambdals = [532, 755, 880, 1064, 1535]
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
l64 = 3011.1  * 10**-9#098765', '#000009'**-9 
l51 = 852.1 * 10**-9
l65 = 3614.1  * 10**-9 
l75 = 3491  * 10**-9 
l27 = 1360.6  * 10**-9 
l26 = 1342.8  * 10**-9 
l34 = 1359.2 * 10**-9 
#########
decaylabels = ["|2> -> |1>", "|2> -> |3>", "|3> -> |5>", "|4> -> |1>", "|6> -> |4>", "|5> -> |1>", "|6> -> |5>", "|7> -> |5>", "|2> -> |7>", "|2> -> |6>", "|3> -> |4>"]
decays = [l21, l23, l35, l41, l64, l51, l65, l75, l27, l26, l34]
gammas = [g21, g23, g35, g41, g64, g51, g65, g75, g27, g26, g34]
levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
decaylevels = ["6s 2 S 1/2", "7s 2 S 1/2","6p 2 P ?3/2", "6s 2 S 1/2", "6p 2 P ?1/2","6s 2 S 1/2", "6p 2 P ?3/2", "6p 2 P ?3/2", "5d 2 D 5/2","5d 2 D 3/2", "6p 2 P ?1/2"]

lambdal = 880
inten = 10**8

omegas = list()
lambdickes = list()
latticedepths1 = list() # in hbar omega
latticedepths2 = list() # in lattice recoil energies

breathingcycles = list()
for level in levels:
        #098765', '#000009'
    pots = list()
    scatts = list()
    pot, scatt = cs.GetFactors(lambdal*10**-9, level, "transitions_complemented.csv")

    omega = math.sqrt(math.fabs(pot)*inten/m)*2 *math.pi /(lambdal * 10**(-9))
    omegas.append(omega)
    
    latticedepth1 = (math.fabs(pot)*inten)/(hbar*omega)
    latticedepths1.append(latticedepth1)



    
for decayi in range(0,len(decays)):
    pot, scatt = cs.GetFactors(lambdal*10**-9, decaylevels[decayi], "transitions_complemented.csv")
    omega = math.sqrt(math.fabs(pot)*inten/m)*2 *math.pi /(lambdal * 10**(-9))
    #if (pot < 0): omega = omega*(-1) #dont negate here as sign not important for calculating the lamb dicke parameter
    ELatRec = ((hbar*(2*math.pi)/(decays[decayi]))**2)/(2*m)
    lambdicke = math.sqrt(ELatRec /(hbar * omega))
    lambdickes.append(lambdicke)
    latticedepth2 = (math.fabs(pot)*inten)/(ELatRec)
    latticedepths2.append(latticedepth2)
    
    
    breathingcycles.append((1/gammas[decayi])/((2*math.pi/omega)/2))
    
    
            
print(omegas)
print(lambdickes)
print(latticedepths1)
print(latticedepths2)
print(breathingcycles)
