# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 11:27:50 2019

@author: J.Wienand
"""

import math
from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
import pandas as pd
import qm
import numpy as np

c = 2.997*10**8
width1 = 32800000
width = 32800000
width2 = 28800000
qma = qm.QM()

lambdas = np.arange(450e-9,850e-9, 0.1e-9)

pot_factors2 = list()
for lambdal in lambdas:
    pfactor = qma.PotPrefactorAntiMagic(lambdal, width, mf=2, pol=1)
    pot_factors2.append(pfactor)
    
pot_factors3 = list()
for lambdal in lambdas:
    pfactor = qma.PotPrefactorAntiMagic(lambdal, width, mf=3, pol=1)
    pot_factors3.append(pfactor)
    
pot_factors22 = list()
for lambdal in lambdas:
    pfactor = qma.PotPrefactorAntiMagic(lambdal, width, mf=2, pol=-1)
    pot_factors22.append(pfactor)
    
pot_factors32 = list()
for lambdal in lambdas:
    pfactor = qma.PotPrefactorAntiMagic(lambdal, width, mf=3, pol=-1)
    pot_factors32.append(pfactor)
    
pot_factors0 = list()
for lambdal in lambdas:
    d2wl=852e-9
    d1wl=894e-9
    pfactor = qma.PotPrefactor(lambdal, width2, d2wl) +qma.PotPrefactor(lambdal, width1, d1wl)
    pot_factors0.append(pfactor) 
    
plt.plot(lambdas*1e9, pot_factors2,label="$\sigma_{+}, m_{F} = 2$")
plt.plot(lambdas*1e9, pot_factors3,label="$\sigma_{+}, m_{F} = 3$")
plt.plot(lambdas*1e9, pot_factors22,label="$\sigma_{-}, m_{F} = 2$")
plt.plot(lambdas*1e9, pot_factors32,label="$\sigma_{-}, m_{F} = 3$")
plt.legend(frameon=False)
#plt.plot(lambdas, pot_factors0)
plt.axhline(0, color="grey")
plt.axvline(870.8, color="grey") #sigma plus
plt.axvline(888.6, color="grey") #sigma minus
plt.xlabel("$\lambda$ (nm)")
plt.ylabel("$u_0$ ($m^2 s$)")
plt.ylim(-0.2e-34, 0.1e-34)

print(qma.PotPrefactorAntiMagic(870.8e-9, width, mf=3, pol=+1)/qma.PotPrefactorAntiMagic(888.6e-9, width, mf=3, pol=-1))
print(qma.ScattPrefactorAntiMagic(870.8e-9, width, mf=3, pol=+1)/qma.ScattPrefactorAntiMagic(888.6e-9, width, mf=3, pol=-1))