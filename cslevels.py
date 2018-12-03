# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 16:31:02 2018

@author: J.Wienand
"""

import math
from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
import pandas as pd
import qm

lattice_laser_wl = 100 #nm

def isNaN(num):
    return num != num

file = 'transitions.csv'
df = pd.read_csv(file, delimiter=";")
df.columns = ['wavelength', 'unc', 'wavenumber', "Int","D","J","Jn",'decaywidth', 'lowerlevel', 'upperlevel',""]
print(df)
levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
for level in levels:
    selection = df.loc[df['lowerlevel'] == level]
    
    
    tot_pot_factor = 0
    c = 0
    for index, row in selection.iterrows():
        if not isNaN(row["decaywidth"]):
            c+=1
            width = float(row["decaywidth"])
            wl = float(row["wavelength"])
            pot_factor = qm.QM.PotPrefactor(lattice_laser_wl*10**-9, width, wl * 10**-10)
            tot_pot_factor += pot_factor
            #print(row["lowerlevel"],wl, width, pot_factor)
    
    print(row["lowerlevel"],tot_pot_factor, c)
            
            
# add scattering rate