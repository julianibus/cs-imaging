# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 07:24:45 2018

@author: J.Wienand
"""

import matplotlib.pyplot as plt  
import numpy as np               
from IPython.core.display import display, HTML 
from arc import *
import math
from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
import pandas as pd 
import os

output = "transitions_complemented.csv"


def isNaN(num):
    return num != num

def getQuantumNumbers(ident):
    try:
    
        spl = ident.split(" ")
        nl = spl[0]
        n = nl[0:-1]
        letter = nl[-1]
        assignment = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6}
        l = assignment[letter]
        jfrac = spl[3]
       
        jfracspl = jfrac.split("/")
        j = float(jfracspl[0][-1])/float(jfracspl[1][0])
       
        return int(n), int(l) ,float(j)
    except:
        return None
    
def DecayWidth(wl,j,jt,D):
    return (8*(3.14159)**2)/(3*8.85419*10**(-12)*1.0545718*10**(-34)*(wl*10**(-10))**3)*(2*j + 1)/(2*jt + 1)*(D**2)*(5.217721092 * 10**(-11))**(2) * (1.602*10**(-19))**2 *1/(math.sqrt(2*j + 1))

atom=Caesium()
    
file = "transitions.csv"
df = pd.read_csv(file, delimiter=";")
df.columns = ['wavelength', 'unc', 'wavenumber', "Int","DQ","D","J","Jn",'decaywidth', 'lowerlevel', 'upperlevel',"",""]


#preselection = df.loc[isNaN(df['decaywidth'])]
#selection = preselection.loc[isNaN(df['D'])]
columns = ['wavelength', 'unc', 'wavenumber', "Int","DQ","D","J","Jn",'decaywidth', 'lowerlevel', 'upperlevel',"",""]
newdf = pd.DataFrame(columns=columns)

for index, row in df.iterrows():
    if not (isNaN(row['decaywidth']) and isNaN(row['D'])):
        newdf.loc[len(newdf)] = row
        continue;
    
    lowerlevel = row["lowerlevel"]
    upperlevel = row["upperlevel"]
    try:
        wl = float(row["wavelength"])
        n,l,j = getQuantumNumbers(lowerlevel)
        nt,lt,jt = getQuantumNumbers(upperlevel)
    except:
        newdf.loc[len(newdf)] = row
        continue;
    #print (row["lowerlevel"], n, l ,j)
    #print (row["upperlevel"], nt, lt ,jt)
    #print(n,l,j,nt,lt,jt)
    DQ = atom.getReducedMatrixElementJ(n,l,j,nt,lt,jt)
    dw = DecayWidth(wl, j, jt,DQ)
    print(row["lowerlevel"],row["upperlevel"],n,l,j,nt,lt,jt,"res",DQ, dw)
    row["DQ"] = DQ
    row["decaywidth"] = dw
    newdf.loc[len(newdf)] = row
    #print (row["DQ"], row["decaywidth"])

print (df)
os.remove(output)
newdf.to_csv(output, sep=";", index=False)
    
    
    