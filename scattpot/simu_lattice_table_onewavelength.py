# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 13:19:13 2018

@author: J.Wienand
"""

import caesium
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
import pandas as pd

#config
lambdals = [532, 755, 880, 1064, 1535]
output = "532.csv"


levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
columns = ['laser wavelength','level', 'pot', 'scatt']
newdf = pd.DataFrame(columns=columns)
for lambdal in lambdals:
    for level in levels:
        pots = list()
        scatts = list()
        pot, scatt = caesium.atom.GetFactors(lambdal*10**-9, level, "transitions_complemented.csv")
        row = {'laser wavelength':lambdal,'level': level, 'pot': pot, 'scatt': scatt}
        newdf.loc[len(newdf)] = row
        
newdf.to_csv(output, sep=";", index=False)
    
    