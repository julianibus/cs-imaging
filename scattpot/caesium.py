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


class atom:
	
    qma = qm.QM()
	
    def __init__(self):
        dummy = 0
		
    def isNaN(self, num):
        return num != num
    
    def GetFactors(self, laser_wl, transition_ident, file):
        df = pd.read_csv(file, delimiter=";")
        df.columns = ['wavelength', 'unc', 'wavenumber', "Int","DQ","D","J","Jn",'decaywidth', 'lowerlevel', 'upperlevel',"hu","ha"]
        #levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
        levels = [transition_ident]
        for level in levels:
            
            
            tot_pot_factor = 0
            tot_scatt_factor = 0
            c = 0
            
            #Resonances upwards
            selection = df.loc[df['lowerlevel'] == level]
            if not (level == "6s 2 S 1/2"):
                for index, row in selection.iterrows():
                    if not self.isNaN(row["decaywidth"]):
                        c+=1
                        width = float(row["decaywidth"])
                        wl = float(row["wavelength"])
                        pot_factor = self.qma.PotPrefactor(laser_wl, width, wl * 10**-10)
                        #print (row["lowerlevel"],row["upperlevel"],pot_factor)
                        scatt_factor = self.qma.ScattPrefactor(laser_wl, width, wl * 10**-10)
                        tot_scatt_factor += scatt_factor
                        tot_pot_factor += pot_factor
            else:
                for index, row in selection.iterrows():
                    potex = 0
                    if not self.isNaN(row["decaywidth"]):
                        if not (row["upperlevel"] == "6p 2 P ?1/2") and not (row["upperlevel"] == "6p 2 P ?3/2"):
                            c+=1
                            width = float(row["decaywidth"])
                            wl = float(row["wavelength"])
                            pot_factor = self.qma.PotPrefactor(laser_wl, width, wl * 10**-10)
                            #print (row["lowerlevel"],row["upperlevel"],pot_factor)
                            scatt_factor = self.qma.ScattPrefactor(laser_wl, width, wl * 10**-10)
                            tot_scatt_factor += scatt_factor
                            tot_pot_factor += pot_factor
                        else:
                            width = float(row["decaywidth"])
                            pot_factor = self.qma.PotPrefactorAntiMagic(laser_wl, width, mf=2, pol=-1)
                            tot_pot_factor += pot_factor/2
                            potex= pot_factor
                #print("ex", potex/tot_pot_factor)
            #print(level, c, "lower levels")
                    #print(row["lowerlevel"],wl, width, pot_factor)
                    
            #Resonances downwards
            selection = df.loc[df['upperlevel'] == level]
            for index, row in selection.iterrows():
                if not self.isNaN(row["decaywidth"]):
                    c+=1
                    width = float(row["decaywidth"])
                    wl = float(row["wavelength"])
                    pot_factor = -1*self.qma.PotPrefactor(laser_wl, width, wl * 10**-10) ## add with minus to account for stark energy shift in opposite direction
                    scatt_factor = self.qma.ScattPrefactor(laser_wl, width, wl * 10**-10)
                    tot_scatt_factor += scatt_factor
                    tot_pot_factor += pot_factor                    
            #print(level, c, "upper levels")
            
            return [tot_pot_factor, tot_scatt_factor]
                
            
# add scattering rate
