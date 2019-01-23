# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 16:35:38 2019

@author: J.Wienand
"""
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
import math
import numpy as np
import matplotlib.pyplot as plt
 

ICs = 7/2
q = 0

def trans_elem(J,Jt, F, mf, Ft, mft):
    q = (mf - mft)
    return (-1)**(Ft - 1 + mf) * math.sqrt(2*F + 1) *wigner_3j(Ft, 1, F, mft , q, -mf) * (-1)**(Ft+J+1+ICs)*math.sqrt((2*Ft+1)*(2*J+1))*wigner_6j(J, Jt, 1, Ft, F, ICs)

def transprobs(J,Jt, F, mf):
    
    choices = list()
    probs = list()
    total = 0
    Fts = np.arange(max(ICs - Jt, F-1), min(ICs + Jt+1,F+1), 1)
    for Fti in Fts:
        mfts = np.arange(max(-Fti, mf - 1), min(Fti, mf + 1) + 1, 1)
        for mfti in mfts:
            total_part = np.square(np.abs(trans_elem(J,Jt, F, mf, Fti, mfti)))
            total = total + total_part
            #print(Fti, mfti, "-> ",total_part)
    
    Fts = np.arange(max(ICs - Jt, F-1), min(ICs + Jt+1,F+1), 1)
    for Fti in Fts:
        mfts = np.arange(max(-Fti, mf - 1), min(Fti, mf + 1) + 1, 1)
        for mfti in mfts:
            selected = np.square(np.abs(trans_elem(J,Jt, F, mf, Fti, mfti)))
            choices.append((Fti, mfti))
            probs.append(selected/total)
    #print(sum(probs))
    return choices, probs

def decay_level_choicer(J,Jt, F, mf):
    ch, prob = transprobs(J,Jt, F, mf)
    choice = np.random.choice(range(0,len(ch)),1, p=prob)
    return ch[choice[0]]


clist = list()
for f in range(0,1000):
    c = decay_level_choicer(5/2,5/2, 3,-1)
    clist.append(c)
    
unique, counts = np.unique(clist, return_counts=True, axis=0)
print(unique, counts)

y_pos = np.arange(len(unique))
 
plt.bar(y_pos, counts, align='center', alpha=0.5)
plt.xticks(y_pos, unique)