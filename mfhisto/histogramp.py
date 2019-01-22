# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 16:35:38 2019

@author: J.Wienand
"""
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
import math
import numpy as np

ICs = 7/2

def trans_elem(J,Jt, F, mf, Ft, mft):
    q = mft - mf
    return (-1)**(Ft - 1 + mf) * math.sqrt(2*F + 1) *wigner_3j(Ft, 1, F, mft , q, -mf) * (-1)**(Ft+J+1+ICs)*math.sqrt((2*Ft+1)*(2*J+1))*wigner_6j(J, Jt, 1, Ft, F, ICs)

def transprob(J,Jt, F, mf, Ft, mft):
    selected = np.square(np.abs(trans_elem(J,Jt, F, mf, Ft, mft)))
     
    total = 0
    Fts = np.arange(ICs - Jt, ICs - Jt, 1)
    for Fti in Fts:
        mfts = np.arange(-Fti, Fti, 1)
        for mfti in mfts:
            total_part = np.square(np.abs(trans_elem(J,Jt, Fti, mf, Ft, mfti)))
            print total_part
            total = total + np.square(np.abs(total_part))
    return selected/total