# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 00:58:01 2019

@author: julian
"""

import numpy as np
from math import *
import math
import scipy.special
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz
from scipy.integrate import odeint
from matplotlib.colors import LogNorm
from scipy.linalg import expm
from numpy import linalg as LA
from scipy.misc import factorial
from scipy.special import *
from numpy.linalg import *
import caesium
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
import pandas as pd
import math

nmax=18

def IntDnormal(eta):
	N = np.arange(0,nmax,1)
	NT = np.arange(0,nmax,1)
	D = np.zeros((nmax, nmax),dtype=complex)
	for n in N:
		for nt in NT:
			D[n,nt] = sqrt(factorial(np.minimum(n, nt))/factorial(np.minimum(n, nt) + np.abs(nt - n)))*(1j*eta *np.cos(math.pi/4))**(np.abs(nt - n))*eval_genlaguerre(np.minimum(n, nt), np.abs(nt - n), (eta*math.cos(math.pi/4))**2)*np.exp(-0.5*(eta*np.cos(math.pi/4))**2)
	
	return D

plt.figure(figsize=(8,2.7))
D1 = IntDnormal(0.15)
D2 = IntDnormal(0.3)
D3 = IntDnormal(1)
Ds = [D1, D2, D3]
etas = [0.15, 0.3, 1]
for i in range(0,3):
    plt.subplot(1,3,i+1)
    plt.imshow(np.power(np.abs(Ds[i]),2),cmap="Reds",aspect="auto")
    plt.gca().text(0.64, 0.90, "$\eta$ = " + str(etas[i]), fontsize=10, transform=plt.gca().transAxes)
    plt.xlabel("n")
    plt.ylabel("n'")
    plt.clim(0,1)
    if i ==2:
        plt.colorbar(fraction=0.056, pad=0.04)
plt.tight_layout()
plt.show()