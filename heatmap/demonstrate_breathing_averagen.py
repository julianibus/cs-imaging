# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 18:31:59 2019

@author: julian
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 09:44:17 2019

@author: J.Wienand
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


import numpy as np
from math import *
import math
import scipy.special
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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
import pandas as pd
import math

c = 2.997*10**8;
hbar = 1.0545718*10**(-34)
m = 132.90 * 1.660539*10**(-27) 

nmax = 30

g23 = 4.05*10**6;
g34 = 6.23*10**6;
g35 = 11.4*10**6;
g41 = 28.6*10**6;
g26 = 0.13*10**6;
g27 = 1.10*10**6;
g75 = 0.78*10**6;
g65 = 0.11*10**6;
g51 = 32.8*10**6;
g64 = 0.91*10**6;
g21 = 1.84*10**6;

gammas = [g21, g23, g35, g41, g64, g51, g65, g75, g27, g26, g34]
gammaslabel = ["g21","g23", "g35", "g41", "g64", "g51", "g65", "g75", "g27", "g26", "g34"]
overlaps = list()
avgoverlaps = list()
totoverlaps = list()

def psix(x,n,omega):
    expression = 1/(np.sqrt(2**n*factorial(n)))*((m*omega)/(np.pi*hbar))**(1/4)*np.exp((-m*omega*x**2)/(2*hbar))*scipy.special.eval_hermite(n, np.sqrt((m*omega)/hbar)*x)
    return expression/LA.norm(expression)

xstep = 0.4 * 10**(-8)
x = np.arange(-0.3*10**(-6), +0.3*10**(-6), xstep)
X,Y = np.meshgrid(x,x)
top = toeplitz(np.concatenate([[-2,1], np.zeros(int(len(x)-2))]))
top[0,len(x)-1] = 1
top[len(x)-1,0] = 1
top = top/(xstep**2)
print(top)

def psit(psip,ta, omegat):
    Ht = (-((hbar**2)/(2*m))*top +0.5*m*omegat*np.abs(omegat) *np.square(np.diag(x)))#-((hbar**2)/(2*m))*top +
    B = (-1j/hbar)*Ht*ta
    #w, P = LA.eig(B)
    #A = np.diag(w)
    #EA = expm(A)
    #Pin = np.linalg.inv(P)
    #EB = np.dot(Pin,np.dot(EA,P))
    EB = expm(B)
    res = np.dot(EB,psip)
	#print (initial_evol, final)
    ns = np.arange(0,30,1)
    overlap = list()
    avgoverlap = 0
    for k in ns:
        ov = np.dot(psix(x,k,2946567),res)
        overlap.append(ov)
    
    overlaps.append(overlap)
    print(((np.asarray(np.nan_to_num(np.square(np.abs(overlap))))* ns)))
    avgoverlap =  np.sum((np.asarray(np.nan_to_num(np.square(np.abs(overlap))))* (ns)))
    totoverlap = np.sum(np.asarray(np.nan_to_num(np.square(np.abs(overlap)))))
    #print (avgoverlap)
    avgoverlaps.append(avgoverlap)
    totoverlaps.append(totoverlap)
    
	
    #print(mat)
    return res

initial = psix(x,2 ,2946567)

tstart = 0*10**(-9)
tstop = 4000*10**(-9)
tstep = 100*10**(-9)

timevol = list()
T = np.arange(tstart, tstop, tstep)
for td in np.arange(tstart, tstop, tstep):
   # print ("hey")
    print (td)
    timevol.append(psit(initial, td,957940.036115241855))
    #print (psit(cpsi, td, 2*10**5))
    #print (tstep*psit(cpsi, td, 2*10**5))



colormap = plt.get_cmap('gist_heat')
plt.gca().set_color_cycle([colormap(k) for k in np.linspace(0, 0.8, 5)])
lines = plt.plot(T/10**-6, np.square(np.abs(overlaps)))
#plt.plot(T/10**-9, totoverlaps, "--b")
plt.figure(figsize=(8,5))
plt.subplot(2,1,1)
l0=plt.plot(np.asarray(T)*1e6, avgoverlaps, color="black", linewidth=1.5, label="$\langle n \\rangle\,(n=2)$")
B=avgoverlaps[16]-2
l1=plt.plot(np.asarray(T)*1e6, 2+B*np.square(np.sin(T*957940.036115241855)), "--r", label="$A sin^2(t \omega) + B$")
omp = 957940.036115241855
avgts = list()
avgns = list()
taus = np.arange(0.1e-6, 5e-6, 0.1e-6)
for tau in taus:
    avgts.append(np.sqrt(2)*tau/np.sqrt(1+4*tau**2*omp**2))
    avgns.append(2 + (2*B*tau**2 *omp**2)/(1+4*tau**2*omp**2))
tau = 0.9e-6
plt.ylabel("$\langle n \\rangle$")
ax1 = plt.gca()
ax2 = plt.gca().twinx()
l2=plt.plot(np.asarray(T)*1e6, np.exp(-T/tau), "b", label="Exp Decay $\\tau=0.5\,\mu s$")
ax2.tick_params(labelcolor='blue')


plt.xlabel("t ($\mu s$)")
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.legend(handles=[l0[0],l1[0],l2[0]],frameon=False)

plt.subplot(2,1,2)
plt.plot(taus*1e6, avgns, color="green",label="$\langle \widetilde{n} \\rangle\,(n=2)$")
plt.xlabel("$\\tau\,(\mu s)$")
plt.ylabel("$\langle  \widetilde{n} \\rangle\,(n=2)$")
plt.tight_layout()
plt.show()