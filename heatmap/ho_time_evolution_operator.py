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

c = 2.997*10**8;
hbar = 1.0545718*10**(-34)
m = 132.90 * 1.660539*10**(-27) 

nmax = 25

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

def psix(x,n,omega):
    expression = 1/(np.sqrt(2**n*factorial(n)))*((m*omega)/(np.pi*hbar))**(1/4)*np.exp((-m*omega*x**2)/(2*hbar))*scipy.special.eval_hermite(n, np.sqrt((m*omega)/hbar)*x)
    return expression/LA.norm(expression)

xstep = 0.5 * 10**(-8)
x = np.arange(-0.4*10**(-6), +0.4*10**(-6), xstep)
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
    ns = np.arange(0,21)
    overlap = list()
    avgoverlap = 0
    for k in ns:
        ov = np.dot(psix(x,k,658860.1901609472),res)
        overlap.append(ov)
    
    overlaps.append(overlap)
    print(((np.asarray(np.nan_to_num(np.square(np.abs(overlap))))* ns)))
    avgoverlap =  np.sum((np.asarray(np.nan_to_num(np.square(np.abs(overlap))))* (ns)))
    #print (avgoverlap)
    avgoverlaps.append(avgoverlap)
    
	
    #print(mat)
    return res

initial = psix(x,0 ,658860.1901609472)

tstart = 0.1*10**(-9)
tstop = 20000*10**(-9)
tstep = 400*10**(-9)

timevol = list()
T = np.arange(tstart, tstop, tstep)
for td in np.arange(tstart, tstop, tstep):
   # print ("hey")
    print (td)
    timevol.append(psit(initial, td,214198.036115241855))
    #print (psit(cpsi, td, 2*10**5))
    #print (tstep*psit(cpsi, td, 2*10**5))

plt.figure(figsize=(7,5))
plt.imshow(np.square(np.abs(timevol)))
plt.axes().set_aspect("auto")
plt.show()
plt.figure(figsize=(7,7))
lines = plt.plot(T, np.square(np.abs(overlaps)))
plt.axes().set_aspect("auto")
plt.legend(lines[:len(lines)], np.arange(0,len(lines),1))
plt.ylabel("Overlap")
plt.xlabel("Time (s)")
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.show()

plt.figure(figsize=(7,5))
plt.plot(T, avgoverlaps, color="black", linewidth=1.5)
plt.ylabel("<n'>")
plt.xlabel("Time (s)")
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.show()