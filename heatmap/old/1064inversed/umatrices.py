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

c = 2.997*10**8;
hbar = 1.0545718*10**(-34)
m = 132.90 * 1.660539*10**(-27) 

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

def uprop(x,xt,t,omegat):
    return ((m*omegat)/(2*np.pi*1j*hbar*np.sin(omegat*t)))**(1/2) * np.exp((1j*m*omegat)/(2*hbar*np.sin(omegat*t))*((x**2 + xt**2)*np.cos(omegat*t) - 2*x*xt))

def psi(x,n,omega):
    return 1/(np.sqrt(2**n*math.factorial(n)))*((m*omega)/(np.pi*hbar))**(1/4)*np.exp((-m*omega*x**2)/(2*hbar))*scipy.special.eval_hermite(n, np.sqrt((m*omega)/hbar)*x)

xstep = 10*10**(-9)
x = np.arange(-10**(-6), +10**(-6), xstep)
X,Y = np.meshgrid(x,x)


#print(x)
#plt.plot(x, psi(x,16,10**5))

#Vpsi = psi(x,3,10**5)

#Mupro = uprop(X,Y,10**(-7),10**5)
#print(Mupro)
#plt.imshow(np.real(Mupro))

#t = np.arange(0,500*10**(-9), 10*10**(-9))

def psit(n,omega,omeganew,t):
    return xstep*np.dot(uprop(X,Y,t,omeganew), psi(x,n,omega))

T = np.arange(80*10**(-9),2000*10**(-9), 10*10**(-9))
timevol= list()
for ti in T:
    vec = psit(3,10**5,5*10**5,ti)
    timevol.append(vec)
    print(vec)

plt.imshow(np.real(timevol),aspect = 1)
plt.axes().set_aspect('equal', 'datalim')
