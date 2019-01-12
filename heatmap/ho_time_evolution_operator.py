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

nmax = 15

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

def psix(x,n,omega):
    return 1/(np.sqrt(2**n*math.factorial(n)))*((m*omega)/(np.pi*hbar))**(1/4)*np.exp((-m*omega*x**2)/(2*hbar))*scipy.special.eval_hermite(n, np.sqrt((m*omega)/hbar)*x)


xstep = 0.1 * 10**(-7)
x = np.arange(-10**(-6), +10**(-6), xstep)
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
    #plt.imshow((np.imag(Ut)),norm = LogNorm())
    res = np.dot(EB,psip)
    #print(mat)
    return res

initial = psix(x,2,200000)

tstart = 0.1*10**(-9)
tstop = 4000*10**(-9)
tstep = 800*10**(-9)

timevol = list()
T = np.arange(tstart, tstop, tstep)
for td in np.arange(tstart, tstop, tstep):
   # print ("hey")
    print (td)
    timevol.append(psit(initial, td, 400000))
    #print (psit(cpsi, td, 2*10**5))
    #print (tstep*psit(cpsi, td, 2*10**5))
    
plt.imshow(np.real(np.square(np.absolute(timevol))))
plt.axes().set_aspect("auto")
plt.show()
