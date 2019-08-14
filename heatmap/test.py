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


c = 2.997*10**8;
hbar = 1.0545718*10**(-34)
m = 132.90 * 1.660539*10**(-27) 

#CONFIGURATION####
nmax = 10
lambdal = 880
inten = 5*10**8
xstep = 0.25  * 10**(-8)
x = np.arange(-0.25*10**(-6), +0.25*10**(-6), xstep)

thetastep = math.pi/50
thetas = np.arange(0, math.pi+thetastep, thetastep)

timeintegration = True
timesN = 40
###################

X,Y = np.meshgrid(x,x)
top = toeplitz(np.concatenate([[-2,1], np.zeros(int(len(x)-2))]))
top[0,len(x)-1] = 1
top[len(x)-1,0] = 1
top = top/(xstep**2)

def IntDr(eta):
    N = np.arange(0,nmax,1)
    NT = np.arange(0,nmax,1)
    Deff = np.zeros((nmax, nmax),dtype=complex)
    for theta in thetas:
    	D = np.zeros((nmax, nmax),dtype=complex)
     
    	for n in N:
    		for nt in NT:
    			D[n,nt] = sqrt(factorial(np.minimum(n, nt))/factorial(np.minimum(n, nt) + np.abs(nt - n)))*(1j*eta *np.cos(math.pi/4))**(np.abs(nt - n))*eval_genlaguerre(np.minimum(n, nt), np.abs(nt - n), (eta*math.cos(math.pi/4))**2)*np.exp(-0.5*(eta*np.cos(math.pi/4))**2)
    	Deff = Deff +  D*np.sin(theta)/2
    return Deff*thetastep



def IntDnormal(eta):
	N = np.arange(0,nmax,1)
	NT = np.arange(0,nmax,1)
	D = np.zeros((nmax, nmax),dtype=complex)
	for n in N:
		for nt in NT:
			D[n,nt] = sqrt(factorial(np.minimum(n, nt))/factorial(np.minimum(n, nt) + np.abs(nt - n)))*(1j*eta *np.cos(math.pi/4))**(np.abs(nt - n))*eval_genlaguerre(np.minimum(n, nt), np.abs(nt - n), (eta*math.cos(math.pi/4))**2)*np.exp(-0.5*(eta*np.cos(math.pi/4))**2)
	
	return D


def eff_decay_time(tau, omeg):
    if (omeg > 0):
        breathing = (2 * math.pi/((omeg)) / 2)
        redped = tau % (2 * math.pi/((omeg)) / 2)
        if (tau > (breathing/2)):
            return breathing - redped
        else:
            return redped
    else:
        return tau      

def IntU(omega, omegat, ti):
    N = np.arange(0,nmax,1)
    NT = np.arange(0,nmax,1)
    timestep = 4*ti/timesN
    times = np.arange(0, 4*ti, timestep)
    
    Ueff = np.zeros((nmax, nmax),dtype=complex)
    
    for tim in times:
        U = np.zeros((nmax, nmax),dtype=complex)
        for nl in N:
            print (nl)
            for ntl in NT:
                U[nl,ntl] = Unn(nl,ntl, omega , omegat,  eff_decay_time(tim, omegat))
        
        Ueff = Ueff + U*np.exp(-tim/ti)/ti
            
    return np.copy(Ueff * timestep)

def IntUnormal(omega, omegat, ti):
    
    if (1==1): 
        U = np.zeros((nmax, nmax),dtype=complex)
        N = np.arange(0,nmax,1)
        NT = np.arange(0,nmax,1)
        for nl in N:
            print (nl)
            for ntl in NT:
                U[nl,ntl] = Unn(nl,ntl, omega , omegat, eff_decay_time(ti, omegat))
        for nl in N:
            U[nl] = U[nl]/LA.norm(U[nl])
            
        return np.copy(U)
    

def Unn(n,nt, omega, omegat, td):
	#print (td,omega, omegat)
	initial = psix(x,n,omega)
	final = psix(x,nt,omega)
	initial_evol = psit(initial, td, omegat)
	#print (initial_evol, final)
	overlap = np.dot(final,initial_evol)
	#print (overlap)
	return overlap

def psit(psip,t, omegat):
    Ht = (-((hbar**2)/(2*m))*top +0.5*m*omegat*np.abs(omegat)*np.square(np.diag(x)))#-((hbar**2)/(2*m))*top +
    B = (-1j/hbar)*Ht*t
    #w, P = LA.eig(B)N = np.arange(0,nmax,1)
    #A = np.diag(w)
    #EA = expm(A)
    #Pin = np.linalg.inv(P)
    #EB = np.dot(Pin,np.dot(EA,P))
    EB = expm(B)
    #plt.imshow((np.imag(Ut)),norm = LogNorm())
    res = np.dot(EB,psip)
    #print(mat)
    return res

def psix(x,n,omega):
    expression = 1/(np.sqrt(2**n*factorial(n)))*((m*omega)/(np.pi*hbar))**(1/4)*np.exp((-m*omega*x**2)/(2*hbar))*scipy.special.eval_hermite(n, np.sqrt((m*omega)/hbar)*x)
    return expression/LA.norm(expression)




plt.figure(1)
A = IntU(10**6, 1.0*10**6, 10**(-6))
plt.imshow(np.square(np.abs(A)))

plt.figure(2)
B = IntUnormal(10**6, 1.0*10**6, 10**(-6))
plt.imshow(np.square(np.abs(B)))
print(np.square(np.abs(A[0])))
print(np.square(np.abs(B[0])))
print(np.sum(np.square(np.abs(A[0]))))
print(np.sum(np.square(np.abs(B[0]))))
print(thetas)