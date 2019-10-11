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
x = np.arange(-0.15*10**(-6), +0.15*10**(-6), xstep)
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
    ns = np.arange(0,10,2)
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
tstep = 10*10**(-9)

timevol = list()
T = np.arange(tstart, tstop, tstep)
for td in np.arange(tstart, tstop, tstep):
   # print ("hey")
    print (td)
    timevol.append(psit(initial, td,957940.036115241855))
    #print (psit(cpsi, td, 2*10**5))
    #print (tstep*psit(cpsi, td, 2*10**5))


plt.figure(figsize=(8,4.9))
plt.subplot(2,1,1)
plt.imshow(np.transpose(np.square(np.abs(timevol))), aspect="auto", cmap="Greys")
plt.xlabel("t ($\mu$s)")
plt.ylabel("x (nm)")
plt.xlim(0,400)
plt.xticks(plt.xticks()[0], [str(t/100) for t in plt.xticks()[0]])
plt.yticks(plt.yticks()[0], [str(t*4- 150) for t in plt.yticks()[0]])
plt.ylim(0, 74)
plt.colorbar(fraction=0.056, pad=0.04, label="Density (arb. u.)")

plt.axvline(x=54, color="red")
plt.axvline(x=0.2*54, color="red")

plt.gca().text(0.14, 0.04, "$\\tau $", fontsize=10, color="red", transform=plt.gca().transAxes)
plt.gca().text(0.03, 0.04, "$0.2 \\tau $", fontsize=10, color="red", transform=plt.gca().transAxes)
plt.gca().text(0.90, 0.85, "$n=2$", fontsize=10, transform=plt.gca().transAxes)

    
plt.subplot(2,3,4)
colormap = plt.get_cmap('gist_heat')
plt.gca().set_color_cycle([colormap(k) for k in np.linspace(0, 0.8, 5)])
lines = plt.plot(T/10**-6, np.square(np.abs(overlaps)))
#plt.plot(T/10**-9, totoverlaps, "--b")
plt.legend(lines[:len(lines)], ["n = 0","n = 2","n = 4","n = 6","n = 8"], loc=0, fontsize=8,frameon=False)
plt.ylabel("$| \langle n\,| \,U(t)\,|\,2 \\rangle |^2$")
plt.xlabel("t ($\mu s$)")
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

#plt.figure(figsize=(7,5))
#plt.plot(T, avgoverlaps, color="black", linewidth=1.5)
#plt.plot(T, 0+2*np.square(np.sin(T*214198.036115241855)), "--r")
#plt.ylabel("<n'>")
#plt.xlabel("Time (s)")
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

#####################################################


nmax = 16
lambdal = 880
inten = 1*10**9 
xstep = 0.2 * 10**(-8)
x = np.arange(-0.15*10**(-6), +0.15*10**(-6), xstep)

thetanormal = 0

thetastep = math.pi/50
thetas = np.arange(0, math.pi+thetastep, thetastep)

timeintegration = False
timesN = 10


###################

cs = caesium.atom()

####### PART 1: CONSTANT ########
c = 2.997*10**8;
hbar = 1.0545718*10**(-34)
m = 132.90 * 1.660539*10**(-27) 


g23 = 4.05*10**6
g34 = 6.23*10**6
g35 = 11.4*10**6
g41 = 28.6*10**6
g26 = 0.13*10**6
g27 = 1.10*10**6
g75 = 0.78*10**6
g65 = 0.11*10**6
g51 = 32.8*10**6
g64 = 0.91*10**6
g21 = 1.84*10**6



l21 = 455.5 * 10**-9 
l23 = 2931.8 * 10**-9 
l35 = 1469.9  * 10**-9 
l41 = 894.3  * 10**-9 
l64 = 3011.1  * 10**-9
l51 = 852.1 * 10**-9
l65 = 3614.1  *10**-9 
l75 = 3491  * 10**-9 
l27 = 1360.6  * 10**-9 
l26 = 1342.8  * 10**-9 
l34 = 1359.2 * 10**-9 

decaylabels = ["|2> -> |1>", "|2> -> |3>", "|3> -> |5>", "|4> -> |1>", "|6> -> |4>", "|5> -> |1>", "|6> -> |5>", "|7> -> |5>", "|2> -> |7>", "|2> -> |6>", "|3> -> |4>"]
decays = [l21, l23, l35, l41, l64, l51, l65, l75, l27, l26, l34]
levels = ["6s 2 S 1/2", "7p 2 P ?3/2", "7s 2 S 1/2", "6p 2 P ?1/2", "6p 2 P ?3/2", "5d 2 D 3/2", "5d 2 D 5/2"]
decaylevels = ["6s 2 S 1/2", "7s 2 S 1/2","6p 2 P ?3/2", "6s 2 S 1/2", "6p 2 P ?1/2","6s 2 S 1/2", "6p 2 P ?3/2", "6p 2 P ?3/2", "5d 2 D 5/2","5d 2 D 3/2", "6p 2 P ?1/2"]
gammas = [g21, g23, g35, g41, g64, g51, g65, g75, g27, g26, g34]
gammaslabel = ["g21","g23", "g35", "g41", "g64", "g51", "g65", "g75", "g27", "g26", "g34"]

omegas = list()
omegas.append(10) #dummy for index 
lambdickes = list()
lambdickes.append(10) #dummy for index 
latticedepths1 = list() # in hbar omega
latticedepths2 = list() # in lattice recoil energies
for level in levels:
        #098765', '#000009'
        pots = list()
        scatts = list()
        pot, scatt = cs.GetFactors(lambdal*10**-9, level, "transitions_complemented.csv")
        #pot = pot*(-1)
        omega = math.sqrt(2*np.absolute(pot)*inten/m)*2 *math.pi /(lambdal * 10**(-9))
        if (pot < 0):
            omegas.append(-omega)
        else:
            omegas.append(omega)
    
        latticedepth1 = (math.fabs(pot)*inten)/(hbar*omega)
        latticedepths1.append(latticedepth1)


    
for decayi in range(0,len(decays)):
    pot, scatt = cs.GetFactors(lambdal*10**-9, decaylevels[decayi], "transitions_complemented.csv")
    pot = pot*(-1)
    omega = math.sqrt(math.fabs(pot)*inten/m)*2 *math.pi /(lambdal * 10**(-9))
    #if (pot < 0): omega = omega*(-1) #dont negate here as sign not important for calculating the lamb dicke parameter
    ELatRec = ((hbar*(2*math.pi)/(decays[decayi]))**2)/(2*m)
    lambdicke = math.sqrt(ELatRec /(hbar * omega))
    lambdickes.append(lambdicke)
    latticedepth2 = (math.fabs(pot)*inten)/(ELatRec)
    latticedepths2.append(latticedepth2)
            
print(omegas)
print(lambdickes)
print(latticedepths1)
print(latticedepths2)

etas = lambdickes





def psix(x,n,omega):
    expression = 1/(np.sqrt(2**n*factorial(n)))*((m*omega)/(np.pi*hbar))**(1/4)*np.exp((-m*omega*x**2)/(2*hbar))*scipy.special.eval_hermite(n, np.sqrt((m*omega)/hbar)*x)
    return expression/LA.norm(expression)


X,Y = np.meshgrid(x,x)
top = toeplitz(np.concatenate([[-2,1], np.zeros(int(len(x)-2))]))
top[0,len(x)-1] = 1
top[len(x)-1,0] = 1
top = top/(xstep**2)

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


def IntDr(eta):
    	N = np.arange(0,nmax,1)
    	NT = np.arange(0,nmax,1)
    	D = np.zeros((nmax, nmax),dtype=complex)
    	for n in N:
    		for nt in NT:
    			D[n,nt] = sqrt(factorial(np.minimum(n, nt))/factorial(np.minimum(n, nt) + np.abs(nt - n)))*(1j*eta *np.cos(math.pi/4))**(np.abs(nt - n))*eval_genlaguerre(np.minimum(n, nt), np.abs(nt - n), (eta*math.cos(math.pi/4))**2)*np.exp(-0.5*(eta*np.cos(math.pi/4))**2)
    	
    	return D

	
def IntD(eta):
    N = np.arange(0,nmax,1)
    NT = np.arange(0,nmax,1)
    D = np.zeros((nmax, nmax),dtype=complex)
    for n in N:
	    for nt in NT:
	        D[n,nt] = sqrt(factorial(np.minimum(n, nt))/factorial(np.minimum(n, nt) + np.abs(nt - n)))*(1j*eta *np.cos(thetanormal))**(np.abs(nt - n))*eval_genlaguerre(np.minimum(n, nt), np.abs(nt - n), (eta*math.cos(thetanormal))**2)*np.exp(-0.5*(eta*np.cos(thetanormal))**2)
    #for n in N:
    #    D[n] = D[n]/LA.norm(D[n]) #fake normalization
	
    return D
	
#tstart = 0.1*10**(-9)
#tstop = 1000*10**(-9)
#tstep = 10*10**(-9)
#T = np.arange(tstart, tstop, tstep)

def Unn(n,nt, omega, omegat, td):
	#print (td,omega, omegat)
	initial = psix(x,n,omega)
	final = psix(x,nt,omega)
	initial_evol = psit(initial, td, omegat)
	#print (initial_evol, final)
	overlap = np.dot(final,initial_evol)
	print (overlap)
	return overlap
	
	
	

def IntU(omega, omegat, ti):
    
    if (timeintegration == False): 
        U = np.zeros((nmax, nmax),dtype=complex)
        N = np.arange(0,nmax,1)
        NT = np.arange(0,nmax,1)
        for nl in N:
            print (nl)
            for ntl in NT:
                U[nl,ntl] = Unn(nl,ntl, omega , omegat, ti)
        #for nl in N:
        #    U[nl] = U[nl]/LA.norm(U[nl])
            
        return np.copy(U)
    else:
        dummy = 0
        
	
#print (omegas[1], omegas[5], 1/g65)
#U65 = IntU(omegas[1], omegas[5], (1/g65)) 
#np.savetxt("U65.csv", np.square(np.abs(U65)), delimiter=",")
    
def eff_decay_time(tau, omeg):
    return tau
    #return np.sqrt(2)*tau/(np.sqrt(1+4*tau**2 * omeg**2))
    #if (omeg > 0):
     #   breathing = (2 * math.pi/((omeg)) / 2)
    #    redped = tau % (2 * math.pi/((omeg)) / 2)
    #    if (tau > (breathing/2)):
    #        return breathing - redped
     #   else:
    #        return redped
    #else:
     #   return tau     
        
plt.subplot(2,3,6)
U51 = IntU(omegas[1], omegas[2], eff_decay_time(1/g21, omegas[2]))
plt.gca().text(0.64, 0.90, "$t = \\tau $", fontsize=10, transform=plt.gca().transAxes)
plt.imshow(np.power(np.abs(U51),2), cmap="Reds",aspect="auto")
plt.colorbar(fraction=0.056, pad=0.04, label="$|\langle n'\,| \,U(t)\,|\,n \\rangle |^2$")
plt.clim(0,1)
plt.ylabel("n'")
plt.xlabel("n")

plt.subplot(2,3,5)
U51 = IntU(omegas[1], omegas[2], eff_decay_time(0.2*1/g21, omegas[2]))
plt.gca().text(0.64, 0.90,  "$t = 0.2 \\tau $", fontsize=10, transform=plt.gca().transAxes)
plt.imshow(np.power(np.abs(U51),2), cmap="Reds",aspect="auto")
plt.ylabel("n'")
plt.xlabel("n")


plt.tight_layout()
plt.show()