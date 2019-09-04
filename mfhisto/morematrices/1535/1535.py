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

#CONFIGURATION####
 
lambdal = 1535
### INTENSITY MULTI ######
xstep = 0.25  * 10**(-8)
x = np.arange(-0.25*10**(-6), +0.25*10**(-6), xstep)
###################
RamanCooling = False
#################
 
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

gtot = 0.201*1/(1/g23+1/g34+1/g41)+0.3678*1/(1/g23+1/g35+1/g51)+0.001969*1/(1/g26+1/g65+1/g51)+0.016289*1/(1/g26+1/g64+1/g41)+0.15449*1/(1/g27+1/g75+1/g51)+0.258427*1/(1/g21)
print("gtot", gtot)



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

collector = list()
collector2 = list()
collectorall = list()



collectorscatt = list()
collectorbluescatt = list()
collectordepth1 = list()
collectordepth2 = list()
intensities = (10**(np.linspace(5,9, num=10)))
icount = 0
for inten in intensities:
    print ("INTENs", inten)
    omegas = list()
    omegas.append(10) #dummy for index 
    scatts = list()
    scatts.append(10) #dummy for index 
    lambdickes = list()
    lambdickes.append(10) #dummy for index 
    latticedepths1 = list() # in hbar omega
    latticedepths2 = list() # in lattice recoil energies
    for level in levels:
        #pots = list()
        #scatts = list()
        pot, scatt = cs.GetFactors(lambdal*10**-9, level, "transitions_complemented.csv")
        omega = math.sqrt(2*np.absolute(pot)*inten/m)*2 *math.pi /(lambdal * 10**(-9))
        if (pot < 0):
            omegas.append(omega)
        else:
            omegas.append(-omega)
        scatts.append(scatt*inten)
        latticedepth1 = (math.fabs(pot)*inten)/(hbar*omega)
        latticedepths1.append(latticedepth1)
    
    
    nmax = int(latticedepths1[0])+2
    print(nmax)
        
    for decayi in range(0,len(decays)):
        pot, scatt = cs.GetFactors(lambdal*10**-9, decaylevels[decayi], "transitions_complemented.csv")
        omega = math.sqrt(2*math.fabs(pot)*inten/m)*2 *math.pi /(lambdal * 10**(-9))
        #if (pot < 0): omega = omega*(-1) #dont negate here as sign not important for calculating the lamb dicke parameter
        ELatRec = ((hbar*(2*math.pi)/(lambdal*10**(-9)))**2)/(2*m)
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
    	        D[n,nt] = sqrt(factorial(np.minimum(n, nt))/factorial(np.minimum(n, nt) + np.abs(nt - n)))*(1j*eta *np.cos(0))**(np.abs(nt - n))*eval_genlaguerre(np.minimum(n, nt), np.abs(nt - n), (eta*math.cos(0))**2)*np.exp(-0.5*(eta*np.cos(0))**2)
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
    	#print (overlap)
    	return overlap
    	
    	
    	
    
    def IntU(omega, omegat, ti):
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
    	
    #print (omegas[1], omegas[5], 1/g65)
    #U65 = IntU(omegas[1], omegas[5], (1/g65)) 
    #np.savetxt("U65.csv", np.square(np.abs(U65)), delimiter=",")
        
    def eff_decay_time(tau, omeg):
        return np.sqrt(2)*tau/(np.sqrt(1+4*tau**2 * omeg**2))
        #if (omeg > 0):
        #    breathing = (2 * math.pi/((omeg)) / 2)
        #    redped = tau % (2 * math.pi/((omeg)) / 2)
        #    if (tau > (breathing/2)):
        #        return breathing - redped
        #    else:
        #        return redped
        #else:
        #    return tau        
    
    print("7/11")
    U65 = IntU(omegas[1], omegas[6], eff_decay_time(1/g65, omegas[6]))
    np.savetxt("U65.csv", np.square(np.abs(np.copy(U65))), delimiter=",")
    print("1/11")
    U21 = IntU(omegas[1], omegas[2], eff_decay_time(1/g21, omegas[2]))
    np.savetxt("U21.csv", np.square(np.abs(np.copy(U21))), delimiter=",")
    print("2/11")
    U23 = IntU(omegas[1], omegas[2], eff_decay_time(1/g23, omegas[2]))
    np.savetxt("U23.csv", np.square(np.abs(np.copy(U23))), delimiter=",")
    print("3/11")
    U35 = IntU(omegas[1], omegas[3], eff_decay_time(1/g35, omegas[3]))
    np.savetxt("U35.csv", np.square(np.abs(np.copy(U35))), delimiter=",")
    print("4/11")
    U41 = IntU(omegas[1], omegas[4], eff_decay_time(1/g41, omegas[4]))  
    np.savetxt("U41.csv", np.square(np.abs(np.copy(U41))), delimiter=",")
    print("5/11")
    U64 = IntU(omegas[1], omegas[6], eff_decay_time(1/g64, omegas[6])) 
    np.savetxt("U64.csv", np.square(np.abs(np.copy(U64))), delimiter=",")
    print("6/11")
    U51 = IntU(omegas[1], omegas[5], eff_decay_time(1/g51, omegas[5]))
    np.savetxt("U51.csv", np.square(np.abs(np.copy(U51))), delimiter=",")
    print("8/11")
    U75 = IntU(omegas[1], omegas[7], eff_decay_time(1/g75, omegas[7]))
    np.savetxt("U75.csv", np.square(np.abs(np.copy(U75))), delimiter=",")
    print("9/11")
    U27 = IntU(omegas[1], omegas[2],eff_decay_time(1/g27, omegas[2]))
    np.savetxt("U27.csv", np.square(np.abs(np.copy(U27))), delimiter=",")
    print("10/11")
    U26 = IntU(omegas[1], omegas[2], eff_decay_time(1/g26, omegas[2]))
    np.savetxt("U26.csv", np.square(np.abs(np.copy(U26))), delimiter=",")
    print("11/11")
    U34 = IntU(omegas[1], omegas[3], eff_decay_time(1/g34, omegas[3]))
    np.savetxt("U34.csv", np.square(np.abs(np.copy(U34))), delimiter=",")
    
    #Blue Imaging
    D12 = IntD(etas[1])
    D21 = IntDr(etas[1])
    D23 = IntDr(etas[2])
    D35 = IntDr(etas[3])
    D41 = IntDr(etas[4])
    D64 = IntDr(etas[5])
    D51 = IntDr(etas[6])
    D65 = IntDr(etas[7])
    D75 = IntDr(etas[8])
    D27 = IntDr(etas[9])
    D26 = IntDr(etas[10])
    D34 = IntDr(etas[11])
    np.savetxt("D12.csv", np.square(np.abs(D12)), delimiter=",")
    np.savetxt("D21.csv", np.square(np.abs(D21)), delimiter=",")
    np.savetxt("D23.csv", np.square(np.abs(D23)), delimiter=",")
    np.savetxt("D35.csv", np.square(np.abs(D35)), delimiter=",")
    np.savetxt("D41.csv", np.square(np.abs(D41)), delimiter=",")
    np.savetxt("D64.csv", np.square(np.abs(D64)), delimiter=",")
    np.savetxt("D51.csv", np.square(np.abs(D51)), delimiter=",")
    np.savetxt("D65.csv", np.square(np.abs(D65)), delimiter=",")
    np.savetxt("D75.csv", np.square(np.abs(D75)), delimiter=",")
    np.savetxt("D27.csv", np.square(np.abs(D27)), delimiter=",")
    np.savetxt("D26.csv", np.square(np.abs(D26)), delimiter=",")
    np.savetxt("D34.csv", np.square(np.abs(D34)), delimiter=",")
    #RedImaging
    D15 = IntD(etas[6])
    np.savetxt("D15.csv", np.square(np.abs(D15)), delimiter=",")
    

    Ra = np.identity(nmax) ##relict of failed implementation of raman sideband cooling matrix
    
    
    ################################
    
    pathA = np.dot(D41,np.dot(U41,np.dot(D34,np.dot(U34,np.dot(D23,np.dot(U23,np.dot(D12,Ra))))))) #(*7p3/2 \[Rule] 7s1/2 \[Rule] \6p1/2\[Rule] 6s1/2*)
    np.savetxt("A.csv", np.square(np.abs(pathA)), delimiter=",")
    pathB = np.dot(D51,np.dot(U51,np.dot(D35,np.dot(U35,np.dot(D23,np.dot(U23,np.dot(D12,Ra)))))))#(*7p3/2 \[Rule] 7s1/2 \[Rule] \6p3/2\[Rule] 6s1/2*)
    np.savetxt("B.csv", np.square(np.abs(pathB)), delimiter=",") 
    pathC = np.dot(D51,np.dot(U51,np.dot(D65,np.dot(U65,np.dot(D26,np.dot(U26,np.dot(D12,Ra)))))))#(*7p3/2 \[Rule] 5d3/2 \[Rule] \6p3/2\[Rule] 6s1/2*)
    np.savetxt("C.csv", np.square(np.abs(pathC)), delimiter=",")
    pathD = np.dot(D41,np.dot(U41,np.dot(D64,np.dot(U64,np.dot(D26,np.dot(U26,np.dot(D12,Ra)))))))#(*7p3/2 \[Rule] 5d3/2 \[Rule] \6p1/2\[Rule] 6s1/2*)
    np.savetxt("D.csv", np.square(np.abs(pathD)), delimiter=",")
    pathE = np.dot(D51,np.dot(U51,np.dot(D75,np.dot(U75,np.dot(D27,np.dot(U27,np.dot(D12,Ra)))))))#(*7p3/2 \[Rule] 5d5/2 \[Rule] \6p1/2\[Rule] 6s1/2*)
    np.savetxt("E.csv", np.square(np.abs(pathE)), delimiter=",")
    pathF = np.dot(D21,np.dot(U21,np.dot(D12,Ra)))#(*7p3/2 \[Rule] 6s1/2*) 
    np.savetxt("F.csv", np.square(np.abs(pathF)), delimiter=",")
    pathR = np.dot(D51,np.dot(U51,np.dot(D15,Ra)))#RED IMAGING 6p3/2 
    np.savetxt("R.csv", np.square(np.abs(pathR)), delimiter=",")

    Pathtot = np.add(0.258427*np.square(np.abs(pathF)),np.add(0.154494*np.square(np.abs(pathE)),np.add(0.016289*np.square(np.abs(pathD)),np.add(0.001969*np.square(np.abs(pathC)),np.add(0.367810*np.square(np.abs(pathB)),0.201004*np.square(np.abs(pathA)))))))
    PathtotCDE = np.add(0.25827/(0.2001004 + 0.3678 + 0.25827)*np.square(np.abs(pathF)),np.add(0.3678/(0.2001004 + 0.3678 + 0.25827)*np.square(np.abs(pathB)),0.201004/(0.2001004 + 0.3678 + 0.25827)*np.square(np.abs(pathA))))
    
    PathtotE = np.add(0.30564774229869446*np.square(np.abs(pathF)),np.add(0*np.square(np.abs(pathE)),np.add(0.019265386644210684*np.square(np.abs(pathD)),np.add(0.002328783000948544*np.square(np.abs(pathC)),np.add(0.43502588982218926*np.square(np.abs(pathB)),0.2377321982339569*np.square(np.abs(pathA)))))))
    
    #PathtotE = np.add(0.258427*np.square(np.abs(pathF)),np.add(0.154494*np.square(np.abs(pathE)),np.add(0.016289*np.square(np.abs(pathD)),np.add(0.001969*np.square(np.abs(pathC)),np.add(0.367810*np.square(np.abs(pathB)),0.201004*np.square(np.abs(pathA)))))))
    
    np.savetxt(str(icount) +"Pathtot.csv", Pathtot, delimiter=",")
    np.savetxt(str(icount) +"PathtotCDE.csv", PathtotCDE, delimiter=",")
    np.savetxt(str(icount) +"PathtotE.csv", PathtotE, delimiter=",") 
    np.savetxt(str(icount) +"PathtotR.csv", np.square(np.abs(pathR)), delimiter=",")
    
    def avgnt(n, pat):
        avg = 0
        if (RamanCooling == True):
            for k in range(0, nmax):
                if (k >= 1):
                    avg = avg + np.square(np.abs(pat[n-1,k]))*k
                else:
                    avg = avg + np.square(np.abs(pat[n-1,k]))*(k+1)
        else:
            for k in range(0, nmax):
                avg = avg + np.square(np.abs(pat[n-1,k]))*(k+1)
        return avg-n
        

    collector.append([avgnt(1, pathA),avgnt(1, pathB),avgnt(1, pathC),avgnt(1, pathD),avgnt(1, pathE),avgnt(1, pathF),avgnt(1, pathR)])
    collector2.append([avgnt(2, pathA),avgnt(2, pathB),avgnt(2, pathC),avgnt(2, pathD),avgnt(2, pathE),avgnt(2, pathF),avgnt(2, pathR)])
    
    collectorallpart = list()
    for z in range(1, nmax):
        collectorallpart.append([avgnt(z, pathA),avgnt(z, pathB),avgnt(z, pathC),avgnt(z, pathD),avgnt(z, pathE),avgnt(z, pathF),avgnt(z, pathR)])
    collectorall.append(collectorallpart)
    
    
    collectorscatt.append((0.254*scatts[1] + 0.254*scatts[2]+0.058*scatts[3]+0.014*scatts[4]+0.029*scatts[5]+0.032*scatts[6]+0.358*scatts[7])*((hbar*(2*math.pi)/(lambdal*10**(-9)))**2)/(2*m)) #mean lattice scattering assuming saturation of the imaging laser)
    collectorbluescatt.append(np.multiply([avgnt(1, pathA),avgnt(1, pathB),avgnt(1, pathC),avgnt(1, pathD),avgnt(1, pathE),avgnt(1, pathF),avgnt(1, pathR)],hbar*omegas[1]*gtot))
    collectordepth1.append(latticedepths1)
    collectordepth2.append(latticedepths2)
    print("latt",collectorscatt[len(collectorscatt)-1])
    print("blue",collectorbluescatt[len(collectorbluescatt)-1])
    
    np.savetxt("intensities.csv", intensities, delimiter=",")
    np.savetxt("heatings.csv", collector, delimiter=",")
    np.savetxt("heating2s.csv", collector2, delimiter=",")
    #np.savetxt("heatingsall.csv", collectorall, delimiter=",")
    np.savetxt("heatings_lattice.csv", collectorscatt, delimiter=",")
    np.savetxt("heatings_blue.csv", collectorbluescatt, delimiter=",")
    np.savetxt("latticedepths1.csv", collectordepth1, delimiter=",")
    np.savetxt("latticedepths2.csv", collectordepth2, delimiter=",")
    
    icount += 1