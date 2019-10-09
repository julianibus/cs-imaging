# -*- coding: utf-8 -*-
import math

c = 2.997*10**8
hbar = 1.054571*10**(-34)

class QM:
    
    #U(x) = PotPrefactor*I(x)
    #assuming large detunings & neglegible saturatio
    def PotPrefactor(self, laser_wl, transition_gamma, transition_wl):
        transition_freq = c* 2*math.pi/transition_wl
        laser_freq = c* 2*math.pi/laser_wl
        return -(3*math.pi*c**2)/(2*(transition_freq**3))*(transition_gamma/(transition_freq-laser_freq)+transition_gamma/(transition_freq+laser_freq))
            
    def PotPrefactorAntiMagic(self, laser_wl, transition_gamma, mf, pol):
        d2freq = c* 2*math.pi/852e-9
        d1freq= c* 2*math.pi/894e-9
        gF = -0.25
        laser_freq = c* 2*math.pi/laser_wl
        return (math.pi*c**2 * transition_gamma)/(2*((d1freq)**3))*((2+pol*gF*mf)/(laser_freq-d2freq)+(1-pol*gF*mf)/(laser_freq-d1freq))     
            

    def ScattPrefactor(self, laser_wl, transition_gamma, transition_wl):
        transition_freq = c*2*math.pi/transition_wl
        laser_freq = c*2*math.pi/laser_wl
        return -(3*math.pi*c**2)/(2*hbar*(transition_freq**3))*(laser_freq/transition_freq)**3 * (transition_gamma/(transition_freq-laser_freq)+transition_gamma/(transition_freq+laser_freq))**2

    def ScattPrefactorAntiMagic(self, laser_wl, transition_gamma, mf, pol):
        d2freq = c* 2*math.pi/852e-9
        d1freq= c* 2*math.pi/894e-9
        laser_freq = c* 2*math.pi/laser_wl
        return (math.pi*c**2)/(2*hbar*(d1freq**3))*((2)/(laser_freq-d2freq)**2+(1)/(laser_freq-d1freq)**2)
            
        
        

#print (QM.PotPrefactor(300*10**-9, 10**6,299*10**-9))
        
    #add scattering rate
