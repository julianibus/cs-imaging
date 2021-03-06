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

    def ScattPrefactor(self, laser_wl, transition_gamma, transition_wl):
        transition_freq = c*2*math.pi/transition_wl
        laser_freq = c*2*math.pi/laser_wl
        return (3*math.pi*c**2)/(2*hbar*(transition_freq**3))*(laser_freq/transition_freq)**3 * (transition_gamma/(transition_freq-laser_freq)+transition_gamma/(transition_freq+laser_freq))**2
		
        
        

#print (QM.PotPrefactor(300*10**-9, 10**6,299*10**-9))
        
    #add scattering rate
