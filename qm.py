# -*- coding: utf-8 -*-
import math

c = 2.997*10**8


class QM:
    
    #U(x) = PotPrefactor*I(x)
    #assuming large detunings & neglegible saturatio
    def PotPrefactor(laser_wl, transition_gamma, transition_wl):
        transition_freq = 2*math.pi/transition_wl
        laser_freq = 2*math.pi/laser_wl
        return -(3*math.pi*c**2)/(2*(transition_freq**3))*(transition_gamma/(transition_freq-laser_freq)-transition_gamma/(transition_freq+laser_freq))

    def ScattPrefactor(laser_wl, transition_gamma, transition_wl):
        transition_freq = 2*math.pi/transition_wl
        laser_freq = 2*math.pi/laser_wl
        return (3*math.pi*c**2)/(2*(transition_freq**3))*(transition_gamma/(transition_freq-laser_freq)-transition_gamma/(transition_freq+laser_freq))**2
		
        
        

#print (QM.PotPrefactor(300*10**-9, 10**6,299*10**-9))
        
    #add scattering rate