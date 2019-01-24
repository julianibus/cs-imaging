# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 16:35:38 2019

@author: J.Wienand
"""
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
import math
import numpy as np
import matplotlib.pyplot as plt
 

ICs = 7/2
q = 0

def trans_elem(J,Jt, F, mf, Ft, mft):
    q = (mf - mft)
    return (-1)**(Ft - 1 + mf) * math.sqrt(2*F + 1) *wigner_3j(Ft, 1, F, mft , q, -mf) * (-1)**(Ft+J+1+ICs)*math.sqrt((2*Ft+1)*(2*J+1))*wigner_6j(J, Jt, 1, Ft, F, ICs)

def transprobs(J,Jt, F, mf):
    
    choices = list()
    probs = list()
    total = 0
    Fts = np.arange(max(ICs - Jt, F-1), min(ICs + Jt+1,F+1)+1, 1)
    for Fti in Fts:
        mfts = np.arange(max(-Fti, mf - 1), min(Fti, mf + 1) + 1, 1)
        for mfti in mfts:
            total_part = np.square(np.abs(trans_elem(J,Jt, F, mf, Fti, mfti)))
            total = total + total_part
            #print(Fti, mfti, "-> ",total_part)
    
    Fts = np.arange(max(ICs - Jt, F-1), min(ICs + Jt+1,F+1)+1, 1)
    for Fti in Fts:
        mfts = np.arange(max(-Fti, mf - 1), min(Fti, mf + 1) + 1, 1)
        for mfti in mfts:
            selected = np.square(np.abs(trans_elem(J,Jt, F, mf, Fti, mfti)))
            choices.append((Fti, mfti))
            probs.append(selected/total)
    #print(sum(probs))
    return choices, probs

def decay_level_choicer(J,Jt, F, mf):
    ch, prob = transprobs(J,Jt, F, mf)
    choice = np.random.choice(range(0,len(ch)),1, p=prob)
    return ch[choice[0]]

def mfhisto(F, mf):
    clist = list()
    
    for f in range(0,50):
        paths = np.arange(0, 6,1)
       # print(paths)
        path = np.random.choice(paths,1, p=[0.201004, 0.367810,0.001969,0.016289,0.15449,0.258427]/np.sum([0.201004, 0.367810,0.001969,0.016289,0.15449,0.258427]))[0]
        try:
            if path == 0:
                c = decay_level_choicer(1/2,1/2, *decay_level_choicer(1/2, 1/2,*decay_level_choicer(3/2,1/2,F,mf)))
            elif path == 1:
                c = decay_level_choicer(3/2,1/2, *decay_level_choicer(1/2, 3/2,*decay_level_choicer(3/2,1/2,F,mf)))
            elif path == 2:
                c = decay_level_choicer(3/2,1/2, *decay_level_choicer(3/2, 3/2,*decay_level_choicer(3/2,3/2,F,mf)))
            elif path == 3:
                c = decay_level_choicer(1/2,1/2, *decay_level_choicer(3/2, 1/2,*decay_level_choicer(3/2,3/2,F,mf)))
            elif path == 4:
                c = decay_level_choicer(3/2,1/2, *decay_level_choicer(5/2, 3/2,*decay_level_choicer(3/2,5/2,F,mf)))
            elif path == 5:
                c = decay_level_choicer(3/2,1/2,F,mf)
                
            clist.append(c)
           # print("yes decay", path)
        except ValueError:
            print("no decay", path)
        
        
        
        
     
    if len(clist) > 0:
        unique, counts = np.unique(clist, return_counts=True, axis=0)
        print(unique, counts)
        
        grounds_posdict = {
                (3,-3): 0,
                (3,-2): 1,
                (3,-1): 2,            
                (3,0): 3,    
                (3,1): 4,    
                (3,2): 5,    
                (3,3): 6,    
                (4,-4): 7,
                (4,-3): 8,
                (4,-2): 9,
                (4,-1): 10,            
                (4,0): 11,    
                (4,1): 12,    
                (4,2): 13,    
                (4,3): 14,    
                (4,4): 15,            
				}
        out_array = np.zeros(16)
        counter = 0
        for un in unique:
            out_array[grounds_posdict[(un[0],un[1])]] = counts[counter]/np.sum(counts)
            counter += 1
        
        
        
    #    y_pos = np.arange(len(unique))
    #    plt.figure(figsize=(20,10))
    #    plt.bar(y_pos, counts, align='center', alpha=0.5)
    #    plt.xticks(y_pos, unique)
    #    plt.show()
        
        print(out_array)
        return out_array
    else:
        return "SHIT"
    
exciteds = [(2,-2),(2,-1),(2,0),(2,1),(2,2),(3,-3),(3,-2),(3,-1),(3,0),(3,1),(3,2),(3,3),(4,-4),(4,-3),(4,-2),(4,-1),(4,0),(4,1),(4,2),(4,3),(4,4),(5,-5),(5,-4),(5,-3),(5,-2),(5,-1),(5,0),(5,1),(5,2),(5,3),(5,4),(5,5)]
grounds = [(3,-3),(3,-2),(3,-1),(3,0),(3,1),(3,2),(3,3),(4,-4),(4,-3),(4,-2),(4,-1),(4,0),(4,1),(4,2),(4,3),(4,4)]
grounds_strings = ["(3,-3)","(3,-2)","(3,-1)","(3,0)","(3,1)","(3,2)","(3,3)","(4,-4)","(4,-3)","(4,-2)","(4,-1)","(4,0)","(4,1)","(4,2)","(4,3)","(4,4)"]
decayb_matrix = list() 

for ex in exciteds:
    row = mfhisto(*ex)
    decayb_matrix.append(row)

decay_matrix = np.transpose(np.asarray(decayb_matrix))


exciteds_posdict = {
        (2,-2): 0,
        (2,-1): 1,            
        (2,0): 2,    
        (2,1): 3,    
        (2,2): 4,    
        (3,-3): 5,
        (3,-2): 6,
        (3,-1): 7,            
        (3,0): 8,    
        (3,1): 9,    
        (3,2): 10,    
        (3,3): 11,    
        (4,-4): 12,
        (4,-3): 13,
        (4,-2): 14,
        (4,-1): 15,            
        (4,0): 16,    
        (4,1): 17,    
        (4,2): 18,    
        (4,3): 19,    
        (4,4): 20,   
        (5,-5): 21,
        (5,-4): 22,
        (5,-3): 23,
        (5,-2): 24,
        (5,-1): 25,            
        (5,0): 26,    
        (5,1): 27,    
        (5,2): 28,    
        (5,3): 29,    
        (5,4): 30, 
        (5,5): 31           
        }


def create_excitation_matrix(mode, deltaf):
    ematrix = list()
    if mode == "pi":
        for gstate in grounds:
            estate = (gstate[0]+deltaf,gstate[1] )
            row = np.zeros(len(exciteds_posdict))
            try:
                estate_index = exciteds_posdict[estate]
                row[estate_index] = 1
            except Exception:
                dummy = 0
            ematrix.append(row)
    elif mode == "sigma+":
        for gstate in grounds:
            estate = (gstate[0]+deltaf,gstate[1]+1)
            row = np.zeros(len(exciteds_posdict))
            try:
                estate_index = exciteds_posdict[estate]
                row[estate_index] = 1
            except Exception:
                dummy = 0
            ematrix.append(row)
    elif mode == "sigma-":
        for gstate in grounds:
            estate = (gstate[0]+deltaf,gstate[1]-1)
            row = np.zeros(len(exciteds_posdict))
            try:
                estate_index = exciteds_posdict[estate]
                row[estate_index] = 1
            except Exception:
                dummy = 0
            
            ematrix.append(row)
    return np.transpose(ematrix)
    
    
initial_state = np.zeros(len(grounds))
initial_state[len(grounds) -1] = 1

############### Pi light ##############

pil = np.dot(decay_matrix,np.dot(create_excitation_matrix("sigma+", 1),initial_state))
print(np.asarray(grounds), pil)
plt.figure(figsize=(10,5))
plt.bar(np.arange(0,16,1), pil, align='center', alpha=0.5)
plt.xticks(np.arange(0,16,1), grounds_strings)
plt.show()
