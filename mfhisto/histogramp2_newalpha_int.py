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
import matplotlib as mpl
import matplotlib.ticker as mtick

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

ICs = 7/2
q = 0

Natoms = 1000


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

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
            total_part = np.square(np.abs(trans_elem(J,Jt, int(F), int(mf), int(Fti), int(mfti))))
            total = total + total_part
            #print(Fti, mfti, "-> ",total_part)
    
    Fts = np.arange(max(ICs - Jt, F-1), min(ICs + Jt+1,F+1)+1, 1)
    for Fti in Fts:
        mfts = np.arange(max(-Fti, mf - 1), min(Fti, mf + 1) + 1, 1)
        for mfti in mfts:
            selected = np.square(np.abs(trans_elem(J,Jt, int(F), int(mf), int(Fti), int(mfti))))
            choices.append((Fti, mfti))
            probs.append(selected/total)
    #print(sum(probs))
    
  #  print (J, Jt, F, mf, choices)
    return choices, probs

def decay_level_choicer(J,Jt, F, mf):
    ch, prob = transprobs(J,Jt, F, mf)
    choice = np.random.choice(range(0,len(ch)),1, p=prob)
    return ch[choice[0]]

def mfhisto(F, mf):
    clist = list()
    
    for f in range(0,Natoms):
        paths = np.arange(0, 6,1)
       # print(paths)
        path = np.random.choice(paths,1, p=[0.201004, 0.367810,0.001969,0.016289,0.15449,0.258427]/np.sum([0.201004, 0.367810,0.001969,0.016289,0.15449,0.258427]))[0]
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
        
        
        
        
     
    if len(clist) > 0:
        unique, counts = np.unique(clist, return_counts=True, axis=0)
        #print(unique, counts)
        
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
exciteds_strings = ["(2,-2)","(2,-1)","(2,0)","(2,1)","(2,2)","(3,-3)","(3,-2)","(3,-1)","(3,0)","(3,1)","(3,2)","(3,3)","(4,-4)","(4,-3)","(4,-2)","(4,-1)","(4,0)","(4,1)","(4,2)","(4,3)","(4,4)","(5,-5)","(5,-4)","(5,-3)","(5,-2)","(5,-1)","(5,0)","(5,1)","(5,2)","(5,3)","(5,4)","(5,5)"]
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

def create_excitation_matrix(mode, deltaf,deltaframan):
    ematrix = list()
    
    if deltaframan == 1:
        nochangestate = 3
    elif deltaframan == -1:
        nochangestate = 4
        
    
    if mode == "pi": 
        for gstate in grounds: 
            estate = (gstate[0]+deltaf,gstate[1] )
            row = np.zeros(len(exciteds_posdict))
            try:
                estate_index = exciteds_posdict[estate]
                if gstate[0] != nochangestate:
                    row[estate_index] = 1
                else:
                    row[estate_index] = 0
            except Exception:
                dummy = 0
            ematrix.append(row)
    elif mode == "sigma+":
        for gstate in grounds:
            estate = (gstate[0]+deltaf,gstate[1]+1)
            row = np.zeros(len(exciteds_posdict))
            try:
                estate_index = exciteds_posdict[estate]
                if gstate[0] != nochangestate:
                    row[estate_index] = 1
                else:
                    row[estate_index] = 0
            except Exception:
                dummy = 0
            ematrix.append(row)
    elif mode == "sigma-":
        for gstate in grounds:
            estate = (gstate[0]+deltaf,gstate[1]-1)
            row = np.zeros(len(exciteds_posdict))
            try:
                estate_index = exciteds_posdict[estate]
                if gstate[0] != nochangestate:
                    row[estate_index] = 1
                else:
                    row[estate_index] = 0
            except Exception:
                dummy = 0
            
            ematrix.append(row)
    return np.transpose(ematrix)

Fcolors ={2: 'red',
          3: 'green',
          4: 'blue',
          5: 'orange'
        
        }

def raman_matrix(mode, deltaf):
    ematrix = list()
    if mode == "pi":
        for gstate in grounds:
            estate = (gstate[0]+deltaf,gstate[1] )
            row = np.zeros(len(grounds_posdict))
            try:
                estate_index = grounds_posdict[estate]
                row[estate_index] = 1
            except Exception:
                gstate_index = grounds_posdict[gstate]
                row[gstate_index] = 1
            ematrix.append(row)
    elif mode == "sigma+":
        for gstate in grounds:
            estate = (gstate[0] +deltaf,gstate[1]+1)
            row = np.zeros(len(grounds_posdict))
            try:
                estate_index = grounds_posdict[estate]
                row[estate_index] = 1
            except Exception:
                gstate_index = grounds_posdict[gstate]
                row[gstate_index] = 1
            ematrix.append(row)
    elif mode == "sigma-":
        for gstate in grounds:
            estate = (gstate[0]+deltaf,gstate[1]-1)
            row = np.zeros(len(grounds_posdict))
            try:
                estate_index = grounds_posdict[estate]
                row[estate_index] = 1
            except Exception:
                gstate_index = grounds_posdict[gstate]
                row[gstate_index] = 1
            
            ematrix.append(row)
    return np.transpose(ematrix)



def barplot(pil, strings,title, xtitle, ytitle):
    plt.figure(figsize=(10,8))
    barlist = plt.bar(np.arange(0,len(strings),1), pil, align='center', alpha=0.7,edgecolor='b')
    for ibar in range(0,len(barlist)):
        try:
            if len(strings)==4:
                barlist[ibar].set_color(Fcolors[int(strings[ibar][-1])])
            else:
                barlist[ibar].set_color(Fcolors[int(strings[ibar][1])])
        except Exception:
            barlist[ibar].set_color("b")
        
    ax = plt.axes()        
    ax.yaxis.grid()
    ax.set_axisbelow(True)
    plt.xticks(np.arange(0,len(strings),1), strings)
    plt.title(title)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.xticks(rotation=90) # horizontal lines
    
def barplotmany(pil, strings,title, xtitle, ytitle):
    plt.figure(figsize=(10,8))
    barlist = plt.bar(np.arange(0,len(strings),1), pil, align='center', alpha=0.7,edgecolor='b')
    for ibar in range(0,len(barlist)):
        try:
            if len(strings)==4:
                barlist[ibar].set_color(Fcolors[int(strings[ibar][-1])])
            else:
                barlist[ibar].set_color(Fcolors[int(strings[ibar][1])])
        except Exception:
            barlist[ibar].set_color("b")
        
    ax = plt.axes()        
    ax.yaxis.grid()
    ax.set_axisbelow(True)
    plt.xticks(np.arange(0,len(strings),1), strings)
    plt.title(title)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.xticks(rotation=90) # horizontal lines




sqrt(-1)

## EXPERIMENT 1: DECAY HISTOGRAMS
shares =list()
maxdecays = list()
for exi in exciteds:
    initial_state = np.zeros(len(exciteds))
    initial_state[exciteds_posdict[exi]] = 1
    out = np.dot(decay_matrix,initial_state) 
    #barplot(out,grounds_strings, "Decay from (F', m$_F$') = "  + str(exi), "Ground State (F, m$_F$)", "Population")
    sum3 = sum(out[0:7])
    sum4 = sum(out[7:16])
    if exi[1] == 0:
        share3 = sum3/(sum3 +sum4)
        shares.append(share3)
    maxdecays.append(max(out))

print (np.asarray(shares), exciteds_strings)
fs = ["F' = 2","F' = 3","F' = 4", "F' = 5"]
barplot(np.asarray(shares), fs, "", "", "Share of Decays into F = 3")
barplot(np.asarray(maxdecays), exciteds_strings, "", "Excited State (F, mF)", "Largest Decay Ratio")
    
## EXPERIMENT 2: Full Cycle

initial_state = np.zeros(len(grounds))
initial_state[11] = 1 #(4,4)

def timeevolution_movie(N, initial_state, repump_mode, repump_deltaf, raman_mode, raman_deltaf):    
    
    tot_probs = list()
    cool_shares = list()
    times = np.arange(0,N,1)
    
    newinitial = initial_state
    plt.clf()
    fig1 = plt.gcf()
    barplot(newinitial, grounds_strings, "","","")
    plt.savefig("-1.png")    
    for n in range(0,N):
        newinitial = np.dot(decay_matrix,np.dot(create_excitation_matrix(repump_mode, repump_deltaf,raman_deltaf),np.dot(raman_matrix(raman_mode,raman_deltaf),newinitial)))
        barplot(newinitial, grounds_strings, "","","")
        plt.savefig(str(n)+".png")
        tot_prob = (np.sum(newinitial))
        sum3 = sum(newinitial[0:7])
        sum4 = sum(newinitial[7:16])
        if raman_deltaf == 1:
            cool_share = sum3/(sum3 +sum4)
        elif raman_deltaf == -1:
            cool_share = sum4/(sum3 +sum4)
        tot_probs.append(tot_prob)
        cool_shares.append(cool_share)
        #print(tot_prob, cool_share)
    #plt.figure(figsize=(10,5))
    #plt.plot(times, tot_probs)
    #plt.figure(figsize=(10,5))
    #plt.plot(times, cool_shares)
    return tot_prob, cool_share

def timeevolution(N, initial_state, repump_mode, repump_deltaf, raman_mode, raman_deltaf):    
    
    tot_probs = list()
    cool_shares = list()
    times = np.arange(0,N,1)
    
    newinitial = initial_state
    
    for n in range(0,N):
        newinitial = np.dot(decay_matrix,np.dot(create_excitation_matrix(repump_mode, repump_deltaf,raman_deltaf),np.dot(raman_matrix(raman_mode,raman_deltaf),newinitial)))
        #barplot(newinitial, grounds_strings, "","","")
        tot_prob = (np.sum(newinitial))
        sum3 = sum(newinitial[0:7])
        sum4 = sum(newinitial[7:16])
        if raman_deltaf == 1:
            cool_share = sum3/(sum3 +sum4)
        elif raman_deltaf == -1:
            cool_share = sum4/(sum3 +sum4)
        tot_probs.append(tot_prob)
        cool_shares.append(cool_share)
        #print(tot_prob, cool_share)
    #plt.figure(figsize=(10,5))
    #plt.plot(times, tot_probs)
    #plt.figure(figsize=(10,5))
    #plt.plot(times, cool_shares)
    return tot_prob, cool_share, tot_probs, cool_shares

def ex2():
    mode_labels = {"sigma+": "$\sigma^{+}$",
                   "sigma-": "$\sigma^{-}$",
                   "pi":"  $\pi$"}
    deltaf_labels = {-1: "-1",
                   0: " 0",
                   1:" 1"}
    label_strings = list()
    eq_tot_probs = list()
    eq_cool_shares = list()
    #for repump_mode in ["pi", "sigma+", "sigma-"]:
    repump_mode = "pi"
    for repump_deltaf in [0,+1]:
                raman_mode = "pi"
            #for raman_mode in ["pi", "sigma+", "sigma-"]:
                raman_mode = "pi"
                for raman_deltaf in [-1,+1]:
                    if repump_mode == "pi" and repump_deltaf == 1 and raman_mode == "pi" and raman_deltaf==-1:
                        eq_tot_prob2, eq_cool_share2, tot_probs2, cool_shares2= timeevolution(10, initial_state, repump_mode, repump_deltaf, raman_mode, raman_deltaf) 
                    eq_tot_prob, eq_cool_share, tot_probs, cool_shares= timeevolution(80, initial_state, repump_mode, repump_deltaf, raman_mode, raman_deltaf)
                    if eq_tot_prob > 0:
                        label = str(deltaf_labels[repump_deltaf]) + ", " + str(deltaf_labels[raman_deltaf]) + ""
                        label_strings.append( label)
                        eq_tot_probs.append(eq_tot_prob)
                        eq_cool_shares.append(eq_cool_share)
    #barplotmany(eq_tot_probs, label_strings, "", "Configuration", "Prob")
    #barplotmany(eq_cool_shares, label_strings, "", "Configuration", "Equilibrium Share of Atoms being Cooled")
    strings = label_strings
    pil = eq_cool_shares
    plt.figure(figsize=(11,5))
    barlist = plt.bar(np.arange(0,len(strings),1), pil, align='center', alpha=0.7,edgecolor='b')
    for ibar in range(0,len(barlist)):
        try:
            if len(strings)==4:
                barlist[ibar].set_color(Fcolors[int(strings[ibar][-1])])
            else:
                barlist[ibar].set_color(Fcolors[int(strings[ibar][1])])
        except Exception:
            barlist[ibar].set_color("grey")
        
    ax = plt.axes()        
    ax.yaxis.grid()
    ax.set_axisbelow(True)
    plt.xticks(np.arange(0,len(strings),1), strings)
    plt.xlabel("Configuration")
    plt.ylabel("$\\alpha_{t\\to\infty}$")
    plt.xticks(rotation=90) # horizontal lines
    plt.tight_layout()
    plt.show()
    plt.figure(figsize=(11,5))
    plt.scatter(np.arange(0,len(cool_shares2),1),cool_shares2)
    plt.show()


## EXPERIMENT 3: TIME EVOLUTION INCLUDING TEMPERATURE
#Loading Matrix
def load_matrix(filename, cutoff):
    totmatrix = np.loadtxt(filename,delimiter=",")[0:cutoff,0:cutoff]
    #plt.matshow(totmatrix)
    rowsums = list()
    columnsums = list()
    for i in range(0,len(totmatrix)):
        rowsum = np.sum(totmatrix[i])
        columnsum = np.sum(totmatrix[:,i])
        #Stabilizer
        if (columnsum) > 1:
            totmatrix[:,i] = totmatrix[:,i]/np.sum(totmatrix[:,i])
            columnsum = np.sum(totmatrix[:,i])
        #print (i, columnsum)  
    return totmatrix
#matrix is loaded from heatmap folder -> first calculate totmatrix with heatpot_final.py


def full_time_evolution(No, ramannn, NRaman, NRepump, deltan, initial_state, initial_state_n, n_matrix, repump_mode, repump_deltaf, raman_mode, raman_deltaf):    
    #off_res_ra is now share of not cooled atoms during raman (n-> n)
    tot_probs = list()
    cool_shares = list()
        
    newinitial = initial_state 
    #print(newinitial)   
    newinitial1_n = initial_state_n  #State from which Raman cools
    newinitial2_n = np.zeros(len(initial_state_n)) #State into which Raman cools
    
    newinitial_n_nocooling = initial_state_n
    newinitial_n_bluephotons = 0
    newinitial_n_nocooling_bluephotons = 0
    newinitial1_n_sums = list()
    newinitial2_n_sums = list()
    newinitial_n_bluephotons_mon = list()
    newinitial_n_nocooling_bluephotons_mon = list()
    newinitial_n_nocooling_sums = list()
    for n in range(0,No):
        #äprint(newinitial)
        newinitial = np.dot(decay_matrix,np.dot(create_excitation_matrix(repump_mode, repump_deltaf, raman_deltaf),np.dot(raman_matrix(raman_mode,raman_deltaf),newinitial)))
        tot_prob = (np.sum(newinitial))
        sum3 = sum(newinitial[0:7])
        sum4 = sum(newinitial[7:16])
        if raman_deltaf == 1:
            cool_share = sum3/(sum3 +sum4)
        elif raman_deltaf == -1:
            cool_share = sum4/(sum3 +sum4)
        tot_probs.append(tot_prob)
        if n > 100:
            cool_shares.append(cool_shares[90])
            cool_share = cool_shares[90]
        else:
            cool_shares.append(cool_share)
         
        #1. COoling
        #for d in np.arange(0, NRaman):
            
            
        for d in np.arange(0, NRaman): 
            newinitial1_n_old = newinitial1_n
            newinitial2_n_old = newinitial2_n
            
            newinitial1_n_cooled = np.zeros(cutoff)
            for j in np.arange(0, len(newinitial1_n)-deltan):
                newinitial1_n_cooled[j] = (1-ramannn)*newinitial1_n[j+deltan] 
            #newinitial1_n_cooled[len(newinitial1_n)-1] = 0
            
            newinitial1_n_cooling = np.zeros(cutoff)
            for j in np.arange(deltan, len(newinitial1_n)):
                newinitial1_n_cooling[j] = (1-ramannn)*newinitial1_n[j]        
            
            newinitial1_n = newinitial1_n - newinitial1_n_cooling - ramannn*newinitial1_n
            newinitial2_n = newinitial2_n + newinitial1_n_cooled + ramannn*newinitial1_n_old
        #2.Heating
        for d in np.arange(0, NRepump):
            newinitial1_n_sum = np.sum(newinitial1_n)
            newinitial2_n_sum = np.sum(newinitial2_n)
            
            newinitial1_n = newinitial1_n + cool_share*np.dot(n_matrix, np.copy(newinitial2_n))
            newinitial2_n = (1-cool_share)*np.dot(n_matrix, np.copy(newinitial2_n))
            
            
            newinitial_n_bluephotons = newinitial_n_bluephotons + 0.25*(newinitial2_n_sum/(newinitial1_n_sum+newinitial2_n_sum))      
        #newinitial1_n = newinitial1_n + cool_share*np.dot(n_matrix, newinitial2_n)- newinitial1_n_cooling - ramannn*newinitial1_n
        #newinitial2_n = newinitial1_n_cooled + (1-cool_share)*np.dot(n_matrix, newinitial2_n) + ramannn*newinitial1_n
        
        newinitial1_n_sum = np.sum(newinitial1_n)
        newinitial2_n_sum = np.sum(newinitial2_n)
        
        newinitial_n_bluephotons_mon.append(newinitial_n_bluephotons)
        newinitial_n_nocooling_bluephotons = newinitial_n_nocooling_bluephotons + 0.25
        newinitial_n_nocooling_bluephotons_mon.append(newinitial_n_nocooling_bluephotons)
        
        newinitial1_n_sums.append(newinitial1_n_sum)
        newinitial2_n_sums.append(newinitial2_n_sum)
        #print(n,newinitial1_n_sum,newinitial2_n_sum,newinitial1_n_sum+ newinitial2_n_sum)
        #print (newinitial1_n_cooling,newinitial1_n_cooled)*
        
        newinitial_n_nocooling    = np.dot(n_matrix, newinitial_n_nocooling)
        newinitial_n_nocooling_sum = np.sum(newinitial_n_nocooling)
        newinitial_n_nocooling_sums.append(newinitial_n_nocooling_sum)
        
    return (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums,newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon)

#creat fname list
sqrt(-1)

measure1s = list()
measure2s = list()
measure3s = list()
measure4s = list()

Nte = 10000

wls = [535, 767,880, 1064]
for wli in range(0, len(wls)):
    wl = wls[wli]
    ifilename = "morematrices/"+str(wl)+"/intensities.csv"
    intlist =  np.loadtxt(ifilename,delimiter=",")
    ldfilename = "morematrices/"+str(wl)+"/latticedepths1.csv"
    ldlist =  np.loadtxt(ldfilename,delimiter=",")[:,0]
    measure1 = list()
    measure2 = list()
    measure3 = list()
    measure4 = list()
    for h in range(0, len(intlist)):
        cutoff = int(ldlist[h])+2  #max n, must matxh size of matrix
        initial_state = np.zeros(len(grounds))
        initial_state[6] = 1 #(4,0)
        initial_state_n = np.zeros(cutoff)
        initial_state_n[0] = 1 #All atoms in ground state
        n_matrix = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"Pathtot.csv",cutoff)
        (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums, newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon) = full_time_evolution(Nte, 0.1, 1,1,2,initial_state, initial_state_n, n_matrix, "sigma+", 1, "pi", -1)
        isel = 0    
        for i in range(0, len(newinitial1_n_sums)):
            left = newinitial1_n_sums[i] + newinitial2_n_sums[i]
            if (left > 0.95):
                isel = i
        print (h, len(intlist), wl, intlist[h], isel, newinitial_n_bluephotons_mon[isel])
        measure1.append(newinitial_n_bluephotons_mon[isel])
        
        initial_state = np.zeros(len(grounds))
        initial_state[6] = 1 #(4,0)
        initial_state_n = np.zeros(cutoff)
        initial_state_n[0] = 1 #All atoms in ground state
        n_matrix = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotE.csv",cutoff)
        (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums, newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon) = full_time_evolution(Nte, 0.1, 1,1,2,initial_state, initial_state_n, n_matrix, "sigma+", 1, "pi", -1)
        isel = 0    
        for i in range(0, len(newinitial1_n_sums)):
            left = newinitial1_n_sums[i] + newinitial2_n_sums[i]
            if (left > 0.95):
                isel = i
        print (intlist[h], isel, newinitial_n_bluephotons_mon[isel])
        measure2.append(newinitial_n_bluephotons_mon[isel])
        
        initial_state = np.zeros(len(grounds))
        initial_state[6] = 1 #(4,0)
        initial_state_n = np.zeros(cutoff)
        initial_state_n[0] = 1 #All atoms in ground state
        n_matrix = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotCDE.csv",cutoff)
        (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums, newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon) = full_time_evolution(Nte, 0.06, 1,1,2,initial_state, initial_state_n, n_matrix, "sigma+", 1, "pi", -1)
        isel = 0    
        for i in range(0, len(newinitial1_n_sums)):
            left = newinitial1_n_sums[i] + newinitial2_n_sums[i]
            if (left > 0.95):
                isel = i
        print (intlist[h], isel, newinitial_n_bluephotons_mon[isel])
        measure3.append(newinitial_n_bluephotons_mon[isel])
        
        initial_state = np.zeros(len(grounds))
        initial_state[6] = 1 #(4,0)
        initial_state_n = np.zeros(cutoff)
        initial_state_n[0] = 1 #All atoms in ground state
        n_matrix = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotR.csv",cutoff)
        (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums, newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon) = full_time_evolution(Nte, 0.06, 1,1,2,initial_state, initial_state_n, n_matrix, "sigma+", 1, "pi", -1)
        isel = 0    
        for i in range(0, len(newinitial1_n_sums)):
            left = newinitial1_n_sums[i] + newinitial2_n_sums[i]
            if (left > 0.95):
                isel = i
        print (intlist[h], isel, newinitial_n_bluephotons_mon[isel])
        measure4.append(newinitial_n_bluephotons_mon[isel]*4)
        
    measure1s.append(measure1)
    measure2s.append(measure2)
    measure3s.append(measure3)
    measure4s.append(measure4)

plt.figure(figsize=(10,8))
for wli in range(0, len(wls)):
    ifilename = "morematrices/"+str(wls[wli])+"/intensities.csv"
    intlist =  np.loadtxt(ifilename,delimiter=",")
    plt.subplot(2,2, wli+1)
    plt.loglog(intlist, measure1s[wli],marker = "o", color="blue", markerfacecolor="blue", markeredgecolor="navy", label="ABCDEF")
    plt.loglog(intlist, measure2s[wli],marker = "o", color="cornflowerblue", markerfacecolor="cornflowerblue", markeredgecolor="royalblue",label="ABCDF")
    plt.loglog(intlist, measure3s[wli],marker = "o", color="mediumpurple", markerfacecolor="mediumpurple", markeredgecolor="rebeccapurple",label="ABF")
    plt.loglog(intlist, 4*np.asarray(measure4s[wli]),marker = "o", color="tomato", markerfacecolor="tomato", markeredgecolor="red",label="R")
    ax = plt.gca()
    plt.text(0.75,0.9, str(wls[wli]) + " nm",transform = ax.transAxes)
    if wli == 0:
        plt.legend(frameon=False)
    plt.xlabel("$I\;(W/m^2)$")
    plt.ylabel("$N_{\gamma}$ (5% Loss)")
    if wli == 0:
        plt.xlim(10**5, 10**10)
    else:
        plt.xlim(10**5, 10**9)
    plt.ylim(0.1,250)

plt.tight_layout()
plt.show()



###################### MIXTURE RED BLUE ############################
measure1s = list()
measure2s = list()
measure3s = list()
measure4s = list()

wls = [535, 767,880, 1064]
c = 0
for wl in wls:
    ifilename = "morematrices/"+str(wl)+"/intensities.csv"
    intlist =  np.loadtxt(ifilename,delimiter=",")
    ldfilename = "morematrices/"+str(wl)+"/latticedepths1.csv"
    ldlist =  np.loadtxt(ldfilename,delimiter=",")[:,0]
    
    
    measure1 = list()
    measure2 = list()
    measure3 = list()
    measure4 = list()
    hs = [10,8,8,7]
    h = hs[c]
    c+=1    
    phis = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]
    for phi in phis:
        if phi >= 0.8:
            Ng = 40000
        elif phi >= 0.9:
            Ng = 100000
        else:
            Ng = 20000
        cutoff = int(ldlist[h])+2  #max n, must matxh size of matrix
        initial_state = np.zeros(len(grounds))
        initial_state[6] = 1 #(4,0)
        initial_state_n = np.zeros(cutoff)
        initial_state_n[0] = 1 #All atoms in ground state
        n_matrix1 = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"Pathtot.csv",cutoff)
        n_matrix2 = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotR.csv",cutoff)
        n_matrix = (1-phi)*n_matrix1 + phi*n_matrix2
        (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums, newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon) = full_time_evolution(Ng, 0.1, 1,1,1,initial_state, initial_state_n, n_matrix, "sigma+", 1, "pi", -1)
        isel = 0    
        for i in range(0, len(newinitial1_n_sums)):
            left = newinitial1_n_sums[i] + newinitial2_n_sums[i]
            if (left > 0.95):
                isel = i
        print (wl, phi, isel, newinitial_n_bluephotons_mon[isel])
        measure1.append(newinitial_n_bluephotons_mon[isel]*(1-phi))
        
        initial_state = np.zeros(len(grounds))
        initial_state[6] = 1 #(4,0)
        initial_state_n = np.zeros(cutoff)
        initial_state_n[0] = 1 #All atoms in ground state
        n_matrix1 = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotE.csv",cutoff)
        n_matrix2 = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotR.csv",cutoff)
        n_matrix = (1-phi)*n_matrix1 + phi*n_matrix2
        (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums, newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon) = full_time_evolution(Ng, 0.1, 1,1,1,initial_state, initial_state_n, n_matrix, "sigma+", 1, "pi", -1)
        isel = 0    
        for i in range(0, len(newinitial1_n_sums)):
            left = newinitial1_n_sums[i] + newinitial2_n_sums[i]
            if (left > 0.95):
                isel = i
        print (intlist[h], isel, newinitial_n_bluephotons_mon[isel])
        measure2.append(newinitial_n_bluephotons_mon[isel]*(1-phi))
        
        initial_state = np.zeros(len(grounds))
        initial_state[6] = 1 #(4,0)
        initial_state_n = np.zeros(cutoff)
        initial_state_n[0] = 1 #All atoms in ground state
        n_matrix1 = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotCDE.csv",cutoff)
        n_matrix2 = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotR.csv",cutoff)
        n_matrix = (1-phi)*n_matrix1 + phi*n_matrix2
        (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums, newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon) = full_time_evolution(Ng, 0.06, 1,1,1,initial_state, initial_state_n, n_matrix, "sigma+", 1, "pi", -1)
        isel = 0    
        for i in range(0, len(newinitial1_n_sums)):
            left = newinitial1_n_sums[i] + newinitial2_n_sums[i]
            if (left > 0.95): 
                isel = i
        print (intlist[h], isel, newinitial_n_bluephotons_mon[isel])
        measure3.append(newinitial_n_bluephotons_mon[isel]*(1-phi))
        
    measure1s.append(measure1)
    measure2s.append(measure2)
    measure3s.append(measure3)
    measure4s.append(measure4)
        
phis = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]
wls = [532, 767,880, 1064]
plt.figure(figsize=(10,8))
for wli in range(0, len(wls)):
    plt.subplot(2,2, wli+1)
    plt.semilogy(phis, measure1s[wli],marker = "o", color="blue", markerfacecolor="blue", markeredgecolor="navy", label="ABCDEF")
    plt.semilogy(phis, measure2s[wli],marker = "o", color="cornflowerblue", markerfacecolor="cornflowerblue", markeredgecolor="royalblue",label="ABCDF")
    plt.semilogy(phis, measure3s[wli],marker = "o", color="mediumpurple", markerfacecolor="mediumpurple", markeredgecolor="rebeccapurple",label="ABF")
    ax = plt.gca()
    plt.text(0.75,0.05, str(wls[wli]) + " nm",transform = ax.transAxes)
    if wli == 0:
        plt.legend(frameon=False)
    plt.xlabel("$\phi$")
    plt.ylabel("$N_{\gamma}$ (5% Loss)")
    plt.ylim(0.1,120)
    plt.xlim(0, 0.90)

plt.tight_layout()
plt.show()


#==== FIXED MICTURE RED BLUE; INDEPENDENCE OF INTENSITY ====
measure1s = list()
measure2s = list()
measure3s = list()
measure4s = list()

Nte = 100000
phi = 0.9
wls = [535, 767,880, 1064]
for wli in range(0, len(wls)):
    wl = wls[wli]
    ifilename = "morematrices/"+str(wl)+"/intensities.csv"
    intlist =  np.loadtxt(ifilename,delimiter=",")
    ldfilename = "morematrices/"+str(wl)+"/latticedepths1.csv"
    ldlist =  np.loadtxt(ldfilename,delimiter=",")[:,0]
    measure1 = list()
    measure2 = list()
    measure3 = list()
    measure4 = list()
    for h in range(0, len(intlist)):
        if h >= 8:
            Nte = 100000
        elif h == 7:
            Nte = 60000
        else:
            Nte = 20000    
        cutoff = int(ldlist[h])+2  #max n, must matxh size of matrix
        initial_state = np.zeros(len(grounds))
        initial_state[6] = 1 #(4,0)
        initial_state_n = np.zeros(cutoff)
        initial_state_n[0] = 1 #All atoms in ground state
        n_matrix1 = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"Pathtot.csv",cutoff)
        n_matrix2 = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotR.csv",cutoff)
        n_matrix = (1-phi)*n_matrix1 + phi*n_matrix2
        (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums, newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon) = full_time_evolution(Nte, 0.1, 1,1,2,initial_state, initial_state_n, n_matrix, "sigma+", 1, "pi", -1)
        isel = 0    
        for i in range(0, len(newinitial1_n_sums)):
            left = newinitial1_n_sums[i] + newinitial2_n_sums[i]
            if (left > 0.95):
                isel = i
        print (wl, intlist[h], isel, newinitial_n_bluephotons_mon[isel])
        measure1.append(newinitial_n_bluephotons_mon[isel]*(1-phi))
        
        initial_state = np.zeros(len(grounds))
        initial_state[6] = 1 #(4,0)
        initial_state_n = np.zeros(cutoff)
        initial_state_n[0] = 1 #All atoms in ground state
        n_matrix1 = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotE.csv",cutoff)
        n_matrix2 = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotR.csv",cutoff)
        n_matrix = (1-phi)*n_matrix1 + phi*n_matrix2
        (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums, newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon) = full_time_evolution(Nte, 0.1, 1,1,2,initial_state, initial_state_n, n_matrix, "sigma+", 1, "pi", -1)
        isel = 0    
        for i in range(0, len(newinitial1_n_sums)):
            left = newinitial1_n_sums[i] + newinitial2_n_sums[i]
            if (left > 0.95):
                isel = i
        print (intlist[h], isel, newinitial_n_bluephotons_mon[isel])
        measure2.append(newinitial_n_bluephotons_mon[isel]*(1-phi))
        
        initial_state = np.zeros(len(grounds))
        initial_state[6] = 1 #(4,0)
        initial_state_n = np.zeros(cutoff)
        initial_state_n[0] = 1 #All atoms in ground state
        n_matrix1 = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotCDE.csv",cutoff)
        n_matrix2 = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotR.csv",cutoff)
        n_matrix = (1-phi)*n_matrix1 + phi*n_matrix2
        (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums, newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon) = full_time_evolution(Nte, 0.06, 1,1,2,initial_state, initial_state_n, n_matrix, "sigma+", 1, "pi", -1)
        isel = 0    
        for i in range(0, len(newinitial1_n_sums)):
            left = newinitial1_n_sums[i] + newinitial2_n_sums[i]
            if (left > 0.95):
                isel = i
        print (intlist[h], isel, newinitial_n_bluephotons_mon[isel])
        measure3.append(newinitial_n_bluephotons_mon[isel]*(1-phi))
        
        initial_state = np.zeros(len(grounds))
        initial_state[6] = 1 #(4,0)
        initial_state_n = np.zeros(cutoff)
        initial_state_n[0] = 1 #All atoms in ground state
        n_matrix = load_matrix("morematrices/"+str(wl)+"/"+ str(h)+"PathtotR.csv",cutoff)
        (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums, newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon) = full_time_evolution(Nte, 0.06, 1,1,2,initial_state, initial_state_n, n_matrix, "sigma+", 1, "pi", -1)
        isel = 0    
        for i in range(0, len(newinitial1_n_sums)):
            left = newinitial1_n_sums[i] + newinitial2_n_sums[i]
            if (left > 0.95):
                isel = i
        print (intlist[h], isel, newinitial_n_bluephotons_mon[isel]*4)
        measure4.append(newinitial_n_bluephotons_mon[isel])
        
    measure1s.append(measure1)
    measure2s.append(measure2)
    measure3s.append(measure3)
    measure4s.append(measure4)

plt.figure(figsize=(10,8))
for wli in range(0, len(wls)):
    ifilename = "morematrices/"+str(wls[wli])+"/intensities.csv"
    intlist =  np.loadtxt(ifilename,delimiter=",")
    plt.subplot(2,2, wli+1)
    if len(measure1s) != 4:
        measure1s = measure1s[1:5]
    plt.loglog(intlist, measure1s[wli],marker = "o", color="blue", markerfacecolor="blue", markeredgecolor="navy", label="ABCDEF")
    plt.loglog(intlist, measure2s[wli],marker = "o", color="cornflowerblue", markerfacecolor="cornflowerblue", markeredgecolor="royalblue",label="ABCDF")
    plt.loglog(intlist, measure3s[wli],marker = "o", color="mediumpurple", markerfacecolor="mediumpurple", markeredgecolor="rebeccapurple",label="ABF")
    plt.loglog(intlist, 4*np.asarray(measure4s[wli]),marker = "o", color="tomato", markerfacecolor="tomato", markeredgecolor="red",label="R")
    ax = plt.gca()
    plt.text(0.75,0.05, str(wls[wli]) + " nm",transform = ax.transAxes)
    if wli == 0:
        plt.legend(frameon=False)
    plt.xlabel("$I\;(W/m^2)$")
    plt.ylabel("$N_{\gamma}$ (5% Loss)")
    plt.ylim(0.1,5000)
    if wli == 0:
        plt.xlim(10**5, 10**10)
    else:
        plt.xlim(10**5, 10**9)

plt.tight_layout()
plt.show()

        
