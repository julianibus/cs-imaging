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

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

ICs = 7/2
q = 0

Natoms = 100

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
    plt.figure(figsize=(10,5))
    barlist = plt.bar(np.arange(0,len(strings),1), pil, align='center', alpha=0.7,edgecolor='b')
    for ibar in range(0,len(barlist)):
        try:
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

## EXPERIMENT 1: DECAY HISTOGRAMS
shares =list()
maxdecays = list()
for exi in exciteds:
    initial_state = np.zeros(len(exciteds))
    initial_state[exciteds_posdict[exi]] = 1
    out = np.dot(decay_matrix,initial_state) 
    #barplot(out,grounds_strings, "Decay from "  + str(exi), "Ground State (F, mF)", "Population")
    sum3 = sum(out[0:7])
    sum4 = sum(out[7:16])
    share3 = sum3/(sum3 +sum4)
    shares.append(share3)
    maxdecays.append(max(out))

print (np.asarray(shares), exciteds_strings)
barplot(np.asarray(shares), exciteds_strings, "", "Excited State (F, mF)", "Share of Decays into F=3")
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
        newinitial = np.dot(decay_matrix,np.dot(create_excitation_matrix(repump_mode, repump_deltaf),np.dot(raman_matrix(raman_mode,raman_deltaf),newinitial)))
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
        newinitial = np.dot(decay_matrix,np.dot(create_excitation_matrix(repump_mode, repump_deltaf),np.dot(raman_matrix(raman_mode,raman_deltaf),newinitial)))
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
    return tot_prob, cool_share

def experiment2():
    mode_labels = {"sigma+": "s+",
                   "sigma-": "s-",
                   "pi":"pi"}
    label_strings = list()
    eq_tot_probs = list()
    eq_cool_shares = list()
    for repump_mode in ["pi", "sigma+", "sigma-"]:
        for repump_deltaf in [0,+1]:
            for raman_mode in ["pi", "sigma+", "sigma-"]:
                for raman_deltaf in [-1,+1]:
                    label = "Re " + mode_labels[repump_mode] + " " +  str(repump_deltaf) + " Ra " + mode_labels[raman_mode] + " " + str(raman_deltaf)
                    print(label)
                    label_strings.append( label)
                    eq_tot_prob, eq_cool_share = timeevolution(60, initial_state, repump_mode, repump_deltaf, raman_mode, raman_deltaf)
                    eq_tot_probs.append(eq_tot_prob)
                    eq_cool_shares.append(eq_cool_share)
    barplot(eq_tot_probs, label_strings, "", "Configuration", "Prob")
    barplot(eq_cool_shares, label_strings, "", "Configuration", "Equilibrium Share of Atoms being Cooled")


## EXPERIMENT 3: TIME EVOLUTION INCLUDING TEMPERATURE
#Loading Matrix
def load_matrix(cutoff):
    totmatrix = np.loadtxt("/home/julian/Projects/cs-imaging/heatmap/PathtotCDE.csv",delimiter=",")[0:cutoff,0:cutoff]
    plt.matshow(totmatrix)
    rowsums = list()
    columnsums = list()
    for i in range(0,len(totmatrix)):
        rowsum = np.sum(totmatrix[i])
        columnsum = np.sum(totmatrix[:,i])
        #Stabilizer
        if (columnsum) > 1:
            totmatrix[:,i] = totmatrix[:,i]/np.sum(totmatrix[:,i])
            columnsum = np.sum(totmatrix[:,i])
        print (i, columnsum)  
    return totmatrix

cutoff = 25  #max n, must matxh size of matrix
initial_state = np.zeros(len(grounds))
initial_state[15] = 1 #(4,0)
initial_state_n = np.zeros(cutoff)
initial_state_n[7] = 1 #All atoms in ground state
n_matrix = load_matrix(cutoff) #matrix is loaded from heatmap folder -> first calculate totmatrix with heatpot_final.py

def full_time_evolution(N, off_res_ra, initial_state, initial_state_n, n_matrix, repump_mode, repump_deltaf, raman_mode, raman_deltaf):    
    
    tot_probs = list()
    cool_shares = list()
        
    newinitial = initial_state 
            
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
    for n in range(0,N):
        #plt.clf()
        #plt.ioff()
        #plt.subplot(211)
        #plt.bar(np.arange(0,len(newinitial),1), newinitial, align='center', alpha=0.7,edgecolor='b')    
        #plt.xticks(np.arange(0,len(grounds_strings),1), grounds_strings)
        #plt.xlabel("State")
        #plt.ylabel("Population")
        #plt.ylim([0,1])
        #plt.xticks(rotation=90)           
        
        #plt.subplot(212)
        #plt.bar(np.arange(0,len(newinitial_n),1), newinitial_n, align='center', alpha=0.7,edgecolor='b')
        
        #plt.xticks(np.arange(0,len(newinitial_n),1), np.arange(0,len(newinitial_n),1))
        #plt.xlabel("State")
        #plt.ylabel("Population")
        #plt.ylim([0,1])
        #plt.xticks(rotation=90)
        
        #plt.savefig("evol"+ '{0:05d}'.format(n) + ".png")
        
        #print ("=====",n)
        newinitial = np.dot(decay_matrix,np.dot(create_excitation_matrix(repump_mode, repump_deltaf),np.dot(raman_matrix(raman_mode,raman_deltaf),newinitial)))
        tot_prob = (np.sum(newinitial))
        sum3 = sum(newinitial[0:7])
        sum4 = sum(newinitial[7:16])
        if raman_deltaf == 1:
            cool_share = sum3/(sum3 +sum4)
        elif raman_deltaf == -1:
            cool_share = sum4/(sum3 +sum4)
        tot_probs.append(tot_prob)
        cool_shares.append(cool_share)
                
        
        #N VECTOR PART
        newinitial1_n_cooled = np.zeros(cutoff)
        for j in range(0, len(newinitial1_n)-1):
            newinitial1_n_cooled[j] = newinitial1_n[j+1] 
        newinitial1_n_cooled[len(newinitial1_n)-1] = 0
        
        newinitial1_n_cooling = np.zeros(cutoff)
        for j in range(1, len(newinitial1_n)):
            newinitial1_n_cooling[j] = newinitial1_n[j] 
        
        newinitial1_n = newinitial1_n + cool_share*np.dot(n_matrix, newinitial2_n) - newinitial1_n_cooling - off_res_ra*newinitial1_n + cool_share*off_res_ra*newinitial1_n 
        newinitial2_n = newinitial1_n_cooled + (1-cool_share)*np.dot(n_matrix, newinitial2_n)+ (1-cool_share)*off_res_ra*newinitial1_n 
        #newinitial_n = np.concatenate([[newinitial1_n],[newinitial2_n]])
        #newinitial2_n = newinitial2_n / np.sum(newinitial_n)
        #newinitial1_n = newinitial1_n / np.sum(newinitial_n)
        
        newinitial1_n_sum = np.sum(newinitial1_n)
        newinitial2_n_sum = np.sum(newinitial2_n)
        newinitial_n_bluephotons = newinitial_n_bluephotons + 0.25*(newinitial2_n_sum/(newinitial1_n_sum+newinitial2_n_sum))
        newinitial_n_bluephotons_mon.append(newinitial_n_bluephotons)
        newinitial_n_nocooling_bluephotons = newinitial_n_nocooling_bluephotons + 0.25
        newinitial_n_nocooling_bluephotons_mon.append(newinitial_n_nocooling_bluephotons)
        newinitial1_n_sums.append(newinitial1_n_sum)
        newinitial2_n_sums.append(newinitial2_n_sum)
        #print(n,newinitial1_n_sum,newinitial2_n_sum,newinitial1_n_sum+ newinitial2_n_sum)
        #print (newinitial1_n_cooling,newinitial1_n_cooled)
        
        newinitial_n_nocooling    = np.dot(n_matrix, newinitial_n_nocooling)
        newinitial_n_nocooling_sum = np.sum(newinitial_n_nocooling)
        newinitial_n_nocooling_sums.append(newinitial_n_nocooling_sum)
        
    return (tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums,newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon)

(tot_probs, cool_shares, newinitial1_n_sums,newinitial2_n_sums,newinitial_n_nocooling_sums, newinitial_n_bluephotons_mon,newinitial_n_nocooling_bluephotons_mon) = full_time_evolution(100, 0, initial_state, initial_state_n, n_matrix, "sigma+", 1, "pi", -1)
plt.clf()
plt.plot(np.arange(0, len(newinitial2_n_sums)), np.asarray(newinitial1_n_sums), "blue")
plt.plot(np.arange(0, len(newinitial2_n_sums)), np.asarray(newinitial2_n_sums), "green")
plt.xlabel("Cycles")
plt.ylabel("Population Probability")
plt.grid(which="both")
plt.xlim((1, len(newinitial2_n_sums)))
plt.ylim((0,1))
plt.savefig("prop-fstates.png")
plt.show()

plt.clf()
plt.plot(np.arange(0, len(newinitial2_n_sums)), np.asarray(newinitial2_n_sums) + np.asarray(newinitial1_n_sums), "black")
plt.plot(np.arange(0, len(newinitial2_n_sums)), np.asarray(newinitial_n_nocooling_sums), "red")
plt.xlabel("Cycles")
plt.ylabel("Share of Remaining Atoms")
plt.grid(which="both")
plt.xlim((1, len(newinitial2_n_sums)))
plt.ylim((0,1))
plt.savefig("loss-cycles.png")
plt.show()

plt.clf()
plt.plot(newinitial_n_bluephotons_mon, np.asarray(newinitial2_n_sums) + np.asarray(newinitial1_n_sums), "black")
plt.plot(newinitial_n_nocooling_bluephotons_mon, np.asarray(newinitial_n_nocooling_sums), "red")
plt.xlabel("Number of Blue Photons")
plt.ylabel("Share of Remaining Atoms")
plt.grid(which="both")
plt.xlim((1, 25))
plt.ylim((0,1))
plt.savefig("loss-bluephotons.png")
plt.show()