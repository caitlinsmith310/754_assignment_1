#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 23:25:42 2020

@author: caitlin
"""


import numpy as np
import matplotlib.pyplot as plt
import random as ran
from pylab import rcParams
from numba import jit
import timeit
start = timeit.default_timer()

iSeed = None#'iSeed'
ran.seed(iSeed)

method="Spanning"


# ------------------------ INPUT -------------------------------------------
if method == "Spanning":
    b=2
if method == "Majority":
    b=3    

NLength = b**7       # Lattice size (NLength x NLength)
prob    = 0.282       # occupation probability
Nren    = 5         # number of renormalization steps

print("original length of square lattice", NLength)

# ---------------- function for renormalization, using majority rule --------
#NOTE: b should be odd number
@jit(nopython=True) #compile to machine code 
def LatticeRenMaj(L,b,lattice):
    #print("next step")
    #print(len(lattice))
    Lnew = int(L/b)
    newlattice = np.empty((Lnew, Lnew))
    maj = int(b*b/2) + 1
    for j in range(Lnew):
        for i in range(Lnew):
            #majority rule
            sum = 0.
            for k in range(b):
                for l in range(b):
                    sum += lattice[b*i+k,b*j+l]
            occ = int(sum/maj)
            newlattice[i,j] = occ
    return newlattice, Lnew

# ---------------- function for renormalization, using top-down spanning rule --------
#renormalization, spanning cluster, top-bottom (do it for b=2 only)
@jit(nopython=True) #compile to machine code 
def LatticeRenSpan(L,b,lattice):
    Lnew = int(L/b)
    newlattice = np.empty((Lnew, Lnew))
    for j in range(Lnew):
        for i in range(Lnew):
            #top down
            sum=0
            column = np.zeros(b)
            row = np.zeros(b)
            for l in range (b): # for each column
                for k in range (b):  #for each row
                    column[l] = lattice[b*i+0,b*j+l]*lattice[b*i+1,b*j+l]  #product of top and bottom of column
                    row[k] = lattice[b*i+k,b*j+0]*lattice[b*i+k,b*j+1]    #product of left and right of row                                    #sum of product of 1st column vs second 
            sum= row[0]+row[1]+column[0]+column[1]                        #if spanning= shoudld be >1
            if sum > 0.5:
                occ=1
            else:
                occ=0
            newlattice[i,j] = occ
    #print(newlattice)
    return newlattice, Lnew


# -------- Initialization of lattice ----------------

Lattice = np.empty((NLength, NLength))

for j in range(NLength):
    for i in range(NLength):
        occ = 1
        r = ran.random() 
        if (r > prob): 
            occ = 0  
        Lattice[i][j] = occ
 

#upper left window of lattice that is shown (size of lattice in last renormalization step)
#Lwindow = int(NLength/b**Nren) 
Lwindow= NLength
#-------- renormalization loop ----------------

for i in range(Nren):
    print(NLength)
    plt.figure(1)
    plt.matshow(Lattice[0:Lwindow,0:Lwindow],vmin=0, vmax=1)
    #plt.title("NLength",i)
    #Latticenew, Lnew = LatticeRenMaj(NLength,b,Lattice) 
    Latticenew, Lnew = LatticeRenSpan(NLength,b,Lattice) 

    Lattice = Latticenew
    NLength = Lnew
    
stop = timeit.default_timer()
print('Execution time:', stop-start, 'seconds')


#R_b(p) for the majority rule, b=3
def rgt_maj(p):
    sum = 0 
    for i in range(5):
        sum += p**(9-i)*(1-p)**i*math.factorial(9)/math.factorial(i)/math.factorial(9-i)
    return sum

def rgt_span_tb(p):             #top bottom spanning
    return 2*p**2 - p**4

def rgt_span_lr(p):             #left-right spanning
    return 2*p**2-p**4

def rgt_span_all(p):            #left-right and top-bottom spanning
    result=(p**4)+(4*((p**3)*(1-p)))+4*(p**2)*((1-p)**2)
    return result


if method=="Majority":
    pvals = np.arange(0,1,0.01)
    plt.figure(10)
    plt.plot(pvals,rgt_maj(pvals))
    plt.plot(pvals,pvals)
    plt.title("RGT for square lattice, b=3, majority rule")
    plt.xlabel("occupation probability $p$")
    plt.ylabel("$R_b(p)$")
    plt.grid(True)
     
if method=="Spanning":  

    pvals = np.arange(0,1,0.001)
    
    plt.figure(10, [8,6])
    plt.plot(pvals,rgt_span_tb(pvals), label="Top-bottom / left-right")
    plt.plot(pvals,rgt_span_all(pvals), label="All spanning rules")
    plt.plot([0.382,0.382],[0,1], linestyle="dotted", color="orangered")
    plt.plot([0.618,0.618],[0,1], linestyle="dotted", color="#1f77b4")
    rcParams['figure.figsize'] = 7, 7
    plt.plot(pvals,pvals, label="$R_b(p)$=p")
    plt.title("RGT for square lattice, b=2, different spanning rules")
    plt.xlabel("Occupation probability $p$")
    plt.ylabel("$R_b(p)$")
    plt.legend()
    plt.grid(True)
    plt.gca().set_aspect('equal', adjustable='box')
    
    
    difference1=abs(pvals-rgt_span_tb(pvals))
    index1=np.where(difference1 <= 0.001)
    index1=index1[0]
    print("ALL:Indices where RGT is very close to pvals at :",index1)
    print("non-trivial fixed point at p=",pvals[index1[3]])
    
    difference2=abs(pvals-rgt_span_all(pvals))
    index2=np.where(difference2 <= 0.001)
    index2=index2[0]
    print("TB only: indices where RGT is very close to pvals at :",index2)
    print("non-trivial fixed point at p=",pvals[index2[3]])
    
    plt.figure(11)
    plt.plot(pvals, difference1, label='Top-bottom/ Left right')
    plt.plot(pvals, difference2, label="All spanning")
    plt.xlabel("Occupation probability $p$")
    plt.ylabel("Absolute difference between P and RGT")
    plt.grid(True)
    plt.title("Difference between RGT and P for top-bottom and all spanning")
    plt.legend()


