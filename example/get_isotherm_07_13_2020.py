# -*- coding: utf-8 -*-
"""
The purpose of this code is to generate the adsorption isotherm (and relative
thermodynamic quantities) of molecules adsorbed wihtin zeolites.

The input is a particle number probability distribution (PNPD) which comes from
the output of a grand canonical transition matrix Monte-Carlo simulation. In
particular, the free energy and advanced sampling toolkit (FEASST) is the software
used to run the simulation.
"""

import numpy as np
import re
import sys
import matplotlib.pyplot as plt
from scipy import constants

NA = constants.N_A
mass_cat = 7*9.58e-21/1000 #kg catalyst

'''
Suboutines and Functions
'''

def parse_colMat(fname):
    
    key_lnz = re.compile('lnz')
    key_beta = re.compile('beta')
    key_Nmax = re.compile('mMax')
    key_Nmin = re.compile('mMin')
    
    keys = [key_lnz, key_beta, key_Nmax, key_Nmin]
    
    main_list = []
    
    f = open(fname, 'r')
    for line in f:
        main_list.append(line.strip())
        
    f.close()
    
    var_list = ['None']*len(keys)
    body = []
    for i,line in enumerate(main_list):
        for j,key in enumerate(keys):
            if key.search(line):
                var_list[j] = float(line.split()[-1])
                
        if line[0] != '#':
            body.append(line.split())
    
            
    colMat_string = np.reshape(body,((int(var_list[-1])-int(var_list[-2])),len(body[0])))
    colMat_float = colMat_string.astype(np.float)
    
                
    return var_list, colMat_float
    
 
def reweight(lnz1,lnz2,N,lnPi_1,Beta):
    mu1 = lnz1 / Beta
    mu2 = lnz2 / Beta
    lnPi_RW = np.zeros(len(lnPi_1))
    #Reweight
    for i,lnPi in enumerate(lnPi_1):
        lnPi_RW[i] = (lnPi+Beta*N[i]*(mu2-mu1))
        
    #Shift new ln Pi to prevent overflow.
    max_lnPi = np.amax(lnPi_RW)
    lnPi_RW[i] -= max_lnPi
    
    #Calculate Pi* (temporary unnormalized Pi)
    Pi_RW = np.exp(lnPi_RW)
    
    #Normalize ln Pi_2
    sum_Pi_RW = np.sum(Pi_RW)
    lnPi_RW -= np.log(sum_Pi_RW)
    
    
        
    return lnPi_RW, lnz2
        
    
def get_pressure(V,Beta,lnPi_RW):
    #NOTE: lnPi must be already seperate based on phases.
    
    #Get Area
    Pi = 0
    for lnPi in lnPi_RW:
        Pi += np.exp(lnPi)
        
    return (-1*lnPi_RW[0]+np.log(Pi))/V/Beta/NA
        

def get_N_avg(N, lnPi):
    N_avg = 0
    for i,n in enumerate(N):
        N_avg += np.exp(lnPi[i])*n
        
    return N_avg

'''
Program Begins Here
'''

var_list_Bulk, colMat_Bulk = parse_colMat(r"TRAPPE_molecules\alkenes_T_300\ethylene\colMat")
var_list_Ads, colMat_Ads = parse_colMat(r"TRAPPE_molecules_ads\alkenes_T_300\ethylene\colMat")

cutoff_Bulk = 248
cutoff_Ads = 100
phase_boundary = 248

#plt.plot(colMat_Bulk[:cutoff_Bulk,0],colMat_Bulk[:cutoff_Bulk,1])
#plt.plot(colMat_Ads[:cutoff_Ads,0],colMat_Ads[:cutoff_Ads,1])

#Generate Isotherm
if True:
    V = 2.7e-26 # Volume of simulation cell in [m^3]
    lnz = -5.700769546847236 # chemical potential
    P_Bulk_list = []
    N_avg_list = []
    for i in np.linspace(-1,0.99685,num=200):
        lnz += i
        lnPi_RW_Bulk, lnz_Bulk = reweight(var_list_Bulk[0],lnz,colMat_Bulk[:cutoff_Bulk,0],colMat_Bulk[:cutoff_Bulk,1],var_list_Bulk[1])
        N_Bulk = get_N_avg(colMat_Bulk[:cutoff_Bulk,0],lnPi_RW_Bulk[:cutoff_Bulk])
        P_Bulk = get_pressure(V,var_list_Bulk[1],lnPi_RW_Bulk[:phase_boundary])/1000.00#Pressure in MPa, Volume give in units of [m^3]
        rho = N_Bulk/V #In units of [#particles/m^3], need MW to convert to mol.
        
        lnPi_RW_Ads, lnz_Ads = reweight(var_list_Ads[0],lnz,colMat_Ads[:cutoff_Ads,0],colMat_Ads[:cutoff_Ads,1],var_list_Ads[1])
        N_Ads = get_N_avg(colMat_Ads[:cutoff_Ads,0],lnPi_RW_Bulk[:cutoff_Ads])/NA/mass_cat
        
        P_Bulk_list.append(P_Bulk)
        N_avg_list.append(N_Ads)
        
if True:        
    plt.plot(P_Bulk_list,N_avg_list,linestyle='None',marker="o")
    plt.xlabel('Pressure [MPa]')
    plt.ylabel('Amount Adsorbed [mol/kg-cat]')




