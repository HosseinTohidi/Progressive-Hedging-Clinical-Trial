#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 19:35:47 2019

@author: hosseintohidi
"""
import matplotlib.pyplot as plt
import numpy as np
import sys
from gurobipy import *
import numpy as np
import pandas as pd
from tabulate import tabulate
import matplotlib.pyplot as plt
from pprint import pprint
import time
import progressbar
from myPackage.myClusters import myClusters
from myPackage.mySimulation import mySimulation
from myPackage.myGplot import myGplot
from myPackage.solve import solve, solveNew, averageCalc, is_inList, solve_for_bound
import networkx as nx
from joblib import Parallel, delayed
from multiprocessing import cpu_count
parallel = Parallel(n_jobs = cpu_count())

# Parameters Definitions 
np.random.seed(123)
M= 1000 #bigM
eps=0.00001
TL=1000 #model time limit

ip=5 #of potential facilities (cluster of facilities)
j= 2 # of drugs to be tested
t=10  #planning horizon
r=1 # of resources
omega=10 #of scenarios


group_size = {5:2, 10:3, 20:5, 100:10}
max_omega = 100
I= np.arange(ip)
J= np.arange(j)
T= np.arange(t)
R= np.arange(r)
Omega= np.arange(omega)
delta = 2 #time from deciding to open till runing clinical trial
Np     =np.random.randint(2,10,size=[5])*100  # number of patients required for each trial if goes to the ned
rho   =np.random.randint(3,size=[5,2]) #resources required for each test
rhoMax=np.random.randint(50,100,size=[10,2]) #max available resource in each center at each time
w     =np.random.random([10,5,30,max_omega])/4 # stochastic term representing the proportion of people who actually arrived
w += 0.75
w = w.round(2)

#w     = np.ones([ip,j,t,omega])
p=np.ones(omega)*1/omega
fp    = np.random.randint(1,100,size=[10,5])*1  #fix cost for oppening centers in 1000 dollar
fpp   = np.random.randint(1,100,size=[10,5,30])*1 #fix cost for opening centers after time zero in 1000 dollar
#changing fpp
for ii in range(10):
    for jj in range(5):
        for tt in range(30):
            fpp[ii,jj,tt] = fp[ii,jj] * 1.1**tt    
g     = np.random.randint(10,100,size=[10,5])*0.1 #maintenance cost for keeping once center open in 1000 dollars
v     = np.random.randint(1,100,size=[10,5])*1000 # capacity of each center 1000
kp    = np.int64(np.random.randint(100,1000,size=[5])*10) # drug profit in 1000 dollars
alpha = np.random.random(size=5)/2
#rounding random numbers to two digit
alpha = alpha.round(2)
alpha = np.zeros(len(J))
l = np.random.randint(3,5,size=j) # study time
Tj=[(j,t) for j in J for t in T[l[j]:]]
rho = rho[:ip,:r]
rhoMax = rhoMax[:ip,:r]
w = w[:ip,:j,:t,:omega]
fp = fp[:ip,:j]
fpp = fpp[:ip,:j,:t]
for tp in T[1:]:
    fpp[:ip,:j,tp] = fpp[:ip,:j,tp-1] *1.4
g = g[:ip,:j]
v = v[:ip,:j]
alpha = alpha[:j]
drugs_actual_p = [0.6,0.65,0.7,0.6,0.5]
N = np.zeros([5,max_omega]) 
k = np.zeros([5,max_omega])
np.random.seed(100)

# simulation starts
print('Simulation Started...')
for drug in range(5):
    f1,f2,f3,s = mySimulation(drugs_actual_p[drug], plot = False)
    cumsum = np.array([f1,f2,f3,s]).cumsum()
    for scenario in range(max_omega):
        rnd = np.random.random()
        for ii in range(4): # 4 possible outcome
            if rnd <= cumsum[0]:
                N[drug, scenario] = Np[drug] * 0.1
                k[drug, scenario] = 0
            elif rnd <= cumsum[1]:
                N[drug, scenario] = Np[drug] * 0.5
                k[drug, scenario] = 0
            elif rnd <= cumsum[2]:
                N[drug, scenario] = Np[drug]
                k[drug, scenario] = 0            
            else:
                N[drug, scenario] = Np[drug]
                k[drug, scenario] = np.round(np.random.triangular(kp[drug]*0.5, kp[drug], kp[drug]*1.1))
print('Simulation ended.')

drugs_actual_p = drugs_actual_p[:j]

kp = kp[:j]
Np = Np[:j]
        
N = N[:j,:omega]
k = k[:j,:omega]        

#not_testing_cost = kp*240 #12 month 20 years worth of patent in 1000 dollars    
not_testing_cost = kp*30 #12 month 20 years worth of patent in 1000 dollars    

####################################   
### Progressive hedging#############
#################################### 
def adjusting_rho(solList):
    rho_dict = {}
    for var in solList[0].keys():
        x_l = [solList[ii][var] for ii in range(len(solList))]
        if 'M[' in var:
            varL= var.split(',')
            i = int(varL[0][2:])
            j = int(varL[1][:-1])
            rho_dict[var] = fp[i,j]/ ( max(x_l) - min(x_l) + 1)
        else:
            j = int(var[2:-1])
            rho_dict[var] = kp[j] / (max(x_l) - min(x_l) + 1)
    return rho_dict


def PH(group_type, rhoCoef, rho_adjustment, rho_dict= {}, plot = False):
    cluster = []
    print('here0')
    if group_type == 'Extensive Form':
        print('here1')

        omega1 = Omega
        first_sol, sol = solve(omega1, I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL, extra_print=True)
        counter = 0
        v_solList = sol
    else:
        print('here2')
        if group_type == 'No Grouping': 
            clNum = -1             
            solList = []
            counter = 0
            # step 2: solve for each scenario
            #print('Proressive Hedging started ...')
            #print('iteration 0')
            #for scenario in Omega:
            #    first_sol, sol = solve(scenario, I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL)
            #    solList.append(first_sol)
            # make the previous lines parallel
            solList2 = parallel(delayed(solve)(scenario, I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL) for scenario in Omega)
            for ij in range(len(solList2)):
                solList.append(solList2[ij][0])
                
            #step 3: average scenario
            averageDict = averageCalc(p,solList)
            
            #step 4
            wDict_s = {}
            if not rho_adjustment:
                for omega in Omega:
                    for var in averageDict.keys():        
                        wDict_s[(var,omega)] = round(rhoCoef* (solList[omega][var]- averageDict[var]),5)
            else:
                rho_dict = adjusting_rho(solList)
                for omega in Omega:
                    for var in averageDict.keys(): 
                        wDict_s[(var,omega)] = round(rho_dict[var] * (solList[omega][var]- averageDict[var]),5)
    
            check = True
            #random numbers for cycle detection
            random_arr1 = np.random.random(size = [len(I), len(J),len(Omega)])
            random_arr2 = np.random.random(size = [len(J),len(Omega)])
            
            hDict = {}
            counter = 0
            fixed_var = {}

            while check and counter <=99:
                # Find LB only fix variables that have been fixed so far in PH which is {} at this point
                #sol_LB = parallel(delayed(solve_for_bound)(omega, I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL, fixed_var) for omega in Omega)
                #obj_LB = sum(sol_LB)
                # find UB
                #chose the first sol and force it to all scenario
                #sol_UB = parallel(delayed(solve_for_bound)(omega, I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL, solList[0]) for omega in Omega)
                #obj_UB = sum(sol_UB)
                
                #print(f'LB = {obj_LB}, UB = {obj_UB}, GAP: {obj_UB - obj_LB} $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
                print(counter)
                
                #step 5
                counter +=1 
                v_solList = []
                #step 6:
                newSolList = []
                #for omega in Omega:
                #    ml, sol, objVal = solveNew(clNum, omega, rhoCoef, wDict_s, averageDict, fixed_var,
                #                               I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL,rho_dict)
                #    newSolList.append(ml)
                #    v_solList.append(sol)
                
                
                
                solList2 = parallel(delayed(solveNew)(clNum, omega, rhoCoef, wDict_s, averageDict, fixed_var,
                                               I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL,rho_dict) for omega in Omega)# Omega)
                for ij in range(len(solList2)):
                    newSolList.append(solList2[ij][0])
                    v_solList.append(solList2[ij][1])
                
                #step 7
                new_averageDict = averageCalc(p, newSolList) 
                
                #step 8  
                new_wDict_s = {}
                if not rho_adjustment:
                    for omega in Omega:
                        for var in averageDict.keys():        
                            new_wDict_s[(var,omega)] = round(wDict_s[(var,omega)] + rhoCoef* (newSolList[omega][var]- new_averageDict[var]),5)
                else:
                    rho_dict = adjusting_rho(newSolList)
                    for omega in Omega:
                        for var in averageDict.keys(): 
                           # print(var)
                            new_wDict_s[(var,omega)] = round(wDict_s[(var,omega)] + rho_dict[var]* (newSolList[omega][var]- new_averageDict[var]),5)
                            
                #cycle detection    
                for var in newSolList[0].keys():
                    if var not in fixed_var.keys():
                        if 'M' in var:
                            varL= var.split(',')
                            i = int(varL[0][2:])
                            j = int(varL[1][:-1])
                            hDict[(counter,var)] = sum([random_arr1[i,j,omega] * wDict_s[(var,omega)]for omega in Omega])
                            
                        elif 'Z' in var:
                            j = int(var[2:-1])
                            hDict[(counter,var)] = sum([random_arr2[j,omega] * wDict_s[(var,omega)]for omega in Omega])
                    
                        # check the convergence of variables
                        if hDict[(counter,var)] in [hDict[(ctr,var)] for ctr in range(1,counter)]:
                            print('warning')
                            # fix that specific var for the next iterations
                            fixed_var[var] = max([newSolList[omega][var] for omega in Omega])     
                            print('**** ', var, fixed_var[var])       
                #step 9 
                gSum = 0
                for omega in Omega:
                    sumP = 0
                    for j in J:
                        sumP += (newSolList[omega]['Z['+str(j)+']'] - new_averageDict['Z['+str(j)+']'])**2
            
                        for i in I:
                            sumP += (newSolList[omega]['M['+str(i)+','+str(j)+']'] - new_averageDict['M['+str(i)+','+str(j)+']'])**2
                    gSum += p[omega]*np.sqrt(sumP)
                
                if gSum < 0.1: 
                    check = False
                wDict_s = new_wDict_s
                averageDict = new_averageDict
                print(f'iteration {counter} , error = {gSum}')
                
            #print(f'Objective function = {objVal}')    
                
        else:   # if need clustering
            # Grouping Scenarios (we have three stochastic parameters.)
            flat_scenarios = []
            for omega in Omega:
                one_scenario = []
                one_scenario.extend(N.transpose()[omega])
                one_scenario.extend(k.transpose()[omega])
                one_scenario.extend(w[:,:,:,omega-1].flatten())
                flat_scenarios.append(one_scenario)
            cluster = myClusters(np.array(flat_scenarios)[:,:20],group_size[len(Omega)], group_type)
            # change the probability of each cluster to sum of the scenario in that cluster:
            prob0 = np.array([sum(p[cluster[ii]]) for ii in sorted(cluster.keys())])
            prob = prob0/sum(prob0)
                
            
            solList = []
            counter = 0
            # step 2: solve for each scenario
            #print('Proressive Hedging started ...')
            #print('iteration 0')
            #for ii in sorted(cluster.keys()):
            #    scenario_list = cluster[ii]
            #    first_sol,sol = solve(scenario_list, I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL)
            #    solList.append(first_sol)
            solList2 = parallel(delayed(solve)(cluster[ii], I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL) for ii in sorted(cluster.keys()))    
            for ij in range(len(solList2)):
                solList.append(solList2[ij][0])
            #step 3: average scenario
            averageDict = averageCalc(prob,solList)
            
            #step 4
            #wList = [[] for omega in Omega] #for step 4
            wDict_s = {}
            if not rho_adjustment:
                for ii in range(len(solList)):
                    for var in averageDict.keys():
                        wDict_s[(var,ii)] = round(rhoCoef* (solList[ii][var]- averageDict[var]),5)
            else:
                rho_dict = adjusting_rho(solList)
                for ii in range(len(solList)):
                    for var in averageDict.keys():
                        wDict_s[(var,ii)] = round(rho_dict[var]* (solList[ii][var]- averageDict[var]),5)
    
            check = True
            #random numbers for cycle detection
            random_arr1 = np.random.random(size = [len(I), len(J),len(solList)])
            random_arr2 = np.random.random(size = [len(J),len(solList)])
            
            hDict = {}
            counter = 0
            fixed_var = {}
            while check and counter <=99:
                #step 5
                counter +=1 
                v_solList = []
                #step 6:
                newSolList = []
                #for ii in sorted(cluster.keys()):
                #    scenario_list = cluster[ii]
                #    ml, sol, objVal = solveNew(ii, scenario_list, rhoCoef, wDict_s, averageDict, fixed_var,
                #                               I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL,rho_dict)
                #    newSolList.append(ml)
                #    v_solList.append(sol)
                
                solList2 = parallel(delayed(solveNew)(ii, cluster[ii], rhoCoef, wDict_s, averageDict, fixed_var,
                                               I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL,rho_dict) for ii in sorted(cluster.keys()))
                for ij in range(len(solList2)):
                    newSolList.append(solList2[ij][0])
                    v_solList.append(solList2[ij][1])
                
                #step 7
                new_averageDict = averageCalc(prob, newSolList) 
                
                #step 8  
                new_wDict_s = {}
                if not rho_adjustment:
                    for ii in range(len(newSolList)):
                        for var in averageDict.keys():
                            new_wDict_s[(var,ii)] = round(wDict_s[(var,ii)] + rhoCoef* (newSolList[ii][var]- new_averageDict[var]),5)
                else:
                    rho_dict = adjusting_rho(newSolList)
                    for ii in range(len(solList)):
                        for var in averageDict.keys():
                            new_wDict_s[(var,ii)] = round(wDict_s[(var,ii)] +  rho_dict[var]* (newSolList[ii][var]- new_averageDict[var]),5)
            
                #cycle detection    
                for var in newSolList[0].keys():
                    if var not in fixed_var.keys():
                        if 'M' in var:
                            varL= var.split(',')
                            i = int(varL[0][2:])
                            j = int(varL[1][:-1])
                            hDict[(counter,var)] = sum([random_arr1[i,j,omega] * wDict_s[(var,omega)]for omega in range(len(solList))])
                            
                        elif 'Z' in var:
                            j = int(var[2:-1])
                            hDict[(counter,var)] = sum([random_arr2[j,omega] * wDict_s[(var,omega)]for omega in range(len(solList))])
                    
                        # check the convergence of variables
                        if hDict[(counter,var)] in [hDict[(ctr,var)] for ctr in range(1,counter)]:
                            print('warning')
                            # fix that specific var for the next iterations
                            fixed_var[var] = max([newSolList[omega][var] for omega in range(len(solList))])     
                            print('**** ', var, fixed_var[var] )
            
                #step 9 
                gSum = 0
                for omega in range(len(solList)):
                    sumP = 0
                    for j in J:
                        sumP += (newSolList[omega]['Z['+str(j)+']'] - new_averageDict['Z['+str(j)+']'])**2
            
                        for i in I:
                            sumP += (newSolList[omega]['M['+str(i)+','+str(j)+']'] - new_averageDict['M['+str(i)+','+str(j)+']'])**2
                    gSum += p[omega]*np.sqrt(sumP)
                
                if gSum < 0.1: 
                    check = False
                wDict_s = new_wDict_s
                averageDict = new_averageDict
                print(f'iteration {counter} , error = {gSum}')
                
        #print(f'Objective function = {objVal}')               
        # visualization:  
    if plot:  #default: False
        for omega in range(len(solList)):
            plt.subplot(len(solList),1,omega+1) 
            G = nx.Graph()
            sol = v_solList[omega]
            for i in range(len(I)*len(J)):
                for t in T:
                    G.add_node(len(T) * i+ t, name = len(T) * i+ t, pos = [2+t, i], nColor = 'Red')           
            for j in J:
                G.add_node('J'+str(j), name = 'J'+str(j), pos = [1, len(I)*j+len(I)/2], nColor = 'Red')
            G.add_node('s', name = 's', pos = [0, len(I)*len(J)//2], nColor = 'Red')
            
            for j in J:
                G.add_edge('s','J'+str(j), edge_color = 'black')   
                for i in I:
                    G.add_edge('J'+str(j), len(T)*(i+len(I)*j), edge_color = 'black')
                    for t in T[:-1]:
                        G.add_edge(len(T)*(i+len(I)*j) +t , len(T)*(i+len(I)*j) +(t+1), edge_color = 'black')
               
            for var in sol.keys():
                if 'Z' in var:
                    j = var[2:-1]
                    G.nodes(data = True)['J'+j]['nColor'] = 'Blue'
                    
                if 'O[' in var:
                    varL = var.split(',')
                    i = int(varL[0][2:])
                    j = int(varL[1])
                    t = int(varL[2])
                    omega = int(varL[3][:-1])
                    G.nodes(data = True)[len(T)*(i+len(I)*j) +t]['nColor'] = 'Green'
                
            myGplot(G, True, False)
        print('warning')
    return cluster, counter, v_solList
        

def extensive_form__obj_calc(sol, Omega):    
    obj1 = sum([fp[i,j] * sol.get('M['+str(i)+','+str(j)+']',0) for i in I for j in J])
    obj6 = sum([(1-sol.get('Z['+str(j)+']',0))*not_testing_cost[j] for j in J])    
    obj2 = sum([p[omega]*fpp[i,j,t]*sol.get('Mp['+str(i)+','+str(j)+','+str(t)+','+str(omega)+']',0) for i in I for j in J for t in T for omega in Omega])    
    obj3 = sum([p[omega]*g[i,j]*sol.get('O['+str(i)+','+str(j)+','+str(t)+','+str(omega)+']',0) for i in I for j in J for t in T for omega in Omega])    
    obj4 = sum([p[omega]*k[j,omega]*(sol.get('D['+str(j)+','+ str(omega)+']',0)- (T[-1]+1)*(1-sol.get('Z['+str(j)+']',0))) for j in J for omega in Omega])    
    obj5 = sum([p[omega]*eps*(sol.get('INV['+str(i)+','+str(j)+','+str(t)+','+str(omega)+']',0)+sol.get('U['+str(i)+','+str(j)+','+str(t)+','+str(omega)+']',0)) for i in I for j in J for t in T for omega in Omega])    
    return obj1+obj2+obj3+obj4+obj5+obj6
    
def obj_recalc(solList, cluster, p):
    Omega = range(len(solList))
    obj1 = sum([fp[i,j] * solList[0].get('M['+str(i)+','+str(j)+']',0) for i in I for j in J])
    obj6 = sum([(1-solList[0].get('Z['+str(j)+']',0))*not_testing_cost[j] for j in J]) 
    obj2 = 0
    obj3 = 0
    obj4 = 0
    obj5 = 0
    if cluster == []:
      cluster = {}
      for omega in Omega:
          cluster[omega] = [omega]   
    for omega in Omega:
        obj2 += sum([p[gamma]*fpp[i,j,t]*solList[omega].get('Mp['+str(i)+','+str(j)+','+str(t)+','+str(gamma)+']',0) for i in I for j in J for t in T for gamma in cluster[omega]])    
        obj3 += sum([p[gamma]*g[i,j]*solList[omega].get('O['+str(i)+','+str(j)+','+str(t)+','+str(gamma)+']',0) for i in I for j in J for t in T for gamma in cluster[omega]])    
        obj4 += sum([p[gamma]*k[j,gamma]*(solList[omega].get('D['+str(j)+','+ str(gamma)+']',0)- (T[-1]+1)*(1-solList[omega].get('Z['+str(j)+']',0))) for j in J for gamma in cluster[omega]])    
        obj5 += sum([p[gamma]*eps*(solList[omega].get('INV['+str(i)+','+str(j)+','+str(t)+','+str(gamma)+']',0)+solList[omega].get('U['+str(i)+','+str(j)+','+str(t)+','+str(gamma)+']',0)) for i in I for j in J for t in T for gamma in cluster[omega]])    
    
    return obj1+obj2+obj3+obj4+obj5+obj6

Grouping_options = ['No Grouping','random', 'min', 'max']
rho_adjustment_options = [True, False]
sol_dict = {}
iteration = 0

def output_finder(solList, cluster, p):
    opened_sites_0, opened_sites_1 = 0, 0
    tested_drugs = 0
    if type(solList)==dict: #answer from extensive form         
        for var in solList.keys():
            if 'M[' in var:
                opened_sites_0 += solList[var]
            elif 'Mp[' in var:
                opened_sites_1 += solList[var]/len(p)
            elif 'Z[' in var: 
                 tested_drugs += solList[var]
     
    else:
        if cluster != []:
            for ii in range(len(solList)):
                for jj in solList[ii].keys():
                    if 'M[' in jj:
                        print(jj,solList[ii][jj])
                        opened_sites_0 += len(cluster[ii])*solList[ii][jj]/len(p)
                    elif 'Mp[' in jj:
                        opened_sites_1 += solList[ii][jj]/len(p)
                    elif 'Z[' in jj: 
                        tested_drugs   += len(cluster[ii])*solList[ii][jj]/len(p)
        else:
            for ii in range(len(solList)):
                for jj in solList[ii].keys():
                    if 'M[' in jj:
                        print(jj,solList[ii][jj])
                        opened_sites_0 += solList[ii][jj]/len(p)
                    elif 'Mp[' in jj:
                        print(jj,solList[ii][jj])
                        opened_sites_1 += solList[ii][jj]/len(p)
                    elif 'Z[' in jj: 
                        print(jj,solList[ii][jj])
                        tested_drugs   += solList[ii][jj]/len(p)
                
    return round(opened_sites_0), round(opened_sites_1), round(tested_drugs)
                    
    
print('**************************************************************')
for rhoCoef in [500]:
    for g_option in ['random']: #['No Grouping','min','max','random'] : #Grouping_options:
        for rho_option in [True]: #rho_adjustment_options:
            print(iteration, rhoCoef,g_option,rho_option)
            iteration += 1
            time0 = time.time()
            rho_dict = {}
            try:
                cluster, counter, solList = PH(g_option, rhoCoef, rho_option, rho_dict,plot = True) # if g_options = 'Extensive Form' the solList will be a dictionary not list
            except:
                cluster = []
                counter = 1000 #max value to say it is not a feasible solution
                solList = []                    
            time1 = time.time()
            if counter < 100:
                obj_val = obj_recalc(solList, cluster,p)
                os0,os1,td = output_finder(solList, cluster, p)   #return total number of sites opened in the first stage, average number of sites opened in the second stage and total number of tested drugs
            else:
                obj_val = np.NaN
                os0,os1,td = np.NaN, np.NaN, np.NaN   
            sol_dict[(g_option, rho_option, rhoCoef)] = (td, os0, os1, obj_val, time1-time0, counter)
            print(f'{time1-time0, counter, obj_val}')
            print('**************************************************************')

# run for extensive form only once
iteration += 1
g_option = 'Extensive Form'
time0 = time.time()
cluster, counter, solList = PH(g_option, rhoCoef, rho_option, rho_dict) # if g_options = 'Extensive Form' the solList will be a dictionary not list
os0, os1, td = output_finder(solList, cluster, p)
time1 = time.time()
obj_val = extensive_form__obj_calc(solList, Omega)
sol_dict[(g_option, rho_option, rhoCoef)] = (td, os0, os1, obj_val, time1-time0, counter)


#store in df to be saved in files
df = pd.DataFrame(sol_dict).T        
df.reset_index(inplace = True)
df.columns = ['grouping method','rho_adjustment','rho_ceof','drugs_tested','first_stage_openSite','second_stage_openSite','objective function','elapsed time', 'No. iterations']
df.to_csv('new_run1.csv')


        

#rho_dict = {}          
#rhoCoef = 1000   
#g_option = 'No Grouping'    
#rho_option = False   
#cluster = []  
#cluster, counter, solList = PH(g_option, rhoCoef, rho_option, rho_dict) # if g_options = 'Extensive Form' the solList will be a dictionary not list
#if g_option == 'Extensive Form':
#    obj_val = extensive_form__obj_calc(solList, Omega)
#else:
#    obj_val = obj_recalc(solList, cluster,p)
#
#
#
#
#for i in range(len(solList)):
#    if i ==0:
#        df_f = pd.DataFrame(solList[i], index=['omega'+str(i)]).T        
#    else:
#        df_temp = pd.DataFrame(solList[i], index=['omega'+str(i)]).T
#
#        df_f = df_f.merge(df_temp, how = 'outer',left_index=True, right_index=True)
#df_new = df_f.reset_index()        
#df_new2 = df_new[df_new['index'].str.contains('M')]
    
        



    



