# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:51:16 2019

@author: atohidi
"""
from gurobipy import *
import numpy as np
import pandas as pd
from tabulate import tabulate
import matplotlib.pyplot as plt
from pprint import pprint
import time

M= 1000 #bigM
eps=0.00001
TL=2000 #model time limit
plotting= True
#defining sets and indecies. 
ip=5 #of potential facilities (cluster of facilities)
j= 3 # of drugs to be tested
t=20  #planning horizon
r=1 # of resources
omega=10 #of scenarios

I= np.arange(ip)
J= np.arange(j)
T= np.arange(t)
R= np.arange(r)
Omega= np.arange(omega)
delta = 2 #time from deciding to open till runing clinical trial
Random = True


if Random:  
    ##parameter random generation
    np.random.seed(25)
    N     =np.random.randint(100,200,size=[j,omega])  # number of patients required for each trial
    rho   =np.random.randint(3,size=[j,r]) #resources required for each test
    rhoMax=np.random.randint(50,100,size=[ip,r]) #max available resource in each center at each time
    #w     =np.random.random([ip,j,t,omega]) # stochastic term representing the proportion of people who actually arrived
    w     =np.ones([ip,j,t,omega])
    p     =np.random.random(size=omega)   # probability of each scenario
    p=np.ones(omega)*1/omega
    p     = p/p.sum() # normalized probabilities
    fp    = np.random.randint(1,100,size=[ip,j])*1000  #fix cost for oppening centers
    fpp   = np.random.randint(1,100,size=[ip,j,t])*1000 #fix cost for opening centers after time zero
    g     = np.random.randint(1,100,size=[ip,j])*100#maintenance cost for keeping once center open 
    theta = np.random.randint(1,10,size=[ip,j,t])*1000 # variable cost (recruiting one patient)
    v     = np.random.randint(1,100,size=[ip,j])*1000 # capacity of each center 
    k     = np.random.randint(100,1000,size=[j,omega])*1000 # drug profit
    l     = np.random.randint(3,5,size=j) # study time
    alpha = np.random.random(size=j)/2
    Tj=[(j,t) for j in J for t in T[l[j]:]]
    #rounding random numbers to two digit
    w = w.round(2)
    alpha = alpha.round(2)
    alpha = np.zeros(len(J))

else: # for parameters given from MSOM
    ip=5 #of potential facilities (cluster of facilities)
    j= 4 # of drugs to be tested
    t=20  #planning horizon
    r=1 # of resources
    omega=10 #of scenarios
    I= np.arange(ip)
    J= np.arange(j)
    T= np.arange(t)
    R= np.arange(r)
    Omega= np.arange(omega)
    N = [2229,1161,5555,564]
    rho = np.ones([j,r])*2
    rhoMax = np.ones([ip,r])*200 
    np.random.seed(25)
    w = np.ones([ip,j,t,omega])
    #  w     =np.random.random([ip,j,t,omega])
    w = w.round(2)
    p = np.random.random(size=omega)   # probability of each scenario
    p = np.ones(omega)*1/omega
    p = p/p.sum() # normalized probabilities
    fp = np.array([[0.96,0.96,0.96,0.96],
                   [1.26,1.26,1.26,1.26],
                   [0.64,0.64,0.64,0.64]])
    theta_temp = np.array([[0.11,0.12,0.01,0.46],
                           [0.16,0.30,0.02,0.83],
                           [0.18,0.39,0.03,1000]])
    theta = np.zeros([len(I),len(J),len(T)])
    for i in I:
        for j in J:
            for t in T:
                theta[i,j,t] = theta_temp[i,j]             
    v = np.ones([len(I),len(J)])*200
    k = [260,800,1195,374]
    l = [3 for j in J]   
    alpha = np.zeros(len(J))

    Tj=[(j,t) for j in J for t in T[int(l[j]):]]
    
def solve(omega1):
    Omega = [omega1]
    model= Model('two stage clinical trial problem')
    
    u  =model.addVars(I,J,T,Omega, name='U', vtype= GRB.INTEGER)
    m  =model.addVars(I,J, name='M', vtype= GRB.BINARY)
    mp =model.addVars(I,J,T,Omega, name='Mp', vtype= GRB.BINARY)
    y  =model.addVars(J,T,Omega, name='Y', vtype= GRB.BINARY)
    d  =model.addVars(J,Omega, name='D', vtype= GRB.INTEGER) # d= Tau
    a  =model.addVars(I,J,T,Omega, name='A', vtype=GRB.INTEGER)
    f  =model.addVars(I,Tj,Omega, name='F', vtype=GRB.INTEGER)
    o =model.addVars(I,J,T,Omega, name='O', vtype= GRB.BINARY)
    c =model.addVars(I,J,T,Omega, name='C', vtype= GRB.BINARY)
    z =model.addVars(J, name='Z', vtype= GRB.BINARY)
    inv=model.addVars(I,J,T,Omega, name='INV', vtype=GRB.INTEGER)
    obj=model.addVar(name='OBJ',vtype=GRB.CONTINUOUS)
    
    model.addConstrs((a[i,j,t,omega]<= w[i,j,t,omega]*u[i,j,t,omega]* (1-alpha[j])**l[j] for i in I for j in J for t in T for omega in Omega),'c1')
    model.addConstrs((f[i,j,t+1,omega]==f[i,j,t,omega]+ a[i,j,t-l[j]+1,omega] for i in I for j in J for t in T[l[j]:-1] for omega in Omega),'c2')
    model.addConstrs((f[i,j,t,omega]== a[i,j,t-l[j],omega] for i in I for j in J for t in T[l[j]:l[j]+1] for omega in Omega),'c3')
    model.addConstrs((inv[i,j,t+1,omega]>= (inv[i,j,t,omega]* (1-alpha[j])-a[i,j,t-l[j],omega]) + w[i,j,t,omega]*u[i,j,t,omega]-M* y[j,t,omega] for i in I for j in J for t in T[l[j]:-1] for omega in Omega),'c4')
    model.addConstrs((inv[i,j,t+1,omega]>= inv[i,j,t,omega]*(1-alpha[j])+ w[i,j,t,omega]*u[i,j,t,omega]  for i in I for j in J for omega in Omega for t in T[0:l[j]]),'c5')
    model.addConstrs((inv[i,j,0,omega]== 0 for i in I for j in J for omega in Omega),'c6')
    model.addConstrs((quicksum(f[i,j,t,omega] for i in I) >= N[j,omega]*y[j,t,omega] for (j,t) in Tj for omega in Omega),'c7')
    model.addConstrs((quicksum(f[i,j,t,omega] for i in I) >= N[j,omega]*z[j] for (j,t) in Tj if t==T[-1] for omega in Omega),'c7_1')
    model.addConstrs((d[j,omega]>= (t+1) * (1-y[j,t,omega]) for t in T for j in J for omega in Omega),'c9')
    model.addConstrs((quicksum(rho[j,r]*inv[i,j,t,omega] for j in J )<= rhoMax[i,r] for i in I for r in R for t in T for omega in Omega),'c10')
    model.addConstrs((u[i,j,t,omega]<= v[i,j]*o[i,j,t,omega] for i in I for j in J for omega in Omega for t in T),'c11')
    model.addConstrs((o[i,j,t,omega]== o[i,j,t-1,omega]+mp[i,j,t-delta,omega]-c[i,j,t,omega] for i in I for j in J for t in T[delta+1:] for omega in Omega),'c12')
    model.addConstrs((o[i,j,delta,omega]==m[i,j] for i in I for j in J for omega in Omega),'c13')
    model.addConstrs((o[i,j,t,omega] ==0 for i in I for j in J for t in T[:delta] for omega in Omega),'c14')
    model.addConstrs((m[i,j]+quicksum( mp[i,j,t,omega] for t in T)<=z[j] for i in I for j in J for omega in Omega),'c15')
    model.addConstrs((c[i,j,t,omega]<= o[i,j,t-1,omega] for i in I for j in J for omega in Omega for t in T[1:]),'c16')
    model.addConstrs((mp[i,j,t,omega]+c[i,j,t,omega]<=1 for i in I for j in J for t in T for omega in Omega),'c18')
    obj1 = quicksum(fp[i,j]*m[i,j] for i in I for j in J)
    obj2 = quicksum(p[omega]*fpp[i,j,t]*mp[i,j,t,omega] for i in I for j in J for t in T for omega in Omega)
    obj3 = quicksum(p[omega]*g[i,j]*o[i,j,t,omega] for i in I for j in J for t in T for omega in Omega)
    obj4 = quicksum(p[omega]*k[j,omega]*d[j,omega] for j in J for omega in Omega)
    obj5 = quicksum(eps*(inv[i,j,t,omega]+u[i,j,t,omega]) for i in I for j in J for t in T for omega in Omega)
    obj  = obj1+obj2+obj3+obj4+obj5
    model.setObjective(obj , GRB.MINIMIZE)
    model.update()
    model.setParam('OutputFlag', 0) #turn off output to console
    model.setParam("TimeLimit", TL)
    model.write('two_stage_closing_decision.lp')
    model.optimize()
    
    Mlist = np.zeros([len(I),len(J)])
    for i in I:
        for j in J:
            Mlist[i,j] = model.getVarByName('M['+str(i) +',' +str(j) +']').x

    return Mlist         


def solveNew(omega1, rhoCoef, wDict, averageList):
    Omega = [omega1]
    model= Model('step6')
    
    u  =model.addVars(I,J,T,Omega, name='U', vtype= GRB.INTEGER)
    m  =model.addVars(I,J, name='M', vtype= GRB.BINARY)
    mp =model.addVars(I,J,T,Omega, name='Mp', vtype= GRB.BINARY)
    y  =model.addVars(J,T,Omega, name='Y', vtype= GRB.BINARY)
    d  =model.addVars(J,Omega, name='D', vtype= GRB.INTEGER) # d= Tau
    a  =model.addVars(I,J,T,Omega, name='A', vtype=GRB.INTEGER)
    f  =model.addVars(I,Tj,Omega, name='F', vtype=GRB.INTEGER)
    o =model.addVars(I,J,T,Omega, name='O', vtype= GRB.BINARY)
    c =model.addVars(I,J,T,Omega, name='C', vtype= GRB.BINARY)
    z =model.addVars(J, name='Z', vtype= GRB.BINARY)
    inv=model.addVars(I,J,T,Omega, name='INV', vtype=GRB.INTEGER)
    obj=model.addVar(name='OBJ',vtype=GRB.CONTINUOUS)
    
    model.addConstrs((a[i,j,t,omega]<= w[i,j,t,omega]*u[i,j,t,omega]* (1-alpha[j])**l[j] for i in I for j in J for t in T for omega in Omega),'c1')
    model.addConstrs((f[i,j,t+1,omega]==f[i,j,t,omega]+ a[i,j,t-l[j]+1,omega] for i in I for j in J for t in T[l[j]:-1] for omega in Omega),'c2')
    model.addConstrs((f[i,j,t,omega]== a[i,j,t-l[j],omega] for i in I for j in J for t in T[l[j]:l[j]+1] for omega in Omega),'c3')
    model.addConstrs((inv[i,j,t+1,omega]>= (inv[i,j,t,omega]* (1-alpha[j])-a[i,j,t-l[j],omega]) + w[i,j,t,omega]*u[i,j,t,omega]-M* y[j,t,omega] for i in I for j in J for t in T[l[j]:-1] for omega in Omega),'c4')
    model.addConstrs((inv[i,j,t+1,omega]>= inv[i,j,t,omega]*(1-alpha[j])+ w[i,j,t,omega]*u[i,j,t,omega]  for i in I for j in J for omega in Omega for t in T[0:l[j]]),'c5')
    model.addConstrs((inv[i,j,0,omega]== 0 for i in I for j in J for omega in Omega),'c6')
    model.addConstrs((quicksum(f[i,j,t,omega] for i in I) >= N[j,omega]*y[j,t,omega] for (j,t) in Tj for omega in Omega),'c7')
    model.addConstrs((quicksum(f[i,j,t,omega] for i in I) >= N[j,omega]*z[j] for (j,t) in Tj if t==T[-1] for omega in Omega),'c7_1')
    model.addConstrs((d[j,omega]>= (t+1) * (1-y[j,t,omega]) for t in T for j in J for omega in Omega),'c9')
    model.addConstrs((quicksum(rho[j,r]*inv[i,j,t,omega] for j in J )<= rhoMax[i,r] for i in I for r in R for t in T for omega in Omega),'c10')
    model.addConstrs((u[i,j,t,omega]<= v[i,j]*o[i,j,t,omega] for i in I for j in J for omega in Omega for t in T),'c11')
    model.addConstrs((o[i,j,t,omega]== o[i,j,t-1,omega]+mp[i,j,t-delta,omega]-c[i,j,t,omega] for i in I for j in J for t in T[delta+1:] for omega in Omega),'c12')
    model.addConstrs((o[i,j,delta,omega]==m[i,j] for i in I for j in J for omega in Omega),'c13')
    model.addConstrs((o[i,j,t,omega] ==0 for i in I for j in J for t in T[:delta] for omega in Omega),'c14')
    model.addConstrs((m[i,j]+quicksum( mp[i,j,t,omega] for t in T)<=z[j] for i in I for j in J for omega in Omega),'c15')
    model.addConstrs((c[i,j,t,omega]<= o[i,j,t-1,omega] for i in I for j in J for omega in Omega for t in T[1:]),'c16')
    model.addConstrs((mp[i,j,t,omega]+c[i,j,t,omega]<=1 for i in I for j in J for t in T for omega in Omega),'c18')
    obj1 = quicksum(fp[i,j]*m[i,j] for i in I for j in J)
    obj2 = quicksum(p[omega]*fpp[i,j,t]*mp[i,j,t,omega] for i in I for j in J for t in T for omega in Omega)
    obj3 = quicksum(p[omega]*g[i,j]*o[i,j,t,omega] for i in I for j in J for t in T for omega in Omega)
    obj4 = quicksum(p[omega]*k[j,omega]*d[j,omega] for j in J for omega in Omega)
    obj5 = quicksum(eps*(inv[i,j,t,omega]+u[i,j,t,omega]) for i in I for j in J for t in T for omega in Omega)
    obj6 = quicksum(wList[omega][i,j]*m[i,j] + rhoCoef/2  * (m[i,j]- averageList[i,j]) * (m[i,j]- averageList[i,j]) for i in I for j in J for omega in Omega)
    obj  = obj1+obj2+obj3+obj4+obj5+obj6
    model.setObjective(obj , GRB.MINIMIZE)
    model.update()
    model.setParam('OutputFlag', 0) #turn off output to console
    model.setParam("TimeLimit", TL)
    model.write('two_stage_closing_decision.lp')
    model.optimize()
    Mlist = np.zeros([len(I),len(J)])
    for i in I:
        for j in J:
            Mlist[i,j] = model.getVarByName('M['+str(i) +',' +str(j) +']').x
    
    sol = {}
    for var in model.getVars():
        if var.x !=0:
            sol.update({var.Varname: var.x}) 
    return Mlist, sol  
 

def averageCalc(p, solList):
    averageList = np.zeros([len(I),len(J)]) #for step 3
    for i in I:
        for j in J:
            averageList[i,j] = sum([p[omega]*solList[omega][i,j] for omega in Omega])   
    return averageList
        
solList = []
counter = 0
rhoCoef = 20000
# step 2: solve for each scenario
for scenario in Omega:
    solList.append(solve(scenario))
    
#step 3: average scenario
averageList = averageCalc(p,solList)

#step 4
wList = [[] for omega in Omega] #for step 4
for omega in Omega:
    wList[omega] = np.zeros([len(I),len(J)])
    for i in I:
        for j in J:
            wList[omega][i,j] = rhoCoef* (solList[omega][i,j] - averageList[i,j])
            
 
# store in df for check:
df0 = pd.DataFrame()
for omega in Omega:
    for i in I:
        df0['omega:'+str(omega)+', i:'+str(i)] = solList[omega][i]  
for i in I:
    df0['average, i:'+str(i)] = averageList[i]  
for omega in Omega:
    for i in I:
        df0['w omega:'+str(omega)+', i:'+str(i)] = wList[omega][i]        


check = True
while check:
    #step 5
    k +=1
    
    #step 6:
    newSolList = []
    for omega in Omega:
        ml,sol = solveNew(omega, rhoCoef, wList, averageList)
        newSolList.append(ml)
    
    #step 7
    newAverageList = averageCalc(p, newSolList) 
    
    #step 8   
    newWList = [[] for omega in Omega] #for step 4
    for omega in Omega:
        newWList[omega] = np.zeros([len(I),len(J)])
        for i in I:
            for j in J:
                newWList[omega][i,j] = wList[omega][i,j]+rhoCoef* (newSolList[omega][i,j] - newAverageList[i,j])
                
    #step 9 
    gSum = 0
    for omega in Omega:
        sumP = 0
        for i in I:
            for j in J:
                sumP += (newSolList[omega][i,j] - newAverageList[i,j])**2
        gSum += p[omega]*np.sqrt(sumP)
    
    if gSum < 0.1: 
        check = False
    else:
        print(gSum)
        wList = newWList 
        solList = newSolList
        averageList = newAverageList
                
                
    