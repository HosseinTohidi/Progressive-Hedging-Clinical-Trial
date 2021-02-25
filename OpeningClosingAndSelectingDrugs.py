#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 12:14:10 2019

@author: hosseintohidi
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
ip=10 #of potential facilities (cluster of facilities)
j= 5 # of drugs to be tested
t=20  #planning horizon
r=1 # of resources
omega=1 #of scenarios

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
    v = np.ones([len(I),len(J)])*1000
    k = [260,800,1195,374]
    l = [3 for j in J]   
    alpha = np.zeros(len(J))

    Tj=[(j,t) for j in J for t in T[int(l[j]):]]
    
Stime=time.time()
## Building the model

model= Model('two stage clinical trial problem')

u  =model.addVars(I,J,T,Omega, name='U', vtype= GRB.INTEGER)
m  =model.addVars(I,J, name='M', vtype= GRB.BINARY)
mp =model.addVars(I,J,T,Omega, name='Mp', vtype= GRB.BINARY)
y  =model.addVars(J,T,Omega, name='Y', vtype= GRB.BINARY)
d  =model.addVars(J,Omega, name='D', vtype= GRB.INTEGER) # d= Tau
a  =model.addVars(I,J,T,Omega, name='A', vtype=GRB.INTEGER)
#f  =model.addVars(I,J,T,Omega, name='F', vtype=GRB.INTEGER)
f  =model.addVars(I,Tj,Omega, name='F', vtype=GRB.INTEGER)
o =model.addVars(I,J,T,Omega, name='O', vtype= GRB.BINARY)
c =model.addVars(I,J,T,Omega, name='C', vtype= GRB.BINARY)
z =model.addVars(J, name='Z', vtype= GRB.BINARY)

inv=model.addVars(I,J,T,Omega, name='INV', vtype=GRB.INTEGER)
obj=model.addVar(name='OBJ',vtype=GRB.CONTINUOUS)

model.addConstrs((a[i,j,t,omega]<= w[i,j,t,omega]*u[i,j,t,omega]* (1-alpha[j])**l[j] for i in I for j in J for t in T for omega in Omega),'c1')

#model.addConstrs((f[i,j,t,omega]<=f[i,j,t-1,omega]+ a[i,j,t-l[j],omega] for i in I for j in J for t in T[l[j]:] for omega in Omega),'c2')
model.addConstrs((f[i,j,t+1,omega]==f[i,j,t,omega]+ a[i,j,t-l[j]+1,omega] for i in I for j in J for t in T[l[j]:-1] for omega in Omega),'c2')

#model.addConstrs((f[i,j,t,omega]== 0 for i in I for j in J for t in T[0:l[j]] for omega in Omega),'c3')
model.addConstrs((f[i,j,t,omega]== a[i,j,t-l[j],omega] for i in I for j in J for t in T[l[j]:l[j]+1] for omega in Omega),'c3')

model.addConstrs((inv[i,j,t+1,omega]>= (inv[i,j,t,omega]* (1-alpha[j])-a[i,j,t-l[j],omega]) + w[i,j,t,omega]*u[i,j,t,omega]-M* y[j,t,omega] for i in I for j in J for t in T[l[j]:-1] for omega in Omega),'c4')

model.addConstrs((inv[i,j,t+1,omega]>= inv[i,j,t,omega]*(1-alpha[j])+ w[i,j,t,omega]*u[i,j,t,omega]  for i in I for j in J for omega in Omega for t in T[0:l[j]]),'c5')

model.addConstrs((inv[i,j,0,omega]== 0 for i in I for j in J for omega in Omega),'c6')

#model.addConstrs((quicksum(f[i,j,t,omega] for i in I) >= N[j]*y[j,t,omega] for j in J for t in T for omega in Omega),'c7')
model.addConstrs((quicksum(f[i,j,t,omega] for i in I) >= N[j,omega]*y[j,t,omega] for (j,t) in Tj for omega in Omega),'c7')

model.addConstrs((quicksum(f[i,j,t,omega] for i in I) >= N[j,omega]*z[j] for (j,t) in Tj if t==T[-1] for omega in Omega),'c7_1')

#model.addConstrs((quicksum(f[i,j,T[-1],omega] for i in I) >= N[j] for j in J for omega in Omega),'c8') #no need to have it
#model.addConstrs((d[j,omega]== quicksum(1-y[j,t,omega] for t in T) for j in J for omega in Omega),'c9') #replaced by other c9 to test performance
#new constraint :: I want to remove I  from objective function and still get zero Inv when study is done
#model.addConstrs(( inv[i,j,t+1,omega]  <= M* (1-y[j,t,omega]) for i in I for j in J for t in T[:-1] for omega in Omega),'C7_1')

model.addConstrs((d[j,omega]>= (t+1) * (1-y[j,t,omega]) for t in T for j in J for omega in Omega),'c9')

model.addConstrs((quicksum(rho[j,r]*inv[i,j,t,omega] for j in J )<= rhoMax[i,r] for i in I for r in R for t in T for omega in Omega),'c10')

model.addConstrs((u[i,j,t,omega]<= v[i,j]*o[i,j,t,omega] for i in I for j in J for omega in Omega for t in T),'c11')

#model.addConstrs((2*o[i,j,t,omega]>= o[i,j,t-1,omega]+mp[i,j,t-delta,omega]-c[i,j,t,omega] for i in I for j in J for t in T[delta+1:] for omega in Omega),'c12')
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
#model.setParam('OutputFlag', 0) #turn off output to console
model.setParam("TimeLimit", TL)
model.write('two_stage_closing_decision.lp')
model.optimize()

Ftime=time.time()
Dtime=Ftime-Stime
#printing outputs and plotting
j=len(J)
t=len(T)
sol={}
mlist=np.zeros([ip,j])  
udict={}
Alist=np.zeros([ip,j,t,omega]) 
fdict={}
dDict={}
Ilist=np.zeros([ip,j,t,omega]) 

flist=np.zeros([ip,j,t,omega])  

sumx=0
appended_data=[]
mpdict={}
odict={}
cdict={}
zdict={}

if model.Status!=3 and plotting == True: #status =3 --> infeasible
    for v in model.getVars():
        if v.x !=0:
            sol.update({v.Varname: v.x})
        if 'U' in v.Varname and v.x !=0:
            udict[v.Varname] = v.x  
            sumx+=v.x
                
        if 'F' in v.Varname :
            name=v.Varname
            nl= name[2:-1].split(',')
            fdict[v.varname]=v.x
            flist[int(nl[0]),int(nl[1]),int(nl[2]),int(nl[3])]= v.x
        if 'A' in v.Varname:
            name=v.Varname
            nl= name[2:-1].split(',')
            Alist[int(nl[0]),int(nl[1]),int(nl[2]),int(nl[3])]= v.x
        if 'I' in v.Varname:
            name=v.Varname
            nl= name[4:-1].split(',')
            Ilist[int(nl[0]),int(nl[1]),int(nl[2]),int(nl[3])]= v.x
            appended_data.append(pd.Series([int(nl[0]),int(nl[1]),int(nl[2]),int(nl[3]),v.x],index=['i','j','t','omega','Inv'] ,dtype='int32' ))
            
        if 'M[' in v.Varname:
            name=v.Varname
            nl= name[2:-1].split(',')
            mlist[int(nl[0])][int(nl[1])]= v.x
        
        if 'Mp' in v.Varname:
            name=v.Varname
            mpdict[name]= v.x
        if 'O' in v.Varname:
            name=v.Varname
            odict[name]= v.x
        if 'C' in v.Varname:
            name=v.Varname
            cdict[name]= v.x    
        if 'D' in v.Varname:
            name=v.Varname
            dDict[v.varname]=v.x
        if 'Z' in v.Varname:
            name=v.Varname
            zdict[v.varname]=v.x    
    appended_data = pd.concat(appended_data, axis=1)
    I_df0=appended_data.transpose()
    
    M_df=pd.DataFrame(mlist,columns=np.arange(j),dtype='int32')            
    print('\n************************************** \n\nThe optimal M[i,j] is as follow : \n ')
    print(tabulate(M_df, headers='keys', tablefmt='psql'))        
    print(f'\n  P.s. # of centers : {ip}, # of drugs: {j} \n\n**************************************')    
    print('tau:\n',dDict)
      #select which center and drug to be plotted
    jlist=[0]
    ilist=[0]
    omegalist=Omega
    
    plt.grid() 
    plt.subplot(3,1,1)
    for j1 in jlist:
        for i1 in ilist:
            for omega1 in omegalist:
                f_df=pd.DataFrame(data=flist[i1,j1,:,omega1])
                plt.plot(f_df, label= f'[{i1},{j1},{omega1}]')
        
    plt.legend()
    plt.xticks(T[::3])
    plt.title(f'Cumulative number of patients in center {ilist} who finished clinical trial {jlist} in scenario {omegalist}')
    enddict={} # to save fdict when t=T
    for i in fdict:
        if i.split(',')[2]== str(t-1):  
           enddict[i]=fdict[i] 
    pprint(enddict)
    plt.subplot(3,1,2)    
    plt.grid() 
        
    for j1 in jlist:
        for i1 in ilist:
            for omega1 in omegalist:
                A_df=pd.DataFrame(data=Alist[i1,j1,:,omega1])
                plt.plot(A_df, label= f'[{i1},{j1},{omega1}]')
    plt.legend()
    plt.xticks(T[::3])
    plt.title(f'Patient arrival to center {ilist} for clinical trial {jlist} in scenario {omegalist}')
    
    plt.subplot(3,1,3)        
    plt.grid() 
    
    for j1 in jlist:
        for i1 in ilist:
            for omega1 in omegalist:
                I_df=pd.DataFrame(data=Ilist[i1,j1,:,omega1])
                plt.plot(I_df, label= f'[{i1},{j1},{omega1}]')
    plt.legend()
    plt.xticks(T[::3])
    plt.title(f'Patient existing in center {ilist} for clinical trial {jlist} in scenario {omegalist}')
    plt.grid() 
    
    s=0
    # check if the number of people for one specific drug is sufficient.
    for i in enddict:
        if i[2:-1].split(',')[1]== '1' and i[2:-1].split(',')[3]== '0':
            s+= enddict[i]
    
    
    df2=I_df0.loc[(I_df0.i==0) & (I_df0.omega==0)]  
    gI= I_df0.groupby(by=['j', 't','omega'],axis=0).sum() 
          
else:
    print(f'model status is {model.Status}')
    if model.Status== 9:
        print(f"time limit of {TL} reached" )
    
    if model.Status== 3:
        print('model is infeasible')
        


    Al=pd.Series([Alist[0,0,t,0] for t in T], dtype='int32' )
    Il=pd.Series([Ilist[0,0,t,0] for t in T] )
    Fl=pd.Series([flist[0,0,t,0] for t in T], dtype='int32' )
    Ul=pd.Series([udict.get( f'U[0,0,{t},0]' ,0 ) for t in T], dtype='int32' )
    wl=pd.Series(w[0,0,:,0].T)
    
    L= pd.DataFrame([ wl, Ul, Al,Il,Fl]).T
    L.columns=['W%','U','A','I','F']
#format all columns but w to int
#for i in L.columns:
#    if i !='W%' or i != "I":
#          L[i]=L[i].astype('int')
        
print(f'Elapsed time: {Dtime}')
print(f'MILP GAP at termination: {model.MIPGap}')
print(f'Best objective found = {model.ObjVal}')

Idf1=pd.DataFrame(data=[Ilist[i,j,:,0] for j in J for i in I]).T

