# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 18:48:54 2019

@author: atohidi
"""

#### extensive form
#        
#model= Model('two stage clinical trial problem')
#
#u  =model.addVars(I,J,T,Omega, name='U', vtype= GRB.INTEGER)
#m  =model.addVars(I,J, name='M', vtype= GRB.BINARY)
#mp =model.addVars(I,J,T,Omega, name='Mp', vtype= GRB.BINARY)
#y  =model.addVars(J,T,Omega, name='Y', vtype= GRB.BINARY)
#d  =model.addVars(J,Omega, name='D', vtype= GRB.INTEGER) # d= Tau
#a  =model.addVars(I,J,T,Omega, name='A', vtype=GRB.INTEGER)
#f  =model.addVars(I,Tj,Omega, name='F', vtype=GRB.INTEGER)
#o =model.addVars(I,J,T,Omega, name='O', vtype= GRB.BINARY)
#c =model.addVars(I,J,T,Omega, name='C', vtype= GRB.BINARY)
#z =model.addVars(J, name='Z', vtype= GRB.BINARY)
#inv=model.addVars(I,J,T,Omega, name='INV', vtype=GRB.INTEGER)
#obj=model.addVar(name='OBJ',vtype=GRB.CONTINUOUS)
#
#model.addConstrs((a[i,j,t,omega]<= w[i,j,t,omega]*u[i,j,t,omega]* (1-alpha[j])**l[j] for i in I for j in J for t in T for omega in Omega),'c1')
#model.addConstrs((f[i,j,t+1,omega]==f[i,j,t,omega]+ a[i,j,t-l[j]+1,omega] for i in I for j in J for t in T[l[j]:-1] for omega in Omega),'c2')
#model.addConstrs((f[i,j,t,omega]== a[i,j,t-l[j],omega] for i in I for j in J for t in T[l[j]:l[j]+1] for omega in Omega),'c3')
#model.addConstrs((inv[i,j,t+1,omega]>= (inv[i,j,t,omega]* (1-alpha[j])-a[i,j,t-l[j],omega]) + w[i,j,t,omega]*u[i,j,t,omega]-M* y[j,t,omega] for i in I for j in J for t in T[l[j]:-1] for omega in Omega),'c4')
#model.addConstrs((inv[i,j,t+1,omega]>= inv[i,j,t,omega]*(1-alpha[j])+ w[i,j,t,omega]*u[i,j,t,omega]  for i in I for j in J for omega in Omega for t in T[0:l[j]]),'c5')
#model.addConstrs((inv[i,j,0,omega]== 0 for i in I for j in J for omega in Omega),'c6')
#model.addConstrs((quicksum(f[i,j,t,omega] for i in I) >= N[j,omega]*y[j,t,omega] for (j,t) in Tj for omega in Omega),'c7')
#model.addConstrs((quicksum(f[i,j,t,omega] for i in I) >= N[j,omega]*z[j] for (j,t) in Tj if t==T[-1] for omega in Omega),'c7_1')
#model.addConstrs((d[j,omega]>= (t+1) * (1-y[j,t,omega]) for t in T for j in J for omega in Omega),'c9')
#model.addConstrs((quicksum(rho[j,r]*inv[i,j,t,omega] for j in J )<= rhoMax[i,r] for i in I for r in R for t in T for omega in Omega),'c10')
#model.addConstrs((u[i,j,t,omega]<= v[i,j]*o[i,j,t,omega] for i in I for j in J for omega in Omega for t in T),'c11')
#model.addConstrs((o[i,j,t,omega]== o[i,j,t-1,omega]+mp[i,j,t-delta,omega]-c[i,j,t,omega] for i in I for j in J for t in T[delta+1:] for omega in Omega),'c12')
#model.addConstrs((o[i,j,delta,omega]==m[i,j] for i in I for j in J for omega in Omega),'c13')
#model.addConstrs((o[i,j,t,omega] ==0 for i in I for j in J for t in T[:delta] for omega in Omega),'c14')
#model.addConstrs((m[i,j]+quicksum( mp[i,j,t,omega] for t in T)<=z[j] for i in I for j in J for omega in Omega),'c15')
#model.addConstrs((c[i,j,t,omega]<= o[i,j,t-1,omega] for i in I for j in J for omega in Omega for t in T[1:]),'c16')
#model.addConstrs((mp[i,j,t,omega]+c[i,j,t,omega]<=1 for i in I for j in J for t in T for omega in Omega),'c18')
#model.addConstrs((y[j,t,omega] <= z[j] for j in J for t in T for omega in Omega),'new cons')
#
#obj1 = quicksum(fp[i,j]*m[i,j] for i in I for j in J)
#obj2 = quicksum(p[omega]*fpp[i,j,t]*mp[i,j,t,omega] for i in I for j in J for t in T for omega in Omega)
#obj3 = quicksum(p[omega]*g[i,j]*o[i,j,t,omega] for i in I for j in J for t in T for omega in Omega)
#obj4 = quicksum(p[omega]*k[j,omega]*(d[j,omega]- (T[-1]+1)*(1-z[j])) for j in J for omega in Omega)
#obj5 = quicksum(eps*(inv[i,j,t,omega]+u[i,j,t,omega]) for i in I for j in J for t in T for omega in Omega)
#obj6 = quicksum((1-z[j])*not_testing_cost[j] for j in J)     # loosing lump sum cost if not testing a drug at all
#obj  = obj1+obj2+obj3+obj4+obj5+obj6
#model.setObjective(obj , GRB.MINIMIZE)
#model.update()
#print('****', model.NumVars,model.NumConstrs)
#
##model.setParam('OutputFlag', 0) #turn off output to console
#model.setParam("TimeLimit", TL)
#model.setParam("MIPGap", 0.06)
#
#model.write('two_stage_closing_decision.lp')
#model.optimize()
#
#
#counter = 0
#counterMP = 0
#for i in model.getVars():
#    if i.x!= 0 and 'M' in i.varName:
#        counter += 1
#        if 'Mp' in i.varname:
#            counterMP+=1
#print(counter-counterMP, counterMP/omega)
#incomp =0
#for i in model.getVars():
#    if i.x!= 0 and 'D' in i.varName:
#        if i.x == t:
#            incomp +=1
#            
#        
#        print(i,i.x)
#
#print(model.objVal)

#Ftime=time.time()
#Dtime=Ftime-Stime
##printing outputs and plotting
#j=len(J)
#t=len(T)
#sol={}
#mlist=np.zeros([ip,j])  
#udict={}
#Alist=np.zeros([ip,j,t,omega]) 
#fdict={}
#dDict={}
#Ilist=np.zeros([ip,j,t,omega]) 
#
#flist=np.zeros([ip,j,t,omega])  
#
#sumx=0
#appended_data=[]
#mpdict={}
#odict={}
#cdict={}
#zdict={}
#
#if model.Status!=3 and plotting == True: #status =3 --> infeasible
#    for v in model.getVars():
#        if v.x !=0:
#            sol.update({v.Varname: v.x})
#        if 'U' in v.Varname and v.x !=0:
#            udict[v.Varname] = v.x  
#            sumx+=v.x
#                
#        if 'F' in v.Varname :
#            name=v.Varname
#            nl= name[2:-1].split(',')
#            fdict[v.varname]=v.x
#            flist[int(nl[0]),int(nl[1]),int(nl[2]),int(nl[3])]= v.x
#        if 'A' in v.Varname:
#            name=v.Varname
#            nl= name[2:-1].split(',')
#            Alist[int(nl[0]),int(nl[1]),int(nl[2]),int(nl[3])]= v.x
#        if 'I' in v.Varname:
#            name=v.Varname
#            nl= name[4:-1].split(',')
#            Ilist[int(nl[0]),int(nl[1]),int(nl[2]),int(nl[3])]= v.x
#            appended_data.append(pd.Series([int(nl[0]),int(nl[1]),int(nl[2]),int(nl[3]),v.x],index=['i','j','t','omega','Inv'] ,dtype='int32' ))
#            
#        if 'M[' in v.Varname:
#            name=v.Varname
#            nl= name[2:-1].split(',')
#            mlist[int(nl[0])][int(nl[1])]= v.x
#        
#        if 'Mp' in v.Varname:
#            name=v.Varname
#            mpdict[name]= v.x
#        if 'O' in v.Varname:
#            name=v.Varname
#            odict[name]= v.x
#        if 'C' in v.Varname:
#            name=v.Varname
#            cdict[name]= v.x    
#        if 'D' in v.Varname:
#            name=v.Varname
#            dDict[v.varname]=v.x
#        if 'Z' in v.Varname:
#            name=v.Varname
#            zdict[v.varname]=v.x    
#    appended_data = pd.concat(appended_data, axis=1)
#    I_df0=appended_data.transpose()
#    
#    M_df=pd.DataFrame(mlist,columns=np.arange(j),dtype='int32')            
#    print('\n************************************** \n\nThe optimal M[i,j] is as follow : \n ')
#    print(tabulate(M_df, headers='keys', tablefmt='psql'))        
#    print(f'\n  P.s. # of centers : {ip}, # of drugs: {j} \n\n**************************************')    
#    print('tau:\n',dDict)
#      #select which center and drug to be plotted
#    jlist=[0]
#    ilist=[0]
#    omegalist=Omega
#    
#    plt.grid() 
#    plt.subplot(3,1,1)
#    for j1 in jlist:
#        for i1 in ilist:
#            for omega1 in omegalist:
#                f_df=pd.DataFrame(data=flist[i1,j1,:,omega1])
#                plt.plot(f_df, label= f'[{i1},{j1},{omega1}]')
#        
#    plt.legend()
#    plt.xticks(T[::3])
#    plt.title(f'Cumulative number of patients in center {ilist} who finished clinical trial {jlist} in scenario {omegalist}')
#    enddict={} # to save fdict when t=T
#    for i in fdict:
#        if i.split(',')[2]== str(t-1):  
#           enddict[i]=fdict[i] 
#    pprint(enddict)
#    plt.subplot(3,1,2)    
#    plt.grid() 
#        
#    for j1 in jlist:
#        for i1 in ilist:
#            for omega1 in omegalist:
#                A_df=pd.DataFrame(data=Alist[i1,j1,:,omega1])
#                plt.plot(A_df, label= f'[{i1},{j1},{omega1}]')
#    plt.legend()
#    plt.xticks(T[::3])
#    plt.title(f'Patient arrival to center {ilist} for clinical trial {jlist} in scenario {omegalist}')
#    
#    plt.subplot(3,1,3)        
#    plt.grid() 
#    
#    for j1 in jlist:
#        for i1 in ilist:
#            for omega1 in omegalist:
#                I_df=pd.DataFrame(data=Ilist[i1,j1,:,omega1])
#                plt.plot(I_df, label= f'[{i1},{j1},{omega1}]')
#    plt.legend()
#    plt.xticks(T[::3])
#    plt.title(f'Patient existing in center {ilist} for clinical trial {jlist} in scenario {omegalist}')
#    plt.grid() 
#    
#    s=0
#    # check if the number of people for one specific drug is sufficient.
#    for i in enddict:
#        if i[2:-1].split(',')[1]== '1' and i[2:-1].split(',')[3]== '0':
#            s+= enddict[i]
#    
#    
#    df2=I_df0.loc[(I_df0.i==0) & (I_df0.omega==0)]  
#    gI= I_df0.groupby(by=['j', 't','omega'],axis=0).sum() 
#          
#else:
#    print(f'model status is {model.Status}')
#    if model.Status== 9:
#        print(f"time limit of {TL} reached" )
#    
#    if model.Status== 3:
#        print('model is infeasible')
#        
#
#
#    Al=pd.Series([Alist[0,0,t,0] for t in T], dtype='int32' )
#    Il=pd.Series([Ilist[0,0,t,0] for t in T] )
#    Fl=pd.Series([flist[0,0,t,0] for t in T], dtype='int32' )
#    Ul=pd.Series([udict.get( f'U[0,0,{t},0]' ,0 ) for t in T], dtype='int32' )
#    wl=pd.Series(w[0,0,:,0].T)
#    
#    L= pd.DataFrame([ wl, Ul, Al,Il,Fl]).T
#    L.columns=['W%','U','A','I','F']
##format all columns but w to int
##for i in L.columns:
##    if i !='W%' or i != "I":
##          L[i]=L[i].astype('int')
#        
#print(f'Elapsed time: {Dtime}')
#print(f'MILP GAP at termination: {model.MIPGap}')
#print(f'Best objective found = {model.ObjVal}')
#
#Idf1=pd.DataFrame(data=[Ilist[i,j,:,0] for j in J for i in I]).T

#print('****', model.NumVars,model.NumConstrs)
#print('incomp =' , incomp)
