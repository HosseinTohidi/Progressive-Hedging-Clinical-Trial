# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:55:30 2019

@author: atohidi
"""

from gurobipy import *
import numpy as np
def solve(omega1, I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL, fixed_vars = {},extra_print= False):
    if type(omega1) != list and type(omega1) != np.ndarray :
        Omega = [omega1]
    else:
        Omega = omega1
    model= Model('two stage clinical trial problem')
    
    u  =model.addVars(I,J,T,Omega, name='U', vtype= GRB.INTEGER) #
    m  =model.addVars(I,J, name='M', vtype= GRB.BINARY)#
    mp =model.addVars(I,J,T,Omega, name='Mp', vtype= GRB.BINARY)#
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
    model.addConstrs((y[j,t,omega] <= z[j] for j in J for t in T for omega in Omega),'new_cons')
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
    obj4 = quicksum(p[omega]*k[j,omega]*(d[j,omega]- (T[-1]+1)*(1-z[j])) for j in J for omega in Omega)
    obj5 = quicksum(eps*(inv[i,j,t,omega]+u[i,j,t,omega]) for i in I for j in J for t in T for omega in Omega)
    obj6 = quicksum((1-z[j])*not_testing_cost[j] for j in J)     # loosing lump sum cost if not testing a drug at all
    
    for var in fixed_vars.keys():  # if fixed_var.key is empty -> nothing happens. 
        if 'M[' in var:
            varL = var.split(',')
            i = int(varL[0][2:])
            j = int(varL[1][:-1])
            model.addConstr((m[i,j] == fixed_vars[var]), 'new cons')
        elif 'Z' in var:
            j = int(var[2:-1])
            model.addConstr((z[j] == fixed_vars[var]),'new cons2')

    obj  = obj1+obj2+obj3+obj4+obj5+obj6
    model.setObjective(obj , GRB.MINIMIZE)
    model.update()
    #model.setParam('OutputFlag', 0) #turn off output to console
    #model.setParam('MIPGap' , 0.01)
    model.setParam("TimeLimit", TL)
    model.write('two_stage_closing_decision.lp')
    if extra_print:
        print(f'***Cons: {model.NumConstrs}, vars: {model.NumVars}')
    model.optimize()
    sol = {}
    for var in model.getVars():
        if var.x !=0:
            sol.update({var.Varname: round(var.x,5)})
    
    first_stage_var = {}
    for var in model.getVars():
        if 'M[' in var.varname or 'Z' in var.varname :
            first_stage_var[var.varname] = round(var.x,5)

    return first_stage_var,sol 


def solveNew(clNum, omega1, rhoCoef, wDict_s, averageDict, fixed_var, I, J, T, R, w,
             alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL, rho_dict = {}):
    if type(omega1) != list and type(omega1) != np.ndarray:
        Omega = [omega1]
    else:
        Omega = omega1
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
    model.addConstrs((y[j,t,omega] <= z[j] for j in J for t in T for omega in Omega),'new_cons')
    model.addConstrs((d[j,omega]>= (t+1) * (1-y[j,t,omega]) for t in T for j in J for omega in Omega),'c9')
    model.addConstrs((quicksum(rho[j,r]*inv[i,j,t,omega] for j in J )<= rhoMax[i,r] for i in I for r in R for t in T for omega in Omega),'c10')
    model.addConstrs((u[i,j,t,omega]<= v[i,j]*o[i,j,t,omega] for i in I for j in J for omega in Omega for t in T),'c11')
    model.addConstrs((o[i,j,t,omega]== o[i,j,t-1,omega]+mp[i,j,t-delta,omega]-c[i,j,t,omega] for i in I for j in J for t in T[delta+1:] for omega in Omega),'c12')
    model.addConstrs((o[i,j,delta,omega]==m[i,j] for i in I for j in J for omega in Omega),'c13')
    model.addConstrs((o[i,j,t,omega] ==0 for i in I for j in J for t in T[:delta] for omega in Omega),'c14')
    model.addConstrs((m[i,j]+quicksum(mp[i,j,t,omega] for t in T)<=z[j] for i in I for j in J for omega in Omega),'c15')
    model.addConstrs((c[i,j,t,omega]<= o[i,j,t-1,omega] for i in I for j in J for omega in Omega for t in T[1:]),'c16')
    model.addConstrs((mp[i,j,t,omega]+c[i,j,t,omega]<=1 for i in I for j in J for t in T for omega in Omega),'c18')
    
    for var in fixed_var.keys():
        if 'M[' in var:
            varL= var.split(',')
            i = int(varL[0][2:])
            j = int(varL[1][:-1])
            model.addConstr(m[i,j] == fixed_var[var],'new_m_cons')
#            if fixed_var[var]==1:
#                model.addConstr(z[j] == 1)
        elif 'Z' in var:
            j = int(var[2:-1])
            model.addConstr(z[j] == fixed_var[var])

    
    obj1 = quicksum(fp[i,j]*m[i,j] for i in I for j in J)
    obj2 = quicksum(p[omega]*fpp[i,j,t]*mp[i,j,t,omega] for i in I for j in J for t in T for omega in Omega)
    obj3 = quicksum(p[omega]*g[i,j]*o[i,j,t,omega] for i in I for j in J for t in T for omega in Omega)
    obj4 = quicksum(p[omega]*k[j,omega]*(d[j,omega]- (T[-1]+1)*(1-z[j])) for j in J for omega in Omega)
    obj5 = quicksum(eps*(inv[i,j,t,omega]+u[i,j,t,omega]) for i in I for j in J for t in T for omega in Omega)
    obj6 = quicksum((1-z[j])*not_testing_cost[j] for j in J)     # loosing lump sum cost if not testing a drug at all
    if rho_dict == {}:
        if clNum != -1:
            obj7 = quicksum(wDict_s[('M['+str(i)+','+str(j)+']',clNum)]*m[i,j] + rhoCoef/2  * 
                                    (m[i,j]- averageDict['M['+str(i)+','+str(j)+']']) * (m[i,j]- averageDict['M['+str(i)+','+str(j)+']']) for i in I for j in J)
            obj8 = quicksum(wDict_s[('Z['+str(j)+']',clNum)]*z[j] + rhoCoef/2  * 
                                    (z[j]- averageDict['Z['+str(j)+']']) * (z[j]- averageDict['Z['+str(j)+']']) for i in I for j in J)
        else:
            obj7 = quicksum(wDict_s[('M['+str(i)+','+str(j)+']',omega)]*m[i,j] + rhoCoef/2  * 
                                (m[i,j]- averageDict['M['+str(i)+','+str(j)+']']) * (m[i,j]- averageDict['M['+str(i)+','+str(j)+']']) for i in I for j in J for omega in Omega)
            obj8 = quicksum(wDict_s[('Z['+str(j)+']',omega)]*z[j] + rhoCoef/2  * 
                                    (z[j]- averageDict['Z['+str(j)+']']) * (z[j]- averageDict['Z['+str(j)+']']) for i in I for j in J for omega in Omega)
    else:
        if clNum != -1:
            obj7 = quicksum(wDict_s[('M['+str(i)+','+str(j)+']',clNum)]*m[i,j] + rho_dict['M['+str(i)+','+str(j)+']'] /2  * 
                                    (m[i,j]- averageDict['M['+str(i)+','+str(j)+']']) * (m[i,j]- averageDict['M['+str(i)+','+str(j)+']']) for i in I for j in J)
            obj8 = quicksum(wDict_s[('Z['+str(j)+']',clNum)]*z[j] + rho_dict['Z['+str(j)+']']/2  * 
                                    (z[j]- averageDict['Z['+str(j)+']']) * (z[j]- averageDict['Z['+str(j)+']']) for i in I for j in J)
        else:
            obj7 = quicksum(wDict_s[('M['+str(i)+','+str(j)+']',omega)]*m[i,j] + rho_dict['M['+str(i)+','+str(j)+']']/2  * 
                                (m[i,j]- averageDict['M['+str(i)+','+str(j)+']']) * (m[i,j]- averageDict['M['+str(i)+','+str(j)+']']) for i in I for j in J for omega in Omega)
            obj8 = quicksum(wDict_s[('Z['+str(j)+']',omega)]*z[j] + rho_dict['Z['+str(j)+']']/2  * 
                                    (z[j]- averageDict['Z['+str(j)+']']) * (z[j]- averageDict['Z['+str(j)+']']) for i in I for j in J for omega in Omega)
    
    obj  = obj1+obj2+obj3+obj4+obj5+obj6+obj7+obj8
    model.setObjective(obj , GRB.MINIMIZE)
    model.update()
    model.setParam('MIPGap' , 0.01)
   # model.setParam('OutputFlag', 0) #turn off output to console
    model.setParam("TimeLimit", TL)
    model.write('solve-new.lp')
    model.optimize()
    if model.status == 3:
        print(f'in solving scenario {omega1} model get infeasible')
        print(fixed_var)
        del(fixed_var[list(fixed_var.keys())[0]])
        return solveNew(clNum, omega1, rhoCoef, wDict_s, averageDict, fixed_var,
                        I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL,rho_dict)
        
    
    else:
        first_stage_var = {}
        for var in model.getVars():
            if 'M[' in var.varname or 'Z' in var.varname :
                first_stage_var[var.varname] = round(var.x,5)
        
        sol = {}
        for var in model.getVars():
            if var.x !=0:
                sol.update({var.Varname: round(var.x,5)}) 
        return first_stage_var, sol, model.objVal 


def averageCalc(p, solList):
    averageDict = {}
    for var in solList[0].keys():
        averageDict[var] = round(sum([p[omega]*solList[omega][var] for omega in range(len(solList))]),5)
    return averageDict

def is_inList(i,H,Hlist):  # to be completed
 return False


def solve_for_bound(omega1, I, J, T, R, w, alpha, not_testing_cost, k, N, eps, g, v, fp, fpp, delta, rho, rhoMax, Np, p, Tj, l, M, TL, fixed_vars = {}):
    ''' This model we are fixing some vars (fixed_vars: Dict) 
    the objective function is multiplied by p for the first stage vars  
    '''
    if type(omega1) != list and type(omega1) != np.ndarray :
        Omega = [omega1]
    else:
        Omega = omega1
    model= Model('UB') 
    
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
    model.addConstrs((y[j,t,omega] <= z[j] for j in J for t in T for omega in Omega),'new_cons')
    model.addConstrs((d[j,omega]>= (t+1) * (1-y[j,t,omega]) for t in T for j in J for omega in Omega),'c9')
    model.addConstrs((quicksum(rho[j,r]*inv[i,j,t,omega] for j in J )<= rhoMax[i,r] for i in I for r in R for t in T for omega in Omega),'c10')
    model.addConstrs((u[i,j,t,omega]<= v[i,j]*o[i,j,t,omega] for i in I for j in J for omega in Omega for t in T),'c11')
    model.addConstrs((o[i,j,t,omega]== o[i,j,t-1,omega]+mp[i,j,t-delta,omega]-c[i,j,t,omega] for i in I for j in J for t in T[delta+1:] for omega in Omega),'c12')
    model.addConstrs((o[i,j,delta,omega]==m[i,j] for i in I for j in J for omega in Omega),'c13')
    model.addConstrs((o[i,j,t,omega] ==0 for i in I for j in J for t in T[:delta] for omega in Omega),'c14')
    model.addConstrs((m[i,j]+quicksum( mp[i,j,t,omega] for t in T)<=z[j] for i in I for j in J for omega in Omega),'c15')
    model.addConstrs((c[i,j,t,omega]<= o[i,j,t-1,omega] for i in I for j in J for omega in Omega for t in T[1:]),'c16')
    model.addConstrs((mp[i,j,t,omega]+c[i,j,t,omega]<=1 for i in I for j in J for t in T for omega in Omega),'c18')
    obj1 = quicksum(p[omega]*fp[i,j]*m[i,j] for i in I for j in J for omega in Omega)
    obj2 = quicksum(p[omega]*fpp[i,j,t]*mp[i,j,t,omega] for i in I for j in J for t in T for omega in Omega)
    obj3 = quicksum(p[omega]*g[i,j]*o[i,j,t,omega] for i in I for j in J for t in T for omega in Omega)
    obj4 = quicksum(p[omega]*k[j,omega]*(d[j,omega]- (T[-1]+1)*(1-z[j])) for j in J for omega in Omega)
    obj5 = quicksum(eps*(inv[i,j,t,omega]+u[i,j,t,omega]) for i in I for j in J for t in T for omega in Omega)
    obj6 = quicksum(p[omega]*(1-z[j])*not_testing_cost[j] for j in J for omega in Omega)     # loosing lump sum cost if not testing a drug at all
    
    for var in fixed_vars.keys():  # if fixed_var.key is empty -> nothing happens. 
        if 'M[' in var:
            varL = var.split(',')
            i = int(varL[0][2:])
            j = int(varL[1][:-1])
            model.addConstr((m[i,j] == fixed_vars[var]), 'new cons')
        elif 'Z' in var:
            j = int(var[2:-1])
            model.addConstr((z[j] == fixed_vars[var]),'new cons2')

    obj  = obj1+obj2+obj3+obj4+obj5+obj6
    model.setObjective(obj , GRB.MINIMIZE)
    model.update()
    model.setParam('OutputFlag', 0) #turn off output to console
    #model.setParam("TimeLimit", TL)
    model.write('UB.lp')
    model.optimize()
    return model.objVal




