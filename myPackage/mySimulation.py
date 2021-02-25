# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 16:00:22 2019

@author: atohidi
"""
import numpy as np
from functools import reduce
from scipy.stats import beta
from scipy.special import beta as beta_func
import statsmodels.stats.proportion as Mn
import operator as op
import time
import progressbar

def mySimulation(p,p_hat = 0.5, etta = 0.95, ppos = 0.6, p_lb = 0.2, N = 100, M = 59, scenario_number = 50, repl_no = 100, plot = False):  
    '''
    p: actual p
    etta: confidence level
    ppos: initial guess for ppos
    p_lb: stopping threshold on posterior prob    
    N: total number of observations for testing a particular drug
    M: total required successful observation needed to conclude success in trial
    scenario_number = 100 # how many instances to be generated in each replication
    repl_no: how many replication
    plot: bolean to turn on or off plotting 
    '''
    N1 = 0.10 * N #first interim point
    N2 = 0.50 * N #second interim point
    # calculate nCr
    def ncr(n, r):
        r = min(r, n-r)
        numer = reduce(op.mul, range(n, n-r, -1), 1)
        denom = reduce(op.mul, range(1, r+1), 1)
        return numer / denom   
    def patient_result(p):
        if np.random.rand() < p:
            return True #success
        else:
            return False #failure   
    def f(k,n,a,b):   #beta-binomial distriution pdf
        return ncr(n,k)* beta_func(k+a,n-k+b)/beta_func(a,b)    
    def ppos_calc(M,N,Alpha,Beta,p=0.5):  #calculate posterior probability
        n= Alpha+Beta #number of observation until the interim analysis
        temp=0
        for k in range(M-Alpha,N-n):
            temp+= f(k,N-n,Alpha+1,Beta+1) 
        rv = beta(Alpha, Beta)    
        pr = 1 - rv.cdf(p) # we need it to be over 0.95
        return (temp,pr)  
    
    def scenario_creation(ppos):
        counter=0
        (succ,fail)= (0,0)
        while ppos>=p_lb and counter<N:
            #print(succ,fail)
            counter+=1 #recruit one more patient
            if patient_result(p): 
                succ += 1
            else: 
                fail += 1           
            if counter == N1: # interim point
                ppos,post_prob = ppos_calc(M,N,succ,fail,p)
                if ppos<p_lb:
                    temp = 'F1'
                    break
            elif counter == N2: # interim point
                ppos,post_prob = ppos_calc(M,N,succ,fail,p)
                if ppos<p_lb:
                    temp = 'F2'
                    break       
            elif counter == N:
                ppos,post_prob = ppos_calc(M,N,succ,fail,p_hat)
                if succ>M:

                    temp= 'S'
                else:
                    temp = 'F3'
        return temp            
    F1=F2=F3=S=0 ; F1list=[]
    F2list=[] ;F3list=[]; Slist=[];f1list=[]
    f2list=[]; f3list=[]; slist=[]
    ci_l_list=[];f1_ci_l = []
    f1_ci_u =[];f2_ci_l = []
    f2_ci_u = []; f3_ci_l = []
    f3_ci_u = [] ;s_ci_l = []; s_ci_u = []
    with progressbar.ProgressBar(max_value=repl_no) as bar:
        time.sleep(0.1)
        for  tt in range(repl_no):
            time.sleep(0.1)
            bar.update(tt)
            f1=f2=f3=s=0
            for i in range(scenario_number):
                result = scenario_creation(ppos)
                if result == 'F1':
                    f1+=1
                elif result =='F2':
                    f2+=1
                elif result == 'F3':
                    f3+=1
                elif result == 'S':
                    s+=1
            F1+=f1
            F2+=f2
            F3+=f3
            S+=s
            count=[F1,F2,F3,S]
            ci = Mn.multinomial_proportions_confint(count,alpha=0.05,method='goodman')
            f1_ci_l.append(ci[0][0]); f1_ci_u.append(ci[0][1]);
            f2_ci_l.append(ci[1][0]); f2_ci_u.append(ci[1][1]);
            f3_ci_l.append(ci[2][0]); f3_ci_u.append(ci[2][1]);
            s_ci_l.append(ci[3][0]);  s_ci_u.append(ci[3][1]);
            f1list.append(f1)
            f2list.append(f2)
            f3list.append(f3)
            slist.append(s)
            F1list.append(F1/((tt+1)* scenario_number))
            F2list.append(F2/((tt+1)* scenario_number))
            F3list.append(F3/((tt+1)* scenario_number))
            Slist.append(S/((tt+1)* scenario_number))
    if plot:
        plt.subplot(2,2,1)
        plt.tight_layout(pad=0.05, w_pad=0.2, h_pad=0.3)
        plt.plot(np.arange(repl_no),F1list)#,label='$p_1$'
        plt.plot(np.arange(repl_no),f1_ci_l ,'--')
        plt.plot(np.arange(repl_no),f1_ci_u , '--')
        plt.text(-4,F1list[0],"$\hat{P}_{f_1}$")
        plt.text(-4,f1_ci_l[0],"$ci_l$")
        plt.text(-4,f1_ci_u[0],"$ci_u$")
        plt.xlabel('replication')
        plt.ylabel('$f_1$')
        
        plt.subplot(2,2,2)
        plt.tight_layout(pad=0.05, w_pad=0.2, h_pad=0.3)
        plt.plot(np.arange(repl_no),F2list,label='F2')
        plt.plot(np.arange(repl_no),f2_ci_l ,'--')
        plt.plot(np.arange(repl_no),f2_ci_u ,'--')
        plt.text(-4,F2list[0],"$\hat{P}_{f_2}$")
        plt.text(-4,f2_ci_l[0],"$ci_l$")
        plt.text(-4,f2_ci_u[0],"$ci_u$")
        plt.xlabel('replication')
        plt.ylabel('$f_2$')
        
        plt.subplot(2,2,3)
        plt.tight_layout(pad=0.05, w_pad=0.2, h_pad=0.3)
        plt.plot(np.arange(repl_no),F3list)
        plt.plot(np.arange(repl_no),f3_ci_l , '--')
        plt.plot(np.arange(repl_no),f3_ci_u , '--')
        plt.text(-4,F3list[0],"$\hat{P}_{f_3}$")
        plt.text(-4,f3_ci_l[0],"$ci_l$")
        plt.text(-4,f3_ci_u[0],"$ci_u$")
        plt.xlabel('replication')
        plt.ylabel('$f_3$')
        plt.subplot(2,2,4)
        plt.tight_layout(pad=0.05, w_pad=0.2, h_pad=0.3)
        plt.plot(np.arange(repl_no),Slist)
        plt.plot(np.arange(repl_no),s_ci_l , '--')
        plt.plot(np.arange(repl_no),s_ci_u , '--')
        plt.xlabel('iteration')
        plt.ylabel('$S$')
        plt.text(-4,Slist[0],"$\hat{P}_s$")
        plt.text(-4,s_ci_l[0],"$ci_l$")
        plt.text(-4,s_ci_u[0],"$ci_u$") 
        plt.show()
    return F1list[-1],F2list[-1],F3list[-1],Slist[-1]

