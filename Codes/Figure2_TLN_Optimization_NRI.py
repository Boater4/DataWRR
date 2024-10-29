# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 10:55:46 2022

@author: hyh
"""
from toolkit import *
import wntr
# from wntr.epanet.toolkit import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.termination import get_termination
from pymoo.optimize import minimize


diameterList = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0] # 1inch=0.0254m
unitCostList = [2.0, 5.0, 8.0 ,11.0, 16.0, 23.0, 32.0, 50.0, 60.0, 90.0, 130.0, 170.0, 300.0, 550.0]

nXF = pd.DataFrame()

nRuns = 1 
nGens = 250 
for num in range(0,nRuns):  
    
    np.random.seed(num) 
    seed = np.random.randint(1, 100) 

    tk = ENepanet()  
    tk.ENopen('networks/TLN.inp', '', '')

    class TLNProblem(ElementwiseProblem):
        def __init__(self):   
            super().__init__(n_var=8,
                             n_obj=2,
                             n_ieq_constr=1,
                             xl=np.array(np.zeros(8)), 
                             xu=np.array(np.full(8,14))) # 14 diameters   
    
        def _evaluate(self, x, out, *args, **kwargs):
    
            totalCost = 0
            pressureViolation = 0
            minPressure = 30.0
            N = 0
            Nmax = 0
            
            juncDemand = list(np.ones(6))
            juncHead = list(np.ones(6))
            juncPressure = list(np.ones(6))
            minHead = list(np.ones(6))
            reservoirDischarge = list(np.ones(1))
            reservoirHead = list(np.ones(1))
        
            
            for i in range(0,8): # 8 Pipes
                dCode = int(np.floor(x[i]))
                pipeDiameter = 25.4*diameterList[dCode] #in. => mm
                tk.ENsetlinkvalue(i+1, 0, pipeDiameter) #EN_DIAMETER=0
                pipeLength = tk.ENgetlinkvalue(i+1, 1) #EN_LENGTH=1
                unitCost = unitCostList[dCode]
                totalCost += unitCost*pipeLength
                
            errorCode = tk.ENsolveH()
            if not errorCode:
                errorCode = 0
            
            for i in range(0,6):  # 6 Junctions
                juncDemand[i] = tk.ENgetnodevalue(i+1, 9)# (Actual demand)EN_DEMAND=9  
                juncHead[i] = tk.ENgetnodevalue(i+1, 10)     # EN_HEAD=10
                juncPressure[i] = tk.ENgetnodevalue(i+1, 11) # EN_PRESSURE=11
                if juncPressure[i]<minPressure:
                    pressureViolation += (minPressure-juncPressure[i])                                               
                juncElevation = tk.ENgetnodevalue(i+1, 0) # EN_ELEVATION=0
                minHead[i] = juncElevation + minPressure
            
            for i in range(0,1): 
                reservoirDischarge[i] = tk.ENgetnodevalue(6+1+i, 9) #6Junc
                reservoirHead[i] = tk.ENgetnodevalue(6+1+i, 10)
            
            # Calculate diameter uniformity coefficient for NRI
            fromNode = list(np.ones(8))
            toNode = list(np.ones(8))            
            for i in range(0,8): 
                fromNode[i] = tk.ENgetlinknodes(i+1)[0]
                toNode[i] = tk.ENgetlinknodes(i+1)[1]
             
            juncUnif = list(np.ones(6))
            for i in range(1,7): 
                nPipe = 0
                maxDiameter = 0 
                sumDiameter = 0
                pipeIndex = list(np.ones(8))
                for j in range(0,8): 
                    if fromNode[j] == i and tk.ENgetlinkvalue(j+1, 0) > 0.01:
                        nPipe += 1
                        pipeIndex[nPipe-1] = 1+j
                    if toNode[j] == i and tk.ENgetlinkvalue(j+1, 0) > 0.01:
                        nPipe += 1
                        pipeIndex[nPipe-1] = 1+j 
                for k in range(0,nPipe):
                    kth_diameter = tk.ENgetlinkvalue(pipeIndex[k], 0)
                    sumDiameter += kth_diameter
                    if kth_diameter > maxDiameter:
                        maxDiameter = kth_diameter
                juncUnif[i-1]=sumDiameter/(nPipe*maxDiameter)
                    
            for i in range(0,6):
                N += juncUnif[i]*juncDemand[i]*(juncHead[i]-minHead[i])
                Nmax -= juncDemand[i]*minHead[i]
            
            for i in range(0,1):
                Nmax += -1*reservoirDischarge[i]*reservoirHead[i]
            networkResilience = N / Nmax
    
    
    
            f1 = totalCost/1000000.0   #$ => $ MM
            f2 = -networkResilience
            g1 = pressureViolation + errorCode 
    
            out["F"] = [f1, f2]
            out["G"] = [g1]
                
    problem = TLNProblem()
    
    algorithm = NSGA2(
        pop_size=100,
        n_offsprings=100, 
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=15), 
        mutation=PM(prob=1/8, eta=7), 
        eliminate_duplicates=True)

    termination = get_termination("n_gen", nGens) 
        
    res = minimize(problem,
                   algorithm,
                   termination,
                   seed,
                   save_history=True,
                   verbose=True)
    
    tk.ENclose()
    
    X = res.X
    F = res.F
    
    # Visualize the Pareto Front
    plt.figure(figsize=(6, 6))
    
    s = plt.scatter(F[:, 0], -F[:, 1], s=30,
                    facecolors='none', edgecolors='lightcoral',label='TLN')
    plt.legend(loc="lower right", fontsize=13)
    
    plt.xlabel("Cost(M$)", fontsize=16)
    plt.ylabel("NRI(-)", fontsize=16)
    
    plt.xlim(0, 4.5)
    plt.xticks(np.linspace(0, 4.5, 10, endpoint=True), fontsize=16)
    plt.ylim(0.1, 1.0)
    plt.yticks(np.linspace(0.1, 1.0, 10, endpoint=True), fontsize=16)
    
    plt.title("Cost VS NRI", fontsize=18)
    plt.show()
    
    X_res = pd.DataFrame(X)
    X_res = X_res.applymap(lambda x: int(np.floor(x)))
    F_res = pd.DataFrame(F)
    F_res = F_res.apply(lambda x:abs(x))
    F_res.columns = ['Cost','NRI']
    XF = pd.concat([X_res,F_res],axis=1)
    # XF.to_csv(f"TLN_NRI_run{num}_gen{nGens}.csv") 
    nXF = pd.concat([nXF,XF],axis=0) 
    nXF = nXF.drop_duplicates(keep=False)

nXF.to_excel(f"TLN_NRI_run{nRuns}_gen{nGens}.xlsx", index=False)    
