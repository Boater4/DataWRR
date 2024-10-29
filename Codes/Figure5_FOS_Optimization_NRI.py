# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 10:55:46 2022

@author: hyh
"""
from toolkit import *
import wntr
import time
import csv
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


diameterList = [16,20.4,26,32.6,40.8,51.4,61.4,73.6,90,102.2,114.6,130.8,147.2,163.6,184,204.6,229.2,257.8,290.6,327.4,368.2,409.2] #英寸,1英寸(in)=0.0254米(m) 
unitCostList = [0.38,0.56,0.88,1.35,2.02,3.21,4.44,6.45,9.59,11.98,14.93,19.61,24.78,30.55,38.71,47.63,59.7,75.61,99.58,126.48,160.29,197.71]
maxPressureList = [55.85,56.6,57.65,58.5,59.76,55.6,53.1,54.5,55,56.83,57.3,58.36,59.1,58.4,57.5,56.7,55.5,56.9,58.1,58.17,58.2,57.1,56.8,53.5,56.6,57.6,57.1,55.35,56.5,56.9,56.6,56.8,56.4,56.3,55.57,55.1]

nXF = pd.DataFrame()

for num in range(0,30): 
    startTime = time.perf_counter()    
    np.random.seed(num) 
    seed = np.random.randint(1, 100) # 双重控制随机，保证可以复现实验

    tk = ENepanet()  # Create EPANET DLL object "tk"
    tk.ENopen('networks/FOS.inp', '', '')
    
    class FOSProblem(ElementwiseProblem):
        def __init__(self):   
            super().__init__(n_var=58,   # FOS has 58 pipes
                             n_obj=2,
                             n_ieq_constr=1,
                             xl=np.array(np.zeros(58)),  
                             xu=np.array(np.full(58,22))) #FOS has 22 diameters 
    
        def _evaluate(self, x, out, *args, **kwargs):
    
            totalCost = 0
            pressureViolation = 0
            velocityViolation = 0
            minPressure = 40.0
            N = 0
            Nmax = 0
            
            juncDemand = list(np.ones(36))
            juncHead = list(np.ones(36))
            juncPressure = list(np.ones(36))
            minHead = list(np.ones(36))
            reservoirDischarge = list(np.ones(1))
            reservoirHead = list(np.ones(1))
            
            for i in range(0,58): # FOS has 58 pipes
                dCode = int(np.floor(x[i]))
                pipeDiameter = diameterList[dCode] 
                # pipeDiameter = 25.4*diameterList[dCode] #in. => mm
                tk.ENsetlinkvalue(i+1, 0, pipeDiameter) #EN_DIAMETER为0
                pipeLength = tk.ENgetlinkvalue(i+1, 1) #EN_LENGTH为1
                unitCost = unitCostList[dCode]
                totalCost += unitCost*pipeLength
    
            
            errorCode = tk.ENsolveH()
            if not errorCode:
                errorCode = 0
            
            for i in range(0,36):  #FOS has 36 junctions.
                juncDemand[i] = tk.ENgetnodevalue(i+1, 9)#/3600    # EN_DEMAND 9 Actual demand 
                juncHead[i] = tk.ENgetnodevalue(i+1, 10)     # EN_HEAD 10 Hydraulic head 
                juncPressure[i] = tk.ENgetnodevalue(i+1, 11) # EN_PRESSURE 11 Pressure 
                
                maxPressure = maxPressureList[i]
                if juncPressure[i] < minPressure:
                    pressureViolation += (minPressure - juncPressure[i])
                elif juncPressure[i] > maxPressure:                                                
                    pressureViolation += (juncPressure[i] - maxPressure)
                juncElevation = tk.ENgetnodevalue(i+1, 0) #EN_ELEVATION 0 Elevation 
                minHead[i] = juncElevation + minPressure
            
            for i in range(0,1): 
                reservoirDischarge[i] = tk.ENgetnodevalue(36+1+i, 9) #36Junc
                reservoirHead[i] = tk.ENgetnodevalue(36+1+i, 10)
            
            for i in range(0,58):
                velocity = tk.ENgetlinkvalue(1+i,9) #EN_VELOCITY 9 Flow velocity  
                if velocity > 1.0:
                    velocityViolation += (velocity - 1.0)
                       
            # Calculate diameter uniformity coefficient  (Revised Version)20231007
            fromNode = list(np.ones(58))
            toNode = list(np.ones(58))     
            for i in range(0,58): 
                fromNode[i] = tk.ENgetlinknodes(i+1)[0]
                toNode[i] = tk.ENgetlinknodes(i+1)[1]
             
            juncUnif = list(np.ones(36))
            for i in range(1,37): 
                nPipe = 0
                maxDiameter = 0 
                sumDiameter = 0
                pipeIndex = list(np.ones(58))   
                for j in range(0,58): 
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
                    
        
            
            for i in range(0,36):
                N += juncUnif[i]*juncDemand[i]*(juncHead[i]-minHead[i])
                Nmax -= juncDemand[i]*minHead[i]
            
            
            for i in range(0,1):
                Nmax += -1*reservoirDischarge[i]*reservoirHead[i]
            NRI = N / Nmax
    
            f1 = totalCost/1000000.0   #$ => $ MM
            f2 = -NRI
            g1 = pressureViolation + velocityViolation + errorCode 
    
            out["F"] = [f1, f2]
            out["G"] = [g1]
                
    problem = FOSProblem()
    
    algorithm = NSGA2(
        pop_size=100,
        n_offsprings=100, #每一代的个体数
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=15), 
        mutation=PM(prob=1/58, eta=7), 
        eliminate_duplicates=True)
    
    termination = get_termination("n_gen", 5000) 

    res = minimize(problem,
                   algorithm,
                   termination,
                   seed,
                   save_history=True,
                   verbose=True)
    
    tk.ENclose()
    
    endTime = time.perf_counter()
    runTime = endTime - startTime
    # 使用 with open 语句写入CSV文件
    with open('OptRes/FOS_NRI_RunTime.csv', 'a', newline='') as file:
        writer = csv.writer(file)
        if num == 0:
            writer.writerow(["Seed", "runTime"])
        writer.writerow([num, f'{runTime:.6f}'])
    
    
    X = res.X
    F = res.F
    
    X_res = pd.DataFrame(X)
    F_res = pd.DataFrame(F)
    F_res[1].apply(lambda x:abs(x))
    XF = pd.concat([X_res,F_res],axis=1)
    XF.to_csv(f"OptRes/FOS_NRI_seed{num}.csv")   #导出决策变量和目标函数
    
    #作出本实验得到的帕累托前锋                                
    plt.figure(figsize=(6, 6))
 
    s = plt.scatter(F[:, 0], -F[:, 1], s=30, facecolors='none', edgecolors='lightcoral',label='FOS')
    plt.legend(loc="lower right",fontsize=13)
    
    plt.xlabel("Cost(M€)", fontsize=16)
    plt.ylabel("NRI(-)", fontsize=16)
    
    plt.xlim(0,1.5)
    plt.xticks(np.linspace(0,1.5,6,endpoint=True))
    # plt.ylim(0.5,1.0)
    # plt.yticks(np.linspace(0.5,1.0,6,endpoint=True))
    plt.yticks(fontsize=16)
    
    plt.title("Cost VS NRI",fontsize=18)
    plt.show()

    nXF = pd.concat([nXF,XF],axis=0) #axis=1沿横轴拼接，=0沿纵轴拼接

nXF.to_csv("opt30ResFOS/FOS1_NRI_opt30.csv")   #导出决策变量和目标函数  
















