# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 10:55:46 2022

@author: hyh
"""
from toolkit import *
import wntr
import time
import csv
import math
import networkx as nx

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


diameterList = [12.0, 16.0, 20.0, 24.0, 30.0, 40.0] #英寸,1英寸(in)=0.0254米(m) 

wn = wntr.network.WaterNetworkModel('networks/HAN.inp') # 用wn格式管径用SI单位m，用EN格式管径用mm
G = wn.get_graph() # directed multigraph    
node_degree = G.degree()
juncDegree = [i[1] for i in node_degree if i[0] != '1']                     

nXF = pd.DataFrame()

for num in range(0,30): 
    startTime = time.perf_counter()    
    np.random.seed(num) 
    seed = np.random.randint(1, 100) # 双重控制随机，保证可以复现实验

    tk = ENepanet()  # 创建EPANET DLL对象tk
    tk.ENopen('networks/HAN.inp', '', '')

    class HANProblem(ElementwiseProblem):
        def __init__(self):   
            super().__init__(n_var=34,
                             n_obj=2,
                             n_ieq_constr=1,
                             xl=np.array(np.zeros(34)), 
                             xu=np.array(np.full(34,6)))    
    
        def _evaluate(self, x, out, *args, **kwargs):
    
            totalCost = 0
            pressureViolation = 0
            minPressure = 30.0
            N = 0
            Nmax = 0
            
            juncDemand = list(np.ones(31))
            juncHead = list(np.ones(31))
            juncPressure = list(np.ones(31))
            minHead = list(np.ones(31))
            reservoirDischarge = list(np.ones(1))
            reservoirHead = list(np.ones(1))
        
            
            for i in range(0,34): # HAN有34条管道,范围[0,34)=[0,33]
                dCode = int(np.floor(x[i]))
                pipeDiameter = 25.4*diameterList[dCode] #in. => mm
                tk.ENsetlinkvalue(i+1, 0, pipeDiameter) #EN_DIAMETER为0
                pipeLength = tk.ENgetlinkvalue(i+1, 1) #EN_LENGTH为1
                unitCost = diameterList[dCode]**1.5 #exponent=1.5
                totalCost += 1.1*unitCost*pipeLength
                
            errorCode = tk.ENsolveH()
            if not errorCode:
                errorCode = 0
            
            for i in range(0,31):  #有31个Junction
                juncDemand[i] = tk.ENgetnodevalue(i+1, 9)#/3600    # EN_DEMAND 9 Actual demand 
                juncHead[i] = tk.ENgetnodevalue(i+1, 10)     # EN_HEAD 10 Hydraulic head 
                juncPressure[i] = tk.ENgetnodevalue(i+1, 11) # EN_PRESSURE 11 Pressure 
                if juncPressure[i]<minPressure:
                    pressureViolation += (minPressure-juncPressure[i])                                               
                juncElevation = tk.ENgetnodevalue(i+1, 0) #EN_ELEVATION 0 Elevation 
                minHead[i] = juncElevation + minPressure
            
            for i in range(0,1): 
                reservoirDischarge[i] = tk.ENgetnodevalue(31+1+i, 9) #31Junc
                reservoirHead[i] = tk.ENgetnodevalue(31+1+i, 10)
            
            # Calculate diameter uniformity coefficient  (Revised Version)20231007
            fromNode = list(np.ones(34))
            toNode = list(np.ones(34))            
            for i in range(0,34): 
                fromNode[i] = tk.ENgetlinknodes(i+1)[0]
                toNode[i] = tk.ENgetlinknodes(i+1)[1]
            

            
            juncUnif = list(np.ones(31))
            juncNeighbor = []
            for i in range(1,32): 
                nPipe = 0
                maxDiameter = 0 
                sumDiameter = 0
                pipeIndex = list(np.ones(34))
                tempNei = []
                for j in range(0,34): 
                    if fromNode[j] == i and tk.ENgetlinkvalue(j+1, 0) > 0.01:
                        nPipe += 1
                        pipeIndex[nPipe-1] = 1+j
                        tempNei.append(toNode[j])
                    if toNode[j] == i and tk.ENgetlinkvalue(j+1, 0) > 0.01:
                        nPipe += 1
                        pipeIndex[nPipe-1] = 1+j 
                        if fromNode[j] != 32:   tempNei.append(fromNode[j]) #排除水源Index=32
                            
                juncNeighbor.append(tempNei)
                
                for k in range(0,nPipe):
                    kth_diameter = tk.ENgetlinkvalue(pipeIndex[k], 0)
                    sumDiameter += kth_diameter
                    if kth_diameter > maxDiameter:
                        maxDiameter = kth_diameter
                juncUnif[i-1]=(sumDiameter/(nPipe*maxDiameter))**(1/juncDegree[i-1])
            
            #print(juncNeighbor)    
            for i in range(0,31):
                tempList = juncNeighbor[i]
                tempUnif = 0
                for t in tempList:
                    tempUnif += juncUnif[t-1]
                tempUnif /= len(tempList)
                juncUnif[i] *= tempUnif
                    

                N += (juncUnif[i])*juncDemand[i]*(juncHead[i]-minHead[i])
                Nmax -= juncDemand[i]*minHead[i]
            
                #print(juncUnif[i]**juncDegree[i])
            for i in range(0,1):
                Nmax += -1*reservoirDischarge[i]*reservoirHead[i]
            networkResilience = N / Nmax
    
    
    
            f1 = totalCost/1000000.0   #$ => $ MM
            f2 = -networkResilience
            g1 = pressureViolation + errorCode 
    
            out["F"] = [f1, f2]
            out["G"] = [g1]
                
    problem = HANProblem()
    
    algorithm = NSGA2(
        pop_size=100,
        n_offsprings=100, #每一代的个体数
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=15), 
        mutation=PM(prob=1/34, eta=7), 
        eliminate_duplicates=True)

    termination = get_termination("n_gen", 500) 
        
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
    with open('OptRes/HAN_NRI2uk_RunTime.csv', 'a', newline='') as file:
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
    
    #导出决策变量和目标函数
    XF.to_csv(f"OptRes/HAN_NRI2uk_seed{num}.csv")           
    
    # Visualize the Pareto Front
    plt.figure(figsize=(6, 6))
    
    s = plt.scatter(F[:, 0], -F[:, 1], s=30,
                    facecolors='none', edgecolors='lightcoral',label='HAN')
    plt.legend(loc="lower right", fontsize=13)
    
    plt.xlabel("Cost(M$)", fontsize=16)
    plt.ylabel("NRI2uk(-)", fontsize=16)
    
    plt.xlim(6, 11.5)
    plt.xticks(np.linspace(6, 12, 7, endpoint=True), fontsize=16)
    # plt.ylim(0.1, 0.4)
    # plt.yticks(np.linspace(0.1, 0.4, 7, endpoint=True), fontsize=16)
    
    plt.title("Cost VS NRI2uk", fontsize=18)
    plt.show()
    
    nXF = pd.concat([nXF,XF],axis=0) #axis=1沿横轴拼接，=0沿纵轴拼接

nXF.to_csv("opt30ResHAN/HAN9_NRI2uk_opt30.csv")   #导出决策变量和目标函数    
