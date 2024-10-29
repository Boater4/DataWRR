# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 13:30:00 2023

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


diameterList = [100.0,125.0,150.0,200.0,250.0,300.0,350.0,
                400.0,450.0,500.0,600.0,700.0,800.0]#英寸,1英寸(in)=0.0254米(m) 
unitCostList = [27.7,38,40.5,55.4,75,92.4,123.1,141.9,169.3,191.5,246,319.6,391.1]
maxPressureList = [54.1,52,53.5,53.2,54.8,50,50.5,51.8,52.7,51.8,
                   29,53.8,53.8,37.8,53,53.5,54.2,54.3,55.2,54.9,
                   55,38.3,53.8,55.4,53.2,52.7,53.4,53.8,54.2,55.2,
                   53.2,54.1,53.3,54.9,50.3,50.3,52.8,54.1,54.1,28.5,
                   29.7,55.9,54.9,54.2,53.7,55.2,55.2,55.4,55.2,55.4,
                   53.7,54.5,54.2,53.8,53.4,53.2,53.9,54,52.8,52.8,
                   54.9,54.9,53.7,54.9,55.5,37.8,54.9,55.5]

nXF = pd.DataFrame()

for num in range(0,30): 
    startTime = time.perf_counter()
    np.random.seed(num) 
    seed = np.random.randint(1, 100) # 双重控制随机，保证可以复现实验

    tk = ENepanet()  # Create EPANET DLL object "tk"
    tk.ENopen('networks/PES.inp', '', '')
    
    class PESProblem(ElementwiseProblem):
        def __init__(self):   
            super().__init__(n_var=99,   # PES has 99 pipes
                             n_obj=2,
                             n_ieq_constr=1,
                             xl=np.array(np.zeros(99)),  
                             xu=np.array(np.full(99,13))) #PES has 13 diameters 
    
        def _evaluate(self, x, out, *args, **kwargs):
    
            totalCost = 0
            pressureViolation = 0
            velocityViolation = 0
            minPressure = 20.0
            N = 0
            Nmax = 0
            
            juncDemand = list(np.ones(68))
            juncHead = list(np.ones(68))
            juncPressure = list(np.ones(68))
            minHead = list(np.ones(68))
            reservoirDischarge = list(np.ones(3))
            reservoirHead = list(np.ones(3))
            
            for i in range(0,99): # PES has 99 pipes
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
            
            for i in range(0,68):  #PES has 68 junctions.
                juncDemand[i] = tk.ENgetnodevalue(i+1, 9)# EN_DEMAND 9 Actual demand 
                juncHead[i] = tk.ENgetnodevalue(i+1, 10)     # EN_HEAD 10 Hydraulic head 
                juncPressure[i] = tk.ENgetnodevalue(i+1, 11) # EN_PRESSURE 11 Pressure 
                
                maxPressure = maxPressureList[i]
                if juncPressure[i] < minPressure:
                    pressureViolation += (minPressure - juncPressure[i])
                elif juncPressure[i] > maxPressure:                                                
                    pressureViolation += (juncPressure[i] - maxPressure)
                juncElevation = tk.ENgetnodevalue(i+1, 0) #EN_ELEVATION 0 Elevation 
                minHead[i] = juncElevation + minPressure
            
            for i in range(0,3): # PES has three reservoirs
                reservoirDischarge[i] = tk.ENgetnodevalue(68+1+i, 9) #68Junc
                reservoirHead[i] = tk.ENgetnodevalue(68+1+i, 10)
            
            for i in range(0,99):
                velocity = tk.ENgetlinkvalue(1+i,9) #EN_VELOCITY 9 Flow velocity  
                if velocity> 2.0:
                    velocityViolation += (velocity- 2.0)
                    
            # for i in range(0,99): 
            #     fromNode[i],toNode[i] = tk.ENgetlinknodes(1+i)
            
            # wn = wntr.network.WaterNetworkModel('networks/PES.inp')
            # for i in range(0,len(wn.pipe_name_list)):
                # pipe_name = wn.pipe_name_list[i]
                # pipe = wn.get_link(pipe_name)
                # var = int(x[i]) #这里已经建立了变量X与NRI的关联
                # pipe.diameter = float(diameterList[var]*0.001) #PES与HAN单位不同，但输入WNTR中都是SI单位
            
            # juncU = []
            # for name, junc in wn.junctions():
                # temp1 = []
                # for pipe_name, pipe in wn.pipes():
                    # if pipe.start_node_name == name and pipe.diameter > 0.01:
                        # temp1.append(pipe_name)
                # n_out = len(temp1)
                # temp2 = []
                # for pipe_name, pipe in wn.pipes():
                    # if pipe.end_node_name == name and pipe.diameter > 0.01:
                        # temp2.append(pipe_name)
                # n_in = len(temp2)
                # n_links = n_out + n_in
                
                # D = []
                # for pipe_name, pipe in wn.pipes():
                    # if pipe_name in temp1:
                        # D.append(pipe.diameter)
                    # if pipe_name in temp2:
                        # D.append(pipe.diameter)
                # Dmax = max(D)
                # Dtot = sum(D)
                # Uniformity = Dtot/(n_links*Dmax)
                # juncU.append(Uniformity) # Give the Uniformity of each node in order       
            
            # print(juncDemand)
            
            # Calculate diameter uniformity coefficient  (Revised Version)20231007
            fromNode = list(np.ones(99))
            toNode = list(np.ones(99))     
            for i in range(0,99): 
                fromNode[i] = tk.ENgetlinknodes(i+1)[0]
                toNode[i] = tk.ENgetlinknodes(i+1)[1]
            
            # print(fromNode,toNode)
            
            juncUnif = list(np.ones(68))
            for i in range(1,69): 
                nPipe = 0
                maxDiameter = 0 
                sumDiameter = 0
                pipeIndex = list(np.ones(99))   
                for j in range(0,99): 
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
                        
            
            
            for i in range(0,68):
                N += juncUnif[i]*juncDemand[i]*(juncHead[i]-minHead[i])
                Nmax -= juncDemand[i]*minHead[i]
            
            
            for i in range(0,3):
                Nmax += -1*reservoirDischarge[i]*reservoirHead[i]
            NRI = N / Nmax
    
            f1 = totalCost/1000000.0   #$ => $ MM
            f2 = -NRI
            g1 = pressureViolation + velocityViolation + errorCode 
    
            out["F"] = [f1, f2]
            out["G"] = [g1]
                
    problem = PESProblem()
    
    algorithm = NSGA2(
        pop_size=200,
        n_offsprings=200, #每一代的个体数
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=15), 
        mutation=PM(prob=1/99, eta=7), 
        eliminate_duplicates=True)
    
    termination = get_termination("n_gen", 10000) #1000

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
    with open('OptRes/PES_NRI_RunTime.csv', 'a', newline='') as file:
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
    XF.to_csv(f"OptRes/PES_NRI_seed{num}.csv")   #导出决策变量和目标函数
    
    #作出本实验得到的帕累托前锋                                
    plt.figure(figsize=(6, 6))
 
    s = plt.scatter(F[:, 0], -F[:, 1], s=30, facecolors='none', edgecolors='lightcoral',label='PES')
    plt.legend(loc="lower right",fontsize=13)
    
    plt.xlabel("Cost(M€)", fontsize=16)
    plt.ylabel("NRI(-)", fontsize=16)
    
    plt.xlim(0,20)
    plt.xticks(np.linspace(0,20,9,endpoint=True))
    plt.ylim(0.2,1.0)
    plt.yticks(np.linspace(0.2,1.0,9,endpoint=True))
    
    plt.title("Cost VS NRI",fontsize=18)
    plt.show()

    nXF = pd.concat([nXF,XF],axis=0) #axis=1沿横轴拼接，=0沿纵轴拼接

nXF.to_csv("opt30Res/PES_NRI_opt30.csv")   #导出决策变量和目标函数  
















