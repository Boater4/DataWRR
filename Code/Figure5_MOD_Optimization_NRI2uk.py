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
maxPressureList = [35.007,34.875,35.797,37.254,38.235,38.545,38.545,38.413,36.321,
                   37.497,38.000,37.112,36.426,37.481,33.243,35.150,34.971,37.906,
                   37.739,36.785,37.188,36.877,37.513,39.295,39.387,39.846,40.175,
                   38.355,38.204,38.403,38.361,38.700,41.239,41.163,40.987,41.800,
                   41.853,41.935,40.935,42.905,43.119,41.833,41.001,40.929,40.726,
                   40.363,40.820,40.794,42.823,41.155,41.668,41.722,35.224,37.377,
                   38.016,38.084,38.365,38.451,37.735,39.016,39.451,39.395,36.549,
                   36.058,36.693,36.282,35.773,35.547,34.799,33.911,33.688,33.436,
                   33.047,32.670,33.065,33.408,33.757,35.895,37.585,37.751,37.687,
                   37.455,38.617,38.046,38.339,39.509,38.888,39.608,38.914,38.800,
                   39.305,38.860,38.571,36.861,37.332,37.395,37.529,37.503,37.761,
                   39.724,40.243,40.840,40.716,40.754,41.123,39.650,40.227,40.203,
                   40.546,40.580,42.183,39.742,40.287,39.576,38.544,43.811,43.905,
                   43.769,43.797,43.480,43.468,42.755,42.500,42.452,42.402,40.740,
                   42.229,42.640,42.083,41.498,40.874,38.134,38.806,38.976,38.940,
                   38.583,39.133,39.443,40.375,35.150,35.396,34.659,34.659,35.051,
                   34.795,36.549,36.890,36.549,38.814,39.183,38.690,38.688,38.481,
                   36.246,36.996,36.964,37.421,37.745,38.615,38.732,39.796,39.131,
                   39.507,38.573,38.235,41.833,41.746,41.616,40.415,38.407,38.451,
                   38.459,38.483,42.743,42.590,42.701,43.017,43.384,43.404,43.306,
                   44.108,43.953,43.366,42.690,42.155,41.674,40.806,41.325,41.271,
                   41.157,40.728,40.732,42.296,40.095,41.111,40.155,39.473,40.061,
                   39.966,39.565,39.796,37.800,38.297,39.469,37.735,38.303,36.621,
                   36.465,37.637,37.262,37.842,38.010,37.200,34.201,34.651,33.502,
                   33.340,39.451,40.580,42.356,40.333,39.403,42.951,42.755,42.434,
                   42.556,42.843,43.460,43.450,36.008,38.816,39.110,39.612,39.642,
                   39.505,41.959,40.087,38.343,39.195,39.329,41.582,41.434,42.590,
                   42.498,42.452,42.446,43.795,43.168,38.204,38.669,37.555,36.487,
                   37.850,37.595,37.727,43.003,35.849,34.957,34.919,34.919,33.949,
                   33.714,33.546,36.745,38.537,37.691,38.289,38.888]

wn = wntr.network.WaterNetworkModel('networks/MOD.inp')
G = wn.get_graph() # directed multigraph    
node_degree = G.degree()
resList = ['269','270','271','272']
juncDegree = [i[1] for i in node_degree if i[0] not in resList] # Reservoir node    

nXF = pd.DataFrame()

for num in range(0,30): 
    startTime = time.perf_counter()
    np.random.seed(num) 
    seed = np.random.randint(1, 100) # 双重控制随机，保证可以复现实验

    tk = ENepanet()  # Create EPANET DLL object "tk"
    tk.ENopen('networks/MOD.inp', '', '')
    
    class MODProblem(ElementwiseProblem):
        def __init__(self):   
            super().__init__(n_var=317,   # MOD has 317 pipes
                             n_obj=2,
                             n_ieq_constr=1,
                             xl=np.array(np.zeros(317)),  
                             xu=np.array(np.full(317,13))) #MOD has 13 diameters 
    
        def _evaluate(self, x, out, *args, **kwargs):
    
            totalCost = 0
            pressureViolation = 0
            velocityViolation = 0
            minPressure = 20.0
            N = 0
            Nmax = 0
            
            juncDemand = list(np.ones(268))
            juncHead = list(np.ones(268))
            juncPressure = list(np.ones(268))
            minHead = list(np.ones(268))
            reservoirDischarge = list(np.ones(4))
            reservoirHead = list(np.ones(4))
            
            for i in range(0,317): # MOD has 317 pipes
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
            
            for i in range(0,268):  #MOD has 268 junctions.
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
            
            for i in range(0,4): # MOD has four reservoirs
                reservoirDischarge[i] = tk.ENgetnodevalue(268+1+i, 9) #268Junc
                reservoirHead[i] = tk.ENgetnodevalue(268+1+i, 10)
            
            for i in range(0,317):
                velocity = tk.ENgetlinkvalue(1+i,9) #EN_VELOCITY 9 Flow velocity  
                if velocity> 2.0:
                    velocityViolation += (velocity- 2.0)
                    
            # for i in range(0,317): 
            #     fromNode[i],toNode[i] = tk.ENgetlinknodes(1+i)
            
            # wn = wntr.network.WaterNetworkModel('networks/MOD.inp')
            
            
            # Calculate diameter uniformity coefficient  (Revised Version)20231007
            fromNode = list(np.ones(317))
            toNode = list(np.ones(317))     
            for i in range(0,317): 
                fromNode[i] = tk.ENgetlinknodes(i+1)[0]
                toNode[i] = tk.ENgetlinknodes(i+1)[1]
            
            # print(fromNode,toNode)
            
            juncUnif = list(np.ones(268))
            juncNeighbor = []
            for i in range(1,269): 
                nPipe = 0
                maxDiameter = 0 
                sumDiameter = 0
                pipeIndex = list(np.ones(317)) 
                tempNei = []
                for j in range(0,317): 
                    if fromNode[j] == i and tk.ENgetlinkvalue(j+1, 0) > 0.01:
                        nPipe += 1
                        pipeIndex[nPipe-1] = 1+j
                        tempNei.append(toNode[j])                           
                    if toNode[j] == i and tk.ENgetlinkvalue(j+1, 0) > 0.01:
                        nPipe += 1
                        pipeIndex[nPipe-1] = 1+j
                        slist = [269,270,271,272]
                        if fromNode[j] not in slist:   tempNei.append(fromNode[j]) #排除水源Index=37
                            
                juncNeighbor.append(tempNei)                
                
                        
                for k in range(0,nPipe):
                    kth_diameter = tk.ENgetlinkvalue(pipeIndex[k], 0)
                    sumDiameter += kth_diameter
                    if kth_diameter > maxDiameter:
                        maxDiameter = kth_diameter
                juncUnif[i-1]=(sumDiameter/(nPipe*maxDiameter))**(1/juncDegree[i-1])
                        
            
            
            for i in range(0,268):
                tempList = juncNeighbor[i]
                tempUnif = 0
                for t in tempList:
                    tempUnif += juncUnif[t-1]
                if tempList:
                    tempUnif /= len(tempList)
                    juncUnif[i] *= tempUnif            
            
                N += juncUnif[i]*juncDemand[i]*(juncHead[i]-minHead[i])
                Nmax -= juncDemand[i]*minHead[i]
            
            
            for i in range(0,4):
                Nmax += -1*reservoirDischarge[i]*reservoirHead[i]
            NRI2uk = N / Nmax
    
            f1 = totalCost/1000000.0   #$ => $ MM
            f2 = -NRI2uk
            g1 = pressureViolation + velocityViolation + errorCode 
    
            out["F"] = [f1, f2]
            out["G"] = [g1]
                
    problem = MODProblem()
    
    algorithm = NSGA2(
        pop_size=200,
        n_offsprings=200, #每一代的个体数
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=15), 
        mutation=PM(prob=1/317, eta=7), 
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
    with open('OptRes/MOD_NRI2uk_RunTime.csv', 'a', newline='') as file:
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
    XF.to_csv(f"OptRes/MOD_NRI2uk_seed{num}.csv")   #导出决策变量和目标函数
    
    #作出本实验得到的帕累托前锋                                
    plt.figure(figsize=(6, 6))
 
    s = plt.scatter(F[:, 0], -F[:, 1], s=30, facecolors='none', edgecolors='lightcoral',label='MOD')
    plt.legend(loc="lower right",fontsize=13)
    
    plt.xlabel("Cost(M€)", fontsize=16)
    plt.ylabel("NRI2uk(-)", fontsize=16)
    
    plt.xlim(0,25)
    plt.xticks(np.linspace(0,25,6,endpoint=True))
    plt.ylim(0.3,1.0)
    plt.yticks(np.linspace(0.3,1.0,8,endpoint=True))
    
    plt.title("Cost VS NRI2uk",fontsize=18)
    plt.show()

    nXF = pd.concat([nXF,XF],axis=0) #axis=1沿横轴拼接，=0沿纵轴拼接

nXF.to_csv("opt30ResMOD/MOD_NRI2uk_opt30.csv")   #导出决策变量和目标函数  
















