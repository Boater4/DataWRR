# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 13:12:35 2022

@author: hyh
"""

import wntr
import wntr.network.controls as controls
import numpy as np
import pandas as pd
import os
import scipy.stats as stats
from itertools import combinations
 
path = '******'
os.chdir(path)

diameterList = [100.0,125.0,150.0,200.0,250.0,300.0,350.0,
                400.0,450.0,500.0,600.0,700.0,800.0]# mm 

files = os.listdir('ndSortResMOD')


for j in files[:]:
    rsm = j.split("_")[1]
    df = pd.read_excel(f"ndSortResMOD/{j}",header=None) #输入xlsx无表头
    print(len(df))
    solution = df.iloc[:,0:317].astype(float).values
    cost = df.iloc[:,317].astype(float).values
    
    sinMRS = []
    cmulMRS = []
    overMRS = []
    
    for i in range(0,len(solution)): 
        wn = wntr.network.WaterNetworkModel("networks/MOD.inp")
        
        solution_i = solution[i]
        diameter_list = [diameterList[int(j)] for j in solution_i]    
    
        k = 0                
        for pipe_name, pipe in wn.pipes():
            pipe.diameter = diameter_list[k]*0.001 #wntr的管径用m作单位
            k += 1
        # wn.write_inpfile(f'all_sol/MOD_AllSol_MRS/MOD_{metric}_Sol{i}.inp', version=2.2)
        #Define the duration and timestep of the delay simulation 
        wn.options.time.duration = 0
        
        # Adopt PDD mode
        wn.options.hydraulic.demand_model = 'PDD'  
        
        wn.options.hydraulic.required_pressure = 20 
        wn.options.hydraulic.minimum_pressure  = 0 
        wn.options.hydraulic.pressure_exponent = 0.5 # default = 0.5 
        
        
        # Single-source network does not consider disconnecting water source pipe 
        MRScore = 0
        cp = 0
        for pipe_name in wn.pipe_name_list[:]: #
            wn.reset_initial_values()
            # StrainMetrics = {}
        
            pipe = wn.get_link(pipe_name)        
            
            #Define the action and condition to create a control
            act = controls.ControlAction(pipe, 'status', 0) # 0 closed
            
            T = '0:00:00' # pipe failure (shutdown)
            cond = controls.SimTimeCondition(wn, '=', T)    
            ctrl = controls.Control(cond, act, name='control0')
        
            # Add controls to close the specific pipe
            wn.add_control('close pipe ' + pipe_name, ctrl)
            
            sim = wntr.sim.EpanetSimulator(wn) # 3.27 EpanetSim
            results = sim.run_sim()
            
            # Generate scores for each node
            pressure = results.node['pressure'].loc[:,wn.junction_name_list]
            score = pressure.copy(deep=True)
            p_min = 20
            for j,p in pressure.iloc[0].items():
                if p <= 0:
                    score.loc[:,j] = 0
                elif 0 < p < p_min:
                    score.loc[:,j] = p/p_min
                elif p > p_min:
                    score.loc[:,j] = 1
        
            # here use Availability to express performance
            demand = results.node['demand'].loc[:,wn.junction_name_list]#.sum(axis=1)
            req_demand = wntr.metrics.expected_demand(wn)
            tot_demand = float(req_demand.sum(axis=1))
            junc_w = req_demand / tot_demand
            # tot_w = junc_w.sum(axis=1)
            
            tot_score = junc_w * score
            TotalScore = float(tot_score.sum(axis=1))
            
            MRScore += TotalScore 
            cp += 1
            # Remove the control
            wn.remove_control('close pipe ' + pipe_name)
        
        MRScore /= cp
        sinMRS.append(MRScore)        
        sMRS = MRScore

        # Calculate critical pipe failure scenarios
        MRScore = 0
        cp = 0        
        keyScene = [('330','331'),('330','335'),('330','336'),
                    ('331','335'),('331','336'),('335','336'),
                    ('331','335','336'),
                    ('330','335','336'), 
                    ('330','331','336'), 
                    ('330','331','335')]
        for FailPipeSet in keyScene:
    
            wn.reset_initial_values()
            T = '0:00:00' # pipe failure (shutdown)
            
            # j = 0
            ctrls = locals() #利用命名空间动态定义变量
            for pipe_name in FailPipeSet: 
                pipe = wn.get_link(pipe_name)           
                #Define the action and condition to create a control
                act = controls.ControlAction(pipe, 'status', 0) # 0 closed
                cond = controls.SimTimeCondition(wn, '=', T)    
                ctrls['ctrl' + str(j)] = controls.Control(cond, act, name=f'close{j}')
                # Add controls to close the specific pipe
                wn.add_control('close pipe ' + pipe_name, ctrls['ctrl' + str(j)])
                # j += 1
            
            sim = wntr.sim.EpanetSimulator(wn) # 3.27 ESim
            results = sim.run_sim()
            # Generate scores for each node
            pressure = results.node['pressure'].loc[:,wn.junction_name_list]
            
            # if len(pressure) == 0:
            #     for pipe_name in FailPipeSet: # Remove the control   
            #         wn.remove_control('close pipe ' + pipe_name)
        
            score = pressure.copy(deep=True)
            p_min = 20
            for j,p in pressure.iloc[0].items():
                if p <= 0:
                    score.loc[:,j] = 0
                elif 0 < p < p_min:
                    score.loc[:,j] = p/p_min
                elif p > p_min:
                    score.loc[:,j] = 1
        
            # here use Availability to express performance
            demand = results.node['demand'].loc[:,wn.junction_name_list]#.sum(axis=1)
            req_demand = wntr.metrics.expected_demand(wn)
            tot_demand = float(req_demand.sum(axis=1))
            junc_w = req_demand / tot_demand
            # tot_w = junc_w.sum(axis=1)
            
            tot_score = junc_w * score
            TotalScore = float(tot_score.sum(axis=1))
            
            MRScore += TotalScore 
            cp += 1
            # Remove the control
            for pipe_name in FailPipeSet: 
                wn.remove_control('close pipe ' + pipe_name)
        
        
        MRScore /= cp
        cmulMRS.append(MRScore)
        cMRS = MRScore
      
        oMRS = np.average(np.array([sMRS,cMRS]))
        overMRS.append(oMRS)
        print(f'{rsm} of {i+1} times.\nsMRS={sMRS:2f},cMRS={cMRS:2f},oMRS={oMRS:2f}.')
    
    
    sinMRS = np.array(sinMRS)
    cmulMRS = np.array(cmulMRS)
    overMRS = np.array(overMRS)
    
         
    df = pd.DataFrame({'Cost':cost, 'sMRS':sinMRS, 'cMRS':cmulMRS, 'oMRS':overMRS})
    df.index = range(0, len(df))
    df.to_excel(f'paper2_MOD_MRSres/{rsm}_MRSres.xlsx')
