# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 20:51:18 2024

@author: 10249
"""

import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

path = '******'
os.chdir(path)

cases = ['HAN','FOS','PES','MOD']
labels = ['NRI', 'NRI$_{2u}$','NRI$_{2uk}$']
labels2 = ['','2u','2uk']
colors = ['darkorange','deeppink','green']
markers = ['o','D','s']
plots = {'HAN':221,'FOS':222,'PES':223,'MOD':224}
 
costDict = {}
for case in cases[:]:
    if case == 'PES' or case == 'MOD': rsmlist = ['NRI','NRI2u','NRI2uk']
    else: rsmlist = ['RI','NRI','MRI','REDU','API','PHRI','WRI','NRI2u','NRI2uk']   
    minCost = 0
    maxCost = 9999
    for rsm in rsmlist[:]:
        df = pd.read_excel(f'paper2_{case}_MRSres/compare/{rsm}_MRSres.xlsx',index_col=0)
        if min(df['Cost']) > minCost:
            minCost = min(df['Cost'])
        if max(df['Cost']) < maxCost:
            maxCost = max(df['Cost'])   
    ΔCost = (maxCost-minCost)/3/4 #取前1/3公共成本区间

    cost5 = []
    for i in range(0,5):
        cost1 = float(f'{minCost+i*ΔCost:.2f}')
        cost5.append(cost1)
        
    costDict[case] = cost5

fig = plt.figure(figsize=(12,12), dpi=600)
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1,wspace=0.3,hspace=0.3)
 
for case in cases[:]:
    r = 0
    rsmlist = ['NRI','NRI2u','NRI2uk']       
    ax1 = plt.subplot(plots[case]) 
    if case == 'HAN':    ax1.set_xlabel("Cost (M$)",fontsize=20)
    else:    ax1.set_xlabel("Cost (M€)",fontsize=20)
    ax1.set_ylabel("Mechanical reliability (-)",fontsize=20)
   
    res = {}
    sol_id = {'NRI':{},'NRI2u':{},'NRI2uk':{}}
    for rsm in rsmlist[:]:    
        df = pd.read_excel(f'paper2_{case}_MRSres/compare/{rsm}_MRSres.xlsx',index_col=0)
        temp_res = {}
        avgMRS = 0
        for i in range(0,5):     
            #cost = minCost*(1+i*0.05) #等成本分成n%档
            cost = costDict[case][i]
            costID = (df.iloc[:,0] - cost).abs().idxmin()

            MRS = df.iloc[costID,-1]
            temp_res[f'MRS{i+1}_{cost:.2f}M'] = MRS
            sol_id[rsm][f'Cost{i+1}'] = f'Sol{costID+1},MRS{MRS:.3f}' #输出ID从1开始            
            avgMRS += MRS
        avgMRS /= 5
        temp_res['X_EqualCost'] = avgMRS
        
        # 按照"Cost"列进行排序
        df_sorted = df.sort_values(by='Cost')

        # 截取"Cost"在决策成本区间内的数据片段
        df_part = df_sorted[(df_sorted['Cost'] >= costDict[case][0]) & (df_sorted['Cost'] <= costDict[case][4])]
        
        costList = df_part.iloc[:,0]
        sinMRS = df_part.iloc[:,1]
        cmulMRS = df_part.iloc[:,2]
        oMRS = df_part.iloc[:,3]
        
        sinMRS1 = []
        costList1 = []
        for i in sinMRS:
            if i < 1: 
                sinMRS1.append(i)
                costList1.append(i)
        spcor0 = stats.spearmanr(costList1, sinMRS1)[0]
        spcor1 = stats.spearmanr(costList, cmulMRS)[0]    
        spcor = (spcor0 + spcor1) / 2
        
        temp_res['SpCor1'] = spcor0       
        temp_res['SpCor2'] = spcor1
        temp_res['Y_Consistency'] = spcor
        
        res[rsm] = temp_res
        
        # 作图
        ax1.plot(costList,oMRS,color=colors[r],
                 label=f"{labels[r]}  ε={avgMRS:.3f}  ρ={spcor:.3f}",
                 linewidth=1.5)#marker=markers[i],markersize=4
        
        lowCost = costDict[case][0]-0.2*(costDict[case][1]-costDict[case][0])
        upCost = costDict[case][4]+0.2*(costDict[case][4]-costDict[case][3])
        plt.xlim(lowCost,upCost)
        plt.xticks(costDict[case],size=16)
        if case == 'HAN':
            plt.ylim(0.61,0.73)
            plt.yticks(np.linspace(0.62,0.72,6,endpoint=True),size=16)
            plt.title("(a)",fontsize=18,x=-0.15,y=1)   
        elif case == 'FOS':
            plt.ylim(0.55,1.05)
            plt.yticks(np.linspace(0.6,1,5,endpoint=True),size=16)
            plt.title("(b)",fontsize=18,x=-0.15,y=1)               
        elif case == 'PES':
            plt.ylim(0.63,0.97)
            plt.yticks(np.linspace(0.65,0.95,7,endpoint=True),size=16)
            plt.title("(c)",fontsize=18,x=-0.15,y=1)  
        elif case == 'MOD':
            plt.ylim(0.67,0.79)
            plt.yticks(np.linspace(0.68,0.78,6,endpoint=True),size=16)
            plt.title("(d)",fontsize=18,x=-0.15,y=1)               
        
        
        plt.legend(loc="lower right",ncol=1,frameon=True,fontsize=14)
        plt.grid(axis="x",c=(217/256,217/256,217/256),alpha=0.8)   #设置网格线,c为颜色
        plt.grid(axis="y",c=(217/256,217/256,217/256),alpha=0.8)   #设置网格线,c为颜色
        
        r += 1
    df_res = pd.DataFrame(res).T
    
    # 以NRI为基准计算ΔECP(%)和ΔTCP(%)
    df_res['ΔX(%)'] = ((df_res['X_EqualCost']-df_res.loc['NRI','X_EqualCost'])/df_res.loc['NRI','X_EqualCost']) * 100
    df_res['ΔY(%)'] = ((df_res['Y_Consistency']-df_res.loc['NRI','Y_Consistency'])/df_res.loc['NRI','Y_Consistency']) * 100

    # 计算ΔECP+ΔTCP(%)
    df_res['ΔX+Y(%)'] = df_res['ΔX(%)'] + df_res['ΔY(%)']
    
    df_res.to_excel(f'paper2_{case}_res/{case}_RSM_ComRes.xlsx')
    df_id = pd.DataFrame(sol_id).T #输出ID从1开始，与帕累托文件行数直接对应
    df_id.to_excel(f'figsManu2/EqualCostID/{case}_EqualSolutionID.xlsx')

plt.savefig("******.svg",dpi=600,bbox_inches='tight')