# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 21:34:15 2024

@author: 10249
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
 
path = '******'
files = os.listdir(path)

rsmlist = ['NRI','NRI2u','NRI2uk','NRI']
colors = ['darkorange','deeppink','green','darkorange']
markers = ['o','D','s','o']
labels = ['NRI', 'NRI$_{2u}$','NRI$_{2uk}$','NRI']
costs = [1.82,3.06,4.31,5.55,6.79]


fig = plt.figure(figsize=(6,6), dpi=600)
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

ax1 = plt.subplot(111)    
ax1.set_xlabel("Cost (M€)",fontsize=24)
ax1.set_ylabel("Surrogate measure (-)",fontsize=24)

for i in range(0,5):
    ax1.axvline(costs[i],ymin=0,ymax=1,c='black',ls='dashed',linewidth=1.2)

plt.axvspan(xmin=costs[0],
            xmax=costs[4],
            facecolor='gray',
            alpha=0.2)  # 绘制垂直于x轴的参考区域

for i in range(0,3):    
    df = pd.read_excel(f'{path}/PES_{rsmlist[i]}_ndSortRes.xlsx',header=None).iloc[:,[-4,-3]].astype(float).values           
    ax1.scatter(df[:, 0], -df[:, 1], s=60, c='none',alpha=0.8, 
                marker=markers[i], edgecolors=colors[i],label=labels[i])

ax1.set_xlim(0.5,19.5)
ax1.set_xticks(np.linspace(2,18,5,endpoint=True))
ax1.set_ylim(0.1,1.1)
ax1.set_yticks(np.linspace(0.2,1,5,endpoint=True))
plt.tick_params(labelsize=22)


ax1.legend(loc="lower right", fontsize=20, frameon=True,ncol=1)
plt.grid(axis="x",c=(217/256,217/256,217/256),alpha=0.8)   #设置网格线,c为颜色
plt.grid(axis="y",c=(217/256,217/256,217/256),alpha=0.8)   #设置网格线,c为颜色

# 补画NRI
for i in range(3,4):    
    df = pd.read_excel(f'{path}/PES_{rsmlist[i]}_ndSortRes.xlsx',header=None).iloc[:,[-4,-3]].astype(float).values           
    ax1.scatter(df[:, 0], -df[:, 1], s=60, c='none',alpha=0.8, 
                marker=markers[i], edgecolors=colors[i],label=labels[i])
 
plt.grid(axis="x",c=(217/256,217/256,217/256),alpha=0.8)   #设置网格线,c为颜色
plt.grid(axis="y",c=(217/256,217/256,217/256),alpha=0.8)   #设置网格线,c为颜色
 
plt.savefig(f"******.svg",dpi=600,bbox_inches='tight')
