# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 20:08:01 2024

@author: 10249
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import wntr
import os

os.chdir('******')

# 管径值
diameters = [100.0,125.0,150.0,200.0,250.0,300.0,350.0,
              400.0,450.0,500.0,600.0,700.0,800.0]# mm 

wn = wntr.network.WaterNetworkModel('PESinp/PES_NRI_Sol362.inp')
diaList = []
for name, pipe in wn.pipes():
    diaList.append(1000*float(f'{pipe.diameter:.3f}'))
NRI = diaList

wn = wntr.network.WaterNetworkModel('PESinp/PES_NRI2u_Sol247.inp')
diaList = []
for name, pipe in wn.pipes():
    diaList.append(1000*float(f'{pipe.diameter:.3f}'))
NRI2u = diaList

wn = wntr.network.WaterNetworkModel('PESinp/PES_NRI2uk_Sol381.inp')
diaList = []
for name, pipe in wn.pipes():
    diaList.append(1000*float(f'{pipe.diameter:.3f}'))
NRI2uk = diaList

# for i in range(99):
#     NRI[i] = diameters[NRI[i]]
#     NRI2u[i] = diameters[NRI2u[i]]
#     NRI2uk[i] = diameters[NRI2uk[i]]

def calculate_offsets(data):
    offset_data = []
    # data = sorted(data)
    interval = 0.2
    # 使用字典将相同数字归为一组
    groups = defaultdict(list)
    for num in data:
        groups[num].append(num)
    for num, group in groups.items():
        count = len(group)
        if count == 1:
            offset_data.append(0)
        elif count > 1:
            templist = list(np.linspace(-interval*(count-1)/2,interval*(count-1)/2, count,endpoint=True))
            for temp in templist:
                offset_data.append(temp)        
    return offset_data

# 指标对应位置
positions = [0, 8, 16]


NRI, NRI2u, NRI2uk = sorted(NRI), sorted(NRI2u), sorted(NRI2uk)
# 计算偏移
NRI_offsets = calculate_offsets(NRI)
NRI2u_offsets = calculate_offsets(NRI2u)
NRI2uk_offsets = calculate_offsets(NRI2uk)

# 绘图
plt.figure(figsize=(6.5, 6), dpi=600)

# 小提琴图
plt.violinplot([NRI, NRI2u, NRI2uk], positions=positions, widths=2, 
               showextrema=False, showmedians=False, showmeans=False, 
               bw_method=0.5)

# 箱线图
plt.boxplot([NRI, NRI2u, NRI2uk], positions=positions, widths=0.9, 
            patch_artist=True, boxprops=dict(facecolor='none'), 
            whiskerprops=dict(color='black'),
            medianprops=dict(color='blue', linewidth=7))

# 散点图
plt.scatter(np.array(positions[0]) + np.array(NRI_offsets), 
            NRI, color='darkorange', label='NRI', alpha=0.5)
plt.scatter(np.array(positions[1]) + np.array(NRI2u_offsets), 
            NRI2u, color='deeppink', label='NRI$_{2u}$', alpha=0.5)
plt.scatter(np.array(positions[2]) + np.array(NRI2uk_offsets), 
            NRI2uk, color='green', label='NRI$_{2uk}$', alpha=0.5)

plt.xticks(positions, ['NRI', 'NRI$_{2u}$', 'NRI$_{2uk}$'],fontsize=25)
plt.ylabel('Diameter (mm)',fontsize=25)
plt.ylim(90,810)
plt.yticks(diameters[:],fontsize=17)
# plt.title('Violin Plot with Box and Scatter')

# plt.legend(loc="lower center", fontsize=16, frameon=True,ncol=3)
plt.grid(axis="x",c=(217/256,217/256,217/256),alpha=0.8)   #设置网格线,c为颜色
plt.grid(axis="y",c=(217/256,217/256,217/256),alpha=0.8)   #设置网格线,c为颜色

plt.tight_layout()
plt.savefig('******/PES_diaDis.svg', dpi=600)  

