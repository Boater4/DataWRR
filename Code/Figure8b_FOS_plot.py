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

wn = wntr.network.WaterNetworkModel('FOSinp/FOS_NRI_Sol93.inp')
diaList = []
for name, pipe in wn.pipes():
    diaList.append(1000*float(f'{pipe.diameter:.4f}'))
NRI = diaList

wn = wntr.network.WaterNetworkModel('FOSinp/FOS_NRI2u_Sol90.inp')
diaList = []
for name, pipe in wn.pipes():
    diaList.append(1000*float(f'{pipe.diameter:.4f}'))
NRI2u = diaList

wn = wntr.network.WaterNetworkModel('FOSinp/FOS_NRI2uk_Sol91.inp')
diaList = []
for name, pipe in wn.pipes():
    diaList.append(1000*float(f'{pipe.diameter:.4f}'))
NRI2uk = diaList

# 管径值
diameters = [16, 20.4, 26, 32.6, 40.8, 51.4, 61.4, 73.6, 90, 102.2, 114.6, 130.8, 147.2, 163.6, 184, 204.6, 229.2, 257.8, 290.6, 327.4, 368.2, 409.2]

    
def calculate_offsets(data):
    offset_data = []
    # data = sorted(data)
    interval = 0.25
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
positions = [0, 6, 12]


NRI, NRI2u, NRI2uk = sorted(NRI), sorted(NRI2u), sorted(NRI2uk)
# 计算偏移
NRI_offsets = calculate_offsets(NRI)
NRI2u_offsets = calculate_offsets(NRI2u)
NRI2uk_offsets = calculate_offsets(NRI2uk)

# 绘图
plt.figure(figsize=(6.5, 6), dpi=600)

# 小提琴图
plt.violinplot([NRI, NRI2u, NRI2uk], positions=positions, widths=1.5, 
               showextrema=False, showmedians=False, showmeans=False, 
               bw_method=0.5)

# 箱线图

plt.boxplot([NRI, NRI2u, NRI2uk], positions=positions, widths=0.7, 
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
plt.ylim(35,235)
# plt.yticks(np.linspace(2,16,15,endpoint=True))
plt.yticks(diameters[4:17],fontsize=17)



#plt.legend(loc="lower center", fontsize=16, frameon=True,ncol=3)
plt.grid(axis="x",c=(217/256,217/256,217/256),alpha=0.8)   #设置网格线,c为颜色
plt.grid(axis="y",c=(217/256,217/256,217/256),alpha=0.8)   #设置网格线,c为颜色

plt.tight_layout()
plt.savefig('******/FOS_diaDis.svg', dpi=600)  # 保存为PNG图像，dpi设置为600

