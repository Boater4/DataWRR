# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 17:48:01 2024

@author: 10249
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import viswaternet as vis
import wntr

path = '******'
os.chdir(path)

flist = ['HAN_NRI_Sol49.inp','HAN_NRI2u_Sol66.inp','HAN_NRI2uk_Sol92.inp',
         'FOS_NRI_Sol93.inp','FOS_NRI2u_Sol90.inp','FOS_NRI2uk_Sol91.inp',
         'PES_NRI_Sol362.inp','PES_NRI2u_Sol247.inp','PES_NRI2uk_Sol381.inp',
         'MOD_NRI_Sol130.inp','MOD_NRI2u_Sol227.inp','MOD_NRI2uk_Sol137.inp']
locations = {'HAN':-0.2,'FOS':-0.23,'PES':-0.17,'MOD':-0.2}
diameters = {'HAN':6,'FOS':25,'PES':13,'MOD':13}
ups = {'HAN':15,'FOS':23,'PES':16,'MOD':15}
units = {'HAN':'in','FOS':'mm','PES':'mm','MOD':'mm'}
points = {'HAN':0,'FOS':1,'PES':0,'MOD':0}
 
for file in flist[:]:
    case = file.split('_')[0]
    
    model = vis.VisWNModel(f'{case}inp/{file}')
    
    model.plot_unique_data(parameter = "diameter", 
                       unit = units[case], 
                       link_width = np.linspace(9,ups[case],diameters[case],endpoint=True),
                       cmap = "Blues", draw_base_legend=False,
                       discrete_legend_loc =(locations[case],0), discrete_legend_title_font_size=30,
                       legend_title = f" Diameter\n    ({units[case]})", 
                       discrete_legend_label_font_size = 30,
                       legend_decimal_places = points[case], 
                       reservoir_shape='s',
                       reservoir_size=650, 
                       reservoir_color='white',
                       reservoir_border_color='deeppink', 
                       reservoir_border_width=3)
    
    plt.savefig(f"******/Fluct_{file.split('.')[0]}.svg",dpi=600,bbox_inches='tight')
