# -*- coding: UTF-8 -*-

# Author: kmihajlo

import os
import shutil

tasks_dir = 'tasks'

if os.path.exists('output'):
 	shutil.rmtree('output')
 	os.makedirs('output')  
else:
    os.makedirs('output') 

if os.path.exists(tasks_dir):
 	shutil.rmtree(tasks_dir)
 	os.makedirs(tasks_dir)  
else:
    os.makedirs(tasks_dir) 



cell_conds = os.listdir('input') 
k1_cc = {'ND_8':84, 'ND_18':83, 'ND_25':83, 'ND_32':79, 'ND_37':76, 'WT_8':83, 'WT_18':83, 'WT_25':84, 'WT_32':84, 'WT_37':80}
k2_cc = {'ND_8':36, 'ND_18':37, 'ND_25':37, 'ND_32':34, 'ND_37':21, 'WT_8':37, 'WT_18':37, 'WT_25':36, 'WT_32':37, 'WT_37':28}

k1s_cc_RA = {'ND_8':[k1_cc['ND_8']-30, k1_cc['ND_8']-15, k1_cc['ND_8'], k1_cc['ND_8']+15, k1_cc['ND_8']+30], 
             'ND_18':[k1_cc['ND_18']-30, k1_cc['ND_18']-15, k1_cc['ND_18'], k1_cc['ND_18']+15, k1_cc['ND_18']+30], 
             'ND_25':[k1_cc['ND_25']-30, k1_cc['ND_25']-15, k1_cc['ND_25'], k1_cc['ND_25']+15, k1_cc['ND_25']+30],  
             'ND_32':[k1_cc['ND_32']-30, k1_cc['ND_32']-15, k1_cc['ND_32'], k1_cc['ND_32']+15, k1_cc['ND_32']+30], 
             'ND_37':[k1_cc['ND_37']-30, k1_cc['ND_37']-15, k1_cc['ND_37'], k1_cc['ND_37']+15, k1_cc['ND_37']+30], 
             'WT_8':[k1_cc['WT_8']-30, k1_cc['WT_8']-15, k1_cc['WT_8'], k1_cc['WT_8']+15, k1_cc['WT_8']+30], 
             'WT_18':[k1_cc['WT_18']-30, k1_cc['WT_18']-15, k1_cc['WT_18'], k1_cc['WT_18']+15, k1_cc['WT_18']+30], 
             'WT_25':[k1_cc['WT_25']-30, k1_cc['WT_25']-15, k1_cc['WT_25'], k1_cc['WT_25']+15, k1_cc['WT_25']+30], 
             'WT_32':[k1_cc['WT_32']-30, k1_cc['WT_32']-15, k1_cc['WT_32'], k1_cc['WT_32']+15, k1_cc['WT_32']+30], 
             'WT_37':[k1_cc['WT_37']-30, k1_cc['WT_37']-15, k1_cc['WT_37'], k1_cc['WT_37']+15, k1_cc['WT_37']+30]}

k2s_cc_RA = {'ND_8':[k2_cc['ND_8']-5, k2_cc['ND_8']-10, k2_cc['ND_8'], k2_cc['ND_8']+5, k2_cc['ND_8']+10], 
             'ND_18':[k2_cc['ND_18']-5, k2_cc['ND_18']-10, k2_cc['ND_18'], k2_cc['ND_18']+5, k2_cc['ND_18']+10], 
             'ND_25':[k2_cc['ND_25']-5, k2_cc['ND_25']-10, k2_cc['ND_25'], k2_cc['ND_25']+5, k2_cc['ND_25']+10], 
             'ND_32':[k2_cc['ND_32']-5, k2_cc['ND_32']-10, k2_cc['ND_32'], k2_cc['ND_32']+5, k2_cc['ND_32']+10], 
             'ND_37':[k2_cc['ND_37']-5, k2_cc['ND_37']-10, k2_cc['ND_37'], k2_cc['ND_37']+5, k2_cc['ND_37']+10], 
             'WT_8':[k2_cc['WT_8']-5, k2_cc['WT_8']-10, k2_cc['WT_8'], k2_cc['WT_8']+5, k2_cc['WT_8']+10],
             'WT_18':[k2_cc['WT_18']-5, k2_cc['WT_18']-10, k2_cc['WT_18'], k2_cc['WT_18']+5, k2_cc['WT_18']+10], 
             'WT_25':[k2_cc['WT_25']-5, k2_cc['WT_25']-10, k2_cc['WT_25'], k2_cc['WT_25']+5, k2_cc['WT_25']+10], 
             'WT_32':[k2_cc['WT_32']-5, k2_cc['WT_32']-10, k2_cc['WT_32'], k2_cc['WT_32']+5, k2_cc['WT_32']+10], 
             'WT_37':[k2_cc['WT_37']-5, k2_cc['WT_37']-10, k2_cc['WT_37'], k2_cc['WT_37']+5, k2_cc['WT_37']+10]}

cell_conds = ['ND_8']
for cell_cond in cell_conds:   
    task_filename = cell_cond + '.txt'
    task_file = open(tasks_dir + '/' + task_filename, 'w')
    for k1 in k1s_cc_RA[cell_cond]:
        for k2 in k2s_cc_RA[cell_cond]:
            task_file.write('python Main_iCell.py' + ' ' + cell_cond + ' ' + str(k1) + ' ' + str(k2) + '\n')
    task_file.close()

for cell_cond in cell_conds: 
    task_filename = cell_cond + '.txt'
    with open(tasks_dir + '/' + task_filename, 'r') as file:
        # Read each line from the file
        for line in file:
            try:
                exec(line)
            except Exception as e:
                print("Error executing line:", e)
            
            
        
        
        
        
        
        
        