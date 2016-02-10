# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 15:59:35 2015

@author: Benjamin
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 12:50:15 2015

@author: Benjamin
"""
import numpy as np
import glob
import matplotlib.pyplot as plt
import sys
import re
import pickle
import os
import collections
import matplotlib
from matplotlib import rcParams


rcParams.update({'figure.autolayout': True})

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 9}

matplotlib.rc('font', **font)

#style.use('grayscale')

def cm2inch(value):
    return value/2.54

def dimensions(a):

    if a =='X':
        return 0,1
    elif a =='Y':
        return 2,3
    else:
        return 4,5
        
delta_dim_1 = 0.2
delta_dim_2 = 0.2

times = np.arange(0,330,30)[:]
mesh_ids = [1,2,3]



def mesh_attributes(time, mesh_id):

    sf=np.loadtxt('../../1_mesh_test_3/osloer_mesh_test_3_dx20_v2/ascii/9.2/SOOT_%i_mesh_%i.txt'%(time, mesh_id), skiprows=2, delimiter=',')
    
    dim1_min=round(np.amin(sf[:,0]),3)
    dim1_max=round(np.amax(sf[:,0]),3)
    
    dim2_min=round(np.amin(sf[:,1]),3)
    dim2_max=round(np.amax(sf[:,1]),3)
    
    dim1=np.arange(dim1_min, dim1_max+0.01, delta_dim_1)
    dim2=np.arange(dim2_min, dim2_max+0.01, delta_dim_2) 
    
    off_dim1 = len(dim1)
    off_dim2 = len(dim2)
    
    return dim1, dim2, off_dim1, off_dim2
    

for time in times:

    consolidated = np.zeros((1000,1000))
    consolidated[:] = np.NAN
    
    dim1_cons = []
    dim2_cons = []
    
    global_offset_dim1 = 500
    global_offset_dim2 = 500

    offset_dim1 = [0]     
    offset_dim2 = [0]    
    dim1_cons = []
    dim2_cons = []    
    
    plt.figure(figsize=(cm2inch(15),cm2inch(10)), dpi=300)
    ax=plt.subplot(111)

    for mesh_id in mesh_ids:
        print '../../1_mesh_test_3/osloer_mesh_test_3_dx20_v2/converted/9.2/SOOT_%i_mesh_%i.csv'%(time, mesh_id)
        
        collect = np.loadtxt('../../1_mesh_test_3/osloer_mesh_test_3_dx20_v2/converted/9.2/SOOT_%i_mesh_%i.csv'%(time, mesh_id), delimiter=',')

        dim1 = offset_dim1, mesh_attributes(time, mesh_id)[0]
        dim2 = offset_dim2, mesh_attributes(time, mesh_id)[1]
                
        dim1_cons = np.append(dim1_cons, dim1[1])
        dim2_cons = np.append(dim2_cons, dim2[1])

       
        consolidated[ 
        
        #rows
        
        global_offset_dim2 + np.amin(dim2[1])*5 
        : 
        global_offset_dim2 + np.amax(dim2[1])*5 + 1
        ,
        
        #cols
        global_offset_dim1 + np.amin(dim1[1])*5 
        : 
        global_offset_dim1 + np.amax(dim1[1])*5 + 1
        
        ] = collect

    
    dim1_resize = []    
    dim2_resize = []
    
    for i, row in enumerate(consolidated):       
        if np.isnan(np.nanmean(row)) == False:
            dim2_resize = np.append(dim2_resize, i)
    
    for i, col in enumerate(consolidated.T):       
        if np.isnan(np.nanmean(col)) == False:
            dim1_resize = np.append(dim1_resize, i)
      
      
    consolidated = consolidated[
    dim2_resize[0] : dim2_resize[-1]
    ,
    dim1_resize[0] : dim1_resize[-1] 
    ]

    dim1_global_min = np.amin(dim1_cons)
    dim2_global_min = np.amin(dim2_cons)
    
    dim1_global_max = np.amax(dim1_cons)
    dim2_global_max = np.amax(dim2_cons)
    
    new_dim1 = np.arange(dim1_global_min, dim1_global_max, 0.2)    
    new_dim2 = np.arange(dim2_global_min, dim2_global_max, 0.2) 
    
    aa = ax.pcolorfast(new_dim1, new_dim2, consolidated, vmin=0, vmax=3)
        
    plt.colorbar(aa)
    plt.axes().set_aspect('equal')
       
    #plt.xlabel('%s (m)'%dimension_1.lower())
    #plt.ylabel('%s (m)'%dimension_2.lower())

    #plt.xticks([-5,0,5,10])


    #plt.minorticks_on()
    #plt.grid(which='major')

    plt.savefig('SOOT_mesh_test_3_%i.png'%time)
    
    plt.close()
#plt.show()

print "\n*** Finished ***"
    