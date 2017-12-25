#!/usr/bin/env python

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
import multiprocessing as mp

CPUs = mp.cpu_count()
print CPUs

rcParams.update({'figure.autolayout': True})

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 9}

matplotlib.rc('font', **font)


def cm2inch(value):
    return value/2.54


def mesh_attributes(time, mesh_id):

    sf=np.loadtxt('../0_slice2ascii/%s_%.6f/%s_%i_mesh_%i.txt'%(specified_location[0], specified_location[1], quantity, time, mesh_id), skiprows=2, delimiter=',')

    delta_dim_1 = abs(np.max(np.diff(sf[:,0])))
    delta_dim_2 = abs(np.max(np.diff(sf[:,1])))

    dim1_min=round(np.amin(sf[:,0]),3)
    dim1_max=round(np.amax(sf[:,0]),3)

    dim2_min=round(np.amin(sf[:,1]),3)
    dim2_max=round(np.amax(sf[:,1]),3)

    dim1=np.arange(dim1_min, dim1_max+0.01, delta_dim_1)
    dim2=np.arange(dim2_min, dim2_max+0.01, delta_dim_2)

    off_dim1 = len(dim1)
    off_dim2 = len(dim2)

    return dim1, dim2, off_dim1, off_dim2, delta_dim_1, delta_dim_2, dimension_1, dimension_2


f = open('data_meshgrid.pckl')
(chid, quantity, specified_location, t_low, id_meshes, plots, dimension_1, dimension_2, dim1, dim2, delta_dim_1, delta_dim_2, geometry, magnitudes, dim1_min, dim1_max, dim2_min, dim2_max) = pickle.load(f)

if not os.path.exists('../2_extinction_grids/%s/%s_%.6f'%(quantity, specified_location[0], specified_location[1])):
    print 'create directory "2_extinction_grids"'
    os.makedirs('../2_extinction_grids/%s/%s_%.6f'%(quantity, specified_location[0], specified_location[1]))

for time in t_low:

    consolidated = np.zeros((5000,5000))
    consolidated[:] = np.NAN

    dim1_cons = []
    dim2_cons = []

    global_offset_dim1 = 1500
    global_offset_dim2 = 1500

    offset_dim1 = [0]
    offset_dim2 = [0]
    dim1_cons = []
    dim2_cons = []

    plt.figure(figsize=(cm2inch(20),cm2inch(15)), dpi=300)
    ax=plt.subplot(111)

    for mesh_id in id_meshes:
        print '%s_%.6f/%s_%i_mesh_%i.csv'%(specified_location[0], specified_location[1], quantity, time, mesh_id)

        collect = np.loadtxt('../1_meshgrid/%s_%.6f/%s_%i_mesh_%i.csv'%(specified_location[0], specified_location[1], quantity, time, mesh_id), delimiter=',')

        dim1 = offset_dim1, mesh_attributes(time, mesh_id)[0]
        dim2 = offset_dim2, mesh_attributes(time, mesh_id)[1]

        dim1_cons = np.append(dim1_cons, dim1[1])
        dim2_cons = np.append(dim2_cons, dim2[1])


        consolidated[

        #rows
        int(global_offset_dim2 + np.amin(dim2[1])/delta_dim_2)
        :
        int(global_offset_dim2 + np.amin(dim2[1])/delta_dim_2 + np.shape(collect)[0])
        ,

        #cols
        int(global_offset_dim1 + np.amin(dim1[1])*(1/delta_dim_1))
        :
        int(global_offset_dim1 + np.amin(dim1[1])*(1/delta_dim_1) + np.shape(collect)[1])

        ] = collect


    dim1_resize = []
    dim2_resize = []

    for i, row in enumerate(consolidated):
        if np.isnan(np.nanmean(row)) == False:
            dim2_resize = np.append(dim2_resize, i)

    for i, col in enumerate(consolidated.T):
        if np.isnan(np.nanmean(col)) == False:
            dim1_resize = np.append(dim1_resize, i)
            
    dim1_resize = dim1_resize.astype(int)
    dim2_resize = dim2_resize.astype(int)   

    consolidated = consolidated[
    dim2_resize[0] : dim2_resize[-1]
    ,
    dim1_resize[0] : dim1_resize[-1]
    ]

    dim1_global_min = np.amin(dim1_cons)
    dim2_global_min = np.amin(dim2_cons)

    dim1_global_max = np.amax(dim1_cons)
    dim2_global_max = np.amax(dim2_cons)

    new_dim1 = np.arange(dim1_global_min, dim1_global_max, delta_dim_1)
    new_dim2 = np.arange(dim2_global_min, dim2_global_max, delta_dim_2)

    header='Room \n dX[m], dY[m] , minX[m] , maxX[m], minY[m], maxY[m] \n  %f  , %f ,  %f  ,  %f  ,  %f ,  %f' \
    %(delta_dim_1, delta_dim_2, dim1_min, dim1_max, dim2_min, dim2_max)

    np.savetxt('../2_extinction_grids/%s/%s_%.6f/t_%.6f.csv'%(quantity, specified_location[0], specified_location[1], time), consolidated, header=header, delimiter=',', comments='')

    if plots == True:
        aa = ax.pcolorfast(new_dim1, new_dim2, consolidated, vmin=0, vmax=5)

        plt.colorbar(aa, label=quantity)
        plt.axes().set_aspect('equal')

        plt.xlabel('%s (m)'%dimension_1.lower())
        plt.ylabel('%s (m)'%dimension_2.lower())

        #plt.minorticks_on()
        #plt.grid(which='major')

        plt.savefig('../2_extinction_grids/%s/%s_%.6f/t_%.6f.pdf'%(quantity, specified_location[0], specified_location[1], time))

        plt.close()

print "\n*** Finished ***"
