# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 15:44:05 2015

@author: Benjamin
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib
from matplotlib import rcParams, gridspec



rcParams.update({'figure.autolayout': True})

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 9}

matplotlib.rc('font', **font)

f = open('store.pckl')
(chid, dimension_1, dimension_2, dim1, dim2, delta_dim_1, delta_dim_2,\
geometry, dim1_min, dim1_max, dim2_min, dim2_max, exits) = pickle.load(f)


from mpltools import style
style.use('grayscale')

def cm2inch(value):
    return value/2.54

start=0
stop=390
step=15

##############################
#clip_dim3s = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
clip_dim3s = [3.0]
#delta_mesh_exits = [4.0, 2.0, 1.0, 0.5, 0.25]
delta_mesh_exits = [1]
##############################

for clip_dim3 in clip_dim3s:

    for delta_mesh_exit in delta_mesh_exits:

        times=np.arange(start,stop+1,step)

        for time in times[11:12]:

        #for time in times:

            test=[]
            plt.close()

            fig  = plt.figure(figsize=(cm2inch(15),cm2inch(12)), dpi=300)
            gs = gridspec.GridSpec(30, 3,bottom=0.18,left=0.18,right=0.88)
            #fig, axes = plt.subplots(figsize=(10, 6), ncols=3)
            #plt.subplots_adjust(wspace = 0.6)
            #plt.suptitle('Edge Factor Grids per Exit', fontsize=14)


            for ax, exit in enumerate(exits):

                # if ax>2:
                #     break

                if ax==1:
                    plt.ylabel('y (m)')

                ax1 = fig.add_subplot(gs[:24,ax])
                plt.xlabel('x (m)')


                mesh_exit_norm = np.loadtxt('../%s/sfgrids/z_%.1f/dx_%.2f/Door_X_%.6f_Y_%.6f/t_%.6f.csv'%(chid, clip_dim3, delta_mesh_exit, exits[exit][0], exits[exit][1], time), delimiter=',', skiprows=3)
                aa = ax1.pcolorfast(dim1, dim2, mesh_exit_norm, cmap='gray_r', vmin=0, vmax=5)

                ax1.text(exits[exit][0]-5, exits[exit][1], 'EXIT %s'%exit, bbox={'facecolor':'w', 'alpha':0.5, 'pad':5})
                ax1.set_aspect('equal')
                ax1.set_xticks(np.arange(dim1_min, dim1_max,5))
                ax1.set_yticks(np.arange(dim2_min+1, dim2_max,5))

            ax2 = fig.add_subplot(gs[29:30,0:6])
            fig.colorbar(aa,ax=ax1,cax=ax2, label=r'$f_{smoke}$', orientation='horizontal')



            plt.savefig('../%s/sfgrids/z_%.1f/dx_%.2f/%i.pdf'%(chid, clip_dim3, delta_mesh_exit, time))

            print '\nSave potential field: %s/sfgrids/z_%.1f/dx_%.2f/%i.pdf'%(chid, clip_dim3, delta_mesh_exit, time)

        #plt.show()

print "\n*** Finished ***"    