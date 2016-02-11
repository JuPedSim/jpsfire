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

def cm2inch(value):
    return value/2.54

f = open('data_sfgrids.pckl')
(chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, jps_path, plots, dimension_1, dimension_2, dim1, dim2, delta_dim_1, delta_dim_2, geometry, magnitudes, dim1_min, dim1_max, dim2_min, dim2_max, exits, delta_smoke_factor_grid) = pickle.load(f)

if plots==True:

    times=np.arange(t_start,t_stop+1,t_step)

    for time in times:

        test=[]
        plt.close()

        fig  = plt.figure(figsize=(cm2inch(21),cm2inch(10)), dpi=300)
        gs = gridspec.GridSpec(15, 4,bottom=0.18,left=0.18,right=0.88)
        #fig, axes = plt.subplots(figsize=(10, 6), ncols=3)
        #plt.subplots_adjust(wspace = 0.6)
        #plt.suptitle('Edge Factor Grids per Exit', fontsize=14)


        for ax, exit in enumerate(exits):

            # if ax>2:
            #     break

            if ax==1:
                plt.ylabel('y (m)')

            ax1 = fig.add_subplot(gs[:12,ax])
            plt.xlabel('x (m)')


            smoke_factor_grid_norm = np.loadtxt('../3_sfgrids/%s_%.2f/dx_%.2f/Door_X_%.6f_Y_%.6f/t_%.6f.csv'%(specified_location[0], specified_location[1], delta_smoke_factor_grid, exits[exit][0], exits[exit][1], time), delimiter=',', skiprows=3)
            aa = ax1.pcolorfast(dim1, dim2, smoke_factor_grid_norm, cmap='coolwarm')#, vmin=0, vmax=10)
            ax1.set_title(exit)
            #ax1.text(exits[exit][0]-5, exits[exit][1], 'EXIT %s'%exit, bbox={'facecolor':'w', 'alpha':0.5, 'pad':5})
            ax1.set_aspect('equal')
            ax1.set_xticks(np.arange(dim1_min, dim1_max,5))
            ax1.set_yticks(np.arange(dim2_min+1, dim2_max,5))
            ax1.minorticks_on()
            ax1.grid(which='major',linestyle='-', lw=1, alpha=0.4)
            ax1.grid(which='minor',linestyle='-', lw=1, alpha=0.4)

        ax2 = fig.add_subplot(gs[14:15,0:len(exits)])
        fig.colorbar(aa,ax=ax1,cax=ax2, label=r'$f_{smoke}$', orientation='horizontal')

        plt.savefig('../3_sfgrids/%s_%.2f/dx_%.2f/sfgrid_%i.pdf'%(specified_location[0], specified_location[1], delta_smoke_factor_grid, time))

        print '\nSave smoke factor grid: ../sfgrids/%s_%.2f/dx_%.2f/sfgrid_%i.pdf'%(specified_location[0], specified_location[1], delta_smoke_factor_grid, time)

        #plt.show()

print "\n*** Finished ***"
