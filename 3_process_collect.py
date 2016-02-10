# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 10:25:22 2015

@author: Benjamin
"""

import glob
import numpy as np
import math
import scipy.ndimage
import matplotlib.pyplot as plt
import scipy.misc
import pickle
import os
import matplotlib
from matplotlib import rcParams, gridspec
import matplotlib.patches as patches
from scipy.interpolate import interp1d

rcParams.update({'figure.autolayout': True})

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 9}

matplotlib.rc('font', **font)


f = open('data_meshgrid.pckl')
(chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, jps_path, plots, dimension_1, dimension_2, dim1, dim2, delta_dim_1, delta_dim_2, geometry, dim1_min, dim1_max, dim2_min, dim2_max, exits) = pickle.load(f)

# Comversion cm to inches
def cm2inch(value):
    return value/2.54

# Comversion m to pixel coordinates
def m_to_pix(x_i, y_i):
    return ((abs(dim1_max-dim1_min) - dim1_max + x_i)/delta_dim_1, (abs(dim2_max-dim2_min) - dim2_max + y_i) / delta_dim_2 )

# main function
def main(convert, exit):

    if not os.path.exists('../3_sfgrids/%s_%.2f/dx_%.2f/Door_X_%.6f_Y_%.6f'%(specified_location[0], specified_location[1], delta_mesh_exit, exits[exit][0], exits[exit][1])):
        os.makedirs('../3_sfgrids/%s_%.2f/dx_%.2f/Door_X_%.6f_Y_%.6f'%(specified_location[0], specified_location[1], delta_mesh_exit, exits[exit][0], exits[exit][1]))

    mesh_exit = np.ones([int(abs(dim2_max-dim2_min)/delta_mesh_exit) ,int(abs(dim1_max-dim1_min)/delta_mesh_exit)])

    x, y = m_to_pix(x_i=exits[exit][0], y_i=exits[exit][1])

    global_D_max = 0

    for a, line in enumerate(mesh_exit):

        for b, col in enumerate(line):

            x0, y0 = m_to_pix(b*delta_mesh_exit+dim1_min, a*delta_mesh_exit+dim2_min)    # virtual agent position
            #print x0, y0

            #### stop of iteration and edge factor storage of predecessor cell
            #### if virtual agent position is identical with exit position
            if x0 == x and y0==y:
                edge_exit = mesh_exit[a,b-1]
                mesh_exit[a,b]=edge_exit
                continue

            #### shifting the x and y position of the exit according to the virtual agent position
            #### as well as the orientation of the exit to avoid to excessive &OBST slicing (Pixels)
            #### Uncomment if needed, else: x_shift, y_shift= 0, 0

            x_shift, y_shift= 0, 0

            # if x0 > x and exits[exit][2]=='y':
            #     x_shift=10
            # elif x0 < x and exits[exit][2]=='y':
            #     x_shift=-10
            #
            # if y0 > y and exits[exit][2]=='x':
            #     y_shift=10
            # elif y0 < y and exits[exit][2]=='x':
            #     y_shift=-10
            #
            # if exit=='C':
            #     y_shift=-4
            #

            x_exit, y_exit = np.linspace(x0, x+x_shift, math.hypot(x+x_shift - x0, y+y_shift - y0)), np.linspace(y0, y+y_shift, math.hypot(x+x_shift - x0, y+y_shift - y0))

            zi_exit = scipy.ndimage.map_coordinates(np.transpose(magnitudes), np.vstack((x_exit,y_exit)))

            #### consideration of visibility straight lines blocked by an &OBST:

            if len(zi_exit)<1:
                continue

            if np.amax(zi_exit)>8:
                edge_exit = 10

            #### visibility straight lines not or only moderately obstructed by smoke:

            else:

                edge_exit = np.amax(zi_exit)*abs(np.trapz(zi_exit))*delta_dim_1*2


                if np.amax(zi_exit) > global_D_max:
                    global_D_max=np.amax(zi_exit)


            #### storage of the edge factor
            mesh_exit[a,b]=edge_exit

            #print x0, y0

            ### This statement yields a representative plot of the line of
            ### sights if wanted
            if plots == True:
                if x0 == 24 and y0 == 4:
                    return a,b, x0, y0, x, y, zi_exit, edge_exit, x_shift, y_shift, time

    print 'Exit ', exit, 'maximum D(l)', global_D_max

    for a, line in enumerate(mesh_exit):
        for b, col in enumerate(line):
            if mesh_exit[a,b] == 10:
                #print exit
                continue
            else:
                #print global_D_max
                mesh_exit[a,b] = mesh_exit[a,b]/global_D_max


    header='Room No. 1 , Exit %s, \n dX[m], dY[m] , minX[m] , maxX[m], minY[m], maxY[m] \n   %f  ,  %f    ,  %f  ,  %f  ,  %f ,  %f' \
    %(exit, delta_mesh_exit, delta_mesh_exit, dim1_min, dim1_max, dim2_min, dim2_max)
    np.savetxt('../3_sfgrids/%s_%.2f/dx_%.2f/Door_X_%.6f_Y_%.6f/t_%6f.csv'\
    %(specified_location[0], specified_location[1], delta_mesh_exit, exits[exit][0], exits[exit][1], float(time)), mesh_exit, header=header, delimiter=',', comments='')
    print 'Write smoke factor grid: ../3_sfgrids/%s_%.2f/dx_%.2f/Door_X_%.6f_Y_%.6f/t_%6f.csv'\
    %(specified_location[0], specified_location[1], delta_mesh_exit, exits[exit][0], exits[exit][1], float(time))

    return a,b, x0, y0, x, y, zi_exit, edge_exit, x_shift, y_shift, time


#==============================================================================
# This is the resolution of the smoke sensor grids: Adjustments may be
# appropriate in fligry geometries, Default: delta_mesh_exit = 1 (m)
delta_mesh_exit = 1.
#==============================================================================

print '\n-----------------------------------------------------------------'
print 'Analysing Slice File at %s=%.2f\n'%(specified_location[0], specified_location[1])


#plt.close()
converted = glob.glob('../2_consolidated/%s_%.2f/%s_*.csv'%(specified_location[0], specified_location[1], quantity))

print '\n-----------------------------------------------------------------'
print 'Generation of potential fields... mesh resolution: %s m\n' %delta_mesh_exit

for convert in converted:

    plt.close()
    print '\n-----------------------------------------------------------------'
    print 'Processing files: %s\n' %convert
    magnitudes = np.loadtxt(convert, delimiter=',')

    time=int(convert[:-4][convert.rfind('_')+1:])

    fig = plt.figure(figsize=(cm2inch(15),cm2inch(9)), dpi=300)

    gs = gridspec.GridSpec(1, 15)

    #plt.suptitle('Time%s'%convert[:-4][convert.find('_'):], fontsize=14)

    ax1 = fig.add_subplot(gs[0,0:7])
    ax2 = fig.add_subplot(gs[0,6:8])
    ax3 = fig.add_subplot(gs[0,8:15])

    for id, exit in enumerate(exits):

        a,b, x0, y0, x, y, zi_exit, edge_exit, x_shift, y_shift, time = main(convert, exit)

        if plots==True:

            ax1.set_xlabel('%s (m)'%dimension_1.lower())
            ax1.set_ylabel('%s (m)'%dimension_2.lower())
            ax1.set_xticks(np.arange(dim1_min, dim1_max, 5))
            ax1.set_yticks(np.arange(dim2_min+1, dim2_max, 5))

            # Plot the line of sights towards each exit
            line_of_sight = [
            [x0*delta_dim_1,
            exits[exit][0]]
            ,
            [y0*delta_dim_2,
            exits[exit][1]],
            ]

            line, = ax1.plot(line_of_sight[0], line_of_sight[1])

            aa = ax1.pcolorfast(dim1, dim2, magnitudes, vmin=0, vmax=2)
            ax1.minorticks_on()
            ax1.axis('image')

            ax1.grid(which='major',linestyle='-', alpha=0.4)
            ax1.grid(which='minor',linestyle='-', alpha=0.4)

            ax3.plot(zi_exit,  lw=2, label='%s: $f_{smoke}$ = %.3f'%(exit, edge_exit))
            first_smoke=np.argmax(zi_exit>0.001)

            ax3.set_xlabel('l (m)')
            labels = ax3.get_xticks()
            labels = (labels*delta_dim_1).astype(int)
            ax3.set_xticklabels(labels)

            ax3.set_ylabel('D (1/m)')
            ax3.set_ylim(0,1)

        cbar = fig.colorbar(aa,ax=ax1,cax=ax2,
            orientation='vertical')


        plt.legend(fontsize=8, loc='upper left')

        plt.savefig('../2_consolidated/%s_%.2f/%s_%s.pdf'%(specified_location[0], specified_location[1], quantity+'_arrows', time ))

print "\n*** Finished ***"

#============storage of the most important variables to store.pckl=============

f = open('data_sfgrids.pckl', 'w')
pickle.dump((chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, jps_path, plots, dimension_1, dimension_2, dim1, dim2, delta_dim_1, delta_dim_2, geometry, dim1_min, dim1_max, dim2_min, dim2_max, exits, delta_mesh_exit), f)

f.close()
