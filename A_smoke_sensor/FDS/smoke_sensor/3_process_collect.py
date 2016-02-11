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
(chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, jps_path, plots, dimension_1, dimension_2, dim1, dim2, delta_dim_1, delta_dim_2, geometry, magnitudes, dim1_min, dim1_max, dim2_min, dim2_max, exits) = pickle.load(f)


# Comversion cm to inches
def cm2inch(value):
    return value/2.54

# Comversion m to pixel coordinates
def m_to_pix(x_i, y_i):
    return ((abs(dim1_max-dim1_min) - dim1_max + x_i)/delta_dim_1, (abs(dim2_max-dim2_min) - dim2_max + y_i) / delta_dim_2 )

# main function
def main(convert, exit):

    if not os.path.exists('../3_sfgrids/%s_%.2f/dx_%.2f/Door_X_%.6f_Y_%.6f'%(specified_location[0], specified_location[1], delta_smoke_factor_grid, exits[exit][0], exits[exit][1])):
        os.makedirs('../3_sfgrids/%s_%.2f/dx_%.2f/Door_X_%.6f_Y_%.6f'%(specified_location[0], specified_location[1], delta_smoke_factor_grid, exits[exit][0], exits[exit][1]))

    smoke_factor_grid = np.ones([int(abs(dim2_max-dim2_min)/delta_smoke_factor_grid) ,int(abs(dim1_max-dim1_min)/delta_smoke_factor_grid)])

    x, y = m_to_pix(x_i=exits[exit][0], y_i=exits[exit][1])

    ### calculation of the global max of the extracted field_qunatity
    ### (without geometry)
    #global_D_max = np.amax(magnitudes-geometry[ :-1 , :-1])
    #print global_D_max

    for a, line in enumerate(smoke_factor_grid):

        for b, col in enumerate(line):

            x0, y0 = m_to_pix(b*delta_smoke_factor_grid+dim1_min, a*delta_smoke_factor_grid+dim2_min)    # virtual agent position
            #print x0, y0

            #### stop of iteration and edge factor storage of predecessor cell
            #### if virtual agent position is identical with exit position
            if x0 == x and y0==y:
                smoke_factor = smoke_factor_grid[a,b-1]
                smoke_factor_grid[a,b]=smoke_factor
                continue

            #### shifting the x and y position of the exit according to the virtual agent position
            #### as well as the orientation of the exit to avoid to excessive &OBST slicing (Pixels)
            #### Uncomment if needed, else: x_shift, y_shift= 0, 0

            x_shift, y_shift= 0, 0

            # if x0 > x and exits[exit][2]=='y':
            #     x_shift=1/delta_dim_1
            # elif x0 < x and exits[exit][2]=='y':
            #     x_shift=-1/delta_dim_1
            #
            # if y0 > y and exits[exit][2]=='x':
            #     y_shift=1/delta_dim_2
            # elif y0 < y and exits[exit][2]=='x':
            #     y_shift=-1/delta_dim_2

            x_exit, y_exit = np.linspace(x0, x+x_shift, math.hypot(x+x_shift - x0, y+y_shift - y0)), np.linspace(y0, y+y_shift, math.hypot(x+x_shift - x0, y+y_shift - y0))

            line_of_sight = scipy.ndimage.map_coordinates(np.transpose(magnitudes), np.vstack((x_exit,y_exit)))


            #### consideration of visibility straight lines whose lengt is lt 1:
            if len(line_of_sight)<1:
                continue

            else:
                #### Computation of the edge factor towards the exit:
                smoke_factor = np.amax(line_of_sight)*abs(np.trapz(line_of_sight, dx=delta_dim_1))

            #### storage of the edge factor
            smoke_factor_grid[a,b]=smoke_factor

            ## This statement yields a representative plot of the line of
            ## sights if wanted. x0 and y0 can be adjusted e.g. for debugging
            if plots == True:
                if x0 == 24 and y0 == 4:
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

                    ax3.plot(line_of_sight,  lw=2, label='%s: $f_{smoke}$ = %.3f'%(exit, smoke_factor))
                    first_smoke=np.argmax(line_of_sight>0.001)

                    ax3.set_xlabel('l (m)')
                    labels = ax3.get_xticks()
                    labels = (labels*delta_dim_1).astype(int)
                    ax3.set_xticklabels(labels)

                    ax3.set_ylabel('D (1/m)')
                    ax3.set_ylim(0,1)

                    cbar = fig.colorbar(aa,ax=ax1,cax=ax2, orientation='vertical')

                    plt.legend(fontsize=8, loc='upper left')

                    plt.savefig('../2_consolidated/%s_%.2f/%s_%s.pdf'%(specified_location[0], specified_location[1], quantity+'_arrows', time ))

    #if np.amax(smoke_factor_grid)>max_Smoke_Factor:
    max_Smoke_Factor = np.amax(smoke_factor_grid)

    print exit, 'maximum Smoke Factor', max_Smoke_Factor

    # Norm of the smoke factor grid multiplied with 10 in order to achieve a
    # factor range from 0..10
    smoke_factor_grid_norm = smoke_factor_grid/max_Smoke_Factor*10


    header='Room No. 1 , Exit %s, \n dX[m], dY[m] , minX[m] , maxX[m], minY[m], maxY[m] \n   %f  ,  %f    ,  %f  ,  %f  ,  %f ,  %f' \
    %(exit, delta_smoke_factor_grid, delta_smoke_factor_grid, dim1_min, dim1_max, dim2_min, dim2_max)
    np.savetxt('../3_sfgrids/%s_%.2f/dx_%.2f/Door_X_%.6f_Y_%.6f/t_%6f.csv'\
    %(specified_location[0], specified_location[1], delta_smoke_factor_grid, exits[exit][0], exits[exit][1], float(time)), smoke_factor_grid_norm, header=header, delimiter=',', comments='')
    print 'Write smoke factor grid: ../3_sfgrids/%s_%.2f/dx_%.2f/Door_X_%.6f_Y_%.6f/t_%6f.csv'\
    %(specified_location[0], specified_location[1], delta_smoke_factor_grid, exits[exit][0], exits[exit][1], float(time))

    return a,b, x0, y0, x, y, line_of_sight, smoke_factor, x_shift, y_shift, time


#==============================================================================
# This is the resolution of the smoke sensor grids: Adjustments may be
# appropriate in fligry geometries, Default: delta_smoke_factor_grid = 1 (m)
delta_smoke_factor_grid = 1.
#==============================================================================

print '\n-----------------------------------------------------------------'
print 'Analysing Slice File at %s=%.2f\n'%(specified_location[0], specified_location[1])

#plt.close()
converted = glob.glob('../2_consolidated/%s_%.2f/%s_*.csv'%(specified_location[0], specified_location[1], quantity))

print '\n-----------------------------------------------------------------'
print 'Generation of potential fields... mesh resolution: %s m\n' %delta_smoke_factor_grid

for convert in converted:

    plt.close()
    print '\n-----------------------------------------------------------------'
    print 'Processing files: %s\n' %convert

    magnitudes = np.loadtxt(convert, delimiter=',')-geometry[:-1,:-1]

    time=int(convert[:-4][convert.rfind('_')+1:])

    fig = plt.figure(figsize=(cm2inch(15),cm2inch(9)), dpi=300)

    gs = gridspec.GridSpec(1, 15)

    #plt.suptitle('Time %i s'%time, fontsize=12)

    ax1 = fig.add_subplot(gs[0,0:7])
    ax2 = fig.add_subplot(gs[0,6:8])
    ax3 = fig.add_subplot(gs[0,8:15])

    for id, exit in enumerate(exits):

        a,b, x0, y0, x, y, line_of_sight, smoke_factor, x_shift, y_shift, time = main(convert, exit)

print "\n*** Finished ***"

#============storage of the most important variables to store.pckl=============

f = open('data_sfgrids.pckl', 'w')
pickle.dump((chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, jps_path, plots, dimension_1, dimension_2, dim1, dim2, delta_dim_1, delta_dim_2, geometry, magnitudes, dim1_min, dim1_max, dim2_min, dim2_max, exits, delta_smoke_factor_grid), f)

f.close()
