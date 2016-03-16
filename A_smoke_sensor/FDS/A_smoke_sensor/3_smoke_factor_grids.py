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

    if not os.path.exists('../3_sfgrids/dx_%.2f/%s_%.2f/Door_X_%.6f_Y_%.6f'%(delta_smoke_factor_grid, specified_location[0], specified_location[1],  exits[exit][0], exits[exit][1])):
        os.makedirs('../3_sfgrids/dx_%.2f/%s_%.2f/Door_X_%.6f_Y_%.6f'%(delta_smoke_factor_grid, specified_location[0], specified_location[1],  exits[exit][0], exits[exit][1]))

    smoke_factor_grid = np.ones([int(abs(dim2_max-dim2_min)/delta_smoke_factor_grid) ,int(abs(dim1_max-dim1_min)/delta_smoke_factor_grid)])

    x, y = m_to_pix(x_i=exits[exit][0], y_i=exits[exit][1])

    ### calculation of the global max of the extracted field_qunatity
    ### (without geometry)
    #global_D_max = np.amax(magnitudes-geometry[ :-1 , :-1])
    #print global_D_max

    for a, line in enumerate(smoke_factor_grid):

        for b, col in enumerate(line):

            x0, y0 = m_to_pix(b*delta_smoke_factor_grid+dim1_min, a*delta_smoke_factor_grid+dim2_min)    # virtual agent position

            #### stop of iteration and edge factor storage of predecessor cell
            #### if virtual agent position is identical with exit position
            if x0 == x and y0==y:
                smoke_factor = smoke_factor_grid[a,b-1]
                smoke_factor_grid[a,b]=smoke_factor
                continue

            x_exit, y_exit = np.linspace(x0, x, math.hypot(x - x0, y - y0)), np.linspace(y0, y, math.hypot(x - x0, y - y0))

            magnitude_along_line_of_sight = scipy.ndimage.map_coordinates(np.transpose(magnitudes), np.vstack((x_exit,y_exit)))

            #### consideration of visibility straight lines whose lengt is lt 1:
            if len(magnitude_along_line_of_sight)<1:
                continue

            else:
                #### Computation of the smoke factor towards the exit:
                smoke_factor = np.nanmax(magnitude_along_line_of_sight)*abs(np.trapz(magnitude_along_line_of_sight, dx=delta_dim_1))

            #### storage of the edge factor
            smoke_factor_grid[a,b]=smoke_factor

    #if np.amax(smoke_factor_grid)>max_Smoke_Factor:
    max_Smoke_Factor = np.amax(smoke_factor_grid)

    #print exit, 'maximum Smoke Factor', max_Smoke_Factor

    # Norm of the smoke factor grid multiplied with 10 in order to achieve a
    # factor range from 0..10
    smoke_factor_grid_norm = smoke_factor_grid/max_Smoke_Factor*10


    header='Room No. 1 , Exit %s, \n dX[m], dY[m] , minX[m] , maxX[m], minY[m], maxY[m] \n   %f  ,  %f    ,  %f  ,  %f  ,  %f ,  %f' \
    %(exit, delta_smoke_factor_grid, delta_smoke_factor_grid, dim1_min, dim1_max, dim2_min, dim2_max)
    np.savetxt('../3_sfgrids/dx_%.2f/%s_%.2f/Door_X_%.6f_Y_%.6f/t_%6f.csv'\
    %(delta_smoke_factor_grid, specified_location[0], specified_location[1],  exits[exit][0], exits[exit][1], float(time)), smoke_factor_grid_norm, header=header, delimiter=',', comments='')
    print 'Write smoke factor grid: ../3_sfgrids/dx_%.2f/%s_%.2f/Door_X_%.6f_Y_%.6f/t_%6f.csv'\
    %(delta_smoke_factor_grid, specified_location[0], specified_location[1],  exits[exit][0], exits[exit][1], float(time))

    return a,b, x0, y0, x, y, magnitude_along_line_of_sight, smoke_factor, time


#==============================================================================
# This is the resolution of the smoke factor grids: Adjustments may be
# appropriate in fligry geometries, Default: delta_smoke_factor_grid = 1 (m)
delta_smoke_factor_grid = 1.
#==============================================================================

print '\n-----------------------------------------------------------------'
print 'Analysing Slice File at %s=%.2f\n'%(specified_location[0], specified_location[1])

#plt.close()
converted = glob.glob('../2_consolidated/%s_%.2f/%s_*.csv'%(specified_location[0], specified_location[1], quantity))

print '\n-----------------------------------------------------------------'
print 'Generation of smoke factor grids... mesh resolution: %s m\n' %delta_smoke_factor_grid


for convert in converted:
    plt.close()
    print '\n-----------------------------------------------------------------'
    print 'Processing files: %s\n' %convert

    magnitudes = np.loadtxt(convert, delimiter=',')

    time=int(convert[:-4][convert.rfind('_')+1:])

    for id, exit in enumerate(exits):

        a,b, x0, y0, x, y, magnitude_along_line_of_sight, smoke_factor, time = main(convert, exit)


if plots==True:

    ## This statement yields a representative plot of the line of
    ## sights if wanted. x0 and y0 can be adjusted e.g. for debugging

    x0 = 3
    y0 = 3
    slicefile = '../2_consolidated/Z_2.25/SOOT_OPTICAL_DENSITY_120.csv'
    sfgrid = '../3_sfgrids/dx_1.00/Z_2.25/Door_X_2.000000_Y_5.500000/t_120.000000.csv'
    exit = 'cross_0'
    time = 120

    x, y = m_to_pix(x_i=exits[exit][0], y_i=exits[exit][1])
    x_exit, y_exit = np.linspace(x0, x, math.hypot(x - x0, y - y0)), np.linspace(y0, y, math.hypot(x - x0, y - y0))

    magnitude_along_line_of_sight = scipy.ndimage.map_coordinates(np.transpose(magnitudes), np.vstack((x_exit,y_exit)))


    fig = plt.figure()

    gs = gridspec.GridSpec(1, 40)

    ax1 = fig.add_subplot(gs[0,0:20])
    ax2 = fig.add_subplot(gs[0,19:22])
    ax3 = fig.add_subplot(gs[0,24:40])
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

    line, = ax1.plot(line_of_sight[0], line_of_sight[1], lw=2)

    slicefile = np.loadtxt(slicefile, delimiter=',')
    aa = ax1.pcolorfast(dim1, dim2, slicefile, cmap='Greys', vmin=0, vmax=1)
    ax1.minorticks_on()
    ax1.axis('image')

    ax1.grid(which='major',linestyle='-', alpha=0.4)
    ax1.grid(which='minor',linestyle='-', alpha=0.4)

    sfgrid = np.loadtxt(sfgrid, skiprows=3, delimiter=',')
    smoke_factor = sfgrid[ x0/delta_smoke_factor_grid, -y0/delta_smoke_factor_grid ]
    #print smoke_factor

    ax3.plot(magnitude_along_line_of_sight,  lw=2, label='%s: $f_{smoke}$ = %.3f'%(exit, smoke_factor))

    ax3.set_xlabel('l (m)')
    labels = ax3.get_xticks()
    labels = (labels*delta_dim_1).astype(int)
    ax3.set_xticklabels(labels)

    ax3.set_ylabel(quantity)
    ax3.set_ylim(0,1)

    cbar = fig.colorbar(aa,ax=ax1,cax=ax2, orientation='vertical')

    plt.legend(loc='upper left')
    plt.grid()

    plt.savefig('../2_consolidated/%s_%.2f/%s_%s.pdf'%(specified_location[0], specified_location[1], quantity+'_debug', time ))


    ## Plot Smoke factor grids - if specified
    times=np.arange(t_start,t_stop+1,t_step)

    for time in times:

        plt.close()

        fig = plt.figure()
        nrows = int(math.ceil( len(exits)**0.5 ))
        ncols = int(math.ceil( len(exits)**0.5 ))
        fig, axes = plt.subplots(nrows+1, ncols)

        gs = gridspec.GridSpec(nrows+1, ncols)

        for i, g in enumerate(gs):

            if i == len(exits):
                cbar_ax = plt.subplot(gs[-1, :])
                fig.colorbar(aa, cax=cbar_ax, label=r'$f_{smoke}$', orientation='horizontal')
                break

            exit = exits.items()[i][0]
            ax = plt.subplot(g)

            plt.xlabel('x (m)')
            plt.ylabel('y (m)')
            plt.tight_layout()
            smoke_factor_grid_norm = np.loadtxt('../3_sfgrids/dx_%.2f/%s_%.2f/Door_X_%.6f_Y_%.6f/t_%.6f.csv'%(delta_smoke_factor_grid, specified_location[0], specified_location[1],  exits[exit][0], exits[exit][1], time), delimiter=',', skiprows=3)
            aa = ax.pcolorfast(dim1, dim2, smoke_factor_grid_norm, cmap='coolwarm', vmin=0, vmax=10)
            ax.set_title(exit)
            ax.set_aspect('equal')
            ax.set_xticks(np.arange(dim1_min, dim1_max,5))
            ax.set_yticks(np.arange(dim2_min+1, dim2_max,5))
            ax.minorticks_on()
            ax.grid(which='major',linestyle='-', lw=1, alpha=0.4)
            ax.grid(which='minor',linestyle='-', lw=1, alpha=0.4)

        plt.savefig('../3_sfgrids/dx_%.2f/%s_%.2f/sfgrid_%i.pdf'%(delta_smoke_factor_grid, specified_location[0], specified_location[1],  time))

        print '\nPlot smoke factor grid: ../sfgrids/dx_%.2f/%s_%.2f/sfgrid_%i.pdf'%(delta_smoke_factor_grid, specified_location[0], specified_location[1],  time)

print "\n*** Finished ***"

#============storage of the most important variables to store.pckl=============

f = open('data_sfgrids.pckl', 'w')
pickle.dump((chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, jps_path, plots, dimension_1, dimension_2, dim1, dim2, delta_dim_1, delta_dim_2, geometry, magnitudes, dim1_min, dim1_max, dim2_min, dim2_max, exits, delta_smoke_factor_grid), f)

f.close()
