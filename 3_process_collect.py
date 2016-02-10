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
(chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, dimension_1, dimension_2, dim1, dim2, delta_dim_1, delta_dim_2, geometry, dim1_min, dim1_max, dim2_min, dim2_max, exits) = pickle.load(f)

#### Import of a function for optical density vs. visibility:
#density_visibility=np.loadtxt('HLKL_1lx.txt', skiprows=1)
#f_Dl = interp1d(density_visibility[:,0], density_visibility[:,1])

# plt.plot(np.arange(0.00,1.2,0.01),f_Dl(np.arange(0.00,1.2,0.01)))
# plt.show()

def cm2inch(value):
    return value/2.54

def m_to_pix(x_i, y_i):
    return ((abs(dim1_max-dim1_min) - dim1_max + x_i)/delta_dim_1, (abs(dim2_max-dim2_min) - dim2_max + y_i) / delta_dim_2 )

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
            #x0=x0 + 5
            #y0=y0 + 5

            #### stop of iteration and edge factor storage of predecessor cell
            #### if virtual agent position is identical with exit position
            if x0 == x and y0==y:
                edge_exit = mesh_exit[a,b-1]
                mesh_exit[a,b]=edge_exit
                continue

            #### shifting the x and y position of the exit according to the virtual agent position
            #### as well as the orientation of the exit to avoid to excessive &OBST slicing (Pixels)

            x_shift, y_shift= 0, 0

            if x0 > x and exits[exit][2]=='y':
                x_shift=10
            elif x0 < x and exits[exit][2]=='y':
                x_shift=-10

            if y0 > y and exits[exit][2]=='x':
                y_shift=10
            elif y0 < y and exits[exit][2]=='x':
                y_shift=-10

            if exit=='C':
                y_shift=-4

            aaa=0
            if exit=='D' or exit=='E' or exit=='F':
                aaa=-10

            x_exit, y_exit = np.linspace(x0+aaa, x+x_shift, math.hypot(x+x_shift - x0, y+y_shift - y0)), np.linspace(y0, y+y_shift, math.hypot(x+x_shift - x0, y+y_shift - y0))


            # if np.isnan(magnitudes[np.amin(y_exit):np.amax(y_exit) , np.amin(x_exit):np.amax(x_exit)]).any()==True:
            #     edge_exit = 1
            #     #print exit, 'NAN!'
            #     zi_exit = scipy.ndimage.map_coordinates(np.transpose(magnitudes), np.vstack((x_exit,y_exit)))
            #
            # else:
            #     if np.amin(y_exit)==np.amax(y_exit):
            #         test=magnitudes[np.amin(y_exit) , np.amin(x_exit):np.amax(x_exit)]
            #         zi_exit=test
            #         #print zi_exit
            #
            #     elif np.amin(x_exit)==np.amax(x_exit):
            #         test=magnitudes[np.amin(y_exit):np.amax(y_exit) , np.amin(x_exit)]
            #         zi_exit=test
            #
            #     else:
            #         test=magnitudes[np.amin(y_exit):np.amax(y_exit) , np.amin(x_exit):np.amax(x_exit)]
            #         zi_exit = scipy.ndimage.map_coordinates(np.transpose(test), np.vstack((x_exit-np.amin(x_exit), y_exit-np.amin(y_exit))))

            #edge_exit = np.amax(zi_exit)*abs(np.trapz(zi_exit))/(len(zi_exit)*delta_dim_1) + 1
            #     #edge_exit = abs(np.trapz(zi_exit))/len(zi_exit) + 1

            zi_exit = scipy.ndimage.map_coordinates(np.transpose(magnitudes), np.vstack((x_exit,y_exit)))
            #edge_exit = np.amax(zi_exit)*abs(np.trapz(zi_exit))/(len(zi_exit)*delta_dim_1) + 1

            # #### consideration if straight line length is higher than visibility throughout smoke:
            # first_smoke=np.argmax(zi_exit>0.001)
            #
            # if f_Dl(np.mean(zi_exit[first_smoke:])) < (len(zi_exit)-first_smoke)*delta_dim_1 and np.amax(zi_exit)<8:
            #     print len(zi_exit)
            #     print 'No visibility!'
            #     edge_exit = 10 + np.amax(zi_exit)*abs(np.trapz(zi_exit))/(len(zi_exit)*delta_dim_1)

            #### consideration of visibility straight lines blocked by an &OBST:

            if len(zi_exit)<1:
                continue

            if np.amax(zi_exit)>8:
                edge_exit = 10

            #### visibility straight lines not or only moderately obstructed by smoke:

            else:

                # print 'min', np.amin(zi_exit)
                # print 'max', np.amax(zi_exit)
                #
                # zi_exit=zi_exit*np.arange(0,len(zi_exit))/len(zi_exit)
                #
                # print 'min_norm', np.amin(zi_exit)
                # print 'max_norm', np.amax(zi_exit)

                #edge_exit = np.amax(zi_exit)*abs(np.trapz(zi_exit))*delta_dim_1 # /(len(zi_exit)*delta_dim_1) # np.amax(zi_exit) * #/(len(zi_exit)*delta_dim_1)# + 1
                #edge_exit = abs(np.trapz(zi_exit))*delta_dim_1 # /(len(zi_exit)*delta_dim_1) # np.amax(zi_exit) * #/(len(zi_exit)*delta_dim_1)# + 1

                edge_exit = np.amax(zi_exit)*abs(np.trapz(zi_exit))*delta_dim_1*2


                if np.amax(zi_exit) > global_D_max:
                    global_D_max=np.amax(zi_exit)

                #if np.amax(zi_exit)>1:
                    #return a,b, x0, y0, x, y, zi_exit, edge_exit, x_shift, y_shift

            #### storage of the edge factor
            mesh_exit[a,b]=edge_exit

            # if x0 == 110 and y0==230:
            #    return a,b, x0, y0, x, y, zi_exit, edge_exit, x_shift, y_shift

    print 'Exit ', exit, 'maximum D(l)', global_D_max

    for a, line in enumerate(mesh_exit):
        for b, col in enumerate(line):
            if mesh_exit[a,b] == 10:
                #print exit
                continue
            else:
                #print global_D_max
                mesh_exit[a,b] = mesh_exit[a,b]/global_D_max


    time=int(convert[:-4][convert.rfind('_')+1:])


    header='Room No. 1 , Exit %s, \n dX[m], dY[m] , minX[m] , maxX[m], minY[m], maxY[m] \n   %f  ,  %f    ,  %f  ,  %f  ,  %f ,  %f' \
    %(exit, delta_mesh_exit, delta_mesh_exit, dim1_min, dim1_max, dim2_min, dim2_max)
    np.savetxt('../3_sfgrids/%s_%.2f/dx_%.2f/Door_X_%.6f_Y_%.6f/t_%6f.csv'\
    %(specified_location[0], specified_location[1], delta_mesh_exit, exits[exit][0], exits[exit][1], float(time)), mesh_exit, header=header, delimiter=',', comments='')
    print 'Write smoke factor grid: ../3_sfgrids/%s_%.2f/dx_%.2f/Door_X_%.6f_Y_%.6f/t_%6f.csv'\
    %(specified_location[0], specified_location[1], delta_mesh_exit, exits[exit][0], exits[exit][1], float(time))

    return a,b, x0, y0, x, y, zi_exit, edge_exit, x_shift, y_shift, time


#==============================================================================

delta_mesh_exit = 1.0

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

    fig = plt.figure(figsize=(cm2inch(15),cm2inch(9)), dpi=300)

    gs = gridspec.GridSpec(1, 15)

    #plt.suptitle('Time%s'%convert[:-4][convert.find('_'):], fontsize=14)

    colors = ['0.1','0.3','0.5','0.7','0.85','1.0']
    ax1 = fig.add_subplot(gs[0,0:7])
    ax2 = fig.add_subplot(gs[0,6:8])
    ax3 = fig.add_subplot(gs[0,8:15])

    for id, exit in enumerate(exits):

        a,b, x0, y0, x, y, zi_exit, edge_exit, x_shift, y_shift, time = main(convert, exit)

        ax1.set_xlabel('%s (m)'%dimension_1.lower())
        ax1.set_ylabel('%s (m)'%dimension_2.lower())
        ax1.set_xticks(np.arange(dim1_min, dim1_max, 5))
        ax1.set_yticks(np.arange(dim2_min+1, dim2_max, 5))


        ax1.add_patch(
        patches.FancyArrow(

        b*delta_mesh_exit+dim1_min+delta_mesh_exit/2.0, \
        a*delta_mesh_exit+dim2_min+delta_mesh_exit/2.0, \
        #-1.0*((b)*delta_mesh_exit+dim1_min-exits[exit][0])+x_shift*delta_dim_1, \

        #-1.0*((a)*delta_mesh_exit+dim2_min-exits[exit][1])+y_shift*delta_dim_1,

        exits[exit][0]-(b*delta_mesh_exit+dim1_min+delta_mesh_exit/2.0)+x_shift/10, \

        exits[exit][1]-(a*delta_mesh_exit+dim2_min+delta_mesh_exit/2.0)+y_shift/10,

        alpha=None, color=colors[id],
        width=0.25,        # Default
        linewidth=0.1,
        head_width=0.75,    # Default: 3 * width
        head_length=0.75,
        ))    # Default: 1.5 * head_width


        aa = ax1.pcolorfast(dim1, dim2, magnitudes, cmap='gray_r', vmin=0, vmax=1)
        ax1.minorticks_on()
        ax1.axis('image')
        #ax1.grid(which='major',linestyle='-', alpha=0.4)
        #ax1.grid(which='minor',linestyle='-', alpha=0.4)

        ax3.plot(zi_exit,  lw=2, color=colors[id], label='Exit %s: $f_{smoke}$ = %.3f'%(exit, edge_exit))
        first_smoke=np.argmax(zi_exit>0.001)
        #ax3.axhline(np.mean(zi_exit[first_smoke:]),color=colors[id])

        ax1.text(13,30, 'Room 3', bbox={'facecolor':'w', 'pad':4})
        ax1.text(13,20, 'Room 1', bbox={'facecolor':'w', 'pad':4})

        ax1.text(13,10, 'Agent \n $P_A :=(x_1|y_1)$', bbox={'facecolor':'w', 'pad':4})
        ax1.text(13,0, 'Exit \n $P_E :=(x_2|y_2)$', bbox={'facecolor':'w', 'pad':4})

        ax3.set_xlabel('l (m)')
        labels = ax3.get_xticks()
        labels = (labels*delta_dim_1).astype(int)
        ax3.set_xticklabels(labels)

        ax3.set_ylabel('D (1/m)')
        ax3.set_ylim(0,1)

    cbar = fig.colorbar(aa,ax=ax1,cax=ax2,
        orientation='vertical')


    #plt.legend(fontsize=9, loc='upper left')

    #plt.savefig('../2_consilidated/%s_%.2f/smoke_%s.pdf'%(specified_location[0], specified_location[1], time ))

print "\n*** Finished ***"
#plt.show()

#============storage of the most important variables to store.pckl=============

f = open('data_sfgrids.pckl', 'w')
pickle.dump((chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, dimension_1, dimension_2, dim1, dim2, delta_dim_1, delta_dim_2, geometry, dim1_min, dim1_max, dim2_min, dim2_max, exits, delta_mesh_exit), f)

f.close()
