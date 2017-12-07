# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Created on Wed Jan 28 10:34:18 2015
@author: Benjamin Schroeder
# Script based on 06.19.08 by bryan.klein@nist.gov
"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec
import re
import time
import argparse
import logging
import glob
import xml.etree.ElementTree as ET
import collections
import math
import scipy.ndimage
import warnings
import sys
sys.path.append("../..")
import fdsreader.slice as fs

__author__ = 'Benjamin Schroeder'
__email__ = 'dev@jupedsim.org'
__copyright__ = '<2009-2017> Forschungszentrum Jülich GmbH. All rights reserved.'
__license__ = 'GNU Lesser General Public License'
__version__ = '0.1'
__status__ = 'Production'

rcParams.update({'figure.autolayout': True})
font = {'family': 'serif',
        'weight': 'normal',
        'size': 9}
matplotlib.rc('font', **font)

basename = os.path.basename(__file__)
logfile = os.path.join(os.path.dirname(__file__), "log_%s.txt" % basename.split(".")[0])
open(logfile, 'w').close()
logging.basicConfig(filename=logfile, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
print("Running <%s>. (Logfile: <%s>)" % (__file__, logfile))

def get_tstop(fds_path, chid):
    fds_file = os.path.join(fds_path, chid+".fds")
    with open(fds_file) as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith("&TIME"):
                time_line = l.strip("/ \n").split("=")
            else:
                continue

            if len(time_line) > 1:
                t_stop = float(time_line[1])
            else:
                logging.critical("could not read t_stop from file <%s>" % fds_file[0])
                sys.exit("ERROR: can not read time from fds-file")
    return t_stop
def get_dims(chid, spec):

    fds_data = open('../'+chid+'.fds')
    fds_lines = fds_data.readlines()

    for l in fds_lines:
        if l.startswith("&MESH"):
            l = l.split('IJK=')[-1]
            l = re.split(r'[,=/]+', l)
            l.remove(' XB'), l.remove('\n')
            l = np.array(l[:], dtype=float)
            if spec == 'coordinates':
                return l[3], l[4], l[5], l[6]
            elif spec == 'grids':
                return (l[4]-l[3])/l[0], (l[6]-l[5])/l[1], (l[8]-l[7])/l[2]
            else:
                return 'Error, wrong specification' #todo: logging
def get_tstep(jps_path):
    geo = glob.glob(jps_path+'*ini*.xml')
    tree = ET.parse(geo[0])
    root = tree.getroot()

    for time in root.iter('A_smoke_sensor'):
            update_time = time.attrib.get('update_time')
    return float(update_time)
def getParserArgs():
    parser = argparse.ArgumentParser(description="Read out of SmokeView slicefiles and conversion to grid")
    parser.add_argument("-q", "--slice_quantity", type=str, default="EXTINCTION COEFFICIENT", help="quantity of the slicefile (default: EXTINCTION COEFFICIENT)", required=False)
    parser.add_argument("-d", "--slice_dim", type=str, default='Z', help="axis of the slicefile (default: Z)", required=False)
    parser.add_argument("-c", "--slice_coord", type=float, default=2.25, help="coordinate on the axis (default: 2.25)", required=False)
    parser.add_argument("-p", "--plot", dest="plot", action='store_true')
    parser.add_argument("-j", "--jps", type=str, default='../../JuPedSim/', help="Path pointing to the JuPedSim simulation directory", required=False)
    parser.add_argument("-s", "--start", type=float, default=0.0, help="start time for fds input (default: 0.0)", required=False)
    parser.add_argument("-e", "--end", type=float, default=get_tstop(fds_path, chid), help="end time for fds input (default: T_END from fds file)", required=False)
    parser.add_argument("-g", "--delta_sfgrid", type=float, default=1.0, help="Resolution of smoke factor grids (default: 1.0)", required=False)

    args = parser.parse_args()
    return args
def cm2inch(value):
    return value/2.54
def isfloat(value):
      try:
        float(value)
        return True
      except ValueError:
        return False
def read_fds_line(fds_entry):

        coordinates = [x.strip() for x in fds_entry.split('XB')][1:]
        coordinates = ''.join(coordinates)
        coordinates = re.split(', |,| |=|/', coordinates)
        for j, cord in enumerate(coordinates):
            if isfloat(cord) == True:
                coordinates[j] = float(cord)
            if cord.isdigit() == True:
                coordinates[j] = float(cord)
            else:
                continue

        coordinates = [x for x in coordinates if type(x) != str]
        coordinates = list(map(float, coordinates))
        # increasing sorting of coordinate pairs:
        coordinates = sorted(coordinates[0:2])+sorted(coordinates[2:4])+sorted(coordinates[4:6])
        return coordinates
def get_fds_geo(chid):

    geometry = np.zeros((len(dim_y),len(dim_x)))

    obsts = []
    holes = []

    fds = open('../'+chid+'.fds', newline=None)
    logging.info('Opening '+chid+'.fds'+' to read out meshes, obstacles and holes')
    fds_entries = fds.readlines()

    for i, fds_entry in enumerate(fds_entries):
        if fds_entry[0:5] == '&OBST':
            obst = read_fds_line(fds_entry)
            if obst[4] < specified_location[1] < obst[5]:
                obst = [(1/dx)*i for i in obst]
                geometry[int(obst[2]-y_max/dy-1): int(obst[3]-y_max/dy-1),
                int(obst[0]-x_max/dx-1): int(obst[1]-x_max/dx-1)][:] = np.nan
                obsts = np.append(obsts,obst)

    obsts = np.reshape(obsts, (-1, 6))
    #np.savetxt('obst.csv', obsts, delimiter=',')

    for i, fds_entry in enumerate(fds_entries):
        if fds_entry[0:5] == '&HOLE':
            hole = read_fds_line(fds_entry)
            if hole[4] < specified_location[1] < hole[5]:
                hole = [(1/dx)*i for i in hole]
                geometry[hole[2]-y_max/dy-1: hole[3]-y_max/dy-1, hole[0]-x_max/dx-1: hole[1]-x_max/dx-1][:] = 0
                holes = np.append(holes, hole)

    holes = np.reshape(holes, (-1, 6))
    #np.savetxt('../1_meshgrid/hole.csv', holes, delimiter=',')

    fds.close()
    #np.savetxt('../geometry.txt', geometry)
    return geometry
def get_exits(jps_path):

    exits_dict = {}
    crossings = []
    transitions = []

    try:
        geo = glob.glob(jps_path+'*eo*.xml')
        tree = ET.parse(geo[0])
        root = tree.getroot()

        for crossing in root.iter('crossing'):
            c_id = 'cross_' + crossing.attrib.get('id')
            for vert in crossing.iter('vertex'):
                px = vert.get('px')
                py = vert.get('py')
                crossings.append([c_id, px,py])

        crossings = np.reshape(crossings, (-1,3))

        for k,cross in enumerate(crossings):
            if k >= len(crossings)-1:
                break
            if crossings[k,0]==crossings[k+1,0]:
                x = (float(crossings[k,1]) + float(crossings[k+1,1]))/2.
                y = (float(crossings[k,2]) + float(crossings[k+1,2]))/2.
                if crossings[k,1]==crossings[k+1,1]:
                    orientation = 'y'
                else:
                    orientation = 'x'
                exits_dict[cross[0]] = x,y,orientation

        for transition in root.iter('transition'):
            t_id = 'trans_' + transition.attrib.get('id')
            for vert in transition.iter('vertex'):
                px = vert.get('px')
                py = vert.get('py')
                transitions.append([t_id, px,py])

        transitions = np.reshape(transitions, (-1,3))

        for l,trans in enumerate(transitions):
            if l >= len(transitions)-1:
                break
            if transitions[l,0]==transitions[l+1,0]:
                x = (float(transitions[l,1]) + float(transitions[l+1,1]))/2.
                y = (float(transitions[l,2]) + float(transitions[l+1,2]))/2.
                if transitions[l,1]==transitions[l+1,1]:
                    orientation = 'y'
                else:
                    orientation = 'x'
                exits_dict[trans[0]] = x,y,orientation
                # format of the exits_dict: exits_dict={"A":(-0.3,2.0,'y')}

        return collections.OrderedDict(sorted(exits_dict.items()))

    except:
        exits_dict = {}
        return logging.info('!!! WARNING No geometry file found - no exits extracted !!!')
def smoke_factor_conv(convert, exit):

    if not os.path.exists('../3_sfgrids/dx_%.2f/%s_%.6f/Door_X_%.6f_Y_%.6f' % (delta_smoke_factor_grid,
                                                                             specified_location[0],
                                                                             specified_location[1],
                                                                             exits[exit][0], exits[exit][1])):
        os.makedirs('../3_sfgrids/dx_%.2f/%s_%.6f/Door_X_%.6f_Y_%.6f' % (delta_smoke_factor_grid, specified_location[0],
                                                                       specified_location[1],  exits[exit][0],
                                                                       exits[exit][1]))

    smoke_factor_grid = np.ones([int(abs(y_max-y_min)/delta_smoke_factor_grid),
                                 int(abs(x_max-x_min)/delta_smoke_factor_grid)])

    x, y = m_to_pix(x_i=exits[exit][0], y_i=exits[exit][1])

    ### calculation of the global max of the extracted field_quantity
    ### (without geometry)
    #global_D_max = np.amax(magnitudes-geometry[ :-1 , :-1])
    #print global_D_max

    for a, line in enumerate(smoke_factor_grid):

        for b, col in enumerate(line):

            x0, y0 = m_to_pix(b*delta_smoke_factor_grid+x_min, a*delta_smoke_factor_grid+y_min)    # virtual agent position

            #### stop of iteration and edge factor storage of predecessor cell
            #### if virtual agent position is identical with exit position
            if x0 == x and y0 == y:
                smoke_factor = smoke_factor_grid[a, b-1]
                smoke_factor_grid[a, b] = smoke_factor
                continue

            x_exit, y_exit = np.linspace(x0, x, math.hypot(x - x0, y - y0)), np.linspace(y0, y, math.hypot(x - x0, y - y0))

            magnitude_along_line_of_sight = scipy.ndimage.map_coordinates(np.transpose(convert), np.vstack((x_exit, y_exit)))

            #### consideration of visibility straight lines whose lengt is lt 1:
            if len(magnitude_along_line_of_sight) < 1:
                continue

            else:
                #### Computation of the smoke factor towards the exit:
                smoke_factor = np.nanmax(magnitude_along_line_of_sight)*abs(np.trapz(magnitude_along_line_of_sight, dx=dx))
                #todo: smoke factor has a unit!!!
            #### storage of the edge factor
            smoke_factor_grid[a, b] = smoke_factor

    #if np.amax(smoke_factor_grid)>max_Smoke_Factor:
    max_Smoke_Factor = np.amax(smoke_factor_grid)

    #print exit, 'maximum Smoke Factor', max_Smoke_Factor

    # Norm of the smoke factor grid multiplied with 10 in order to achieve a
    # factor range from 0..10
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        smoke_factor_grid_norm = smoke_factor_grid/max_Smoke_Factor*10


    header='Room No. 1 , Exit %s, \n dX[m], dY[m] , minX[m] , maxX[m], minY[m], maxY[m] \n   %f  ,  %f    ,  %f  ,  %f  ,  %f ,  %f' \
    %(exit, delta_smoke_factor_grid, delta_smoke_factor_grid, x_min, x_max, y_min, y_max)
    np.savetxt('../3_sfgrids/dx_%.2f/%s_%.6f/Door_X_%.6f_Y_%.6f/t_%.0f.000000.csv'\
    %(delta_smoke_factor_grid, specified_location[0], specified_location[1],  exits[exit][0], exits[exit][1], time), smoke_factor_grid_norm, header=header, delimiter=',', comments='')
    logging.info('Write smoke factor grid: ../3_sfgrids/dx_%.2f/%s_%.6f/Door_X_%.6f_Y_%.6f/t_%.0f.000000.csv'\
    %(delta_smoke_factor_grid, specified_location[0], specified_location[1],  exits[exit][0], exits[exit][1], time))

    return a,b, x0, y0, x, y, magnitude_along_line_of_sight, smoke_factor, time
def m_to_pix(x_i, y_i):
    # Conversion m to pixel coordinates
    return (abs(x_max-x_min) - x_max + x_i)/dx, (abs(y_max-y_min) - y_max + y_i)/dy

# Path pointing to the fire simulation directory
script_dir = os.path.dirname(os.path.realpath(__file__))
fds_path = os.path.join(script_dir, "..")
logging.info("script path: <%s>" % script_dir)
logging.info("fds path: <%s>" % fds_path)

# FDS chid
fds_file = glob.glob('../*.fds')
chid = (fds_file[0].strip('.fds')).strip('../')

# Parse parameters
cmdl_args = getParserArgs()

jps_path = cmdl_args.jps
logging.info("jps_path: %s" % jps_path)

quantity = cmdl_args.slice_quantity
quantity = quantity.replace(' ', '_')
logging.info("Quantity: %s" % quantity)

specified_location = (cmdl_args.slice_dim, cmdl_args.slice_coord) # todo: is this actually necessary? maybe also read from fds
logging.info("Specified location: %s %s" % specified_location)

t_start = cmdl_args.start
logging.info("t_start: %.1f" % t_start)
t_stop = cmdl_args.end
logging.info("t_stop: %.1f" % t_stop)
t_step = get_tstep(jps_path)
logging.info("t_step: %.1f" % t_step)

plots = cmdl_args.plot
logging.info("Plot: On" if plots else "Plot: Off")

# Get spatial extend
x_min, x_max, y_min, y_max = get_dims(chid, 'coordinates')
logging.info('FDS coordinates are xmin = %.2f, xmax = %.2f, ymin = %.2f, ymax = %.2f' % (x_min, x_max, y_min, y_max))

# Grid resolution in x, y and z
dx, dy, dz = get_dims(chid, 'grids')
logging.info('FDS grid resolution is dx = %.2f, dy = %.2f, dz = %.2f' % (dx, dy, dz))

dim_x = np.arange(x_min, x_max+0.01, dx)
dim_y = np.arange(y_min, y_max+0.01, dy)

# ==============================================================================
# Readout of obstacles and holes from .fds file
# ==============================================================================

geometry = get_fds_geo(chid)

# ==============================================================================
# Readout of slicefiles from .smv file using fdsreader
# ==============================================================================

# locate smokeview file
root_dir = "../../FDS"
smv_fn = fs.scanDirectory(root_dir)
logging.info("smv file found: %s" % smv_fn)

# parse smokeview file for slice information
sc = fs.readSliceInfos(os.path.join(root_dir, smv_fn))
# print all found information
# sc.print()

# read in meshes
meshes = fs.readMeshes(os.path.join(root_dir, smv_fn))

# select matching slice
slice_label = 'ext_coef_C0.9H0.1'
sid = -1
for iis in range(len(sc.slices)):
    if sc[iis].label == slice_label:
        sid = iis
        logging.info("found matching slice")
        break

if sid == -1:
    logging.info("no slice matching label: {}".format(slice_label))
    sys.exit()

slice = sc[sid]

# read in time information
slice.readAllTimes(root_dir)

# read in slice data
slice.readTimeSelection(root_dir, dt=t_step, average_dt=1)

# map data on mesh
slice.mapData(meshes)

# get max value
max_coefficient = 0
for it in range(0, slice.times.size):
    cmax = np.max(slice.sd[it])
    max_coefficient = max(cmax, max_coefficient)

if not os.path.exists('../slicedata'): os.makedirs('../slicedata')
logging.info('create directory <slicedata>')

Z_directory = os.path.join('../slicedata', '%s_%.2f' % (specified_location[0], specified_location[1]))

if not os.path.exists(Z_directory):
    os.makedirs(Z_directory)
logging.info("create directory <%s>" % Z_directory)

# plot slice data
for it in range(0, slice.times.size):  # todo: is it possible to get defined times? (@Lukas) warum geo unterschiedlich
    np.savetxt('../slicedata/%s_%.2f/%s_%s.txt' % (specified_location[0], specified_location[1], quantity, slice.times[it]), slice.sd[it])

    if plots == True:
        collect = slice.sd[it]
        plt.imshow(slice.sd[it], cmap='coolwarm', vmax=max_coefficient,
                   origin='lower', extent=slice.sm.extent)
        plt.title("time = {:.2f}".format(slice.times[it]))
        plt.colorbar(label="{} [{}]".format(slice.quantity, slice.units))
        plt.savefig("../slicedata/%s_%.2f/single_slice_%s.pdf" % (specified_location[0], specified_location[1], slice.times[it]))
        plt.clf()
        plt.imshow(geometry, cmap='coolwarm', vmax=max_coefficient,
                   origin='lower', extent=slice.sm.extent)
        plt.title("time = {:.2f}".format(slice.times[it]))
        plt.colorbar(label="{} [{}]".format(slice.quantity, slice.units))
        plt.savefig("../slicedata/%s_%.2f/single_slice_geo_%s.pdf" % (specified_location[0], specified_location[1], slice.times[it]))
        plt.clf()

# todo: plot geometry too!

# ==============================================================================
# Readout of crossings and transitions from jps geometry (saved as dictionary)
# ==============================================================================

exits = get_exits(jps_path)

# ==============================================================================
# Calculation of smoke factor grids
# Default resolution of smoke factor grids: delta_smoke_factor_grid = 1 (m)
# ==============================================================================

delta_smoke_factor_grid = cmdl_args.delta_sfgrid
logging.info('Generation of smoke factor grids with mesh resolution: %s [m]' % delta_smoke_factor_grid)

for it in range(0, slice.times.size):

    convert = slice.sd[it]
    time = slice.times[it]

    logging.info('-----------------------------------------------------------------')
    logging.info('Processing files for time: %s' % time)

    for id, exit in enumerate(exits):

        a, b, x0, y0, x, y, magnitude_along_line_of_sight, smoke_factor, time = smoke_factor_conv(convert, exit)

        if plots == True:

            # plot line of sights

            ## ==== adjust here for debugging =====

            x_0 = 12.5  # todo: parse both
            y_0 = 5.5

            #slicefile = '../slicedata/Z_2.25/OPTICAL_DENSITY_%i.txt' %(time)
            sfgrid = '../3_sfgrids/dx_1.00/Z_2.250000/Door_X_25.000000_Y_6.000000/t_%.f.000000.csv' % time
            exit = 'trans_0' #e.g. exit with smoke

            ### ==== automatic part ======

            x0 = x_0 / dx
            y0 = y_0 / dx

            x, y = m_to_pix(x_i=exits[exit][0], y_i=exits[exit][1])
            x_exit, y_exit = np.linspace(x0, x, math.hypot(x - x0, y - y0)), np.linspace(y0, y, math.hypot(x - x0, y - y0))

            magnitude_along_line_of_sight = scipy.ndimage.map_coordinates(np.transpose(convert), np.vstack((x_exit,y_exit)))

            fig = plt.figure()

            gs = gridspec.GridSpec(1, 40)

            ax1 = fig.add_subplot(gs[0,:20])
            ax2 = fig.add_subplot(gs[0,15:22])
            ax3 = fig.add_subplot(gs[0,22:])
            ax1.set_xlabel('x [m]')
            ax1.set_ylabel('y [m]')
            ax1.set_xticks(np.arange(x_min, x_max+1, 5))
            ax1.set_yticks(np.arange(y_min, y_max+1, 5))

            # Plot the line of sights towards each exit
            ax1.plot([x_0, exits[exit][0]], [y_0, exits[exit][1]], lw=1.5, ls=':', color='red', label='line of sight')
            ax1.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05), ncol=2, frameon=False)
            aa = ax1.pcolorfast(dim_x, dim_y, convert, cmap='Greys', vmin=0, vmax=2)
            ax1.minorticks_on()
            ax1.axis('image')

            ax1.grid(which='major', lw=0.5, alpha=0.5)
            ax1.grid(which='minor', lw=0.5, alpha=0.5)

            sfgrid = np.loadtxt(sfgrid, skiprows=3, delimiter=',')

            for i in range(int((1/dx)/2)):
                sfgrid = np.kron(sfgrid, [[1, 1], [1, 1]])

            smoke_factor = sfgrid[int(y0/delta_smoke_factor_grid), int(-x0/delta_smoke_factor_grid)]

            ax3.plot(magnitude_along_line_of_sight,  lw=1.5, color='black', label='%s: $f_{smoke}$ = %.2f' % (exit, smoke_factor))

            ax3.set_xlabel('l [m]')
            labels = ax3.get_xticks()
            labels = (labels*dx).astype(int)
            ax3.set_xticklabels(labels)

            ax3.set_ylabel('Extinction Coefficient')
            ax3.set_ylim(0, 2)

            fig.colorbar(aa, ax=ax1, cax=ax2, orientation='vertical')

            ax3.legend(loc='upper left')
            ax3.grid(ls='--', lw=0.5)

            plt.savefig('../3_sfgrids/dx_%.2f/%s_%.6f/%s_%s_%.f.pdf' % (delta_smoke_factor_grid, specified_location[0],
                                                                    specified_location[1], quantity+'_debug', exit, time))
            plt.close()

            #########################
            # Plot Smoke factor grids
            #########################

            # todo: only trans
            fig = plt.figure()
            nrows = int(math.ceil(len(exits)**0.5))
            print(nrows)
            ncols = int(math.ceil(len(exits)**0.5))

            #fig, axes = plt.subplots(nrows, ncols)
            gs = gridspec.GridSpec(nrows+1, ncols)

            for i, g in enumerate(gs):

                if i == len(exits):
                    cbar_ax = plt.subplot(gs[-1,:])
                    fig.colorbar(aa, cax=cbar_ax, label=r'$f_{smoke}$', orientation='horizontal')
                    break

                exit = list(exits.items())[i][0]
                ax = plt.subplot(g)
                plt.xlabel('x [m]')
                plt.ylabel('y [m]')
                plt.tight_layout()
                smoke_factor_grid_norm = np.loadtxt('../3_sfgrids/dx_%.2f/%s_%.6f/Door_X_%.6f_Y_%.6f/t_%.f.000000.csv'%
                                                    (delta_smoke_factor_grid, specified_location[0], specified_location[1],
                                                     exits[exit][0], exits[exit][1], time), delimiter=',', skiprows=3)

                if time == 0.25: #todo: workaround as starter slice is filled with NAN's!!!
                    continue
                aa = ax.pcolorfast(dim_x, dim_y, smoke_factor_grid_norm, cmap='jet', vmin=0, vmax=10)
                ax.set_title(exit)
                ax.set_aspect('equal')
                ax.set_xticks(np.arange(x_min, x_max+1, 5))
                ax.set_yticks(np.arange(y_min, y_max+1, 5))
                ax.minorticks_on()
                ax.grid(which='major', linestyle='-', lw=0.5, alpha=0.4)
                ax.grid(which='minor', linestyle='-', lw=0.5, alpha=0.4)

            plt.savefig('../3_sfgrids/dx_%.2f/%s_%.6f/sfgrids_%i.pdf'%(delta_smoke_factor_grid, specified_location[0], specified_location[1], time))
            plt.close()
            logging.info('Plot smoke factor grid: ../sfgrids/dx_%.2f/%s_%.2f/sfgrid_%i.pdf'%(delta_smoke_factor_grid, specified_location[0], specified_location[1],  time))

#dimension_1 = x
#dimension_2 = y
#dim1 = dim_x
#dim2 = dim_y
#dim1_min = x_min
#dim1_max = x_max
#dim2_min = y_min
#dim2_max = y_max

logging.info("*** Finished ***")
