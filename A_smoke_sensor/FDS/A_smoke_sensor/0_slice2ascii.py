# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Created on Wed Jan 28 10:34:18 2015
@author: Benjamin Schroeder
# Script based on 06.19.08 by bryan.klein@nist.gov
"""
__author__ = 'Benjamin Schroeder'
__email__ = 'dev@jupedsim.org'
__copyright__ = '<2009-2017> Forschungszentrum Jülich GmbH. All rights reserved.'
__license__ = 'GNU Lesser General Public License'
__version__ = '0.1'
__status__ = 'Production'


import subprocess
import csv
import os
import numpy as np
import shutil
import re
import time
import sys
import pickle
import argparse
import logging
import glob
import xml.etree.ElementTree as ET

basename = os.path.basename(__file__)
logfile =  os.path.join(os.path.dirname(__file__), "log_%s.txt" % basename.split(".")[0])

open(logfile, 'w').close()
logging.basicConfig(filename=logfile, level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')
print("Running <%s>. (Logfile: <%s>)" % (__file__, logfile))

def get_tstop(fds_path, chid):
    fds_file = os.path.join(fds_path, chid+".fds")
    logging.info("enter get_tstop with file <%s>" % fds_file)
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

def get_gridres(chid):

    fds_data = open('../'+chid+'.fds')
    fds_lines = fds_data.readlines()

    for l in fds_lines:
        if l.startswith("&MESH"):
            l = l.split('IJK=')[-1]
            l = re.split(r'[,=/]+', l)
            return (float(l[5])-float(l[4]))/float(l[0]), (float(l[7])-float(l[6]))/float(l[1]), (float(l[9])-float(l[8]))/float(l[2])

def get_tstep(jps_path):
    geo = glob.glob(jps_path+'*ini*.xml')
    tree = ET.parse(geo[0])
    root = tree.getroot()

    for time in root.iter('A_smoke_sensor'):
            update_time = time.attrib.get('update_time')
    return float(update_time)

def getParserArgs():
    parser = argparse.ArgumentParser(description="Automation of NIST's fds2ascii. Conversion of FDS-data (slicefiles) to ascii format")
    parser.add_argument("-q", "--slice_quantity", type=str, default = "OPTICAL DENSITY", help="quantity of the slicefile (default: OPTICAL DENSITY)", required=False)
    parser.add_argument("-d", "--slice_dim", type=str, default='Z', help="axis of the slicefile (default: Z)", required=False)
    parser.add_argument("-c", "--slice_coord", type=float, default=2.25, help="coordinate on the axis (default: 2.25)", required=False)
    parser.add_argument("-p", "--plot", dest="plot", action='store_true')
    parser.add_argument("-f", "--fds", type=str, default='/Applications/FDS/FDS6/bin/fds2ascii', help="Absolute path to fds2ascii", required=False)
    parser.add_argument("-j", "--jps", type=str, default='../../JuPedSim/', help="Path pointing to the JuPedSim simulation directory", required=False)
    parser.add_argument("-s", "--start", type=float, default=0, help="start time for fds input (default: 0)", required=False)
    parser.add_argument("-e", "--end", type=float, default=get_tstop(fds_path, chid), help="end time for fds input (default: T_END from fds file)", required=False)
    
    args = parser.parse_args()
    return args

# Path pointing to the fire simulation directory
script_dir = os.path.dirname(os.path.realpath(__file__))
fds_path = os.path.join(script_dir, "..")
logging.info('Processing FDS slice quantities and locations defined in smoke_sensor.sh')
logging.info("script path: <%s>" % script_dir)
logging.info("fds path: <%s>" % fds_path)

# FDS CHID:
chid_fds = glob.glob('../*.fds')
chid = (chid_fds[0].strip('.fds')).strip('../')

# Parse parameters
cmdl_args = getParserArgs()

quantity = cmdl_args.slice_quantity
logging.info("Quantity: " + quantity)

location = (cmdl_args.slice_dim, cmdl_args.slice_coord)
logging.info("location: " +  ", ".join(map(str, location)))

fds2ascii_exec = cmdl_args.fds
logging.info("fds2ascii_exec: %s" % fds2ascii_exec)

jps_path = cmdl_args.jps
logging.info("jps_path: %s" % jps_path)

plots = cmdl_args.plot
logging.info("Plot: On" if plots else "Plot: Off")

# Grid resolution in x, y and z
dx, dy, dz = get_gridres(chid)
logging.info('Grid resolution is dx = %.2f, dy = %.2f, dz = %.2f' % (dx, dy, dz))

# Parameters to be tunneled to fds2ascii:
# types: (SLCF file = 2)
data_type=2
# data extend: (all data = 1)
extend=1
# domain size: (not limited ='n')
domain_size="n"
# Interpolation time sequences; adjust arrange values and/or window
t_start= cmdl_args.start
logging.info("t_start: %.2f", t_start)
t_stop = cmdl_args.end
logging.info("t_stop: %.2f", t_stop)
t_step=get_tstep(jps_path) #at the moment read in from ini
logging.info("t_step: %s" % t_step)

#interpolation duration prompted by fds2ascii
t_window=1

#===============================================================================
# BEGINNING OF THE AUTOMATIC PART - NO CHANGES NEEDED
if not os.path.exists(fds2ascii_exec):
    logging.critical("<%s> does not exist. exits." %fds2ascii_exec)
    sys.exit("<%s> does not exist. exits." %fds2ascii_exec)

#===============================================================================
# config_file_name = os.path.join("config_fds2ascii.csv")
config_file_name = os.path.join(script_dir, "config_fds2ascii.csv")

specified_location = (location[0].upper(), location[1])

t_low=np.arange(t_start,t_stop+1,t_step)
t_up=[i+t_window for i in t_low]
number_slices=1

smv_slcfs = []
id_slices = []
id_meshes = []
meshes = []

#==============================================================================
# Readout of FDS file:

fds = open('../'+chid+'.fds', newline=None)
fds_entries = fds.readlines()

for i, fds_entry in enumerate(fds_entries):

    # Detection of MESH entries:
    if fds_entry[0:5]=='&MESH':
        fds_entry = fds_entry.split('XB')[-1]
        fds_entry = re.split(r'[,=/]+', fds_entry)
        fds_entry = np.array(fds_entry[1:-1], dtype=float)
        meshes = np.append(meshes, fds_entry)

meshes = np.reshape(meshes, (-1, 6))

#==============================================================================
# Readout of SMV file:
smv = open('../'+chid+'.smv', newline=None)
logging.info('open ../'+chid+'.smv')
smv_entries = smv.readlines()

for i, smv_entry in enumerate(smv_entries):

    # Detection of SLCF entries:
    if smv_entry[0:4]=='SLCF':
        smv_slcfs=np.append(smv_slcfs, re.split('&|\n', smv_entry) +\
        re.split('&|\n', smv_entries[i+2]))

smv_slcfs = np.reshape(smv_slcfs, (-1,5))

for i, smv_slcf in enumerate(smv_slcfs):
    mesh = int(smv_slcf[0].split()[1])
    coordinates = smv_slcf[1].split()

    if coordinates[0] == coordinates[1]:
            location = ('X' , round( meshes[mesh-1, 0] + int(coordinates[0]) * dx, 2))
    elif coordinates[2] == coordinates[3]:
            location = ('Y' , round( meshes[mesh-1, 2]+ int(coordinates[2]) * dy, 2))
    else:
        location = ('Z' , round(meshes[mesh-1, 4] + int(coordinates[4]) * dz, 2))

    if quantity in smv_slcf[-2] and specified_location == location:
        id_slices = np.append(id_slices, i+1)
        id_meshes = np.append(id_meshes, mesh)

if len(id_slices)==0:
    sys.exit('Mismatch between specified SLICE locations and/or quantities! -> EXIT')

#==============================================================================

quantity = quantity.replace(' ', '_')

if os.path.exists(os.path.join('../0_slice2ascii')):
    shutil.rmtree(os.path.join('../0_slice2ascii'))


logging.info('create directory <0_slice2ascii>')
os.makedirs('../0_slice2ascii')
Z_directory = os.path.join('../0_slice2ascii', '%s_%.2f'%(specified_location[0], specified_location[1]))
if not os.path.exists(Z_directory):
    os.makedirs(Z_directory)
    logging.info("create directory <%s>" % Z_directory)

if os.path.exists('data_slice2ascii.pckl'):
    os.remove('data_slice2ascii.pckl')

#==============================================================================

#for l, id_mesh in enumerate(id_meshes):
for k, id_slice in enumerate(id_slices):

    data=[]

    for a,b in enumerate(t_low):
        data.extend([fds_path, chid, data_type, extend, domain_size,\
        t_low[a], t_up[a], number_slices, id_slice, '0_slice2ascii/%s_%.2f/%s_%i_mesh_%i.txt'%(specified_location[0], specified_location[1], quantity, t_low[a], id_meshes[k])])

    data = np.reshape(data, (len(t_low),-1))


    file=open(config_file_name, 'w')
    for line in data:
        for element in line:
            file.write(element+ ', ')
        file.write('\n')
    file.close()

    #Open Config File Object

    try:
        fh=open(config_file_name, newline=None)
    except:
        logging.info("!!! The Config File %s does not exist or the path defined in the script is incorrect."% config_file_name)
        #exit()

    #Read file with csv module.
    data_array = csv.reader(fh)
    #Convert into List object
    config_lists = [list(sublist) for sublist in data_array]
    logging.info("%d lines read in from %s"% (len(config_lists), config_file_name))
    def extract_ascii_data(input_data):
        # change to working directory, the first string in each row of the config file.
        os.chdir(input_data[0])

        logging.info("Working Directory: %s"% os.getcwd())
        # Open fds2ascii and create a seperate output and input stream.
        proc = subprocess.Popen(fds2ascii_exec,
                                shell=True,
                                stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                universal_newlines=True
        )
        # Loop through input fields (columns) for the line passed in
        #from the fds2ascii_Config_File.csv file.
        logging.info("Input: %s"% " ".join(map(str, input_data[1:])))
        for field in input_data[1:]:
            field_d = field.strip()
            if field_d:
                # logging.debug("field: <%s> " %field)
                proc.stdin.write('%s\n' % field)
            else:
                # logging.debug("not field: <%s> " %field)
                pass

        proc.communicate("\n")[0]


    for data_row in config_lists:
        extract_ascii_data(data_row)
        time.sleep(0.5)

pckl_file = os.path.join('A_smoke_sensor', 'data_slice2ascii.pckl')
f = open(pckl_file, 'wb')
#print ((chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, jps_path, plots))
pickle.dump((chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, jps_path, plots, dx), f)
f.close()

logging.info("*** Finished ***")
