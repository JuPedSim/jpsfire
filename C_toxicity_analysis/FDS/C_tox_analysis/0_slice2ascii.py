# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Created on Wed Jan 28 10:34:18 2015

@author: Benjamin Schroeder
"""
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

parser = argparse.ArgumentParser()
parser.add_argument("slice_quantity", type=str, help="quantity of the slicefile")
parser.add_argument("slice_dim", type=str, help="axis of the slicefile")
parser.add_argument("slice_coord", type=float, help="coordinate on the axis" )

#===================CHANGE PARAMETERS AS NECEASSARY=============================
# location of your local fds2ascii binary
fds2ascii_path = '/Applications/FDS/FDS6/bin/fds2ascii'
# Path pointing to the fire simulation directory
fds_path = os.path.join(os.getcwd()+'/../')

### FDS CHID:
chid = 'toxicity_analysis'

### SLICE QUANTITY and location:
### check if multiple quantities have been submitted from tox_analysis.sh
### check if multiple locations have been submitted from tox_analysis.sh
try:
    cmdl_args = parser.parse_args()
    quantity = cmdl_args.slice_quantity
    location = (cmdl_args.slice_dim, cmdl_args.slice_coord)

    print '\nProcessing FDS slice quantities and locations defined in tox_analysis.sh\n'
except:
    ### if not use single quantity defined explicitly in 0_slice2ascii.py
    ### if not use single loation defined explicitly in 0_slice2ascii.py
    print '\nProcessing single FDS slice quantity and location defined in 0_slice2ascii.py\n'
    ### change quantity here - if necessary
    quantity = 'CARBON DIOXIDE VOLUME FRACTION'
    ### change dimension and location of the slicefile here - if necessary
    location = ('Z', 2.0)

# Grid resolution in x, y and z
dx=0.25
dy=0.25
dz=0.25

# Parameters to be tunneled to fds2ascii:
# types: (slice = 2)
data_type=2
# data extend: (all = 1)
extend=1
# domain size: (not limited ='n')
domain_size="n"
# Interpolation times sequences; adjust arange values and/or window
t_start=0
t_stop=120
t_step=20
t_window=1            # interpolation duration prompted by fds2ascii

# Do you want to have plots produced? May be computaionally intensive depending
# on your FDS simulation extend!
plots = True

#===============================================================================

# BEGINNING OF THE AUTOMATIC PART - NO CHANGES NEEDED

#===============================================================================
config_file_name = os.path.join("config_fds2ascii.csv")

specified_location = (location[0].upper(), location[1])

t_low=np.arange(t_start,t_stop+1,t_step)
t_up=[i+t_window for i in t_low]

number_slices=1

output=quantity+'_'

smv_slcfs = []
id_slices = []
id_meshes = []
meshes = []

#==============================================================================
# Readout of FDS file:

fds = open('../'+chid+'.fds', 'U')
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

smv = open('../'+chid+'.smv', 'U')
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

if os.path.exists(os.path.join('../0_slice2ascii')):
    shutil.rmtree(os.path.join('../0_slice2ascii'))

if not os.path.exists('../0_slice2ascii'):
    print 'create directory "0_slice2ascii"'
    os.makedirs('../0_slice2ascii')

if not os.path.exists(os.path.join('../0_slice2ascii/%s_%.2f'%(specified_location[0], specified_location[1]))):
    os.makedirs(os.path.join('../0_slice2ascii/%s_%.2f'%(specified_location[0], specified_location[1])))

try:
    os.remove('data_slice2ascii.pckl')
except OSError:
    pass

#==============================================================================

for k, id_slice in enumerate(id_slices):

    data=[]

    for a,b in enumerate(t_low):
        data.extend([fds_path, chid, data_type, extend, domain_size,\
        t_low[a], t_up[a], number_slices, id_slice, '0_slice2ascii/%s_%.2f/%s_mesh_%i.txt'%(specified_location[0], specified_location[1], output+str(t_low[a]), id_meshes[k])])

    data = np.reshape(data, (len(t_low),-1))


    file=open(config_file_name, 'w')
    for line in data:
        for element in line:
            file.write(element+ ', ')
        file.write('\n')
    file.close()

    #Open Config File Object

    try:
        fh=open(config_file_name, 'U')
    except:
        print"!!! The Config File "+config_file_name+\
        " does not exist or the path defined in the script is incorrect. !!!"
        #exit()

    #Read file with csv module.
    data_array = csv.reader(fh)
    #Convert into List object
    config_lists = [list(sublist) for sublist in data_array]
    print str(len(config_lists))+" lines read in from "+config_file_name+"\n"

    def extract_ascii_data(input_data):
        # change to working directory, the first string in each row of the config file.
        #print input_data[0]
        os.chdir(input_data[0])

        print "Working Directory:",os.getcwd()
        # Open fds2ascii and create a seperate output and input stream.
        proc = subprocess.Popen(os.path.join(fds2ascii_path),
                               shell=True,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               )
        # Loop through input fields (columns) for the line passed in
        #from the fds2ascii_Config_File.csv file.
        for field in input_data[1:]:
            if field != '':
                print "Input:",field
                proc.stdin.write('%s\n' % field)
            else:
                pass
        proc.communicate()[0]

    # Pass each line in the config file to the extraction function.
    for data_row in config_lists:
        extract_ascii_data(data_row)
        time.sleep(0.5)

f = open('C_tox_analysis/data_slice2ascii.pckl', 'w+')
pickle.dump((chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, plots), f)
f.close()

print "*** Finished ***"
# Script based on 06.19.08 by bryan.klein@nist.gov
