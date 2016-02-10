# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 10:34:18 2015

@author: Benjamin Schroeder
"""
import subprocess
import csv
import os
import numpy as np
import glob
import shutil
import re
import time

### Directory, that consists the FDS data:
simulations = glob.glob(os.path.join('../osloer_mesh_test_1_dx20'))

### FDS CHID:
chid = 'osloer_mesh_test_1_dx20'
#quantity = 'SOOT EXTINCTION COEFFICIENT'
quantity = 'TEMPERATURE'
specified_location = ('x'.upper(), 9.2)
# Grid resolution in x, y and z
dx=0.2
dy=0.2
dz=0.2

# Parameters for the slice conversion:

path = os.path.join("../%s"%chid)
config_file_name = os.path.join("../%s/config_fds2ascii.csv"%chid)

# types: slice = 2
data_type=2

# data extend: all = 1
extend=1

# domain size: not limited ='n'
domain_size="n"

# Interpolation times sequences; adjust arange values and/or window
start=0
stop=300
step=10
window=5            # interpolation duration prompted by fds2ascii

t_low=np.arange(start,stop+1,step)
t_up=[i+window for i in t_low]

number_slices=1

output="TEMP_"

#==============================================================================

for simulation in simulations:

    smv_slcfs = []
    id_slices = []
    id_meshes = []


#==============================================================================
    # Readout of SMV file:

    smv = open(simulation+'/'+chid+'.smv', 'U')
    smv_entries = smv.readlines()[-300:]

    for i, smv_entry in enumerate(smv_entries):
        if smv_entry[0:4]=='SLCF':
            smv_slcfs=np.append(smv_slcfs, re.split('&|\n', smv_entry) +\
            re.split('&|\n', smv_entries[i+2]))

    smv_slcfs = np.reshape(smv_slcfs, (-1,5))

    for i, smv_slcf in enumerate(smv_slcfs):
        coordinates = smv_slcf[1].split()
        if coordinates[0] == coordinates[1]:
                location = ('X' , round(int(coordinates[0]) * dx, 2))
        elif coordinates[2] == coordinates[3]:
                location = ('Y' , round(int(coordinates[2]) * dy, 2))
        else:
                location = ('Z' , round(int(coordinates[4]) * dz, 2))

        if quantity in smv_slcf[-2] and location == specified_location:
            id_slices = np.append(id_slices, i+1)

            mesh = smv_slcf[0].split()[1]
            id_meshes = np.append(id_meshes, int(mesh))

        else:
            print 'geht nicht'


#==============================================================================


    if os.path.exists(os.path.join('../%s/ascii'%chid)):
        shutil.rmtree(os.path.join('../%s/ascii'%chid))


    if not os.path.exists('%s/ascii'%simulation):
        print 'create directory "ascii"'
        os.makedirs('%s/ascii/'%simulation)


    if os.path.exists(os.path.join('../%s/ascii/%s'%(chid, str(specified_location[0])+'_'+str(specified_location[1])))):
        continue
    else:
        os.makedirs(os.path.join('../%s/ascii/%s'%(chid, str(specified_location[0])+'_'+str(specified_location[1]))))


    for k, id_slice in enumerate(id_slices):


        data=[]

        for a,b in enumerate(t_low):
            data.extend([path, chid, data_type, extend, domain_size,\
            t_low[a], t_up[a], number_slices, id_slice, 'ascii/%s/'\
            %(str(specified_location[0])+'_'+str(specified_location[1]))+output+str(t_low[a])+'_mesh_%i'\
            %id_meshes[k]+'.txt'])

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
            proc = subprocess.Popen(os.path.join('/Applications/FDS/FDS6/bin/fds2ascii'),
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

    print "*** Finished ***"
    # Script based on 06.19.08 by bryan.klein@nist.gov
