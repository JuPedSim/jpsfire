# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 10:34:18 2015

@author: Benjamin
"""
import subprocess
import csv
import os
import numpy as np
import glob
import shutil

simulations = glob.glob(os.path.join('../offenbach'))
for simulation in simulations:
    print simulation

    if not os.path.exists('%s/ascii'%simulation):
        print 'create directory "ascii"'
        os.makedirs('%s/ascii/'%simulation)


#==============================================================================
    # Parameters for the slice conversion:

    chid = simulation[simulation.index('/')+1:]

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
    stop=400
    step=60
    window=1            # interpolation duration prompted by fds2ascii

    t_low=np.arange(start,stop+1,step)
    t_up=[i+window for i in t_low]

    number_slices=1
    slices_dict = {16:"2.0", 14:"2.2", 12:"2.4", 10:"2.6", 6:"2.8", 2:"3.0"}

    id_slices=[16,14,12,10,6,2][4:5]    # SOOT_DENSITY
    #id_slices=[39, 80, 121, 162, 203]    # SOOT_DENSITY
    #id_slices=[24, 42, 106, 147, 188]    # TEMPERATURE

    if os.path.exists(os.path.join('../%s/ascii'%chid)):
        shutil.rmtree(os.path.join('../%s/ascii'%chid))

    for id_slice in id_slices:
        output="SOOT_"
    #==============================================================================

        os.makedirs(os.path.join('../%s/ascii/%s'%(chid, slices_dict[id_slice])))

        data=[]
        for a,b in enumerate(t_low):
            data.extend([path, chid, data_type, extend, domain_size,\
            t_low[a], t_up[a], number_slices, id_slice, 'ascii/%s/'%slices_dict[id_slice]+output+str(t_low[a])+'.txt'])

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
            print"!!! The Config File "+config_file_name+" does not exist or the path defined in the script is incorrect. !!!"
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
            # Loop through input fields (columns) for the line passed in from the fds2ascii_Config_File.csv file.
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

    print "*** Finished ***"
    # Script Written on 06.19.08 by bryan.klein@nist.gov
    # For usage instructions, see: http://code.google.com/p/fds-smv/wiki/fds2asciiBatchProcess
