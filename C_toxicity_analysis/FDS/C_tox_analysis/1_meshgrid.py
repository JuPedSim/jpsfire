# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Created on Mon Jan 19 12:50:15 2015

@author: Benjamin
"""
import numpy as np
import glob
import matplotlib.pyplot as plt
import sys
import re
import pickle
import os
import collections
import matplotlib
from matplotlib import rcParams
import xml.etree.ElementTree as ET

rcParams.update({'figure.autolayout': True})

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 9}

matplotlib.rc('font', **font)


def cm2inch(value):
    return value/2.54

def dimensions(a):

    if a =='X':
        return 0,1
    elif a =='Y':
        return 2,3
    else:
        return 4,5

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def read_fds_line():

    coordinates=[x.strip() for x in lines.split('XB')][1:]
    coordinates=''.join(coordinates)
    coordinates=re.split(', |,| |=|/',coordinates)
    for j, cord in enumerate(coordinates):
        if isfloat(cord) == True:
            coordinates[j] = float(cord)
        if cord.isdigit() == True:
            coordinates[j] = float(cord)
        else:
            continue

    coordinates=[ x for x in coordinates if type(x)<>str]
    coordinates=map(float, coordinates)
    #increasiing sorting of coordinate pairs:
    coordinates=sorted(coordinates[0:2])+sorted(coordinates[2:4])+sorted(coordinates[4:6])
    return coordinates

f = open('data_slice2ascii.pckl')
(chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, plots) = pickle.load(f)

if not os.path.exists('../1_meshgrid/%s_%.2f'%(specified_location[0], specified_location[1])):
    print 'create directory "1_meshgrid"'
    os.makedirs('../1_meshgrid/%s_%.2f'%(specified_location[0], specified_location[1]))

try:
    os.remove('data_meshgrid.pckl')
except OSError:
    pass

if plots==True:
    if not os.path.exists('../slicefiles/%s_%.2f'%(specified_location[0], specified_location[1])):
        print 'create directory "slicefiles"'
        os.makedirs('../slicefiles/%s_%.2f'%(specified_location[0], specified_location[1]))


fds=glob.glob('../%s.fds'%chid)
#===== if necessary, adjust path and/or name of the converted asci slices======

slicefiles=glob.glob('../0_slice2ascii/%s_%.2f/*'%(specified_location[0], specified_location[1]))

if len(slicefiles)==0:
    sys.exit('CAUTION: directory comprises no converted slicefiles!')
if len(fds)>1:
    sys.exit('CAUTION: directory comprises more than 1 FDS file')
if len(fds)==0:
    sys.exit('CAUTION: directory comprises no FDS file!')

with open(slicefiles[0]) as sl:
    c = sl.readlines()[0:1]
    c=''.join(c)

sl.close()

dims=['X', 'Y', 'Z']

dimension_1=c[1]

print '\nInitialisation...'
print '\nFDS-File: ',fds[0]
print '\nDimension 1:', dimension_1
dimension_2=c[3]
print 'Dimension 2:', dimension_2

if dimension_1 in dims: dims.remove(dimension_1)
if dimension_2 in dims: dims.remove(dimension_2)

dimension_3=dims[0]
print 'Dimension 3:', dimension_3
print ''

print 'projection level geometry', dimension_3, '=', specified_location[1]

dim1_1, dim1_2 = dimensions(dimension_1)
dim2_1, dim2_2 = dimensions(dimension_2)
dim3_1, dim3_2 = dimensions(dimension_3)

obsts=[]
holes=[]
magnitudes=[]

for slicefile in slicefiles:
#for slicefile in slicefiles[6:7]:

    plt.close()
    sf=np.loadtxt(slicefile, skiprows=2, delimiter=',')
    print '\n-----------------------------------------------------------------'
    print 'Read files:', slicefile, '\n'

    dim1=sf[:,0]
    dim2=sf[:,1]
    magnitudes=sf[:,2]

    delta_dim_1=round(abs(dim1[1]-dim1[0]),3)
    print 'd%s'%dimension_1, '=', delta_dim_1
    for i,x in enumerate(dim2):
        if x<>dim2[i-1]:
            delta_dim_2=round(abs(dim2[i]-dim2[i-1]),3)
    print 'd%s'%dimension_2, '=', delta_dim_2

    dim1_min=round(np.amin(sf[:,0]),3)
    dim1_max=round(np.amax(sf[:,0]),3)

    dim2_min=round(np.amin(sf[:,1]),3)
    dim2_max=round(np.amax(sf[:,1]),3)

    print '\nDimension 1: %f (m) ... %f (m)' %(dim1_min, dim1_max)
    print 'Dimension 2: %f (m) ... %f (m)' %(dim2_min, dim2_max)

    dim1=np.arange(dim1_min, dim1_max+0.01, delta_dim_1)
    dim2=np.arange(dim2_min, dim2_max+0.01, delta_dim_2)

    XX, YY = np.meshgrid(dim1, dim2)
    magnitudes = np.reshape(magnitudes, (len(dim2),len(dim1)))

#==========extraction of &OBST + &HOLE lines ==================================
    geometry = np.zeros((len(dim2),len(dim1)))
    with open(fds[0]) as fd:
        content = fd.readlines()

        for lines in content:
            if lines[0:5]=='&OBST':
                obst=read_fds_line()
                if obst[dim3_1] < specified_location[1] < obst[dim3_2]:
                    obst=[(1/delta_dim_1)*i for i in obst]
                    geometry[obst[dim2_1]-dim2_max/delta_dim_2-1 : obst[dim2_2]-dim2_max/delta_dim_2-1, \
                    obst[dim1_1]-dim1_max/delta_dim_1-1 : obst[dim1_2]-dim1_max/delta_dim_1-1] [:]=10

                    if slicefile == slicefiles[0]:
                        obsts=np.append(obsts,obst)

        if slicefile == slicefiles[0]:
            obsts=np.reshape(obsts,(-1,6))
            np.savetxt('../1_meshgrid/obst.csv', obsts, delimiter=',')

        for lines in content:
            if lines[0:5]=='&HOLE':
                hole=read_fds_line()
                if hole[dim3_1] < specified_location[1] < hole[dim3_2]:
                   hole=[(1/delta_dim_1)*i for i in hole]
                   geometry[hole[dim2_1]-dim2_max/delta_dim_2-1 : hole[dim2_2]-dim2_max/delta_dim_2-1, \
                   hole[dim1_1]-dim1_max/delta_dim_1-1 : hole[dim1_2]-dim1_max/delta_dim_1-1] [:]=0

                   if slicefile == slicefiles[0]:
                       holes=np.append(holes,hole)

        if slicefile == slicefiles[0]:
            holes=np.reshape(holes,(-1,6))
            np.savetxt('../1_meshgrid/hole.csv', holes, delimiter=',')

        fd.close()

    collect = magnitudes + geometry

    #collect = collect[:-1,:-1]

    np.savetxt('../1_meshgrid/%s_%.2f/%s.csv'%(specified_location[0], specified_location[1], slicefile[slicefile.rfind('/')+1:-4]), collect, delimiter=',')
    print '\nWrite files: ../1_meshgrid/%s_%.2f/%s.csv'%(specified_location[0], specified_location[1], slicefile[slicefile.rfind('/')+1:-4])

    if plots==True:
        plt.figure(figsize=(cm2inch(15),cm2inch(10)), dpi=300)
        ax=plt.subplot(111)
        aa = ax.pcolorfast(dim1, dim2, collect, vmin=0, vmax=0.00001)
        plt.colorbar(aa, label=quantity)

        plt.xlabel('%s (m)'%dimension_1.lower())
        plt.ylabel('%s (m)'%dimension_2.lower())

        plt.axes().set_aspect('equal')
        plt.minorticks_on()
        plt.grid(which='major')

        plt.savefig('../slicefiles/%s_%.2f/%s.pdf'%(specified_location[0],specified_location[1], slicefile[slicefile.rfind('/')+1:-4]))


print "\n*** Finished ***"


#============storage of the most important variables to store.pckl=============

f = open('data_meshgrid.pckl', 'w')
pickle.dump((chid, quantity, specified_location, t_start, t_stop, t_step, id_meshes, plots, dimension_1, dimension_2, dim1, dim2, delta_dim_1, delta_dim_2, geometry, magnitudes, dim1_min, dim1_max, dim2_min, dim2_max), f)

f.close()
