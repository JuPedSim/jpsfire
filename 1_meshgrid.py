# -*- coding: utf-8 -*-
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
from mpltools import style
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 9}

matplotlib.rc('font', **font)

style.use('grayscale')

def cm2inch(value):
    return value/2.54

def dimensions(a):

    if a =='X':
        return 0,1
    elif a =='Y':
        return 2,3
    else:
        return 4,5
        
def read_fds_line():
    
    coordinates=[x.strip() for x in lines.split('XB')][1:]
    coordinates=''.join(coordinates)
    coordinates=re.split(', |,| |=|/',coordinates)
    coordinates=[ x for x in coordinates if '.' in x][0:6]
    coordinates=map(float, coordinates)
    #increasiing sorting of coordinate pairs:
    coordinates=sorted(coordinates[0:2])+sorted(coordinates[2:4])+sorted(coordinates[4:6])
    return coordinates

#===============================================================#
                                                                #
clip_dim3=2.8         #Dimension 3 zum Schnitt der Geometrie    #
                                                                #
#===============================================================#

simulations = glob.glob('../offenbach')
for simulation in simulations:
    chid = simulation[simulation.index('/')+1:]

    if not os.path.exists('../%s/converted/%s'%(chid, clip_dim3)):
        print 'create directory "converted"'
        os.makedirs('../%s/converted/%s'%(chid, clip_dim3))
    if not os.path.exists('../%s/slicefiles/%s'%(chid, clip_dim3)):
        print 'create directory "slicefiles"'
        os.makedirs('../%s/slicefiles/%s'%(chid, clip_dim3))

    
    fds=glob.glob('../%s/*.fds'%chid)
    #===== if necessary, adjust path and/or name of the converted asci slices======
    slicefiles=glob.glob('../%s/ascii/%s/*'%(chid, clip_dim3))
    
    if len(slicefiles)==0:
        sys.exit('CAUTION: directory comprises no converted slicefiles!')
    if len(fds)>1:
        sys.exit('CAUTION: directory comprises more than 1 FDS file')
    if len(fds)==0:
        sys.exit('CAUTION: directory comprises no FDS file!')
        
    with open(slicefiles[0]) as f:
        c = f.readlines()[0:1]    
        c=''.join(c)
    
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

    print 'projection level geometry', dimension_3, '=', clip_dim3
    
    dim1_1, dim1_2 = dimensions(dimension_1)                        
    dim2_1, dim2_2 = dimensions(dimension_2)
    dim3_1, dim3_2 = dimensions(dimension_3)
    
    obsts=[]
    holes=[]
    magnitudes=[]
    
    #for slicefile in slicefiles:
    for slicefile in slicefiles[6:7]:
            
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
        with open(fds[0]) as f:
            content = f.readlines()
        
            for lines in content:
                if lines[0:5]=='&OBST':
                    obst=read_fds_line()
                    if obst[dim3_1] < clip_dim3 < obst[dim3_2]:     
                        obst=[(1/delta_dim_1)*i for i in obst]
                        geometry[obst[dim2_1]-dim2_max/delta_dim_2-1 : obst[dim2_2]-dim2_max/delta_dim_2-1, \
                        obst[dim1_1]-dim1_max/delta_dim_1-1 : obst[dim1_2]-dim1_max/delta_dim_1-1] [:]=10
                        
                        if slicefile == slicefiles[0]:
                            obsts=np.append(obsts,obst)
            
            if slicefile == slicefiles[0]:
                obsts=np.reshape(obsts,(-1,6))
                np.savetxt('../%s/converted/obst.csv'%chid, obsts, delimiter=',')
    
            for lines in content:
                if lines[0:5]=='&HOLE':
                    hole=read_fds_line()
                    if hole[dim3_1] < clip_dim3 < hole[dim3_2]:
                       hole=[(1/delta_dim_1)*i for i in hole]
                       geometry[hole[dim2_1]-dim2_max/delta_dim_2-1 : hole[dim2_2]-dim2_max/delta_dim_2-1, \
                       hole[dim1_1]-dim1_max/delta_dim_1-1 : hole[dim1_2]-dim1_max/delta_dim_1-1] [:]=0
                       
                       if slicefile == slicefiles[0]:
                           holes=np.append(holes,hole)
            
            if slicefile == slicefiles[0]:
                holes=np.reshape(holes,(-1,6))
                np.savetxt('../%s/converted/hole.csv'%chid, holes, delimiter=',')
                
        collect = magnitudes + geometry
        #collect = collect[:-1,:-1] 
            
        np.savetxt('../%s/converted/%s/%s.csv'%(chid, clip_dim3, slicefile[slicefile.rfind('/')+1:-4]), collect, delimiter=',')
        print '\nWrite files: ../%s/converted/%s/%s.csv'%(chid, clip_dim3, slicefile[slicefile.rfind('/')+1:-4])
        
        #plt.imshow(collect)
#        plt.contourf(dim1, dim2, collect)   

        plt.figure(figsize=(cm2inch(15),cm2inch(10)), dpi=300)
        #plt.figure(figsize=(cm2inch(30),cm2inch(20)), dpi=300)
        ax=plt.subplot(111)
        aa = ax.pcolorfast(dim1, dim2, collect, vmin=0, vmax=1, cmap='gray_r')
        plt.colorbar(aa)
    
        plt.xlabel('%s (m)'%dimension_1.lower())
        plt.ylabel('%s (m)'%dimension_2.lower())

        plt.xticks([-5,0,5,10])

        plt.axes().set_aspect('equal')
        #plt.minorticks_on()
        #plt.grid(which='major')

        plt.savefig('../%s/slicefiles/%s/slice%s.pdf'%(chid, clip_dim3, slicefile[slicefile.rfind('/')+1:-4]))
        #plt.show()

    print "\n*** Finished ***"
    
#============storage of the most important variables to store.pckl=============

exits_dict ={"A":(-0.3,2.0,'y')}#, "C":(5.0,24.3,'x'), "B":(-0.3,22.0,'y')}#, "D":(-0.3,27.6,'y'), "E":(-4.3,-7.0,'y'), "F":(-2.3,34.3,'x')}
exits = collections.OrderedDict(sorted(exits_dict.items()))
   

f = open('store.pckl', 'w')
pickle.dump((chid, dimension_1, dimension_2, dim1, dim2, delta_dim_1, delta_dim_2,\
geometry, dim1_min, dim1_max, dim2_min, dim2_max, exits), f)

f.close() 