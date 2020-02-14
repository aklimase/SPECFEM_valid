#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 17:16:20 2019

@author: aklimase
"""
import numpy as np
import cPickle as pickle
#import veldefs
import pandas as pd
from pyproj import Proj

workingdir= '/Users/aklimase/Documents/Cascadia_subduction/MESH/'
#my input file
velocity = workingdir + 'casc1.6_velmdl.txt'
#output tomo file for specfem
out = workingdir + 'tomo_elasticEug_150.xyz'
#spacing of model in meters
spacing = 500

#boundaries of model in utm (or whatever coordinates match your mesh)
#x0 = 400200.0

#eug no water
x0 = 450200.0
x1 = 600200.0
y0 = 4800300.0
y1 = 4950300.0
z0 = -60000.0
z1 = 0.0

#calculate number of points for the tomo file header
nx = int(abs((x1-x0)/spacing + 1))
ny = int(abs((y1-y0)/spacing + 1))
nz = int(abs((z1-z0)/spacing + 1))

print nx, ny, nz

#velocity boundaries for the header
#vpmin = 1488.0
#vpmax = 8483.0
#vsmin = 1.0
#vsmax = 4825.0

#Eug 150
vpmin = 1158.0
vpmax = 8454.0
vsmin = 607.0
vsmax = 4809.0

vpminkm = vpmin/1000.
vpmaxkm = vpmax/1000.
vsminkm = vsmin/1000.
vsmaxkm = vsmax/1000.

#use Tom Brocher's rho formula
rhomin = round(1000.*(1.6612*vpminkm - 0.4721*(vpminkm**2.) + 0.0671*(vpminkm**3.) - 0.0043*(vpminkm**4.) + 0.000106*(vpminkm**5.)),4)
rhomax = round(1000.*(1.6612*vpmaxkm - 0.4721*(vpmaxkm**2.) + 0.0671*(vpmaxkm**3.) - 0.0043*(vpmaxkm**4.) + 0.000106*(vpmaxkm**5.)),4)

#write the tomo file header in specfem format
input_file = open(velocity,'r')
outfile = open(out, 'w')
header1 = str(x0) + '\t' + str(y0) + '\t'+ str(z0)+ '\t'+ str(x1) + '\t'+ str(y1) + '\t'+ str(z1) + '\n'
header2 = str(spacing) + '\t' +str(spacing) + '\t'+ str(spacing) + '\n'
header3 = str(nx) + '\t' +str(ny) + '\t'+ str(nz)+ '\n'
header4 = str(vpmin) + '\t' +str(vpmax) + '\t'+ str(vsmin) + '\t' + str(vsmax) + '\t' +str(rhomin) + '\t'+ str(rhomax)+ '\n'
outfile.write(header1)
outfile.write(header2)
outfile.write(header3)
outfile.write(header4)

count = 0

vsmin = 3000
vsmax = 3000
vpmin = 3000
vpmax = 3000
while(1):
    vslist = []
    for lines in range(500000):
        l = input_file.readline() 
        x,y,z,vp,vs = l.split()
        #if the 
        if float(x) >= float(x0) and float(x) <= x1 and float(y) >= float(y0) and float(y) <= float(y1):
        	#set water velocity from 0 to min value
            if float(vs) < float(vsmin):
                vsmin = vs
            elif float(vs) > float(vsmax):
                 vsmax = vs
            if float(vp) < float(vpmin):
                vpmin = vp
            elif float(vp) > float(vpmax):
                vpmax = vp
            
            if float(vs) == 0.:
                vs = 1.0
#                print 'water'
    
            vp = float(vp)/1000.
            #make all of the depths negative
            z = -1*float(z)
            rho  = round(1000.*(1.6612*vp - 0.4721*(vp**2.) + 0.0671*(vp**3.) - 0.0043*(vp**4.) + 0.000106*(vp**5.)),4)
    
            outfile.write(str(float(x)*1.0) + '\t' + str(float(y)*1.0) +'\t' + str(float(z)*1.0) +'\t'+ str(vp*1000.) +'\t'+ str(float(vs)*1.0)+'\t' + str(rho) + '\n')
    count += 500000
    print 'number of lines read: ',count

#usually i have to run this after all the lines are read in and an error is thrown
outfile.close()
    
    
print vsmin, vsmax
print vpmin, vpmax
    
