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

workingdir= '/Users/aklimase/Documents/Cascadia_subduction/'
velocity = workingdir + 'casc1.6_velmdl.txt'
out = workingdir + 'tomo_elasticsmall.xyz'
spacing = 500

lon0, lon1, lat0,lat1 = -122, -129, 40.2, 50,  
myProj = Proj("+proj=utm +zone=10, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
x0,y0= myProj(lon0, lat0)
x1,y1 = myProj(lon1, lat1)
x0 = int(x0)
x1 = int(x1)
y0 = int(y0)
y1 = int(y1)
z0,z1 = -60000.0,0.0

#x0 = -10800.0
#x1 = 643200.0
#y0 = 4467300.0
#y1 = 5540300.0

x0 = 200200.0
x1 = 643200.0
y0 = 4600300.0
y1 = 5200300.0

nx = int(abs((x1-x0)/spacing + 1))
ny = int(abs((y1-y0)/spacing + 1))
nz = int(abs((z1-z0)/spacing + 1))

print nx, ny, nz
print 'number of lines:', nx*ny*nz

vpmin = 1107.0
vpmax = 8502.0
vsmin = 1.0
vsmax = 4834.0

vpminkm = float(vpmin)/1000.
vpmaxkm = float(vpmax)/1000.


rhomin = round(1000.*(1.6612*vpminkm - 0.4721*(vpminkm**2.) + 0.0671*(vpminkm**3.) - 0.0043*(vpminkm**4.) + 0.000106*(vpminkm**5.)),4)
rhomax = round(1000.*(1.6612*vpmaxkm - 0.4721*(vpmaxkm**2.) + 0.0671*(vpmaxkm**3.) - 0.0043*(vpmaxkm**4.) + 0.000106*(vpmaxkm**5.)),4)


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

#x0 = 0.0
#x1 = 0.0
#y0 = 4500000.0
#y1 = 5000000.0
#rho0 = 2
#rho1 = 2
vsmin = 3000
vsmax = 3000
vpmin = 3000
vpmax = 3000
while(1):
#    xlist = []
#    ylist = []
#    rholist = []
    vslist = []
#    vplist = []
    for lines in range(500000):
        l = input_file.readline()
        
        x,y,z,vp,vs = l.split()
        
        if float(x) >= float(x0) and float(x) <= float(x1) and float(y) >= float(y0) and float(y)<=float(y1):
    #        xlist.append(int(x))
    #        ylist.append(int(y))
            vslist.append(int(vs))
    #        vplist.append(int(vp))
            if float(vs) == 0.:
                vs = 1.0
    
            vp = float(vp)/1000.
            z = -1*float(z)
            rho  = round(1000.*(1.6612*vp - 0.4721*(vp**2.) + 0.0671*(vp**3.) - 0.0043*(vp**4.) + 0.000106*(vp**5.)),4)
    #        rholist.append(rho)
            outfile.write(str(float(x)*1.0) + '\t' + str(float(y)*1.0) +'\t' + str(float(z)*1.0) +'\t'+ str(vp*1000.) +'\t'+ str(float(vs)*1.0)+'\t' + str(rho) + '\n')
    count += 500000
    print 'number of lines read: ',count
#    if min(xlist) <= x0:
#        x0 = min(xlist)
#    if max(xlist) >= x1:
#        x1 = max(xlist)
#    if min(ylist) <= y0:
#        y0 = min(ylist)
#    if max(ylist) >= y1:
#        y1 = max(ylist)
#    if max(rholist) >= rho1:
#        rho1= max(rholist)
#    if min(rholist) <= rho0:
#        rho0= min(rholist)
#    if min(vslist) <= vsmin:
#        vsmin = min(vslist)
#    if max(vslist) >= vsmax:
#        vsmax = max(vslist)
#    if max(vplist) >= vpmax:
#        vpmax= max(vplist)
#    if min(vplist) <= vpmin:
#        vpmin =  min(vplist)
#        
        
        
#    user_input = raw_input('Type STOP to quit, otherwise press the Enter/Return key ')
#    if user_input == 'STOP':
#        break
#        out_file.close()
#print x0,x1,y0,y1, rho0, rho1
print vsmin
print vsmax
outfile.close()
    
    
    
    
    
    