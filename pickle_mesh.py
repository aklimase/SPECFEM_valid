#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 12:06:09 2019

@author: aklimase
"""

import numpy as np
import cPickle as pickle
import veldefs
import glob
import pandas as pd

workingdir= '/Users/aklimase/Documents/Cascadia_subduction/Mesh/'
#velocity = workingdir + 'casc_1000.txt'
velocity = workingdir + 'casc1.6_velmdl.txt'
#data = np.genfromtxt(velocity,usecols = 0)

        
        
def make_pickle(test,i):
    
    print len(test)
    print len(test[0])
    
    x_list = []
    y_list = []
    z_list = []
    vp_list = []
    vs_list = []
    
    for j in range(len(test)):
#        for k in range(len(test[j])):
        x_list.append(float(test[j][0]))
        y_list.append(float(test[j][1]))
        z_list.append(float(test[j][2]))
        vp_list.append(float(test[j][3]))
        vs_list.append(float(test[j][4]))

    database = veldefs.velmod(x_list,y_list,z_list,vp_list,vs_list)
    pickle_out = open(workingdir + "/pckl_zslice/velocity_model_" + str(i).zfill(2) + ".pckl","wb")
    pickle.dump(database, pickle_out)
    #pickle.dump(database, pickle_out_all)
    pickle_out.close()





#f = open(velocity)
#namecount = 0
#
#c = 0
#t = []
#z_cut = 0.
#with open(velocity) as f:
#    for i in xrange(300000000):
#        f.next()
#    for line in f:
#        ###################################################
#        l = line.split()
#        z = float(l[2])
##        print l
##        print z_cut
##        print len(t)
#        if z == z_cut:
#            print l
#            #if c%10 == 0:#every 10th line
#            t.append(l)
##        else:
#    make_pickle(t, z_cut)
##            t = []
##            c = 0
##            namecount +=1
##            z_cut = z
#
#f.close()



#pfiles = [(workingdir + "pckl/velocity_model_" + str(i) + ".pckl") for i in range(30)]
pfiles = glob.glob(workingdir + "pckl_zslice/velocity_model_*")
print pfiles
pfiles = sorted(pfiles)

p0 = pickle.load(open(pfiles[0], 'rb'))
d = {'x': p0.x, 'y':p0.y, 'z':p0.z,'vs':p0.vs,'vp':p0.vp,'rho':p0.rho}
df = pd.DataFrame(data=d)


for i in range(1,len(pfiles)):#leave off first
    p_temp = pickle.load(open(pfiles[i], 'rb'))
    print i

    d_temp = {'x': p_temp.x, 'y':p_temp.y, 'z':p_temp.z,'vs':p_temp.vs,'vp':p_temp.vp,'rho':p_temp.rho}

    df_temp = pd.DataFrame(data=d_temp)
    
    df = pd.concat([df, df_temp], ignore_index=True)

database = veldefs.velmod(df.x,df.y,df.z,df.vp,df.vs)
pickle_out_all = open(workingdir + "pckl_zslice/velocity_model_all.pckl","wb")
pickle.dump(database, pickle_out_all)
pickle_out_all.close()
