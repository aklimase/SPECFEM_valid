#python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 17:16:20 2019

@author: aklimase
"""
import numpy as np
import cPickle as pickle
#import pickle
import veldefs
from veldefs import velmod
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import glob
import pandas as pd
from numpy import ones,r_,c_
from scipy.stats import binned_statistic_dd

workingdir= '/Users/aklimase/Documents/Cascadia_subduction/Mesh/'

pfiles = glob.glob(workingdir + "pckl_zslice/velocity_model_*")
pfiles = sorted(pfiles)


#read in first file to create a dataframe to concatenate with
f = open(pfiles[0], 'rb')
p0 = pickle.load(f)
d = {'x': p0.x, 'y':p0.y, 'z':p0.z,'vs':p0.vs,'vp':p0.vp,'rho':p0.rho}
df = pd.DataFrame(data=d)
f.close()

#loop takes a while to run
for i in range(1,len(pfiles)):#leave off first
	f = open(pfiles[i], 'rb')
	p_temp = pickle.load(f)
	print i
	f.close()
	d_temp = {'x': p_temp.x, 'y':p_temp.y, 'z':p_temp.z,'vs':p_temp.vs,'vp':p_temp.vp,'rho':p_temp.rho}

	df_temp = pd.DataFrame(data=d_temp)
    
	df = pd.concat([df, df_temp], ignore_index=True)

p = veldefs.velmod(df.x,df.y,df.z,df.vp,df.vs)

##plot the x cross sections, every 60 files
##grid and create a color mesh
##change vp to vs is needed
#xcross = sorted(list(set(p.x)))
#for i in range(0,len(xcross),60):  
#    print i
#    x_here = xcross[i]
#
#    ind = np.where(p.x == x_here)
#
#    xlist = np.asarray([p.x[i] for i in ind])[0]
#    ylist = np.asarray([p.y[i] for i in ind])[0]
#    zlist = np.asarray([p.z[i] for i in ind])[0]
#    vplist = [p.vp[i] for i in ind]
#    a = np.asarray(vplist)
#    a = a[0]
#    
#    
#    y_edges = np.arange(np.min(ylist),np.max(ylist), 1000)
#    z_edges = np.arange(np.min(zlist),np.max(zlist), 500)
#    
#    sample = c_[ylist, zlist]
#    med = np.median(vplist)
#    bindims = [y_edges, z_edges]
#    statistic_y,bin_edges,binnumber=binned_statistic_dd(sample, a, statistic= 'median', bins=bindims)
#         
#    cmap = mpl.cm.get_cmap('viridis')
#    normalize = mpl.colors.Normalize(vmin=1000, vmax=max(a))
#    colors = [cmap(normalize(value)) for value in a]
#    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
#    s_m.set_array([])       
#
#    fig = plt.figure(figsize = (12,5))
#    plt.subplot(111)
#    plt.tick_params(axis='both', which='major', labelsize=14)
#    plt.tick_params(axis='both', which='both', length = 3, width = 1)
#    plt.gca().invert_yaxis()
#
#    X,Y = np.meshgrid(y_edges, z_edges)
#    plt.pcolormesh(X,Y, statistic_y.T, cmap = cmap, norm = normalize)
#    
#    plt.xlabel('y coordinate (m)', fontsize = 20)
#    plt.ylabel('z coordinate (m)', fontsize = 20)
#    plt.title('x =' + str(x_here), fontsize = 22)
#
#    fig.subplots_adjust(right=0.84, bottom = 0.14)
#    cbar_ax = fig.add_axes([0.86, 0.12, 0.02, 0.75])
#    cbar = plt.colorbar(s_m, cax=cbar_ax)
#    cbar.set_label(ur"Vp (m/s)", fontsize = 20)#"$azimuth$ (\u00b0)"
#    cbar.ax.tick_params(labelsize = 18)
#    plt.savefig(workingdir + 'yz_cross_Vp/x_' + str(x_here) + '.png')
#    plt.close()
#    
#


#plot y cross sections every 80 files
ycross = sorted(list(set(p.y)))
for i in range(0,len(ycross), 80):  
    y_here = ycross[i]

    ind = np.where(p.y == y_here)
    
    xlist = np.asarray([p.x[i] for i in ind])[0]
    ylist = np.asarray([p.y[i] for i in ind])[0]
    zlist = np.asarray([p.z[i] for i in ind])[0]
    vplist = [p.vp[i] for i in ind]
    a = np.asarray(vplist)
    a = a[0]
    
    x_edges = np.arange(np.min(xlist),np.max(xlist), 1000)
    z_edges = np.arange(np.min(zlist),np.max(zlist), 500)
    
    sample = c_[xlist, zlist]
    med = np.median(vplist)
    bindims = [x_edges, z_edges]
    statistic_y,bin_edges,binnumber=binned_statistic_dd(sample, a, statistic= 'median', bins=bindims)

    cmap = mpl.cm.get_cmap('viridis')
    normalize = mpl.colors.Normalize(vmin=1000, vmax=max(a))
    colors = [cmap(normalize(value)) for value in a]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])       

    fig = plt.figure(figsize = (18,4))
    plt.subplot(111)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tick_params(axis='both', which='both', length = 3, width = 1)
    plt.gca().invert_yaxis()

    X,Y = np.meshgrid(x_edges, z_edges)
    plt.pcolormesh(X,Y, statistic_y.T, cmap = cmap, norm = normalize)

    plt.xlabel('x coordinate (m)', fontsize = 20)
    plt.ylabel('z coordinate (m)', fontsize = 20)
    plt.title('y =' + str(y_here), fontsize = 22)

    fig.subplots_adjust(right=0.84, bottom = 0.14)
    cbar_ax = fig.add_axes([0.86, 0.12, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"Vp (m/s)", fontsize = 20)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)
    plt.savefig(workingdir + 'xz_cross_Vp_long/y_' + str(y_here) + '.png')
    plt.close()
    




flist = glob.glob(workingdir + "pckl_zslice/velocity_*.pckl")

for i in range(0,len(flist),5):  
    
    fname =flist[i]
    f = open(fname,'rb')
    p = pickle.load(f)
    f.close()
    
#    p = pickle.load(open(flist[i],"r"))
    d = flist[i].split('_')[-1].split('.')[0]
    d = float(d)
    print d

    xlist = np.asarray(p.x)
    ylist = np.asarray(p.y)
    zlist = np.asarray(p.z)
    vslist = p.vs
    a = np.asarray(vslist)
    
    x_edges = np.arange(np.min(xlist),np.max(xlist), 2000)
    y_edges = np.arange(np.min(ylist),np.max(ylist), 2000)
    
    sample = c_[xlist, ylist]
    med = np.median(vslist)
    bindims = [x_edges, y_edges]
    statistic_y,bin_edges,binnumber=binned_statistic_dd(sample, a, statistic= 'median', bins=bindims)

    cmap = mpl.cm.get_cmap('viridis')
    normalize = mpl.colors.Normalize(vmin=100, vmax=5000)
    colors = [cmap(normalize(value)) for value in a]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])       

    fig = plt.figure(figsize = (10,12))
    plt.subplot(111)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tick_params(axis='both', which='both', length = 3, width = 1)
    
    X,Y = np.meshgrid(x_edges, y_edges)
    plt.pcolormesh(X,Y, statistic_y.T, cmap = cmap, norm = normalize)

    plt.xlabel('x coordinate (m)', fontsize = 20)
    plt.ylabel('y coordinate (m)', fontsize = 20)
    plt.title('z =' + str(d/1000.) + ' km', fontsize = 22)
    plt.tight_layout()

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"Vs (m/s)", fontsize = 16)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 12)
    plt.savefig(workingdir + 'depth_cross_Vs/z_' + str(d/1000.) + '_km.png')
    plt.close()
    print 'done with : ', d
