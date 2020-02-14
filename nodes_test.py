#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 13:57:55 2019

@author: aklimase
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import numpy as np
import pyproj
import cPickle as pickle
from numpy import ones,r_,c_
from scipy.stats import binned_statistic_dd
import veldefs
from shapely.geometry.polygon import Polygon

stamen_terrain = cimgt.Stamen('terrain-background')

x0 = 450200.0
x1 = 600200.0
y0 = 4800300.0
y1 = 4950300.0


#plot validation stations and event in lat, lon (convert from x,y in utm)
velfile = '/Users/aklimase/Documents/Cascadia_subduction/Mesh/pckl_zslice/velocity_model_0.0.pckl'

mpl.rcParams['font.size'] = 18

workingdir = '/Users/aklimase/Documents/Cascadia_subduction/validation/'

stnsf = '/Users/aklimase/Documents/Cascadia_subduction/validation/PNSN_Network.txt'

#convert latitude and longitude to mesh coordinates
def latlon_to_mesh(lat, lon):
    p = pyproj.Proj(proj='utm', zone=10, ellps='WGS84')
    x, y = p(lon, lat)#lon, lat
    return x, y

#convert mesh coordinates (utm 10) to lat, lon to plot on the m43ap
def mesh_to_latlon(x,y):
    p = pyproj.Proj(proj='utm', zone=10, ellps='WGS84')
    lon, lat = p(x,y, inverse = True)
    return lat, lon

mesh_to_latlon(x0,y0)
mesh_to_latlon(x1,y1)

#read in station file
s = np.genfromtxt(stnsf, comments = '#', delimiter = '\t',names = True, dtype = None, encoding = None)
testlon = -123.086667
testlat = 44.051944

source_lat = 44.090 #N
source_lon = -122.831 #°W 


#all in utm 10
for i in range(len(s['lat'])):
    x,y = latlon_to_mesh(s['lat'][i], s['lon'][i])
    print s['stnm'][i],s['net'][i], round(y,4),round(x,4), 0.0, 0.0

x,y = latlon_to_mesh(s['lat'], s['lon'])
x_source, y_source = latlon_to_mesh(source_lat, source_lon)


fig = plt.figure(figsize = (10,10))#6,10
ax = plt.axes(projection=stamen_terrain.crs)
#ax.set_extent([-124.01, -121.99, 43.3, 44.8])
ax.set_extent([-123.615, -121.735, 43.354, 44.699])


ax.add_image(stamen_terrain, 10, alpha = 0.5)#, cmap = 'gray')
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='both', length = 5, width = 1)

xticks = [-130, -128, -126, -124, -123, -122.5, -122, -120, -118]
yticks = [40, 42, 43.5, 44, 44.5, 46, 48, 50]
ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels = True)

transform = ccrs.Geodetic()._as_mpl_transform(ax)


    
plt.scatter(s['lon'], s['lat'], s = 40, marker = '^',c = ['red'],transform=ccrs.Geodetic(),label = 'stations', zorder = 1000)
for i in range(len(s['stnm'])):
#    ax.text(s['lon'][i]+0.01, s['lat'][i]+0.01, str(s['stnm'][i]), fontsize=14, transform=ccrs.Geodetic())
    an = plt.annotate(str(s['stnm'][i]), xy = (s['lon'][i]+0.01, s['lat'][i]+0.01),xycoords=transform, xytext = (2, 2), textcoords = 'offset points')
    an.draggable()
    
plt.scatter(testlon, testlat, s = 30,transform=ccrs.Geodetic())
#ax.text(testlon-0.05, testlat-0.03, 'Eugene', fontsize=14, transform=ccrs.Geodetic())
an = plt.annotate('Eugene', xy = (testlon-0.05, testlat-0.03),xycoords=transform, xytext = (2, 2), textcoords = 'offset points')
an.draggable()

scale_bar(ax, 50)


plt.scatter(source_lon, source_lat, s = 50, marker = '*',c = 'black',transform=ccrs.Geodetic(),label = 'source', zorder = 1000)

plt.legend()
plt.savefig(workingdir + 'validiation_ev.png')
plt.show()



#%%

#plot events and stuff
def scale_bar(ax, length=None, location=(0.8, 0.05), linewidth=3):
    """
    ax is the axes to draw the scalebar on.
    length is the length of the scalebar in km.
    location is center of the scalebar in axis coordinates.
    (ie. 0.5 is the middle of the plot)
    linewidth is the thickness of the scalebar.
    """
    #Get the limits of the axis in lat long
    llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.Geodetic())
    #Make tmc horizontally centred on the middle of the map,
    #vertically at scale bar location
    sbllx = (llx1 + llx0) / 2
    sblly = lly0 + (lly1 - lly0) * location[1]
    tmc = ccrs.TransverseMercator(sbllx, sblly)
    #Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(tmc)
    #Turn the specified scalebar location into coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    #Calculate a scale bar length if none has been given
    #(Theres probably a more pythonic way of rounding the number but this works)
    if not length: 
        length = (x1 - x0) / 5000 #in km
        ndim = int(np.floor(np.log10(length))) #number of digits in number
        length = round(length, -ndim) #round to 1sf
        #Returns numbers starting with the list
        def scale_number(x):
            if str(x)[0] in ['1', '2', '5']: return int(x)        
            else: return scale_number(x - 10 ** ndim)
        length = scale_number(length) 

    #Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx - length * 500, sbx + length * 500]
    #Plot the scalebar
    ax.plot(bar_xs, [sby, sby], transform=tmc, color='k', linewidth=linewidth)
    #Plot the scalebar label
    ax.text(sbx, sby, str(length) + ' km', transform=tmc,
            horizontalalignment='center', verticalalignment='bottom')
    



#%%
#open velocity model file and plot the stations and event on top

f = open(velfile,'rb')
p = pickle.load(f)
f.close()
xlist = p.x[700000:1500000]
ylist = p.y[700000:1500000]
zlist = p.z[700000:1500000]
vslist = p.vs[700000:1500000]
vplist = p.vp[700000:1500000]

a = np.asarray(vslist)
vslat, vslon = mesh_to_latlon(xlist, ylist)

fig = plt.figure(figsize = (10,10))#6,10
ax = plt.axes(projection=stamen_terrain.crs)
ax.set_extent([-124.01, -121.99, 43.3, 44.8])
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='both', length = 5, width = 1)

xticks = [-130, -128, -126, -124, -122, -120, -118]
yticks = [40, 42, 44, 46, 48, 50]
ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels = True)

x_edges = np.arange(np.min(vslon),np.max(vslon), 0.01)
y_edges = np.arange(np.min(vslat),np.max(vslat), 0.01)
    
sample = c_[vslon, vslat]
med = np.median(vslist)
bindims = [x_edges, y_edges]
statistic_y,bin_edges,binnumber=binned_statistic_dd(sample, a, statistic= 'median', bins=bindims)

cmap = mpl.cm.get_cmap('viridis')
normalize = mpl.colors.Normalize(vmin=100, vmax=5000)
colors = [cmap(normalize(value)) for value in a]
s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
s_m.set_array([])   
   
X,Y = np.meshgrid(x_edges, y_edges)
plt.pcolormesh(X,Y, statistic_y.T, cmap = cmap, norm = normalize, transform=ccrs.PlateCarree())

plt.scatter(s['lon'], s['lat'], s = 20, marker = '^',c = ['red'],transform=ccrs.PlateCarree(),label = 'stations', zorder = 1000)
for i in range(len(s['stnm'])):
    ax.text(s['lon'][i]+0.01, s['lat'][i]+0.01, str(s['stnm'][i]), fontsize=14, transform=ccrs.PlateCarree())
    
plt.scatter(testlon, testlat, s = 30,transform=ccrs.PlateCarree())
ax.text(testlon+0.01, testlat+0.01, 'Eugene', fontsize=14, transform=ccrs.PlateCarree())

plt.scatter(source_lon, source_lat, s = 30, marker = '*',c = 'black',transform=ccrs.PlateCarree(),label = 'source', zorder = 1000)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.86, 0.12, 0.02, 0.75])
cbar = plt.colorbar(s_m, cax=cbar_ax)
cbar.set_label(ur"Vs (m/s)", fontsize = 20)#"$azimuth$ (\u00b0)"
cbar.ax.tick_params(labelsize = 18)


plt.legend()
plt.savefig(workingdir + 'validiation_ev_velmodel.png')



#%%

#cross sections
y1 = 4947300
y2 = 4507300
y3 = 5307300

lat1, lon1 = mesh_to_latlon(-126,y1)
lat2, lon2 = mesh_to_latlon(-126,y2)
lat3, lon3 = mesh_to_latlon(-126,y3)


vertices = [(-122, 40.2),(-129, 40.2),(-129,50),(-122,50),(-122, 40.2)]

#polygon = Polygon([[40.2,122],[40.2,129],[50,122],[50,129]], True,transform=ccrs.PlateCarree())

fig = plt.figure(figsize = (8,10))#6,10
ax = plt.axes(projection=stamen_terrain.crs)
#ax.set_extent([-124.01, -121.99, 43.3, 44.8])
ax.set_extent([-130.1, -117.9, 39.9, 50.1])

ax.add_image(stamen_terrain, 8, alpha = 0.5)#, cmap = 'gray')
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='both', length = 5, width = 1)

poly = Polygon(vertices)
x,y = poly.exterior.xy
ax.plot(x,y, color = 'red', transform=ccrs.Geodetic(), label = 'model extent', lw = 2)

ax.plot([-129,-122],[lat2,lat2],transform=ccrs.Geodetic(), lw = 2, label = 'y = 4507300, cross section b')
ax.plot([-129,-122],[lat1,lat1],transform=ccrs.Geodetic(), lw = 2, label = 'y = 4947300, cross section a')
ax.plot([-129,-122],[lat3,lat3],transform=ccrs.Geodetic(), lw = 2, label = 'y = 5307300, cross section c')

        
xticks = [-132, -130, -128, -126, -124, -122, -120, -118, -116]
yticks = [38, 40, 42, 44, 46, 48, 50, 52]
ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels = True)


scale_bar(ax, 100)

plt.legend()
plt.savefig(workingdir + 'ycrosslines.png')
plt.show()


