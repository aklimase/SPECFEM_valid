n#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 10:13:53 2019

@author: aklimase
"""

import obspy
from obspy import read
from obspy.core.utcdatetime import UTCDateTime
import pandas as pd
import numpy as np
import glob
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
client = Client("IRIS")
import os.path as path
import dread
from obspy.signal.filter import bandpass
import seaborn as sns
import scipy.signal as signal
from mtspec import mtspec

#Bad stations = bro, fris b030
plt.style.use("classic")

sns.set_context("poster")
sns.set_style('whitegrid')

#mpl.rcParams['font.size'] = 20
sns.set(font_scale=1.5)

#%%
#change these variables for different model runs


examplename = '150km_tripling'#'500m_mesh_150km'# '500m_mesh_150km_dep4km'#
freqmax = 1.4
pre_filt = (0.08,  0.09, 0.7*freqmax, freqmax)
samplingrate = 0.004

#plotchan = 'E'
#plotchanspecfem = 'X'
#plotchan = 'N'
#plotchanspecfem = 'Y'
plotchan = 'Z'
plotchanspecfem = 'Z'
seistype = 'VEL'
spectype = 'semv'

#examplename = 'simple_tomo_2_Talapas'
#freqmax = 0.25
#samplingrate = 0.005
#pre_filt = (0.08,  0.09, 0.7*freqmax, freqmax)

output_dir = '/Users/aklimase/software/specfem3d/my_examples/' + examplename
samplingrate = 0.004
workingdir = '/Users/aklimase/Documents/Cascadia_subduction/validation/runs/'

#%%

start_t = '2015-07-04T15:42:18.000Z'

### Directories to be used in pulling data ##
## Need a list of earthquakes and stations
stndata = pd.read_csv('/Users/aklimase/Documents/Cascadia_subduction/validation/PNSN_Network_noclip.txt', delim_whitespace = True, header = 0)
net = np.array(stndata["net"])
stnm = np.array(stndata["#stnm"])
chan = np.array(stndata["chan"])
loc = np.array(['*']* len(net))
lat = np.array(stndata["lat"])
lon = np.array(stndata["lon"])
#%% 
       
#plot filtered data
test = glob.glob(output_dir + '/OUTPUT_FILES/*CX' + plotchanspecfem + '.' + spectype)


 #plot each specfem station on its on plot
for i in range(len(test)):
    name = path.basename(test[i])[:-4]
    print name

    trace2 = np.genfromtxt(test[i])
    trspecfem = trace2.T[1]
    tspecfem = trace2.T[0]
    stats = obspy.core.trace.Stats()
    stats.sampling_rate = 1.0/samplingrate
    stats.npts = len(trspecfem)
    trt = obspy.core.trace.Trace(data = trspecfem, header = stats)
    trt = bandpass(trt, 0.1, 0.7*freqmax, trt.stats.sampling_rate, corners=2, zerophase=True)

    fig, axs = plt.subplots(1,2, sharex=True,figsize = (10,2))

    axs[1].plot(tspecfem,trt)
    axs[0].plot(tspecfem,trspecfem)
    axs[1].set_title('specfem filtered')
    axs[0].set_title('specfem unfiltered')
    plt.tight_layout()

    plt.xlim([-5,45])
    plt.savefig(workingdir + examplename + '/'+ seistype + '/specfem_out_' + name  + 'png')
    plt.close()

#plot all filtered seismograms


#%%
#order stations by distance to the source
source_lat = 44.090 #N
source_lon = -122.831 #°W 
source_dep = 4.0 #km

x = []
for i in range(len(test)):
    name_specfem = (test[i].split('/')[-1])[0:-4]
    nets, stas, chans, ns =  name_specfem.split('.')
    print name_specfem 
    for j in range(len(stnm)):
        if stnm[j] == stas and chan[j][2] == plotchan and chans[2] == plotchanspecfem:
            print stnm[j], stas, chans, ns
            dist = dread.compute_rrup(source_lon, source_lat, source_dep, lon[j], lat[j], 0.0)
            x.append([stnm[j], test[i], dist])

station_dist = sorted(x, key=lambda x:x[2])

n = len(station_dist)

##%%
##make a n x 2 figure of each station with data on left and specfem on right
#n = len(station_dist)
#fig, (axs) = plt.subplots(n,2, sharex=True,figsize = (16,10))
#
#for i in range(len(station_dist)):
#    print station_dist[i][0]
#    specfem =  station_dist[i][1]
#    r = round(station_dist[i][2],2)
#    trace2 = np.genfromtxt(specfem)
#    trspecfem = trace2.T[1]
#    tspecfem = trace2.T[0]
#    stats = obspy.core.trace.Stats()
#    stats.sampling_rate = 1.0/samplingrate
#    stats.npts = len(trspecfem)
#    trt = obspy.core.trace.Trace(data = trspecfem, header = stats)
#    trt = bandpass(trt, 0.1, 0.7*freqmax, trt.stats.sampling_rate, corners=2, zerophase=True)
#    
#    stn = station_dist[i][0]
#    
#    #########################
#    data = glob.glob(workingdir + examplename + '/' + seistype + '/' + '*' + stn + '*' + plotchan + '.sac')[0]
#    
#    st = read(data)
#    tr = st[0]
#    t = tr.times(reftime=UTCDateTime(start_t))
#    
#    axs[i][0].plot(t,tr.data)
#    axs[i][1].plot(tspecfem,trt)
#    
#    axs[i][0].set_xlim([-5,45])
#    axs[i][1].set_xlim([-5,45])
#
#    
#    limdata = max(abs(tr.data))
#    limspecfem = max(abs(trt))
#    lim = max(limdata, limspecfem)
#
#    
#    axs[i][0].set_ylim([-1*lim, lim])
#    axs[i][1].set_ylim([-1*lim, lim])
#
#    axs[i][0].locator_params(tight=True, nbins=4, axes = 'y')
#    axs[i][1].locator_params(tight=True, nbins=4, axes = 'y')
#
#        
#    axs[i][0].annotate(stn + ' data ' + str(r) + ' km', xy = (0.5,0.8), xycoords = 'axes fraction')
#    axs[i][1].annotate(stn + ' specfem ' +  str(r) + ' km', xy = (0.5,0.8), xycoords = 'axes fraction')
#
#    fig.text(0.5, 0.04, 'time after event (s)', ha='center', fontsize = 20)
#    fig.text(0.04, 0.5, 'm/s', va='center', rotation='vertical', fontsize = 20)
#        
#    plt.savefig(workingdir + examplename + '/compare_all_corrected_sharey_' + plotchan + 'comp.png')
    
    
#%%
 #plot each  station on same axes
fig, (axs) = plt.subplots(n,figsize = (11,10), sharex = True)
for i in range(len(station_dist)):
    print station_dist[i][0]
    specfem =  station_dist[i][1]
    
    r = round(station_dist[i][2],2)

    trace2 = np.genfromtxt(specfem)
    trspecfem = trace2.T[1]
    tspecfem = trace2.T[0]
    stats = obspy.core.trace.Stats()
    stats.sampling_rate = 1.0/samplingrate
    stats.npts = len(trspecfem)
    trt = obspy.core.trace.Trace(data = trspecfem, header = stats)
    trt = bandpass(trt, 0.1, 0.7*freqmax, trt.stats.sampling_rate, corners=2, zerophase=True)
    
    stn = station_dist[i][0]
    
    #########################
    data = glob.glob(workingdir + examplename + '/' + seistype + '/' + '*' + stn + '*' + plotchan + '.sac')[0]
    
    st = read(data)
    tr = st[0]
    t = tr.times(reftime=UTCDateTime(start_t))
    
    axs[i].plot(t,tr.data, c = 'green', label = 'data')
    axs[i].plot(tspecfem,trt,  c = 'blue',label  = 'specfem')
    
    axs[i].set_xlim([-5,45])

    limdata = max(abs(tr.data))
    limspecfem = max(abs(trt))
    lim = max(limdata, limspecfem)
    
    axs[i].set_ylim([-1*lim, lim])
    axs[i].set_ylim([-1*lim, lim])

    axs[i].locator_params(tight=True, nbins=4, axes = 'y')
#    axs[i].annotate(stn + ' data ' + str(r) + ' km', xy = (0.5,0.8), xycoords = 'axes fraction')
    axs[i].set_title(stn + ' data ' + str(r) + ' km', fontsize = 16)
    
    fig.text(0.45, 0.02, 'time after event (s)', ha='center', fontsize = 20)
    fig.text(0.0, 0.5, 'm/s', va='center', rotation='vertical', fontsize = 20)
#    axs[0].legend(loc = 'upper right')
#    plt.legend(bbox_to_anchor=(0, 1), loc='upper left', ncol=1)
    axs[0].legend(bbox_to_anchor=(1.0, 1), loc='upper left')
    
#    plt.tight_layout()
    
plt.subplots_adjust(right = 0.8) 
plt.savefig(workingdir + examplename + '/compare_sameax_sharey_' +seistype + '_' + plotchan + 'comp.png')

    
plt.close('all')


#%%
##make spectrograms
#
##st.spectrogram(log=True, title='BW.RJOB ' + str(st[0].stats.starttime))
#
#fig, (axs) = plt.subplots(n,2, sharex=True,figsize = (16,10))
#for i in range(len(station_dist)):
#    print station_dist[i][0]
#    specfem =  station_dist[i][1]
#    
#    r = round(station_dist[i][2],2)
#
#    trace2 = np.genfromtxt(specfem)
#    trspecfem = trace2.T[1]
#    tspecfem = trace2.T[0]
#    stats = obspy.core.trace.Stats()
#    stats.sampling_rate = 1.0/samplingrate
#    stats.npts = len(trspecfem)
#    trt = obspy.core.trace.Trace(data = trspecfem, header = stats)
#    samplerate = trt.stats.sampling_rate
#    trt = bandpass(trt, 0.1, 0.7*freqmax, trt.stats.sampling_rate, corners=2, zerophase=True)
#    
#    stn = station_dist[i][0]
#    
#    #########################
#    data = glob.glob(workingdir + examplename + '/' + seistype + '/' + '*' + stn + '*' + plotchan + '.sac')[0]
##    
#    st = read(data)
#    tr = st[0]
#    t = tr.times(reftime=UTCDateTime(start_t))
#    
#    
#    tr.spectrogram(log=False, title='test' + str(tr.stats.starttime), axes = axs[i][0])
#    axs[i][0].set_ylim([0,1])
#
#    trt = obspy.core.trace.Trace(data = trspecfem, header = stats)
#    trt.spectrogram(log=False, title='test' + str(trt.stats.starttime), axes = axs[i][1])
#    axs[i][1].set_ylim([0,1])
#    
#    axs[i][0].locator_params(tight=True, nbins=4, axes = 'y')
#    axs[i][1].locator_params(tight=True, nbins=4, axes = 'y')
#    
#    axs[i][0].annotate(stn + ' data ' + str(r) + ' km', xy = (0.5,0.8), xycoords = 'axes fraction', color = 'white')
#    axs[i][1].annotate(stn + ' specfem ' +  str(r) + ' km', xy = (0.5,0.8), xycoords = 'axes fraction',color = 'white')
#    
#    fig.text(0.5, 0.04, 'time after event (s)', ha='center', fontsize = 20)
#    fig.text(0.04, 0.5, 'frequency', va='center', rotation='vertical', fontsize = 20)
#
#    plt.savefig(workingdir + examplename + '/spectrogram_' + plotchan + 'comp.png')
#


#%%
 #plot spectra on same axes
#fig, (axs) = plt.subplots(n,figsize = (10,10), sharex = True)
#for i in range(len(station_dist)):
#    print station_dist[i][0]
#    r = round(station_dist[i][2],2)
#    stn = station_dist[i][0]
#    
#    #specfem
#    specfem =  station_dist[i][1]
#    trace2 = np.genfromtxt(specfem)
#    trspecfem = trace2.T[1]
#    tspecfem = trace2.T[0]
#    stats = obspy.core.trace.Stats()
#    stats.sampling_rate = 1.0/samplingrate
#    stats.npts = len(trspecfem)
#    trt = obspy.core.trace.Trace(data = trspecfem, header = stats)
##    print trt.data
#    filt = bandpass(trt.data, 0.1, 0.7*freqmax, trt.stats.sampling_rate, corners=4, zerophase=True)
##    print trt
#    
#
##    filt = obspy.signal.filter.bandpass(data=trt.data, freqmin=0.1, freqmax= 0.5, df = 1.0/samplingrate, corners=2, zerophase=False)
#    
#    spec1, freq1 , jack, fstat, dof =  mtspec(data = filt, delta = samplingrate, time_bandwidth = 4, number_of_tapers=7, quadratic = True, statistics = True)
#    
#    
##    spec1, freq1, jackknife, _, _ = mtspec(samplingrate, delta=samplingrate, time_bandwidth=3.5,number_of_tapers=5, nfft=312, statistics=True)
#
#    
#    #########################
#    #data
#    data = glob.glob(workingdir + examplename + '/' + seistype + '/' + '*' + stn + '*' + plotchan + '.sac')[0]
#    
#    st = read(data)
#    tr = st[0]
#    t = tr.times(reftime=UTCDateTime(start_t))
#    # Calculate the spectral estimation.
#    spec2, freq2, jackknife, _, _ = mtspec(data=tr.data, delta=samplingrate, time_bandwidth = 4, number_of_tapers=7, quadratic = True, statistics = True)
#    
#    #plot both
#    axs[i].plot(freq1, spec1, c = 'blue', label = 'data')
#    axs[i].plot(freq2, spec2,  c = 'orange', label  = 'specfem')
#    
#    axs[i].set_xlim([0,2])
##
##    limdata = max(abs(tr.data))
##    limspecfem = max(abs(trt))
##    lim = max(limdata, limspecfem)
##    
##    axs[i].set_ylim([-1*lim, lim])
##    axs[i].set_ylim([-1*lim, lim])
##
##    axs[i].locator_params(tight=True, nbins=4, axes = 'y')
#    
#    fig.text(0.5, 0.04, 'frequency', ha='center', fontsize = 20)
#    fig.text(0.04, 0.5, 'm/s', va='center', rotation='vertical', fontsize = 20)
#    axs[0].legend(loc = 'upper right')
#        
#    plt.savefig(workingdir + examplename + '/PSD_' + plotchan + 'comp.png')
#
#    
#plt.close('all')
