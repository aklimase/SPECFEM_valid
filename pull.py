#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 16:25:43 2019

@author: aklimase
"""
import obspy
from obspy import read
from obspy import Stream
#from obspy.io.segy.core import _read_segy
from obspy.core.utcdatetime import UTCDateTime
import pandas as pd
import numpy as np
import glob
import os
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt
client = Client("IRIS")
import waveforms as wf
import os.path as path



##          Script to get synthetic seismograms from IRIS DMC                ##
##                        Gabriel Ferragut                                   ##
##                      November Dec 1st 2018                                ##



## Directories to be used in pulling data ##
workingdir = '/Users/aklimase/Documents/Cascadia_subduction/valid_seismo/'

# Need a list of earthquakes and stations
stndata = pd.read_csv(workingdir + 'PNSN_Network.txt', delim_whitespace = True, header = 0)

net = np.array(stndata["net"])
stnm = np.array(stndata["#stnm"])
chan = np.array(stndata["chan"])
loc = np.array(['*']* len(net))

# Start and end date/time for resp file:
starttime = '2015-01-01T00:00.000'
endtime = '2015-12-31T23:59.000'


## USGS Events from Northern Pacific Rim with P travel time gap (175 events)
#evdata_Ptt = pd.read_csv('/Users/aklimase/Desktop/test/N_PacificRim_Eq_Catalog_P-gap.txt', delim_whitespace = True)
## USGS Events from Northern Pacific Rim, no travel times (~ 4000 events) 
evdata = pd.read_csv(workingdir + 'cat.txt', delim_whitespace = True)

evnm = np.array(evdata["id"])
evlon = np.array(evdata["lon"])
evlat = np.array(evdata["lat"])
evDepth = np.array(evdata["depth"])
evMag = np.array(evdata["mag"])
tGap = np.array(evdata["gap"])
start_t = np.array(evdata["starttime"])

    
for i in range(len(net)):

    print("Working on station " + str(stnm[i]))

    time = UTCDateTime(start_t[0])
    bulk = [(net[i], stnm[i], '*', chan[i], time-40, time+120)]
    for ii in range(len(evnm)):
        st = client.get_waveforms_bulk(bulk)
        print st
        
        st.write(workingdir + 'uncorrected/' + str(net[i]) + "." + 
                 str(str(stnm[i]) + "." + chan[i] + "." + str(evnm[ii]) + ".sac"), format="sac")  
        
        print("Got a seismogram! Count: " + str(ii))
        
        
        
#make response files in response directory for each combo of networks and stations
response_path = workingdir + 'respfiles'

####################################################
#comment out if already have resp files
#makes a response file for each station and channel
for i in range(len(stnm)):
    respfile = response_path + '/' + net[i] + '.' + stnm[i] + '.' + chan[i] + '.resp'
    network, station = net[i], stnm[i]
    location= loc[i]
    channel = chan[i]
    print network, station
    wf.download_response(network,station,location,channel,starttime,endtime,respfile)


#print(st.sort()) 
f = glob.glob(workingdir + 'uncorrected/*.sac')
for i in range(len(f)):
    st = read(f[i])
    fig = st.plot()
    name = (f[i].split('/')[-1])[0:-4]#.split('.')[0:4]
    print name
    fig.savefig(workingdir + 'uncorrected/' + name + '.png')
        
#read in all uncorrected sac files and loop through and for each look in response directory, then remove response and save to directory of corrected files
icorr_path = workingdir + 'corrected'

#read in cut data, rmean and rtrend, find .resp file, correct and add to icorr dir
for i in range(0,len(f)):#in this case event paths are all sac files
    base = path.basename(f[i])
    network, station, channel, n, ext = base.split('.')
#    station = stnm[i]
#    channel = chan[i]
    #find response file
    respfile = response_path + '/' + network + '.' + station + '.' + channel + '.resp'
    #first rmean and rtrend
    stream = read(f[i])
    st = stream
#    st = stream
    tr = stream[0]
    ny = tr.stats.sampling_rate*0.5
    print ny
#    ny = 0.125
    ####################
    prefilt = (0.0,0.01,0.7*ny,ny)

    tr.detrend(type = 'demean')#removes mean of data
    tr.detrend(type = 'simple')#rtrend linear from first and last samples
    #rewrite to a sac file
    tr.write('temp.sac', format = 'sac')
    sacfile = 'temp.sac'
    icorr_sacfile = icorr_path +'/'+ base

    wf.remove_response(sacfile,respfile,icorr_sacfile, prefilt,'VEL')#, plot=True)#prefilt values for HH
    
    
f = glob.glob(workingdir + 'corrected/*.sac')
for i in range(len(f)):
    st = read(f[i])
    name = (f[i].split('/')[-1])[0:-4]#.split('.')[0:4]
    print name
    fig = st.plot()
    fig.savefig(workingdir + 'corrected/' + name + '.png')    
    
    
#f = glob.glob(workingdir + 'corrected/*.sac')
#for i in range(len(f)):
#    st = read(f[i])
#    name = (f[i].split('/')[-1])[0:-4]#.split('.')[0:4]
#    print name
##    st.filter("bandpass", freqmin=.001, freqmax=0.125,corners=2, zerophase=False)  
##    st.detrend(type = 'demean')#removes mean of data
##    st.detrend(type = 'simple')#rtrend linear from first and last samples
#    st.filter('lowpass', freq=.125, corners=2, zerophase=True)
#    fig = st.plot(type='relative')
#    fig.savefig(workingdir + 'lowpass/' + name + '.png')    
