o#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 17:11:16 2020

@author: aklimase
"""

#pull seismograms, remove response, and filter for the given specfem run
from obspy.core.utcdatetime import UTCDateTime
import pandas as pd
import numpy as np
from obspy.clients.fdsn import Client
import matplotlib.pyplot as plt
client = Client("IRIS")



#change these variables for different model runs
examplename = '150km_tripling'#'500m_mesh_150km'# '500m_mesh_150km_dep4km'#
freqmax = 1.4
pre_filt = (0.08,  0.09, 0.7*freqmax, freqmax)
seistype = 'VEL'

workingdir = '/Users/aklimase/Documents/Cascadia_subduction/validation/runs/'

start_t = '2015-07-04T15:42:18.000Z'

## Directories to be used in pulling data ##
# Need a list of earthquakes and stations
stndata = pd.read_csv('/Users/aklimase/Documents/Cascadia_subduction/validation/PNSN_Network_noclip.txt', delim_whitespace = True, header = 0)
net = np.array(stndata["net"])
stnm = np.array(stndata["#stnm"])
chan = np.array(stndata["chan"])
loc = np.array(['*']* len(net))
lat = np.array(stndata["lat"])
lon = np.array(stndata["lon"])

t1 = UTCDateTime(start_t)
t2 = t1 + (45)#delta in seconds

fdsn_client = Client('IRIS')
# Fetch waveform from IRIS FDSN web service into a ObsPy stream object
# and automatically attach correct response
for i in range(0,len(net)):
    st = fdsn_client.get_waveforms(network=net[i], station=stnm[i], location='*',
                                   channel=chan[i], starttime=t1, endtime=t2,
                                   attach_response=True)
    
    tr = st[0]
    data_unfilt = tr.data
    tr.detrend(type = 'linear')
    tr.remove_response(output=seistype, pre_filt=pre_filt)
#    tr.remove_response(output='DISP', pre_filt=pre_filt)


    
    t = tr.times(reftime=UTCDateTime(start_t))
    data = tr.data

    icorr_sacfile = workingdir + examplename + '/'+ seistype + '/' + net[i] + '.' + stnm[i] + '.' + chan[i] + '.sac'
    tr.write(icorr_sacfile, format = 'sac')
    
    fig, axs = plt.subplots(1,2, sharex=True,figsize = (10,2))
    axs[0].plot(t, data_unfilt)
    axs[1].plot(t,data)
    axs[0].set_title('data instrument response corrected')
    axs[1].set_title('data filtered')
    plt.tight_layout()
    plt.xlim([-5,45])
    plt.savefig(workingdir + examplename + '/'+ seistype + '/data_' + stnm[i]  + '.' + chan[i] +'.png')
    plt.close()