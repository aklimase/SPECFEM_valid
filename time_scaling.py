#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 12:17:02 2020

@author: aklimase
"""
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

#elements, proc, solver time s, freq

s_500m = [10800000, 32, 101*60, 1]
s_500m_triple = [4000000, 32, 28*60, 1.4]

full_2km =  [5267970, 32, 50*60, 0.25]

full1_5km = [337980, 1, 1279, 0.06]
full32_5km = [337980, 32, 152, 0.06]
full2_5km = [337980, 2,1342, 0.06]
full8_5km = [337980, 8, 463, 0.06]



#
##elements per proc
#x = [(full2_5km[0]), (s_500m_triple[0]), (full_2km[0]), (s_500m[0])]
##solver time
#y = [(full2_5km[2]), (s_500m_triple[2]), (full_2km[2]), (s_500m[2])]
#fig = plt.figure()
#ax = plt.gca()
#ax.scatter((full1_5km[0]),(full1_5km[2]), c = 'b', label = '5km 1 proc')
#ax.scatter((full2_5km[0]),(full2_5km[2]), c = 'g', label = '5km 32 proc')
#ax.scatter((full_2km[0]),(full_2km[2]), c = 'violet', label = '2km 32 proc')
#ax.scatter((s_500m[0]),(s_500m[2]), c = 'r', label = '500m 32 proc')
#ax.scatter((s_500m_triple[0]),(s_500m_triple[2]), c = 'orange', label = '500m triple layer 32 proc')
#ax.plot(x, y, zorder = 0, ls = '--', label = '32 proc')
##ax.set_yscale('log')
##ax.set_xscale('log')
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0d'))
#plt.xlabel('elements')
#plt.ylabel('solver time (s)')
#ax.set_ylim([1,7000])
#plt.legend(loc = 'upper left')
#plt.title('solver time vs. spectral elements')
#plt.show()
#
#fig = plt.figure()
#ax = plt.gca()
#ax.scatter((full1_5km[0]),(full1_5km[2]), c = 'b', label = '5km 1 proc')
#ax.scatter((full32_5km[0]),(full32_5km[2]), c = 'g', label = '5km 32 proc')
#ax.scatter((full_2km[0]),(full_2km[2]), c = 'violet', label = '2km 32 proc')
#ax.scatter((s_500m[0]),(s_500m[2]), c = 'r', label = '500m 32 proc')
#ax.scatter((s_500m_triple[0]),(s_500m_triple[2]), c = 'orange', label = '500m triple layer 32 proc')
#ax.plot(x, y, zorder = 0, ls = '--', label = '32 proc')
#ax.set_yscale('log')
#ax.set_xscale('log')
#plt.xlabel('elements')
#plt.ylabel('solver time (s)')
#ax.set_ylim([1,7000])
#plt.legend(loc = 'upper left')
#plt.title('solver time vs. spectral elements')
#plt.show()
#
#
#

#elements and solver time
fig = plt.figure()
ax = plt.gca()
#ax.scatter((full1_5km[0])/(full1_5km[1]),(full1_5km[2]), c = 'b', label = '5km 1 proc')
ax.scatter((full32_5km[0])/(full32_5km[1]),(full32_5km[2]), c = 'g', label = '5km 32 proc')
ax.scatter((full_2km[0])/(full_2km[1]),(full_2km[2]), c = 'violet', label = '2km 32 proc')
ax.scatter((s_500m[0])/(s_500m[1]),(s_500m[2]), c = 'r', label = '500m 32 proc')
ax.scatter((s_500m_triple[0])/(s_500m_triple[1]),(s_500m_triple[2]), c = 'orange', label = '500m triple layer 32 proc')

ax.scatter((full2_5km[0])/(full2_5km[1]),(full2_5km[2]), c = 'g', marker = '*', label = '5km 2 proc')
ax.scatter((full2_5km[0])/(full8_5km[1]),(full8_5km[2]), c = 'g', marker = '^', label = '5km 8 proc')

ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('elements per process')
plt.ylabel('solver time (s)')
#ax.set_ylim([10,120])
plt.legend(loc = 'upper left')
plt.title('solver time vs. spectral elements per process')
plt.show()



#proc and solver time
fig = plt.figure()
ax = plt.gca()
ax.scatter((full1_5km[1]),(full1_5km[2]), c = 'b', label = '5km 1 proc')
ax.scatter((full2_5km[1]),(full2_5km[2]), c = 'g', label = '5km 2 proc')
ax.scatter((full8_5km[1]),(full8_5km[2]), c = 'orange', label = '5km 8 proc')
ax.scatter((full32_5km[1]),(full32_5km[2]), c = 'r', label = '5km 32 proc')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('number of processes')
plt.ylabel('solver time (s)')
#ax.set_ylim([1,120])
plt.legend(loc = 'lower left')
plt.title('solver time vs. number of proc')
plt.show()

