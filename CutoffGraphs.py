# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 20:52:12 2015

@author: Travis
"""

from numpy import *
from CTL import *
import datamanager as datamgr
import matplotlib.pyplot as plt
import os
import errno

###
# Functions
###
# Helper function for fullAnalysis to allow str(None) to return ''
def xstr(s):
    if s is None:
        return ''
    return str(s)

# Perform the full analysis, with graphs saved to workDir
# N is window size
# lowPassFilter is expected to be a list/array of [order, cutoff],
#               where cutoff is between 0 to 1, as digital filter
# Set sampleNumCol to 1 if they first column of the data file is simply the data index
def fullAnalysis(recordList, workDir='', lowPassFilter=None, sampleNumCol=0, embedDim=20, N=50):
    indSlope = []   
    # Ensure working directory exists
    # Note that data file must be contained in the working directory
    # ex. if recordList[0] = '31' then data file is '31.txt' in working directory
    try:
        os.makedirs(workDir[:len(workDir)-1])
    except OSError:
        if not os.path.isdir(workDir[:len(workDir)-1]):
            raise
    
    for record in recordList:
        datamgr.processData(record, workDir, sampleNumCol) # Ensure hdf5 file for record exists
        datafp = datamgr.accessData(record, workDir)
        data = datafp[record] # numpy array to be used    
        rows, cols = data.shape    
        
        for col in xrange(0, cols):
            # Low Pass Filter
            order, cutoff = None, None
            if lowPassFilter is not None:
                order = lowPassFilter[0]
                cutoff = lowPassFilter[1]
                b, a = signal.butter(order, cutoff, btype='lowpass', analog=False, output='ba')
                data[:, col] = signal.filtfilt(b, a, data[:, col])
                
                # For saving images with cutoff + order info
                order = 'Order' + str(order)
                cutoff = 'Cutoff' + str(cutoff)

            clist=[]; slist=[]; vlist=[]; aclist = []
            scores = [0,0,0]
            P = len(data) / N # window size
            r = embedDim # embedding dimension
            for k in range(2*N-1):
                r = embedDim; 
                S = ssa(centre( data[int(k*P/2):int((k+2)*P/2), col] ), r) #int(P/2)
                L = (S/S[0])**2
                idc = array(range(1,r+1)); m=sum(L[:r])
                clist.append( dot(L[:r],idc)/m )
                slist.append(S[0]**2/P)
                variance = var(data[int(k*P/2):int((k+2)*P/2), col])
                vlist.append(variance)
                aclist.append(dot(data[int(k*P/2):int((k+2)*P/2)-10, col],data[int(k*P/2)+10:int((k+2)*P/2), col])/variance)
                
            # alist will be big, therefore create a dataset inside
            # the temporary datafp to store alist
            alist = datafp.require_dataset('time', (len(data), len(data)), dtype=int32)
            alist = linspace(0, len(data), len(clist))
            
            pc = slope(alist,clist,0.05)
            ps = slope(alist,slist,0.05)
            pv = slope(alist,vlist,0.05)
            pac = slope(alist,aclist,0.05)
            if pc[0]<0 and abs(pc[1]/pc[0])<1: scores[0]=scores[0]+1
            if ps[0]>0 and abs(ps[1]/ps[0])<1: scores[1]=scores[1]+1
            if pv[0]>0 and abs(pv[1]/pv[0])<1: scores[2]=scores[2]+1
            plt.clf() # Free up memory  
            indSlope.append([alist, alist*pc[0]+pc[2], alist*ps[0]+ps[2], alist*pv[0]+pv[2]])
        datafp.close()
    return indSlope

###
# Runtime
###
indSlopes = []
recordList = ["AirTemp"]
order = 10
lowPassFilter = [None, [order, 0.05], [order, 0.1], [order, 0.25], [order, 0.5], [order, 0.75]]
for lpf in lowPassFilter:
    indSlopes.append(fullAnalysis(recordList, 'cutoffgraphs/', lpf))

lowPassFilter[0] = [0, "No Filter"]
label = []
for ix in range(len(indSlopes)):
    label.append(r"cutoff = %s" % lowPassFilter[ix][1])

workDir = 'cutoffgraphs/'
colors = ['r', 'g', 'b', 'y', 'c', 'k']
plt.figure(figsize=(22,6), dpi=150); 
plt.subplot(131)

#plt.plot(alist,clist,"bo"); 
for ix, indSlope in enumerate(indSlopes):
    plt.plot(indSlope[0][0], indSlope[0][1], colors[ix]);
plt.title("$\\nu$ vs $s$", fontsize=26)
plt.xlabel("$s$", fontsize=26)
plt.ylabel("$\\nu$", fontsize=26)
#plt.ylim(1,4)
#plt.savefig(workDir + record + str(col+1) + xstr(order) + xstr(cutoff) + 'nuvss.png', dpi=150)
#plt.show()

plt.subplot(132)
#plt.figure(figsize=(6,6)); 
#plt.plot(alist,slist,"bo"); 
for ix, indSlope in enumerate(indSlopes):
    plt.plot(indSlope[0][0], indSlope[0][2], colors[ix])
plt.title("$\\lambda_1$ vs $s$", fontsize=26)
plt.xlabel("$s$", fontsize=26)
plt.ylabel("$\\lambda_1$", fontsize=26)
#plt.ylim(1,3)
#plt.savefig(workDir + record + str(col+1) + xstr(order) + xstr(cutoff) + 'lamvss.png', dpi=150)
#plt.show()

plt.subplot(133)
#plt.figure(figsize=(6,6)); 
#plt.plot(alist,vlist,"bo"); 
for ix, indSlope in enumerate(indSlopes):
    plt.plot(indSlope[0][0], indSlope[0][3], colors[ix]);
plt.title("$\\sigma^2$ vs $s$", fontsize=26)
plt.xlabel("$s$", fontsize=26)
plt.ylabel("$\\sigma^2$", fontsize=26)
plt.legend(label, loc='upper right')
#plt.ylim(1,3)
#plt.savefig(workDir + record + str(col+1) + xstr(order) + xstr(cutoff) + 'sigvss.png', dpi=300)
plt.savefig(workDir + 'indicatorslopes', dpi=300)
plt.show()