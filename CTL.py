# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 20:53:32 2015

@author: nakamura
"""
from numpy import *
import scipy.linalg as linalg
import scipy.stats as stats
import scipy.signal as signal
import datamanager as datamgr
import matplotlib.pyplot as plt
import os



def centre(x):
    nt = len(x); mean = sum(x)/nt
    return x - ones(nt)*mean

def detrend1(x):# Linear trend
    nlist = range(len(x))
    nx = nlist*x; nn = nlist*nlist; v = nn.mean() - nlist.mean()**2
    b1 = (nx.mean() - nlist.mean() * x.mean()) / v
    b0 = x.mean() - b1 * nlist.mean() 
    return x - b0 + b1*nlist

def detrendMV(x, k): #k-moving average
    return 1

def ssa(x, r): # x is nt x 1 vector, r is embedding dimension
    #nt = len(x); mean = sum(x)/nt
    #x = x - ones(nt)*mean
    H = linalg.hankel(x, zeros(r)) #Construct Hankel matrix
    return linalg.svd(H, compute_uv=False)

def slope(x, y, alpha):# returns slope and margin of error
    x = array(x); y = array(y); n = len(x)
    xy = x * y; xx = x * x; v = xx.mean() - x.mean()**2
    b1 = (xy.mean() - x.mean() * y.mean()) / v
    b0 = y.mean() - b1 * x.mean()
    s2 = 1./n * sum([(y[i] - b0 - b1 * x[i])**2 for i in xrange(n)])
    c = -1 * stats.t.ppf(alpha/2.,n-2)
    db1 = c * sqrt( s2 / ((n-2) * v) )
    return [b1, db1, b0] # slope, uncertainty in slope, intercept
    
def slk(R, s, A, x0, dt, tmax):  
    #R:birth rate, s:standard dev., A:interaction matrix, func of time
    d = len(x0) 
    xt = array([x0]) 
    eye = ones(d)
    #Solving stochastic Lotka-Volterra equation 
    for n in range(0, int(tmax/dt)): 
        x = xt[n]
        xnew = (eye+random.normal(0,s,d)*dt)*x+dt*R*x*(eye-dot(A(n*dt),x))
        xt = vstack((xt, array([xnew])))
    return xt # xt is (tmax/dt) x d
    
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
    """
    Perform the full ssa analysis, with graphs saved to workDir
    N is window size
    lowPassFilter is expected to be a list/array of [order, cutoff],
                  where cutoff is between 0 to 1, as digital filter
    Set sampleNumCol to 1 if they first column of the data file is simply the data index
    Ensure working directory exists
    Note that data file must be contained in the working directory
    ex. if recordList[0] = '31' then data file is '31.txt' in working directory
    """
    if not os.path.exists(workDir[:len(workDir)-1]):
        os.makedirs(workDir[:len(workDir)-1])
    
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
            
            # Graphs
            plt.figure(figsize=(12,3)); 
            steps = 1
            if len(data) >= 100e3:
                steps = 10
            plt.plot(data[::steps, col], linewidth=0.5); 
            plt.xlabel('Sample Number', fontsize=20)
            plt.ylabel('$x$', fontsize=26)
            plt.savefig(workDir + record + str(col+1) + xstr(order) + xstr(cutoff) + '.png', dpi=150)
            plt.show()
            plt.clf() # Free up memory
            
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
            outstr = "{0:6.3f}+/-{1:6.3f} ({2:4.2f}), {3:6.3f}+/-{4:6.3f} ({5:4.2f}), {6:6.3f}+/-{7:6.3f} ({8:4.2f})"
            print(outstr.format(pc[0], pc[1], abs(pc[1]/pc[0]), ps[0], ps[1], abs(ps[1]/ps[0]), pv[0], pv[1], abs(pv[1]/pv[0]) ))
            
            print scores
            
            plt.figure(figsize=(22,6), dpi=150); 
            plt.subplot(131)
            
            plt.plot(alist,clist,"bo"); 
            plt.plot(alist,alist*pc[0]+pc[2],"g");
            plt.title("$\\nu$ vs $s$", fontsize=26)
            plt.xlabel("$s$", fontsize=26)
            plt.ylabel("$\\nu$", fontsize=26)
            #plt.ylim(1,4)
            #plt.savefig(workDir + record + str(col+1) + xstr(order) + xstr(cutoff) + 'nuvss.png', dpi=150)
            #plt.show()
            
            plt.subplot(132)
            #plt.figure(figsize=(6,6)); 
            plt.plot(alist,slist,"bo"); 
            plt.plot(alist,alist*ps[0]+ps[2],"g");
            plt.title("$\\lambda_1$ vs $s$", fontsize=26)
            plt.xlabel("$s$", fontsize=26)
            plt.ylabel("$\\lambda_1$", fontsize=26)
            #plt.ylim(1,3)
            #plt.savefig(workDir + record + str(col+1) + xstr(order) + xstr(cutoff) + 'lamvss.png', dpi=150)
            #plt.show()
            
            plt.subplot(133)
            #plt.figure(figsize=(6,6)); 
            plt.plot(alist,vlist,"bo"); 
            plt.plot(alist,alist*pv[0]+pv[2],"g");
            plt.title("$\\sigma^2$ vs $s$", fontsize=26)
            plt.xlabel("$s$", fontsize=26)
            plt.ylabel("$\\sigma^2$", fontsize=26)
            #plt.ylim(1,3)
            #plt.savefig(workDir + record + str(col+1) + xstr(order) + xstr(cutoff) + 'sigvss.png', dpi=300)
            plt.savefig(workDir + record + str(col+1) + xstr(order) + xstr(cutoff) + 'indicators500.png', dpi=150)
            plt.show()
            
            """
            plt.figure(figsize=(5,4)); 
            plt.plot(alist,aclist,"bo"); 
            plt.plot(alist,alist*pac[0]+pac[2],"g");
            plt.title("Autocorrelation vs $s$", fontsize=26)
            plt.xlabel("$s$", fontsize=26)
            plt.ylabel("Autocorrelation", fontsize=26)
            #plt.ylim(1,3)
            #plt.savefig(workDir + record + str(col+1) + order + cutoff + 'sigvss.png', dpi=150)
            plt.show()
            """
            plt.clf() # Free up memory    
        datafp.close()