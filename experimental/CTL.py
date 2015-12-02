# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 20:53:32 2015

@author: Travis
"""
from numpy import *
import scipy.linalg as scilinalg
import scipy.stats as stats
import scipy.signal as signal
import matplotlib.pyplot as plt
import pandas as pd


"""
SSA Functions + Helpers
"""
def findMean(x):
    return sum(x) / len(x), len(x)

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

def ssa(x, r, fullMatr=True, comp_uv=False): # x is nt x 1 vector, r is embedding dimension
    #nt = len(x); mean = sum(x)/nt
    #x = x - ones(nt)*mean
    H = scilinalg.hankel(x, zeros(r)) #Construct Hankel matrix
    return scilinalg.svd(H, full_matrices=fullMatr, compute_uv=comp_uv)

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

"""
Analysis Functions
"""
def nanStartEnd(dfArray):
    """
    Find the start and end of data in array with nans
    """
    start = 0
    end = len(dfArray)
    for idx in range(1, len(dfArray)):
        if isnan(dfArray[idx-1]) and not isnan(dfArray[idx]) and start == 0:
            start = idx
        
        if isnan(dfArray[idx]) and not isnan(dfArray[idx-1]):
            end = idx # not x-1 as array slice will stop at x
    return start, end
    
def xstr(s):
    """
    Return empty string if s doesn't exist
    """
    if s is None:
        return ''
    return str(s)

def analysis(df, title, workDir='', winSize=50, lowPassFilter=None):
    """
    Perform analysis, with graphs saved to workDir

    lowPassFilter is expected to be a list/array of [order, cutoff],
        where cutoff is between 0 to 1, as digital filter
    """
    for col in df:
        start, end = nanStartEnd(df.ix[:, col])
        data = df.ix[start:end, col].copy()
            
        # Low Pass Filter
        order, cutoff = None, None
        if lowPassFilter is not None:
            order = lowPassFilter[0]
            cutoff = lowPassFilter[1]
            b, a = signal.butter(order, cutoff, btype='lowpass', analog=False, output='ba')
            data[:] = signal.filtfilt(b, a, data[:])
            
            # For saving images with cutoff + order info
            order = '_' + str(order) + 'order'
            cutoff = '_' + str(cutoff) + 'order'
        
        # Graphs
        graphTitle = workDir + title + '_' + str(col) + 'col_' + str(winSize) + 'winSize' + xstr(order) + xstr(cutoff)
        fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(15,15))
        data[:].plot(ax=axes[0]); 
        axes[0].set_ylabel('$x$', fontsize=26)
        axes[0].set_title("Data", fontsize=26)
        
        clist=[]; slist=[]; vlist=[]; aclist = []
        scores = [0,0,0]
        N = int(len(data) / winSize) 
        for k in range(2*N-1):
            r = N; 
            S = ssa(centre( data[int(k*winSize/2):int((k+2)*winSize/2)] ), r) #int(P/2)
            L = (S/S[0])**2
            idc = array(range(1,r+1)); m=sum(L[:r])
            clist.append( dot(L[:r],idc)/m )
            slist.append(S[0]**2/winSize)
            variance = var(data[int(k*winSize/2):int((k+2)*winSize/2)])
            vlist.append(variance)
            aclist.append(dot(data[int(k*winSize/2):int((k+2)*winSize/2)-10],data[int(k*winSize/2)+10:int((k+2)*winSize/2)])/variance)
        
        alist = linspace(0, len(data), len(clist))
        
        # Generate date index for plotting
        xlist = []
        for x in alist:
            x = int(x)
            while x >= len(data.index):
                x = x - 1 # if x is over by one
            xlist.append(data.index[x])
        
        pc = slope(alist,clist,0.05)
        ps = slope(alist,slist,0.05)
        pv = slope(alist,vlist,0.05)
        pac = slope(alist,aclist,0.05)
        if pc[0]<0 and abs(pc[1]/pc[0])<1: scores[0]=scores[0]+1
        if ps[0]>0 and abs(ps[1]/ps[0])<1: scores[1]=scores[1]+1
        if pv[0]>0 and abs(pv[1]/pv[0])<1: scores[2]=scores[2]+1
            
        pcSlope = "{0:6.3f}+/-{1:6.3f} ({2:4.2f})".format(pc[0], pc[1], abs(pc[1]/pc[0]))
        psSlope = "{0:6.3f}+/-{1:6.3f} ({2:4.2f})".format(ps[0], ps[1], abs(ps[1]/ps[0]))
        pvSlope = "{0:6.3f}+/-{1:6.3f} ({2:4.2f})".format(pv[0], pv[1], abs(pv[1]/pv[0]))
        
        print(pcSlope + ', ' + psSlope + ', ' + pvSlope)
        
        print scores
        
        axes[1].plot(xlist,clist,"bo");
        axes[1].plot(xlist,alist*pc[0]+pc[2],"g",)
        axes[1].set_title("$\\nu$ vs $s$: " + pcSlope, fontsize=26)
        axes[1].set_ylabel("$nu$", fontsize=26)
        axes[1].grid(True)

        axes[2].plot(xlist,slist,"bo");
        axes[2].plot(xlist,alist*ps[0]+ps[2],"g")
        axes[2].set_title("$\\lambda_1$ vs $s$: " + psSlope, fontsize=26)
        axes[2].set_ylabel("$\\lambda_1$", fontsize=26)
        axes[2].grid(True)
        
        axes[3].plot(xlist,vlist,"bo")
        axes[3].plot(xlist,alist*pv[0]+pv[2],"g")
        axes[3].set_title("$\\sigma^2$ vs $s$: " + pvSlope, fontsize=26)
        axes[3].set_ylabel("$\\sigma^2$", fontsize=26)
        axes[3].grid(True)
        
        fig.autofmt_xdate()
        plt.show()
        fig.savefig(graphTitle + '_indicators500.png', dpi=300)
        
        """
        plt.figure(figsize=(5,4)); 
        plt.plot(alist,aclist,"bo"); 
        plt.plot(alist,alist*pac[0]+pac[2],"g");
        plt.title("Autocorrelation vs $s$", fontsize=26)
        plt.xlabel("$s$", fontsize=26)
        plt.ylabel("Autocorrelation", fontsize=26)
        #plt.ylim(1,3)
        #plt.savefig(graphTitle + '_sigvss.png', dpi=150)
        plt.show()
        """
        plt.clf() # Free up memory    

"""
Interpolation Functions
"""
def ssaReconstruction(df, start=0, end=None, numEigen=6, winSize=50):
    """
    Reconstruct a signal with SSA using the specified number of eignevectors
    """
    # To prevent overwriting the caller's array
    df = df.copy()
    dfArray = df.values
    
    if end == None:
        end = len(dfArray)
    shortdfArray = dfArray[start:end]
    
    # Calculate the eigenvectors
    N = int(len(shortdfArray) / winSize)
    for k in range(2*N-1):
        mean, length = findMean(shortdfArray[int(k*winSize/2):int((k+2)*winSize/2)])
        U, S, V = ssa(centre(shortdfArray[int(k*winSize/2):int((k+2)*winSize/2)]), N, False, True)
        reducedS = zeros(len(S))
        for idx in range(numEigen):
            reducedS[idx] = S[idx]
        S = diag(S)
        reconstr = dot(U, dot(S, V))
        dfArray[int(k*winSize/2)+start:int((k+2)*winSize/2)+start] = reconstr[:, 0] + ones(length) * mean
    return df
    
def ssaInterpolation(df, numEigen=6, winSize=50, minError=0.01):
    """
    Use SSA to interpolate a signal, however don't touch invalid values before
    the first valid data point, and the last valid data point
    """
    # To prevent overwriting the caller's array
    df = df.copy()
    
    for col in df:
        dfArray = df.ix[:, col].values
        start, end = nanStartEnd(dfArray)
        
        error = minError + 1.0 # ensure loop is hit first time
        while error > minError:
            oldArray = df.ix[:, col].values.copy()
            reconstArray = ssaReconstruction(df.ix[:, col], start, end, numEigen, winSize)
            dfArray[start:end] = reconstArray[start:end]
            error = ((linalg.norm(oldArray[start:end]) - linalg.norm(dfArray[start:end])) / 
                    linalg.norm(oldArray[start:end]))

            
    return df

def nanAvgFill(df):
    """
    df is assumed to be a signal which has invalid (nan) points within it
    
    For a given invalid point, the average value from valid points before, and after,
    will be used to linearly interpolate the invalid point
    """
    # To prevent overwriting the caller's array
    df = df.copy()    
    
    for col in df:
        # Mark where all invalid entries start and end
        # Except for all NaN values before the first valid data point
        dfArray = df.ix[:, col].values # for speed
        dfIndex = df.index.values.astype('datetime64[s]')
        
        nanidxs = findNans(dfArray)
        
        # Interpolate NaN data points linearly with averages from valid points
        for idx in range(1, len(nanidxs)-1):
            currNanStart = nanidxs[idx][0]     
            currNanEnd = nanidxs[idx][1]
            prevNanEnd = nanidxs[idx-1][1]
            nextNanStart = nanidxs[idx+1][0]
    
            avgStart = average(dfArray[prevNanEnd:currNanStart])
            avgEnd = average(dfArray[currNanEnd:nextNanStart])
            
            timeStart = dfIndex[currNanStart-1]
            timeEnd = dfIndex[currNanEnd]
            deltaTime = pd.to_timedelta(timeEnd - timeStart) / timedelta64(1,'s')
            
            slope = (avgEnd - avgStart) / deltaTime # per second
            
            deltaTime = pd.to_timedelta(dfIndex[currNanStart] - dfIndex[currNanStart-1]) / timedelta64(1,'s')
            dfArray[currNanStart] = avgStart + deltaTime * slope 
            for idy in range(1, len(dfArray[currNanStart:currNanEnd])):
                deltaTime = pd.to_timedelta(dfIndex[currNanStart+idy] - dfIndex[currNanStart+idy-1]) / timedelta64(1,'s')
                df.ix[currNanStart+idy, col] = dfArray[currNanStart+idy-1] + deltaTime * slope 
    return df

def nanHybridFill(df):
    """
    df is assumed to be a signal which has invalid (nan) points within it
    
    First, for a given invalid point, this function will determine if it is a long
    string of invalid points (< 3 in a row), and will then use pandas to do a simple
    linear interpolation to fill the points
    
    For a long string of invalid points (> 3), the average value before and after the
    string will be calculated and be used for invalid points. This method exists to
    explore the results of several different methods for interpolation
    """
    # To prevent overwriting the caller's array
    df = df.copy()    
    
    for col in df:
        # Mark where all invalid entries start and end
        # Except for all NaN values before the first valid data point
        dfArray = df.ix[:, col].values # for speed and overwrite protection
        dfIndex = df.index.values.astype('datetime64[s]')
        
        nanidxs = findNans(dfArray)
        
        # Linearly interpolate small groups of nans inbetween data
        newidxs = []
        interp = df[col].interpolate()
        for idx in range(len(nanidxs)):
            start = nanidxs[idx][0]
            end = nanidxs[idx][1]
            
            if start == 0 or end == len(dfArray):
                newidxs.append(nanidxs[idx])
                continue # skip beginning and ending NaNs
            
            if (end - start) < 3:
                df.ix[start:end, col] = interp[start:end]
            else:
                newidxs.append(nanidxs[idx])
        
        nanidxs = newidxs
        
        # Interpolate NaN data points linearly with averages from valid points
        for idx in range(1, len(nanidxs)-1):
            currNanStart = nanidxs[idx][0]     
            currNanEnd = nanidxs[idx][1]
            prevNanEnd = nanidxs[idx-1][1]
            nextNanStart = nanidxs[idx+1][0]
    
            avgStart = average(dfArray[prevNanEnd:currNanStart])
            avgEnd = average(dfArray[currNanEnd:nextNanStart])
            
            timeStart = dfIndex[currNanStart-1]
            timeEnd = dfIndex[currNanEnd]
            deltaTime = pd.to_timedelta(timeEnd - timeStart) / timedelta64(1,'s')
            
            slope = (avgEnd - avgStart) / deltaTime # per second
            
            deltaTime = pd.to_timedelta(dfIndex[currNanStart] - dfIndex[currNanStart-1]) / timedelta64(1,'s')
            dfArray[currNanStart] = avgStart + deltaTime * slope 
            for idy in range(1, len(dfArray[currNanStart:currNanEnd])):
                deltaTime = pd.to_timedelta(dfIndex[currNanStart+idy] - dfIndex[currNanStart+idy-1]) / timedelta64(1,'s')
                dfArray[currNanStart+idy] = dfArray[currNanStart+idy-1] + deltaTime * slope 
    return df

def findNans(npArray):
    """
    This will return the indexes of all the nan values in an array
    """
    nanidxs = []
    start = None
    for x in range(len(npArray)):
        # If current entry is NaN, and previous entry is valid
        if isnan(npArray[x]) and not isnan(npArray[x-1]):
            start = x
        
        # If current entry is valid, and previous entry is NaN
        # Then add to the start and end index to validEntries
        if isnan(npArray[x-1]) and not isnan(npArray[x]):
            if start == None: # Special case; first point is NaN
                nanidxs.append([0, x])
            else:
                nanidxs.append([start, x])
                
    # Special case to catch last set of NaNs
    if start not in nanidxs:
        nanidxs.append([start, len(npArray)])
    
    return nanidxs