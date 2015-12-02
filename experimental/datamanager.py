# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 17:54:36 2015

@author: Travis
"""

import os
from sys import exit
from shutil import copy2
import h5py
import numpy as np
import pandas as pd

###
# File Vars
###
hdf5_ext = ".hdf5"

###
# Functions
###

def floatOrInt(fp):     
    """ 
    Read first line of a file to figure out if stored data
    is a float or integer
    """
    isInt = 1
    
    pos = fp.tell()
    line = fp.readline()
    if line.find('.'):
        isInt = 0
    fp.seek(pos)
    
    if isInt:
        return np.dtype('int32')
    else:
        return np.dtype('float64')

def readBlock(fp, size=65536):
    """
    Helper function to read a file in blocks
    """
    while True:
        block = fp.read(size)
        if not block: break
        yield block

   
def countRowCol(filename):
    """
    Fast method to count rows and cols in large file
    Assumed white space between each column
    Return as (rows, cols)     
    """
    # Rows
    with open(filename, 'r') as fp:
        rows = sum(blk.count("\n") for blk in readBlock(fp)) + 1
    
    # Cols
    with open(filename, 'r') as fp:
        line = fp.readline()
        cols = len((line.strip()).split())
    
    return rows, cols

def processData(record, workDir='', timeCol=0):
    """
    Read in a text data file, and store into an hdf5 file for later access
    
    Set timeCol = 1 if the first column of the text file are timestamps in the form of:
    yyyy-mm-ddThh:mm:ss
    """
    path_txt = workDir + record + ".txt"
    path_hdf5 = workDir + record + hdf5_ext
    tDataType = "S" + str(len('yyyy-mm-ddThh:mm:ss'))
    
    if not os.path.isfile(path_hdf5):
        # Ensure record.txt exists before running code
        try: 
            fp = open(path_txt, 'r')
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
            exit(1)
        
        # Count number of lines and columns in data file for allocation
        rows, cols = countRowCol(path_txt)
        
        # Allocate hdf5 dataset
        dataType = floatOrInt(fp)
        datafp = h5py.File(path_hdf5, 'a')
        data = datafp.require_dataset(record, (rows, cols - timeCol), dtype=dataType)
        if timeCol:
            time = datafp.require_dataset(record + 'time', (rows, 1), dtype=tDataType)
        
         # Read data file into numpy array as integer
        BUF_SIZE = int(1e6) # ~1 MB
        readBuf = fp.readlines(BUF_SIZE)
        readpos = np.uint32(0) # Position in read file
        datapos = np.uint32(0) # Position in dataset
        
        # Create the write buffer
        # Ensure that BUF_SIZE / 4 is int
        try:    
            writeBuf = np.zeros(shape=(len(readBuf), cols - timeCol), dtype=dataType)
            timeWriteBuf = np.zeros(shape=(len(readBuf), 1), dtype=tDataType)
        except TypeError:
            print "BUF_SIZE/4 is not an int"
            fp.close()
            datafp.close()
            exit(1)
            
        while readBuf != []:
            # Put each line into the write buffer
            for line in readBuf:
                line = line.split()
                if timeCol:
                    timeWriteBuf[readpos-datapos] = line.pop(0)
                line = np.array(line, dtype=dataType)
                for idcol, sample in enumerate(line):
                    try:
                        writeBuf[readpos-datapos][idcol] = sample
                    except IndexError:
                        print "Index error on writeBuf"
                        fp.close()
                        datafp.close()
                        exit(1)
                readpos += 1
            
            # Write the buffer into the hdf5
            try:
                data[datapos:readpos] = writeBuf[:len(readBuf)]  
                if timeCol:
                    time[datapos:readpos] = timeWriteBuf[:len(readBuf)]
            except IndexError:
                print "Index error on dataset and writeBuf"
                fp.close()
                datafp.close()
                exit(1)
            datapos = readpos
            readBuf = fp.readlines(BUF_SIZE)
                    
        fp.close()
        datafp.close()
        
def accessData(record, workDir=''):
    """
    Access the stored hdf5 data file, and pass it back as a
    temporary file to allow it to be written without affecting
    the original file
    Returns the record numpy array from the temp file
    """
    path_hdf5 = workDir + record + hdf5_ext
    path_hdf5_temp = "temp" + hdf5_ext
    copy2(path_hdf5, path_hdf5_temp)
    return h5py.File(path_hdf5_temp, 'a')
    
def saveDataFrame(df, filename, workDir=''):
    """
    Save the provided dataframe, as filename, to workDir
    """
    with h5py.File(os.path.join(workDir, filename + '.hdf5'), 'a') as writefp:
        dataArray = df.values
        timeArray = df.index.values  
        tDataType = "S" + str(len(str(timeArray[0])))
        
        rows, cols = df.shape
        data = writefp.require_dataset(filename, (rows, cols), dtype=df[0].dtype, data=dataArray)
        time = writefp.require_dataset(filename + 'time', (rows, 1), dtype=tDataType)
        
        for idx in range(rows):
            time[idx] = str(timeArray[idx])
        writefp.flush()
    
def h5pyToDataFrame(record, workDir='', resample=False, deltaTime='5S'):
    """
    Convert an h5py (hdf) array with timestamps into a dataframe, and return it    
    
    Assumed that first column is earliest time, and last column is
    the latest time
    """
    with h5py.File(os.path.join(workDir, record + '.hdf5'), 'r') as readfp:
        origData = readfp[record]
        origTime = readfp[record + 'time']
        
        # Create a pandas DatetimeIndex with the timestamps
        time = pd.to_datetime(origTime[:len(origTime), 0])
        
        # Original Data
        origdf = pd.DataFrame(origData[:len(origTime)], index=time)
        
        if resample == True:
            # Create the new dataframe with NaN values
            index = pd.date_range(time[0], time[-1], freq=deltaTime)
            df = pd.DataFrame(index=index, columns=origdf.columns)
            
            # Merge them together
            df = origdf.combine_first(df)
        else:
            df = origdf
            
    return df