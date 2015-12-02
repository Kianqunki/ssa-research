# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 17:54:36 2015

@author: Travis
"""

from os import path
from sys import exit
from itertools import islice
import h5py
from shutil import copy2
import numpy as np

###
# Global Vars
###
hdf5_ext = ".hdf5"
script_path = path.dirname(__file__)

###
# Functions
###
def floatOrInt(fp):
    """
    Read first line to figure out if data is float or int
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
    
# Helper function to read a file in blocks
def readBlock(fp, size=65536):
    """
    Helper function to read a file in blocks
    """
    while True:
        block = fp.read(size)
        if not block: break
        yield block

# Fast method to count rows and cols in large file
# Assumed white space between each column
# Return as (rows, cols)        
def countRowCol(filename):
    """
    Fast method to count rows and cols in large file
    Assumed white space between each column
    Return as (rows, cols)     
    """
    # Rows
    with open(filename, 'r') as fp:
        rows = sum(blk.count("\n") for blk in readBlock(fp)) + 1
#    with open(filename, 'rb') as fp:
#        bufgen = takewhile(lambda x: x, (fp.read(1024*1024) for _ in repeat(None)))
#        rows = sum(buf.count(b'\n') for buf in bufgen if buf)
    
    # Cols
    with open(filename, 'r') as fp:
        line = fp.readline()
        cols = len((line.strip()).split())
    
    return rows, cols

# Read in a text data file from PhysioNet
# and store into an hdf5 file for later access
#
# sampleNumberCol = 1 indicates that the first column of the
# data file represents the sample number, used for PhysioNet data 
def processData(recordnum, data_path='data/', sampleNumberCol=0):
    """
    Read in a text data file and store into an hdf5 file for later access
    
    sampleNumberCol = 1 indicates that the first column of the
    data file represents the sample number, used for PhysioNet data 
    """
    path_txt = script_path + data_path + recordnum + ".txt"
    path_hdf5 = script_path + data_path + recordnum + ".hdf5"
    
    if not path.isfile(path_hdf5):
        # Ensure recordnum.txt exists before running code
        try: 
            fp = open(path_txt, 'r')
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
            exit(1)
        
        # Count number of lines and columns in data file for allocation
        rows, cols = countRowCol(path_txt)
        cols -= sampleNumberCol # remove sample index from cols count for PhysioNet
        
        # Allocate hdf5 dataset
        # cols-1 as we don't store the index of a sample
        dataType = floatOrInt(fp)
        datafp = h5py.File(path_hdf5, 'a')
        data = datafp.require_dataset(recordnum, (rows, cols), dtype=dataType)
        
        # Read data file into numpy array as integer
        BUF_SIZE = int(1e6) # 1 MB
        readBuf = fp.readlines(BUF_SIZE)
        readpos = np.uint32(0) # Position in read file
        datapos = np.uint32(0) # Position in dataset
        
        # Create the write buffer
        # Ensure that BUF_SIZE / 4 is int
        try:    
            writeBuf = np.zeros(shape=(len(readBuf),cols), dtype=dataType)
        except TypeError:
            print "BUF_SIZE/4 is not an int"
            fp.close()
            datafp.close()
            exit(1)
            
        while readBuf != []:
            # Put each line into the write buffer
            for line in readBuf:      
                line = np.array(line.split(), dtype=dataType)
                for idcol, sample in enumerate(islice(line, sampleNumberCol, None)):
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
            except IndexError:
                print "Index error on dataset and writeBuf"
                fp.close()
                datafp.close()
                exit(1)
            datapos = readpos
            readBuf = fp.readlines(BUF_SIZE)
        fp.close()
        datafp.close()
        
# Access the stored hdf5 data file, and pass it back as a
# temporary file to allow it to be written without affecting
# the original file
# Returns the recordnum numpy array from the temp file
def accessData(recordnum, data_path='data/'):
    """
    Access the stored hdf5 data file, and pass it back as a
    temporary file to allow it to be written without affecting
    the original file
    Returns the recordnum numpy array from the temp file
    """
    path_hdf5 = script_path + data_path + recordnum + hdf5_ext
    path_hdf5_temp = script_path + "temporaryData" + hdf5_ext
    copy2(path_hdf5, path_hdf5_temp)
    return h5py.File(path_hdf5_temp, 'a')
    