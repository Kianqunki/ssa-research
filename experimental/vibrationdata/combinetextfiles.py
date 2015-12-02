# -*- coding: utf-8 -*-
"""
Created on Fri Jul 03 18:04:20 2015

@author: Travis
"""

"""
Expected format is
1,1-Feb-15,02:21:21,6038.2,5640.4,5262.0,3.52,8.82
samplenumber,day-mth-yr,hour:min:sec,data0,data1 etc.

All files in the script directory will be ignored
ex: combinetextfiles.py isn't read
"""

import os
from datetime import datetime

data = []

# read data entries into a python list
for dirpath, dirnames, filenames in os.walk('.'):   
    for filename in filenames:
        # ignore script path and bad data folder
        if dirpath != '.':
            with open(os.path.join(dirpath, filename), 'r') as fp:
                # Remove blank/whitespaced lines
                lines = (line.rstrip() for line in fp)
                lines = (line for line in lines if line)
                
                buf = []
                for line in lines:
                    line = line.split(',')

                    # deal with entries split onto multiple lines
                    if len(line) < 8:
                        buf = buf + line
                    
                    # remove sample index + reformat timestamp
                    # an expected example string is
                    # 1,1-Feb-15,02:21:21,6038.2,5640.4,5262.0,3.52,8.82
                    # line = ['1', '1-Feb-15', '02:21:21', '6038.2', etc.]
                    # time = '1-Feb-15 02:21:21'
                    if len(line) == 8:
                        for idx, item in enumerate(line):
                            if item in ['', 'g0.?', '.']:
                                line[idx] = 'nan' # set invalid entries to NaN
                        
                        # Convert datetimes to ISO format
                        time = datetime.strptime(line[1] + ' ' + line[2], "%d-%b-%y %H:%M:%S")
                        line = line[2:]
                        line[0] = datetime.isoformat(time)
                        
                        # Save data to our array
                        data.append(line)
                        buf = []
                        
# sort data by earliest date to latest date
data = sorted(data, key=lambda row: row[0])

# count the different delta ts for inspection with variable explorer in Spyder
dtimes = {} # empty dictionary
for idx, entry in enumerate(data):
    if idx != 0:
        d2 = datetime.strptime(data[idx][0], "%Y-%m-%dT%H:%M:%S")
        d1 = datetime.strptime(data[idx - 1][0], "%Y-%m-%dT%H:%M:%S")
        dtime = str(d2 - d1)
        
        if dtime in dtimes:
            dtimes[dtime] += 1
        else:
            dtimes[dtime] = 1
            
              
# write data entries into central text file
with open('viadata.txt', 'wb') as writefp:
    for entry in data:
        writefp.write(' '.join(entry) + '\n')