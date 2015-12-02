# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 20:52:12 2015

@author: Travis
"""

from numpy import *
import CTL

###
# Runtime
###
recordList = ["31"]
order = 10
#lowPassFilter = [None, [order, 0.1], [order, 0.25], [order, 0.5]]
lowPassFilter = [[order, 0.05]]
for lpf in lowPassFilter:
    CTL.fullAnalysis(recordList, 'physionetdata/', lpf, 1)