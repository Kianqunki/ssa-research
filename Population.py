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
# Had to lower embedding dimension and N due to small datasets
embeddingDim = 10
N = 2

recordList = ["SockeyeSalmon", "ElephantSeal", "ChumSalmon", "Beaver"]
order = 2
lowPassFilter = [None, [order, 0.1], [order, 0.25], [order, 0.5]]
for lpf in lowPassFilter:
    CTL.fullAnalysis(recordList, 'populationdata/', lpf, 0, embeddingDim, N)