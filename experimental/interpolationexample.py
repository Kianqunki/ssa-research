# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 20:52:12 2015

@author: Travis
"""

from CTL import *
import datamanager as datamgr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import h5py
import os


###
# Runtime
###
workDir = 'testdata/'
record = 'testcase'

numEigen = 6
winSize = 50
minError = 1e-10

datamgr.processData(record, workDir, 1)
df = datamgr.h5pyToDataFrame(record, workDir, resample=True, deltaTime='1T')

title = 'AverageInterpolation'
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(15,15))
df[0].plot(ax=axes[0]); axes[0].set_title('Original Data')
df = nanAvgFill(df)
df[0].plot(ax=axes[1]); axes[1].set_title('Missing Points Interpolation w/ Avg')
df = ssaInterpolation(df, numEigen, winSize, minError)
df[0].plot(ax=axes[2]); axes[2].set_title('SSA Interpolation w/ eigens=' + str(numEigen)
                            + ' winSize=' + str(winSize))
plt.show()
fig.savefig(workDir + record + title + '.png', dpi=300)
datamgr.saveDataFrame(df, title + '_' + str(winSize) + 'winSize_' + str(numEigen) + 'eigen', workDir)

title = 'HybridInterpolation'
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(15,15))
df = datamgr.h5pyToDataFrame(record, workDir, resample=True, deltaTime='1T')
df[0].plot(ax=axes[0]); axes[0].set_title('Original Data')
df = nanHybridFill(df)
df[0].plot(ax=axes[1]); axes[1].set_title('Missing Points Interpolation w/ Avg Hybrid')
df = ssaInterpolation(df, numEigen, winSize, minError)
df[0].plot(ax=axes[2]); axes[2].set_title('SSA Interpolation w/ eigens=' + str(numEigen)
                            + ' winSize=' + str(winSize))
plt.show()
fig.savefig(workDir + record + title + '.png', dpi=300)
datamgr.saveDataFrame(df, title + '_' + str(winSize) + 'winSize_' + str(numEigen) + 'eigen', workDir)

title = 'LinearInterpolation'
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(15,15))
minError = 0.01
datamgr.processData(record, workDir, 1)
df = datamgr.h5pyToDataFrame(record, workDir, resample=True, deltaTime='1T')
df[0].plot(ax=axes[0]); axes[0].set_title('Original Data')
df = df.interpolate()
df[0].plot(ax=axes[1]); axes[1].set_title('Missing Points Interpolation w/ Linear')
df = ssaInterpolation(df, numEigen, winSize, minError)
df[0].plot(ax=axes[2]); axes[2].set_title('SSA Interpolation w/ eigens=' + str(numEigen)
                            + ' winSize=' + str(winSize))
plt.show()
fig.savefig(workDir + record + title + '.png', dpi=300)
datamgr.saveDataFrame(df, title + '_' + str(winSize) + 'winSize_' + str(numEigen) + 'eigen', workDir)