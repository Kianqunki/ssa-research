# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 20:52:12 2015

@author: Travis
"""

from CTL import *
import datamanager as datamgr
#import matplotlib # For westgrid
#matplotlib.use('Agg') # For westgrid
import matplotlib.pyplot as plt


###
# Runtime
###
workDir = 'vibrationdata/'
record = 'viadata'

"""
numEigens = [6, 7, 8, 9, 10, 11, 12, 13]
winSizes = [500, 1000, 5000]
ssaWinSizes = [500, 1000, 5000, 10000, 50000]
"""
numEigens = [6]
winSizes = [1000]
ssaWinSizes = [5000, 10000]
minError = 1e-4

datamgr.processData(record, workDir, 1)
origdf = datamgr.h5pyToDataFrame(record, workDir, resample=True, deltaTime='1T')

workDir = workDir + 'results/' # save everything in results
for numEigen in numEigens:
    for winSize in winSizes:
        title = 'AverageInterpolation'
        df = origdf.copy()
        df_filled = nanAvgFill(df)
        df_ssa = ssaInterpolation(df_filled, numEigen, winSize, minError)
        """
        for col in df:
            fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(15,15))
            df[col].plot(ax=axes[0]); axes[0].set_title('Original Data')
            df_filled[col].plot(ax=axes[1]); axes[1].set_title('Missing Points Interpolation w/ Average')
            df_ssa[col].plot(ax=axes[2]); axes[2].set_title('SSA Interpolation w/ eigens=' + str(numEigen)
                                    + ' winSize=' + str(winSize))
            #plt.show()
            fig.savefig(workDir + title + '_' + str(col) + 'col_' + str(numEigen) + 'numEigen_' + str(winSize) + 'winSize' + '.png', dpi=300)
        """
        #datamgr.saveDataFrame(df, title + '_' + str(numEigen) + 'numEigen_' + str(winSize) + 'winSize', workDir)
        for ssaWinSize in ssaWinSizes:
            analysis(df_ssa, title + '_' + str(numEigen) + 'numEigen', workDir, ssaWinSize)
        
        title = 'HybridInterpolation'
        df = origdf.copy()
        df_filled = nanHybridFill(df)
        df_ssa = ssaInterpolation(df_filled, numEigen, winSize, minError)
        """
        for col in df:
            fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(15,15))
            df[col].plot(ax=axes[0]); axes[0].set_title('Original Data')
            df_filled[col].plot(ax=axes[1]); axes[1].set_title('Missing Points Interpolation w/ Hybrid')
            df_ssa[col].plot(ax=axes[2]); axes[2].set_title('SSA Interpolation w/ eigens=' + str(numEigen)
                                    + ' winSize=' + str(winSize))
            #plt.show()
            fig.savefig(workDir + title + '_' + str(col) + 'col_' + str(numEigen) + 'numEigen_' + str(winSize) + 'winSize' + '.png', dpi=300)
        """
        #datamgr.saveDataFrame(df, title + '_' + str(numEigen) + 'numEigen_' + str(winSize) + 'winSize', workDir)
        for ssaWinSize in ssaWinSizes:
            analysis(df_ssa, title + '_' + str(numEigen) + 'numEigen', workDir, ssaWinSize)
        
        title = 'LinearInterpolation'
        df = origdf.copy()
        df_filled = df.interpolate()
        df_ssa = ssaInterpolation(df_filled, numEigen, winSize, minError)
        """
        for col in df:
            fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(15,15))
            df[col].plot(ax=axes[0]); axes[0].set_title('Original Data')
            df_filled[col].plot(ax=axes[1]); axes[1].set_title('Missing Points Interpolation w/ Pandas Linear')
            df_filled[col].plot(ax=axes[2]); axes[2].set_title('SSA Interpolation w/ eigens=' + str(numEigen)
                                    + ' winSize=' + str(winSize))
            #plt.show()
            fig.savefig(workDir + title + '_' + str(col) + 'col_' + str(numEigen) + 'numEigen_' + str(winSize) + 'winSize' + '.png', dpi=300)
        """
        #datamgr.saveDataFrame(df, title + '_' + str(numEigen) + 'numEigen_' + str(winSize) + 'winSize', workDir)
        for ssaWinSize in ssaWinSizes:        
            analysis(df_ssa, title + '_' + str(numEigen) + 'numEigen', workDir, ssaWinSize)

