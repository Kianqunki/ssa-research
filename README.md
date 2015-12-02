# Singular Spectrum Analysis (SSA)
These Python scripts are used to perform singular spectrum analysis on various signals retrieved from the internet to
predict critical transitions in a time series. During a critical transition the dimensionality (v here) will decrease
as the system is driven towards a critical point. During the decrease in dimensionality, variance and autocorrelation
will experience increases however these indicators are not known to work well for time series with a small sample size
or large noise.  
  
These scripts were created for research Takashi Nakamura was performing on SSA, and their results were shown at Montreal's High
Performance Computing Symposium 2015. More research on these methods have yet to be done.

To skim down on the size of the repository, the original vibration text files and PhysioNet data (> 1 GB!) were excluded.

## Results
Here are some sample graphs produced, check data sources readme to see where the data came from.

v = dimensionality
lambda_one = largest eigenvalue
sigma^2 = variance
s = sample index

SSA on a heart EKG looking at sudden cardiac death (data from PhysioNet)
![PhysioNet](/physionetdata/311indicators500.png)

SSA on Air Temperature
![AirTemp](/climatedata/AirTemp1indicators500.png)

SSA on Dow Jones index
![Financial](/financialdata/DowJonesClose1Indicators500.png)

Running a low-pass filter over the data before the SSA was also investigated, with the resulting graph below.
![lpf](/cutoffgraphs/indicatorslopes.png)

## How to run
To perform the SSA analysis run one of the following:
'Climate.py', 'financial.py' ('PhysioNet.py' doesn't have required data for reasons above)

To see the effects of adding a low-pass filter for noise, run 'CutoffGraphs.py'

'CTL.py' contains the mathematical functions used in the scripts.
'datamanager.py' contains functions to manage/parse the data in a very quick manner.

## Experimental
A folder named experimental contains modified scripts which add support for timestamps in the form of 'yyyy-mm-ddThh:mm:ss',
and add support for Pandas dataframes. This was seperated as the scripts currently run on the provided vibration examples, but
currently don't run on the previous examples used. 
