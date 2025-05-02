import xarray as xr  # Library for working with labeled multi-dimensional arrays
import numpy as np  # Library for numerical computations
import matplotlib.pyplot as plt  # Library for plotting
import glob  # Library for file handling (not used in this snippet)
import requests  # Library for making HTTP requests (used to fetch data)
import os  # Library for interacting with the operating system
import pandas as pd
import pickle # Library for saving python objects
import seaborn as sns
import lmoments3 as lm
from lmoments3 import distr
from scipy.stats import genextreme
import os

def array_mean(series):
    return series.values.mean()

def calc_rp_values(data, type = 'empirical'):
  """
  Calculates return period values for precipitation data.

  Args:
    data: A pandas Series containing DAILY precipitation values.

  Returns:
    A pandas DataFrame containing return period values for durations and return periods.
  """

  #Check Data Format, make Pandas Series
  if isinstance(data, np.ndarray):
    data = pd.Series(data)

  return_periods = [2, 5, 10, 25, 50, 100]
  duration = ['24-hr', '2-day', '3-day', '4-day', '7-day', '10-day',
       '20-day', '30-day', '45-day', '60-day']

  #Initialize DataFrame to return threshold values
  loca_thresholds_cp = pd.DataFrame(columns=return_periods, index=duration)
  num_years = int(data.size/365)

  for d in duration:
    if d.split('-')[1]=='day':
      window_length = int(d.split('-')[0])
    elif d.split('-')[1]=='hr':
      window_length = int(int(d.split('-')[0])/24)
    r_sum = data.rolling(window=window_length).sum()
    #Sort values largest to smallest
    r_sum = r_sum.sort_values(ascending=False, na_position='last')

    for rp in return_periods:
      n_days = rp * 365.25  # Convert return period in years to number of days
      exceedance_prob = 365.25 / n_days  # Probability of exceeding this threshold in a given year
      if type == 'empirical':
        # Compute threshold using percentile (inverse of exceedance probability) using first n-values of precipitation given n-year long record
        threshold = np.percentile(r_sum[:num_years], 100 - (exceedance_prob * 100))
      elif type == 'lmom':
        #Compute threshold according to L-moments fitted GEV distribution
        fit_paras = distr.gev.lmom_fit(r_sum[:num_years].fillna(0).values) #Fit the distribution to the data
        fitted_gev = distr.gev(**fit_paras)
        threshold = fitted_gev.ppf(1 - exceedance_prob)
      elif type == 'mle':
        #Compute threshold according to MLE fitted GEV distribution
        fit_paras_mle = genextreme.fit(r_sum[:num_years].fillna(0).values) #Fit the distribution to the data
        threshold = genextreme.ppf(1 - exceedance_prob, *fit_paras_mle)
      # Store the computed threshold
      loca_thresholds_cp.loc[d,rp] = threshold

  return loca_thresholds_cp

def processDS_Precip(ds, locs):
  #Daily (24-hr) Precipitation
  precip = ds['pr']

  #Slice xr for station location
  precip = precip.sel(
      lon=locs[0],
      lat=locs[1],
      method='nearest'
  )
  # Convert precipitation from kg m⁻² s⁻¹ to inches per day:
  precip = np.array(precip*86400*0.0393701)
  # Replace negative values with NaN (negative values indicate missing data)
  precip[precip<0.] = np.nan

  return precip

def annual_totals(precip):
  daily_data = pd.Series(precip, index=pd.date_range('1950-01-01', periods=len(precip)))

  annual_totals = daily_data.resample('YE').sum()
  return annual_totals
