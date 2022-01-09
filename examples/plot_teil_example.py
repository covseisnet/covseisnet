#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa
"""
Beamforming and Network Response Function: Teil earthquake example
=========================================================


This example shows how to calculate the network response function and do beamforming on real seismic data acquired on the Teil earthquake.

Four covariance matrix spectral width plots are shown to show the effects or preprocessing:
plot #1 - spectral width
plot #2 - network response function
"""

import covseisnet as csn
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import numpy as np


# # download data near Teil with RESIF Seismic data portal
# client = Client("RESIF")
# signal_duration_sec = 600
# t = UTCDateTime("2019-11-11T10:50:00.00")
# list_stations = [
#     "BANN",
#     "OGCB",
#     "OGCC",
#     "OGDF",
#     "SAUF",
# ]
# stream = csn.arraystream.ArrayStream()
# for sta in list_stations:
#     st = client.get_waveforms("FR", sta, "00", "HHZ", t, t + signal_duration_sec)
#     stream.append(st[0])


# from file (temp)
file = "C:/sync/AlpArray/teil/input/mseed/2019-11-11T1050HHZ.mseed"
stream = csn.arraystream.read(file)
# remove distant stations
for tr in stream.select(station="BALS"):
    stream.remove(tr)
for tr in stream.select(station="BSTF"):
    stream.remove(tr)
for tr in stream.select(station="GRN"):
    stream.remove(tr)
for tr in stream.select(station="MLYF"):
    stream.remove(tr)
for tr in stream.select(station="OGCN"):
    stream.remove(tr)
for tr in stream.select(station="OGLC"):
    stream.remove(tr)
for tr in stream.select(station="OGS2"):
    stream.remove(tr)
for tr in stream.select(station="OGS3"):
    stream.remove(tr)
for tr in stream.select(station="OGVG"):
    stream.remove(tr)


# Settings
# --------
# Paths
traveltime_filepath = "./data"

# Read data
# ---------
## frequency limits for filtering (depend on the target signal)
band = [0.1, 50.]
low_pass = 0.1
high_pass = 50

duration_sec =  10 * 60 # length of stream to be proceseed, in seconds
sampling_rate = 100 # resampling frequency

window_duration_sec = 35 #orig was 50
average = 6 #orig was 40
overlap=0.5
preproc_spectral_secs = window_duration_sec * average * overlap

# Correlation Smoothing
sigma = 30 

# grid extent
lon_min = 4.15000
lon_max = 5.69301
lat_min = 43.9300
lat_max = 44.9843
depth_min = -1.3
depth_max = 58.7

                
## merge traces to have one trace per station
stream.merge(method=1,fill_value='interpolate',interpolation_samples=-1)
      
## set start time
start = UTCDateTime(stream[0].stats.starttime)



## synchronize traces in the stream
stream = stream.synchronize(start, duration_sec, method="linear")
nsta = len(stream)



### filtering
stream.detrend(type='demean')
stream.detrend(type='linear')
stream.filter(type='bandpass', freqmin=band[0], freqmax=band[1])  
        
## Preprocess : (can try different pre-processing)
stream.preprocess(domain='spectral', method='onebit', window_duration_sec=preproc_spectral_secs)
stream.trim(start, start+duration_sec) #preprocessing can add artifacts to ends




# Calculate coherence
# -------------------


times, frequencies, covariances = csn.covariancematrix.calculate(
    stream, window_duration_sec, average
)




# Spectral width
spectral_width = covariances.coherence(kind="spectral_width")
   

# Eigenvector decomposition
covariance_1st = covariances.eigenvectors(covariance=True, rank=0)

# Extract cross-correlations from the covariance matrix filtered by the 1st eigenvector
correlation =  csn.correlationmatrix.cross_correlation(covariance_1st, sampling_rate)
   
# load traveltime grids from file
teil_traveltimes = csn.traveltime.TravelTime(stream, traveltime_filepath)

# Initiate beam object and set geographical extent of grid
nwin = correlation.nwin() # number of time windows
teil_beam = csn.beam.Beam(nwin, teil_traveltimes)
teil_beam.set_extent(lon_min, lon_max, lat_min, lat_max, depth_min, depth_max)

# Loop through all windows and calculate likelihood and nrf
for i in range(0, nwin):
    print('Processing window', i+1, 'of', nwin)

    correl = correlation[i]

    # Filter correlation
    correl = correl.bandpass(low_pass, high_pass, sampling_rate)
    
    # Smooth correlation
    correl = correl.hilbert_envelope()
    correl = correl.smooth(sigma=sigma) #default sigma is 5
           
    teil_beam.calculate_likelihood(correl.T, sampling_rate, i)
    teil_beam.calculate_nrf(i) 
    
    beam_max = teil_beam.max_likelihood(i)
    print("Maximum likelihood of", round(beam_max[3],3), "at", round(beam_max[0],4), "\N{DEGREE SIGN},", round(beam_max[1],4), "\N{DEGREE SIGN},", round(beam_max[2],1), "km")



# plot nrf and spectral width
fig, ax = plt.subplots(2, constrained_layout=True, figsize=(10,8)) #stretched plot

ax[0].plot(np.linspace(0, 10, nwin), teil_beam.nrf, "k")
ax[0].set_title("NRF")
ax[0].set_ylabel("Network response function, unnormalized")
ax[0].set_xlabel("Minutes")
ax[0].set_xlim(0,10)

img = ax[1].imshow(spectral_width.T, origin='lower', cmap='jet_r', interpolation='none', extent=[0, duration_sec/60, 0, sampling_rate], aspect='auto')
ax[1].set_ylim([0.1, stream[0].stats.sampling_rate / 2])
ax[1].set_xlabel("Minutes")
ax[1].set_ylabel("Frequency (Hz)")
ax[1].set_yscale("log")
ax[1].set_title("Spectral Width")
plt.colorbar(img).set_label("Covariance matrix spectral width")