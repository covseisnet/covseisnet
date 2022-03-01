#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa
"""
Beamforming and Network Response Function: Teil earthquake example
=========================================================


This example shows how to calculate the network response function and likelihood location on seismic data of the Teil earthquake.

The plot, from top to bottom, shows
plot #1 - network response function
plot #2 - spectral width
plot #3 - waveform from the vertical channel of station OGDF

"""

import covseisnet as csn
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import RegularPolygon
from obspy import UTCDateTime, read_inventory
from obspy.clients.fdsn import Client
import numpy as np

# download seismic data near Teil using the RESIF Seismic data portal
client = Client("RESIF")
duration_sec = 30 * 60 # length of stream to be processed, in seconds
t = UTCDateTime("2019-11-11T10:30:00.00")
list_stations = [
    "BANN",
    "OGCB",
    "OGCC",
    "OGDF",
    "SAUF",
]
stream = csn.arraystream.ArrayStream()
for sta in list_stations:
    st = client.get_waveforms("FR", sta, "00", "HHZ", t, t + duration_sec)
    stream.append(st[0])

# choose waveform from channel HHZ of station OGDF for plotting
stream_plot = stream.select(station="OGDF", channel="HHZ")
stream_plot.detrend(type='demean')
stream_plot.detrend(type='linear')
band = [2, 50.]
stream_plot.filter(type='bandpass', freqmin=band[0], freqmax=band[1])
trace_plot = stream_plot[0].data


# Settings
# --------
traveltime_filepath = "./data/Teil" # path to traveltime grids
metadata_filepath = "./data/Teil/*.xml"

# frequency limits for filtering (depends on the target signal)
band = [0.1, 50.]
low_pass = 0.1
high_pass = 50

sampling_rate = 100 # resampling frequency
window_duration_sec = 35
average = 6
overlap=0.5
preproc_spectral_secs = window_duration_sec * average * overlap
sigma = 30 # Correlation Smoothing

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

## filtering
stream.detrend(type='demean')
stream.detrend(type='linear')
stream.filter(type='bandpass', freqmin=band[0], freqmax=band[1])  
        
## Preprocess
stream.preprocess(domain='spectral', method='onebit', window_duration_sec=preproc_spectral_secs)
# stream.trim(start, start+duration_sec) #preprocessing can add artifacts to ends




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
           
    teil_beam.calculate_likelihood(correl, sampling_rate, i)
    teil_beam.calculate_nrf(i) 
    
    beam_max = teil_beam.max_likelihood(i)
    print("Maximum likelihood of", round(beam_max[3],3), "at", round(beam_max[0],4), "\N{DEGREE SIGN},", round(beam_max[1],4), "\N{DEGREE SIGN},", round(beam_max[2],1), "km")



# plot nrf and spectral width and waveforms from closest station
fig, ax = plt.subplots(3, constrained_layout=True, figsize=(10,8)) #stretched plot

ax[0].plot(np.linspace(0, duration_sec/60, nwin), teil_beam.nrf, "k")
ax[0].set_title("NRF")
ax[0].set_ylabel("Network response function, unnormalized")
ax[0].set_xlim(0,duration_sec/60)

img = ax[1].imshow(spectral_width.T, origin='lower', cmap='jet_r', interpolation='none', extent=[0, duration_sec/60, 0, sampling_rate], aspect='auto')
ax[1].set_ylim([0.1, stream[0].stats.sampling_rate / 2])
ax[1].set_ylabel("Frequency (Hz)")
ax[1].set_yscale("log")
ax[1].set_title("Spectral Width")
plt.colorbar(img).set_label("Covariance matrix spectral width")

ax[2].plot(np.linspace(0,duration_sec/60, len(trace_plot)), trace_plot, "k")
ax[2].set_title("Vertical channel of station OGDF")
ax[2].set_xlabel("Minutes")
ax[2].set_ylabel("OGDF.HHZ (counts)")
ax[2].set_xlim([0, duration_sec/60])

# create dictionary of station metadata from station xml
inv = read_inventory(metadata_filepath)
net = {"lat":[],"lon":[]}
for tr in stream:
    inv_sel = inv.select(station=tr.stats.station)
    net["lat"].append(inv_sel[0][0].latitude)
    net["lon"].append(inv_sel[0][0].longitude)
    
# Project max likelihood to a plane along each axis
volume = teil_beam.likelihood[25,:,:,:]
yz_slice = np.max(volume, axis=0)
xz_slice = np.max(volume, axis=1)
xy_slice = np.max(volume, axis=2)
   
# colorbar_min = np.median(volume)
colorbar_min = 0.9*np.max(volume)
colorbar_max = np.max(volume)

fig, ax = plt.subplots(1, constrained_layout=True, figsize=(10,5), dpi=200)
img_xy = ax.imshow(xy_slice.T, interpolation='none', origin = 'lower',  cmap = 'turbo', aspect="auto", extent=[lon_min,lon_max,lat_min,lat_max], norm=colors.LogNorm(vmin=colorbar_min,vmax=colorbar_max))
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.set_title("Likelihood location, map view")
for x,y in zip(net["lon"],net["lat"]):
    triangle = RegularPolygon((x,y), facecolor = 'white', edgecolor ='black', numVertices=3, radius = 0.02)
    ax.add_patch(triangle)
plt.colorbar(img_xy).set_label("Likelihood")

fig, ax = plt.subplots(1, constrained_layout=True, figsize=(9.6,3), dpi=200)
img_xz = ax.imshow(xz_slice.T, interpolation='none', origin = 'upper',  cmap = 'turbo', aspect="auto", extent=[lon_min,lon_max, depth_max, depth_min], norm=colors.LogNorm(vmin=colorbar_min,vmax=colorbar_max))
ax.set_xlabel("Longitude")
ax.set_ylabel("Depth (km)")
ax.set_title("Likelihood location, depth view")
plt.colorbar(img_xz).set_label("Likelihood")