#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa
"""
Beamforming and Network Response Function Example: Tremor at the Piton de la Fournaise
=========================================================
This example shows how to calculate the network response function and 
likelihood location on seismic data associated with the eruption at Piton de la
Fournaise on December 21, 2021.

In this example we select the 3rd window for plotting the likelihood function.

The first plot shows,from top to bottom, 1) waveforms from the vertical channel
of station FOR near the eruption, 2) the network response function, and 3) the 
spectral width.

The second plot shows the normalized likelihood function plotted onto shaded 
topography in map view. The third and fourth plots show the same thing but in
depth view longitudinally and latitudinally respectively.
"""

import covseisnet as csn
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import matplotlib.colors as mcolors
from matplotlib.colors import LightSource
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
import numpy as np
import rasterio
import scipy.interpolate
# import elevation

# # download seismic data near using the RESIF Seismic data portal
client = Client("RESIF")
duration_sec = 60 * 60 # length of stream to be processed, in seconds
t = UTCDateTime("2021-12-21T23:40:00.00")
list_stations = [
    "BON",
    "BOR",
    "CRA",
    "CSS",
    "DSO",
    "ENO",
    "FJS",
    "FOR",
    "FRE",
    "GBS",
    "GPN",
    "GPS",
    "HDL",
    "HIM",
    "PHR",
    "PRA",
    # "PVD", #not available publicly
    # "RVA", #not available publicly
    "SNE",
]
# stream = csn.arraystream.ArrayStream()
# for sta in list_stations:
#     st = client.get_waveforms("PF", sta, "*", "*HZ", t, t + duration_sec)
#     stream.append(st[0])

# download metadata
inv = client.get_stations(network="PF", location="*", channel="*HZ", starttime=t, endtime=t + duration_sec)


# from file (temp)
# file = "/mnt/c/sync/volcanoes/PdF/input/mseed/2021-12-21T23.mseed"
file = "C:/sync/volcanoes/PdF/input/mseed/2021-12-21T23-1.mseed"
stream = csn.arraystream.read(file)
stream.trim(t, t+duration_sec)
# for tr in stream.select(station="BOR"): #BOR ok
#     stream.remove(tr)
for tr in stream.select(station="RVA"): # causes issues if you include RVA
    stream.remove(tr)
for tr in stream.select(station="PVD"): # also causes issue if PVD incl
    stream.remove(tr)

# decimate data to 20Hz as we are only looking at low frequency tremors
stream.decimate(5)
# stream.resample(20)



# Settings
# --------
traveltime_filepath = "./data/PdF" # path to traveltime grids


# frequency limits for filtering (depends on the target signal)
low_pass = 0.5
high_pass = 10.0

# optimized for VT earthquakes
# window_duration_sec = 12
# average = 20

# optimized for tremors
window_duration_sec = 60 
average = 30

overlap=0.5
sampling_rate = stream[0].stats.sampling_rate # sampling rate, assumes all streams have the same
preproc_spectral_secs = window_duration_sec * average * overlap
sigma = 30 # Correlation Smoothing

# grid extent Longitude: 55.67° to 55.81° (145 points), Latitude: -21.3° to -21.2° (110 points)
lon_min = 55.67 
lon_max = 55.81
lat_min = -21.3
lat_max = -21.2
depth_min = -2.6 # edifice summit (km)
depth_max = 3.0 # max depth of traveltime grids (km) 





                
## merge traces to have one trace per station
stream.merge(method=1,fill_value='interpolate',interpolation_samples=-1)
      
## synchronize traces in the stream
stream = stream.synchronize(t, duration_sec, method="linear")

## filtering
stream.detrend(type='demean')
stream.detrend(type='linear')
stream.filter(type='bandpass', freqmin=low_pass, freqmax=high_pass, zerophase=True)  

# choose waveform from channel HHZ of station OGDF for plotting
trace_plot = stream.select(station="FOR", channel="HHZ")[0].data
# stream_plot.detrend(type='demean')
# stream_plot.detrend(type='linear')
# stream_plot.filter(type='bandpass', freqmin=low_pass, freqmax=high_pass)
# trace_plot = stream_plot[0].data

        
## Preprocess
stream.preprocess(domain='spectral', method='onebit', window_duration_sec=preproc_spectral_secs)
stream.trim(t, t+duration_sec) #preprocessing can add artifacts to ends




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
PdF_traveltimes = csn.traveltime.TravelTime(stream, traveltime_filepath)

# Initiate beam object and set geographical extent of grid
nwin = correlation.nwin() # number of time windows
PdF_beam = csn.beam.Beam(nwin, PdF_traveltimes)
PdF_beam.set_extent(lon_min, lon_max, lat_min, lat_max, depth_min, depth_max)

# Loop through all windows and calculate likelihood and nrf
for i in range(0, nwin):
    print('Processing window', i+1, 'of', nwin)

    correl = correlation[i]

    # Filter correlation
    correl = correl.bandpass(low_pass, high_pass, sampling_rate)
    
    # Smooth correlation
    correl = correl.hilbert_envelope()
    correl = correl.smooth(sigma=sigma) #default sigma is 5
           
    PdF_beam.calculate_likelihood(correl, sampling_rate, i)
    PdF_beam.calculate_nrf(i) 
    
    beam_max = PdF_beam.max_likelihood(i)
    print("Maximum likelihood of", round(beam_max[3],3), "at", round(beam_max[0],4), "\N{DEGREE SIGN},", round(beam_max[1],4), "\N{DEGREE SIGN},", round(beam_max[2],1), "km")



# plot nrf and spectral width and waveforms from closest station
fig, ax = plt.subplots(3, constrained_layout=True, figsize=(10,8)) #stretched plot

ax[0].plot(np.linspace(0,duration_sec/60, len(trace_plot)), trace_plot, "k")
ax[0].set_title("Vertical channel of station FOR")
ax[0].set_xlabel("Minutes")
ax[0].set_ylabel("FOR.HHZ (counts)")
ax[0].set_xlim([0, duration_sec/60])

ax[1].plot(np.linspace(0, duration_sec/60, nwin), PdF_beam.nrf, "k")
ax[1].set_title("NRF")
ax[1].set_ylabel("Network response function, unnormalized")
ax[1].set_xlim(0,duration_sec/60)

img = ax[2].imshow(spectral_width.T, origin='lower', cmap='jet_r', interpolation='none', extent=[0, duration_sec/60, 0, sampling_rate], aspect='auto')
ax[2].set_ylim([0.1, stream[0].stats.sampling_rate / 2])
ax[2].set_ylabel("Frequency (Hz)")
ax[2].set_yscale("log")
ax[2].set_title("Spectral Width")
plt.colorbar(img).set_label("Covariance matrix spectral width")



for ii in range(0,6):
    # Choose the 3nd window for plotting likelihood
    likelihood_xyz = PdF_beam.likelihood[ii,:,:,:]
    
    
    # Take slices at point of max likelihood
    i_max,j_max,k_max = np.unravel_index(likelihood_xyz.argmax(), likelihood_xyz.shape)
    likelihood_xy = likelihood_xyz[:, :, k_max]
    likelihood_xz = likelihood_xyz[:, j_max]
    likelihood_yz = likelihood_xyz[i_max]
    
    # Normalize likelihood between 0 and 1
    likelihood_xy =  (likelihood_xy - likelihood_xy.min()) / (likelihood_xy.max()-likelihood_xy.min())
    likelihood_xz =  (likelihood_xz - likelihood_xz.min()) / (likelihood_xz.max()-likelihood_xz.min())
    likelihood_yz =  (likelihood_yz - likelihood_yz.min()) / (likelihood_yz.max()-likelihood_yz.min())
    
    
    
    # Download DEM and interpolate to grid
    # dem_path = './DEM.tif'
    # elevation.clip(bounds=(lon_min, lat_min, lon_max, lat_max), output=dem_path, product='SRTM3')
    
    dem_path = 'C:/sync/volcanoes/PdF/input/DEM.tif'
    dem=rasterio.open(dem_path) #open dem file
    
    #data
    dem1=dem.read(1)
    dem1 = np.where(dem1 == -32768, 0, dem1)
    
    
    # old dem grid
    x, y = np.mgrid[0:119:120j, 0:167:168j]
    
    # new dem grid, with dimensions matching our traveltime grid
    x2, y2 = np.mgrid[0:119:110j, 0:167:145j]
    
    # interpolate onto the new grid
    dem2 = scipy.interpolate.griddata((x.ravel(), y.ravel()), dem1.ravel(), (x2, y2), method='linear')
    dem_x = -1*dem2[i_max,:]/1000 #dem along xz slice, convert to km
    dem_y = np.flip(-1*dem2[:,j_max]/1000) #dem along yz slice, convert to km
    
    
    # create dictionary of station metadata from station xml
    net = {"lat":[],"lon":[]}
    for tr in stream:
        inv_sel = inv.select(station=tr.stats.station)
        net["lat"].append(inv_sel[0][0].latitude)
        net["lon"].append(inv_sel[0][0].longitude)
    
    
    
    # create custom discrete colormap for likelihood
    low = 4
    levels = 12
    n_colours = 16 #number of discrete colours to split colourbar into
    custom_cmap = plt.cm.get_cmap('RdYlBu_r')(np.linspace(0, 1, levels))
    for i in range(low):
        custom_cmap[i, :] = [1, 1, 1, 1]
    for i in range(1, levels):
        custom_cmap[i, -1] = np.sqrt(i / levels)
    custom_cmap = mcolors.LinearSegmentedColormap.from_list('RdYlBu_r', custom_cmap, N=n_colours)
    
    # create custome discrete colormap for topography
    levels = 12
    custom_cmap_dem = plt.cm.get_cmap('Greys')(np.linspace(0.2, 0.5, levels))
    custom_cmap_dem = mcolors.LinearSegmentedColormap.from_list('Greys', custom_cmap_dem)
    
    # prepare shaded topography
    # create light source object.
    ls = LightSource(azdeg=315, altdeg=45)
    # shade data, creating an rgb array.
    rgb = ls.shade(dem2, custom_cmap_dem)
    # rgb = ls.shade(dem2, custom_cmap_dem, vert_exag=100, dx=0.122, dy=0.122, fraction=1, )
    
    # plot map view
    fig, ax = plt.subplots(1, constrained_layout=True, figsize=(10,5), dpi=100)
    img_xy = ax.imshow(likelihood_xy.T, interpolation='none', origin = 'lower',  cmap = custom_cmap, aspect="auto", extent=[lon_min,lon_max,lat_min,lat_max], vmin=0.5)
    img_xy_dem = ax.imshow(rgb, interpolation='none', alpha=0.35, cmap=plt.get_cmap('Greys'), aspect="auto", extent=[lon_min,lon_max,lat_min,lat_max])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Likelihood location, map view")
    for x,y in zip(net["lon"],net["lat"]):
        triangle = RegularPolygon((x,y), facecolor = 'white', edgecolor ='black', numVertices=3, radius = 0.002)
        ax.add_patch(triangle)
    plt.colorbar(img_xy).set_label("Likelihood")
    # plt.savefig("xy.png")

# plot depth view
fig, ax = plt.subplots(1, constrained_layout=True, figsize=(9.6,3), dpi=100)
img_xz = ax.imshow(likelihood_xz.T, interpolation='none', origin = 'upper',  cmap = custom_cmap, aspect="auto", extent=[lon_min,lon_max, depth_max, depth_min], vmin=0.5)
ax.fill_between(np.linspace(lon_min, lon_max, len(dem_x)), depth_min, dem_x,facecolor='w', edgecolor='k', lw=.4) #crop out data above topo
ax.set_xlabel("Longitude")
ax.set_ylabel("Depth (km)")
ax.set_title("Likelihood location, depth view")
plt.colorbar(img_xz).set_label("Likelihood")
# plt.savefig("xz.png")


fig, ax = plt.subplots(1, constrained_layout=True, figsize=(9.6,3), dpi=100)
img_yz = ax.imshow(likelihood_yz.T, interpolation='none', origin = 'upper',  cmap = custom_cmap, aspect="auto", extent=[lat_min,lat_max, depth_max, depth_min], vmin=0.5)
ax.fill_between(np.linspace(lat_min, lat_max, len(dem_y)), depth_min, dem_y,facecolor='w', edgecolor='k', lw=.4) #crop out data above topo
ax.set_xlabel("Latitude")
ax.set_ylabel("Depth (km)")
ax.set_title("Likelihood location, depth view")
plt.colorbar(img_yz).set_label("Likelihood")
# plt.savefig("yz.png")
