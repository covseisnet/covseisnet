#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa
"""
Network covariance matrix spectral width: natural case with preprocessing
=========================================================


This example shows how to apply preprocessing in addition to the steps done in the previous examples. The spectral width of the network covariance matrix is calculated for 24 hours of real seismic data acquired on 2010.10.14 by the Piton de la Fournaise seismic network.  This dataset contains seismovolcanic signals linked to a pre-eruptive seismic swarm and to the beginning of a co-eruptive tremor starting at 15:20 (GMT).

Four covariance matrix spectral width plots are shown to show the effects or preprocessing:
plot #1 - no preprocessing
plot #2 - preprocessing using one-bit spectral whitening
plot #3 - preprocessing using one-bit temporal normalization
plot #4 - preprocessing using both smooth spectral whitening and smooth temporal normalization
"""

import covseisnet as csn
import obspy
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# # download data from the YA Undervolc seismic network with RESIF Seismic data portal
client = Client("RESIF")
signal_duration_sec = 24 * 3600
t = UTCDateTime("2010-10-14T00:00:00.00")
list_stations = [
    "SNE",
    "FOR",
    "UV01",
    "UV02",
    "UV03",
    "UV04",
    "UV05",
    "UV06",
    "UV07",
    "UV08",
    "UV09",
    "UV10",
    "UV11",
    "UV12",
    "UV13",
    "UV14",
    "UV15",
]
stream = csn.arraystream.ArrayStream()
for sta in list_stations:
    st = client.get_waveforms("YA", sta, "00", "HHZ", t, t + signal_duration_sec)
    stream.append(st[0])

raw_stream = stream.copy()

# create plot #1 - no preprocessing

# downsample data to 25 Hz
stream.decimate(4)

# synchronize data
stream = stream.synchronize(start=t, duration_sec=signal_duration_sec, method="linear")

# calculate covariance from stream
window_duration_sec = 20
average = 15
times, frequencies, covariances = csn.covariancematrix.calculate(
    stream, window_duration_sec, average
)

# calculate spectral width
spectral_width = covariances.coherence(kind="spectral_width")

# show network covariance matrix spectral width
fig, ax = plt.subplots(1, constrained_layout=True)
img = ax.pcolormesh(times / 3600, frequencies, spectral_width.T, rasterized=True)
ax.set_ylim([0, stream[0].stats.sampling_rate / 2])
ax.set_xlabel("2010.10.14 (hours)")
ax.set_ylabel("Frequency (Hz)")
plt.colorbar(img).set_label("No preprocessing")


# create plot #2 - preprocess data with one-bit spectral whitening
stream = raw_stream.copy()

# preprocess using default setting of one-bit spectral whitening
stream.preprocess()

# downsample data to 25 Hz
stream.decimate(4)

# synchronize data
stream = stream.synchronize(start=t, duration_sec=signal_duration_sec, method="linear")


# calculate covariance from stream
times, frequencies, covariances = csn.covariancematrix.calculate(
    stream, window_duration_sec, average
)

# calculate spectral width
spectral_width = covariances.coherence(kind="spectral_width")

# show network covariance matrix spectral width
fig, ax = plt.subplots(1, constrained_layout=True)
img = ax.pcolormesh(times / 3600, frequencies, spectral_width.T, rasterized=True)
ax.set_ylim([0, stream[0].stats.sampling_rate / 2])
ax.set_xlabel("2010.10.14 (hours)")
ax.set_ylabel("Frequency (Hz)")
plt.colorbar(img).set_label("One-bit spectral whitening")


# create plot #3 - preprocess data with one-bit temporal normalization
stream = raw_stream.copy()

# preprocess using default setting of one-bit temporal normalization
stream.preprocess(domain='temporal', method='onebit')

# downsample data to 25 Hz
stream.decimate(4)

# synchronize data
stream = stream.synchronize(start=t, duration_sec=signal_duration_sec, method="linear")


# calculate covariance from stream
times, frequencies, covariances = csn.covariancematrix.calculate(
    stream, window_duration_sec, average
)

# calculate spectral width
spectral_width = covariances.coherence(kind="spectral_width")

# show network covariance matrix spectral width
fig, ax = plt.subplots(1, constrained_layout=True)
img = ax.pcolormesh(times / 3600, frequencies, spectral_width.T, rasterized=True)
ax.set_ylim([0, stream[0].stats.sampling_rate / 2])
ax.set_xlabel("2010.10.14 (hours)")
ax.set_ylabel("Frequency (Hz)")
plt.colorbar(img).set_label("One-bit temporal normalization")


# create plot #4 - preprocess data with smooth spectral whitening and temporal normalization
stream = raw_stream.copy()

# preprocess using default setting of one-bit temporal normalization
stream.preprocess(domain='spectral', method='smooth')
stream.preprocess(domain='temporal', method='smooth')

# downsample data to 25 Hz
stream.decimate(4)

# synchronize data
stream = stream.synchronize(start=t, duration_sec=signal_duration_sec, method="linear")


# calculate covariance from stream
times, frequencies, covariances = csn.covariancematrix.calculate(
    stream, window_duration_sec, average
)

# calculate spectral width
spectral_width = covariances.coherence(kind="spectral_width")

# show network covariance matrix spectral width
fig, ax = plt.subplots(1, constrained_layout=True)
img = ax.pcolormesh(times / 3600, frequencies, spectral_width.T, rasterized=True)
ax.set_ylim([0, stream[0].stats.sampling_rate / 2])
ax.set_xlabel("2010.10.14 (hours)")
ax.set_ylabel("Frequency (Hz)")
plt.colorbar(img).set_label("Smooth spectral and temporal preprocessing")