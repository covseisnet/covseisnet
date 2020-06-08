#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa
"""
Network covariance matrix spectral width: natural case
=========================================================


This example shows how to calculate the spectral width of the network covariance matrix. It makes use of real seismic data acquired the 2010.10.14 by the Piton de la Fournaise seismic network. This day of data contains seismovolcanic signals linked to the pre-eruptive seismic swarm and to the beginning of the co-eruptive tremor starting at 15:20 (GMT). This basic example does not apply any synchronization or pre-processing.

Considering a seismic network with 21 stations (using only the vertical component), the covariance is of dimensions 21 times 21. The following example use a Fourier estimation window of 50 seconds and is estimated over 40 consecutive windows. The overlapping between the Fourier windows is 50%, and the step between two consecutive averaging window is 50%.

The lower the covariance matrix spectral width (dark blue), the more coherent the seismic wavefield. Note that the sampling rate of the seismic station is 20 Hz, the dark band at 10 Hz therefore corresponds to the Nyquist frequency (and should be disregarded).
"""

import covseisnet as csn
import obspy
import matplotlib.pyplot as plt


# read ObsPy's example stream
file_seismograms = '/path/to/data'
stream = obspy.read(file_seismograms)

# calculate covariance from stream
window_duration_sec = 50
average = 40
times, frequencies, covariances = csn.covariancematrix.calculate(
    stream, window_duration_sec, average
)

# calculate spectral width
spectral_width = covariances.coherence(kind="spectral_width")

# show network covariance matrix spectral width
fig, ax = plt.subplots(1, constrained_layout=True)
img = ax.pcolormesh(times, frequencies, spectral_width.T, rasterized=True)
ax.set_ylim([0, stream[0].stats.sampling_rate / 2])
ax.set_xlabel("Times (seconds)")
ax.set_ylabel("Frequency (Hz)")
plt.colorbar(img).set_label("Covariance matrix spectral width")
