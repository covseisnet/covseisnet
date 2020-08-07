#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa
"""
Single-station covariance matrix spectral width
===============================================


This example shows how to calculate the spectral width of the interchannel spectral covariance matrix. It makes use of the obspy example trace available when installing obspy. This basic example does not apply any synchronization or pre-processing.

Considering a single seismic station with three channels (NS, EW, and Z), the covariance is of dimensions 3 times 3. The following example use a Fourier estimation window of 1 second and is estimated over 5 consecutive windows. The overlapping between the Fourier windows is 50%, and the step between two consecutive averaging window is 50%.

The lower the covariance matrix spectral width (yellow), the more polarized the seismic wavefield (into a given direction). Note that the sampling rate of the seismic station is 100 Hz, the dark band at 50 Hz therefore corresponds to the Nyquist frequency (and should be disregarded).
"""

import covseisnet as csn
import obspy
import matplotlib.pyplot as plt


# read ObsPy's example stream
stream = obspy.read()

# calculate covariance from stream
window_duration_sec = 1.0
average = 5
times, frequencies, covariances = csn.covariancematrix.calculate(
    stream, window_duration_sec, average
)

# calculate spectral width
spectral_width = covariances.coherence(kind="spectral_width")

# show covariance at first time window and first frequency
fig, ax = plt.subplots(1, constrained_layout=True)
img = ax.pcolormesh(
    times, frequencies, spectral_width.T, rasterized=True, cmap="viridis_r"
)
ax.set_ylim([0, stream[0].stats.sampling_rate / 2])
ax.set_xlabel("Times (seconds)")
ax.set_ylabel("Frequency (Hz)")
plt.colorbar(img).set_label("Covariance matrix spectral width")
