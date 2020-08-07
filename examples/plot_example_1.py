#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa
"""
Single-station covariance matrix
================================


This example shows how to calculate the interchannel spectral covariance matrix. It makes use of the obspy example trace available when installing obspy. This basic example does not apply any synchronization or pre-processing.

Considering a single seismic station with three channels (NS, EW, and Z), the covariance is of dimensions 3 times 3. The following example use a Fourier estimation window of 1 second and is estimated over 5 consecutive windows.
"""

import covseisnet as csn
import obspy
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import MaxNLocator

# read ObsPy's example stream
stream = obspy.read()
channels = [s.stats.channel for s in stream]

# calculate covariance from stream
window_duration_sec = 1.0
average = 5
times, frequencies, covariances = csn.covariancematrix.calculate(
    stream, window_duration_sec, average
)

# show covariance from first window and first frequency
covariance_show = np.abs(covariances[0, 0])
fig, ax = plt.subplots(1, constrained_layout=True)
img = ax.imshow(covariance_show, origin="lower", cmap="viridis_r")
ax.set_xticks(range(len(channels)))
ax.set_xticklabels(channels)
ax.set_yticks(range(len(channels)))
ax.set_yticklabels(channels)
ax.set_title("Single-station multiple channels covariance")
plt.colorbar(img, shrink=0.6).set_label("Covariance modulus")
