#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Covariance matrix in spectral domain."""

import covnet as cn
import obspy
import matplotlib.pyplot as plt

# import numpy as np

# read ObsPy's example stream
stream = obspy.read()

# calculate covariance from stream
window_duration_sec = 1.0
average = 5
times, frequencies, covariances = cn.covariance.calculate(
    stream, window_duration_sec, average
)

# calculate spectral width
spectral_width = covariances.coherence(kind="spectral_width")
print(spectral_width)

# show covariance at first time window and first frequency
fig, ax = plt.subplots(1, constrained_layout=True)
img = ax.pcolormesh(times, frequencies, spectral_width.T)
ax.set_ylim([0, stream[0].stats.sampling_rate / 2])
ax.set_xlabel("Times (seconds)")
ax.set_ylabel("Frequency (Hz)")
plt.colorbar(img).set_label("Covariance matrix spectral width")
plt.show()
