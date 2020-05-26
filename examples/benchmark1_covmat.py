#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Covariance matrix in spectral domain."""

import covseisnet as csn
import obspy
import matplotlib.pyplot as plt
import numpy as np

# read ObsPy's example stream
stream = obspy.read()

# calculate covariance from stream
window_duration_sec = 1.0
average = 5
times, frequencies, covariances = csn.covariancematrix.calculate(
    stream, window_duration_sec, average
)

# show covariance at first time window and first frequency
fig, ax = plt.subplots(1)
ax.imshow(np.abs(covariances[0, 0]))
plt.show()
