#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa
"""
Cross correlation coefficients 
=========================================================


This example shows how to use cross-correlation coefficients 


"""

import numpy as np
import covseisnet as csn
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# # download data from the YA Undervolc seismic network with RESIF Seismic data portal
client = Client("RESIF")
signal_duration_sec = 24 * 3600
t = UTCDateTime("2010-10-14T00:00:00.00")

streamZ = client.get_waveforms("YA", "UV*", "00", "HHZ", t, t + signal_duration_sec)
streamN = client.get_waveforms("YA", "UV*", "00", "HHN", t, t + signal_duration_sec)
streamE = client.get_waveforms("YA", "UV*", "00", "HHE", t, t + signal_duration_sec)

streamZ = csn.arraystream.ArrayStream(streamZ)
streamN = csn.arraystream.ArrayStream(streamN)
streamE = csn.arraystream.ArrayStream(streamE)


# Downsample stream
streamZ.decimate(4)
streamN.decimate(4)
streamE.decimate(4)

# Merge traces in stream
streamZ.merge(method=1, fill_value="interpolate", interpolation_samples=-1)
streamN.merge(method=1, fill_value="interpolate", interpolation_samples=-1)
streamE.merge(method=1, fill_value="interpolate", interpolation_samples=-1)

# Synchronize traces in stream
streamZ = streamZ.synchronize(
    start=t, duration_sec=signal_duration_sec, method="linear"
)
streamN = streamN.synchronize(
    start=t, duration_sec=signal_duration_sec, method="linear"
)
streamE = streamE.synchronize(
    start=t, duration_sec=signal_duration_sec, method="linear"
)


# Filter data with bandpass between 1Hz and 10Hz
fmin = 1
fmax = 10
streamZ.filter("bandpass", freqmin=fmin, freqmax=fmax, corners=4, zerophase=True)
streamN.filter("bandpass", freqmin=fmin, freqmax=fmax, corners=4, zerophase=True)
streamE.filter("bandpass", freqmin=fmin, freqmax=fmax, corners=4, zerophase=True)


# Preprocess data
streamZ.preprocess()
streamN.preprocess()
streamE.preprocess()

# Run cross-correlation
ccc = csn.crosscorrelation.cross_correlation(streamZ, streamN, streamE)

## average over component pairs
ccc_mean = ccc.mean(axis=1)

## average over all stations
ccc_network = ccc_mean.mean(axis=1)


### plot cc function averaged over all stations
plt.figure()
plt.plot(np.linspace(0, 24, len(ccc_network)), ccc_network, "k")
plt.title("Average over all stations")
plt.ylabel("Correlation coefficient")
plt.xlabel("14.10.2010 (hours)")
