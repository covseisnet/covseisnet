#!/usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa
"""
Cross correlation network response function 
=========================================================


This example shows how to calculate network response function of a seismic network. Data from the Piton de la Fournaise seismic network is used to illustrate how cross correlation coefficients relate to seismovolcanic activity. A 24 hour period of data from October 14, 2010 is used. This data contains seismovolcanic signals linked to a pre-eruptive seismic swarm and to the beginning of a co-eruptive tremor starting at 15:20 (GMT).


"""

import numpy as np
import matplotlib.pyplot as plt
import covseisnet as csn
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# Download data from the YA Undervolc seismic network with RESIF Seismic data portal
client = Client("RESIF")
signal_duration_sec = 24 * 3600
t = UTCDateTime("2010-10-14T00:00:00.00")

# Download stream for each component
list_stations = [
    "FJS",
    "FLR",
    "FOR",
    "HDL",
    "RVL",
    "SNE",
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

# Downsample stream
stream.decimate(4)

# Merge traces in stream
stream.merge(method=1, fill_value="interpolate", interpolation_samples=-1)

# Synchronize traces in stream
stream = stream.synchronize(t, signal_duration_sec)

# Filter data with bandpass between 1Hz and 10Hz
band = [1.0, 10.0]
stream.filter("bandpass", freqmin=band[0], freqmax=band[1])

# Preprocess data
stream.preprocess()

# Compute the network response function, specify the location of the per-station travel time grids
NRF = csn.crosscorrelation.nrf(stream, T_path="./data")

# Plot the network response function
plt.figure()
plt.plot(np.linspace(0, 24, len(NRF)), NRF, "k")
plt.title("Seismovolcanic activity on October 14, 2010")
plt.ylabel("Network response function, unnormalized")
plt.xlabel("14.10.2010 (hours)")
plt.xlim(0, 24)
