#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compute cross-correlation in the time domain between component pairs for each station. 
A correlation coefficient for consecutive windows is computed and an average taken of
consecutive computations of this correlation coefficient. If there is a seismic source 
stable in time, then the cross-correlation will be stable in time too and the average 
of the correlation coefficient will be higher. 



Todo
----
"""

import numpy as np
import scipy
from scipy import signal
from skimage.util.shape import view_as_windows


def average_windowing(streams, window_duration=200, overlap=0.5):

    ncomp = len(streams)
    nsta = len(streams[0])
    fs = streams[0][0].stats.sampling_rate
    npts = streams[0][0].stats.npts

    # Cut traces into windows
    window_npts = int(window_duration * fs)
    step = int(window_npts * overlap)
    shape = (nsta, window_npts)

    # check sampling rate is the same for all streams
    if not all(
        fs == st[i].stats.sampling_rate for st in streams for i in range(len(st))
    ):
        raise ValueError("Not all streams have the same sampling rate.")

    # check every component has the same number of streams
    if not all(nsta == len(streams[i]) for i in range(len(streams))):
        raise ValueError(
            "Not all streams have the same length (ie. same number of stations)"
        )

    # TODO: check window_duration, overlap, and lagtime are acceptable

    # Extract data into numpy arrays
    # new method
    data = [np.zeros((nsta, npts)) for x in range(ncomp)]
    dataslide = []

    for i in range(0, ncomp):

        for j in range(nsta):
            data[i][j] = streams[i][j].data

        # perform windowing
        dataslide.append(view_as_windows(data[i], shape, step=step)[0])

    return dataslide


def nrf(stream, window_duration=200, overlap=0.5, lag_time=50, gaussian_window=2):

    ncomp = len(stream)
    nsta = len(stream)
    fs = stream[0].stats.sampling_rate

    dataslide = average_windowing([stream], window_duration, overlap)
    dataslide = dataslide[0]

    nslide = dataslide.shape[0]
    nwin = dataslide.shape[2]
    ncc = nwin * 2 - 1
    center = (ncc - 1) // 2 + 1
    lag_npts = int(lag_time * fs)
    npairs = int(ncomp * (ncomp - 1) / 2)

    # Method 3 - NRF
    cc = np.zeros([nslide, npairs, lag_npts * 2])

    for i in range(0, nslide):
        n = 0
        while n < npairs:
            for a in range(0, nsta):
                for b in range(a + 1, nsta):

                    crosscorr = signal.correlate(
                        dataslide[i, b], dataslide[i, a], mode="full", method="fft"
                    )
                    cc[i, n] = crosscorr[center - lag_npts : center + lag_npts]
                    n = n + 1

    gaussian_window_npts = gaussian_window * fs
    c = np.zeros((nslide, npairs, lag_npts * 2))
    env = np.zeros([nslide, npairs, lag_npts * 2])
    for i in range(0, nslide):
        for j in range(0, npairs):
            c[i, j] = np.abs(signal.hilbert(cc[i][j]))
            env[i, j] = scipy.ndimage.filters.gaussian_filter1d(
                c[i, j], gaussian_window_npts, mode="constant"
            )

    # TODO: figure out best way to store travel time files

    # # get list of stations in stream
    # list_avail = [stream[i].stats.station for i in range(nsta)]
    # list_stations = [
    #     "FJS",
    #     "FLR",
    #     "FOR",
    #     "HDL",
    #     "RVL",
    #     "SNE",
    #     "UV01",
    #     "UV02",
    #     "UV03",
    #     "UV04",
    #     "UV05",
    #     "UV06",
    #     "UV07",
    #     "UV08",
    #     "UV09",
    #     "UV10",
    #     "UV11",
    #     "UV12",
    #     "UV13",
    #     "UV14",
    #     "UV15",
    # ]

    # for tr in range(0,len(stream)):
    #     list_avail.append(stream[tr].stats.station)

    T = np.load("travel_time_PdF.npy")
    nlon = T.shape[1]
    nlat = T.shape[2]
    ndep = T.shape[3]

    # T = []
    # for i,sta in enumerate(list_stations):
    #     boolean = sta in list_avail
    #     if boolean == True:
    #       T.append(T_stations[i])
    # T = np.array((T))

    delta_t = np.zeros((npairs, nlon, nlat, ndep))

    n = 0
    while n < npairs:
        for a in range(0, nsta):
            for b in range(a + 1, nsta):
                delta_t[n] = T[a] - T[b]
                n = n + 1

    delta_index = delta_t * fs

    delta_int = np.round(delta_index)
    delta_int = delta_int.astype(np.int64)

    center = (env.shape[2] - 1) // 2 + 1
    R = np.zeros([nslide, nlon, nlat, ndep])

    for i in range(0, nslide):
        for j in range(0, npairs):
            R[i] = R[i] + env[i, j, center - delta_int[j]]

    NRF = np.zeros((nslide))

    ## Take reference as average of (R_max - R_min) during tremor
    start_tremor = int(nslide / 24 * 15) + 1
    Rdiff_tremor = []
    for i in range(start_tremor, len(NRF)):
        Rdiff_tremor.append(R[i].max() - R[i].min())

    R_ref = np.mean(Rdiff_tremor)

    for i in range(0, nslide):
        NRF[i] = R[i].max() - R[i].min()
    NRF = 100 * NRF / R_ref

    return NRF


# Method: Compute inter-components cross-correlation coefficients
def cross_correlation(
    *streams, window_duration=200, overlap=0.5, lag_time=50, consec_win=6
):
    r"""Calculate cross correlation coefficients from three streams.

    Arguments
    --------
    streamZ: :class:`ArrayStream`
        The input data stream for the Z component.
        
    streamN: :class:`ArrayStream`
        The input data stream for the N component.

    streamE: :class:`ArrayStream`
        The input data stream for the E component.
        
    window_duration: int
        Width of the sliding window in seconds. If unspecified, a default of 200 seconds is used.

    overlap: int
        Amount of overlap between windows in seconds. If unspecified, a default of 1 second is used.
        
    lagtime: int
        Amount of lag in seconds. If unspecified, a default of 50 seconds is used.

    consec_win
        Number of consecutive windows used to average the correlation coefficients. If unspecified, a default of 6 is used.
        
    Returns
    -------
    The cross-correlation coefficient matrix.
    """

    # TODO check that more than one stream provided
    ncomp = len(streams)
    nsta = len(streams[0])
    fs = streams[0][0].stats.sampling_rate

    dataslide = average_windowing(streams, window_duration, overlap)

    nslide = dataslide[0].shape[0]
    nwin = dataslide[0].shape[2]
    ncc = nwin * 2 - 1
    center = (ncc - 1) // 2 + 1
    lag_npts = int(lag_time * fs)
    npairs = int(ncomp * (ncomp - 1) / 2)

    cc = {}

    # initialize dictionary of pairs
    for i in range(0, ncomp):
        for j in range(i + 1, ncomp):
            cc[i, j] = np.zeros([nslide, nsta, ncc])

    for i in range(0, nslide):
        for j in range(0, nsta):
            for pair in cc:
                cc[pair][i, j] = signal.correlate(
                    dataslide[pair[0]][i, j],
                    dataslide[pair[1]][i, j],
                    mode="full",
                    method="fft",
                )

    ccx = np.concatenate(list(cc.values()), axis=1)
    ccx = ccx.reshape(nslide, npairs, nsta, ncc)
    ccx = ccx[:, :, :, center - lag_npts : center + lag_npts]

    delta = np.arange(1, consec_win)
    sizelim = nslide - consec_win

    ccc = np.zeros([sizelim, npairs, nsta])
    slide_loop = np.linspace(consec_win, nslide - 1, sizelim)

    ### Compute correlation coefficient and average over consecutive windows
    for n in range(0, npairs):

        for i, c in enumerate(slide_loop):
            c = int(c)
            for k in range(0, nsta):

                tmp = np.zeros((consec_win))
                for j in delta:
                    tmp[j - 1] = np.corrcoef(ccx[c, n, k], ccx[c - j, n, k])[0, 1]
                ccc[i, n, k] = tmp.mean()

    return ccc
