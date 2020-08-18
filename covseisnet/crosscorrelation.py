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
from scipy import signal
from skimage.util.shape import view_as_windows

def cross_correlation(streamZ, streamN, streamE, window_duration=200, overlap=1, lag_time = 50):
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


    Returns
    -------
    The cross-correlation coefficient matrix.
    """
    
    # TODO: check sampling rate is the same for each component
    # TODO: check stream lengths are the same for each component
    # TODO: check window_duration, overlap, and lagtime are acceptable
    
    nsta = len(streamZ)
    fs = streamZ[0].stats.sampling_rate
    npts = streamZ[0].stats.npts

    
    # Extract data into numpy arrays   
    data_Z = np.zeros((nsta,npts))
    for i in range(nsta):
        data_Z[i] = streamZ[i].data
                
    data_N = np.zeros((nsta,npts))
    for i in range(nsta):
        data_N[i] = streamN[i].data
        
    data_E = np.zeros((nsta,npts))
    for i in range(nsta):
        data_E[i] = streamE[i].data

    # Cut traces into windows
    window_npts = int(window_duration * fs)
    step = int(window_npts*overlap)
    shape = (nsta, window_npts)

    dataslideZ = view_as_windows(data_Z, shape, step = step)
    dataslideZ = dataslideZ[0]
    dataslideN = view_as_windows(data_N, shape, step = step)
    dataslideN = dataslideN[0]
    dataslideE = view_as_windows(data_E, shape, step = step)
    dataslideE = dataslideE[0]
    

    nslide = dataslideZ.shape[0]
    nwin = dataslideZ.shape[2]
    ncc = nwin*2-1
    center = (ncc - 1) // 2 + 1
    lag_npts = int(lag_time*fs)
    
    
    ## Compute inter-components cross-correlations
    
    ccZN = np.zeros([nslide,nsta,ncc])
    ccZE = np.zeros([nslide,nsta,ncc])
    ccNE = np.zeros([nslide,nsta,ncc])
    
    for i in range(0,nslide):
        for j in range(0,nsta):
            ccZN[i,j] = signal.correlate(dataslideZ[i,j],dataslideN[i,j],mode="full",method="fft")
            ccZE[i,j] = signal.correlate(dataslideZ[i,j],dataslideE[i,j],mode="full",method="fft")
            ccNE[i,j] = signal.correlate(dataslideN[i,j],dataslideE[i,j],mode="full",method="fft")
    
    npairs=3  
    cc3 = np.concatenate((ccZN,ccZE,ccNE),axis=1)
    cc3 = cc3.reshape(nslide,npairs,nsta,ncc)
    cc3 =  cc3[:,:,:,center-lag_npts:center+lag_npts]
       
    x = 6
    delta = np.arange(1,x)
    sizelim = nslide - x
    
    ccc = np.zeros([sizelim,npairs,nsta])
    slide_loop = np.linspace(x,nslide-1,sizelim)
    
    
    ### Compute correlation coefficient and average over consecutive windows
    for n in range(0,npairs):
    
        for i,c in enumerate(slide_loop):
            c=int(c)
            for k in range(0,nsta):
              
                tmp = np.zeros((x))
                for j in delta:
                    tmp[j-1] = np.corrcoef(cc3[c,n,k],cc3[c-j,n,k])[0,1]
                ccc[i,n,k] = tmp.mean()
   
    return ccc
               