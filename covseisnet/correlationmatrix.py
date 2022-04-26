#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Correlation matrix in spectral domain.

"""
import numpy as np

from scipy.signal import butter, filtfilt, hilbert
from scipy.ndimage import gaussian_filter1d


class CorrelationMatrix(np.ndarray):
    r"""Correlation Matrix.
    
    This class is a subclass of :class:`numpy.ndarray`.
    """

    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        return obj

    def nwin(self):
        """ Returns the number of windows in the correlation matrix.
        
        Returns
        -------
        int
            The number of windows in the correlation matrix.
        
        """
        return self.shape[0]

    def hilbert_envelope(self, **kwargs):
        """Apply the Hilbert transform to the correlation matrix. Uses
        :func:`~scipy.signal.hilbert`
             
        """
        return np.abs(hilbert(self, axis=0, **kwargs)).view(CorrelationMatrix)

    def smooth(self, sigma, **kwargs):
        """Apply a 1-D Gaussian filter to the correlation matrix. Uses
        :func:`~scipy.ndimage.gaussian_filter1d`.
        
        Parameters
        ----------
        
        sigma: float
            Standard deviation for Gaussian kernel
                   
        """
        return gaussian_filter1d(self, sigma, axis=0, **kwargs).view(CorrelationMatrix)

    def bandpass(self, low_cut, high_cut, sampling_rate, **kwargs):
        """Apply a Butterworth bandpass filter to the correlation matrix. Uses
        :func:`~scipy.signal.butter` and :func:`~scipy.signal.filtfilt`.
        
        Parameters
        ----------
        
        low_cut: float
            Pass band low corner frequency.
            
        high_cut: float
            Pass band high corner frequency.
            
        sampling_rate: float
            Sampling rate in Hz.        
            
        """
        # calculate the Nyquist frequency
        nyquist = 0.5 * sampling_rate

        # design filter
        order = 4
        low = low_cut / nyquist
        high = high_cut / nyquist - 0.1
        b, a = butter(order, [low, high], btype="band", **kwargs)

        filtered = np.zeros(self.shape)
        for i in range(self.shape[1]):
            filtered[:, i] = filtfilt(b, a, self[:, i], **kwargs)
        return filtered.view(CorrelationMatrix)


def cross_correlation(covariance, sampling_rate):
    """Extract correlation in time domain from the given covariance matrix.
    
    Parameters
    ----------
    
    covariance: :class:`covseisnet.covariancematrix.CovarianceMatrix` object
        A covariance matrix already computed.
        
    sampling_rate: float
        Sampling rate in Hz.
        

    Returns
    -------
    :class:`~tuple`
    
        :class:`~numpy.ndarray` The lag time between stations.
    
        :class:`~covseisnet.correlationmatrix.CorrelationMatrix` The correlation matrix.
        
    """

    # Extract upper triangular
    covariance = covariance.triu(k=1)

    # Inverse Fourier transform
    correlation = np.real(np.fft.fftshift(np.fft.ifft(covariance, axis=-2), axes=-2))

    # Calculate lags
    n_lags = correlation.shape[-2]
    n_lag_symm = (n_lags - 1) // 2
    lags = np.arange(-n_lag_symm, n_lag_symm + 1) / sampling_rate

    correlation_unshifted = correlation.view(CorrelationMatrix)
    return lags, correlation_unshifted
