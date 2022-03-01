#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Correlation matrix in spectral domain.

"""
import numpy as np

from scipy.signal import butter, filtfilt, hilbert
from scipy.ndimage import gaussian_filter1d

def cross_correlation(covariance, sampling_rate):
    """Extract correlation in time domain from the given covariance matrix.
    Arguments:
    ----------

    """

    # Extract upper triangular
    covariance = covariance.triu(k=1)

    # Inverse Fourier transform
    correlation = np.real(np.fft.fftshift(
        np.fft.ifft(covariance, axis=-2), axes=-2))

    # Calculate lags (note: omit from code for now as plotting functionality not provided)
    # n_lags = correlation.shape[-2]
    # n_lag_symm = (n_lags - 1) // 2
    # lags = np.arange(-n_lag_symm, n_lag_symm + 1) / sampling_rate

    # return lags, correlation.view(CorrelationMatrix)
    return correlation.view(CorrelationMatrix)

class CorrelationMatrix(np.ndarray):
    """Correlation Matrix.

    """

    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        return obj

    def nwin(self):
        return self.shape[0]

    def hilbert_envelope(self):
        return np.abs(hilbert(self, axis=0)).view(CorrelationMatrix)

    def smooth(self, sigma):
        return gaussian_filter1d(self, sigma, axis=0).view(CorrelationMatrix)

    def bandpass(self, low_cut, high_cut, sampling_rate):
        """Filter the signal within bandwidth."""

        # calculate the Nyquist frequency
        nyquist = 0.5 * sampling_rate
   
        # design filter
        order = 4
        low = low_cut / nyquist
        high = high_cut / nyquist -0.1
        b, a = butter(order, [low, high], btype='band')


        filtered = np.zeros(self.shape)
        for i in range(self.shape[1]):
            filtered[:, i] = filtfilt(b, a, self[:, i])
        return filtered.view(CorrelationMatrix)