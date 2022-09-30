.. _guide_correlation:

Correlation in the time domain
==============================

Given a covariance matrix, a corresponding correlation matrix is computed by extracting the upper triangular matrix and performing an inverse fourier transform. A set of lag times is also provided as output so that the amount of shift in the correlation envelopes can be seen.

The following methods are also provided to filter and smoothen the correlation matrix:
* bandpass filter :meth:`~covseisnet.correlationmatrix.CorrelationMatrix.bandpass`
* 1-D gaussian filter :meth:`~covseisnet.correlationmatrix.CorrelationMatrix.smooth`
* hilbert envelope :meth:`~covseisnet.correlationmatrix.CorrelationMatrix.hilbert_envelope`