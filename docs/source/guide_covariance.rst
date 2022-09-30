.. _guide_covariance:

Network covariance matrix analysis
==================================

One of the goal of the package is to provide detection strategies based on the properties of the spectral covariance matrix of the array seismic data. The spectral network covariance matrix is the Fourier transform of the time-domain (local) inter-station cross-correlation matrix. The :class:`~covseisnet.covariancematrix.CovarianceMatrix` object is based on a :class:`numpy.ndarray` with additional covariance-based analysis tools. One should never instanciate a :class:`~covseisnet.covariancematrix.CovarianceMatrix` object directly, but calculate it from an :class:`~covseisnet.arraystream.ArrayStream` object (or an obspy's :class:`~obspy.core.stream.Stream` directly if the user ensure that the traces are already synchronous and pre-processed) with the :class:`covseisnet.covariancematrix.calculate` function.

The shape of a :class:`~covseisnet.covariancematrix.CovarianceMatrix` object calculated from :math:`N` traces is at least :math:`N \times N`. Depending on the averaging size and frequency content, the covariance matrix can be of shape

- ``(n_sta, n_sta)`` if a single frequency and time sample is given.

- ``(n_freq, n_sta, n_sta)`` for a single time sample and ``n_freq`` frequency points.

- ``(n_times, n_freq, n_sta, n_sta)`` for ``n_times`` and ``n_freq`` dimensions.

All the methods defined in the the :class:`~arrayprocessing.covariance.CovarianceMatrix` class are performed on the flattened array with the private method :class:`arrayprocessing.covariance.CovarianceMatrix._flat`, which allows to obtain as many :math:`N \times N` covariance matrices as time and frequency samples.

1. Covariance matrix estimation from an :class:`~covseisnet.arraystream.ArrayStream` object with :func:`covseisnet.covariancematrix.calculate`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The function :func:`covseisnet.covariancematrix.calculate` allows to compute the spectral network covariance matrix from a synchronous stream object (a manually synchronized obspy's :class:`~obspy.core.stream.Stream` object, or an :class:`~covseisnet.arraystream.ArrayStream` from this package). This function makes use of the :func:`covseisnet.covariancematrix.stft` function to calculate the Fourier spectra, and of :func:`covseisnet.covariancematrix.xcov` to build the covariance matrix.

In order to estimate the covariance matrix, two main parameters are of importance: (1) the `window_duration_sec` which defines the duration of Fourier spectral segments; and (2) `average` which is the number of consecutive Fourier segments to average in order to estimate the covariance. By default, the Fourier spectral segments and the final averaged window are both overlaped by 50%.

2. Wavefield coherence with :meth:`~covseisnet.covariancematrix.CovarianceMatrix.coherence`
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The spatial coherence is a well-defined concept for continuous wavefields. Indeed, it is related to the number of coefficients required to decompose the observed wavefield onto a basis. When the wavefield is recorded at discrete spatial locations (seismic stations), the concept of wavefield coherence can still be assessed from the covariance matrix eigenstructure. In particular, we provide the covariance matrix spectral width for assessing the spatial coherence.

The **covariance matrix spectral width** is a real positive scalar number which measures the width of the network covariance matrix eigenvalues distribution. This measurement can be represented in a time and frequency diagram.

..
	The **Shannon entropy** provides a measurement of the degree of information present in a given dataset. Extended to the case of discrete operators by Van Neumann, it can also be calculated from the eigenvalue distribution.

This coherence measurement is delivered by the :meth:`~covseisnet.covariancematrix.CovarianceMatrix.coherence` method (see the documentation of the method for more details). Please visit the :ref:`examples` page for examples.


3. direct eigenvalue assessment with :meth:`~covseisnet.covariancematrix.CovarianceMatrix.eigenvalues`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Other measurments of the wavefield coherence can be also implemented manually by the user from the eigenvalue distribution (for instance polarization analysis in the case of 3-component single-station data). The method :meth:`covseisnet.covariancematrix.CovarianceMatrix.eigenvalues` allows to directly extract the eigenvalues of the covariance matrices collected at different times and frequencies.

4. direct eigenvector assessment with :meth:`~covseisnet.covariancematrix.CovarianceMatrix.eigenvectors`
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Many source-separation algorithms are based on the eigenvectors of the network covariance matrix. We therefore provide a method :meth:`covseisnet.covariancematrix.CovarianceMatrix.eigenvectors` to access it from the covariance matrix.

.. _user-guide-traveltime: