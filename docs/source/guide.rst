.. _guide:

User's Guide
============


The `covseisnet` package provides tools for array signal processing, with a focus on data from seismic networks. The central analyzed mathematical construction is the array covariance matrix (sometimes called network covariance matrix). The core signal detection algorithms are based on the analysis of the eigenvalues of this matrix. Eigenvector decomposition provides a basis for a blind source separation. In addition, the covariance matrix can be used as input for classical array processing tools such as beamforming and inter-station cross-correlations. Covseisnet objects are mostly inherited from obspy (seismology-oriented python package) and `numpy`.

The code repository is hosted on GitHub at https://github.com/covseisnet/covseisnet

Roadmap
-------

A brief summary of the evolution of the 'covseisnet' package:

* **version 0.4.1**: introduced the core detection algorithms from the array covariance matrix, where two main objects are defined:

    - :class:`~covseisnet.arraystream.ArrayStream` – a synchronous obspy stream collected from an array of seismic stations (see the :ref:`user-guide-arraystream` section below, or see the object documentation).

    - :class:`~covseisnet.covariancematrix.CovarianceMatrix` – a numpy array with covariance-analysis methods (see the :ref:`user-guide-covariancematrix` section below, or see the object documentation).

* **version 0.5.2**: introduced beamforming and the ability to create more advanced products such as likelihood and network response functions: 

    - :class:`~covseisnet.traveltime.TravelTime` – a list-like object containing traveltime grids loaded from file (see the :ref:`user-guide-traveltime` section below, or see the object documentation).

    - :class:`~covseisnet.correlationmatrix.CorrelationMatrix` – a numpy array with methods to extract correlation in time domain from a given covariance matrix (see the :ref:`user-guide-beam` section below, or see the object documentation).

    - :class:`~covseisnet.beam.Beam` – a class containing methods for calculating the likelihood location and network response function from a correlation matrix (see the :ref:`user-guide-beam` section below, or see the object documentation).

.. _user-guide-arraystream:

Dealing with array seismic data
-------------------------------


Because we do array seismic data analysis, we need synchronized seismic traces. The whole package is based on obspy's :class:`~obspy.core.stream.Stream` object in order to use the seismic and signal-processing tools therein defined. Nevertheless, obspy's :class:`~obspy.core.stream.Stream` allow to gather different traces from different seismic stations with different sampling properties (sampling rate, start time, duration...). Yet, in order to do array signal processing, we need to have synchronous seismic traces, and to make sure that any pre-processing is applied to the whole array seismic data in a similar fashion.

We therefore provide an :class:`~covseisnet.arraystream.ArrayStream` class, which primary goal is to synchronize and pre-process a collection of seismic traces collected at different seismic stations. This class directly inherits from the obspy's :class:`~obspy.core.stream.Stream` class, with four additional methods.

1. Traces synchronization with :meth:`~covseisnet.arraystream.ArrayStream.synchronize`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The :meth:`covseisnet.arraystream.ArrayStream.synchronize` method allows to trim the seismic traces on similar starting and ending dates (and thus, duration) and similar sampling rate. In addition, the method can perform sub-sampling interpolation if the traces are time-shifted below the sampling rate.

2. Traces preprocessing with :meth:`~covseisnet.arraystream.ArrayStream.preprocess`
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The :meth:`covseisnet.arraystream.ArrayStream.preprocess` method provides a great diveristy of pre-processing in the spectral and temporal domains.

3. Traces trimming with :meth:`~covseisnet.arraystream.ArrayStream.cut`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The :meth:`covseisnet.arraystream.ArrayStream.cut` method is a wrapper for the :meth:`~obspy.core.stream.Stream.trim` method; the only difference is that it can work with date strings instead of :class:`~obspy.core.utcdatetime.UTCDateTime` objects.

4. Array seismic data time vector with :meth:`~covseisnet.arraystream.ArrayStream.times`
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In a :class:`~covseisnet.arraystream.ArrayStream` instance, the time vectors of each individual traces is supposed to be the same after the synchronization. Note that all the array operations performed by other classes and methods of the package consider that the traces are synchronous. Therefore, there is only a single time vector that should be considered for all traces. The
:meth:`covseisnet.arraystream.ArrayStream.times` method is a wrapper for the :meth:`obspy.core.trace.Trace.times` method, where only the first seismic station (by default) time vector is considered. This method can return the times in different formats, please check the documentation for more details.

.. _user-guide-covariancematrix:

Network covariance matrix analysis
----------------------------------

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
