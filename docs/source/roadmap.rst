.. _roadmap:

Roadmap
=======

A brief summary of the evolution of the `covseisnet` package:

* **version 0.4.1**: introduced the core detection algorithms from the array covariance matrix, where two main objects are defined:

    - :class:`~covseisnet.arraystream.ArrayStream` – a synchronous obspy stream collected from an array of seismic stations (see the :ref:`user-guide-arraystream` section below, or see the object documentation).

    - :class:`~covseisnet.covariancematrix.CovarianceMatrix` – a numpy array with covariance-analysis methods (see the :ref:`user-guide-covariancematrix` section below, or see the object documentation).

* **version 0.5.3/1.0.0**: introduced beamforming and the ability to create more advanced products such as likelihood and network response functions. Three additional classes were defined: 

    - :class:`~covseisnet.traveltime.TravelTime` – a list-like object containing traveltime grids loaded from file (see the :ref:`user-guide-traveltime` section below, or see the object documentation).

    - :class:`~covseisnet.correlationmatrix.CorrelationMatrix` – a numpy array with methods to extract correlation in time domain from a given covariance matrix (see the :ref:`user-guide-correlationmatrix` section below, or see the object documentation).

    - :class:`~covseisnet.beam.Beam` – a class containing methods for calculating the likelihood location and network response function from a correlation matrix (see the :ref:`user-guide-beam` section below, or see the object documentation).

.. _user-guide-arraystream: