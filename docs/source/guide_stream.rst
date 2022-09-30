.. _guide_deal:

Dealing with array seismic data
===============================

Because we do array seismic data analysis, we need synchronized seismic traces. The whole package is based on obspy's :class:`~obspy.core.stream.Stream` object in order to use the seismic and signal-processing tools therein defined. Nevertheless, obspy's :class:`~obspy.core.stream.Stream` allow to gather different traces from different seismic stations with different sampling properties (sampling rate, start time, duration...). Yet, in order to do array signal processing, we need to have synchronous seismic traces, and to make sure that any pre-processing is applied to the whole array seismic data in a similar fashion.

The ArrayStream
+++++++++++++++

We therefore provide an :class:`~covseisnet.arraystream.ArrayStream` class, which primary goal is to synchronize and pre-process a collection of seismic traces collected at different seismic stations. This class directly inherits from the obspy's :class:`~obspy.core.stream.Stream` class, with four additional methods.

Traces synchronization
++++++++++++++++++++++

The :meth:`covseisnet.arraystream.ArrayStream.synchronize` method allows to trim the seismic traces on similar starting and ending dates (and thus, duration) and similar sampling rate. In addition, the method can perform sub-sampling interpolation if the traces are time-shifted below the sampling rate.

Traces preprocessing
++++++++++++++++++++

The :meth:`covseisnet.arraystream.ArrayStream.preprocess` method provides a great diveristy of pre-processing in the spectral and temporal domains.

Traces trimming
+++++++++++++++

The :meth:`covseisnet.arraystream.ArrayStream.cut` method is a wrapper for the :meth:`~obspy.core.stream.Stream.trim` method; the only difference is that it can work with date strings instead of :class:`~obspy.core.utcdatetime.UTCDateTime` objects.

Array seismic data time vector
++++++++++++++++++++++++++++++

In a :class:`~covseisnet.arraystream.ArrayStream` instance, the time vectors of each individual traces is supposed to be the same after the synchronization. Note that all the array operations performed by other classes and methods of the package consider that the traces are synchronous. Therefore, there is only a single time vector that should be considered for all traces. The
:meth:`covseisnet.arraystream.ArrayStream.times` method is a wrapper for the :meth:`obspy.core.trace.Trace.times` method, where only the first seismic station (by default) time vector is considered. This method can return the times in different formats, please check the documentation for more details.