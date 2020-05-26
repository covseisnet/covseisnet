#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Read and pre-process seismic data.

Todo
----
- Implement the new :meth:`~covseisnet.data.ArrayStream.synchronize` method.
- Implement a :meth:`covseisnet.data.ArrayStream.check_synchronicity` method.
"""

import obspy
import numpy as np

from scipy import signal


class ArrayStream(obspy.core.stream.Stream):
    """List-like object of multiple ObsPy trace objects (synchroneous).

    This class is a subclass of the :class:`obspy.core.stream.Stream` class,
    with additional methods pre-processing methods. The main idea is to
    gather traces with an exact same amount of samples per traces in order
    to perform array-procesing methods onto it. The synchronization of the
    different traces objects is not automatic. There are two options for
    synchronizing the stream:

    (1) use the :meth:`~covseisnet.data.ArrayStream.synchronize`

        >>> import covseisnet as cn
        >>> stream = cn.data.read()
        >>> stream.synchronize()

    (2) perform a manual synchronization before turning the stream
        into an :class:`covseisnet.ArrayStream` object with

        >>> import covseisnet as cn
        >>> stream = obspy.read()
        >>> # manually synchronize
        >>> stream = cn.data.ArrayStream(stream)

    Note
    -----
    All the original methods of :class:`obspy.core.stream.Stream` objects
    remain available in :class:`~covseisnet.data.ArrayStream` objects. For more
    information, please visit the ObsPy documentation at
    https://examples.obspy.org.
    """

    def __init__(self, *args, **kwargs):
        r"""Subclassing."""
        super(ArrayStream, self).__init__(*args, **kwargs)

    def cut(self, starttime, endtime, **kwargs):
        r"""Cut (trim) seismic traces between given start and end times.

        This function is a wrapper to the :meth:`obspy.core.stream.Stream.trim`
        method, but works directly with datetimes in :class:`str` format.
        The function uses the native ObsPy
        :class:`~obspy.core.utcdatetime.UTCDateTime` class in order
        to convert the datetimes from :class:`str` into
        :class:`obspy.core.utcdatetime.UTCDateTime` format.

        Parameters
        ----------
        starttime : str
            The starting date time. Can be in any format that the
            :class:`~obspy.core.utcdatetime.UTCDateTime` can parse.
            For instance ``starttime="2010-01-01 10:34"``.

        endtime : str
            The ending date time, with the same constrain that the
            ``starttime`` parameter.

        Keyword arguments
        -----------------
        **kwargs: dict, optional
            Additional keyword arguments are passed to the
            :meth:`obspy.core.stream.Stream.trim` method.

        Example
        -------

        >>> stream = cn.data.read()
        >>> stream.cut('2009-08-24 00:20:05', '2009-08-24 00:20:12')
        >>> print(stream)
        3 Trace(s) in Stream:
        BW.RJOB..EHZ | 2009-08-24T00:20:05.000000Z... | 100.0 Hz, 701 samples
        BW.RJOB..EHN | 2009-08-24T00:20:05.000000Z... | 100.0 Hz, 701 samples
        BW.RJOB..EHE | 2009-08-24T00:20:05.000000Z... | 100.0 Hz, 701 samples
        """
        starttime = obspy.UTCDateTime(starttime)
        endtime = obspy.UTCDateTime(endtime)
        self.trim(starttime, endtime, **kwargs)

    def normalize(
        self, method="onebit", smooth_length=None, smooth_order=None, epsilon=1e-10
    ):
        r"""Normalize the seismic traces in temporal domain.

        Considering :math:`x_i(t)` being the seismic trace :math:`x_i(t)`, the
        normalized trace :math:`\tilde{x}_i(t)` is obtained with

        .. math::
            \tilde{x}_i(t) = \frac{x_i(t)}{Fx_i(t) + \epsilon}

        where :math:`Fx` is a characteristic of the trace :math:`x` that
        depends on the ``method`` argument, and :math:`\epsilon > 0` is a
        regularization value to avoid division by 0, set by the ``epsilon``
        keyword argument.

        Keyword arguments
        -----------------
        method : str, optional
            Must be one of "onebit" (default), "mad", or "smooth".

            - "onebit" compress the seismic trace into a series of 0 and 1.
              In this case, :math:`F` is defined as :math:`Fx(t) = |x(t)|`.

            - "mad" normalize each trace by its median absolute deviation.
              In this case, :math:`F` delivers a scalar value defined as
              :math:`Fx(t) = \text{MAD}x(t) =
              \text{median}(|x(t) - \langle x(t)\rangle|)`, where
              :math:`\langle x(t)\rangle)` is the signal's average.

            - "smooth" normalize each trace by a smooth version of its
              envelope. In this case, :math:`F` is obtained from the
              signal's Hilbert envelope.

        smooth_length: int, optional
            If the ``method`` keyword argument is set to "smooth", the
            normalization is performed with the smoothed trace envelopes,
            calculated over a sliding window of `smooth_length` samples.
            Note that this parameter is not consider if ``method`` is not
            set to "smooth". (`None` by default)

        smooth_order: int, optional
            If the ``method`` keyword argument is set to "smooth", the
            normalization is performed with the smoothed trace envelopes.
            The smoothing order is set by the ``smooth_order`` parameter.
            Note that this parameter is not consider if ``method`` is not
            set to "smooth". (`None` by default)

        epsilon: float, optional
            Regularization parameter in division, set to ``1e-10`` by default.

        """
        if method == "onebit":
            for trace in self:
                trace.data = trace.data / (np.abs(trace.data) + epsilon)

        if method == "smooth":
            for trace in self:
                trace_env_smooth = signal.savgol_filter(
                    np.abs(trace.data), smooth_length, smooth_order
                )
                trace.data = trace.data / (trace_env_smooth + epsilon)

        if method == "demad":
            for trace in self:
                trace.data /= signal.mad(trace.data) + epsilon

    def synchronize(
        self,
        sampling_rate=20.0,
        method="linear",
        start="2010-01-01",
        npts=24 * 3600 * 20,
    ):
        r"""Synchronize seismic traces into the same times.

        So far, this function uses the
        :meth:`obspy.core.trace.Trace.interpolate` method in order to
        synchronize the traces. This is not yet optimal (see to-do list).

        Keyword arguments
        -----------------
        sampling_rate : float, optional
            The desired final sampling rate in Hz. Default to 20 Hz.

        method: str, optional
            Interpolation method. Default to "linear".

        start: str, optional
            Start date of the interpolation. Default to "2010-01-01".

        npts: int, default
            Number of samples per traces. Default to 1,728,000 samples, the
            total number of samples on a single day at 20 Hz.


        Warning
        -------
        This function is deprecated, and wil be replaced in near future.
        """
        start = obspy.UTCDateTime(start)
        self.interpolate(sampling_rate, method, start, npts)

    def times(self, station_index=0, **kwargs):
        r"""Common time vector of the ArrayStream.

        Because the :class:`~covseisnet.data.ArrayStream` is supposed to handle
        traces with exact same number of samples (and sampling frequency), the
        time vector of each traces is supposed to be the same. This function
        only returns the times of one of the traces, accessible from the
        :meth:`obspy.core.trace.Trace.times` method.

        Keyword arguments
        -----------------

        station_index: int, optional
            The trace index to extract the time vector from. This has no
            influence on the returned time vector is the traces have indeed
            the same sampling, otherwise, you should consider synchronizing
            the traces first. By default, the first trace is considered.

        **kwargs: dict, optional
            Additional keyword arguments are directly passed to the
            :meth:`obspy.core.trace.Trace.times` (for instance,
            ``type="matplotlib"`` allows to recover matplotlib timestamps
            provided by the :func:`matplotlib.dates.date2num` function.

        Returns
        -------
        :class:`numpy.ndarray` or :class:`list`.
            An array of timestamps in a :class:`numpy.ndarray` or in a
            :class:`list`.

        Note
        ----
        If the times are not synchroneous, the time vector will only correspond
        to the times of the trace indexed with ``station_index``. The user
        should ensure that the traces are synchroneous first.

        Tip
        ---
        In order to extract times in matplotlib format, you can set the
        ``type`` parameter of the
        :meth:`~obspy.core.trace.Trace.times` method such as

        >>> import covseisnet as cn
        >>> st = cn.data.read()
        >>> st.times(type='matplotlib')
        array([ 733643.01392361,  733643.01392373,  733643.01392384, ...,
        733643.01427049,  733643.0142706 ,  733643.01427072])
        """
        return self[0].times(**kwargs)

    def whiten(self, segment_duration_sec, method="pure", smooth=11):
        """Spectral normalization of the traces.

        Parameters
        ----------
        segment_duration_sec : float
            Duration of the segments for Fourier transformation.

        Keyword arguments
        -----------------
        method : str
            ``"pure"`` or ``"smooth"``. Wheter to consider the division with
            direct Fourier transform modulus, or a smooth version.

        smooth : int
            Smoothing window length in points.

        """
        # Define method
        if method == "pure":
            whiten_method = signal.phase
        elif method == "smooth":
            whiten_method = signal.detrend_spectrum

        # Initialize for waitbar
        # waitbar = logtable.waitbar('Whiten', len(self))
        fft_size = int(segment_duration_sec * self[0].stats.sampling_rate)
        duration = self[0].times()[-1]

        # Whiten
        for index, trace in enumerate(self):
            # waitbar.progress(index)
            data = trace.data
            _, _, data_fft = signal.stft(data, nperseg=fft_size)
            data_fft = whiten_method(data_fft, smooth=smooth)
            _, data = signal.istft(data_fft, nperseg=fft_size)
            trace.data = data

        # Trim
        self.cut(
            pad=True,
            fill_value=0,
            starttime=self[0].stats.starttime,
            endtime=self[0].stats.starttime + duration,
        )


def read(pathname_or_url=None, **kwargs):
    """Read seismic waveforms files into an ArrayStream object.

    This function uses the :func:`obspy.core.stream.read` function to read
    the streams. A detailed list of arguments and options are available at
    https://docs.obspy.org. This function opens either one or multiple
    waveform files given via file name or URL using the ``pathname_or_url``
    attribute. The format of the waveform file will be automatically detected
    if not given. See the `Supported Formats` section in
    the :func:`obspy.core.stream.read` function.

    This function returns an :class:`~covseisnet.data.ArrayStream` object, an
    object directly inherited from the :class:`obspy.core.stream.Stream`
    object.

    Keyword arguments
    -----------------
    pathname_or_url: str or io.BytesIO or None
        String containing a file name or a URL or a open file-like object.
        Wildcards are allowed for a file name. If this attribute is omitted,
        an example :class:`~covseisnet.data.ArrayStream` object will be returned.

    Other parameters
    ----------------
    **kwargs: dict
        Other parameters are passed to the :func:`obspy.core.stream.read`
        directly.

    Returns
    -------
    :class:`~covseisnet.data.ArrayStream`
        An :class:`~covseisnet.data.ArrayStream` object.

    Example
    -------

    In most cases a filename is specified as the only argument to
    :func:`obspy.core.stream.read`. For a quick start you may omit all
    arguments and ObsPy will create and return a basic example seismogram.
    Further usages of this function can be seen in the ObsPy documentation.

    >>> import covseisnet as cn
    >>> stream = cn.data.read()
    >>> print(stream)
    3 Trace(s) in Stream:
    BW.RJOB..EHZ | 2009-08-24T00:20:03.000000Z - ... | 100.0 Hz, 3000 samples
    BW.RJOB..EHN | 2009-08-24T00:20:03.000000Z - ... | 100.0 Hz, 3000 samples
    BW.RJOB..EHE | 2009-08-24T00:20:03.000000Z - ... | 100.0 Hz, 3000 samples

    .. rubric:: _`Further Examples`

    Example waveform files may be retrieved via https://examples.obspy.org.
    """
    stream = obspy.read(pathname_or_url, **kwargs)
    stream = ArrayStream(stream)
    return stream
