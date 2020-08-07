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

from functools import partial
from scipy import signal, stats


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

    def preprocess(self, domain="spectral", **kwargs):
        r"""Pre-process each trace in temporal or spectral domain."""
        kwargs.setdefault("epsilon", 1e-10)
        if domain == "spectral":
            whiten(self, **kwargs)
        elif domain == "temporal":
            normalize(self, **kwargs)
        else:
            raise ValueError(
                "Invalid preprocessing domain {} - please specify 'spectral' or 'temporal'".format(
                    domain
                )
            )
        pass

    def synchronize(
        self, start="2010-01-01T00:00:00.00", duration_sec=24 * 3600, method="linear"
    ):
        r"""Synchronize seismic traces into the same times.

        This function uses the
        :meth:`obspy.core.stream.Stream.trim` method in order to cut all
        traces of the Stream object to given start and end time with
        padding if necessary, and the
        :meth:`obspy.core.trace.Trace.interpolate` method in order to
        synchronize the traces. Then data are linearly interpolated.

        Keyword arguments
        -----------------
        start: str, default
            Start date of the interpolation.
            Default to "2010-01-01T00:00:00.00".

        duration_sec: int, default
            Duration of the data traces in second. Default to 86,400 seconds,
            the total number of seconds on a single day.

        method: str, default
            Interpolation method. Default to "linear".
        """
        # Duplicate stream
        stream_i = self.copy()
        start = obspy.UTCDateTime(start)
        end = start + duration_sec
        sampling = self[0].stats.sampling_rate
        npts = int(sampling * duration_sec)
        stream_i.trim(start, end, pad=True, fill_value=0, nearest_sample=False)

        for tr, tr_i in zip(self, stream_i):
            t = tr.times() / duration_sec + tr.stats.starttime.matplotlib_date
            shift = tr_i.stats.starttime - start
            tr_i.interpolate(
                sampling, method, start=start, npts=npts, time_shift=-shift
            )
            ti = tr_i.times() / duration_sec + tr_i.stats.starttime.matplotlib_date
            tr_i.data = np.interp(ti, t, tr.data)

        return stream_i

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
        an example :class:`~covseisnet.data.ArrayStream` object will be
        returned.

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


def whiten(
    stream,
    method="onebit",
    window_duration_sec=2,
    smooth_length=11,
    smooth_order=1,
    epsilon=1e-10,
):
    r"""Normalize in the spectral domain."""
    if method == "onebit":
        whiten_method = phase
    elif method == "smooth":
        whiten_method = partial(
            detrend_spectrum, smooth=smooth_length, order=smooth_order, epsilon=epsilon
        )
    else:
        raise ValueError("Unknown method {}".format(method))
    r"""Whiten traces in the spectral domain."""
    fft_size = int(window_duration_sec * stream[0].stats.sampling_rate)
    for index, trace in enumerate(stream):
        data = trace.data
        _, _, data_fft = signal.stft(data, nperseg=fft_size)
        data_fft = whiten_method(data_fft)
        _, data = signal.istft(data_fft, nperseg=fft_size)
        trace.data = data
    pass


def detrend_spectrum(x, smooth=None, order=None, epsilon=1e-10):
    r"""Smooth modulus spectrum.

    Arugments
    --------
    x: :class:`np.ndarray`
        The spectra to detrend. Must be of shape `(n_frequencies, n_times)`.

    smooth: int
        Smoothing window size in points.

    order: int
        Smoothing order. Please check the :func:`savitzky_golay` function
        for more details.

    Keyword arguments
    -----------------
    epsilon: float, optional
        A regularizer for avoiding zero division.

    Returns
    -------
    The spectrum divided by the smooth modulus spectrum.
    """
    n_frequencies, n_times = x.shape
    for t in range(n_times):
        x_smooth = signal.savgol_filter(np.abs(x[:, t]), smooth, order)
        x[:, t] /= x_smooth + epsilon
    return x


def normalize(stream, method="onebit", smooth_length=11, smooth_order=1, epsilon=1e-10):
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


    smooth_order: int, optional
        If the ``method`` keyword argument is set to "smooth", the
        normalization is performed with the smoothed trace envelopes.
        The smoothing order is set by the ``smooth_order`` parameter.


    epsilon: float, optional
        Regularization parameter in division, set to ``1e-10`` by default.

    """
    if method == "onebit":
        for trace in stream:
            trace.data = trace.data / (np.abs(trace.data) + epsilon)

    elif method == "smooth":
        for trace in stream:
            trace_env_smooth = signal.savgol_filter(
                np.abs(trace.data), smooth_length, smooth_order
            )
            trace.data = trace.data / (trace_env_smooth + epsilon)

    elif method == "mad":
        for trace in stream:
            trace.data = trace.data / (
                stats.median_absolute_deviation(trace.data) + epsilon
            )

    else:
        raise ValueError("Unknown method {}".format(method))


def phase(x):
    r"""Complex phase extraction.

    Given a complex number (or complex-valued array)
    :math:`x = r e^{\imath \phi}`, where :math:`r` is the complex modulus
    and :math:`phi` the complex phase, the function returns the unitary-modulus
    complex number such as

    .. math::

              \tilde{x} = e^{\imath \phi}

    Arguments
    ---------
    x: :class:`np.ndarray`
        The complex-valued data to extract the complex phase from.
    """
    return np.exp(1j * np.angle(x))
